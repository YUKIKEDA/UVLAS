#include "../support/allocation_audit.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <uvlas/numerics/dense_solve.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;

using Complex = std::complex<double>;
constexpr auto operation = core::Operation::SolveLinear;
Context context() {
    return {{core::ModelFingerprint{}, core::ModelGeneration{1}, core::GeometryGeneration{2},
             core::ChartGeneration{3}, core::SparsityGeneration{4}, core::NumericalGeneration{5},
             core::ScaleGeneration{6}, core::PrecisionGeneration{7}, core::Precision::Binary64,
             core::WakePolicyGeneration{8}, core::TimeIntervalGeneration{9},
             core::BackendGeneration{10}},
            core::TrialId{11},
            core::ScopeId{12},
            core::ScopeId{13},
            0.25};
}
template <StorageScalar T>
class NoChecks final : public ConsumerEvaluator<T> {
public:
    [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
        return SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<T>&, const PrecisionBuffer<T>&, const Context&,
           const ConsumerCheck& check, core::CertificateContext c) const override {
        FAIL_CHECK("algebraic solve must not call an absent physical requirement");
        return core::ErrorCertificate::unevaluated(c, check.budget);
    }
};

template <StorageScalar T, std::size_t N>
DenseMatrix<T> matrix(core::Index rows, core::Index columns, const std::array<T, N>& values,
                      core::ResourcePlan& plan) {
    auto value = DenseMatrix<T>::copy({rows, columns, rows}, values, plan, operation);
    REQUIRE(value);
    return std::move(*value);
}
template <StorageScalar T>
const DenseMatrix<T>& high(const PrecisionBuffer<T>& values) {
    if (auto* b = std::get_if<Binary64<T>>(&values.storage()))
        return b->values;
    return std::get<TwoWord<T>>(values.storage()).hi;
}
template <StorageScalar T>
T low(const PrecisionBuffer<T>& values, std::size_t i) {
    if (auto* b = std::get_if<TwoWord<T>>(&values.storage()))
        return b->lo.values()[i];
    return {};
}
struct Setup {
    core::ResourcePlan plan{1 << 24};
    std::array<double, 3> s{1, 2, 4};
    std::array<UnitDimensions, 3> units{};
    [[nodiscard]] Scales scales(std::size_t n) const {
        return {std::span{s}.first(n), std::span{s}.first(n), std::span{units}.first(n),
                std::span{units}.first(n), context()};
    }
    [[nodiscard]] Limits limits(core::Index n, core::Index k = 1, core::Index r = 0) const {
        return {n, k, r, 5, EstimateAdmission{1 << 24}};
    }
    [[nodiscard]] Accuracy accuracy(double tolerance = 1e-12) const {
        return {Norm::Infinity,         tolerance, 0, tolerance, core::ObservationId{20},
                core::ObservationId{21}};
    }
};
} // namespace
TEST_CASE("Dense solve has only its declared owning output allocations",
          "[property][N4-05][N4-06][N4-09]") {
    Setup f;
    const NoChecks<double> evaluator;
    auto a = matrix(1, 1, std::array{3.}, f.plan);
    auto b = matrix(1, 2, std::array{1., 2.}, f.plan);
    auto rhs = PrecisionBuffer<double>::binary64(std::move(b), core::NumericalGeneration{30});
    REQUIRE(rhs);
    auto scales = f.scales(1);
    auto limits = f.limits(1, 2);
    const DenseOperator op{a, context()};
    auto prepared = DensePrepared::prepare(op, scales, limits, f.plan);
    REQUIRE(prepared);
    auto initial_hi = matrix(1, 2, std::array{1., 2.}, f.plan);
    auto initial_lo = matrix(1, 2, std::array{std::ldexp(1., -60), std::ldexp(1., -59)}, f.plan);
    auto retained_rhs = PrecisionBuffer<double>::two_word(
        std::move(initial_hi), std::move(initial_lo), core::NumericalGeneration{30});
    REQUIRE(retained_rhs);
    for (const bool initially_two_word : {false, true})
        for (const bool projection : {false, true})
            for (const bool strict : {false, true}) {
                const Precision<double> precision{
                    projection
                        ? ReturnedPrecision{Binary64Projection{core::NumericalGeneration{31},
                                                               core::NumericalGeneration{32}}}
                        : ReturnedPrecision{RetainedPrecision{core::NumericalGeneration{31}}},
                    {},
                    evaluator};
                test_support::AllocationAudit audit{test_support::strict_kernel_backends,
                                                    test_support::AuditThreadDomain::CallingThread};
                auto solved = solve_linear(op,
                                           Rhs{initially_two_word ? *retained_rhs : *rhs, context(),
                                               core::RhsGeneration{40}},
                                           scales, precision, f.accuracy(strict ? 1e-30 : 1e-12),
                                           limits, *prepared, f.plan);
                const auto observed = audit.finish();
                // Six metadata arrays + two candidate planes, and one additional projected plane.
                // Kernel allocation is detected as an excess, not inferred from an allocation-free
                // public API.
                CHECK(observed.cpp_allocation_attempts == (projection ? 9 : 8));
                CHECK(observed.status == test_support::AllocationAuditStatus::AllocationObserved);
                if (projection && strict)
                    REQUIRE_FALSE(solved);
                else
                    REQUIRE(solved);
            }
}
TEST_CASE("Complex rank-update solve has only its declared owning output allocations",
          "[property][N4-05][N4-06][N4-09]") {
    Setup f;
    const NoChecks<Complex> evaluator;
    auto a = matrix(1, 1, std::array{2.}, f.plan);
    auto u = matrix(1, 1, std::array{Complex{0, 1. / 16}}, f.plan);
    auto v = matrix(1, 1, std::array{Complex{1, 0}}, f.plan);
    auto b = matrix(1, 1, std::array{Complex{31. / 16, 33. / 16}}, f.plan);
    auto rhs = PrecisionBuffer<Complex>::binary64(std::move(b), core::NumericalGeneration{30});
    REQUIRE(rhs);
    const auto scales = f.scales(1);
    const auto limits = f.limits(1, 1, 1);
    const ColumnUpdatedOperator op{a, u, v, context()};
    auto prepared = UpdatedPrepared::prepare(op, scales, limits, f.plan);
    REQUIRE(prepared);
    for (const bool projection : {false, true}) {
        const Precision<Complex> precision{
            projection ? ReturnedPrecision{Binary64Projection{core::NumericalGeneration{31},
                                                              core::NumericalGeneration{32}}}
                       : ReturnedPrecision{RetainedPrecision{core::NumericalGeneration{31}}},
            {},
            evaluator};
        test_support::AllocationAudit audit{test_support::strict_kernel_backends,
                                            test_support::AuditThreadDomain::CallingThread};
        auto solved = solve_linear(op, Rhs{*rhs, context(), core::RhsGeneration{40}}, scales,
                                   precision, f.accuracy(), limits, *prepared, f.plan);
        const auto observed = audit.finish();
        REQUIRE(solved);
        CHECK(observed.cpp_allocation_attempts == (projection ? 9 : 8));
        CHECK(observed.status == test_support::AllocationAuditStatus::AllocationObserved);
    }
}
