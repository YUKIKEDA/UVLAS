#include "../support/allocation_audit.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <span>
#include <utility>
#include <uvlas/numerics/gmres_solve.hpp>

namespace {
using namespace uvlas;
using namespace numerics;
using namespace numerics::linear;
constexpr auto operation = core::Operation::SolveLinear;
[[nodiscard]] Context context() {
    return {{core::ModelFingerprint{}, core::ModelGeneration{1}, core::GeometryGeneration{2},
             core::ChartGeneration{3}, core::SparsityGeneration{4}, core::NumericalGeneration{5},
             core::ScaleGeneration{6}, core::PrecisionGeneration{7}, core::Precision::Binary64,
             core::WakePolicyGeneration{8}, core::TimeIntervalGeneration{9},
             core::BackendGeneration{10}},
            core::TrialId{11},
            core::ScopeId{12},
            core::ScopeId{13},
            .25};
}
// The bounded scalar fixture needs no external backend workspace or allocation.
class ScalarOperator final : public GmresOperator {
public:
    [[nodiscard]] core::Index order() const noexcept override {
        return 1;
    }
    [[nodiscard]] Context context() const noexcept override {
        return ::context();
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        return {};
    }
    [[nodiscard]] core::Result<void> apply_word(std::span<const double> x, std::span<double> value,
                                                std::span<double> tail,
                                                std::span<double> error) const override {
        value[0] = 3 * x[0];
        tail[0] = std::fma(3., x[0], -value[0]);
        error[0] = 0.;
        return {};
    }
    [[nodiscard]] core::Result<void> absolute_action_lower(std::span<const double> weights,
                                                           std::span<double> value) const override {
        const double product = 3 * weights[0];
        const double remainder = std::fma(3., weights[0], -product);
        value[0] = remainder < 0 ? std::nextafter(product, 0.) : product;
        return {};
    }
};
class Identity final : public FixedPreconditioner {
public:
    [[nodiscard]] core::Index order() const noexcept override {
        return 1;
    }
    [[nodiscard]] FixedPreconditionerIdentity identity() const noexcept override {
        return {::context(), ::context().cache};
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        return {};
    }
    [[nodiscard]] core::Result<void> apply(std::span<const double> x,
                                           std::span<double> y) const override {
        y[0] = x[0];
        return {};
    }
};
class NoChecks final : public ConsumerEvaluator<double> {
public:
    [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
        return SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<double>&, const PrecisionBuffer<double>&, const Context&,
           const ConsumerCheck& check, core::CertificateContext expected) const override {
        FAIL_CHECK("No physical observation was requested in this allocation fixture");
        return core::ErrorCertificate::unevaluated(expected, check.budget);
    }
};
[[nodiscard]] DenseReal scalar(double value, core::ResourcePlan& plan) {
    const auto input = std::array{value};
    auto result = DenseReal::copy({1, 1, 1}, input, plan, operation);
    REQUIRE(result);
    return std::move(*result);
}
} // namespace

TEST_CASE("GMRES allocates only declared owning outputs for every retained and projection path",
          "[property][N4-05][N4-06][N4-09]") {
    core::ResourcePlan plan{1 << 20};
    const auto binary =
        PrecisionBuffer<double>::binary64(scalar(1., plan), core::NumericalGeneration{30});
    const auto words = PrecisionBuffer<double>::two_word(
        scalar(1., plan), scalar(std::ldexp(1., -60), plan), core::NumericalGeneration{30});
    REQUIRE(binary);
    REQUIRE(words);
    const std::array scales{1.};
    const std::array<UnitDimensions, 1> units{};
    const Scales physical{scales, scales, units, units, context()};
    const ScalarOperator op;
    const Identity preconditioner;
    const NoChecks evaluator;
    const GmresLimits limits{1, 1, EstimateAdmission{0}};
    auto prepared = GmresPrepared::prepare(op, preconditioner, physical, limits, plan);
    REQUIRE(prepared);
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
                const double tolerance = strict ? 1.e-30 : 1.e-12;
                const Accuracy accuracy{
                    Norm::Infinity,         tolerance, 0., tolerance, core::ObservationId{20},
                    core::ObservationId{21}};
                const Rhs<double> rhs{initially_two_word ? *words : *binary, context(),
                                      core::RhsGeneration{40}};
                test_support::AllocationAudit audit{test_support::strict_kernel_backends,
                                                    test_support::AuditThreadDomain::CallingThread};
                const auto result = solve_linear(op, preconditioner, rhs, physical, precision,
                                                 accuracy, limits, *prepared, plan);
                const auto observed = audit.finish();
                CAPTURE(initially_two_word, projection, strict);
                if (!result) {
                    CAPTURE(static_cast<unsigned>(result.error().failure().code),
                            static_cast<unsigned>(result.error().diagnostics().stop),
                            result.error().diagnostics().operator_actions);
                }
                // Source inspection identifies six metadata arrays and two candidate planes.
                // Projection adds one plane. The public-call count catches excess allocation;
                // it is not a separate direct measurement of a zero-allocation inner interval.
                CHECK(observed.cpp_allocation_attempts == (projection ? 9 : 8));
                CHECK(observed.status == test_support::AllocationAuditStatus::AllocationObserved);
                if (projection && strict) {
                    REQUIRE_FALSE(result);
                    CHECK(result.error().failure().code == core::FailureCode::PrecisionUnresolved);
                } else {
                    REQUIRE(result);
                    CHECK(result->diagnostics().operator_actions > 0);
                    const auto retained = initially_two_word || strict ? core::Precision::TwoWord
                                                                       : core::Precision::Binary64;
                    const auto returned = projection ? core::Precision::Binary64 : retained;
                    CHECK(result->solution().values().precision() == returned);
                    CHECK(result->solution()
                              .certificate()
                              .retained_residual()
                              .context()
                              .representation == retained);
                    CHECK(result->solution().certificate().residual().context().representation ==
                          returned);
                }
            }
}
