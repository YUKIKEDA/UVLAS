#include "compensated.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <uvlas/numerics/gmres_solve.hpp>
namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;
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

class Blocks final : public GmresOperator {
public:
    std::array<double, 16> a{4, 1, 2, 0, 0, 3, 0, 1, 1, 0, 5, -1, 0, 2, 0, 4};
    [[nodiscard]] core::Index order() const noexcept override {
        return 4;
    }
    [[nodiscard]] Context context() const noexcept override {
        return ::context();
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        return {};
    }
    [[nodiscard]] core::Result<void> apply_word(std::span<const double> x, std::span<double> out,
                                                std::span<double> tail,
                                                std::span<double> error) const override {
        for (std::size_t i = 0; i < 4; ++i) {
            uvlas::numerics::detail::CompensatedAccumulator sum;
            for (std::size_t j = 0; j < 4; ++j) {
                auto r = sum.product(a[i * 4 + j], x[j], 1);
                if (!r)
                    return r;
            }
            auto value = sum.result();
            if (!value)
                return std::unexpected(value.error());
            out[i] = value->value;
            error[i] = value->error;
            tail[i] = 0;
        }
        return {};
    }
    [[nodiscard]] core::Result<void> absolute_action_lower(std::span<const double> w,
                                                           std::span<double> out) const override {
        // Directed-down products/sums are independent of the Krylov search.
        for (std::size_t i = 0; i < 4; ++i) {
            double sum = 0;
            for (std::size_t j = 0; j < 4; ++j) {
                const double p = std::abs(a[i * 4 + j]) * w[j];
                sum = std::nextafter(sum + std::nextafter(p, 0.0), 0.0);
            }
            out[i] = sum;
        }
        return {};
    }
};
class BlockP final : public FixedPreconditioner {
public:
    [[nodiscard]] core::Index order() const noexcept override {
        return 4;
    }
    [[nodiscard]] FixedPreconditionerIdentity identity() const noexcept override {
        return {::context(), ::context().cache};
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        return {};
    }
    [[nodiscard]] core::Result<void> apply(std::span<const double> r,
                                           std::span<double> x) const override {
        x[1] = r[1] / 3;
        x[0] = (r[0] - x[1]) / 4;
        x[3] = r[3] / 4;
        x[2] = (r[2] + x[3]) / 5;
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
           const ConsumerCheck&, core::CertificateContext) const override {
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::SolveLinear, {}});
    }
};
} // namespace
TEST_CASE("Restarted GMRES retains every explicit coupling block", "[N4-05][N4-06]") {
    Blocks const a;
    BlockP const p;
    NoChecks const consumer;
    core::ResourcePlan inputs{1 << 20}, workspace{1 << 20}, output{1 << 20};
    const std::array<double, 4> scales_data{1, 2, 4, 8};
    const std::array<UnitDimensions, 4> units{};
    const Scales scales{scales_data, scales_data, units, units, context()};
    const GmresLimits limits{4, 1, EstimateAdmission{1 << 20}, 2, 512};
    auto prepared = GmresPrepared::prepare(a, p, scales, limits, workspace);
    REQUIRE(prepared);
    const std::array<double, 4> b{8, -7, 17, -8};
    auto dense = DenseReal::copy({4, 1, 4}, b, inputs, core::Operation::SolveLinear);
    REQUIRE(dense);
    auto rhs = PrecisionBuffer<double>::binary64(std::move(*dense), core::NumericalGeneration{1});
    REQUIRE(rhs);
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{2}}, {}, consumer};
    const Accuracy accuracy{Norm::Infinity,        1e-13, 0, 1e-13, core::ObservationId{1},
                            core::ObservationId{2}};
    auto solved = solve_linear(a, p, Rhs<double>{*rhs, context(), core::RhsGeneration{1}}, scales,
                               precision, accuracy, limits, *prepared, output);
    REQUIRE(solved);
    const auto* binary = std::get_if<Binary64<double>>(&solved->solution().values().storage());
    REQUIRE(binary);
    const std::array<double, 4> expected{1, -2, 3, -1};
    for (std::size_t i = 0; i < 4; ++i)
        CHECK(std::abs(binary->values.values()[i] - expected[i]) < 1e-12);
    CHECK(solved->diagnostics().restarts > 0);
    CHECK(solved->diagnostics().operator_actions <= 512);
    CHECK(solved->solution().certificate().refinements() == 0);
}
TEST_CASE("GMRES default resources are independent of numerical accuracy", "[N4-05][N4-06]") {
    const GmresLimits defaults{4, 1, EstimateAdmission{1024}};
    CHECK(defaults.restart == 64);
    CHECK(defaults.maximum_operator_actions == 512);
}
TEST_CASE("Inconsistent singular breakdown cannot publish a Hessenberg solution",
          "[N4-05][N4-06]") {
    Blocks a;
    a.a.fill(0);
    BlockP const p;
    NoChecks const consumer;
    core::ResourcePlan inputs{1 << 20}, workspace{1 << 20}, output{1 << 20};
    const std::array<double, 4> scale{1, 1, 1, 1};
    const std::array<UnitDimensions, 4> units{};
    const Scales scales{scale, scale, units, units, context()};
    const GmresLimits limits{4, 1, EstimateAdmission{1024}, 4, 512};
    auto prepared = GmresPrepared::prepare(a, p, scales, limits, workspace);
    REQUIRE(prepared);
    const std::array<double, 4> b{1, 0, 0, 0};
    auto matrix = DenseReal::copy({4, 1, 4}, b, inputs, core::Operation::SolveLinear);
    REQUIRE(matrix);
    auto rhs = PrecisionBuffer<double>::binary64(std::move(*matrix), core::NumericalGeneration{1});
    REQUIRE(rhs);
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{2}}, {}, consumer};
    const Accuracy accuracy{
        Norm::L2, 1e-12, 0, 1e-12, core::ObservationId{1}, core::ObservationId{2}};
    auto result = solve_linear(a, p, Rhs<double>{*rhs, context(), core::RhsGeneration{1}}, scales,
                               precision, accuracy, limits, *prepared, output);
    REQUIRE_FALSE(result);
    CHECK(result.error().failure().code == core::FailureCode::LinearUnconverged);
    CHECK(result.error().diagnostics().stop == GmresStopReason::KrylovBreakdown);
    CHECK(result.error().diagnostics().operator_actions == 1);
    CHECK(result.error().diagnostics().restarts == 0);
}
TEST_CASE("Zero right hand side qualifies without numerical callbacks", "[N4-05][N4-06]") {
    Blocks const a;
    BlockP const p;
    NoChecks const consumer;
    core::ResourcePlan inputs{1 << 20}, workspace{1 << 20}, output{1 << 20};
    const std::array<double, 4> scale{1, 1, 1, 1};
    const std::array<UnitDimensions, 4> units{};
    const Scales scales{scale, scale, units, units, context()};
    const GmresLimits limits{4, 1, EstimateAdmission{1024}, 4, 0};
    auto prepared = GmresPrepared::prepare(a, p, scales, limits, workspace);
    REQUIRE(prepared);
    const std::array<double, 4> b{};
    auto matrix = DenseReal::copy({4, 1, 4}, b, inputs, core::Operation::SolveLinear);
    REQUIRE(matrix);
    auto rhs = PrecisionBuffer<double>::binary64(std::move(*matrix), core::NumericalGeneration{1});
    REQUIRE(rhs);
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{2}}, {}, consumer};
    const Accuracy accuracy{Norm::L2, 0, 0, 0, core::ObservationId{1}, core::ObservationId{2}};
    auto result = solve_linear(a, p, Rhs<double>{*rhs, context(), core::RhsGeneration{1}}, scales,
                               precision, accuracy, limits, *prepared, output);
    REQUIRE(result);
    CHECK(result->diagnostics().operator_actions == 0);
    CHECK(result->diagnostics().preconditioner_actions == 0);
}
