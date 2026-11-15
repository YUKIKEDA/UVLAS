#include "linear_candidate_detail.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <span>
#include <variant>

namespace {
using namespace uvlas;
using namespace numerics::linear;
using namespace numerics::linear::detail;
[[nodiscard]] Accuracy accuracy(double absolute = 0) {
    return {Norm::Infinity, absolute, 0., .1, core::ObservationId{1}, core::ObservationId{2}};
}
} // namespace

TEST_CASE("Refinement distinguishes usable residuals from finite certificate metrics",
          "[N4-05][N4-06]") {
    const std::array<double, 1> nonzero{1.}, zero{}, unit{1.};
    const auto unmet = assess_candidate_column_metrics<double>(1, nonzero, zero, unit, nonzero, {},
                                                               unit, accuracy());
    REQUIRE(unmet);
    CHECK(std::holds_alternative<UnrepresentableMetrics>(*unmet));
    CHECK_FALSE(candidate_meets_accuracy(*unmet, accuracy()));
    const auto exact =
        assess_candidate_column_metrics<double>(1, zero, zero, unit, nonzero, {}, unit, accuracy());
    REQUIRE(exact);
    CHECK(candidate_meets_accuracy(*exact, accuracy()));
    const auto finite_unmet = assess_candidate_column_metrics<double>(
        1, nonzero, zero, unit, nonzero, {}, unit, accuracy(.25));
    REQUIRE(finite_unmet);
    CHECK(std::holds_alternative<ColumnMetrics>(*finite_unmet));
    CHECK_FALSE(candidate_meets_accuracy(*finite_unmet, accuracy(.25)));
    const auto empty =
        assess_candidate_column_metrics<double>(0, {}, {}, {}, {}, {}, {}, accuracy());
    REQUIRE(empty);
    CHECK(candidate_meets_accuracy(*empty, accuracy()));
}

TEST_CASE("An early unrepresented ratio cannot hide a later invalid raw row",
          "[negative][N4-05][N4-06][N4-08]") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (unsigned fault = 0; fault < 7; ++fault) {
        std::array<double, 2> residual{1., 0.}, error{}, denominator{1., 1.}, rhs{1., 1.}, low{},
            scales{1., 1.};
        auto expected = core::FailureCode::Nonfinite;
        switch (fault) {
        case 0:
            residual[1] = nan;
            break;
        case 1:
            error[1] = std::numeric_limits<double>::infinity();
            expected = core::FailureCode::PrecisionUnresolved;
            break;
        case 2:
            error[1] = -1.;
            expected = core::FailureCode::InvalidInput;
            break;
        case 3:
            denominator[1] = nan;
            break;
        case 4:
            rhs[1] = nan;
            break;
        case 5:
            low[1] = nan;
            break;
        default:
            scales[1] = 0.;
            expected = core::FailureCode::InvalidInput;
            break;
        }
        const auto rejected = assess_candidate_column_metrics<double>(
            2, residual, error, denominator, rhs, low, scales, accuracy());
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == expected);
    }
    using Complex = std::complex<double>;
    const std::array<Complex, 2> residual{Complex{1., 0.}, Complex{0., nan}}, rhs{1., 1.};
    const std::array<double, 2> zero{}, unit{1., 1.};
    const auto imaginary = assess_candidate_column_metrics<Complex>(2, residual, zero, unit, rhs,
                                                                    {}, unit, accuracy());
    REQUIRE_FALSE(imaginary);
    CHECK(imaginary.error().code == core::FailureCode::Nonfinite);
}

TEST_CASE("Default refinement history is unqualified and promotion preserves finite half progress",
          "[property][N4-05][N4-06]") {
    const RefinementProgress prepared;
    CHECK_FALSE(candidate_meets_accuracy(prepared.current, accuracy()));
    CHECK_FALSE(candidate_meets_accuracy(prepared.previous, accuracy()));
    const CandidateMetrics before{ColumnMetrics{8., .5, 0.}};
    const CandidateMetrics halved{ColumnMetrics{4., .25, 0.}};
    const CandidateMetrics stalled{ColumnMetrics{5., .3, 0.}};
    const CandidateMetrics achieved{ColumnMetrics{.5, .05, 0.}};
    CHECK_FALSE(needs_precision_promotion(halved, before, accuracy()));
    CHECK(needs_precision_promotion(stalled, before, accuracy()));
    CHECK_FALSE(needs_precision_promotion(achieved, before, accuracy()));
    CHECK(needs_precision_promotion(prepared.current, before, accuracy()));
    CHECK(needs_precision_promotion(stalled, prepared.previous, accuracy()));
    CHECK_FALSE(needs_precision_promotion(achieved, prepared.previous, accuracy()));
    const CandidateMetrics residual_met{ColumnMetrics{.5, .3, 0.}};
    CHECK(needs_precision_promotion(residual_met, before, accuracy()));
}
