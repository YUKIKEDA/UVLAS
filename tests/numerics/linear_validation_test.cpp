#include "linear_solve_detail.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>

namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;
namespace checks = uvlas::numerics::linear::detail;
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
class Consumer final : public ConsumerEvaluator<double> {
public:
    SupportedPrecision support = SupportedPrecision::Binary64AndTwoWord;
    [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
        return support;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<double>&, const PrecisionBuffer<double>&, const Context&,
           const ConsumerCheck& request, core::CertificateContext expected) const override {
        FAIL_CHECK("preflight must not execute a numerical consumer");
        return core::ErrorCertificate::unevaluated(expected, request.budget);
    }
};
struct Inputs {
    core::ResourcePlan plan{4096};
    std::array<double, 1> unknown{1};
    std::array<double, 1> residual{2};
    std::array<UnitDimensions, 1> units{UnitDimensions{}};
    Consumer consumer;
    [[nodiscard]] Scales scales() const {
        return {unknown, residual, units, units, context()};
    }
    [[nodiscard]] Accuracy accuracy() const {
        return {Norm::Infinity,         1.e-12, 0, 1.e-12, core::ObservationId{20},
                core::ObservationId{21}};
    }
    [[nodiscard]] Limits limits() const {
        return {1, 1, 1, 5, EstimateAdmission{4096}};
    }
    PrecisionBuffer<double> buffer(bool two_word = false) {
        const std::array high{1.0};
        auto hi = DenseReal::copy({1, 1, 1}, high, plan, core::Operation::SolveLinear);
        REQUIRE(hi);
        if (two_word) {
            const std::array low{std::ldexp(1.0, -54)};
            auto lo = DenseReal::copy({1, 1, 1}, low, plan, core::Operation::SolveLinear);
            REQUIRE(lo);
            auto value = PrecisionBuffer<double>::two_word(std::move(*hi), std::move(*lo),
                                                           core::NumericalGeneration{30});
            REQUIRE(value);
            return std::move(*value);
        }
        auto value =
            PrecisionBuffer<double>::binary64(std::move(*hi), core::NumericalGeneration{30});
        REQUIRE(value);
        return std::move(*value);
    }
};
constexpr std::array requirement{
    ConsumerCheck{CheckKind::Observation, core::ObservationId{40}, core::ScopeId{41}, 1.e-6}};
} // namespace

TEST_CASE("Changed RHS identity never bypasses any matrix cache field",
          "[negative][N4-03][N4-05]") {
    using Mutate = void (*)(Context&);
    const std::array<Mutate, 16> mutations{
        [](Context& c) { c.cache.fingerprint.bytes[0] = std::byte{1}; },
        [](Context& c) { c.cache.model = core::ModelGeneration{99}; },
        [](Context& c) { c.cache.geometry = core::GeometryGeneration{99}; },
        [](Context& c) { c.cache.chart = core::ChartGeneration{99}; },
        [](Context& c) { c.cache.sparsity = core::SparsityGeneration{99}; },
        [](Context& c) { c.cache.matrix = core::NumericalGeneration{99}; },
        [](Context& c) { c.cache.scales = core::ScaleGeneration{99}; },
        [](Context& c) { c.cache.precision = core::PrecisionGeneration{99}; },
        [](Context& c) { c.cache.representation = core::Precision::TwoWord; },
        [](Context& c) { c.cache.wake = core::WakePolicyGeneration{99}; },
        [](Context& c) { c.cache.interval = core::TimeIntervalGeneration{99}; },
        [](Context& c) { c.cache.backend = core::BackendGeneration{99}; },
        [](Context& c) { c.trial = core::TrialId{99}; },
        [](Context& c) { c.unknown_layout = core::ScopeId{99}; },
        [](Context& c) { c.residual_layout = core::ScopeId{99}; },
        [](Context& c) { c.time = 0.5; }};
    Inputs input;
    const auto values = input.buffer();
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{50}}, {}, input.consumer};
    REQUIRE(checks::check_solve_inputs(
        1, 0, context(), Rhs<double>{values, context(), core::RhsGeneration{999}}, input.scales(),
        precision, input.accuracy(), input.limits()));
    for (const auto mutate : mutations) {
        auto changed = context();
        mutate(changed);
        const auto result = checks::check_solve_inputs(
            1, 0, context(), Rhs<double>{values, changed, core::RhsGeneration{999}}, input.scales(),
            precision, input.accuracy(), input.limits());
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::StaleGeneration);
    }
}

TEST_CASE("Physical scale metadata is finite positive complete and bound to its context",
          "[negative][N4-05]") {
    Inputs input;
    const auto values = input.buffer();
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{50}}, {}, input.consumer};
    for (int fault = 0; fault < 7; ++fault) {
        input.unknown[0] = 1;
        input.residual[0] = 2;
        auto scales = input.scales();
        if (fault == 0) {
            input.unknown[0] = 0;
        }
        if (fault == 1) {
            input.residual[0] = -1;
        }
        if (fault == 2) {
            input.residual[0] = std::numeric_limits<double>::infinity();
        }
        if (fault == 3) {
            scales.unknown_units = {};
        }
        if (fault == 4) {
            scales.residual_units = {};
        }
        if (fault == 5) {
            scales.context.cache.scales = core::ScaleGeneration{99};
        }
        if (fault == 6) {
            scales.context.time = std::numeric_limits<double>::quiet_NaN();
        }
        const auto result = checks::check_solve_inputs(
            1, 0, context(), Rhs<double>{values, context(), core::RhsGeneration{1}}, scales,
            precision, input.accuracy(), input.limits());
        CHECK_FALSE(result);
    }
}

TEST_CASE("Known unsupported consumers fail preflight without evaluating a candidate",
          "[negative][N4-06]") {
    Inputs input;
    const auto words = input.buffer(true);
    const auto plain = input.buffer();
    input.consumer.support = SupportedPrecision::Binary64Only;
    const Precision<double> retained{RetainedPrecision{core::NumericalGeneration{50}}, requirement,
                                     input.consumer};
    const auto unsupported = checks::check_solve_inputs(
        1, 0, context(), Rhs<double>{words, context(), core::RhsGeneration{1}}, input.scales(),
        retained, input.accuracy(), input.limits());
    REQUIRE_FALSE(unsupported);
    CHECK(unsupported.error().code == core::FailureCode::PrecisionUnsupported);
    input.consumer.support = static_cast<SupportedPrecision>(255);
    const auto invalid = checks::check_solve_inputs(
        1, 0, context(), Rhs<double>{plain, context(), core::RhsGeneration{1}}, input.scales(),
        retained, input.accuracy(), input.limits());
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == core::FailureCode::PrecisionUnsupported);
}

TEST_CASE("Capacities use inclusive boundaries and preserve the combined refinement limit",
          "[negative][N4-05][N4-09]") {
    const Inputs input;
    REQUIRE(checks::check_capacity_inputs(1, 1, 1, input.limits()));
    for (int fault = 0; fault < 7; ++fault) {
        auto limits = input.limits();
        core::Index rows = 1, columns = 1, rank = 1;
        if (fault == 0) {
            rows = 2;
        }
        if (fault == 1) {
            columns = 2;
        }
        if (fault == 2) {
            rank = 2;
        }
        if (fault == 3) {
            limits.maximum_order = -1;
        }
        if (fault == 4) {
            limits.maximum_rhs_columns = -1;
        }
        if (fault == 5) {
            limits.maximum_update_rank = -1;
        }
        if (fault == 6) {
            limits.maximum_refinements = 6;
        }
        const auto result = checks::check_capacity_inputs(rows, columns, rank, limits);
        REQUIRE_FALSE(result);
        CHECK(result.error().code ==
              (fault < 3 ? core::FailureCode::ResourceExceeded : core::FailureCode::InvalidInput));
    }
    CHECK(checks::check_capacity_inputs(0, 0, 0, Limits{0, 0, 0, 0, EstimateAdmission{0}}));
}

TEST_CASE("A projected TwoWord source still needs a consumer that can assess its retained words",
          "[negative][N4-06]") {
    Inputs input;
    const auto words = input.buffer(true);
    input.consumer.support = SupportedPrecision::Binary64Only;
    const Precision<double> projection{
        Binary64Projection{core::NumericalGeneration{50}, core::NumericalGeneration{51}},
        requirement, input.consumer};
    const auto result = checks::check_solve_inputs(
        1, 0, context(), Rhs<double>{words, context(), core::RhsGeneration{1}}, input.scales(),
        projection, input.accuracy(), input.limits());
    REQUIRE_FALSE(result);
    CHECK(result.error().code == core::FailureCode::PrecisionUnsupported);
}

TEST_CASE("Failed policy emplacement is preserved or rejected before numerical work",
          "[negative][N4-06][N4-08]") {
    struct ThrowingGeneration {
        operator core::NumericalGeneration() const {
            throw std::runtime_error("policy construction failed");
        }
    };
    Inputs input;
    const auto rhs = input.buffer();
    ReturnedPrecision returned{RetainedPrecision{core::NumericalGeneration{50}}};
    bool threw = false;
    try {
        returned.emplace<Binary64Projection>(ThrowingGeneration{}, core::NumericalGeneration{51});
    } catch (const std::runtime_error&) {
        threw = true;
    }
    REQUIRE(threw);
    const Precision<double> precision{returned, {}, input.consumer};
    const auto result = checks::check_solve_inputs(
        1, 0, context(), Rhs<double>{rhs, context(), core::RhsGeneration{1}}, input.scales(),
        precision, input.accuracy(), input.limits());
    if (returned.valueless_by_exception()) {
        REQUIRE_FALSE(result);
        CHECK(result.error().code == core::FailureCode::PrecisionUnsupported);
    } else {
        // The standard library may retain the old alternative as a stronger guarantee.
        REQUIRE(std::holds_alternative<RetainedPrecision>(returned));
        CHECK(result);
    }
}
