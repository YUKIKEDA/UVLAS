#include "backend_memory_detail.hpp"
#include "linear_solve_detail.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <variant>

namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;
namespace factory = uvlas::numerics::linear::detail;

Context system_context() {
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
Accuracy request() {
    return {Norm::Infinity, 1.e-12, 0, 1.e-12, core::ObservationId{20}, core::ObservationId{21}};
}
BackendAdmission admission() {
    auto value = admit_backend_memory(factory::BackendPlanFactory::make(16), EstimateAdmission{32});
    REQUIRE(value);
    return *value;
}

struct Fixture {
    core::ResourcePlan plan{65536};
    std::array<double, 1> unknown{2};
    std::array<double, 1> residual{3};
    std::array<UnitDimensions, 1> unknown_units{UnitDimensions{1, 0, 0, 0, 0, 0, 0}};
    std::array<UnitDimensions, 1> residual_units{UnitDimensions{1, 1, -2, 0, 0, 0, 0}};
    Scales scales() const {
        return {unknown, residual, unknown_units, residual_units, system_context()};
    }
    factory::MetadataBuffers metadata(std::span<const ConsumerCheck> checks = {}) {
        auto value = factory::SolutionFactory::prepare_metadata(1, 1, scales(), checks, plan);
        REQUIRE(value);
        value->retained_columns.values()[0] = ColumnMetrics{0.25, 1.e-14, 0};
        value->returned_columns.values()[0] = ColumnMetrics{0.5, 2.e-14, 0};
        return std::move(*value);
    }
    PrecisionBuffer<double> buffer(core::NumericalGeneration generation, bool two_word = false) {
        const std::array high{4.0};
        auto hi = DenseReal::copy({1, 1, 1}, high, plan, core::Operation::SolveLinear);
        REQUIRE(hi);
        if (two_word) {
            const std::array low{std::ldexp(1.0, -54)};
            auto lo = DenseReal::copy({1, 1, 1}, low, plan, core::Operation::SolveLinear);
            REQUIRE(lo);
            auto value =
                PrecisionBuffer<double>::two_word(std::move(*hi), std::move(*lo), generation);
            REQUIRE(value);
            return std::move(*value);
        }
        auto value = PrecisionBuffer<double>::binary64(std::move(*hi), generation);
        REQUIRE(value);
        return std::move(*value);
    }
    factory::Finalization finalization(bool two_word = false) const {
        return {system_context(),
                core::RhsGeneration{30},
                core::NumericalGeneration{31},
                request(),
                admission(),
                admission(),
                2,
                two_word ? RetentionPath::InitiallyTwoWord : RetentionPath::Binary64,
                ReciprocalOneNormEstimate{0.1}};
    }
};

enum class Fault {
    None,
    Observation,
    Scope,
    Generation,
    Representation,
    Budget,
    Exceeded,
    Unevaluated,
    Throw
};
class Evaluator final : public ConsumerEvaluator<double> {
public:
    Fault fault = Fault::None;
    SupportedPrecision support = SupportedPrecision::Binary64AndTwoWord;
    unsigned fault_call = 1;
    mutable unsigned calls = 0;
    mutable bool saw_retained_low = false;
    mutable bool saw_projection = false;
    SupportedPrecision supported_precision() const noexcept override {
        return support;
    }
    core::Result<core::ErrorCertificate> assess(const PrecisionBuffer<double>& retained,
                                                const PrecisionBuffer<double>& candidate,
                                                const Context& context, const ConsumerCheck& check,
                                                core::CertificateContext expected) const override {
        ++calls;
        const auto active_fault = calls == fault_call ? fault : Fault::None;
        CHECK(context == system_context());
        CHECK(expected.numerical_generation == candidate.generation());
        if (const auto* words = std::get_if<TwoWord<double>>(&retained.storage())) {
            saw_retained_low = words->lo.values()[0] == std::ldexp(1.0, -54);
        }
        saw_projection = retained.generation() != candidate.generation();
        if (active_fault == Fault::Throw) {
            throw std::runtime_error("consumer failed");
        }
        if (active_fault == Fault::Observation) {
            expected.observation = core::ObservationId{999};
        }
        if (active_fault == Fault::Scope) {
            expected.scope = core::ScopeId{999};
        }
        if (active_fault == Fault::Generation) {
            expected.numerical_generation = core::NumericalGeneration{999};
        }
        if (active_fault == Fault::Representation) {
            expected.representation = core::Precision::TwoWord;
        }
        const double budget = active_fault == Fault::Budget ? check.budget * 2 : check.budget;
        if (active_fault == Fault::Unevaluated) {
            return core::ErrorCertificate::unevaluated(expected, budget);
        }
        return core::ErrorCertificate::assess(
            expected, budget, active_fault == Fault::Exceeded ? budget * 2 : budget / 2);
    }
};
constexpr std::array checks{
    ConsumerCheck{CheckKind::Observation, core::ObservationId{40}, core::ScopeId{41}, 1.e-6},
    ConsumerCheck{CheckKind::Constraint, core::ObservationId{42}, core::ScopeId{43}, 2.e-6},
    ConsumerCheck{CheckKind::Reaction, core::ObservationId{44}, core::ScopeId{45}, 3.e-6}};
} // namespace

TEST_CASE("Result owns physical metadata and every necessary consumer proof",
          "[unit][verification][N4-05][N4-06]") {
    Fixture fixture;
    Evaluator evaluator;
    auto metadata = fixture.metadata(checks);
    fixture.unknown[0] = 100;
    fixture.residual[0] = 200;
    auto value = factory::SolutionFactory::finalize(
        fixture.buffer(core::NumericalGeneration{50}), std::optional<PrecisionBuffer<double>>{},
        Precision<double>{RetainedPrecision{core::NumericalGeneration{50}}, checks, evaluator},
        std::move(metadata), fixture.finalization());
    REQUIRE(value);
    CHECK(evaluator.calls == 3);
    CHECK_FALSE(evaluator.saw_projection);
    const auto& certificate = value->certificate();
    CHECK(certificate.context() == system_context());
    CHECK(certificate.rhs_generation() == core::RhsGeneration{30});
    CHECK(certificate.rhs_buffer_generation() == core::NumericalGeneration{31});
    CHECK(certificate.unknown_scales()[0] == 2);
    CHECK(certificate.residual_scales()[0] == 3);
    CHECK(certificate.unknown_units()[0][0] == 1);
    CHECK(certificate.residual_units()[0][2] == -2);
    CHECK(certificate.consumer_count() == 3);
    CHECK(certificate.columns()[0].normalized_original_residual_upper == 0.25);
    CHECK(certificate.refinements() == 2);
    CHECK(certificate.preparation_admission().admitted_payload_bytes() == 16);
    for (core::Index i = 0; i < 3; ++i) {
        const auto proof = certificate.consumer_proof(i);
        REQUIRE(proof);
        CHECK(proof->request.kind == checks[static_cast<std::size_t>(i)].kind);
        CHECK(proof->assessment.context().numerical_generation == value->values().generation());
        CHECK(std::holds_alternative<core::Achieved>(proof->assessment.assessment()));
    }
    CHECK_FALSE(certificate.consumer_proof(-1));
    CHECK_FALSE(certificate.consumer_proof(3));
}

TEST_CASE("Consumer proofs must match exact identity budget representation and achievement",
          "[negative][N4-05][N4-06]") {
    for (const auto fault :
         {Fault::Observation, Fault::Scope, Fault::Generation, Fault::Representation, Fault::Budget,
          Fault::Exceeded, Fault::Unevaluated, Fault::Throw}) {
        Fixture fixture;
        Evaluator evaluator;
        evaluator.fault = fault;
        const auto result = factory::SolutionFactory::finalize(
            fixture.buffer(core::NumericalGeneration{50}), std::optional<PrecisionBuffer<double>>{},
            Precision<double>{RetainedPrecision{core::NumericalGeneration{50}}, checks, evaluator},
            fixture.metadata(checks), fixture.finalization());
        REQUIRE_FALSE(result);
        CHECK(evaluator.calls == 1);
        CHECK(fixture.unknown[0] == 2);
        CHECK(fixture.residual[0] == 3);
    }
}

TEST_CASE("Projection certifies the returned candidate while retaining source diagnostics",
          "[verification][negative][N4-06]") {
    Fixture fixture;
    Evaluator evaluator;
    auto result = factory::SolutionFactory::finalize(
        fixture.buffer(core::NumericalGeneration{50}, true),
        std::optional{fixture.buffer(core::NumericalGeneration{51})},
        Precision<double>{
            Binary64Projection{core::NumericalGeneration{50}, core::NumericalGeneration{51}},
            checks, evaluator},
        fixture.metadata(checks), fixture.finalization(true));
    REQUIRE(result);
    CHECK(evaluator.saw_retained_low);
    CHECK(evaluator.saw_projection);
    CHECK(result->values().precision() == core::Precision::Binary64);
    const auto& certificate = result->certificate();
    CHECK(certificate.retained_residual().context().numerical_generation ==
          core::NumericalGeneration{50});
    CHECK(certificate.retained_residual().context().representation == core::Precision::TwoWord);
    CHECK(certificate.residual().context().numerical_generation == core::NumericalGeneration{51});
    CHECK(certificate.residual().context().representation == core::Precision::Binary64);
    CHECK(certificate.retained_columns()[0].normalized_original_residual_upper == 0.25);
    CHECK(certificate.columns()[0].normalized_original_residual_upper == 0.5);
}

TEST_CASE("No algebraic failure or uncomputed column is published as an owning solution",
          "[negative][N4-05][N4-06]") {
    for (int fault = 0; fault < 5; ++fault) {
        Fixture fixture;
        Evaluator evaluator;
        auto metadata = fixture.metadata(checks);
        if (fault == 0) {
            metadata.retained_columns.values()[0].normalized_original_residual_upper = 2;
        }
        if (fault == 1) {
            metadata.retained_columns.values()[0].componentwise_backward_error_upper = 1;
        }
        if (fault == 2) {
            metadata.returned_columns.values()[0].normalized_original_residual_upper = 2;
        }
        if (fault == 3) {
            metadata.returned_columns.values()[0].componentwise_backward_error_upper = 1;
        }
        if (fault == 4) {
            auto uncomputed = factory::SolutionFactory::prepare_metadata(1, 1, fixture.scales(),
                                                                         checks, fixture.plan);
            REQUIRE(uncomputed);
            metadata = std::move(*uncomputed);
        }
        const auto result = factory::SolutionFactory::finalize(
            fixture.buffer(core::NumericalGeneration{50}, true),
            std::optional{fixture.buffer(core::NumericalGeneration{51})},
            Precision<double>{
                Binary64Projection{core::NumericalGeneration{50}, core::NumericalGeneration{51}},
                checks, evaluator},
            std::move(metadata), fixture.finalization(true));
        CHECK_FALSE(result);
        CHECK(evaluator.calls == 0);
    }
}

TEST_CASE("Result preparation validates before allocation and cannot exceed its payload budget",
          "[negative][N4-09]") {
    Fixture fixture;
    core::ResourcePlan none{0};
    CHECK_FALSE(factory::SolutionFactory::prepare_metadata(1, 1, fixture.scales(), checks, none));
    CHECK(none.allocated_bytes() == 0);
    fixture.residual[0] = 0;
    auto invalid = factory::SolutionFactory::prepare_metadata(1, 1, fixture.scales(), checks, none);
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == core::FailureCode::InvalidInput);
    CHECK(none.allocated_bytes() == 0);
}

TEST_CASE("Returned representation generation and consumer precision are enforced",
          "[negative][N4-06]") {
    for (int fault = 0; fault < 4; ++fault) {
        Fixture fixture;
        Evaluator evaluator;
        if (fault == 3) {
            evaluator.support = SupportedPrecision::Binary64Only;
        }
        const auto retained_generation =
            fault == 0 ? core::NumericalGeneration{99} : core::NumericalGeneration{50};
        const auto target_generation =
            fault == 1 ? core::NumericalGeneration{50} : core::NumericalGeneration{51};
        auto finalized = fixture.finalization(true);
        if (fault == 2) {
            finalized.refinements = 6;
        }
        const auto result = factory::SolutionFactory::finalize(
            fixture.buffer(core::NumericalGeneration{50}, true),
            std::optional{fixture.buffer(target_generation)},
            Precision<double>{Binary64Projection{retained_generation, target_generation}, checks,
                              evaluator},
            fixture.metadata(checks), finalized);
        CHECK_FALSE(result);
        CHECK(evaluator.calls == 0);
    }
}

TEST_CASE("A failed final reaction check cannot be replaced by earlier observation success",
          "[negative][N4-05][N4-06]") {
    Fixture fixture;
    Evaluator evaluator;
    evaluator.fault = Fault::Exceeded;
    evaluator.fault_call = 3;
    const auto result = factory::SolutionFactory::finalize(
        fixture.buffer(core::NumericalGeneration{50}), std::optional<PrecisionBuffer<double>>{},
        Precision<double>{RetainedPrecision{core::NumericalGeneration{50}}, checks, evaluator},
        fixture.metadata(checks), fixture.finalization());
    CHECK_FALSE(result);
    CHECK(evaluator.calls == 3);
}
