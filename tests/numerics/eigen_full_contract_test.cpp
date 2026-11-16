#include <catch2/catch_test_macros.hpp>

#include <array>
#include <bit>
#include <complex>
#include <concepts>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <uvlas/numerics/eigen_full.hpp>

namespace {
using namespace uvlas;
namespace eig = numerics::eigen;
namespace lin = numerics::linear;
using Complex = std::complex<double>;
constexpr auto operation = core::Operation::ExtractEigen;
constexpr eig::StateDiagnosticContext reporting{core::AcceptedStateId{31},
                                                core::AcceptedStateId{32}};

template <class T>
concept TemporaryView = requires(T&& value) { std::move(value).view(); };
template <class T>
concept TemporaryProof = requires(T&& value) { std::move(value).basis_proof(); };
static_assert(!TemporaryView<eig::FullSpectrum>);
static_assert(!TemporaryProof<eig::FullSpectrum>);
static_assert(!std::copy_constructible<eig::FullSpectrum>);
static_assert(!std::is_same_v<eig::FullSpectrum, eig::FullPairs>);
static_assert(!std::is_same_v<eig::FullSpectrum, eig::FullCluster>);

[[nodiscard]] eig::Context identity() {
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
[[nodiscard]] eig::FullRequest requirement() {
    return {{{1.e-12, core::ObservationId{20}, core::ScopeId{21}},
             {1.e-12, core::ObservationId{22}, core::ScopeId{23}},
             eig::NotRequested{}},
            lin::RetainedPrecision{core::NumericalGeneration{40}},
            {16, 16, {32, 32, 0, 4}, lin::EstimateAdmission{std::uint64_t{1024} * 1024}},
            {},
            nullptr,
            false};
}
struct ScaleInput {
    std::array<double, 2> components{1., 1.};
    std::array<lin::UnitDimensions, 2> units{};
    std::array<lin::UnitDimensions, 4> coefficient_units{};
    double eigenvalue = 1;
    [[nodiscard]] eig::FullScales view() const {
        return {components, eigenvalue, units, {}, coefficient_units};
    }
};
template <numerics::StorageScalar Scalar>
[[nodiscard]] numerics::PrecisionBuffer<Scalar>
input(numerics::DenseShape shape, std::span<const Scalar> values, core::ResourcePlan& plan) {
    auto matrix = numerics::DenseMatrix<Scalar>::copy(shape, values, plan, operation);
    if (!matrix)
        throw std::runtime_error("Invalid contract-test matrix");
    auto result = numerics::PrecisionBuffer<Scalar>::binary64(std::move(*matrix),
                                                              core::NumericalGeneration{5});
    if (!result)
        throw std::runtime_error("Invalid contract-test precision buffer");
    return std::move(*result);
}
[[nodiscard]] eig::Result<eig::FullSpectrum> owned_diagonal() {
    core::ResourcePlan inputs{1024}, work{std::uint64_t{1024} * 1024};
    const auto original = input<double>({2, 2, 2}, std::array{2., 0., 0., 5.}, inputs);
    const ScaleInput scales;
    return eig::full_spectrum(original, identity(), reporting, scales.view(), requirement(), work);
}
class Consumer final : public eig::FullConsumer {
public:
    mutable unsigned capability_calls = 0;
    mutable unsigned evaluations = 0;
    mutable std::array<std::optional<core::NumericalGeneration>, 8> generations{};
    unsigned refuse_call = 0;
    unsigned failure_call = 0;
    unsigned throw_call = 0;
    bool wrong_generation = false;
    bool wrong_budget = false;
    bool unevaluated = false;
    [[nodiscard]] numerics::SupportedPrecision supported_precision() const noexcept override {
        ++capability_calls;
        return numerics::SupportedPrecision::Binary64Only;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    evaluate(eig::FullPacketView packet, eig::MetricRequirement request) const override {
        ++evaluations;
        if (evaluations <= generations.size())
            generations[evaluations - 1] = packet.generation;
        if (evaluations == throw_call)
            throw std::runtime_error("Contract-test consumer failure");
        if (evaluations == failure_call) {
            core::FailureDetails details;
            details.parent = core::AcceptedStateId{901};
            details.last_accepted = core::AcceptedStateId{902};
            return std::unexpected(
                core::Failure{core::FailureCode::Cancelled, core::Operation::SolveLinear, details});
        }
        const core::CertificateContext context{request.observation, request.scope,
                                               wrong_generation ? core::NumericalGeneration{999}
                                                                : packet.generation,
                                               core::Precision::Binary64};
        const double budget = wrong_budget ? 2 * request.budget : request.budget;
        if (unevaluated)
            return core::ErrorCertificate::unevaluated(context, budget);
        return core::ErrorCertificate::assess(context, budget,
                                              evaluations == refuse_call ? 2 * request.budget : 0.);
    }
};
constexpr std::array checks{
    eig::MetricRequirement{1.e-12, core::ObservationId{70}, core::ScopeId{71}},
    eig::MetricRequirement{1.e-12, core::ObservationId{72}, core::ScopeId{73}}};
} // namespace

TEST_CASE("Full results outlive every borrowed input and safely expose moved-from status",
          "[N4-07][N4-08][property]") {
    auto result = owned_diagonal();
    REQUIRE(result);
    auto moved = std::move(*result);
    CHECK(moved.active());
    CHECK_FALSE(result->active());
    CHECK_FALSE(result->view().has_value());
    CHECK(result->returned_precision() == nullptr);
    CHECK(result->basis_proof() == nullptr);
    CHECK(result->diagonal_mapping().empty());
    CHECK_FALSE(result->consumer_proof(0));
    const auto actual = moved.view();
    if (!actual.has_value()) {
        FAIL("Active packet lost its owned values");
        return;
    }
    REQUIRE(actual->eigenvalues.size() == 2);
    CHECK(std::abs(actual->eigenvalues[0] - Complex{2, 0}) <= 1.e-12);
    CHECK(std::abs(actual->eigenvalues[1] - Complex{5, 0}) <= 5.e-12);
    CHECK(actual->state.parent == reporting.parent);
    CHECK(actual->state.last_accepted == reporting.last_accepted);
    CHECK(actual->generation == core::NumericalGeneration{40});
}

TEST_CASE("Empty full packets perform no callback including capability metadata",
          "[N4-07][N4-02][negative]") {
    core::ResourcePlan inputs{0}, work{0};
    const auto original = input<double>({0, 0, 0}, {}, inputs);
    Consumer consumer;
    auto request = requirement();
    request.limits = {0, 0, {0, 0, 0, 0}, lin::RequireProvenBound{0}};
    request.evaluator = &consumer;
    auto result =
        eig::full_spectrum(original, identity(), reporting, {{}, 1., {}, {}, {}}, request, work);
    REQUIRE(result);
    CHECK(consumer.capability_calls == 0);
    CHECK(consumer.evaluations == 0);
    CHECK(work.allocated_bytes() == 0);
    const auto actual = result->view();
    if (!actual.has_value()) {
        FAIL("Empty successful packet must remain active");
        return;
    }
    CHECK(actual->eigenvalues.empty());
    CHECK(result->basis_proof() != nullptr);
    request.consumers = checks;
    const auto refused =
        eig::full_spectrum(original, identity(), reporting, {{}, 1., {}, {}, {}}, request, work);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().cause.code == core::FailureCode::PrecisionUnresolved);
    CHECK(consumer.capability_calls == 0);
    CHECK(consumer.evaluations == 0);
}

TEST_CASE("Full preflight rejects stale units scales budgets and unsupported words before payload",
          "[N4-07][N4-06][negative]") {
    core::ResourcePlan inputs{1024};
    const auto original = input<double>({2, 2, 2}, std::array{2., 0., 0., 5.}, inputs);
    for (unsigned mode = 0; mode < 8; ++mode) {
        auto context = identity();
        auto request = requirement();
        ScaleInput scales;
        core::ResourcePlan work{0};
        core::FailureCode expected = core::FailureCode::InvalidInput;
        switch (mode) {
        case 0:
            context.cache.matrix = core::NumericalGeneration{999};
            expected = core::FailureCode::StaleGeneration;
            break;
        case 1:
            scales.components[0] = 0;
            break;
        case 2:
            scales.components[1] = std::numeric_limits<double>::quiet_NaN();
            expected = core::FailureCode::Nonfinite;
            break;
        case 3:
            scales.coefficient_units[0][2] = 1;
            break;
        case 4:
            request.basis.residual.budget = -1;
            break;
        case 5:
            request.basis.orthogonality.budget = std::numeric_limits<double>::infinity();
            expected = core::FailureCode::Nonfinite;
            break;
        case 6:
            context.time = std::numeric_limits<double>::quiet_NaN();
            expected = core::FailureCode::Nonfinite;
            break;
        default:
            request.require_physical_coverage = true;
            expected = core::FailureCode::EigenCoverageUnresolved;
            break;
        }
        const auto result =
            eig::full_spectrum(original, context, reporting, scales.view(), request, work);
        REQUIRE_FALSE(result);
        CHECK(result.error().cause.code == expected);
        CHECK(result.error().cause.operation == operation);
        CHECK(result.error().cause.details.parent == reporting.parent);
        CHECK(work.allocated_bytes() == 0);
    }
    auto hi = numerics::DenseReal::copy({2, 2, 2}, std::array{2., 0., 0., 5.}, inputs, operation);
    auto lo =
        numerics::DenseReal::copy({2, 2, 2}, std::array{0., 0x1p-80, 0., 0.}, inputs, operation);
    REQUIRE(hi);
    REQUIRE(lo);
    const auto retained = numerics::PrecisionBuffer<double>::two_word(
        std::move(*hi), std::move(*lo), core::NumericalGeneration{5});
    REQUIRE(retained);
    const ScaleInput scales;
    core::ResourcePlan work{0};
    const auto result =
        eig::full_spectrum(*retained, identity(), reporting, scales.view(), requirement(), work);
    REQUIRE_FALSE(result);
    CHECK(result.error().cause.code == core::FailureCode::PrecisionUnsupported);
    CHECK(work.allocated_bytes() == 0);
    const auto& words = std::get<numerics::TwoWord<double>>(retained->storage());
    CHECK(words.lo.values()[1] == 0x1p-80);
}

TEST_CASE("Full action budgets admit equality and reject before one excess application",
          "[N4-07][N4-02][negative]") {
    core::ResourcePlan inputs{1024};
    const auto original = input<double>({2, 2, 2}, std::array{2., 0., 0., 5.}, inputs);
    const ScaleInput scales;
    auto request = requirement();
    request.limits.actions.forward = 2;
    request.limits.actions.attempts = 1;
    core::ResourcePlan enough{std::uint64_t{1024} * 1024};
    const auto accepted =
        eig::full_spectrum(original, identity(), reporting, scales.view(), request, enough);
    REQUIRE(accepted);
    CHECK(accepted->diagnostics().forward_applications == 2);
    CHECK(accepted->diagnostics().extraction_attempts == 1);
    for (unsigned mode = 0; mode < 2; ++mode) {
        auto limited = request;
        if (mode == 0)
            limited.limits.actions.forward = 1;
        else
            limited.limits.actions.attempts = 0;
        core::ResourcePlan work{std::uint64_t{1024} * 1024};
        const auto result =
            eig::full_spectrum(original, identity(), reporting, scales.view(), limited, work);
        REQUIRE_FALSE(result);
        CHECK(result.error().cause.code == core::FailureCode::ResourceExceeded);
        CHECK(result.error().diagnostics.forward_applications.value_or(0) == (mode == 0 ? 1 : 0));
        CHECK(accepted->diagnostics().forward_applications == 2);
    }
}

TEST_CASE("Full explicit projection requalifies consumers for both actual generations",
          "[N4-07][N4-06][negative]") {
    core::ResourcePlan inputs{1024}, work{std::uint64_t{1024} * 1024};
    const auto original = input<double>({2, 2, 2}, std::array{2., 0., 0., 5.}, inputs);
    const ScaleInput scales;
    Consumer consumer;
    auto request = requirement();
    request.returned =
        lin::Binary64Projection{core::NumericalGeneration{40}, core::NumericalGeneration{41}};
    request.consumers = checks;
    request.evaluator = &consumer;
    request.limits.actions.forward = 4;
    const auto accepted =
        eig::full_spectrum(original, identity(), reporting, scales.view(), request, work);
    REQUIRE(accepted);
    CHECK(consumer.evaluations == 4);
    CHECK(consumer.generations[0] == core::NumericalGeneration{40});
    CHECK(consumer.generations[1] == core::NumericalGeneration{40});
    CHECK(consumer.generations[2] == core::NumericalGeneration{41});
    CHECK(consumer.generations[3] == core::NumericalGeneration{41});
    CHECK(accepted->diagnostics().forward_applications == 4);
    CHECK(accepted->consumer_count() == 2);
    const auto proof = accepted->consumer_proof(1);
    REQUIRE(proof);
    CHECK(proof->context().numerical_generation == core::NumericalGeneration{41});
    CHECK_FALSE(accepted->consumer_proof(-1));
    CHECK_FALSE(accepted->consumer_proof(2));
    for (unsigned mode = 0; mode < 6; ++mode) {
        Consumer bad;
        bad.refuse_call = mode == 0 ? 4 : 0;
        bad.wrong_generation = mode == 1;
        bad.wrong_budget = mode == 2;
        bad.unevaluated = mode == 3;
        bad.throw_call = mode == 4 ? 1 : 0;
        bad.failure_call = mode == 5 ? 1 : 0;
        request.evaluator = &bad;
        core::ResourcePlan retry{std::uint64_t{1024} * 1024};
        const auto result =
            eig::full_spectrum(original, identity(), reporting, scales.view(), request, retry);
        REQUIRE_FALSE(result);
        CHECK(result.error().cause.code ==
              (mode == 5 ? core::FailureCode::Cancelled : core::FailureCode::PrecisionUnresolved));
        CHECK(result.error().cause.operation ==
              (mode == 5 ? core::Operation::SolveLinear : operation));
        CHECK(result.error().cause.details.parent ==
              (mode == 5 ? std::optional{core::AcceptedStateId{901}} : reporting.parent));
        if (mode == 5)
            CHECK(result.error().cause.details.last_accepted == core::AcceptedStateId{902});
        CHECK(accepted->diagnostics().forward_applications == 4);
    }
    request.evaluator = &consumer;
    request.limits.actions.forward = 3;
    core::ResourcePlan shortage{std::uint64_t{1024} * 1024};
    const auto refused =
        eig::full_spectrum(original, identity(), reporting, scales.view(), request, shortage);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().cause.code == core::FailureCode::ResourceExceeded);
    CHECK(refused.error().diagnostics.forward_applications == 3);
}

TEST_CASE("Full selections validate every index and are bound to the original owning packet",
          "[N4-07][N4-08][negative]") {
    const auto first = owned_diagonal();
    const auto second = owned_diagonal();
    REQUIRE(first);
    REQUIRE(second);
    core::ResourcePlan selections{1024};
    const auto valid = eig::FullSelection::make(*first, std::array<core::Index, 1>{0}, selections);
    REQUIRE(valid);
    const std::array invalid{std::array<core::Index, 2>{0, 0}, std::array<core::Index, 2>{-1, 0},
                             std::array<core::Index, 2>{0, 2}};
    for (const auto& indices : invalid) {
        core::ResourcePlan none{0};
        const auto result = eig::FullSelection::make(*first, indices, none);
        REQUIRE_FALSE(result);
        CHECK(result.error().cause.code == core::FailureCode::InvalidInput);
        CHECK(none.allocated_bytes() == 0);
    }
    const eig::FullClusterRequest request{requirement(), {1.e-12, 1.e-12}, false};
    core::ResourcePlan none{0};
    const auto wrong_owner = eig::full_cluster(*second, *valid, request, none);
    REQUIRE_FALSE(wrong_owner);
    CHECK(wrong_owner.error().cause.code == core::FailureCode::StaleGeneration);
    CHECK(none.allocated_bytes() == 0);
}

TEST_CASE("Unrequested full consumers never call capability metadata", "[N4-07][N4-10][negative]") {
    core::ResourcePlan inputs{1024}, work{std::uint64_t{1024} * 1024};
    const auto original = input<double>({2, 2, 2}, std::array{2., 0., 0., 5.}, inputs);
    Consumer consumer;
    auto request = requirement();
    request.evaluator = &consumer;
    const ScaleInput scales;
    const auto result =
        eig::full_spectrum(original, identity(), reporting, scales.view(), request, work);
    REQUIRE(result);
    CHECK(consumer.capability_calls == 0);
    CHECK(consumer.evaluations == 0);
}
