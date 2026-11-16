#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <uvlas/numerics/eigen_generalized.hpp>

namespace {
using namespace uvlas;
namespace eig = numerics::eigen;
using Exact = boost::multiprecision::cpp_rational;
constexpr auto op = core::Operation::ExtractEigen;
[[nodiscard]] eig::Context context() {
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
[[nodiscard]] eig::MetricRequirement criterion(std::uint64_t id) {
    return {1.e-12, core::ObservationId{id}, core::ScopeId{id + 100}};
}
[[nodiscard]] eig::GeneralizedRequest
request(std::span<const double> scales, std::span<const numerics::linear::UnitDimensions> units,
        double lambda = 1, double mass = 1) {
    return {{{scales, lambda}, mass, units, {}, {}, context()},
            {criterion(20), criterion(21), criterion(22), eig::NotRequested{}},
            {core::Precision::Binary64,
             numerics::linear::RetainedPrecision{core::NumericalGeneration{40}},
             {},
             nullptr},
            {100, 100, {200, 0, 0, 1}, numerics::linear::EstimateAdmission{10000000}},
            {core::AcceptedStateId{30}, core::AcceptedStateId{31}},
            false};
}
template <std::size_t Count>
[[nodiscard]] numerics::DenseReal matrix(core::Index n, const std::array<double, Count>& values,
                                         core::ResourcePlan& plan) {
    auto result = numerics::DenseReal::copy({n, n, n}, values, plan, op);
    REQUIRE(result);
    return std::move(*result);
}
[[nodiscard]] double upper(const core::ErrorCertificate& proof) {
    const auto* achieved = std::get_if<core::Achieved>(&proof.assessment());
    REQUIRE(achieved != nullptr);
    return achieved->upper_bound;
}
void exact_original_check(const eig::GeneralizedSpectrum& result) {
    const auto n = static_cast<std::size_t>(result.basis().shape().rows);
    const auto phi = result.basis().values(), values = result.eigenvalues();
    for (std::size_t column = 0; column < n; ++column) {
        Exact norm = 0, maximum_residual = 0, vector_norm = 0;
        std::size_t pivot = 0;
        for (std::size_t row = 0; row < n; ++row) {
            Exact k = 0, m = 0;
            for (std::size_t j = 0; j < n; ++j) {
                k += Exact{result.stiffness()
                               .values()[j * static_cast<std::size_t>(
                                                 result.stiffness().shape().leading_dimension) +
                                         row]} *
                     Exact{phi[column * n + j]};
                m += Exact{result.mass().values()[j * static_cast<std::size_t>(
                                                          result.mass().shape().leading_dimension) +
                                                  row]} *
                     Exact{phi[column * n + j]};
            }
            const Exact residual =
                abs(Exact{result.component_scales()[row]} * (k - Exact{values[column]} * m) /
                    (Exact{result.mass_scale()} * Exact{result.eigenvalue_scale()}));
            const Exact coordinate =
                abs(Exact{phi[column * n + row]} / Exact{result.component_scales()[row]});
            if (residual > maximum_residual)
                maximum_residual = residual;
            if (coordinate > vector_norm) {
                vector_norm = coordinate;
                pivot = row;
            }
            norm += Exact{phi[column * n + row]} * m / Exact{result.mass_scale()};
        }
        REQUIRE(vector_norm > 0);
        CHECK(phi[column * n + pivot] >= 0);
        const Exact normalized = maximum_residual / vector_norm;
        CHECK(normalized <= Exact{1.e-12});
        CHECK(normalized <= Exact{upper(result.metrics().residual)});
        CHECK(abs(norm - 1) <= Exact{upper(result.metrics().normalization)});
        Exact gram_row_sum = 0;
        for (std::size_t other = 0; other < n; ++other) {
            Exact gram = 0;
            for (std::size_t row = 0; row < n; ++row)
                for (std::size_t j = 0; j < n; ++j)
                    gram += Exact{phi[column * n + row]} *
                            Exact{result.mass()
                                      .values()[j * static_cast<std::size_t>(
                                                        result.mass().shape().leading_dimension) +
                                                row]} *
                            Exact{phi[other * n + j]} / Exact{result.mass_scale()};
            gram_row_sum += abs(gram - (column == other ? 1 : 0));
        }
        CHECK(gram_row_sum <= Exact{1.e-12});
        CHECK(gram_row_sum <= Exact{upper(result.metrics().orthogonality)});
    }
}
class Consumer final : public eig::GeneralizedConsumer {
public:
    mutable int calls = 0;
    int reject_on = -1;
    bool wrong_generation = false;
    bool throw_exception = false;
    bool allocation_failure = false;
    bool callback_failure = false;
    mutable bool saw_projection = false;
    core::Result<core::ErrorCertificate> assess(const eig::GeneralizedSpectrum& retained,
                                                const eig::GeneralizedSpectrum& returned,
                                                const numerics::linear::ConsumerCheck& check,
                                                core::CertificateContext expected) const override {
        ++calls;
        saw_projection = retained.generation() != returned.generation();
        if (allocation_failure)
            throw std::bad_alloc{};
        if (throw_exception)
            throw std::runtime_error("consumer");
        if (callback_failure)
            return std::unexpected(
                core::Failure{core::FailureCode::Cancelled, core::Operation::SolveLinear});
        if (wrong_generation)
            expected.numerical_generation = core::NumericalGeneration{999};
        return core::ErrorCertificate::assess(expected, check.budget,
                                              calls == reject_on ? 2 * check.budget : 0);
    }
};
} // namespace
TEST_CASE("Generalized full extraction preserves signed roots and original mass geometry",
          "[N4-07][verification]") {
    core::ResourcePlan input(10000), plan(1000000);
    const auto k = matrix(2, std::array{-2., -2., -2., 10.}, input);
    const auto m = matrix(2, std::array{2., 2., 2., 5.}, input);
    const std::array d{1., 1.};
    const std::array<numerics::linear::UnitDimensions, 2> units{};
    auto r = request(d, units);
    auto mc = context();
    mc.cache.matrix = core::NumericalGeneration{99};
    auto result = eig::symmetric_generalized_spectrum({k, context()}, {m, mc}, r, plan);
    REQUIRE(result);
    CHECK(result->eigenvalues().size() == 2);
    CHECK(std::abs(result->eigenvalues()[0] + 1) < 1.e-12);
    CHECK(std::abs(result->eigenvalues()[1] - 4) < 4.e-12);
    CHECK(result->mass_context().cache.matrix == core::NumericalGeneration{99});
    CHECK(result->state_context().parent == core::AcceptedStateId{30});
    CHECK(result->mass_qualification_shift() > 0);
    CHECK(result->backend_iteration_limit() == 60);
    CHECK_FALSE(result->diagnostics().backend_iterations.has_value());
    CHECK_FALSE(result->physical_coverage_qualified());
    exact_original_check(*result);
    const auto phi = result->basis().values();
    CHECK(std::abs(phi[1]) < 1.e-12);
    CHECK(std::abs(phi[2] + phi[3]) < 1.e-12);
    CHECK(phi[0] > 0);
    CHECK(phi[3] > 0);
}
TEST_CASE("Nonbinary scales and units are preserved and certified on actual physical vectors",
          "[N4-07][verification]") {
    core::ResourcePlan input(10000), plan(1000000);
    const auto k = matrix(2, std::array{-2., -2., -2., 10.}, input);
    const auto m = matrix(2, std::array{2., 2., 2., 5.}, input);
    const std::array d{3., .2};
    std::array<numerics::linear::UnitDimensions, 2> units{};
    units[0][0] = 1;
    units[1][2] = -1;
    auto r = request(d, units, 7, 3);
    r.scales.eigenvalue_units[2] = -2;
    r.scales.mass_units[0] = 2;
    auto result = eig::symmetric_generalized_spectrum({k, context()}, {m, context()}, r, plan);
    REQUIRE(result);
    CHECK(std::abs(result->eigenvalues()[0] + 1) < 1.e-12);
    CHECK(std::abs(result->eigenvalues()[1] - 4) < 4.e-12);
    exact_original_check(*result);
    const auto phi = result->basis().values();
    CHECK(std::abs(phi[0] - std::sqrt(1.5)) < 1.e-12);
    CHECK(std::abs(phi[2] + 1) < 1.e-12);
    CHECK(std::abs(phi[3] - 1) < 1.e-12);
    CHECK(result->component_units()[0] == units[0]);
    CHECK(result->component_scales()[1] == .2);
    CHECK(result->eigenvalue_units()[2] == -2);
    CHECK(result->mass_units()[0] == 2);
}
TEST_CASE("Tiny generalized eigenvalues keep relative accuracy sign and an exact zero",
          "[N4-07][negative][verification]") {
    constexpr double delta = 0x1p-60;
    core::ResourcePlan input(10000), plan(1000000);
    const auto k = matrix(3, std::array{-delta, 0., 0., 0., 0., 0., 0., 0., 16 * delta}, input);
    const auto m = matrix(3, std::array{1., 0., 0., 0., 4., 0., 0., 0., 16.}, input);
    const std::array d{1., 1., 1.};
    const std::array<numerics::linear::UnitDimensions, 3> units{};
    auto result = eig::symmetric_generalized_spectrum({k, context()}, {m, context()},
                                                      request(d, units), plan);
    REQUIRE(result);
    CHECK(result->eigenvalues()[0] < 0);
    CHECK(result->eigenvalues()[1] == 0);
    CHECK(result->eigenvalues()[2] > 0);
    CHECK(std::abs(result->eigenvalues()[0] / (-delta) - 1) < 1.e-12);
    CHECK(std::abs(result->eigenvalues()[2] / delta - 1) < 1.e-12);
    exact_original_check(*result);
}
TEST_CASE("Empty generalized extraction has zero proofs and no backend work", "[N4-07][negative]") {
    const numerics::DenseReal k, m;
    core::ResourcePlan plan(0);
    auto r = request({}, {});
    r.limits = {0, 0, {0, 0, 0, 0}, numerics::linear::RequireProvenBound{0}};
    auto result = eig::symmetric_generalized_spectrum({k, context()}, {m, context()}, r, plan);
    REQUIRE(result);
    CHECK(result->eigenvalues().empty());
    CHECK(upper(result->metrics().residual) == 0);
    CHECK(upper(result->metrics().normalization) == 0);
    CHECK(upper(result->metrics().orthogonality) == 0);
    CHECK(result->diagnostics().extraction_attempts == 0);
    CHECK(result->backend_admission().admitted_payload_bytes() == 0);
    CHECK(plan.allocated_bytes() == 0);
}
TEST_CASE("Original mass qualification and complete symmetry cannot be replaced by scaled LLT",
          "[N4-07][negative]") {
    const std::array d{1., 1.};
    const std::array<numerics::linear::UnitDimensions, 2> units{};
    core::ResourcePlan input(10000);
    const auto k = matrix(2, std::array{1., 0., 0., -1.}, input);
    for (const auto values : {std::array{1., 0., 0., 0.}, std::array{1., 0., 0., -1.}}) {
        const auto m = matrix(2, values, input);
        core::ResourcePlan plan(1000000);
        const auto result = eig::symmetric_generalized_spectrum({k, context()}, {m, context()},
                                                                request(d, units), plan);
        REQUIRE_FALSE(result);
        CHECK(result.error().cause.code == core::FailureCode::RankUnresolved);
        CHECK(result.error().cause.operation == op);
    }
    const auto m = matrix(2, std::array{1., 0., .25, 1.}, input);
    core::ResourcePlan plan(1000000);
    const auto result = eig::symmetric_generalized_spectrum({k, context()}, {m, context()},
                                                            request(d, units), plan);
    REQUIRE_FALSE(result);
    CHECK(result.error().cause.code == core::FailureCode::InvalidInput);
}
TEST_CASE(
    "Generalized entry rejects stale unsupported and physical coverage requests before allocation",
    "[N4-07][negative]") {
    core::ResourcePlan input(10000);
    const auto k = matrix(1, std::array{2.}, input), m = matrix(1, std::array{1.}, input);
    const std::array d{1.};
    const std::array<numerics::linear::UnitDimensions, 1> units{};
    for (int mode = 0; mode < 4; ++mode) {
        auto r = request(d, units);
        auto mc = context();
        core::ResourcePlan plan(0);
        if (mode == 0)
            mc.trial = core::TrialId{777};
        if (mode == 1)
            r.precision.input = core::Precision::TwoWord;
        if (mode == 2)
            r.require_physical_coverage = true;
        if (mode == 3)
            r.limits.maximum_results = 0;
        const auto result = eig::symmetric_generalized_spectrum({k, context()}, {m, mc}, r, plan);
        REQUIRE_FALSE(result);
        const std::array expected{
            core::FailureCode::StaleGeneration, core::FailureCode::PrecisionUnsupported,
            core::FailureCode::EigenCoverageUnresolved, core::FailureCode::ResourceExceeded};
        CHECK(result.error().cause.code == expected[static_cast<std::size_t>(mode)]);
        CHECK(result.error().cause.details.parent == core::AcceptedStateId{30});
        CHECK(plan.allocated_bytes() == 0);
    }
}
TEST_CASE("Owned generalized results survive moved and destroyed original sources",
          "[N4-07][property]") {
    const std::array d{1.};
    const std::array<numerics::linear::UnitDimensions, 1> units{};
    core::ResourcePlan input(10000), plan(1000000);
    auto k = matrix(1, std::array{2.}, input), m = matrix(1, std::array{4.}, input);
    auto result = eig::symmetric_generalized_spectrum({k, context()}, {m, context()},
                                                      request(d, units), plan);
    REQUIRE(result);
    k.values()[0] = 999;
    m.values()[0] = 777;
    auto moved = std::move(*result);
    CHECK(moved.eigenvalues()[0] == .5);
    CHECK(moved.stiffness().values()[0] == 2);
    CHECK(moved.mass().values()[0] == 4);
    CHECK(moved.basis().values()[0] == .5);
}
TEST_CASE("Projection requalifies its distinct generation and every required consumer",
          "[N4-07][N4-06][negative]") {
    core::ResourcePlan input(10000), plan(1000000);
    const auto k = matrix(1, std::array{2.}, input), m = matrix(1, std::array{4.}, input);
    const std::array d{1.};
    const std::array<numerics::linear::UnitDimensions, 1> units{};
    const std::array checks{
        numerics::linear::ConsumerCheck{numerics::linear::CheckKind::Observation,
                                        core::ObservationId{80}, core::ScopeId{81}, 1.e-12},
        numerics::linear::ConsumerCheck{numerics::linear::CheckKind::Reaction,
                                        core::ObservationId{82}, core::ScopeId{83}, 1.e-12}};
    Consumer consumer;
    auto r = request(d, units);
    r.precision.returned = numerics::linear::Binary64Projection{core::NumericalGeneration{40},
                                                                core::NumericalGeneration{41}};
    r.precision.required_checks = checks;
    r.precision.evaluator = &consumer;
    const auto result =
        eig::symmetric_generalized_spectrum({k, context()}, {m, context()}, r, plan);
    REQUIRE(result);
    CHECK(consumer.calls == 2);
    CHECK(consumer.saw_projection);
    CHECK(result->generation() == core::NumericalGeneration{41});
    CHECK(result->metrics().residual.context().numerical_generation ==
          core::NumericalGeneration{41});
    CHECK(result->consumer_count() == 2);
    CHECK(result->consumer(0));
    CHECK_FALSE(result->consumer(2));
    for (int mode = 0; mode < 4; ++mode) {
        core::ResourcePlan retry(1000000);
        Consumer bad;
        bad.reject_on = mode == 0 ? 2 : -1;
        bad.wrong_generation = mode == 1;
        bad.throw_exception = mode == 2;
        bad.callback_failure = mode == 3;
        r.precision.evaluator = &bad;
        const auto refused =
            eig::symmetric_generalized_spectrum({k, context()}, {m, context()}, r, retry);
        REQUIRE_FALSE(refused);
        CHECK(refused.error().cause.code ==
              (mode == 3 ? core::FailureCode::Cancelled : core::FailureCode::PrecisionUnresolved));
        CHECK(refused.error().cause.operation == (mode == 3 ? core::Operation::SolveLinear : op));
        CHECK(result->generation() == core::NumericalGeneration{41});
        CHECK(result->eigenvalues()[0] == .5);
        CHECK(k.values()[0] == 2);
        CHECK(m.values()[0] == 4);
    }
}
TEST_CASE("Generalized resource bounds and nonfinite data fail without mutating input",
          "[N4-07][N4-02][negative]") {
    core::ResourcePlan input(10000);
    auto k = matrix(1, std::array{2.}, input), m = matrix(1, std::array{4.}, input);
    const std::array d{1.};
    const std::array<numerics::linear::UnitDimensions, 1> units{};
    for (int mode = 0; mode < 4; ++mode) {
        core::ResourcePlan plan(mode == 0 ? 0 : 1000000);
        auto r = request(d, units);
        if (mode == 1)
            r.limits.backend_memory = numerics::linear::EstimateAdmission{0};
        if (mode == 2)
            r.limits.backend_memory = numerics::linear::RequireProvenBound{1000000};
        if (mode == 3)
            k.values()[0] = std::numeric_limits<double>::quiet_NaN();
        const auto result =
            eig::symmetric_generalized_spectrum({k, context()}, {m, context()}, r, plan);
        REQUIRE_FALSE(result);
        const std::array expected{core::FailureCode::ResourceExceeded,
                                  core::FailureCode::ResourceExceeded,
                                  core::FailureCode::Unsupported, core::FailureCode::Nonfinite};
        CHECK(result.error().cause.code == expected[static_cast<std::size_t>(mode)]);
        CHECK(m.values()[0] == 4);
    }
}

namespace {
template <class T>
concept RvalueValues = requires(T value) { std::move(value).eigenvalues(); };
template <class T>
concept RvalueBasis = requires(T value) { std::move(value).basis(); };
template <class T>
concept RvalueMass = requires(T value) { std::move(value).mass(); };
template <class T>
concept RvalueStiffness = requires(T value) { std::move(value).stiffness(); };
template <class T>
concept RvalueMetrics = requires(T value) { std::move(value).metrics(); };
static_assert(!RvalueValues<eig::GeneralizedSpectrum>);
static_assert(!RvalueBasis<eig::GeneralizedSpectrum>);
static_assert(!RvalueMass<eig::GeneralizedSpectrum>);
static_assert(!RvalueStiffness<eig::GeneralizedSpectrum>);
static_assert(!RvalueMetrics<eig::GeneralizedSpectrum>);
} // namespace
TEST_CASE("Generalized values outlive their source plan and all original matrices",
          "[N4-07][property]") {
    auto retained = [] {
        core::ResourcePlan input(10000), plan(1000000);
        const auto k = matrix(1, std::array{2.}, input), m = matrix(1, std::array{4.}, input);
        const std::array d{1.};
        const std::array<numerics::linear::UnitDimensions, 1> units{};
        return eig::symmetric_generalized_spectrum({k, context()}, {m, context()},
                                                   request(d, units), plan);
    }();
    REQUIRE(retained);
    CHECK(retained->eigenvalues()[0] == .5);
    CHECK(retained->stiffness().values()[0] == 2);
    CHECK(retained->mass().values()[0] == 4);
    CHECK(retained->component_scales()[0] == 1);
    exact_original_check(*retained);
}
TEST_CASE("Empty problems do not call a required consumer or invent its achievement",
          "[N4-07][N4-06][negative]") {
    const numerics::DenseReal k, m;
    core::ResourcePlan plan(0);
    Consumer evaluator;
    auto r = request({}, {});
    const std::array checks{
        numerics::linear::ConsumerCheck{numerics::linear::CheckKind::Observation,
                                        core::ObservationId{80}, core::ScopeId{81}, 1.e-12}};
    r.precision.required_checks = checks;
    r.precision.evaluator = &evaluator;
    const auto result =
        eig::symmetric_generalized_spectrum({k, context()}, {m, context()}, r, plan);
    REQUIRE_FALSE(result);
    CHECK(result.error().cause.code == core::FailureCode::PrecisionUnresolved);
    CHECK(evaluator.calls == 0);
    CHECK(plan.allocated_bytes() == 0);
}
TEST_CASE("Generalized extraction ignores finite padding and preserves distinct source strides",
          "[N4-07][verification]") {
    core::ResourcePlan input(10000), plan(1000000);
    const std::array kval{-2., -2., 12345., -2., 10., -54321.};
    const std::array mval{2., 2., 999., 888., 2., 5., -777., -666.};
    auto k = numerics::DenseReal::copy({2, 2, 3}, kval, input, op);
    auto m = numerics::DenseReal::copy({2, 2, 4}, mval, input, op);
    REQUIRE(k);
    REQUIRE(m);
    const std::array d{3., .2};
    const std::array<numerics::linear::UnitDimensions, 2> units{};
    const auto result = eig::symmetric_generalized_spectrum({*k, context()}, {*m, context()},
                                                            request(d, units, 7, 3), plan);
    REQUIRE(result);
    CHECK(std::abs(result->eigenvalues()[0] + 1) < 1.e-12);
    CHECK(std::abs(result->eigenvalues()[1] - 4) < 4.e-12);
    CHECK(result->stiffness().shape().leading_dimension == 3);
    CHECK(result->mass().shape().leading_dimension == 4);
    CHECK(result->stiffness().values()[2] == 12345.);
    exact_original_check(*result);
}
TEST_CASE("Invalid generalized budgets scales and matrix dimensions fail at entry",
          "[N4-07][negative]") {
    core::ResourcePlan input(10000);
    const auto k = matrix(1, std::array{2.}, input), m = matrix(1, std::array{4.}, input);
    const std::array good{1.};
    const std::array<numerics::linear::UnitDimensions, 1> units{};
    for (int mode = 0; mode < 8; ++mode) {
        core::ResourcePlan plan(0);
        auto r = request(good, units);
        auto mc = context();
        const std::array bad{0.};
        if (mode == 0)
            r.accuracy.residual.budget = -1;
        if (mode == 1)
            r.accuracy.normalization.budget = std::numeric_limits<double>::quiet_NaN();
        if (mode == 2)
            r.scales.metric.components = {};
        if (mode == 3)
            r.scales.metric.components = bad;
        if (mode == 4)
            r.scales.mass_scale = std::numeric_limits<double>::infinity();
        if (mode == 5)
            mc.time = std::numeric_limits<double>::quiet_NaN();
        if (mode == 6)
            r.precision.returned = numerics::linear::Binary64Projection{
                core::NumericalGeneration{40}, core::NumericalGeneration{40}};
        if (mode == 7)
            r.accuracy.residual.budget = 1.e-9;
        const auto result = eig::symmetric_generalized_spectrum({k, context()}, {m, mc}, r, plan);
        REQUIRE_FALSE(result);
        const std::array expected{
            core::FailureCode::InvalidInput,    core::FailureCode::Nonfinite,
            core::FailureCode::InvalidInput,    core::FailureCode::InvalidInput,
            core::FailureCode::Nonfinite,       core::FailureCode::Nonfinite,
            core::FailureCode::StaleGeneration, core::FailureCode::InvalidInput};
        CHECK(result.error().cause.code == expected[static_cast<std::size_t>(mode)]);
        CHECK(plan.allocated_bytes() == 0);
    }
    const auto non_square = numerics::DenseReal::copy({1, 2, 1}, std::array{2., 3.}, input, op);
    REQUIRE(non_square);
    core::ResourcePlan plan(0);
    const auto result = eig::symmetric_generalized_spectrum(
        {*non_square, context()}, {m, context()}, request(good, units), plan);
    REQUIRE_FALSE(result);
    CHECK(result.error().cause.code == core::FailureCode::InvalidInput);
}
TEST_CASE("Original mass positivity is not inferred from a numerically positive backend pivot",
          "[N4-07][negative]") {
    core::ResourcePlan input(10000), plan(1000000);
    const auto k = matrix(2, std::array{1., 0., 0., 1.}, input);
    const auto m = matrix(2, std::array{1., 1., 1., 1. + 0x1p-52}, input);
    const std::array d{1., 1.};
    const std::array<numerics::linear::UnitDimensions, 2> units{};
    const auto result = eig::symmetric_generalized_spectrum({k, context()}, {m, context()},
                                                            request(d, units), plan);
    REQUIRE_FALSE(result);
    CHECK(result.error().cause.code == core::FailureCode::RankUnresolved);
    CHECK(result.error().diagnostics.extraction_attempts == 0);
}
TEST_CASE(
    "Generalized original actions and extraction attempts are charged before every application",
    "[N4-07][N4-02][negative]") {
    core::ResourcePlan input(10000);
    const auto k = matrix(2, std::array{-2., -2., -2., 10.}, input);
    const auto m = matrix(2, std::array{2., 2., 2., 5.}, input);
    const std::array d{1., 1.};
    const std::array<numerics::linear::UnitDimensions, 2> units{};
    for (const bool projection : {false, true}) {
        const std::uint64_t necessary = projection ? 4 : 2;
        for (const bool sufficient : {false, true}) {
            core::ResourcePlan plan(1000000);
            auto r = request(d, units);
            r.limits.actions.forward = sufficient ? necessary : necessary - 1;
            if (projection)
                r.precision.returned = numerics::linear::Binary64Projection{
                    core::NumericalGeneration{40}, core::NumericalGeneration{41}};
            const auto result =
                eig::symmetric_generalized_spectrum({k, context()}, {m, context()}, r, plan);
            if (sufficient) {
                REQUIRE(result);
                CHECK(result->diagnostics().forward_applications == necessary);
                CHECK(result->diagnostics().adjoint_applications == 0);
                CHECK(result->diagnostics().inverse_applications == 0);
                CHECK(result->diagnostics().extraction_attempts == 1);
                CHECK(result->limits().actions.forward == necessary);
            } else {
                REQUIRE_FALSE(result);
                CHECK(result.error().cause.code == core::FailureCode::ResourceExceeded);
                CHECK(result.error().diagnostics.forward_applications == necessary - 1);
                const auto& resources = result.error().cause.details.resources;
                if (!resources.has_value()) {
                    FAIL("Missing action charge diagnostic");
                    return;
                }
                CHECK(resources->observed == necessary - 1);
                CHECK(resources->required == necessary);
            }
        }
    }
    core::ResourcePlan plan(1000000);
    auto r = request(d, units);
    r.limits.actions.attempts = 0;
    const auto refused =
        eig::symmetric_generalized_spectrum({k, context()}, {m, context()}, r, plan);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().cause.code == core::FailureCode::ResourceExceeded);
    CHECK(refused.error().diagnostics.extraction_attempts == 0);
    CHECK(refused.error().diagnostics.forward_applications == 0);
}
TEST_CASE("Generalized product timing leaves all numerical bytes unchanged when enabled",
          "[N4-07][N4-09][property]") {
    core::ResourcePlan input(10000), first_plan(1000000), second_plan(1000000);
    const auto k = matrix(2, std::array{-2., -2., -2., 10.}, input);
    const auto m = matrix(2, std::array{2., 2., 2., 5.}, input);
    const std::array d{3., .2};
    const std::array<numerics::linear::UnitDimensions, 2> units{};
    auto r = request(d, units, 7, 3);
    const auto first =
        eig::symmetric_generalized_spectrum({k, context()}, {m, context()}, r, first_plan);
    REQUIRE(first);
    r.measurement = eig::GeneralizedMeasurement::Enabled;
    const auto second =
        eig::symmetric_generalized_spectrum({k, context()}, {m, context()}, r, second_plan);
    REQUIRE(second);
    for (std::size_t i = 0; i < 2; ++i)
        CHECK(std::bit_cast<std::uint64_t>(first->eigenvalues()[i]) ==
              std::bit_cast<std::uint64_t>(second->eigenvalues()[i]));
    for (std::size_t i = 0; i < 4; ++i)
        CHECK(std::bit_cast<std::uint64_t>(first->basis().values()[i]) ==
              std::bit_cast<std::uint64_t>(second->basis().values()[i]));
    CHECK(std::bit_cast<std::uint64_t>(upper(first->metrics().residual)) ==
          std::bit_cast<std::uint64_t>(upper(second->metrics().residual)));
    CHECK(std::bit_cast<std::uint64_t>(upper(first->metrics().normalization)) ==
          std::bit_cast<std::uint64_t>(upper(second->metrics().normalization)));
    CHECK(std::bit_cast<std::uint64_t>(upper(first->metrics().orthogonality)) ==
          std::bit_cast<std::uint64_t>(upper(second->metrics().orthogonality)));
    CHECK(first_plan.allocated_bytes() == second_plan.allocated_bytes());
    CHECK_FALSE(first->timings().enabled);
    CHECK(second->timings().enabled);
    for (const auto sample : first->timings().stages) {
        CHECK(sample.calls == 0);
        CHECK(sample.elapsed_nanoseconds == 0);
    }
    for (std::size_t stage = 0; stage < 4; ++stage)
        CHECK(second->timings().stages[stage].calls == 1);
    CHECK(second->timings().stages[4].calls == 0);
}
TEST_CASE("Both original matrices validate every logical and padding component",
          "[N4-07][negative]") {
    core::ResourcePlan input(10000);
    const std::array d{1., 1.};
    const std::array<numerics::linear::UnitDimensions, 2> units{};
    for (const bool mass_side : {false, true}) {
        auto k = matrix(2, std::array{2., 0., 0., 3.}, input);
        auto m = matrix(2, std::array{1., 0., 0., 1.}, input);
        auto& changed = mass_side ? m : k;
        changed.values()[2] = .125;
        core::ResourcePlan first(0);
        const auto asymmetric = eig::symmetric_generalized_spectrum({k, context()}, {m, context()},
                                                                    request(d, units), first);
        REQUIRE_FALSE(asymmetric);
        CHECK(asymmetric.error().cause.code == core::FailureCode::InvalidInput);
        changed.values()[2] = std::numeric_limits<double>::infinity();
        core::ResourcePlan second(0);
        const auto nonfinite = eig::symmetric_generalized_spectrum({k, context()}, {m, context()},
                                                                   request(d, units), second);
        REQUIRE_FALSE(nonfinite);
        CHECK(nonfinite.error().cause.code == core::FailureCode::Nonfinite);
    }
    const auto m = matrix(1, std::array{1.}, input);
    const std::array padded{2., 0.};
    auto k = numerics::DenseReal::copy({1, 1, 2}, padded, input, op);
    REQUIRE(k);
    k->values()[1] = std::numeric_limits<double>::quiet_NaN();
    const std::array one{1.};
    const std::array<numerics::linear::UnitDimensions, 1> one_unit{};
    core::ResourcePlan plan(0);
    const auto refused = eig::symmetric_generalized_spectrum({*k, context()}, {m, context()},
                                                             request(one, one_unit), plan);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().cause.code == core::FailureCode::Nonfinite);
}
TEST_CASE("Generalized consumer allocation failures preserve a value failure and old results",
          "[N4-07][N4-06][negative]") {
    core::ResourcePlan input(10000), good_plan(1000000), bad_plan(1000000);
    const auto k = matrix(1, std::array{2.}, input), m = matrix(1, std::array{4.}, input);
    const std::array d{1.};
    const std::array<numerics::linear::UnitDimensions, 1> units{};
    const auto accepted = eig::symmetric_generalized_spectrum({k, context()}, {m, context()},
                                                              request(d, units), good_plan);
    REQUIRE(accepted);
    Consumer evaluator;
    evaluator.allocation_failure = true;
    auto r = request(d, units);
    const std::array checks{
        numerics::linear::ConsumerCheck{numerics::linear::CheckKind::Observation,
                                        core::ObservationId{80}, core::ScopeId{81}, 1.e-12}};
    r.precision.required_checks = checks;
    r.precision.evaluator = &evaluator;
    const auto refused =
        eig::symmetric_generalized_spectrum({k, context()}, {m, context()}, r, bad_plan);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().cause.code == core::FailureCode::ResourceExceeded);
    CHECK(refused.error().cause.details.parent == core::AcceptedStateId{30});
    CHECK(accepted->eigenvalues()[0] == .5);
    CHECK(k.values()[0] == 2);
    CHECK(m.values()[0] == 4);
}
