#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <functional>
#include <limits>
#include <new>
#include <span>
#include <stdexcept>
#include <utility>
#include <uvlas/numerics/gmres_solve.hpp>
#include <variant>

namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;
constexpr auto operation = core::Operation::SolveLinear;

Context identity() {
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
core::Failure inner_failure() {
    core::FailureDetails details;
    details.model = core::ModelGeneration{101};
    details.trial = core::TrialId{102};
    details.numerical_generation = core::NumericalGeneration{103};
    details.problem = core::ObservationId{104};
    details.resources = core::ResourceUsage{core::ResourceKind::Bytes, 107, 109};
    return {core::FailureCode::Cancelled, core::Operation::EvaluateAero, details};
}
core::Failure precision_failure() {
    auto failure = inner_failure();
    failure.code = core::FailureCode::PrecisionUnresolved;
    return failure;
}

enum class Fault : std::uint8_t {
    None,
    Returned,
    Allocation,
    Nonfinite,
    NegativeBound,
    ChangeIdentity,
    OmitValue,
    OmitTail,
    OmitError
};

// All fixtures are bounded diagonal systems. FMA gives the exact product remainder
// in this normal-exponent range, independently of the solver's accumulation code.
class Diagonal final : public GmresOperator {
public:
    core::Index n = 2;
    std::array<double, 2> diagonal{1, 1};
    mutable Context state = identity();
    Fault word_fault = Fault::None;
    Fault absolute_fault = Fault::None;
    double fixed_error_upper = 0;
    bool identity_center_in_tail = false;
    mutable bool nonzero_error_returned = false;
    bool validation_failure = false;
    mutable std::uint64_t forwards = 0;
    mutable std::uint64_t absolutes = 0;
    std::uint64_t unresolved_forward = 0;
    std::uint64_t unresolved_absolute = 0;
    mutable std::uint64_t validations = 0;
    [[nodiscard]] core::Index order() const noexcept override {
        return n;
    }
    [[nodiscard]] Context context() const noexcept override {
        return state;
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        ++validations;
        if (validation_failure) {
            return std::unexpected(inner_failure());
        }
        return {};
    }
    [[nodiscard]] core::Result<void> apply_word(std::span<const double> x,
                                                std::span<double> estimate, std::span<double> tail,
                                                std::span<double> error) const override {
        ++forwards;
        if (forwards == unresolved_forward) {
            return std::unexpected(precision_failure());
        }
        if (word_fault == Fault::Returned) {
            return std::unexpected(inner_failure());
        }
        if (word_fault == Fault::Allocation) {
            throw std::bad_alloc{};
        }
        for (std::size_t i = 0; i < x.size(); ++i) {
            const double product = diagonal[i] * x[i];
            if (word_fault != Fault::OmitValue || i + 1 != x.size()) {
                estimate[i] = identity_center_in_tail ? 0 : product;
            }
            if (word_fault != Fault::OmitTail || i + 1 != x.size()) {
                // This alternate representation is enabled only for the exact A=I fixture.
                tail[i] = identity_center_in_tail ? product : std::fma(diagonal[i], x[i], -product);
            }
            if (word_fault != Fault::OmitError || i + 1 != x.size()) {
                error[i] = fixed_error_upper;
                nonzero_error_returned = nonzero_error_returned || error[i] > 0;
            }
        }
        if (!estimate.empty() && word_fault == Fault::Nonfinite) {
            estimate[0] = std::numeric_limits<double>::infinity();
        }
        if (!error.empty() && word_fault == Fault::NegativeBound) {
            error[0] = -1;
        }
        if (word_fault == Fault::ChangeIdentity) {
            state.cache.matrix = core::NumericalGeneration{999};
        }
        return {};
    }
    [[nodiscard]] core::Result<void> absolute_action_lower(std::span<const double> weights,
                                                           std::span<double> lower) const override {
        ++absolutes;
        if (absolutes == unresolved_absolute) {
            return std::unexpected(precision_failure());
        }
        if (absolute_fault == Fault::Returned) {
            return std::unexpected(inner_failure());
        }
        for (std::size_t i = 0; i < weights.size(); ++i) {
            const double product = std::abs(diagonal[i]) * weights[i];
            const double remainder = std::fma(std::abs(diagonal[i]), weights[i], -product);
            if (absolute_fault != Fault::OmitValue || i + 1 != weights.size()) {
                lower[i] = remainder < 0 ? std::nextafter(product, 0.0) : product;
            }
        }
        if (!lower.empty() && absolute_fault == Fault::NegativeBound) {
            lower[0] = -1;
        }
        return {};
    }
};
class Preconditioner final : public FixedPreconditioner {
public:
    core::Index n = 2;
    mutable FixedPreconditionerIdentity state{::identity(), ::identity().cache};
    Fault fault = Fault::None;
    bool zero = false;
    bool validation_failure = false;
    mutable std::uint64_t calls = 0;
    mutable std::uint64_t validations = 0;
    [[nodiscard]] core::Index order() const noexcept override {
        return n;
    }
    [[nodiscard]] FixedPreconditionerIdentity identity() const noexcept override {
        return state;
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        ++validations;
        if (validation_failure) {
            return std::unexpected(inner_failure());
        }
        return {};
    }
    [[nodiscard]] core::Result<void> apply(std::span<const double> residual,
                                           std::span<double> correction) const override {
        ++calls;
        if (fault == Fault::Returned) {
            return std::unexpected(inner_failure());
        }
        for (std::size_t i = 0; i < residual.size(); ++i) {
            if (fault != Fault::OmitValue || i + 1 != residual.size()) {
                correction[i] = zero ? 0 : residual[i];
            }
        }
        if (!correction.empty() && fault == Fault::Nonfinite) {
            correction[0] = std::numeric_limits<double>::quiet_NaN();
        }
        if (fault == Fault::ChangeIdentity) {
            state.factors.backend = core::BackendGeneration{999};
        }
        return {};
    }
};
class Consumer final : public ConsumerEvaluator<double> {
public:
    std::function<void()> hook;
    bool reject_second = false;
    mutable unsigned calls = 0;
    mutable bool saw_low = false;
    [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
        return SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<double>& retained, const PrecisionBuffer<double>& candidate,
           const Context&, const ConsumerCheck& check,
           core::CertificateContext expected) const override {
        ++calls;
        if (hook) {
            hook();
        }
        if (const auto* words = std::get_if<TwoWord<double>>(&retained.storage())) {
            saw_low = words->lo.values()[0] != 0;
        }
        CHECK(expected.numerical_generation == candidate.generation());
        return core::ErrorCertificate::assess(
            expected, check.budget,
            reject_second && check.observation == core::ObservationId{52} ? 2 * check.budget : 0);
    }
};
constexpr std::array requirements{
    ConsumerCheck{CheckKind::Constraint, core::ObservationId{50}, core::ScopeId{51}, 1.e-12},
    ConsumerCheck{CheckKind::Reaction, core::ObservationId{52}, core::ScopeId{53}, 1.e-12}};
struct Fixture {
    Diagonal a;
    Preconditioner p;
    Consumer consumer;
    core::Index restart_width = 64;
    core::ResourcePlan input_plan{65536};
    core::ResourcePlan workspace_plan{1 << 20};
    core::ResourcePlan output_plan{1 << 20};
    std::array<double, 2> unknown{1, 1}, residual{1, 1};
    std::array<UnitDimensions, 2> unknown_units{}, residual_units{};
    explicit Fixture(core::Index n = 2) {
        a.n = n;
        p.n = n;
    }
    [[nodiscard]] Scales scales() const {
        const auto size = static_cast<std::size_t>(a.n);
        return {std::span{unknown}.first(size), std::span{residual}.first(size),
                std::span{unknown_units}.first(size), std::span{residual_units}.first(size),
                ::identity()};
    }
    [[nodiscard]] GmresLimits limits(std::uint64_t quota = 512) const {
        return {a.n, 2, EstimateAdmission{1 << 20}, restart_width, quota};
    }
    [[nodiscard]] Accuracy accuracy(double budget = 1.e-12) const {
        return {Norm::Infinity,         budget, 0, budget, core::ObservationId{20},
                core::ObservationId{21}};
    }
    [[nodiscard]] PrecisionBuffer<double> rhs(std::span<const double> hi, core::Index columns = 1,
                                              std::span<const double> lo = {}) {
        auto high = DenseReal::copy({a.n, columns, a.n}, hi, input_plan, operation);
        REQUIRE(high);
        if (!lo.empty()) {
            auto low = DenseReal::copy({a.n, columns, a.n}, lo, input_plan, operation);
            REQUIRE(low);
            auto result = PrecisionBuffer<double>::two_word(std::move(*high), std::move(*low),
                                                            core::NumericalGeneration{30});
            REQUIRE(result);
            return std::move(*result);
        }
        auto result =
            PrecisionBuffer<double>::binary64(std::move(*high), core::NumericalGeneration{30});
        REQUIRE(result);
        return std::move(*result);
    }
    [[nodiscard]] GmresPrepared prepare() {
        auto result = GmresPrepared::prepare(a, p, scales(), limits(), workspace_plan);
        if (!result) {
            FAIL("GMRES preparation failed");
        }
        return std::move(result).value();
    }
    [[nodiscard]] GmresResult<GmresSolution>
    solve(const PrecisionBuffer<double>& values, GmresPrepared& prepared, std::uint64_t quota = 512,
          std::span<const ConsumerCheck> checks = {}, bool projection = false,
          double budget = 1.e-12, Norm norm = Norm::Infinity) {
        const Precision<double> precision{
            projection ? ReturnedPrecision{Binary64Projection{core::NumericalGeneration{40},
                                                              core::NumericalGeneration{41}}}
                       : ReturnedPrecision{RetainedPrecision{core::NumericalGeneration{40}}},
            checks, consumer};
        auto requested_accuracy = accuracy(budget);
        requested_accuracy.norm = norm;
        return solve_linear(a, p, Rhs<double>{values, ::identity(), core::RhsGeneration{1}},
                            scales(), precision, requested_accuracy, limits(quota), prepared,
                            output_plan);
    }
};
[[nodiscard]] std::span<const double> values(const GmresSolution& result) {
    const auto* binary = std::get_if<Binary64<double>>(&result.solution().values().storage());
    if (binary == nullptr) {
        throw std::runtime_error("Expected binary64 fixture result");
    }
    return binary->values.values();
}
void check_exact_values(const GmresSolution& result, std::span<const double> expected) {
    using Exact = boost::multiprecision::cpp_rational;
    const auto& storage = result.solution().values().storage();
    if (const auto* binary = std::get_if<Binary64<double>>(&storage)) {
        REQUIRE(binary->values.values().size() == expected.size());
        for (std::size_t i = 0; i < expected.size(); ++i) {
            CHECK(Exact{binary->values.values()[i]} == Exact{expected[i]});
        }
        return;
    }
    const auto* words = std::get_if<TwoWord<double>>(&storage);
    if (words == nullptr) {
        throw std::runtime_error("Missing retained representation");
    }
    REQUIRE(words->hi.values().size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        CHECK(Exact{words->hi.values()[i]} + Exact{words->lo.values()[i]} == Exact{expected[i]});
    }
}
void check_counts(const GmresDiagnostics& d, const Fixture& f) {
    CHECK(d.forward_actions == f.a.forwards);
    CHECK(d.absolute_actions == f.a.absolutes);
    CHECK(d.preconditioner_actions == f.p.calls);
    CHECK(d.operator_actions == d.forward_actions + d.absolute_actions);
    CHECK(d.total_action_callbacks == d.operator_actions + d.preconditioner_actions);
}
void reset_counts(Fixture& f) {
    f.a.forwards = 0;
    f.a.absolutes = 0;
    f.p.calls = 0;
}
} // namespace

TEST_CASE("GMRES charges started actions before failure and owns every invocation diagnostic",
          "[negative][N4-05][N4-08]") {
    for (unsigned site = 0; site < 3; ++site) {
        Fixture f{1};
        const auto rhs = f.rhs(std::array{1.0});
        auto prepared = f.prepare();
        if (site == 0) {
            f.a.word_fault = Fault::Returned;
        }
        if (site == 1) {
            f.a.absolute_fault = Fault::Returned;
        }
        if (site == 2) {
            f.p.fault = Fault::Returned;
        }
        const auto failed = f.solve(rhs, prepared);
        REQUIRE_FALSE(failed);
        CHECK(failed.error().failure() == inner_failure());
        const auto snapshot = failed.error().diagnostics();
        check_counts(snapshot, f);
        CHECK(snapshot.total_action_callbacks > 0);
        f.a.word_fault = Fault::None;
        f.a.absolute_fault = Fault::None;
        f.p.fault = Fault::None;
        reset_counts(f);
        const auto success = f.solve(rhs, prepared);
        REQUIRE(success);
        CHECK(failed.error().failure() == inner_failure());
        CHECK(failed.error().diagnostics().total_action_callbacks ==
              snapshot.total_action_callbacks);
        CHECK(failed.error().diagnostics().operator_actions == snapshot.operator_actions);
    }
}

TEST_CASE("GMRES quota is shared across RHS and rejects callbacks outside its remaining budget",
          "[negative][N4-05][N4-09]") {
    Fixture f;
    const auto single = f.rhs(std::array{1.0, 0.0});
    auto prepared = f.prepare();
    const auto control = f.solve(single, prepared);
    REQUIRE(control);
    const auto one_column_cost = control->diagnostics().operator_actions;
    REQUIRE(one_column_cost > 0);
    reset_counts(f);
    const auto multiple = f.rhs(std::array{1.0, 0.0, 0.0, 1.0}, 2);
    const auto limited = f.solve(multiple, prepared, one_column_cost);
    REQUIRE_FALSE(limited);
    CHECK(limited.error().failure().code == core::FailureCode::LinearUnconverged);
    CHECK(limited.error().diagnostics().stop == GmresStopReason::OperatorLimit);
    CHECK(limited.error().diagnostics().operator_actions <= one_column_cost);
    check_counts(limited.error().diagnostics(), f);
    for (const auto quota : {std::uint64_t{0}, std::uint64_t{1}}) {
        reset_counts(f);
        const auto rejected = f.solve(single, prepared, quota);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().failure().code == core::FailureCode::LinearUnconverged);
        CHECK(rejected.error().diagnostics().operator_actions == quota);
        check_counts(rejected.error().diagnostics(), f);
    }
}

TEST_CASE("GMRES freezes the complete preconditioner identity and checks actions on return",
          "[negative][N4-03][N4-05]") {
    for (unsigned mutation = 0; mutation < 4; ++mutation) {
        Fixture f{1};
        const auto rhs = f.rhs(std::array{1.0});
        auto prepared = f.prepare();
        if (mutation == 0) {
            f.p.state.factors.matrix = core::NumericalGeneration{999};
        }
        if (mutation == 1) {
            f.p.state.factors.backend = core::BackendGeneration{999};
        }
        if (mutation == 2) {
            f.p.fault = Fault::ChangeIdentity;
        }
        if (mutation == 3) {
            f.a.word_fault = Fault::ChangeIdentity;
        }
        const auto stale = f.solve(rhs, prepared);
        REQUIRE_FALSE(stale);
        CHECK(stale.error().failure().code == core::FailureCode::StaleGeneration);
        check_counts(stale.error().diagnostics(), f);
    }
}

TEST_CASE("GMRES refuses malformed numerical actions and preserves callback exception failure",
          "[negative][N4-05][N4-08]") {
    for (unsigned mutation = 0; mutation < 5; ++mutation) {
        Fixture f{1};
        const auto rhs = f.rhs(std::array{1.0});
        auto prepared = f.prepare();
        if (mutation == 0) {
            f.a.word_fault = Fault::Nonfinite;
        }
        if (mutation == 1) {
            f.a.word_fault = Fault::NegativeBound;
        }
        if (mutation == 2) {
            f.a.absolute_fault = Fault::NegativeBound;
        }
        if (mutation == 3) {
            f.p.fault = Fault::Nonfinite;
        }
        if (mutation == 4) {
            f.a.word_fault = Fault::Allocation;
        }
        const auto failed = f.solve(rhs, prepared);
        REQUIRE_FALSE(failed);
        constexpr std::array expected{core::FailureCode::Nonfinite, core::FailureCode::InvalidInput,
                                      core::FailureCode::InvalidInput, core::FailureCode::Nonfinite,
                                      core::FailureCode::ResourceExceeded};
        CHECK(failed.error().failure().code == expected[mutation]);
        check_counts(failed.error().diagnostics(), f);
    }
}

TEST_CASE("Zero left preconditioning cannot qualify a nonzero original residual",
          "[negative][N4-05]") {
    Fixture f{1};
    const auto rhs = f.rhs(std::array{1.0});
    f.p.zero = true;
    auto prepared = f.prepare();
    const auto failed = f.solve(rhs, prepared);
    REQUIRE_FALSE(failed);
    CHECK(failed.error().failure().code == core::FailureCode::LinearUnconverged);
    check_counts(failed.error().diagnostics(), f);
}

TEST_CASE("GMRES values and certificates survive all input and prepared owner destruction",
          "[integration][N4-05][N4-06]") {
    const auto result = [] {
        Fixture f{1};
        const auto rhs = f.rhs(std::array{2.0});
        auto prepared = f.prepare();
        return f.solve(rhs, prepared);
    }();
    REQUIRE(result);
    CHECK(std::abs(values(*result)[0] - 2) <= 1.e-12);
    CHECK(result->solution().certificate().unknown_scales()[0] == 1);
    CHECK(result->solution().certificate().refinements() == 0);
    CHECK(result->diagnostics().stop == GmresStopReason::Qualified);
    CHECK(result->diagnostics().restarts == 0);
}

TEST_CASE("GMRES consumer reentry fails without changing the outer owning result",
          "[negative][N4-05][N4-06]") {
    Fixture f{1};
    const auto rhs = f.rhs(std::array{2.0});
    auto prepared = f.prepare();
    bool refused = false;
    f.consumer.hook = [&] {
        const auto inner = f.solve(rhs, prepared);
        refused = !inner && inner.error().failure().code == core::FailureCode::InvalidInput;
    };
    const auto result = f.solve(rhs, prepared, 512, std::span{requirements}.first(1));
    REQUIRE(result);
    CHECK(refused);
    CHECK(f.consumer.calls == 1);
    CHECK(std::abs(values(*result)[0] - 2) <= 1.e-12);
    f.consumer.hook = {};
    f.consumer.reject_second = true;
    const auto rejected = f.solve(rhs, prepared, 512, requirements);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().failure().code == core::FailureCode::PrecisionUnresolved);
    CHECK(f.consumer.calls == 3);
    CHECK(std::abs(values(*result)[0] - 2) <= 1.e-12);
}

TEST_CASE("GMRES keeps a TwoWord RHS for real consumers and rejects an uncertified projection",
          "[verification][negative][N4-05][N4-06]") {
    Fixture f{1};
    const double low = std::ldexp(1.0, -80);
    const auto rhs = f.rhs(std::array{1.0}, 1, std::array{low});
    auto prepared = f.prepare();
    const auto kept = f.solve(rhs, prepared, 512, std::span{requirements}.first(1), false, 1.e-30);
    REQUIRE(kept);
    const auto* pair = std::get_if<TwoWord<double>>(&kept->solution().values().storage());
    if (pair == nullptr) {
        FAIL_CHECK("The nonzero low word must remain owned");
        return;
    }
    CHECK(pair->hi.values()[0] == 1);
    CHECK(pair->lo.values()[0] == low);
    CHECK(f.consumer.saw_low);
    CHECK(kept->solution().certificate().refinements() == 0);
    const auto projected = f.solve(rhs, prepared, 512, {}, true, 1.e-30);
    REQUIRE_FALSE(projected);
    CHECK(projected.error().failure().code == core::FailureCode::PrecisionUnresolved);
    CHECK(pair->lo.values()[0] == low);
}

TEST_CASE("GMRES preallocation enforces strict admission scale contents and move ownership",
          "[negative][N4-03][N4-09]") {
    Fixture f{1};
    const auto rhs = f.rhs(std::array{1.0});
    auto prepared = f.prepare();
    auto strict = f.limits();
    strict.backend_memory = RequireProvenBound{1 << 20};
    const auto no_bound = GmresPrepared::prepare(f.a, f.p, f.scales(), strict, f.workspace_plan);
    REQUIRE_FALSE(no_bound);
    CHECK(no_bound.error().failure().code == core::FailureCode::Unsupported);
    f.unknown[0] = 2;
    const auto changed = f.solve(rhs, prepared);
    REQUIRE_FALSE(changed);
    CHECK(changed.error().failure().code == core::FailureCode::StaleGeneration);
    f.unknown[0] = 1;
    core::ResourcePlan none{0};
    const auto no_output = solve_linear(
        f.a, f.p, Rhs<double>{rhs, ::identity(), core::RhsGeneration{1}}, f.scales(),
        Precision<double>{RetainedPrecision{core::NumericalGeneration{40}}, {}, f.consumer},
        f.accuracy(), f.limits(), prepared, none);
    REQUIRE_FALSE(no_output);
    CHECK(no_output.error().failure().code == core::FailureCode::ResourceExceeded);
    CHECK(none.allocated_bytes() == 0);
    auto moved = std::move(prepared);
    // NOLINTNEXTLINE(bugprone-use-after-move): inactive query is the documented move contract.
    CHECK_FALSE(prepared.active());
    CHECK(moved.active());
    CHECK_FALSE(f.solve(rhs, prepared));
    CHECK(f.solve(rhs, moved));
}

TEST_CASE("GMRES includes nonzero producer error bounds in original qualification",
          "[negative][N4-05][N4-06]") {
    Fixture f{1};
    const auto rhs = f.rhs(std::array{1.0});
    auto prepared = f.prepare();
    // A=I is exact, but the producer deliberately offers only a valid loose enclosure.
    // Discarding that enclosure would falsely qualify the exact candidate immediately.
    f.a.fixed_error_upper = 1.e-6;
    const auto unresolved = f.solve(rhs, prepared, 24);
    REQUIRE_FALSE(unresolved);
    CHECK(f.a.nonzero_error_returned);
    // The exact center residual becomes zero after the identity cycle, while
    // its nonzero producer enclosure still fails qualification. P*r then has
    // zero norm, which is a deterministic inexact-breakdown refusal, not success.
    CHECK(unresolved.error().failure().code == core::FailureCode::LinearUnconverged);
    CHECK(unresolved.error().diagnostics().stop == GmresStopReason::KrylovBreakdown);
    CHECK(unresolved.error().diagnostics().operator_actions <= 24);
    check_counts(unresolved.error().diagnostics(), f);
}

TEST_CASE("GMRES preserves frozen-state failure before numerical callbacks",
          "[negative][N4-03][N4-08]") {
    for (const bool operator_failure : {false, true}) {
        Fixture f{1};
        const auto rhs = f.rhs(std::array{1.0});
        auto prepared = f.prepare();
        if (operator_failure) {
            f.a.validation_failure = true;
        } else {
            f.p.validation_failure = true;
        }
        const auto rejected = f.solve(rhs, prepared);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().failure() == inner_failure());
        CHECK(rejected.error().diagnostics().total_action_callbacks == 0);
        check_counts(rejected.error().diagnostics(), f);
    }
}

TEST_CASE("GMRES empty operators require no numerical callbacks or frozen-state evaluation",
          "[negative][N4-05][N4-09]") {
    Fixture f{0};
    f.a.validation_failure = true;
    f.p.validation_failure = true;
    const auto rhs = f.rhs({}, 2);
    auto prepared = f.prepare();
    const auto result = f.solve(rhs, prepared, 0);
    REQUIRE(result);
    CHECK(result->solution().values().shape().rows == 0);
    CHECK(result->solution().values().shape().columns == 2);
    CHECK(result->solution().certificate().columns().size() == 2);
    CHECK(f.a.validations == 0);
    CHECK(f.p.validations == 0);
    CHECK(result->diagnostics().total_action_callbacks == 0);
}

TEST_CASE("GMRES rejects partially unwritten callback output after successful scratch reuse",
          "[negative][N4-05][N4-08]") {
    for (unsigned site = 0; site < 5; ++site) {
        Fixture f;
        const auto rhs = f.rhs(std::array{1.0, 2.0});
        auto prepared = f.prepare();
        const auto previous = f.solve(rhs, prepared);
        REQUIRE(previous);
        reset_counts(f);
        // Only the final component is omitted. Previous finite scratch must not
        // masquerade as a freshly evaluated estimate, bound or preconditioned value.
        if (site == 0) {
            f.a.word_fault = Fault::OmitValue;
        }
        if (site == 1) {
            f.a.word_fault = Fault::OmitTail;
        }
        if (site == 2) {
            f.a.word_fault = Fault::OmitError;
        }
        if (site == 3) {
            f.a.absolute_fault = Fault::OmitValue;
        }
        if (site == 4) {
            f.p.fault = Fault::OmitValue;
        }
        const auto rejected = f.solve(rhs, prepared);
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().failure().code == core::FailureCode::Nonfinite);
        check_counts(rejected.error().diagnostics(), f);
        CHECK(std::abs(values(*previous)[0] - 1) <= 1.e-12);
        CHECK(std::abs(values(*previous)[1] - 2) <= 1.e-12);
    }
}

TEST_CASE("GMRES Krylov directions consume a legitimate action center stored entirely in its tail",
          "[verification][negative][N4-05][N4-06]") {
    Fixture f{1};
    f.a.identity_center_in_tail = true;
    const auto rhs = f.rhs(std::array{2.0});
    auto prepared = f.prepare();
    const auto result = f.solve(rhs, prepared);
    REQUIRE(result);
    CHECK(std::abs(values(*result)[0] - 2) <= 1.e-12);
    CHECK(result->diagnostics().stop == GmresStopReason::Qualified);
    CHECK(result->diagnostics().restarts == 0);
    check_counts(result->diagnostics(), f);
}

TEST_CASE("GMRES exact solutions remain reachable when initial metric ratios are unrepresentable",
          "[verification][negative][N4-05][N4-06]") {
    // 2^-1063 is a fixed positive binary64 allowance. 1/allowance overflows,
    // although the final exact residual divided by it is representable zero.
    const double tiny_allowance = std::ldexp(1.0, -1063);
    const std::array budgets{0.0, tiny_allowance, 1.e-12};
    for (const auto norm : {Norm::Infinity, Norm::L2}) {
        for (const double rhs_value : {0.0, 1.0}) {
            for (const double budget : budgets) {
                Fixture f{1};
                const std::array input{rhs_value};
                const auto rhs = f.rhs(input);
                auto prepared = f.prepare();
                const auto result = f.solve(rhs, prepared, 512, {}, false, budget, norm);
                REQUIRE(result);
                check_exact_values(*result, input);
                const auto metric = result->solution().certificate().columns()[0];
                CHECK(metric.normalized_original_residual_upper == 0);
                CHECK(metric.componentwise_backward_error_upper == 0);
                CHECK(result->diagnostics().stop == GmresStopReason::Qualified);
                if (rhs_value == 0) {
                    CHECK(result->diagnostics().total_action_callbacks == 0);
                } else {
                    CHECK(result->diagnostics().operator_actions > 0);
                }
                check_counts(result->diagnostics(), f);
            }
        }
    }
}

TEST_CASE("GMRES continues beyond unrepresentable intermediate ratios until its real action limit",
          "[negative][N4-05][N4-06][N4-09]") {
    for (const auto norm : {Norm::Infinity, Norm::L2}) {
        Fixture f;
        f.a.diagonal = {1, 2};
        f.restart_width = 1;
        const auto rhs = f.rhs(std::array{1.0, 1.0});
        auto prepared = f.prepare();
        // A one-direction cycle cannot solve this two-eigenvalue RHS exactly.
        // It uses one Arnoldi action and two qualification actions; the next
        // Arnoldi action reaches quota four. An intermediate ratio/0 must not
        // stop earlier as PrecisionUnresolved or manufacture a qualified metric.
        const auto result = f.solve(rhs, prepared, 4, {}, false, 0, norm);
        REQUIRE_FALSE(result);
        CHECK(result.error().failure().code == core::FailureCode::LinearUnconverged);
        CHECK(result.error().diagnostics().stop == GmresStopReason::OperatorLimit);
        CHECK(result.error().diagnostics().operator_actions == 4);
        CHECK(result.error().diagnostics().restarts > 0);
        check_counts(result.error().diagnostics(), f);
    }
}

TEST_CASE("GMRES zero allowances never qualify nonzero producer uncertainty",
          "[negative][N4-05][N4-06]") {
    for (const auto norm : {Norm::Infinity, Norm::L2}) {
        Fixture f{1};
        f.a.fixed_error_upper = 1.e-6;
        const auto rhs = f.rhs(std::array{1.0});
        auto prepared = f.prepare();
        const auto result = f.solve(rhs, prepared, 24, {}, false, 0, norm);
        REQUIRE_FALSE(result);
        CHECK(result.error().failure().code == core::FailureCode::LinearUnconverged);
        CHECK(result.error().diagnostics().stop == GmresStopReason::KrylovBreakdown);
        CHECK(f.a.nonzero_error_returned);
        check_counts(result.error().diagnostics(), f);
    }
}

TEST_CASE("GMRES preserves genuine precision failures from initial and later numerical callbacks",
          "[negative][N4-05][N4-06][N4-08]") {
    for (unsigned site = 0; site < 4; ++site) {
        Fixture f;
        f.a.diagonal = {1, 2};
        f.restart_width = 1;
        if (site < 2) {
            f.a.unresolved_forward = site + 1;
        } else {
            f.a.unresolved_absolute = site - 1;
        }
        const auto rhs = f.rhs(std::array{1.0, 1.0});
        auto prepared = f.prepare();
        const auto result = f.solve(rhs, prepared, 24);
        REQUIRE_FALSE(result);
        CHECK(result.error().failure() == precision_failure());
        CHECK(result.error().diagnostics().stop == GmresStopReason::ActionFailure);
        if (site < 2) {
            CHECK(f.a.forwards == site + 1);
        } else {
            CHECK(f.a.absolutes == site - 1);
        }
        check_counts(result.error().diagnostics(), f);
    }
}

TEST_CASE("GMRES cannot replace a genuinely unrepresentable original residual with an unmet marker",
          "[negative][N4-05][N4-06]") {
    Fixture f{1};
    const std::array huge{std::numeric_limits<double>::max()};
    // Each input word is finite; their initial residual center cannot be
    // evaluated by this finite-word envelope. No finite residual exists for
    // a Krylov update, unlike a representable nonzero residual divided by tau=0.
    const auto rhs = f.rhs(huge, 1, huge);
    auto prepared = f.prepare();
    const auto result = f.solve(rhs, prepared);
    REQUIRE_FALSE(result);
    CHECK(result.error().failure().code == core::FailureCode::PrecisionUnresolved);
    CHECK(result.error().diagnostics().stop == GmresStopReason::QualificationFailure);
    CHECK(result.error().diagnostics().total_action_callbacks == 0);
    check_counts(result.error().diagnostics(), f);
}
