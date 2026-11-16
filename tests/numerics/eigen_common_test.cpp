#include <catch2/catch_test_macros.hpp>

#include <array>
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <uvlas/numerics/eigen_common.hpp>
#include <variant>

namespace {
using namespace uvlas;
namespace eig = numerics::eigen;
using Complex = std::complex<double>;

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
[[nodiscard]] eig::StateDiagnosticContext state() {
    return {core::AcceptedStateId{20}, core::AcceptedStateId{21}};
}
[[nodiscard]] eig::AttemptDiagnostics diagnostics() {
    return {24, 6, 7, 2, 30, 31, 32};
}

template <numerics::StorageScalar Scalar>
[[nodiscard]] Scalar scalar(double real, double imaginary) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return real;
    } else {
        return {real, imaginary};
    }
}
template <numerics::StorageScalar Scalar, std::size_t Size>
void same_words(const std::array<Scalar, Size>& actual, const std::array<Scalar, Size>& expected) {
    for (std::size_t i = 0; i < Size; ++i) {
        if constexpr (std::is_same_v<Scalar, double>) {
            CHECK(std::bit_cast<std::uint64_t>(actual[i]) ==
                  std::bit_cast<std::uint64_t>(expected[i]));
        } else {
            CHECK(std::bit_cast<std::uint64_t>(actual[i].real()) ==
                  std::bit_cast<std::uint64_t>(expected[i].real()));
            CHECK(std::bit_cast<std::uint64_t>(actual[i].imag()) ==
                  std::bit_cast<std::uint64_t>(expected[i].imag()));
        }
    }
}

template <numerics::StorageScalar Scalar>
void verify_views() {
    using View = eig::WordView<Scalar>;
    STATIC_REQUIRE_FALSE(std::is_default_constructible_v<View>);
    STATIC_REQUIRE(std::is_copy_constructible_v<View>);
    const std::array high{scalar<Scalar>(1, 2), scalar<Scalar>(0, -3), scalar<Scalar>(4, 0)};
    const std::array low{scalar<Scalar>(-1, -2), scalar<Scalar>(std::ldexp(1., -80), 0),
                         scalar<Scalar>(0, std::ldexp(1., -70))};
    const auto binary = View::make(high, {}, core::Precision::Binary64);
    REQUIRE(binary);
    CHECK(binary->size() == 3);
    CHECK(binary->precision() == core::Precision::Binary64);
    CHECK(binary->high().data() == high.data());
    CHECK(binary->low().empty());
    const auto words = View::make(high, low, core::Precision::TwoWord);
    REQUIRE(words);
    CHECK(words->size() == 3);
    CHECK(words->precision() == core::Precision::TwoWord);
    CHECK(words->high().data() == high.data());
    CHECK(words->low().data() == low.data());
    STATIC_REQUIRE(std::is_same_v<decltype(words->high()[0]), const Scalar&>);
    STATIC_REQUIRE(std::is_same_v<decltype(words->low()[0]), const Scalar&>);
    const auto copied = *words;
    CHECK(copied.high().data() == high.data());
    CHECK(copied.low().data() == low.data());
    for (std::size_t i = 0; i < high.size(); ++i) {
        CHECK(copied.high()[i] == high[i]);
        CHECK(copied.low()[i] == low[i]);
    }
    for (const auto precision : {core::Precision::Binary64, core::Precision::TwoWord}) {
        const auto empty = View::make({}, {}, precision);
        REQUIRE(empty);
        CHECK(empty->size() == 0);
        CHECK(empty->high().empty());
        CHECK(empty->low().empty());
        CHECK(empty->precision() == precision);
    }
    STATIC_REQUIRE(noexcept(View::make({}, {}, core::Precision::Binary64)));
}

template <numerics::StorageScalar Scalar>
void verify_invalid_views() {
    using View = eig::WordView<Scalar>;
    const std::array high{scalar<Scalar>(1, 2), scalar<Scalar>(3, 4)};
    const std::array low{scalar<Scalar>(0, 1)};
    const auto unexpected_low = View::make(high, low, core::Precision::Binary64);
    REQUIRE_FALSE(unexpected_low);
    CHECK(unexpected_low.error().code == core::FailureCode::InvalidInput);
    CHECK(unexpected_low.error().operation == core::Operation::ExtractEigen);
    const auto unequal_words = View::make(high, low, core::Precision::TwoWord);
    REQUIRE_FALSE(unequal_words);
    CHECK(unequal_words.error().code == core::FailureCode::InvalidInput);
    const auto missing_words = View::make(high, {}, core::Precision::TwoWord);
    REQUIRE_FALSE(missing_words);
    CHECK(missing_words.error().code == core::FailureCode::InvalidInput);
    const auto unknown = View::make(high, {}, static_cast<core::Precision>(255));
    REQUIRE_FALSE(unknown);
    CHECK(unknown.error().code == core::FailureCode::PrecisionUnsupported);
}

template <numerics::StorageScalar Scalar>
void verify_nonfinite_words() {
    using View = eig::WordView<Scalar>;
    const std::array bad{std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity()};
    for (const auto value : bad) {
        for (const bool corrupt_low : {false, true}) {
            for (const bool imaginary : {false, true}) {
                if constexpr (std::is_same_v<Scalar, double>) {
                    if (imaginary) {
                        continue; // A real scalar has one component, not an imaginary lane.
                    }
                }
                for (std::size_t i = 0; i < 3; ++i) {
                    std::array high{scalar<Scalar>(1, 2), scalar<Scalar>(3, 4),
                                    scalar<Scalar>(5, 6)};
                    std::array low{scalar<Scalar>(0, 1), scalar<Scalar>(0, 2),
                                   scalar<Scalar>(0, 3)};
                    auto& target = corrupt_low ? low[i] : high[i];
                    target = imaginary ? scalar<Scalar>(1, value) : scalar<Scalar>(value, 1);
                    const auto high_before = high, low_before = low;
                    const auto rejected = View::make(high, low, core::Precision::TwoWord);
                    REQUIRE_FALSE(rejected);
                    CHECK(rejected.error().code == core::FailureCode::Nonfinite);
                    CHECK(rejected.error().operation == core::Operation::ExtractEigen);
                    same_words(high, high_before);
                    same_words(low, low_before);
                }
            }
        }
    }
}

void check_completed_context(const eig::EigenFailure& failure) {
    CHECK(failure.cause.details.model == core::ModelGeneration{1});
    CHECK(failure.cause.details.trial == core::TrialId{11});
    CHECK(failure.cause.details.numerical_generation == core::NumericalGeneration{5});
    CHECK(failure.cause.details.parent == core::AcceptedStateId{20});
    CHECK(failure.cause.details.last_accepted == core::AcceptedStateId{21});
}

struct EigenCallback {
    eig::Result<int> operator()() const {
        return 1;
    }
};
struct CoreCallback {
    core::Result<int> operator()() const {
        return 1;
    }
};
template <class Function>
concept EigenGuardCallback = requires(Function function) {
    eig::guard_entry<int>(context(), state(), diagnostics(), function);
};
} // namespace

TEST_CASE("Eigen word views retain both borrowed planes and accept typed emptiness",
          "[unit][property][N4-02][N4-06][N4-07]") {
    verify_views<double>();
    verify_views<Complex>();
}
TEST_CASE("Eigen word views reject incompatible lengths and unknown representations",
          "[negative][N4-02][N4-06][N4-07]") {
    verify_invalid_views<double>();
    verify_invalid_views<Complex>();
}
TEST_CASE("Every eigen word plane component rejects NaN and both infinities without mutation",
          "[negative][N4-06][N4-07]") {
    verify_nonfinite_words<double>();
    verify_nonfinite_words<Complex>();
}
TEST_CASE("Optional eigen backward requirements preserve absence and the exact requested budget",
          "[unit][N4-07][N4-08]") {
    const eig::BackwardRequirement absent;
    CHECK(std::holds_alternative<eig::NotRequested>(absent));
    const eig::BackwardRequirement requested{
        eig::MetricRequirement{-0.0, core::ObservationId{41}, core::ScopeId{42}}};
    const auto* metric = std::get_if<eig::MetricRequirement>(&requested);
    if (metric == nullptr) {
        FAIL_CHECK("An explicit backward requirement must remain requested");
        return;
    }
    CHECK(std::bit_cast<std::uint64_t>(metric->budget) == std::bit_cast<std::uint64_t>(-0.0));
    CHECK(metric->observation == core::ObservationId{41});
    CHECK(metric->scope == core::ScopeId{42});
}
TEST_CASE("Eigen failure completion fills only known missing reporting identities",
          "[unit][negative][N4-07][N4-08]") {
    const core::Failure empty{core::FailureCode::PrecisionUnresolved, core::Operation::Recover};
    const auto completed = eig::complete_failure(empty, context(), state(), diagnostics());
    CHECK(completed.cause.code == core::FailureCode::PrecisionUnresolved);
    CHECK(completed.cause.operation == core::Operation::Recover);
    check_completed_context(completed);
    CHECK(completed.diagnostics == diagnostics());
    CHECK_FALSE(completed.cause.details.resources.has_value());
    const auto unknown_parent = eig::complete_failure(empty, context(), {});
    CHECK_FALSE(unknown_parent.cause.details.parent.has_value());
    CHECK_FALSE(unknown_parent.cause.details.last_accepted.has_value());
    CHECK(unknown_parent.diagnostics == eig::AttemptDiagnostics{});

    core::FailureDetails supplied{};
    supplied.model = core::ModelGeneration{101};
    supplied.parent = core::AcceptedStateId{102};
    supplied.trial = core::TrialId{103};
    supplied.numerical_generation = core::NumericalGeneration{104};
    supplied.component = core::ComponentId{105};
    supplied.problem = core::ObservationId{106};
    supplied.resources = core::ResourceUsage{core::ResourceKind::Bytes, 108, 107};
    supplied.last_accepted = core::AcceptedStateId{109};
    const core::Failure inner{core::FailureCode::ResourceExceeded, core::Operation::EvaluateAero,
                              supplied};
    const auto preserved = eig::complete_failure(inner, context(), state(), diagnostics());
    CHECK(preserved.cause == inner);
    CHECK(preserved.diagnostics == diagnostics());
    CHECK(inner.details == supplied);
}
TEST_CASE("Eigen failure counters distinguish unreported values from a measured zero",
          "[unit][property][N4-07][N4-08]") {
    using Counter = std::optional<std::uint64_t> eig::AttemptDiagnostics::*;
    // Each named field is independent; ordinal iteration could miss a swapped field.
    const std::array<Counter, 7> counters{&eig::AttemptDiagnostics::requested_count,
                                          &eig::AttemptDiagnostics::converged_count,
                                          &eig::AttemptDiagnostics::backend_iterations,
                                          &eig::AttemptDiagnostics::extraction_attempts,
                                          &eig::AttemptDiagnostics::forward_applications,
                                          &eig::AttemptDiagnostics::adjoint_applications,
                                          &eig::AttemptDiagnostics::inverse_applications};
    const auto fallback = diagnostics();
    for (const auto counter : counters) {
        eig::AttemptDiagnostics original{};
        original.*counter = 0;
        const eig::EigenFailure inner{
            core::Failure{core::FailureCode::Cancelled, core::Operation::ExtractEigen}, original};
        auto expected = fallback;
        expected.*counter = 0;
        const auto completed = eig::complete_failure(inner, context(), state(), fallback);
        CHECK(completed.diagnostics == expected);
        CHECK(inner.diagnostics == original);
    }
    auto later = diagnostics();
    const auto owned = eig::complete_failure(
        core::Failure{core::FailureCode::Cancelled, core::Operation::ExtractEigen}, context(),
        state(), later);
    later = {};
    CHECK(owned.diagnostics == diagnostics());
    CHECK(later == eig::AttemptDiagnostics{});
    STATIC_REQUIRE(std::is_nothrow_copy_constructible_v<eig::EigenFailure>);
    STATIC_REQUIRE(std::is_nothrow_move_constructible_v<eig::EigenFailure>);
}
TEST_CASE("Eigen entry guards accept only eigen results and retain success or original failure",
          "[unit][negative][N4-07][N4-08]") {
    STATIC_REQUIRE(EigenGuardCallback<EigenCallback>);
    STATIC_REQUIRE_FALSE(EigenGuardCallback<CoreCallback>);
    const auto success = eig::guard_entry<int>(context(), state(), diagnostics(),
                                               []() -> eig::Result<int> { return 42; });
    REQUIRE(success);
    CHECK(*success == 42);
    const auto empty =
        eig::guard_entry<void>(context(), state(), {}, []() -> eig::Result<void> { return {}; });
    REQUIRE(empty);
    const core::Failure inner{core::FailureCode::StaleGeneration, core::Operation::Commit};
    const eig::EigenFailure original{inner, {0, 0, {}, {}, 0, {}, {}}};
    const auto failed =
        eig::guard_entry<int>(context(), state(), diagnostics(),
                              [&]() -> eig::Result<int> { return std::unexpected(original); });
    REQUIRE_FALSE(failed);
    CHECK(failed.error().cause.code == inner.code);
    CHECK(failed.error().cause.operation == inner.operation);
    check_completed_context(failed.error());
    auto expected = diagnostics();
    expected.requested_count = 0;
    expected.converged_count = 0;
    expected.forward_applications = 0;
    CHECK(failed.error().diagnostics == expected);
    CHECK(original.cause == inner);
}
TEST_CASE("Eigen entry guards convert allocation and backend exceptions into owned failures",
          "[negative][N4-07][N4-08]") {
    for (const bool allocation : {false, true}) {
        const auto failure =
            eig::guard_entry<int>(context(), state(), diagnostics(), [&]() -> eig::Result<int> {
                if (allocation)
                    throw std::bad_alloc{};
                throw std::runtime_error("backend failed");
            });
        REQUIRE_FALSE(failure);
        CHECK(failure.error().cause.code == (allocation ? core::FailureCode::ResourceExceeded
                                                        : core::FailureCode::FactorizationFailed));
        CHECK(failure.error().cause.operation == core::Operation::ExtractEigen);
        check_completed_context(failure.error());
        CHECK(failure.error().diagnostics == diagnostics());
    }
    const auto nonstandard = eig::guard_entry<void>(context(), state(), diagnostics(),
                                                    []() -> eig::Result<void> { throw 17; });
    REQUIRE_FALSE(nonstandard);
    CHECK(nonstandard.error().cause.code == core::FailureCode::FactorizationFailed);
    CHECK(nonstandard.error().diagnostics == diagnostics());
    STATIC_REQUIRE(noexcept(
        eig::guard_entry<int>(std::declval<const eig::Context&>(), eig::StateDiagnosticContext{},
                              std::declval<const eig::AttemptDiagnostics&>(), EigenCallback{})));
}
