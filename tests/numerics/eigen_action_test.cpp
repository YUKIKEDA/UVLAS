#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <new>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <uvlas/numerics/eigen_action.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::eigen;
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
enum class Fault : std::uint8_t {
    None,
    High,
    Low,
    Error,
    Negative,
    Imaginary,
    Return,
    Throw,
    Allocate
};
template <StorageScalar Scalar>
class Provider final : public ActionProvider<Scalar> {
public:
    core::Index n = 2;
    mutable ActionStamp identity{context(), core::NumericalGeneration{17}, {0, 0}};
    mutable std::uint64_t calls = 0;
    mutable std::uint64_t validations = 0;
    mutable std::array<std::uint64_t, 3> kinds{};
    Fault fault = Fault::None;
    bool change_stamp = false;
    bool change_validation_stamp = false;
    SupportedPrecision support = SupportedPrecision::Binary64AndTwoWord;
    std::function<void()> hook;
    std::function<void()> metadata_hook;
    int metadata_site = -1;
    core::Failure inner = [] {
        core::FailureDetails details;
        details.numerical_generation = core::NumericalGeneration{99};
        details.resources = core::ResourceUsage{core::ResourceKind::Bytes, 101, 37};
        return core::Failure{core::FailureCode::Cancelled, core::Operation::EvaluateAero, details};
    }();
    [[nodiscard]] core::Index order() const noexcept override {
        if (metadata_site == 0 && metadata_hook)
            metadata_hook();
        return n;
    }
    [[nodiscard]] ActionStamp stamp() const noexcept override {
        if (metadata_site == 1 && metadata_hook)
            metadata_hook();
        return identity;
    }
    [[nodiscard]] SupportedPrecision supported_precision(ActionKind) const noexcept override {
        if (metadata_site == 2 && metadata_hook)
            metadata_hook();
        return support;
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        ++validations;
        if (change_validation_stamp)
            identity.transformation = core::NumericalGeneration{18};
        return {};
    }
    [[nodiscard]] core::Result<void> apply(ActionKind kind, WordView<Scalar> input,
                                           std::span<Scalar> high, std::span<Scalar> low,
                                           std::span<double> error) const override {
        ++calls;
        ++kinds[static_cast<std::size_t>(kind)];
        if (hook)
            hook();
        if (fault == Fault::Return)
            return std::unexpected(inner);
        if (fault == Fault::Throw)
            throw std::runtime_error("provider failure");
        if (fault == Fault::Allocate)
            throw std::bad_alloc{};
        if (change_stamp)
            identity.transformation = core::NumericalGeneration{19};
        const Scalar factor = kind == ActionKind::Inverse ? Scalar{0.5} : Scalar{2};
        for (std::size_t i = 0; i < input.size(); ++i) {
            if (fault != Fault::High || i + 1 != input.size())
                high[i] = factor * input.high()[i];
            if (fault != Fault::Low || i + 1 != input.size())
                low[i] = input.low().empty() ? Scalar{} : factor * input.low()[i];
            if (fault != Fault::Error || i + 1 != input.size())
                error[i] = 0;
        }
        if (fault == Fault::Negative)
            error.back() = -1;
        if constexpr (std::same_as<Scalar, std::complex<double>>) {
            if (fault == Fault::Imaginary)
                high.back().imag(std::numeric_limits<double>::infinity());
        }
        return {};
    }
};
constexpr ActionLimits ample{12, 12, 12, 3};
} // namespace

TEST_CASE("Eigen actions retain both numerical planes and independent callback counts",
          "[unit][N4-06][N4-07][N4-09]") {
    const Provider<double> provider;
    auto session = ActionSession<double>::make(provider, ample);
    REQUIRE(session);
    REQUIRE(session->start_attempt());
    const std::array high{2.0, -3.0}, low{0x1p-80, -0x1p-82};
    const auto input = WordView<double>::make(high, low, core::Precision::TwoWord);
    REQUIRE(input);
    std::array<double, 2> out{}, tail{}, error{};
    for (const auto kind : {ActionKind::Forward, ActionKind::Adjoint, ActionKind::Inverse}) {
        REQUIRE(session->apply(provider, kind, *input, out, tail, error));
        if (kind == ActionKind::Inverse) {
            CHECK(out == std::array{1.0, -1.5});
            CHECK(tail == std::array{0x1p-81, -0x1p-83});
        } else {
            CHECK(out == std::array{4.0, -6.0});
            CHECK(tail == std::array{0x1p-79, -0x1p-81});
        }
        CHECK(error == std::array{0.0, 0.0});
    }
    const auto counts = session->diagnostics();
    CHECK(counts.extraction_attempts == 1);
    CHECK(counts.forward_applications == 1);
    CHECK(counts.adjoint_applications == 1);
    CHECK(counts.inverse_applications == 1);
    CHECK_FALSE(counts.backend_iterations);
    CHECK_FALSE(counts.converged_count);
    CHECK(high == std::array{2.0, -3.0});
    CHECK(low == std::array{0x1p-80, -0x1p-82});
    CHECK_FALSE(session->first_failure());
}

TEST_CASE("Complex eigen actions consume real imaginary and trailing components",
          "[unit][N4-06][N4-07]") {
    using Complex = std::complex<double>;
    Provider<Complex> provider;
    auto session = ActionSession<Complex>::make(provider, ample);
    REQUIRE(session);
    const std::array high{Complex{1, 2}, Complex{-3, 4}};
    const std::array low{Complex{0x1p-80, -0x1p-81}, Complex{-0x1p-82, 0x1p-83}};
    const auto input = WordView<Complex>::make(high, low, core::Precision::TwoWord);
    REQUIRE(input);
    std::array<Complex, 2> out{}, tail{};
    std::array<double, 2> error{};
    REQUIRE(session->apply(provider, ActionKind::Adjoint, *input, out, tail, error));
    CHECK(out == std::array{Complex{2, 4}, Complex{-6, 8}});
    CHECK(tail == std::array{Complex{0x1p-79, -0x1p-80}, Complex{-0x1p-81, 0x1p-82}});
    provider.fault = Fault::Imaginary;
    const auto failed = session->apply(provider, ActionKind::Adjoint, *input, out, tail, error);
    REQUIRE_FALSE(failed);
    CHECK(failed.error().cause.code == core::FailureCode::Nonfinite);
    CHECK(failed.error().diagnostics.adjoint_applications == 2);
}

TEST_CASE("Every eigen action quota rejects before the excess callback and latches failure",
          "[negative][N4-07][N4-09]") {
    for (const auto kind : {ActionKind::Forward, ActionKind::Adjoint, ActionKind::Inverse}) {
        const Provider<double> provider;
        auto session = ActionSession<double>::make(provider, {1, 1, 1, 1});
        REQUIRE(session);
        const std::array high{1.0, 2.0};
        const auto input = WordView<double>::make(high, {}, core::Precision::Binary64);
        REQUIRE(input);
        std::array<double, 2> out{}, tail{}, error{};
        REQUIRE(session->apply(provider, kind, *input, out, tail, error));
        const auto failed = session->apply(provider, kind, *input, out, tail, error);
        REQUIRE_FALSE(failed);
        CHECK(failed.error().cause.code == core::FailureCode::ResourceExceeded);
        CHECK(failed.error().cause.details.resources ==
              core::ResourceUsage{core::ResourceKind::OperatorApplications, 2, 1});
        CHECK(provider.calls == 1);
        const auto other = session->apply(provider, ActionKind::Forward, *input, out, tail, error);
        REQUIRE_FALSE(other);
        CHECK(other.error() == failed.error());
        CHECK(provider.calls == 1);
        CHECK(out == (kind == ActionKind::Inverse ? std::array{0.5, 1.0} : std::array{2.0, 4.0}));
    }
    const Provider<double> provider;
    auto attempts = ActionSession<double>::make(provider, {0, 0, 0, 1});
    REQUIRE(attempts);
    REQUIRE(attempts->start_attempt());
    const auto no_retry = attempts->start_attempt();
    REQUIRE_FALSE(no_retry);
    CHECK(no_retry.error().diagnostics.extraction_attempts == 1);
    CHECK(no_retry.error().cause.details.resources ==
          core::ResourceUsage{core::ResourceKind::Iterations, 2, 1});
    auto zero = ActionSession<double>::make(provider, {0, 0, 0, 0});
    REQUIRE(zero);
    CHECK_FALSE(zero->start_attempt());
}

TEST_CASE("Eigen actions poison every output plane before entry and check bound signs",
          "[negative][N4-06][N4-07]") {
    for (const auto fault : {Fault::High, Fault::Low, Fault::Error, Fault::Negative}) {
        Provider<double> provider;
        auto session = ActionSession<double>::make(provider, ample);
        REQUIRE(session);
        provider.fault = fault;
        const std::array high{1.0, 2.0};
        const auto input = WordView<double>::make(high, {}, core::Precision::Binary64);
        REQUIRE(input);
        std::array<double, 2> out{9, 9}, tail{9, 9}, error{9, 9};
        const auto failed = session->apply(provider, ActionKind::Forward, *input, out, tail, error);
        REQUIRE_FALSE(failed);
        CHECK(failed.error().cause.code == (fault == Fault::Negative
                                                ? core::FailureCode::InvalidInput
                                                : core::FailureCode::Nonfinite));
        CHECK(failed.error().diagnostics.forward_applications == 1);
    }
}

TEST_CASE("Eigen actions reject input alias stale state and borrowed-value mutations before charge",
          "[negative][N4-07][N4-08]") {
    for (int site = 0; site < 9; ++site) {
        const Provider<double> provider;
        auto session = ActionSession<double>::make(provider, ample);
        REQUIRE(session);
        std::array high{1.0, 2.0}, input_low{0x1p-80, 0x1p-81};
        const auto input = WordView<double>::make(high, input_low, core::Precision::TwoWord);
        REQUIRE(input);
        std::array<double, 2> out{}, tail{}, error{};
        auto out_span = std::span{out};
        auto tail_span = std::span{tail};
        auto error_span = std::span{error};
        if (site == 0)
            out_span = high;
        if (site == 1)
            tail_span = high;
        if (site == 2)
            error_span = high;
        if (site == 3)
            out_span = input_low;
        if (site == 4)
            tail_span = input_low;
        if (site == 5)
            error_span = input_low;
        if (site == 6)
            tail_span = out;
        if (site == 7)
            error_span = out;
        if (site == 8)
            error_span = tail;
        const auto failed =
            session->apply(provider, ActionKind::Forward, *input, out_span, tail_span, error_span);
        REQUIRE_FALSE(failed);
        CHECK(failed.error().cause.code == core::FailureCode::InvalidInput);
        CHECK(provider.calls == 0);
    }
    const Provider<double> provider;
    auto session = ActionSession<double>::make(provider, ample);
    REQUIRE(session);
    std::array high{1.0, 2.0};
    const auto input = WordView<double>::make(high, {}, core::Precision::Binary64);
    REQUIRE(input);
    high[1] = std::numeric_limits<double>::infinity();
    std::array<double, 2> out{}, tail{}, error{};
    const auto bad = session->apply(provider, ActionKind::Forward, *input, out, tail, error);
    REQUIRE_FALSE(bad);
    CHECK(bad.error().cause.code == core::FailureCode::Nonfinite);
    CHECK(provider.calls == 0);
}

TEST_CASE("Eigen stamps time and transformation are checked before and after actions",
          "[negative][N4-07][N4-08]") {
    for (int site = 0; site < 4; ++site) {
        Provider<double> provider;
        if (site == 0)
            provider.identity.original.time = std::numeric_limits<double>::infinity();
        if (site == 1)
            provider.identity.shift.imag(std::numeric_limits<double>::quiet_NaN());
        if (site == 2)
            provider.change_validation_stamp = true;
        auto session = ActionSession<double>::make(provider, ample);
        if (site != 3) {
            REQUIRE_FALSE(session);
            CHECK(session.error().cause.code ==
                  (site == 2 ? core::FailureCode::StaleGeneration : core::FailureCode::Nonfinite));
        } else {
            REQUIRE(session);
            provider.change_stamp = true;
            const std::array high{1.0, 2.0};
            const auto input = WordView<double>::make(high, {}, core::Precision::Binary64);
            REQUIRE(input);
            std::array<double, 2> out{}, tail{}, error{};
            const auto failed =
                session->apply(provider, ActionKind::Forward, *input, out, tail, error);
            REQUIRE_FALSE(failed);
            CHECK(failed.error().cause.code == core::FailureCode::StaleGeneration);
            CHECK(provider.calls == 1);
        }
    }
}

TEST_CASE("Eigen callback failure and exceptions preserve first owned diagnostics",
          "[negative][N4-07][N4-08]") {
    for (const auto fault : {Fault::Return, Fault::Throw, Fault::Allocate}) {
        Provider<double> provider;
        auto session = ActionSession<double>::make(
            provider, ample, {core::AcceptedStateId{71}, core::AcceptedStateId{72}});
        REQUIRE(session);
        provider.fault = fault;
        const std::array high{1.0, 2.0};
        const auto input = WordView<double>::make(high, {}, core::Precision::Binary64);
        REQUIRE(input);
        std::array<double, 2> out{}, tail{}, error{};
        const auto failed = session->apply(provider, ActionKind::Inverse, *input, out, tail, error);
        REQUIRE_FALSE(failed);
        const auto& cause = failed.error().cause;
        if (fault == Fault::Return) {
            CHECK(cause.code == provider.inner.code);
            CHECK(cause.operation == provider.inner.operation);
            CHECK(cause.details.numerical_generation == core::NumericalGeneration{99});
            CHECK(cause.details.resources == provider.inner.details.resources);
        } else {
            CHECK(cause.code == (fault == Fault::Throw ? core::FailureCode::FactorizationFailed
                                                       : core::FailureCode::ResourceExceeded));
        }
        CHECK(cause.details.model == core::ModelGeneration{1});
        CHECK(cause.details.trial == core::TrialId{11});
        CHECK(cause.details.parent == core::AcceptedStateId{71});
        CHECK(cause.details.last_accepted == core::AcceptedStateId{72});
        CHECK(failed.error().diagnostics.inverse_applications == 1);
        CHECK(session->first_failure() == failed.error());
        CHECK_FALSE(session->start_attempt());
        CHECK(session->first_failure() == failed.error());
    }
}

TEST_CASE("Eigen sessions reject reentry and unsupported precision without nested producer calls",
          "[negative][N4-06][N4-07]") {
    Provider<double> provider;
    auto session = ActionSession<double>::make(provider, ample);
    REQUIRE(session);
    const std::array high{1.0, 2.0};
    const auto input = WordView<double>::make(high, {}, core::Precision::Binary64);
    REQUIRE(input);
    std::array<double, 2> out{}, tail{}, error{};
    provider.hook = [&] {
        CHECK_FALSE(session->apply(provider, ActionKind::Forward, *input, out, tail, error));
    };
    const auto reentry = session->apply(provider, ActionKind::Forward, *input, out, tail, error);
    REQUIRE_FALSE(reentry);
    CHECK(reentry.error().cause.code == core::FailureCode::InvalidInput);
    CHECK(provider.calls == 1);
    CHECK(reentry.error().diagnostics.forward_applications == 1);
    provider.hook = {};
    provider.support = SupportedPrecision::Binary64Only;
    auto other = ActionSession<double>::make(provider, ample);
    REQUIRE(other);
    const std::array low{0x1p-80, 0.0};
    const auto two = WordView<double>::make(high, low, core::Precision::TwoWord);
    REQUIRE(two);
    const auto unsupported = other->apply(provider, ActionKind::Forward, *two, out, tail, error);
    REQUIRE_FALSE(unsupported);
    CHECK(unsupported.error().cause.code == core::FailureCode::PrecisionUnsupported);
    CHECK(unsupported.error().cause.operation == core::Operation::ExtractEigen);
    CHECK(provider.calls == 1);
}

TEST_CASE("Empty eigen sessions have no callbacks while moves do not duplicate quota",
          "[unit][negative][N4-07][N4-09]") {
    Provider<double> provider;
    provider.n = 0;
    auto session = ActionSession<double>::make(provider, {0, 0, 0, 0});
    REQUIRE(session);
    const auto input = WordView<double>::make({}, {}, core::Precision::TwoWord);
    REQUIRE(input);
    CHECK(session->apply(provider, ActionKind::Inverse, *input, {}, {}, {}));
    CHECK(provider.validations == 0);
    CHECK(provider.calls == 0);
    CHECK(session->diagnostics().inverse_applications == 0);
    auto moved = std::move(*session);
    CHECK_FALSE(session->active());
    CHECK_FALSE(session->apply(provider, ActionKind::Forward, *input, {}, {}, {}));
    CHECK(moved.apply(provider, ActionKind::Forward, *input, {}, {}, {}));
    CHECK_FALSE(moved.first_failure());
}

TEST_CASE("Eigen owned-payload planning checks full padded shapes and sums independently",
          "[unit][negative][N4-02][N4-09]") {
    const std::array extents{PayloadExtent{2, 3, 4, 8}, PayloadExtent{3, 2, 3, 16}};
    const auto bytes = checked_payload_bytes(extents);
    REQUIRE(bytes);
    CHECK(*bytes == 192);
    CHECK(checked_payload_bytes({}) == 0);
    const auto maximum = std::numeric_limits<core::Index>::max();
    for (const auto extent :
         {PayloadExtent{0, maximum, maximum, 8}, PayloadExtent{maximum, 0, maximum, 8}}) {
        const std::array arrays{extent};
        CHECK(checked_payload_bytes(arrays) == 0);
    }
    for (const auto extent : {PayloadExtent{-1, 0, 0, 8}, PayloadExtent{0, -1, 0, 8},
                              PayloadExtent{2, 1, 1, 8}, PayloadExtent{0, 0, 0, 0}}) {
        const std::array arrays{extent};
        const auto bad = checked_payload_bytes(arrays);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().code == core::FailureCode::InvalidInput);
    }
    for (const auto extent : {PayloadExtent{2, maximum, 2, 8}, PayloadExtent{1, maximum, 1, 8}}) {
        const std::array arrays{extent};
        const auto overflow = checked_payload_bytes(arrays);
        REQUIRE_FALSE(overflow);
        CHECK(overflow.error().code == core::FailureCode::ResourceExceeded);
    }
    const std::array sum{PayloadExtent{1, 1, 1, std::numeric_limits<std::size_t>::max()},
                         PayloadExtent{1, 1, 1, 1}};
    const auto overflow = checked_payload_bytes(sum);
    REQUIRE_FALSE(overflow);
    CHECK(overflow.error().code == core::FailureCode::ResourceExceeded);
}
static_assert(!std::is_copy_constructible_v<ActionSession<double>>);
static_assert(!std::is_move_assignable_v<ActionSession<double>>);

TEST_CASE("Eigen action entry rejects stale and malformed requests without numerical work",
          "[negative][N4-07][N4-08][N4-09]") {
    for (int site = 0; site < 8; ++site) {
        Provider<double> provider;
        auto session = ActionSession<double>::make(provider, ample);
        REQUIRE(session);
        const std::array high{1.0, 2.0};
        const auto input = WordView<double>::make(high, {}, core::Precision::Binary64);
        REQUIRE(input);
        std::array<double, 2> out{}, tail{}, error{};
        if (site == 0)
            provider.identity.original.cache.model = core::ModelGeneration{91};
        if (site == 1)
            provider.identity.original.trial = core::TrialId{92};
        if (site == 2)
            provider.identity.transformation = core::NumericalGeneration{93};
        if (site == 3)
            provider.identity.shift = {4, -5};
        if (site == 4)
            provider.n = 1;
        if (site == 7)
            provider.support = static_cast<SupportedPrecision>(255);
        const auto kind = site == 6 ? static_cast<ActionKind>(255) : ActionKind::Forward;
        const auto high_output =
            site == 5 ? std::span<double>{out}.first(1) : std::span<double>{out};
        const auto failed = session->apply(provider, kind, *input, high_output, tail, error);
        REQUIRE_FALSE(failed);
        const auto expected = site < 5    ? core::FailureCode::StaleGeneration
                              : site == 7 ? core::FailureCode::PrecisionUnsupported
                                          : core::FailureCode::InvalidInput;
        CHECK(failed.error().cause.code == expected);
        CHECK(provider.calls == 0);
        CHECK(failed.error().diagnostics.forward_applications == 0);
        CHECK(out == std::array{0.0, 0.0});
    }
}

TEST_CASE("Eigen metadata callbacks cannot reenter an active extraction session",
          "[negative][N4-07][N4-08]") {
    for (int site = 0; site < 3; ++site) {
        Provider<double> provider;
        auto session = ActionSession<double>::make(provider, ample);
        REQUIRE(session);
        provider.metadata_site = site;
        provider.metadata_hook = [&] { CHECK_FALSE(session->start_attempt()); };
        const std::array high{1.0, 2.0};
        const auto input = WordView<double>::make(high, {}, core::Precision::Binary64);
        REQUIRE(input);
        std::array<double, 2> out{}, tail{}, error{};
        const auto failed = session->apply(provider, ActionKind::Forward, *input, out, tail, error);
        REQUIRE_FALSE(failed);
        CHECK(failed.error().cause.code == core::FailureCode::InvalidInput);
        CHECK(failed.error().diagnostics.forward_applications == 0);
        CHECK(failed.error().diagnostics.extraction_attempts == 0);
        CHECK(provider.calls == 0);
        CHECK(out == std::array{0.0, 0.0});
        CHECK(session->first_failure() == failed.error());
    }
}

TEST_CASE("Eigen original coefficient representation is validated independently of action input",
          "[negative][N4-06][N4-07]") {
    for (const core::Index n : {0, 2}) {
        Provider<double> provider;
        provider.n = n;
        provider.identity.original.cache.representation = static_cast<core::Precision>(255);
        const auto unknown = ActionSession<double>::make(provider, ample);
        REQUIRE_FALSE(unknown);
        CHECK(unknown.error().cause.code == core::FailureCode::PrecisionUnsupported);
        CHECK(unknown.error().cause.operation == core::Operation::ExtractEigen);
        CHECK(provider.validations == 0);
        CHECK(provider.calls == 0);
    }
    Provider<double> provider;
    provider.identity.original.cache.representation = core::Precision::TwoWord;
    const auto known = ActionSession<double>::make(provider, ample);
    REQUIRE(known);
    CHECK(provider.validations == 1);
    CHECK(provider.calls == 0);
}

TEST_CASE("Eigen capability metadata state changes reject before charging an action",
          "[negative][N4-07][N4-08][N4-09]") {
    Provider<double> provider;
    auto session = ActionSession<double>::make(provider, ample);
    REQUIRE(session);
    provider.metadata_site = 2;
    provider.metadata_hook = [&] {
        provider.identity.transformation = core::NumericalGeneration{94};
    };
    const std::array high{1.0, 2.0};
    const auto input = WordView<double>::make(high, {}, core::Precision::Binary64);
    REQUIRE(input);
    std::array<double, 2> out{}, tail{}, error{};
    const auto failed = session->apply(provider, ActionKind::Forward, *input, out, tail, error);
    REQUIRE_FALSE(failed);
    CHECK(failed.error().cause.code == core::FailureCode::StaleGeneration);
    CHECK(failed.error().cause.operation == core::Operation::ExtractEigen);
    CHECK(failed.error().diagnostics.forward_applications == 0);
    CHECK(provider.calls == 0);
    CHECK(out == std::array{0.0, 0.0});
    CHECK(session->first_failure() == failed.error());
}
