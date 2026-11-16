#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <utility>
#include <uvlas/core/checked_dimensions.hpp>
#include <uvlas/numerics/eigen_action.hpp>

namespace uvlas::numerics::eigen {
namespace {
constexpr auto operation = core::Operation::ExtractEigen;
using core::FailureCode;

template <StorageScalar Scalar>
bool finite(Scalar value) noexcept {
    if constexpr (std::same_as<Scalar, double>)
        return std::isfinite(value);
    else
        return std::isfinite(value.real()) && std::isfinite(value.imag());
}
template <class A, class B>
bool overlap(std::span<A> a, std::span<B> b) noexcept {
    if (a.empty() || b.empty())
        return false;
    // std::less supplies a total pointer order without subtracting unrelated pointers.
    const std::less<const void*> before;
    return before(a.data(), b.data() + b.size()) && before(b.data(), a.data() + a.size());
}
bool known_kind(ActionKind kind) noexcept {
    return kind == ActionKind::Forward || kind == ActionKind::Adjoint ||
           kind == ActionKind::Inverse;
}
} // namespace

template <StorageScalar Scalar>
ActionSession<Scalar>::ActionSession(core::Index order, ActionStamp stamp, ActionLimits limits,
                                     StateDiagnosticContext state) noexcept
    : order_(order), stamp_(stamp), limits_(limits), state_(state) {
    diagnostics_.extraction_attempts = 0;
    diagnostics_.forward_applications = 0;
    diagnostics_.adjoint_applications = 0;
    diagnostics_.inverse_applications = 0;
}

template <StorageScalar Scalar>
Result<ActionSession<Scalar>> ActionSession<Scalar>::make(const ActionProvider<Scalar>& provider,
                                                          ActionLimits limits,
                                                          StateDiagnosticContext state) noexcept {
    const auto order = provider.order();
    const auto stamp = provider.stamp();
    ActionSession session{order, stamp, limits, state};
    if (order < 0 || !std::isfinite(stamp.original.time) || !finite(stamp.shift)) {
        auto result = session.fail(order < 0 ? FailureCode::InvalidInput : FailureCode::Nonfinite);
        return std::unexpected(result.error());
    }
    if (stamp.original.cache.representation != core::Precision::Binary64 &&
        stamp.original.cache.representation != core::Precision::TwoWord) {
        auto result = session.fail(FailureCode::PrecisionUnsupported);
        return std::unexpected(result.error());
    }
    if (order != 0) {
        auto valid = core::guard_entry(operation, FailureCode::FactorizationFailed, {},
                                       [&] { return provider.validate_frozen_state(); });
        if (!valid) {
            auto result = session.remember(valid.error());
            return std::unexpected(result.error());
        }
        if (provider.order() != order || provider.stamp() != stamp) {
            auto result = session.fail(FailureCode::StaleGeneration);
            return std::unexpected(result.error());
        }
    }
    return session;
}

template <StorageScalar Scalar>
ActionSession<Scalar>::ActionSession(ActionSession&& other) noexcept
    : order_(other.order_),
      stamp_(other.stamp_),
      limits_(other.limits_),
      state_(other.state_),
      diagnostics_(other.diagnostics_),
      first_(other.first_),
      active_(std::exchange(other.active_, false)) {}

template <StorageScalar Scalar>
Result<void> ActionSession<Scalar>::remember(core::Failure cause) noexcept {
    if (!first_)
        first_ = complete_failure(cause, stamp_.original, state_, diagnostics_);
    return std::unexpected(*first_);
}
template <StorageScalar Scalar>
Result<void> ActionSession<Scalar>::fail(FailureCode code) noexcept {
    return remember(core::Failure{code, operation});
}
template <StorageScalar Scalar>
Result<void> ActionSession<Scalar>::start_attempt() noexcept {
    if (first_)
        return std::unexpected(*first_);
    if (!active_ || busy_)
        return fail(FailureCode::InvalidInput);
    auto& count = diagnostics_.extraction_attempts;
    // These four counts are always known: the session owns every start and retry.
    if (!count || *count == std::numeric_limits<std::uint64_t>::max())
        return fail(FailureCode::ResourceExceeded);
    if (*count >= limits_.attempts) {
        core::FailureDetails details;
        details.resources = core::ResourceUsage{core::ResourceKind::Iterations, *count + 1, *count};
        return remember(core::Failure{FailureCode::ResourceExceeded, operation, details});
    }
    ++*count;
    return {};
}

template <StorageScalar Scalar>
Result<void> ActionSession<Scalar>::apply(const ActionProvider<Scalar>& provider, ActionKind kind,
                                          WordView<Scalar> input, std::span<Scalar> high,
                                          std::span<Scalar> low,
                                          std::span<double> error_upper) noexcept {
    if (first_)
        return std::unexpected(*first_);
    if (!active_ || busy_)
        return fail(FailureCode::InvalidInput);
    busy_ = true;
    struct ResetBusy {
        bool& busy;
        ~ResetBusy() {
            busy = false;
        }
    };
    const ResetBusy reset{busy_};
    // Even metadata getters belong to user code and can attempt synchronous reentry.
    const auto check_stamp = [&]() -> Result<void> {
        const auto current_order = provider.order();
        if (first_)
            return std::unexpected(*first_);
        if (current_order != order_)
            return fail(FailureCode::StaleGeneration);
        const auto current_stamp = provider.stamp();
        if (first_)
            return std::unexpected(*first_);
        if (current_stamp != stamp_)
            return fail(FailureCode::StaleGeneration);
        return {};
    };
    if (!known_kind(kind) || std::cmp_not_equal(input.size(), order_) ||
        high.size() != input.size() || low.size() != input.size() ||
        error_upper.size() != input.size())
        return fail(FailureCode::InvalidInput);
    // A previously validated borrowed view can outlive caller writes; validate again here.
    const auto checked = WordView<Scalar>::make(input.high(), input.low(), input.precision());
    if (!checked)
        return remember(checked.error());
    if (overlap(high, low) || overlap(high, error_upper) || overlap(low, error_upper) ||
        overlap(input.high(), high) || overlap(input.high(), low) ||
        overlap(input.high(), error_upper) || overlap(input.low(), high) ||
        overlap(input.low(), low) || overlap(input.low(), error_upper))
        return fail(FailureCode::InvalidInput);
    auto before = check_stamp();
    if (!before)
        return before;
    if (order_ == 0)
        return {};
    const auto capability = provider.supported_precision(kind);
    if (first_)
        return std::unexpected(*first_);
    const auto supported = require_precision(input.precision(), capability);
    if (!supported) {
        // This is our representation preflight, not a failure reported by an operator.
        return remember(
            core::Failure{supported.error().code, operation, supported.error().details});
    }
    // Capability metadata is user code too; detect its state change before any action charge.
    auto before_callback = check_stamp();
    if (!before_callback)
        return before_callback;
    auto* count = &diagnostics_.forward_applications;
    std::uint64_t limit = limits_.forward;
    if (kind == ActionKind::Adjoint) {
        count = &diagnostics_.adjoint_applications;
        limit = limits_.adjoint;
    } else if (kind == ActionKind::Inverse) {
        count = &diagnostics_.inverse_applications;
        limit = limits_.inverse;
    }
    if (!*count || **count == std::numeric_limits<std::uint64_t>::max())
        return fail(FailureCode::ResourceExceeded);
    if (**count >= limit) {
        core::FailureDetails details;
        details.resources =
            core::ResourceUsage{core::ResourceKind::OperatorApplications, **count + 1, **count};
        return remember(core::Failure{FailureCode::ResourceExceeded, operation, details});
    }
    ++**count;
    // Poison all scratch before entry: an omitted component cannot reuse old finite data.
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    const auto sentinel = [nan] {
        if constexpr (std::same_as<Scalar, double>)
            return nan;
        else
            return Scalar{nan, nan};
    }();
    std::fill(high.begin(), high.end(), sentinel);
    std::fill(low.begin(), low.end(), sentinel);
    std::fill(error_upper.begin(), error_upper.end(), nan);
    auto result = core::guard_entry(operation, FailureCode::FactorizationFailed, {}, [&] {
        return provider.apply(kind, input, high, low, error_upper);
    });
    // A reentrant failure occurs inside the producer before its outer return.
    if (first_)
        return std::unexpected(*first_);
    if (!result)
        return remember(result.error());
    auto after = check_stamp();
    if (!after)
        return after;
    for (std::size_t i = 0; i < high.size(); ++i) {
        if (!finite(high[i]) || !finite(low[i]) || !std::isfinite(error_upper[i]))
            return fail(FailureCode::Nonfinite);
        if (error_upper[i] < 0)
            return fail(FailureCode::InvalidInput);
    }
    return {};
}

template <StorageScalar Scalar>
AttemptDiagnostics ActionSession<Scalar>::diagnostics() const noexcept {
    return diagnostics_;
}
template <StorageScalar Scalar>
std::optional<EigenFailure> ActionSession<Scalar>::first_failure() const noexcept {
    return first_;
}
template <StorageScalar Scalar>
bool ActionSession<Scalar>::active() const noexcept {
    return active_;
}

core::Result<std::uint64_t> checked_payload_bytes(std::span<const PayloadExtent> arrays) noexcept {
    std::uint64_t total = 0;
    for (const auto& extent : arrays) {
        if (extent.rows < 0 || extent.columns < 0 || extent.leading_dimension < extent.rows ||
            extent.leading_dimension < 0 || extent.element_bytes == 0)
            return std::unexpected(core::Failure{FailureCode::InvalidInput, operation});
        const std::array factors{extent.rows == 0 ? core::Index{0} : extent.leading_dimension,
                                 extent.columns};
        auto count = core::checked_extent_product(factors, operation);
        if (!count)
            return std::unexpected(count.error());
        auto bytes = core::checked_byte_count(*count, extent.element_bytes, operation);
        if (!bytes)
            return std::unexpected(bytes.error());
        if (std::cmp_greater(*bytes, std::numeric_limits<std::uint64_t>::max() - total))
            return std::unexpected(core::Failure{FailureCode::ResourceExceeded, operation});
        total += static_cast<std::uint64_t>(*bytes);
    }
    return total;
}

template class ActionSession<double>;
template class ActionSession<std::complex<double>>;
} // namespace uvlas::numerics::eigen
