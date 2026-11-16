#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <functional>
#include <limits>
#include <new>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>
#include <uvlas/numerics/linear_solve.hpp>
#include <variant>

namespace uvlas::numerics::eigen {

/// Operator, layout, scale, policy, model, trial and physical-time identity.
using Context = linear::Context;

/// Optional reporting identities; these do not validate an AcceptedState's ancestry.
struct StateDiagnosticContext {
    std::optional<core::AcceptedStateId> parent;        ///< Explicitly supplied accepted parent.
    std::optional<core::AcceptedStateId> last_accepted; ///< Last state that remains available.
};

/// Owned attempt counters. Missing means unavailable or not performed, never a measured zero.
struct AttemptDiagnostics {
    std::optional<std::uint64_t> requested_count;      ///< Number of requested eigenvalues.
    std::optional<std::uint64_t> converged_count;      ///< Known backend convergence count.
    std::optional<std::uint64_t> backend_iterations;   ///< Counter reported by the actual backend.
    std::optional<std::uint64_t> extraction_attempts;  ///< Started extraction attempts.
    std::optional<std::uint64_t> forward_applications; ///< Actual original-operator calls.
    std::optional<std::uint64_t> adjoint_applications; ///< Actual adjoint-operator calls.
    std::optional<std::uint64_t> inverse_applications; ///< Actual shifted-inverse calls.
    /// Compare both the availability and value of every counter.
    bool operator==(const AttemptDiagnostics&) const = default;
};

/// Allocation-free failure preserving its original cause and owned attempt diagnostics.
struct EigenFailure {
    core::Failure cause; ///< Existing failure code, operation and known context are preserved.
    AttemptDiagnostics diagnostics; ///< Counters describing this attempt, not a later retry.
    /// Compare all owned context and counters.
    bool operator==(const EigenFailure&) const = default;
};

/// An owned eigen result or an allocation-free failure; success does not imply physical coverage.
template <class T>
using Result = std::expected<T, EigenFailure>;

/// One explicit error budget and the quantity/domain to which it applies.
struct MetricRequirement {
    double budget; ///< Nonnegative finite error budget, in the metric's documented scale.
    core::ObservationId observation; ///< Quantity evaluated against this budget.
    core::ScopeId scope;             ///< Exact domain covered by the requirement.
};
/// This metric is not requested; it must not be reported as an achieved zero.
struct NotRequested {};
/// Optional extra componentwise criterion, independent of the mandatory norm residual.
using BackwardRequirement = std::variant<NotRequested, MetricRequirement>;

/// Supply only missing cause fields/counters from the explicit outer context.
/// Existing callback diagnostics, operation, resources and identifiers are never replaced.
[[nodiscard]] EigenFailure complete_failure(core::Failure cause, const Context& context,
                                            StateDiagnosticContext state,
                                            AttemptDiagnostics diagnostics = {}) noexcept;
/// Complete an existing eigen failure without replacing its known counters.
[[nodiscard]] EigenFailure complete_failure(EigenFailure failure, const Context& context,
                                            StateDiagnosticContext state,
                                            const AttemptDiagnostics& fallback) noexcept;

/// Convert allocation/backend exceptions to an owned failure at a synchronous entry boundary.
/// A returned failure retains its operation and known diagnostics; context only fills gaps.
/// The callable returns eigen::Result<T>; adapt a core::Result with complete_failure first.
template <class T, class Function>
    requires std::same_as<std::invoke_result_t<Function>, Result<T>>
[[nodiscard]] Result<T> guard_entry(const Context& context, StateDiagnosticContext state,
                                    const AttemptDiagnostics& diagnostics,
                                    Function&& function) noexcept {
    try {
        auto result = std::invoke(std::forward<Function>(function));
        if (!result) {
            return std::unexpected(
                complete_failure(std::move(result.error()), context, state, diagnostics));
        }
        return result;
    } catch (const std::bad_alloc&) {
        return std::unexpected(complete_failure(
            core::Failure{core::FailureCode::ResourceExceeded, core::Operation::ExtractEigen},
            context, state, diagnostics));
    } catch (...) {
        return std::unexpected(complete_failure(
            core::Failure{core::FailureCode::FactorizationFailed, core::Operation::ExtractEigen},
            context, state, diagnostics));
    }
}

/// Validated, synchronous read-only numerical words; the caller keeps both spans alive.
/// No ownership or allocation is introduced, and a two-word low plane is never discarded.
template <StorageScalar Scalar>
class WordView {
public:
    /// Validate the representation, equal word lengths and all finite real/imaginary parts.
    /// Binary64 requires an empty low span. Unknown representations are unsupported.
    [[nodiscard]] static core::Result<WordView> make(std::span<const Scalar> high,
                                                     std::span<const Scalar> low,
                                                     core::Precision representation) noexcept {
        const auto fail = [](core::FailureCode code) {
            return std::unexpected(core::Failure{code, core::Operation::ExtractEigen});
        };
        if (representation != core::Precision::Binary64 &&
            representation != core::Precision::TwoWord) {
            return fail(core::FailureCode::PrecisionUnsupported);
        }
        if (std::cmp_greater(high.size(), std::numeric_limits<core::Index>::max()) ||
            (representation == core::Precision::Binary64 && !low.empty()) ||
            (representation == core::Precision::TwoWord && high.size() != low.size())) {
            return fail(core::FailureCode::InvalidInput);
        }
        const auto finite = [](Scalar value) {
            if constexpr (std::same_as<Scalar, double>) {
                return std::isfinite(value);
            } else {
                return std::isfinite(value.real()) && std::isfinite(value.imag());
            }
        };
        for (const auto value : high) {
            if (!finite(value)) {
                return fail(core::FailureCode::Nonfinite);
            }
        }
        for (const auto value : low) {
            if (!finite(value)) {
                return fail(core::FailureCode::Nonfinite);
            }
        }
        return WordView{high, low, representation};
    }
    /// Read the borrowed leading plane during the synchronous call.
    [[nodiscard]] std::span<const Scalar> high() const noexcept {
        return high_;
    }
    /// Read the borrowed trailing plane; empty only for binary64 or an empty vector.
    [[nodiscard]] std::span<const Scalar> low() const noexcept {
        return low_;
    }
    /// Arithmetic representation of these exact words.
    [[nodiscard]] core::Precision precision() const noexcept {
        return precision_;
    }
    /// Logical number of entries; the view has no padding or leading dimension.
    [[nodiscard]] std::size_t size() const noexcept {
        return high_.size();
    }

private:
    WordView(std::span<const Scalar> high, std::span<const Scalar> low,
             core::Precision representation) noexcept
        : high_(high), low_(low), precision_(representation) {}
    std::span<const Scalar> high_;
    std::span<const Scalar> low_;
    core::Precision precision_;
};

} // namespace uvlas::numerics::eigen
