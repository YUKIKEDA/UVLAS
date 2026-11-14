#pragma once

#include <cstdint>
#include <expected>
#include <functional>
#include <new>
#include <optional>
#include <type_traits>
#include <utility>
#include <uvlas/core/identifier.hpp>
#include <variant>

namespace uvlas::core {

/// Reason an operation failed. Inspect named values; do not persist their ordinals.
enum class FailureCode : std::uint8_t {
    InvalidInput,            ///< Input violates the operation's contract.
    Unsupported,             ///< Requested capability is unsupported.
    GeometryInvalid,         ///< Geometry is invalid for the requested operation.
    RankUnresolved,          ///< Required rank could not be established.
    ChartSingular,           ///< Coordinate chart is singular.
    Nonfinite,               ///< A required finite value is NaN or infinity.
    FactorizationFailed,     ///< Required factorization failed.
    LinearUnconverged,       ///< Linear solve did not meet its convergence criteria.
    NonlinearUnconverged,    ///< Nonlinear solve did not meet its convergence criteria.
    EigenCoverageUnresolved, ///< Required eigenvalue coverage could not be certified.
    PrecisionUnsupported,    ///< Requested numerical representation is unsupported.
    PrecisionUnresolved,     ///< Required accuracy could not be certified.
    ResourceExceeded,        ///< Resource limit or allocation failure prevented completion.
    Cancelled,               ///< Operation was cancelled.
    StaleGeneration,         ///< Supplied identity refers to an incompatible generation.
    CheckpointCorrupt,       ///< Checkpoint integrity validation failed.
    CheckpointIncompatible   ///< Checkpoint is incompatible with the receiving context.
};

/// Entry operation to which a failure belongs, including certificate validation.
enum class Operation : std::uint8_t {
    PrepareModel,         ///< Model preparation.
    EvaluateStructure,    ///< Structural evaluation.
    EvaluateAero,         ///< Aerodynamic evaluation.
    Initialize,           ///< Initial-state construction.
    EvaluateStep,         ///< Time-step evaluation.
    SolveLinear,          ///< Linear-system solve.
    AdvanceTrial,         ///< Trial-state advancement.
    Commit,               ///< Acceptance of a trial state that passed the required checks.
    LinearizeEquilibrium, ///< Equilibrium linearization.
    SearchFlutter,        ///< Flutter search.
    Recover,              ///< Checkpoint recovery.
    CertifyError          ///< Error-certificate input validation.
};

/// Unit shared by the required and observed counts in ResourceUsage.
enum class ResourceKind : std::uint8_t {
    Bytes,                ///< Memory measured in bytes.
    Iterations,           ///< Count of iterations.
    OperatorApplications, ///< Count of operator applications.
    Nanoseconds           ///< Duration measured in nanoseconds.
};

/// Resource diagnostic with explicit units; omit the entire value if unknown.
struct ResourceUsage {
    /// Unit for both counts; counts are not interchangeable across resource kinds.
    ResourceKind kind;
    /// Required amount determined by the reporting operation.
    std::uint64_t required;
    /// Amount measured by the reporting operation.
    std::uint64_t observed;
    /// Store supplied counts in the same unit; does not estimate resource demand.
    constexpr ResourceUsage(ResourceKind unit, std::uint64_t required_count,
                            std::uint64_t observed_count) noexcept
        : kind(unit), required(required_count), observed(observed_count) {}
    /// Compare unit and both counts exactly.
    bool operator==(const ResourceUsage&) const = default;
};

/// Owned identity of the object implicated in a failure; the alternative keeps its domain.
using ProblemId =
    std::variant<DofId, PanelId, WakeVertexId, WakeFaceId, LayerId, ObservationId, ScopeId>;

/// Owned, allocation-free diagnostic context. A missing field is unknown, never ID zero.
/// The reporting operation supplies these identities; this value does not validate them.
struct FailureDetails {
    /// Model generation involved in the operation, when available.
    std::optional<ModelGeneration> model;
    /// Accepted parent from which the operation started, when available.
    std::optional<AcceptedStateId> parent;
    /// Trial involved in the operation, when available.
    std::optional<TrialId> trial;
    /// Numerical generation used by the failed evaluation, when available.
    std::optional<NumericalGeneration> numerical_generation;
    /// Component containing the problematic object, when identifiable.
    std::optional<ComponentId> component;
    /// Problematic object with its identity domain preserved, when identifiable.
    std::optional<ProblemId> problem;
    /// Known resource demand and measurement in one unit, when available.
    std::optional<ResourceUsage> resources;
    /// Last accepted state that remains available after failure, when available.
    std::optional<AcceptedStateId> last_accepted;
    /// Compare all fields, including whether optional information is present.
    bool operator==(const FailureDetails&) const = default;
};

/// Owned diagnostic safe to construct and copy while handling allocation failure.
/// No field borrows exception text or temporary calculation storage.
struct Failure {
    /// Failure category, distinct from the operation where it occurred.
    FailureCode code;
    /// Entry operation responsible for this diagnostic.
    Operation operation;
    /// Available identities and resource measurements, stored by value.
    FailureDetails details;

    /// Copy supplied context without allocation; absent details remain unknown.
    Failure(FailureCode failure_code, Operation process, FailureDetails context = {}) noexcept
        : code(failure_code), operation(process), details(context) {}
    /// Compare category, operation and all diagnostic context exactly.
    bool operator==(const Failure&) const = default;
};

/// An owned result or an explicit Failure; an empty success must not hide failure.
/// @tparam T Success value type, or void for a successful operation without a value.
template <class T>
using Result = std::expected<T, Failure>;

namespace detail {
template <class T>
inline constexpr bool is_result = false;
template <class T>
inline constexpr bool is_result<Result<T>> = true;
} // namespace detail

// No allocating diagnostic is constructed while translating an allocation failure.
static_assert(std::is_nothrow_copy_constructible_v<Failure>);

/// Invoke an action once and return exceptions as Failure values without throwing.
/// @param operation Entry operation used only for translated exceptions.
/// @param exception_code Caller-selected category for exceptions other than bad_alloc.
/// @param details Context copied into translated failures; borrowed only during this call.
/// @param action Callable returning Result<T>, including Result<void>.
/// @return The action's result unchanged, or ResourceExceeded for bad_alloc, or
/// exception_code for any other exception. Existing returned failures keep their context.
/// @note Diagnostic construction does not allocate. The action may allocate or throw.
/// This guard does not validate inputs or undo changes made by action. An action
/// must leave any previously accepted state unchanged on failure; compute changes
/// in a separate candidate and publish them only after success.
template <class F>
    requires detail::is_result<std::invoke_result_t<F>>
[[nodiscard]] auto guard_entry(Operation operation, FailureCode exception_code,
                               const FailureDetails& details, F&& action) noexcept
    -> std::invoke_result_t<F> {
    try {
        return std::invoke(std::forward<F>(action));
    } catch (const std::bad_alloc&) {
        // N4-08: allocation failure has a fixed category regardless of the entry's fallback.
        // Copy only pre-existing value context; allocating exception text could fail again.
        return std::unexpected(Failure{FailureCode::ResourceExceeded, operation, details});
    } catch (...) {
        // The entry chooses the operation-specific code; existing expected failures
        // returned by action are preserved, not reclassified by this guard.
        return std::unexpected(Failure{exception_code, operation, details});
    }
}

} // namespace uvlas::core
