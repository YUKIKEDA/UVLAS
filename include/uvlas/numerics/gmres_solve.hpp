#pragma once

#include <cstdint>
#include <expected>
#include <memory>
#include <span>
#include <utility>
#include <uvlas/numerics/linear_solve.hpp>

namespace uvlas::numerics::linear {

/// Why one GMRES invocation stopped; the common Failure retains its precise category.
enum class GmresStopReason : std::uint8_t {
    NotStarted,          ///< Input/preparation failed before numerical action.
    Qualified,           ///< Every original-system and necessary consumer check passed.
    OperatorLimit,       ///< Original-operator action budget prevented further work.
    KrylovBreakdown,     ///< The available Krylov direction cannot qualify the original system.
    Stagnated,           ///< A complete fixed cycle left both retained words bitwise unchanged.
    ActionFailure,       ///< An operator/preconditioner callback or fixed-state check failed.
    ArithmeticFailure,   ///< Search-space arithmetic could not produce a finite next candidate.
    QualificationFailure ///< Required retained/projected/consumer qualification failed.
};

/// Value-owned counts for the complete call, summed across all RHS columns.
struct GmresDiagnostics {
    std::uint64_t operator_actions = 0; ///< All started forward and absolute-lower actions.
    std::uint64_t forward_actions = 0;  ///< Started A*word callbacks, including qualification.
    std::uint64_t absolute_actions = 0; ///< Started |A|*weights lower-bound callbacks.
    std::uint64_t preconditioner_actions = 0; ///< Started fixed left-preconditioner callbacks.
    std::uint64_t total_action_callbacks = 0; ///< Forward + absolute + preconditioner callbacks.
    std::uint64_t arnoldi_steps = 0; ///< Completed Arnoldi steps, never LU refinement count.
    std::uint64_t restarts = 0;      ///< Started cycles after the initial Krylov cycle.
    GmresStopReason stop = GmresStopReason::NotStarted; ///< Actual termination reason.
};

namespace detail {
class GmresFactory;
struct GmresData;
} // namespace detail

/// Own the original common Failure and the counters of this failed invocation.
/// Later preparation/workspace reuse cannot alter either value.
class GmresFailure {
public:
    /// Preserve ordinary value semantics for allocation-free failure reporting.
    GmresFailure(const GmresFailure&) noexcept = default;
    /// Copy both owned diagnostic values.
    GmresFailure& operator=(const GmresFailure&) noexcept = default;
    /// Move both owned diagnostic values.
    GmresFailure(GmresFailure&&) noexcept = default;
    /// Replace both owned diagnostic values.
    GmresFailure& operator=(GmresFailure&&) noexcept = default;
    /// Release value-only diagnostics.
    ~GmresFailure() = default;
    /// Read the common category, operation and inner context without relabelling it.
    [[nodiscard]] const core::Failure& failure() const& noexcept {
        return failure_;
    }
    /// A temporary cannot provide a diagnostic borrow.
    [[nodiscard]] const core::Failure& failure() const&& = delete;
    /// Copy the exact counters, including calls that returned or threw a failure.
    [[nodiscard]] GmresDiagnostics diagnostics() const noexcept {
        return diagnostics_;
    }

private:
    friend class detail::GmresFactory;
    GmresFailure(core::Failure failure, GmresDiagnostics diagnostics) noexcept
        : failure_(failure), diagnostics_(diagnostics) {}
    core::Failure failure_;
    GmresDiagnostics diagnostics_;
};

/// GMRES keeps its own diagnostics on failure instead of replacing core::Failure fields.
template <class Value>
using GmresResult = std::expected<Value, GmresFailure>;

/// Matrix-free access to one immutable finite real linear operator in physical coordinates.
/// Implementations own or otherwise freeze their coefficient state and workspace lifetime.
/// Numerical callbacks are synchronous and allocate nothing; Python callbacks are unsupported.
/// Outputs are borrowed only for the call. Every bound refers to the same original A.
class GmresOperator {
public:
    /// Release through the interface.
    virtual ~GmresOperator() = default;
    /// Nonnegative square-system order; metadata access does not consume an action.
    [[nodiscard]] virtual core::Index order() const noexcept = 0;
    /// Complete model/chart/matrix/scale/policy identity, unchanged throughout a solve.
    /// Any change of the represented action must change its adopted identity.
    [[nodiscard]] virtual Context context() const noexcept = 0;
    /// Validate frozen finite coefficients and producer-owned state before numerical work.
    /// No action is charged for validation; empty systems do not invoke this callback.
    [[nodiscard]] virtual core::Result<void> validate_frozen_state() const = 0;
    /// Evaluate one binary64 word plane: abs(A*x-estimate_hi-estimate_lo)_i <= error_upper_i.
    /// Signed low words preserve exact-product remainders; providers may return zero with a valid
    /// error bound. All spans have order() entries. Return unavailable bounds as
    /// PrecisionUnresolved. Error upper bounds are finite/nonnegative. One started call consumes
    /// one operator action. GMRES invokes this separately for nonzero high/low planes and never
    /// drops a low word.
    [[nodiscard]] virtual core::Result<void> apply_word(std::span<const double> x,
                                                        std::span<double> estimate_hi,
                                                        std::span<double> estimate_lo,
                                                        std::span<double> error_upper) const = 0;
    /// Return finite 0<=lower_i<=sum_j abs(A_ij)*weights_j for nonnegative finite weights.
    /// GMRES supplies component lower bounds of abs(actual x_hi+x_lo), so monotonicity
    /// gives a valid standard backward-error denominator lower bound.
    /// A JVP alone is insufficient to implement this callback. Unavailable bounds reject.
    /// One started call consumes one operator action, including retained/projection checks.
    [[nodiscard]] virtual core::Result<void>
    absolute_action_lower(std::span<const double> weights, std::span<double> lower) const = 0;
};

/// Bind a stationary linear preconditioner to its target and actual frozen factors/policies.
struct FixedPreconditionerIdentity {
    Context target; ///< Original system, evaluation time and physical layouts being solved.
    core::CacheKey
        factors; ///< Full factor identity, including numeric, precision and backend state.
    /// Compare every identity field, not the callback address or only its dimensions.
    bool operator==(const FixedPreconditionerIdentity&) const = default;
};

/// One stationary linear left-preconditioning action, fixed for a complete solve.
/// Implementers must freeze mathematical state; const alone is not that guarantee.
/// Scratch may change internally only with producer-enforced exclusive access.
/// The adapter owner separately accounts all retained/backend/workspace resources.
class FixedPreconditioner {
public:
    /// Release through the interface.
    virtual ~FixedPreconditioner() = default;
    /// Input/output dimension of the residual-to-unknown correction map.
    [[nodiscard]] virtual core::Index order() const noexcept = 0;
    /// Full target and factor/policy identities, checked before and after every action.
    [[nodiscard]] virtual FixedPreconditionerIdentity identity() const noexcept = 0;
    /// Check frozen finite state before numerical work; not called for empty systems.
    [[nodiscard]] virtual core::Result<void> validate_frozen_state() const = 0;
    /// Apply the same linear map to an original-residual vector, returning physical unknowns.
    /// Spans have order() entries; they do not alias. No allocation is permitted.
    /// Each started call is counted separately from the original-operator action quota.
    [[nodiscard]] virtual core::Result<void> apply(std::span<const double> residual,
                                                   std::span<double> correction) const = 0;
};

/// Explicit capacities, operator quota and independently scoped memory admission.
struct GmresLimits {
    /// Set size capacities and memory mode; defaults 64/512 are resources, not tolerances.
    GmresLimits(core::Index order, core::Index rhs_columns, BackendMemoryLimit memory,
                core::Index restart_width = 64, std::uint64_t operator_actions = 512) noexcept
        : maximum_order(order),
          maximum_rhs_columns(rhs_columns),
          restart(restart_width),
          maximum_operator_actions(operator_actions),
          backend_memory(memory) {}
    core::Index maximum_order;       ///< Nonnegative prepared order capacity.
    core::Index maximum_rhs_columns; ///< Nonnegative RHS capacity.
    core::Index restart; ///< Positive Krylov width; preparation clips storage to the system order.
    std::uint64_t maximum_operator_actions; ///< Shared across all RHS and qualification/projection.
    BackendMemoryLimit backend_memory; ///< GMRES-owned scope; opaque adapter memory is separate.
};

/// Own a common qualified solution and distinct GMRES iteration diagnostics.
class GmresSolution {
public:
    /// Transfer owning values; invalidate former borrows.
    GmresSolution(GmresSolution&&) noexcept = default;
    /// Replace owning values; invalidate former destination borrows.
    GmresSolution& operator=(GmresSolution&&) noexcept = default;
    /// Implicit copying is forbidden.
    GmresSolution(const GmresSolution&) = delete;
    /// Copy assignment is forbidden.
    GmresSolution& operator=(const GmresSolution&) = delete;
    /// Release owned values and certificate independently of prepared scratch.
    ~GmresSolution() = default;
    /// Read the exact retained/projected owning result and common certificate.
    [[nodiscard]] const Solution<double>& solution() const& noexcept {
        return solution_;
    }
    /// A temporary cannot provide a solution borrow.
    [[nodiscard]] const Solution<double>& solution() const&& = delete;
    /// Copy per-call GMRES counters and stop reason; common refinements remains LU-specific.
    [[nodiscard]] GmresDiagnostics diagnostics() const noexcept {
        return diagnostics_;
    }

private:
    friend class detail::GmresFactory;
    GmresSolution(Solution<double> solution, GmresDiagnostics diagnostics) noexcept
        : solution_(std::move(solution)), diagnostics_(diagnostics) {}
    Solution<double> solution_;
    GmresDiagnostics diagnostics_;
};

/// Own reusable Krylov arrays and scale/unit/identity snapshots, never callback borrows.
class GmresPrepared {
public:
    /// Validate immutable adapters and preallocate solver-owned payload before iteration.
    /// A/P objects are borrowed only during this call; solve supplies live adapters again.
    [[nodiscard]] static GmresResult<GmresPrepared>
    prepare(const GmresOperator& original, const FixedPreconditioner& preconditioner,
            const Scales& scales, const GmresLimits& limits,
            core::ResourcePlan& workspace_plan) noexcept;
    /// Transfer preparation ownership; the source becomes inactive.
    GmresPrepared(GmresPrepared&&) noexcept;
    /// Replace preparation ownership.
    GmresPrepared& operator=(GmresPrepared&&) noexcept;
    /// Copying is forbidden.
    GmresPrepared(const GmresPrepared&) = delete;
    /// Copy assignment is forbidden.
    GmresPrepared& operator=(const GmresPrepared&) = delete;
    /// Release scratch independently of every previously returned result.
    ~GmresPrepared();
    /// Whether this object still owns its preparation.
    [[nodiscard]] bool active() const noexcept {
        return data_ != nullptr;
    }

private:
    friend GmresResult<GmresSolution> solve_linear(const GmresOperator&, const FixedPreconditioner&,
                                                   const Rhs<double>&, const Scales&,
                                                   const Precision<double>&, const Accuracy&,
                                                   const GmresLimits&, GmresPrepared&,
                                                   core::ResourcePlan&) noexcept;
    explicit GmresPrepared(std::unique_ptr<detail::GmresData> data) noexcept;
    std::unique_ptr<detail::GmresData> data_;
};

/// Restarted real GMRES with fixed left preconditioning and zero initial correction.
/// The original operator is used for both retained and requested projection qualification.
/// Every started A-word/absolute action is charged before invocation; high/low planes,
/// restarts and final checks share one quota across all columns. No action is issued after it.
/// P calls and Arnoldi/restarts have separate bounded counters. Input parents never change.
/// Unsupported/failed callbacks preserve their common Failure; unmet original criteria return
/// LinearUnconverged, whereas an unmet explicit output projection is PrecisionUnresolved.
/// Output allocations precede iteration; no successful or failed result borrows workspace.
[[nodiscard]] GmresResult<GmresSolution>
solve_linear(const GmresOperator& original, const FixedPreconditioner& preconditioner,
             const Rhs<double>& rhs, const Scales& scales, const Precision<double>& precision,
             const Accuracy& accuracy, const GmresLimits& limits, GmresPrepared& prepared,
             core::ResourcePlan& result_plan) noexcept;

} // namespace uvlas::numerics::linear
