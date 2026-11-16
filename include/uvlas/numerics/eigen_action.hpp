#pragma once

#include <complex>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <uvlas/numerics/backend_memory.hpp>
#include <uvlas/numerics/eigen_common.hpp>

namespace uvlas::numerics::eigen {

/// Mathematical action requested from one frozen original problem.
enum class ActionKind : std::uint8_t {
    Forward, ///< Original operator in the declared coordinates.
    Adjoint, ///< Conjugate transpose under the same declared inner product.
    Inverse  ///< Inverse of the declared shifted operator, with its own error bound.
};

/// Complete original identity and the fixed spectral transformation, shared by all actions.
struct ActionStamp {
    /// Supply the original identity and fixed transformation explicitly.
    ActionStamp(Context context, core::NumericalGeneration generation,
                std::complex<double> physical_shift) noexcept
        : original(context), transformation(generation), shift(physical_shift) {}
    Context original; ///< Original model, coefficients, layout, scales, time and backend policy.
    core::NumericalGeneration transformation; ///< Generation of shift/factor/transform policy.
    std::complex<double> shift;               ///< Physical shift; zero for an unshifted policy.
    /// Compare every identity field and the physical shift.
    bool operator==(const ActionStamp&) const = default;
};

/// Explicit budgets for one complete extraction, including retries and qualification.
struct ActionLimits {
    std::uint64_t forward;  ///< Maximum started original-operator callbacks.
    std::uint64_t adjoint;  ///< Maximum started adjoint callbacks.
    std::uint64_t inverse;  ///< Maximum started inverse callbacks.
    std::uint64_t attempts; ///< Maximum started extraction attempts.
};

/// Synchronous immutable operator access; implementations freeze coefficients and policy.
/// Calls allocate nothing and must fill every output entry. A/adjoint/inverse all refer
/// to the same stamp. The owner separately accounts opaque retained/workspace memory.
/// const and a matching stamp do not prove that external coefficients stayed unchanged.
template <StorageScalar Scalar>
class ActionProvider {
public:
    /// Release an implementation through this interface.
    virtual ~ActionProvider() = default;
    /// Logical square order, without padding.
    [[nodiscard]] virtual core::Index order() const noexcept = 0;
    /// Frozen original problem and inverse transformation identity; no numerical work.
    [[nodiscard]] virtual ActionStamp stamp() const noexcept = 0;
    /// Input representation actually consumed for this action; no implicit low-word loss.
    [[nodiscard]] virtual SupportedPrecision
    supported_precision(ActionKind kind) const noexcept = 0;
    /// Check finite immutable state before use; no hidden iteration or unbounded work.
    [[nodiscard]] virtual core::Result<void> validate_frozen_state() const = 0;
    /// Return hi/tail and component bounds satisfying |action(x)-hi-tail| <= error_upper.
    /// Every span has order() entries, inputs and outputs do not overlap, and bounds
    /// are real, finite and nonnegative. Both input planes must be consumed. Outputs
    /// are scratch and may change on failure; no pointer/view may escape this call.
    [[nodiscard]] virtual core::Result<void> apply(ActionKind kind, WordView<Scalar> input,
                                                   std::span<Scalar> high, std::span<Scalar> low,
                                                   std::span<double> error_upper) const = 0;
};

/// Allocation-free per-extraction quota, fixed stamp and first-failure ownership.
/// Use exclusively from one thread. Reentrant calls fail and poison this session;
/// subsequent adapter calls return the first failure without invoking the producer.
/// Each extraction, including its final qualification, shares this single session.
template <StorageScalar Scalar>
class ActionSession {
public:
    /// Snapshot identities and check finite state. Empty order skips state validation.
    /// No provider is retained; each apply receives a live provider for this stamp.
    /// Reject negative order, nonfinite time/shift and unknown original-coefficient precision
    /// before producer validation. Known coefficient precision is separate from each input view.
    [[nodiscard]] static Result<ActionSession> make(const ActionProvider<Scalar>& provider,
                                                    ActionLimits limits,
                                                    StateDiagnosticContext state = {}) noexcept;
    /// Transfer counters and fixed identity, leaving the source inactive.
    /// Neither object may be moved while an action is in progress.
    ActionSession(ActionSession&& other) noexcept;
    /// Implicit copies would duplicate quota and are forbidden.
    ActionSession(const ActionSession&) = delete;
    /// Assignment is forbidden, preserving each extraction's first failure and budget.
    ActionSession& operator=(const ActionSession&) = delete;
    /// Move assignment is forbidden; create a new extraction session explicitly.
    ActionSession& operator=(ActionSession&&) = delete;
    /// Release value-only state independently of the provider and returned packets.
    ~ActionSession() = default;
    /// Charge before starting an extraction attempt, retaining previous attempts' cost.
    [[nodiscard]] Result<void> start_attempt() noexcept;
    /// Validate spans and current stamp, charge one action, then invoke the provider.
    /// Unknown kind/shape/alias rejects before charging. Callback failures and invalid
    /// output still consume their started call. Empty vectors invoke no callback and
    /// consume no action. Inputs must remain alive and unchanged during this call.
    [[nodiscard]] Result<void> apply(const ActionProvider<Scalar>& provider, ActionKind kind,
                                     WordView<Scalar> input, std::span<Scalar> high,
                                     std::span<Scalar> low, std::span<double> error_upper) noexcept;
    /// Copy known started counts. Backend iterations/convergence are not measured here.
    [[nodiscard]] AttemptDiagnostics diagnostics() const noexcept;
    /// Copy the first owned failure, if any; later calls cannot replace its diagnostics.
    [[nodiscard]] std::optional<EigenFailure> first_failure() const noexcept;
    /// Whether this object still owns its extraction quota.
    [[nodiscard]] bool active() const noexcept;

private:
    ActionSession(core::Index order, ActionStamp stamp, ActionLimits limits,
                  StateDiagnosticContext state) noexcept;
    [[nodiscard]] Result<void> remember(core::Failure cause) noexcept;
    [[nodiscard]] Result<void> fail(core::FailureCode code) noexcept;
    core::Index order_;
    ActionStamp stamp_;
    ActionLimits limits_;
    StateDiagnosticContext state_;
    AttemptDiagnostics diagnostics_;
    std::optional<EigenFailure> first_;
    bool active_ = true;
    bool busy_ = false;
};

/// One prospective column-major allocation, including its leading-dimension padding.
struct PayloadExtent {
    core::Index rows;              ///< Nonnegative logical rows.
    core::Index columns;           ///< Nonnegative logical columns.
    core::Index leading_dimension; ///< At least rows, nonnegative even for empty arrays.
    std::size_t element_bytes;     ///< Positive storage size of one element.
};
/// Sum explicit simultaneous owned payload without allocating or modifying ResourcePlan.
/// Negative shapes, invalid leading dimensions and zero element sizes are InvalidInput;
/// count/byte/sum overflow is ResourceExceeded. Empty shapes consume zero bytes.
/// This is neither a backend scratch estimate nor an RSS/whole-process upper bound.
[[nodiscard]] core::Result<std::uint64_t>
checked_payload_bytes(std::span<const PayloadExtent> arrays) noexcept;

} // namespace uvlas::numerics::eigen
