#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <uvlas/core/state.hpp>
#include <uvlas/numerics/backend_memory.hpp>
#include <uvlas/numerics/precision.hpp>
#include <variant>

namespace uvlas::numerics::linear {
namespace detail {
class SolutionFactory;
}

/// Complete factorization identity plus the evaluated trial and physical component layouts.
/// Layout identities name component order, SI units and physical-coordinate meaning in the model.
struct Context {
    core::CacheKey
        cache; ///< All adopted model/policy/matrix generations and model content identity.
    core::TrialId trial;           ///< Trial at which the operator and constraints were evaluated.
    core::ScopeId unknown_layout;  ///< Physical unknown component ordering and units.
    core::ScopeId residual_layout; ///< Original residual component ordering and units.
    double time;                   ///< Finite physical evaluation time in SI seconds.
    /// Compare all identities and evaluation time; no pointer/shape shortcut is used.
    bool operator==(const Context&) const = default;
};

/// Right-hand side retaining both numerical words, borrowed for a synchronous solve only.
/// Its RHS generation is distinct from the numerical generation of the actual buffer.
template <StorageScalar Scalar>
class Rhs {
public:
    /// Bind a live owning buffer and its explicit operator/trial/RHS identity.
    Rhs(const PrecisionBuffer<Scalar>& values, Context context, core::RhsGeneration rhs) noexcept
        : values_(values), context_(context), rhs_(rhs) {}
    /// A temporary buffer cannot provide a retained borrow.
    Rhs(PrecisionBuffer<Scalar>&&, Context, core::RhsGeneration) = delete;
    /// A const temporary buffer cannot provide a retained borrow.
    Rhs(const PrecisionBuffer<Scalar>&&, Context, core::RhsGeneration) = delete;
    /// Read the exact source while the original buffer remains alive and unmoved.
    [[nodiscard]] const PrecisionBuffer<Scalar>& values() const noexcept {
        return values_;
    }
    /// Copy all operator/trial/layout metadata expected by this RHS.
    [[nodiscard]] Context context() const noexcept {
        return context_;
    }
    /// Request always includes every matrix-key field even when only the RHS changed.
    [[nodiscard]] core::CacheRequest cache_request() const noexcept {
        return {context_.cache, rhs_};
    }

private:
    const PrecisionBuffer<Scalar>& values_;
    Context context_;
    core::RhsGeneration rhs_;
};

/// Exponents of the seven SI base dimensions in the model's fixed component order.
using UnitDimensions = std::array<std::int8_t, 7>;
/// Model-defined physical scales and units, broadcast across every RHS column.
/// Both value spans are positive finite and every span has the operator's row count.
/// Inputs are borrowed only for a synchronous prepare/solve call; results copy their metadata.
struct Scales {
    /// Bind model-owned component metadata for a synchronous preparation or solve.
    Scales(std::span<const double> x, std::span<const double> r, std::span<const UnitDimensions> xu,
           std::span<const UnitDimensions> ru, Context identity) noexcept
        : unknown(x), residual(r), unknown_units(xu), residual_units(ru), context(identity) {}
    std::span<const double> unknown; ///< Positive physical coordinate scales, fixed before solving.
    std::span<const double>
        residual; ///< Positive original-residual scales, including zero RHS cases.
    std::span<const UnitDimensions> unknown_units;  ///< Physical units of each unknown component.
    std::span<const UnitDimensions> residual_units; ///< Physical units of each residual component.
    Context context; ///< Scale and layout generations must match the operator.
};
/// Norm applied independently to each RHS's normalized residual components.
enum class Norm : std::uint8_t {
    Infinity, ///< Maximum component magnitude.
    L2        ///< Euclidean norm, evaluated without squaring extreme physical values.
};
/// Fixed accuracy requirements; no scale or threshold is chosen from the computed solution.
/// For each RHS, tau_i=atol*s_i+rtol*abs(actual b_i), then norm(r_i/tau_i)<=1.
/// Zero tau_i requires exactly zero residual. Conservative error bounds may be unresolved.
struct Accuracy {
    Norm norm;                 ///< Per-column residual norm; all columns must independently pass.
    double absolute_tolerance; ///< Nonnegative finite multiplier of model residual scales.
    double relative_tolerance; ///< Nonnegative finite multiplier of the actual RHS magnitude.
    double componentwise_backward_tolerance; ///< Separate nonnegative finite backward-error budget.
    core::ObservationId residual_observation; ///< Aggregate worst-column residual requirement.
    core::ObservationId backward_observation; ///< Separate componentwise-backward requirement.
};

/// Keep the solver's actual representation; a TwoWord result is never silently rounded.
struct RetainedPrecision {
    core::NumericalGeneration
        generation; ///< Generation assigned to the exact returned retained values.
};
/// Explicit request for a checked binary64 projection after retained-solution qualification.
struct Binary64Projection {
    core::NumericalGeneration
        retained_generation; ///< Generation of retained values being projected.
    core::NumericalGeneration
        target_generation; ///< Distinct generation for the exact projected values.
};
/// A retained result or an explicitly requested projection; no unused target-generation field.
using ReturnedPrecision = std::variant<RetainedPrecision, Binary64Projection>;
/// Distinguish physical outputs from constraint and reaction checks.
enum class CheckKind : std::uint8_t {
    Observation, ///< Physical output or other required numerical observation.
    Constraint,  ///< Independent constraint-position/velocity requirement.
    Reaction     ///< Independent reaction magnitude/sign/work requirement.
};
/// One necessary consumer requirement, using the evaluator's declared physical scale.
struct ConsumerCheck {
    CheckKind kind;                  ///< Why this independent requirement is needed.
    core::ObservationId observation; ///< Quantity actually evaluated.
    core::ScopeId scope;             ///< Domain covered by this requirement.
    double budget; ///< Nonnegative finite permitted error in the evaluator's explicit scale.
};
/// Synchronous C++ evaluator of every requested physical check on the exact candidate.
/// It must use both numerical words when declaring TwoWord support; Python callbacks
/// do not belong in numerical iterations. Caller owns the evaluator during the solve.
template <StorageScalar Scalar>
class ConsumerEvaluator {
public:
    /// Release through this interface.
    virtual ~ConsumerEvaluator() = default;
    /// Representation capability of the implemented numerical consumer.
    [[nodiscard]] virtual SupportedPrecision supported_precision() const noexcept = 0;
    /// Compute a fresh bound from the exact retained and consumed candidates.
    /// Both references denote the same values when no projection was requested.
    /// The returned certificate must match expected_context and request.budget exactly.
    /// A precomputed proof for another candidate must not be relabelled or reused.
    /// Unevaluated/Exceeded results cannot qualify a solution. Exceptions become Failure.
    [[nodiscard]] virtual core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<Scalar>& retained_source, const PrecisionBuffer<Scalar>& candidate,
           const Context& system, const ConsumerCheck& request,
           core::CertificateContext expected_context) const = 0;
};
/// Explicit returned representation and complete necessary consumer checks.
/// An empty requirement span requests algebraic qualification only; no physical acceptance
/// is inferred from it. All borrowed requirements/evaluator live through the synchronous call.
template <StorageScalar Scalar>
struct Precision {
    ReturnedPrecision returned; ///< Requested actual return representation and generation(s).
    std::span<const ConsumerCheck>
        required_checks; ///< Every necessary physical/constraint/reaction check.
    const ConsumerEvaluator<Scalar>&
        evaluator; ///< Evaluates the actual retained/projected candidate.
};
/// Explicit capacities and refinement limit. Refinements across both representations total at
/// most 5. Backend admission is separate from ResourcePlan's cumulative owned-array payload
/// accounting.
struct Limits {
    core::Index maximum_order; ///< Nonnegative maximum matrix order admitted before decomposition.
    core::Index
        maximum_rhs_columns; ///< Prepared RHS capacity; exceeding it requires new preparation.
    core::Index maximum_update_rank; ///< Nonnegative maximum complex update rank.
    unsigned maximum_refinements;    ///< 0..5; promotion does not reset the count.
    BackendMemoryLimit
        backend_memory; ///< Explicit forecast/proven-bound request, readmitted per solve.
};

/// No condition estimate was evaluated for this operator.
struct NotEvaluated {};
/// An uncertified estimate of reciprocal one-norm condition, never a rank threshold or bound.
struct ReciprocalOneNormEstimate {
    double value; ///< Finite nonnegative estimated reciprocal kappa_1.
};
/// Diagnostic estimate state; neither alternative asserts a guaranteed condition upper bound.
using ConditionDiagnostic = std::variant<NotEvaluated, ReciprocalOneNormEstimate>;
/// Conservative numerical metrics for one independently checked RHS column.
struct ColumnMetrics {
    double normalized_original_residual_upper =
        0; ///< Upper bound in the requested tau normalization.
    double componentwise_backward_error_upper =
        0; ///< Upper bound using actual retained-value magnitudes.
    double scaled_correction_upper =
        0; ///< Last correction in model-defined unknown component scales.
};
/// Owned necessary-check result, kept distinct from algebraic residual/backward certification.
struct ConsumerProof {
    ConsumerCheck request;             ///< Exact requirement actually checked.
    core::ErrorCertificate assessment; ///< Bound, context, representation and achieved state.
};
/// How the returned numerical words were produced.
enum class RetentionPath : std::uint8_t {
    Binary64,         ///< Binary64 retention met the requested criteria.
    InitiallyTwoWord, ///< Both RHS words required two-word retention from the start.
    PromotedTwoWord   ///< Binary64 correction stagnated and retention was promoted.
};

/// Owning qualification metadata bound to one exact solver candidate, not workspace storage.
class Certificate {
public:
    /// Transfer owned diagnostics; source borrows are invalidated.
    Certificate(Certificate&&) noexcept = default;
    /// Replace owned diagnostics.
    Certificate& operator=(Certificate&&) noexcept = default;
    /// Implicit copying of owned diagnostic arrays is forbidden.
    Certificate(const Certificate&) = delete;
    /// Copy assignment is forbidden.
    Certificate& operator=(const Certificate&) = delete;
    /// Release metadata storage.
    ~Certificate() = default;
    /// Full cached operator identity and physical evaluation context.
    [[nodiscard]] Context context() const noexcept {
        return context_;
    }
    /// Explicit RHS generation, distinct from matrix/factor identity.
    [[nodiscard]] core::RhsGeneration rhs_generation() const noexcept {
        return rhs_;
    }
    /// Actual numerical generation of the RHS buffer consumed by the solve.
    [[nodiscard]] core::NumericalGeneration rhs_buffer_generation() const noexcept {
        return rhs_buffer_;
    }
    /// Requested fixed numerical accuracy requirements.
    [[nodiscard]] Accuracy accuracy() const noexcept {
        return accuracy_;
    }
    /// Borrow the retained candidate's original-residual proof.
    [[nodiscard]] const core::ErrorCertificate& retained_residual() const& noexcept {
        return retained_residual_;
    }
    /// Borrow its independently evaluated backward-error proof.
    [[nodiscard]] const core::ErrorCertificate& retained_backward() const& noexcept {
        return retained_backward_;
    }
    /// Borrow per-RHS retained diagnostics, distinct from the returned projection.
    [[nodiscard]] std::span<const ColumnMetrics> retained_columns() const& noexcept {
        return retained_columns_.values();
    }
    /// Read owned physical coordinate scales after input/model buffers have expired.
    [[nodiscard]] std::span<const double> unknown_scales() const& noexcept {
        return unknown_scales_.values();
    }
    /// Read owned original-residual physical scales.
    [[nodiscard]] std::span<const double> residual_scales() const& noexcept {
        return residual_scales_.values();
    }
    /// Read owned unknown SI dimension exponents in component order.
    [[nodiscard]] std::span<const UnitDimensions> unknown_units() const& noexcept {
        return unknown_units_.values();
    }
    /// Read owned residual SI dimension exponents in component order.
    [[nodiscard]] std::span<const UnitDimensions> residual_units() const& noexcept {
        return residual_units_.values();
    }
    /// Borrow the exact returned candidate's aggregate original-residual proof.
    [[nodiscard]] const core::ErrorCertificate& residual() const& noexcept {
        return residual_;
    }
    /// Borrow independent backward-error proof for the same candidate.
    [[nodiscard]] const core::ErrorCertificate& backward() const& noexcept {
        return backward_;
    }
    /// Borrow per-RHS diagnostics for this owned result.
    [[nodiscard]] std::span<const ColumnMetrics> columns() const& noexcept {
        return columns_.values();
    }
    /// A temporary cannot provide this borrowed diagnostic view.
    [[nodiscard]] const core::ErrorCertificate& retained_residual() const&& = delete;
    /// A temporary cannot provide this borrowed diagnostic view.
    [[nodiscard]] const core::ErrorCertificate& retained_backward() const&& = delete;
    /// A temporary cannot provide this borrowed diagnostic view.
    [[nodiscard]] const core::ErrorCertificate& residual() const&& = delete;
    /// A temporary cannot provide this borrowed diagnostic view.
    [[nodiscard]] const core::ErrorCertificate& backward() const&& = delete;
    /// A temporary cannot provide this borrowed diagnostic view.
    [[nodiscard]] std::span<const ColumnMetrics> retained_columns() const&& = delete;
    /// A temporary cannot provide this borrowed diagnostic view.
    [[nodiscard]] std::span<const ColumnMetrics> columns() const&& = delete;
    /// A temporary cannot provide this borrowed diagnostic view.
    [[nodiscard]] std::span<const double> unknown_scales() const&& = delete;
    /// A temporary cannot provide this borrowed diagnostic view.
    [[nodiscard]] std::span<const double> residual_scales() const&& = delete;
    /// A temporary cannot provide this borrowed diagnostic view.
    [[nodiscard]] std::span<const UnitDimensions> unknown_units() const&& = delete;
    /// A temporary cannot provide this borrowed diagnostic view.
    [[nodiscard]] std::span<const UnitDimensions> residual_units() const&& = delete;
    /// Number of required physical check proofs owned by this result.
    [[nodiscard]] core::Index consumer_count() const noexcept;
    /// Copy one immutable proof, or InvalidInput for a negative/out-of-range index.
    [[nodiscard]] core::Result<ConsumerProof> consumer_proof(core::Index index) const noexcept;
    /// Actual admission performed before the cached backend was prepared.
    [[nodiscard]] BackendAdmission preparation_admission() const noexcept {
        return prepared_admission_;
    }
    /// New comparison against this solve call's requested backend mode/limit.
    [[nodiscard]] BackendAdmission solve_admission() const noexcept {
        return solve_admission_;
    }
    /// Actual number of attempted corrections, never reset by precision promotion.
    [[nodiscard]] unsigned refinements() const noexcept {
        return refinements_;
    }
    /// Retention path, distinct from the requested returned precision.
    [[nodiscard]] RetentionPath retention_path() const noexcept {
        return retention_;
    }
    /// Operator-specific heuristic reciprocal-condition diagnostic; no guarantee is inferred.
    [[nodiscard]] ConditionDiagnostic condition() const noexcept {
        return condition_;
    }

private:
    friend class detail::SolutionFactory;
    Certificate(Context, core::RhsGeneration, core::NumericalGeneration, Accuracy,
                core::ErrorCertificate, core::ErrorCertificate, core::ErrorCertificate,
                core::ErrorCertificate, core::PlannedBuffer<ColumnMetrics>,
                core::PlannedBuffer<ColumnMetrics>, core::PlannedBuffer<double>,
                core::PlannedBuffer<double>, core::PlannedBuffer<UnitDimensions>,
                core::PlannedBuffer<UnitDimensions>,
                core::PlannedBuffer<std::optional<ConsumerProof>>, BackendAdmission,
                BackendAdmission, unsigned, RetentionPath, ConditionDiagnostic) noexcept;
    Context context_;
    core::RhsGeneration rhs_;
    core::NumericalGeneration rhs_buffer_;
    Accuracy accuracy_;
    core::ErrorCertificate retained_residual_;
    core::ErrorCertificate retained_backward_;
    core::ErrorCertificate residual_;
    core::ErrorCertificate backward_;
    core::PlannedBuffer<ColumnMetrics> retained_columns_;
    core::PlannedBuffer<ColumnMetrics> columns_;
    core::PlannedBuffer<double> unknown_scales_;
    core::PlannedBuffer<double> residual_scales_;
    core::PlannedBuffer<UnitDimensions> unknown_units_;
    core::PlannedBuffer<UnitDimensions> residual_units_;
    core::PlannedBuffer<std::optional<ConsumerProof>> consumers_;
    BackendAdmission prepared_admission_;
    BackendAdmission solve_admission_;
    unsigned refinements_;
    RetentionPath retention_;
    ConditionDiagnostic condition_;
};

/// Owns the exact qualified retained/projected numerical values and their certificate.
/// Its immutable values survive subsequent prepared-workspace reuse.
template <StorageScalar Scalar>
class Solution {
public:
    /// Transfer ownership and invalidate source borrows.
    Solution(Solution&&) noexcept = default;
    /// Replace owned values and diagnostics.
    Solution& operator=(Solution&&) noexcept = default;
    /// Copying is forbidden.
    Solution(const Solution&) = delete;
    /// Copy assignment is forbidden.
    Solution& operator=(const Solution&) = delete;
    /// Release values independently of the prepared solver's lifetime.
    ~Solution() = default;
    /// Read the exact returned representation until move/assignment/destruction.
    [[nodiscard]] const PrecisionBuffer<Scalar>& values() const& noexcept {
        return values_;
    }
    /// A temporary cannot provide borrowed values.
    [[nodiscard]] const PrecisionBuffer<Scalar>& values() const&& = delete;
    /// Read this candidate's immutable qualification until move/assignment/destruction.
    [[nodiscard]] const Certificate& certificate() const& noexcept {
        return certificate_;
    }
    /// A temporary cannot provide borrowed diagnostics.
    [[nodiscard]] const Certificate& certificate() const&& = delete;

private:
    friend class detail::SolutionFactory;
    Solution(PrecisionBuffer<Scalar> values, Certificate certificate) noexcept
        : values_(std::move(values)), certificate_(std::move(certificate)) {}
    PrecisionBuffer<Scalar> values_;
    Certificate certificate_;
};

/// Borrowed nonsymmetric real operator. Matrix must remain live during prepare/solve.
struct DenseOperator {
    /// Bind live original coefficients and the full cached evaluation identity.
    DenseOperator(const DenseReal& a, Context identity) noexcept : values(a), context(identity) {}
    const DenseReal& values; ///< Original binary64 coefficients, never the overwritten LU storage.
    Context context;         ///< Complete cache identity, with Binary64 cached factor arithmetic.
};
/// Borrowed A=A0+U*V^T operator, explicitly without conjugation.
/// Effective binary64 coefficients are defined by fixed-rank-order coefficient evaluation.
/// Original residual and backward error use those same rounded coefficients, not staged U(V^Tx).
struct ColumnUpdatedOperator {
    /// Bind live real base and both unconjugated complex update factors.
    ColumnUpdatedOperator(const DenseReal& a, const DenseComplex& left, const DenseComplex& right,
                          Context identity) noexcept
        : base(a), u(left), v(right), context(identity) {}
    const DenseReal& base; ///< Original real n-by-n base.
    const DenseComplex& u; ///< Complex n-by-r update factor.
    const DenseComplex& v; ///< Complex n-by-r transpose factor, not an adjoint factor.
    Context context;       ///< Identity covers base, both update factors, policies and layouts.
};
class DensePrepared;
class UpdatedPrepared;
/// Solve using prepared factors/workspaces; all output allocation occurs before the inner kernel.
/// Inputs and earlier results stay unchanged. Generation/content mismatches are StaleGeneration.
/// result_plan accounts owning output buffers; no returned array borrows prepared scratch.
[[nodiscard]] core::Result<Solution<double>> solve_linear(const DenseOperator&, const Rhs<double>&,
                                                          const Scales&, const Precision<double>&,
                                                          const Accuracy&, const Limits&,
                                                          DensePrepared&,
                                                          core::ResourcePlan& result_plan) noexcept;
/// Solve the original effective complex system using real LU and its prepared small complex Schur.
/// Woodbury is an approximate inverse action; qualification always reevaluates original
/// coefficients.
[[nodiscard]] core::Result<Solution<std::complex<double>>>
solve_linear(const ColumnUpdatedOperator&, const Rhs<std::complex<double>>&, const Scales&,
             const Precision<std::complex<double>>&, const Accuracy&, const Limits&,
             UpdatedPrepared&, core::ResourcePlan& result_plan) noexcept;

} // namespace uvlas::numerics::linear
