#pragma once

#include <span>
#include <uvlas/core/instrumentation.hpp>
#include <uvlas/numerics/eigen_action.hpp>
#include <uvlas/numerics/eigen_metrics.hpp>

namespace uvlas::numerics::eigen {
namespace detail {
class GeneralizedFactory;
}

/// Borrow a live original real matrix and its complete physical/numerical identity.
/// The source remains alive, unmoved and unchanged throughout the synchronous call.
class GeneralizedOperator {
public:
    /// Bind a matrix; stiffness and mass keep their own numerical generations.
    GeneralizedOperator(const DenseReal& values, Context identity) noexcept
        : values_(values), identity_(identity) {}
    /// Temporaries cannot supply borrowed coefficients.
    GeneralizedOperator(DenseReal&&, Context) = delete;
    /// Const temporaries cannot supply borrowed coefficients.
    GeneralizedOperator(const DenseReal&&, Context) = delete;
    /// Read the original immutable matrix during the call.
    [[nodiscard]] const DenseReal& values() const noexcept {
        return values_;
    }
    /// Copy all source identities without changing the matrix generation.
    [[nodiscard]] Context context() const noexcept {
        return identity_;
    }

private:
    const DenseReal& values_;
    Context identity_;
};

/// Model-defined component scales and units; no scale is chosen from computed modes.
/// Coefficient dimensions follow K_ij = [m0*s_lambda]/[phi_i*phi_j] and
/// M_ij = [m0]/[phi_i*phi_j]. Dense storage has no independent unit metadata
/// to cross-check; the caller declares that its coefficients obey these units.
struct GeneralizedScales {
    /// Supply every fixed physical scale, unit and identity explicitly.
    GeneralizedScales(MetricScales metric_scales, double quadratic_scale,
                      std::span<const linear::UnitDimensions> units,
                      linear::UnitDimensions lambda_units, linear::UnitDimensions quadratic_units,
                      Context identity) noexcept
        : metric(metric_scales),
          mass_scale(quadratic_scale),
          component_units(units),
          eigenvalue_units(lambda_units),
          mass_units(quadratic_units),
          context(identity) {}
    MetricScales metric; ///< Fixed positive D and eigenvalue scale s_lambda.
    double mass_scale;   ///< Positive scale of phi^T M phi, including coordinate units.
    std::span<const linear::UnitDimensions> component_units; ///< One SI exponent vector per row.
    linear::UnitDimensions eigenvalue_units; ///< SI exponent vector of physical lambda.
    linear::UnitDimensions mass_units;       ///< SI exponent vector of phi^T M phi.
    Context context; ///< Matches stiffness identity; mass may have its own matrix generation.
};

/// Explicit finite capacities; external estimates are separate from owned payload accounting.
struct GeneralizedLimits {
    /// Select finite capacities and an explicit external admission mode.
    GeneralizedLimits(core::Index order, core::Index results, ActionLimits applications,
                      linear::BackendMemoryLimit memory) noexcept
        : maximum_order(order),
          maximum_results(results),
          actions(applications),
          backend_memory(memory) {}
    core::Index maximum_order;   ///< Maximum admitted order, including zero.
    core::Index maximum_results; ///< Capacity for every eigenvalue; no clipping occurs.
    /// One forward application evaluates the original pair (K*phi,M*phi) for one column.
    /// Projection requalification is charged again. Adjoint/inverse are unused; attempts
    /// bound candidate extraction calls independently of Eigen internal iterations.
    ActionLimits actions;
    linear::BackendMemoryLimit backend_memory; ///< Explicit forecast or proven-bound policy.
};

/// Algebraic extraction leaves physical mode coverage unresolved.
struct GeneralizedPhysicalCoverageUnresolved {};

class GeneralizedSpectrum;
/// Synchronous evaluator of a necessary check on the exact physical returned packet.
class GeneralizedConsumer {
public:
    /// Destroy through the evaluator interface.
    virtual ~GeneralizedConsumer() = default;
    /// Evaluate fresh retained and returned packets, including projection when requested.
    /// Certificates must match expected context and request budget exactly; exceptions are
    /// failures.
    [[nodiscard]] virtual core::Result<core::ErrorCertificate>
    assess(const GeneralizedSpectrum& retained, const GeneralizedSpectrum& returned,
           const linear::ConsumerCheck& request, core::CertificateContext expected) const = 0;
};

/// Explicit returned representation, numerical generations and all necessary consumers.
struct GeneralizedPrecision {
    /// Supply the exact input representation, return policy and complete consumer set.
    GeneralizedPrecision(core::Precision representation, linear::ReturnedPrecision returned_values,
                         std::span<const linear::ConsumerCheck> checks,
                         const GeneralizedConsumer* consumer) noexcept
        : input(representation),
          returned(returned_values),
          required_checks(checks),
          evaluator(consumer) {}
    core::Precision input; ///< Currently Binary64 only; TwoWord is rejected before computation.
    linear::ReturnedPrecision
        returned; ///< Retained or explicit Binary64 projection, never implicit.
    std::span<const linear::ConsumerCheck> required_checks; ///< All checks must pass.
    const GeneralizedConsumer*
        evaluator; ///< Required when checks are nonempty; synchronous borrow.
};

/// Coarse wall-time measurements can be disabled without reading the clock.
enum class GeneralizedMeasurement : std::uint8_t {
    Disabled, ///< No clock reads; all timing counters remain zero.
    Enabled   ///< Record product stages without changing numerical policy.
};
/// Nonoverlapping stages measured on the product extraction path.
enum class GeneralizedStage : std::uint8_t {
    Preparation,       ///< Entry validation and owned/scratch allocation.
    MassQualification, ///< Canonical original mass and rigorous qualification.
    Backend,           ///< Scaled candidate construction, LLT and generalized eigensolver.
    Qualification,     ///< Original actions and actual returned-value certification.
    Consumer,          ///< Necessary consumer evaluations.
    Count              ///< Sentinel, not a measurable stage.
};
/// Owned fixed-size measurement counters, separate from numerical certificates.
struct GeneralizedTimings {
    bool enabled = false; ///< Disabled counters are zero and require no clock reads.
    std::array<core::StageTiming, static_cast<std::size_t>(GeneralizedStage::Count)>
        stages{}; ///< Completed intervals.
};

/// Every budget and reporting identity supplied by the caller before extraction.
struct GeneralizedRequest {
    /// Construct a complete request with no implicit numerical or resource policy.
    GeneralizedRequest(
        GeneralizedScales physical_scales, GeneralizedMetricRequest requirements,
        GeneralizedPrecision precision_policy, GeneralizedLimits capacities,
        StateDiagnosticContext identities, bool physical_coverage,
        GeneralizedMeasurement measurements = GeneralizedMeasurement::Disabled) noexcept
        : scales(physical_scales),
          accuracy(requirements),
          precision(precision_policy),
          limits(capacities),
          state(identities),
          require_physical_coverage(physical_coverage),
          measurement(measurements) {}
    GeneralizedScales scales; ///< Fixed model scales and SI units.
    GeneralizedMetricRequest
        accuracy;                   ///< Independent residual, mass-normalization and Gram budgets.
    GeneralizedPrecision precision; ///< Input and actual returned representations/generations.
    GeneralizedLimits limits;       ///< Owned dimensions and external memory admission.
    StateDiagnosticContext state;   ///< Optional reporting IDs, not physical ancestry validation.
    bool require_physical_coverage; ///< True rejects: algebraic extraction cannot prove physical
                                    ///< coverage.
    GeneralizedMeasurement measurement; ///< Optional coarse timing; never changes numerical policy.
};

/// Owned full real generalized spectrum qualified against immutable original K and M.
/// Success proves a full algebraic basis, not physical mode coverage or stability.
/// Results remain valid after source matrices and evaluator are destroyed.
class GeneralizedSpectrum {
public:
    /// Transfer independent owned values and proofs.
    GeneralizedSpectrum(GeneralizedSpectrum&&) noexcept = default;
    /// Replace this result; existing borrows expire.
    GeneralizedSpectrum& operator=(GeneralizedSpectrum&&) noexcept = default;
    /// Implicit ownership copies are forbidden.
    GeneralizedSpectrum(const GeneralizedSpectrum&) = delete;
    /// Implicit ownership assignment is forbidden.
    GeneralizedSpectrum& operator=(const GeneralizedSpectrum&) = delete;
    /// Release arrays without depending on a prepared backend or callback.
    ~GeneralizedSpectrum() = default;
    /// Full signed physical eigenvalues, ascending; zero and negative modes are preserved.
    [[nodiscard]] std::span<const double> eigenvalues() const& noexcept {
        return values_.values();
    }
    /// No borrowed view from a temporary result.
    [[nodiscard]] std::span<const double> eigenvalues() const&& = delete;
    /// Physical column-major M-normalized full basis, n by n.
    [[nodiscard]] const DenseReal& basis() const& noexcept {
        return basis_;
    }
    /// No borrowed matrix from a temporary result.
    [[nodiscard]] const DenseReal& basis() const&& = delete;
    /// Full original stiffness snapshot associated with the proof.
    [[nodiscard]] const DenseReal& stiffness() const& noexcept {
        return stiffness_;
    }
    /// Temporaries cannot supply borrowed result storage.
    [[nodiscard]] const DenseReal& stiffness() const&& = delete;
    /// Full original mass snapshot associated with its positive-definite proof.
    [[nodiscard]] const DenseReal& mass() const& noexcept {
        return mass_;
    }
    /// Temporaries cannot supply borrowed result storage.
    [[nodiscard]] const DenseReal& mass() const&& = delete;
    /// Source stiffness context; mass retains its distinct numerical generation.
    [[nodiscard]] Context stiffness_context() const noexcept {
        return stiffness_context_;
    }
    /// Original mass context, without generation relabelling.
    [[nodiscard]] Context mass_context() const noexcept {
        return mass_context_;
    }
    /// Actual returned values' generation; projected results have their distinct target generation.
    [[nodiscard]] core::NumericalGeneration generation() const noexcept {
        return generation_;
    }
    /// Owned evaluated residual, normalization and full Gram certificates.
    [[nodiscard]] const GeneralizedMetrics& metrics() const& noexcept {
        return metrics_;
    }
    /// Temporaries cannot supply borrowed result storage.
    [[nodiscard]] const GeneralizedMetrics& metrics() const&& = delete;
    /// Original-M Rump qualification shift, zero for an empty matrix.
    [[nodiscard]] double mass_qualification_shift() const noexcept {
        return mass_shift_;
    }
    /// Model's fixed diagonal coordinate scales, copied into this packet.
    [[nodiscard]] std::span<const double> component_scales() const& noexcept {
        return scales_.values();
    }
    /// Temporaries cannot supply borrowed result storage.
    [[nodiscard]] std::span<const double> component_scales() const&& = delete;
    /// Fixed physical eigenvalue scale.
    [[nodiscard]] double eigenvalue_scale() const noexcept {
        return eigenvalue_scale_;
    }
    /// Fixed scale of the original mass quadratic form.
    [[nodiscard]] double mass_scale() const noexcept {
        return mass_scale_;
    }
    /// Fixed SI component units in the copied physical row order.
    [[nodiscard]] std::span<const linear::UnitDimensions> component_units() const& noexcept {
        return units_.values();
    }
    /// Temporaries cannot supply borrowed result storage.
    [[nodiscard]] std::span<const linear::UnitDimensions> component_units() const&& = delete;
    /// SI dimensions of lambda, preserved without frequency/square-root conversion.
    [[nodiscard]] linear::UnitDimensions eigenvalue_units() const noexcept {
        return eigenvalue_units_;
    }
    /// SI dimensions of the mass quadratic form.
    [[nodiscard]] linear::UnitDimensions mass_units() const noexcept {
        return mass_units_;
    }
    /// Number of necessary consumer certificates in this packet.
    [[nodiscard]] std::size_t consumer_count() const noexcept {
        return consumers_.size();
    }
    /// Borrow the certificate at a valid index; invalid indices return InvalidInput.
    [[nodiscard]] core::Result<core::ErrorCertificate> consumer(std::size_t index) const noexcept;
    /// Explicitly unresolved physical coverage; this entry never creates a physical qualification.
    [[nodiscard]] bool physical_coverage_qualified() const noexcept {
        return false;
    }
    /// Explicit unresolved physical-coverage state owned by this result.
    [[nodiscard]] GeneralizedPhysicalCoverageUnresolved physical_coverage() const noexcept {
        return coverage_;
    }
    /// Caller-supplied reporting identities, without ancestry inference.
    [[nodiscard]] StateDiagnosticContext state_context() const noexcept {
        return state_;
    }
    /// Retained/projected policy and both generations when projection was requested.
    [[nodiscard]] linear::ReturnedPrecision returned_precision() const noexcept {
        return returned_;
    }
    /// Fixed internal eigensolver stopping rule, not measured iterations.
    [[nodiscard]] core::Index backend_iteration_limit() const noexcept {
        return basis_.shape().rows * 30;
    }
    /// Extraction diagnostics; backend internal iteration count remains unavailable.
    [[nodiscard]] AttemptDiagnostics diagnostics() const noexcept {
        return diagnostics_;
    }
    /// Copy completed product-stage timings; disabled recording performs no clock reads.
    [[nodiscard]] GeneralizedTimings timings() const noexcept {
        return timings_;
    }
    /// Copy every requested capacity and action budget used for this extraction.
    [[nodiscard]] GeneralizedLimits limits() const noexcept {
        return limits_;
    }
    /// Forecast admitted before external allocations; not an RSS bound.
    [[nodiscard]] linear::BackendAdmission backend_admission() const noexcept {
        return admission_;
    }

private:
    friend class detail::GeneralizedFactory;
    GeneralizedSpectrum(DenseReal stiffness, DenseReal mass, DenseReal basis,
                        core::PlannedBuffer<double> values, core::PlannedBuffer<double> scales,
                        core::PlannedBuffer<linear::UnitDimensions> units,
                        core::PlannedBuffer<std::optional<core::ErrorCertificate>> consumers,
                        Context stiffness_context, Context mass_context,
                        core::NumericalGeneration generation, GeneralizedMetrics metrics,
                        double shift, double eigenvalue_scale, double mass_scale,
                        linear::UnitDimensions eigenvalue_units, linear::UnitDimensions mass_units,
                        AttemptDiagnostics diagnostics, linear::BackendAdmission admission,
                        StateDiagnosticContext state, linear::ReturnedPrecision returned,
                        GeneralizedLimits limits) noexcept;
    DenseReal stiffness_, mass_, basis_;
    core::PlannedBuffer<double> values_, scales_;
    core::PlannedBuffer<linear::UnitDimensions> units_;
    core::PlannedBuffer<std::optional<core::ErrorCertificate>> consumers_;
    Context stiffness_context_, mass_context_;
    core::NumericalGeneration generation_;
    GeneralizedMetrics metrics_;
    double mass_shift_, eigenvalue_scale_, mass_scale_;
    linear::UnitDimensions eigenvalue_units_, mass_units_;
    AttemptDiagnostics diagnostics_;
    linear::BackendAdmission admission_;
    GeneralizedLimits limits_;
    StateDiagnosticContext state_;
    linear::ReturnedPrecision returned_;
    GeneralizedPhysicalCoverageUnresolved coverage_;
    GeneralizedTimings timings_;
};

/// Extract every eigenvalue and a mass-orthonormal basis of K phi = lambda M phi.
/// Zero order performs no backend or consumer call; a nonempty necessary-consumer set
/// on an empty problem returns PrecisionUnresolved rather than fabricating certificates.
/// K/M must be exactly symmetric, share physical identities, and have independent live borrows.
/// M needs an original-matrix positive-definite qualification; K may be indefinite or singular.
/// plan cumulatively accounts owned/scratch payload. Input, prior results and parents stay
/// unchanged. Invalid shapes/budgets yield InvalidInput, unknown generations StaleGeneration,
/// unsupported precision PrecisionUnsupported, unqualified mass RankUnresolved, unmet numerical
/// proofs PrecisionUnresolved, and capacity/allocation failures ResourceExceeded. No exception
/// escapes.
[[nodiscard]] Result<GeneralizedSpectrum>
symmetric_generalized_spectrum(const GeneralizedOperator& stiffness,
                               const GeneralizedOperator& mass, const GeneralizedRequest& request,
                               core::ResourcePlan& plan) noexcept;
} // namespace uvlas::numerics::eigen
