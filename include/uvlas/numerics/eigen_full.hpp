#pragma once

#include <complex>
#include <memory>
#include <optional>
#include <span>
#include <uvlas/core/instrumentation.hpp>
#include <uvlas/numerics/eigen_action.hpp>
#include <uvlas/numerics/eigen_metrics.hpp>

namespace uvlas::numerics::eigen {
namespace detail {
struct FullData;
struct FullIdentity;
class FullFactory;
} // namespace detail

/// Fixed physical component/eigenvalue scales and SI dimensions of every original coefficient.
struct FullScales {
    /// Bind immutable model scales and all coefficient units for this synchronous call.
    FullScales(std::span<const double> d, double scale,
               std::span<const linear::UnitDimensions> units, linear::UnitDimensions lambda_units,
               std::span<const linear::UnitDimensions> a_units) noexcept
        : components(d),
          eigenvalue_scale(scale),
          component_units(units),
          eigenvalue_units(lambda_units),
          coefficient_units(a_units) {}

    std::span<const double> components; ///< Positive finite D diagonal, fixed before extraction.
    double eigenvalue_scale; ///< Positive finite s_lambda, never selected from computed values.
    std::span<const linear::UnitDimensions> component_units; ///< One SI exponent vector per row.
    linear::UnitDimensions eigenvalue_units; ///< Dimensions of each physical eigenvalue.
    std::span<const linear::UnitDimensions> coefficient_units; ///< Packed column-major n*n units.
};
/// Independent dimensions and finite budgets for one complete full/derived extraction.
struct FullLimits {
    /// Explicit independent order/result capacities, shared action quotas and memory admission.
    FullLimits(core::Index order, core::Index results, ActionLimits budget,
               linear::BackendMemoryLimit memory) noexcept
        : maximum_order(order), maximum_results(results), actions(budget), backend_memory(memory) {}
    core::Index maximum_order;   ///< Nonnegative admissible order, including empty systems.
    core::Index maximum_results; ///< Maximum returned roots/pairs/subspace columns; never clipped.
    ActionLimits
        actions; ///< Includes every original/adjoint qualification and explicit projection.
    linear::BackendMemoryLimit backend_memory; ///< Forecast admission is not an allocation/RSS cap.
};
/// Policy for matching eigenvalues in dimensionless coordinates, never a separation proof.
struct FullAssociation {
    /// Explicit finite absolute and relative association policy; no hidden separation tolerance.
    FullAssociation(double a, double b) noexcept : absolute(a), relative(b) {}

    double absolute; ///< Finite nonnegative absolute association width.
    double relative; ///< Finite nonnegative relative association width.
};
/// Candidate generator of the original complete packet; derived packets retain its identity.
enum class FullBackend : std::uint8_t {
    RealSchurComplexified, ///< Real Schur followed by unitary conversion of real 2x2 blocks.
    ComplexSchur           ///< Complex Schur of an originally complex matrix.
};
/// The algebraic result being inspected; it does not identify a physical mode.
enum class FullKind : std::uint8_t {
    Spectrum, ///< All eigenvalues and a complete complex Schur factorization.
    Pairs,    ///< Selected physical right and optionally left eigenvectors.
    Cluster   ///< Selected invariant right and optionally independent left subspaces.
};
/// Numerical extraction does not establish physical model coverage.
enum class FullPhysicalCoverage : std::uint8_t { Unresolved };
/// Distinct algebraic success types; none can manufacture physical-coverage qualification.
struct FullSpectrumQualified {};
/// Selected original eigenpairs, separate from full Schur qualification.
struct FullPairsQualified {};
/// A requested invariant subspace, separate from the complete spectrum.
struct FullClusterQualified {};

/// Read-only synchronous view of the exact packet evaluated by a necessary consumer.
struct FullPacketView {
    FullKind kind; ///< Controls whether right/left are dimensionless Q or physical eigenvectors.
    bool left_requested; ///< True even for an empty requested left basis.
    core::Index rows;    ///< Original dimension.
    core::Index columns; ///< Returned pair/subspace count.
    std::span<const std::complex<double>> eigenvalues; ///< Sorted physical lambda values.
    std::span<const std::complex<double>> right; ///< Packed Q for Spectrum/Cluster, phi for Pairs.
    std::span<const std::complex<double>> left;  ///< Requested left Q/psi; empty when unrequested.
    std::span<const std::complex<double>> triangular; ///< Packed Schur T, empty for Pairs.
    std::span<const std::complex<double>>
        left_triangular; ///< Independent adjoint Schur T for a left cluster.
    std::span<const double>
        triangular_error_upper; ///< Enclosure of physical lambda/scale on T diagonal.
    std::span<const double>
        left_triangular_error_upper;      ///< Corresponding adjoint triangular enclosure.
    Context context;                      ///< Original model/layout/scale/policy identity.
    StateDiagnosticContext state;         ///< Reporting IDs, not a validation of state ancestry.
    core::NumericalGeneration generation; ///< Actual retained or projected packet being evaluated.
};
/// C++ consumer of the actual numerical packet; must not mutate/move borrowed owners or reenter.
class FullConsumer {
public:
    virtual ~FullConsumer() = default; ///< Destroy through the interface.
    /// Declare which actual representations the consumer evaluates.
    [[nodiscard]] virtual SupportedPrecision supported_precision() const noexcept = 0;
    /// Evaluate one exact observation/scope/budget on this packet, without retaining any view.
    [[nodiscard]] virtual core::Result<core::ErrorCertificate>
    evaluate(FullPacketView packet, MetricRequirement requirement) const = 0;
};
/// Whether extraction records wall-clock phase times.
enum class FullTimingMode : std::uint8_t {
    Disabled, ///< No clock reads or measurements.
    Enabled   ///< Record three named phases without changing numerical work.
};
/// Owned extraction phase measurements, separate from physical acceptance.
struct FullTimings {
    core::StageTiming preparation; ///< Validation, allocation and candidate matrix preparation.
    core::StageTiming backend; ///< Schur, unitary exchange and triangular candidate construction.
    core::StageTiming qualification; ///< Original-system metrics and necessary consumers.
};
/// Output, consumer and capacity controls for this complete-method family.
struct FullOutputRequest {
    /// Bind explicit output policy; optional recording defaults to disabled.
    FullOutputRequest(linear::ReturnedPrecision precision, FullLimits capacity,
                      std::span<const MetricRequirement> checks, const FullConsumer* consumer,
                      bool physical = false,
                      FullTimingMode recording = FullTimingMode::Disabled) noexcept
        : returned(precision),
          limits(capacity),
          consumers(checks),
          evaluator(consumer),
          require_physical_coverage(physical),
          timing(recording) {}
    linear::ReturnedPrecision returned; ///< Retained generation or distinct explicit projection.
    FullLimits limits;                  ///< No unlimited implicit size or action quota.
    std::span<const MetricRequirement> consumers; ///< Every necessary observation must qualify.
    const FullConsumer* evaluator;  ///< Synchronous evaluator; required for nonempty checks.
    bool require_physical_coverage; ///< Reject: this layer cannot qualify physical model coverage.
    FullTimingMode timing;          ///< Changes measurements only, never numerical policy.
};
/// Schur-basis accuracy criteria and independent output controls.
struct FullRequest : FullOutputRequest {
    /// Construct a fully specified request without implicit error budgets.
    FullRequest(BasisMetricRequest accuracy, linear::ReturnedPrecision precision,
                FullLimits capacity, std::span<const MetricRequirement> checks,
                const FullConsumer* consumer, bool physical = false,
                FullTimingMode recording = FullTimingMode::Disabled) noexcept
        : FullOutputRequest(precision, capacity, checks, consumer, physical, recording),
          basis(accuracy) {}
    /// Copy only common output controls for a pair request, without slicing its basis criteria.
    [[nodiscard]] FullOutputRequest output() const noexcept {
        return {returned, limits, consumers, evaluator, require_physical_coverage, timing};
    }
    BasisMetricRequest
        basis; ///< Whole original residual, complete Gram and optional backward test.
};
/// Requested pair checks and optional left vectors; all selected pairs must qualify.
struct FullPairRequest {
    /// Specify result controls, pair budgets and explicit left/independence requests.
    FullPairRequest(FullOutputRequest output, PairMetricRequest metric, FullAssociation policy,
                    bool request_left, bool full_basis) noexcept
        : common(output),
          pair(metric),
          association(policy),
          left(request_left),
          require_full_basis(full_basis) {}

    FullOutputRequest common; ///< Representation, capacities and necessary consumer checks.
    PairMetricRequest pair;   ///< Actual right and requested left residual/normalization budgets.
    FullAssociation association; ///< Explicit one-to-one root association policy.
    bool left;                   ///< Construct and qualify independent left eigenvectors when true.
    bool
        require_full_basis; ///< Require all n pairs and a certified left-inverse independence test.
};
/// Requested cluster, including an independently constructed adjoint subspace when needed.
struct FullClusterRequest {
    /// Specify subspace budgets, matching policy and the independent left-space request.
    FullClusterRequest(FullRequest output, FullAssociation policy, bool request_left) noexcept
        : common(output), association(policy), left(request_left) {}

    FullRequest common;          ///< Actual subspace budgets, generation, resources and consumers.
    FullAssociation association; ///< Matching policy for an independently extracted left space.
    bool left;                   ///< Never reuse the right subspace as its left counterpart.
};
/// Diagnostic sensitivity of an approximate left/right pair, never a rigorous forward-error bound.
struct FullSensitivity {
    /// Store a finite diagnostic estimate or an explicitly unresolved denominator.
    FullSensitivity(std::optional<double> value, bool unresolved) noexcept
        : estimate(value), condition_unresolved(unresolved) {}

    std::optional<double> estimate; ///< Missing when denominator is not positively established.
    bool condition_unresolved; ///< Distinguishes an unresolved sensitivity from an unrequested one.
};
/// Owned achieved pair proofs, with optional left proof and sensitivity only when requested.
struct FullPairProof {
    PairMetrics right; ///< Proof evaluated from actual physical right vector and original A.
    std::optional<PairMetrics> left;            ///< Proof for separately constructed physical psi.
    std::optional<FullSensitivity> sensitivity; ///< Diagnostic estimate, not an error certificate.
};

template <class Qualification>
class FullPacket;
using FullSpectrum = FullPacket<FullSpectrumQualified>; ///< Owned all-n Schur packet.
using FullPairs = FullPacket<FullPairsQualified>;       ///< Owned selected qualified pairs.
using FullCluster = FullPacket<FullClusterQualified>;   ///< Owned qualified invariant subspace.

/// Owned packet-bound index set. Its private identity survives owner destruction/address reuse.
class FullSelection {
public:
    /// Validate uniqueness/range and bind to exactly this full packet and generation.
    [[nodiscard]] static Result<FullSelection> make(const FullSpectrum& packet,
                                                    std::span<const core::Index> indices,
                                                    core::ResourcePlan& plan) noexcept;
    FullSelection(FullSelection&&) noexcept = default; ///< Transfer ownership and identity.
    FullSelection& operator=(FullSelection&&) noexcept = default; ///< Replace owned selection.
    FullSelection(const FullSelection&) = delete;                 ///< No implicit array copies.
    FullSelection& operator=(const FullSelection&) = delete;      ///< No implicit array copies.
    ~FullSelection() = default; ///< Release indices and immutable identity token.
private:
    friend class detail::FullFactory;
    FullSelection(std::shared_ptr<const detail::FullIdentity> identity,
                  core::PlannedBuffer<core::Index> indices) noexcept;
    std::shared_ptr<const detail::FullIdentity> identity_;
    core::PlannedBuffer<core::Index> indices_;
};

/// Own all numerical results, original coefficient/scale/unit snapshots and achieved proofs.
/// No getter exposes mutable storage. Borrows expire on owner move, assignment or destruction.
/// Calls using this packet must not concurrently move/destroy it or invoke a mutating alias.
template <class Qualification>
class FullPacket {
public:
    FullPacket(FullPacket&&) noexcept; ///< Transfer ownership, leaving an inactive source.
    FullPacket& operator=(FullPacket&&) noexcept; ///< Replace ownership and invalidate old borrows.
    FullPacket(const FullPacket&) = delete;       ///< No implicit deep copies.
    FullPacket& operator=(const FullPacket&) = delete; ///< No implicit deep copies.
    ~FullPacket(); ///< Release independently of every input, producer and prior packet.
    [[nodiscard]] bool active() const noexcept; ///< Whether this object still owns its packet.
    [[nodiscard]] std::optional<FullPacketView>
    view() const& noexcept; ///< Borrow the exact returned packet; missing for an inactive
                            ///< moved-from object.
    [[nodiscard]] std::optional<FullPacketView>
    view() const&& = delete; ///< Temporary views cannot escape.
    [[nodiscard]] AttemptDiagnostics
    diagnostics() const noexcept; ///< Known counts, unknown QR steps.
    /// Copy measured phase times; missing when recording was disabled.
    [[nodiscard]] std::optional<FullTimings> timings() const noexcept;
    /// Exact recorder clock invocations in this extraction; excludes clocks inside user code.
    [[nodiscard]] std::uint64_t clock_reads() const noexcept;
    /// Copy the strict basis independence bound only when explicitly requested.
    [[nodiscard]] std::optional<double> independence_upper() const noexcept;
    [[nodiscard]] FullPhysicalCoverage physical_coverage() const noexcept; ///< Always unresolved.
    /// Read the exact retained/projection request; null for an inactive moved-from object.
    [[nodiscard]] const linear::ReturnedPrecision* returned_precision() const& noexcept;
    /// Temporaries cannot supply a borrowed precision request.
    [[nodiscard]] const linear::ReturnedPrecision* returned_precision() const&& = delete;
    /// Borrow owned fixed scales/units; missing only for an inactive object.
    [[nodiscard]] std::optional<FullScales> scales() const& noexcept;
    [[nodiscard]] std::optional<FullScales> scales() const&& = delete;
    /// Borrow original packed physical A; real inputs have exactly zero imaginary parts.
    [[nodiscard]] std::span<const std::complex<double>> original_coefficients() const& noexcept;
    [[nodiscard]] std::span<const std::complex<double>> original_coefficients() const&& = delete;
    /// Copy explicit order/result/action capacities, or no value for an inactive object.
    [[nodiscard]] std::optional<FullLimits> limits() const noexcept;
    /// Copy factory admission, which is neither actual allocation nor process RSS.
    [[nodiscard]] std::optional<linear::BackendAdmission> backend_admission() const noexcept;
    /// Copy a derived-root association policy; absent for complete spectra.
    [[nodiscard]] std::optional<FullAssociation> association_policy() const noexcept;
    /// Original complete-packet candidate generator; missing for an inactive object.
    [[nodiscard]] std::optional<FullBackend> source_backend() const noexcept;
    /// Configured QR ceiling 30*n for a complete spectrum, not observed iterations.
    [[nodiscard]] std::optional<core::Index> schur_iteration_limit() const noexcept;
    /// Map each sorted eigenvalue to its actual T diagonal index; empty only for Pairs.
    [[nodiscard]] std::span<const core::Index> diagonal_mapping() const& noexcept;
    [[nodiscard]] std::span<const core::Index> diagonal_mapping() const&& = delete;
    /// Whole right-basis proof; missing for pair packets, never an invented zero certificate.
    [[nodiscard]] const BasisMetrics* basis_proof() const& noexcept;
    [[nodiscard]] const BasisMetrics* basis_proof() const&& = delete;
    /// Independently evaluated left-basis proof, when a cluster requested it.
    [[nodiscard]] const BasisMetrics* left_basis_proof() const& noexcept;
    [[nodiscard]] const BasisMetrics* left_basis_proof() const&& = delete;
    /// Copy one pair proof with checked index; invalid/unrequested indices are failures.
    [[nodiscard]] Result<FullPairProof> pair_proof(core::Index index) const noexcept;
    /// Number of necessary consumer proofs; retained and returned phases are checked separately.
    [[nodiscard]] std::size_t consumer_count() const noexcept;
    /// Copy one achieved proof of the actual returned generation; invalid index is a failure.
    [[nodiscard]] Result<core::ErrorCertificate> consumer_proof(core::Index index) const noexcept;

private:
    friend class detail::FullFactory;
    explicit FullPacket(std::unique_ptr<detail::FullData> data) noexcept;
    std::unique_ptr<detail::FullData> data_;
};

/// Extract and qualify all n roots and complex Schur Q/T from immutable binary64 A.
/// Supports real/complex inputs. Reject unknown/TwoWord coefficients before backend entry.
/// Original generation must equal context.cache.matrix. Physical units must satisfy
/// [A_ij]=[lambda][phi_i]/[phi_j]. Counts and all payload/backend admission are checked
/// before numerical iteration. No input/old packet is changed on any failure.
template <StorageScalar Scalar>
[[nodiscard]] Result<FullSpectrum> full_spectrum(const PrecisionBuffer<Scalar>& original,
                                                 Context context, StateDiagnosticContext state,
                                                 FullScales scales, const FullRequest& request,
                                                 core::ResourcePlan& plan) noexcept;
/// Extract individually qualified pairs associated with the packet-bound selected roots.
[[nodiscard]] Result<FullPairs> full_pairs(const FullSpectrum& packet,
                                           const FullSelection& selection,
                                           const FullPairRequest& request,
                                           core::ResourcePlan& plan) noexcept;
/// Unitary Schur exchange of a packet-bound explicit set, followed by original-system
/// qualification.
[[nodiscard]] Result<FullCluster> full_cluster(const FullSpectrum& packet,
                                               const FullSelection& selection,
                                               const FullClusterRequest& request,
                                               core::ResourcePlan& plan) noexcept;

} // namespace uvlas::numerics::eigen
