#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <uvlas/core/instrumentation.hpp>
#include <uvlas/numerics/eigen_action.hpp>
#include <uvlas/numerics/eigen_metrics.hpp>
#include <variant>

namespace uvlas::numerics::eigen {

/// Right original problem or conjugate-transposed original problem.
enum class ArnoldiSide : std::uint8_t {
    Right, ///< Original eigen equation.
    Left   ///< Conjugate-transposed original eigen equation.
};
/// Request right quantities alone or both independently qualified sides.
enum class ArnoldiSides : std::uint8_t {
    RightOnly, ///< Do not invoke or store left quantities.
    Both       ///< Independently extract and associate both sides.
};
/// Deterministic values are identical with timing enabled or disabled.
enum class ArnoldiTimingMode : std::uint8_t {
    Disabled, ///< Do not read a clock.
    Enabled   ///< Save actual phase durations.
};
/// Physical mode coverage is outside this algebraic extraction's qualification.
struct ArnoldiPhysicalCoverageUnresolved {};
/// No clocks were read for this extraction.
struct ArnoldiUnmeasured {};
/// Inclusive wall-clock phase totals, without a performance acceptance claim.
struct ArnoldiTimings {
    core::StageTiming preparation; ///< Input/result allocation and snapshot copying.
    core::StageTiming krylov; ///< Native search and callbacks, including inverse linear proofs.
    core::StageTiming qualification; ///< Final original-system and consumer checks.
    core::StageTiming subspace;      ///< Small Schur, reordering and lifting, when requested.
};
/// Either explicit absence of measurements or owned actual phase samples.
using ArnoldiMeasurements = std::variant<ArnoldiUnmeasured, ArnoldiTimings>;

/// Fixed physical component scales, units and their exact problem identity.
struct ArnoldiScales {
    MetricScales metric; ///< Positive component D and positive eigenvalue scale.
    std::span<const linear::UnitDimensions> component_units; ///< One SI exponent vector per row.
    linear::UnitDimensions eigenvalue_units; ///< SI exponents of the physical eigenvalue.
    Context context; ///< Must match the frozen original provider in every field.
};
/// Explicit association policy, distinct from a forward eigenvalue-error guarantee.
struct ArnoldiAssociation {
    double absolute; ///< Nonnegative finite distance allowance in dimensionless eigenvalues.
    double relative; ///< Nonnegative finite multiplier of max(abs(mu1),abs(mu2)).
};
/// Pair qualification uses the same requirements independently on each requested side.
struct ArnoldiPairRequest {
    PairMetricRequest accuracy; ///< Original norm residual, normalization and optional beta.
};
/// Invariant-space qualification; it never substitutes individual Ritz-vector residuals.
struct ArnoldiSubspaceRequest {
    BasisMetricRequest accuracy; ///< Original matrix residual, full Gram and optional beta.
};
/// Only metrics belonging to the requested output are supplied.
using ArnoldiQualification = std::variant<ArnoldiPairRequest, ArnoldiSubspaceRequest>;
/// Extract the largest magnitudes of the original dimensionless operator.
struct ArnoldiLargestMagnitude {};

/// Synchronous access to actual shifted inverse values and their linear qualification.
/// Each provider evaluates a complex physical inverse, including for a real original A.
/// Right is (A-sigma I)^-1; Left is (A^*-conj(sigma) I)^-1. Each side's frozen stamp
/// names its actual factors, side policy and physical shift. No provider is retained.
/// absolute_error_upper bounds the exact inverse output by producer contract; a small
/// linear residual alone is never reinterpreted as that forward-error bound.
class ArnoldiInverseAccess {
public:
    /// Release through the inverse interface.
    virtual ~ArnoldiInverseAccess() = default;
    /// Borrow the frozen physical inverse provider for the requested side.
    [[nodiscard]] virtual const ActionProvider<std::complex<double>>&
    provider(ArnoldiSide side) const = 0;
    /// Assess the exact actual physical RHS and inverse words from this call.
    /// Match expected context and budget bitwise. Called synchronously after a successful
    /// inverse action; exceptions and unachieved/mismatched certificates reject extraction.
    [[nodiscard]] virtual core::Result<core::ErrorCertificate>
    qualify(ArnoldiSide side, WordView<std::complex<double>> rhs,
            WordView<std::complex<double>> solution, std::span<const double> inverse_error_upper,
            const ActionStamp& actual_stamp, MetricRequirement requirement,
            core::CertificateContext expected) const = 0;
};
/// Complex physical shift with mandatory inverse-output bounds and linear certificates.
struct ArnoldiShiftInvert {
    std::complex<double> shift; ///< Physical sigma; the left physical inverse uses its conjugate.
    const ArnoldiInverseAccess& inverse; ///< Synchronous inverse provider/evaluator lifetime.
    MetricRequirement linear_accuracy;   ///< Required shifted original-system linear proof.
    core::NumericalGeneration
        first_inverse_generation; ///< First caller-reserved inverse-candidate ID.
    // The caller reserves limits.actions.inverse consecutive unique IDs, disjoint
    // from other numerical fields. Right/left share that range and every started
    // inverse burns one ordinal, even on failure. This is not a factor/shift stamp.
};
/// No implicit spectral transform, preconditioning or small-order fallback is allowed.
using ArnoldiSelection = std::variant<ArnoldiLargestMagnitude, ArnoldiShiftInvert>;

/// Read-only packed logical matrix; capacity outside values is never exposed or qualified.
/// The spectrum owns its storage. Moving or destroying that owner invalidates this borrow.
struct ArnoldiMatrixView {
    DenseShape shape; ///< Logical dimensions and packed leading dimension.
    std::span<const std::complex<double>> values; ///< Exactly rows*columns actual entries.
};

/// Read-only eigenvector matrix in physical coordinates; each column belongs to one value.
struct ArnoldiPairView {
    ArnoldiMatrixView vectors; ///< Packed n-by-k actual returned physical vectors.
};
/// Actual dimensionless invariant basis and its k-by-k triangular restriction.
struct ArnoldiSubspaceView {
    ArnoldiMatrixView basis;      ///< n-by-k orthonormal dimensionless columns.
    ArnoldiMatrixView triangular; ///< k-by-k upper triangular restriction of the original action.
    std::span<const core::Index>
        diagonal_columns; ///< Sorted public value index to this T diagonal.
};
/// Output kind remains explicit even for an empty range.
using ArnoldiVectorView = std::variant<ArnoldiPairView, ArnoldiSubspaceView>;
/// Unrequested left quantities are not stored as empty achieved arrays.
using ArnoldiLeftView = std::variant<NotRequested, ArnoldiPairView, ArnoldiSubspaceView>;
/// A synchronous borrow of actual retained or returned quantities; no scratch owner escapes.
struct ArnoldiCandidateView {
    ArnoldiMatrixView eigenvalues;        ///< k-by-1 physical original eigenvalues.
    ArnoldiVectorView right;              ///< Actual right vectors or invariant space.
    ArnoldiLeftView left;                 ///< Actual left quantities if requested.
    Context context;                      ///< Original problem and layout identity.
    StateDiagnosticContext state;         ///< Explicit report-only accepted-state identities.
    core::NumericalGeneration generation; ///< Identity of these exact candidate values.
    core::Precision representation; ///< Present implementation returns explicitly Binary64 values.
};
/// Synchronous physical/constraint/reaction evaluator of the actual result being consumed.
class ArnoldiConsumerEvaluator {
public:
    /// Release through the consumer interface.
    virtual ~ArnoldiConsumerEvaluator() = default;
    /// Declare the numerical representations actually consumed.
    [[nodiscard]] virtual SupportedPrecision supported_precision() const noexcept = 0;
    /// Re-evaluate this requirement on the exact retained and returned candidates.
    /// A certificate for another observation, scope, budget, generation or precision rejects.
    [[nodiscard]] virtual core::Result<core::ErrorCertificate>
    assess(const ArnoldiCandidateView& retained, const ArnoldiCandidateView& candidate,
           const linear::ConsumerCheck& check, core::CertificateContext expected) const = 0;
};

/// Complete explicit extraction request. Input arrays remain alive and unchanged until
/// the preparation snapshot finishes. Request values, D, units, seed and checks are copied.
/// The provider/evaluator remain alive, and the provider frozen, for the whole call.
template <StorageScalar Scalar>
struct ArnoldiRequest {
    core::Index nev; ///< Minimum requested count; real conjugate completion can return nev+1.
    core::Index ncv; ///< Explicit Krylov width; must satisfy nev+2<=ncv<=n, with no clipping.
    WordView<Scalar> initial; ///< Explicit finite nonzero Binary64 native-coordinate seed.
    double backend_tolerance; ///< Positive finite candidate-extraction tolerance, not acceptance.
    std::uint64_t
        backend_iterations;     ///< Upper limit passed to Spectra, independent of its counter.
    ArnoldiSelection selection; ///< Original-largest-magnitude or an explicit complex shift.
    ArnoldiQualification qualification; ///< Pair or invariant-space requirements.
    ArnoldiSides sides; ///< Left and right are qualified independently when Both is requested.
    ArnoldiAssociation association; ///< Fixed matching policy for sides and selected groups.
    linear::ReturnedPrecision
        returned;                 ///< Retained or explicit projection with distinct generations.
    StateDiagnosticContext state; ///< Optional accepted parent/last accepted reporting identities.
    std::span<const linear::ConsumerCheck> required_checks; ///< Every necessary consumer check.
    const ArnoldiConsumerEvaluator& evaluator;              ///< Called only for supplied checks.
    ArnoldiTimingMode timing;                               ///< Disabled performs no clock reads.
};
/// Explicit allocation, search and action limits; no unlimited defaults are supplied.
struct ArnoldiLimits {
    core::Index maximum_order;    ///< Maximum original square order.
    core::Index maximum_returned; ///< Capacity including real conjugate completion.
    core::Index maximum_ncv;      ///< Maximum active Krylov storage width.
    std::uint64_t
        maximum_backend_iterations; ///< Capacity of the requested backend iteration limit.
    ActionLimits
        actions; ///< All extraction retries and original/adjoint/inverse calls share limits.
    linear::BackendMemoryLimit
        backend_memory; ///< Estimate admission or genuinely proven-bound requirement.
};
/// The pinned backend does not expose an actual restart-count getter.
struct ArnoldiRestartsNotExposed {};
/// Fixed owned extraction counters and explicitly unexposed backend information.
struct ArnoldiDiagnostics {
    AttemptDiagnostics attempts; ///< Counts across both sides and final qualification.
    ArnoldiRestartsNotExposed actual_restarts; ///< Not a measured zero or a hard restart cap.
    core::Index requested_nev;                 ///< Caller request, never silently clipped.
    core::Index requested_ncv;                 ///< Caller request, never silently clipped.
    core::Index returned_count; ///< Actual count including any qualified conjugate completion.
    core::Index right_krylov_dimension; ///< Actual factorization dimension, not allocated ncv.
    std::optional<core::Index>
        left_krylov_dimension; ///< Absent when the left side was not requested.
    std::uint64_t
        backend_reported_operations; ///< Backend counter, including internal adapter calls.
    double initial_norm;             ///< Actual native-coordinate Euclidean seed norm.
    std::array<std::byte, 32>
        initial_sha256; ///< Exact seed bits and representation, with domain/shape.
};

/// This packet qualifies the explicitly requested eigenvalue subset, not full coverage.
struct ArnoldiRequestedSubsetQualified {};
/// This packet qualifies the explicitly requested invariant space, not an eigenbasis.
struct ArnoldiRequestedSubspaceQualified {};
/// Explicit algebraic qualification carried independently of matrix views.
using ArnoldiQualificationState =
    std::variant<ArnoldiRequestedSubsetQualified, ArnoldiRequestedSubspaceQualified>;
/// A finite approximate sensitivity from actual normalized left/right vectors.
struct ArnoldiConditionEstimate {
    double value; ///< Diagnostic estimate, not a forward-error bound.
};
/// A positive finite sensitivity denominator could not be established.
struct ArnoldiConditionUnresolved {};
/// Sensitivity is absent when not requested, unresolved, or explicitly an estimate.
using ArnoldiSensitivity =
    std::variant<NotRequested, ArnoldiConditionUnresolved, ArnoldiConditionEstimate>;
/// Original scalar domain, distinct from a complex transformed native search.
enum class ArnoldiDomain : std::uint8_t {
    Real,   ///< Real original or native scalars.
    Complex ///< Complex original or native scalars.
};
/// Owned metadata of the actual pinned numerical backends.
struct ArnoldiBackend {
    std::array<unsigned, 3> spectra_version; ///< Actual Spectra major/minor/patch.
    std::array<unsigned, 3> eigen_version;   ///< Actual Eigen major/minor/patch.
    ArnoldiDomain original_domain;           ///< Scalar domain of the original provider.
    ArnoldiDomain native_domain;             ///< Native real or complex candidate search.
};
/// Owned shift policy/stamps, without retaining the inverse provider/evaluator.
struct ArnoldiShiftMetadata {
    std::complex<double> physical_shift; ///< Original sigma, not its dimensionless value.
    ActionStamp right_inverse;           ///< Actual right inverse factors/policy identity.
    std::variant<NotRequested, ActionStamp>
        left_inverse;                  ///< Actual left inverse identity, if requested.
    MetricRequirement linear_accuracy; ///< Required inverse linear qualification.
    core::NumericalGeneration
        first_inverse_generation;           ///< Caller-reserved first inverse-candidate ID.
    std::uint64_t inverse_generation_count; ///< Size of the validated reserved range.
};
/// An owned original or shift-near extraction scope.
using ArnoldiSelectionMetadata = std::variant<ArnoldiLargestMagnitude, ArnoldiShiftMetadata>;
/// Complete explicit native search controls and qualification request.
struct ArnoldiSettings {
    core::Index nev;                       ///< Minimum original requested count.
    core::Index ncv;                       ///< Explicit native Krylov width.
    double backend_tolerance;              ///< Candidate tolerance actually passed to Spectra.
    std::uint64_t backend_iteration_limit; ///< Limit actually passed to Spectra.
    ArnoldiSides sides;                    ///< Requested original sides.
    ArnoldiQualification qualification;    ///< Original mandatory and optional metric budgets.
    ArnoldiTimingMode timing;              ///< Explicit measurement policy.
};
/// Synchronous metadata view. All referenced arrays are owned by the spectrum.
struct ArnoldiMetadataView {
    Context context; ///< Original complete problem identity.
    ActionStamp
        original_stamp; ///< Complete frozen original action, including its transform policy.
    StateDiagnosticContext state;       ///< Report-only accepted-state identities.
    std::span<const double> components; ///< Original fixed physical coordinate scales.
    std::span<const linear::UnitDimensions> component_units; ///< Original component SI exponents.
    double eigenvalue_scale;                 ///< Original physical eigenvalue scale.
    linear::UnitDimensions eigenvalue_units; ///< Original eigenvalue SI exponents.
    ArnoldiMatrixView
        initial_vector; ///< Exact supplied seed, including zero imaginary lanes of real data.
    ArnoldiAssociation association;     ///< Fixed dimensionless association policy.
    ArnoldiSelectionMetadata selection; ///< Actual transform/region and inverse stamps.
    ArnoldiBackend backend;             ///< Actual fixed numerical versions and domains.
    linear::ReturnedPrecision returned; ///< Actual retention/projection generations.
    ArnoldiSettings settings;           ///< Actual search and qualification controls.
    ArnoldiLimits limits;               ///< Every explicit resource and action limit.
    std::size_t required_check_count;   ///< All necessary consumers are owned; read each by value.
};

/// One actual shifted-inverse call's independently owned numerical guarantees.
struct ArnoldiInverseAssessment {
    ArnoldiSide side;                    ///< Physical right inverse or its shifted adjoint inverse.
    core::NumericalGeneration candidate; ///< Unique actual output words of this call.
    double
        native_output_error_upper; ///< Conservative infinity bound after fixed D/s_lambda scaling.
    core::ErrorCertificate linear; ///< Separate shifted-original linear residual qualification.
};

namespace detail {
struct ArnoldiData;
struct ArnoldiWorkspaceData;
class ArnoldiWorkspaceAccess;
} // namespace detail
/// Caller-owned, nonreentrant provider I/O and qualification scratch.
/// Prepare once for an explicit order/capacity. All provider metadata and actions
/// run under the same exclusive lease; separate workspaces remain independent.
class ArnoldiWorkspace {
public:
    /// Allocate actual real/complex action planes and bounded qualification scratch.
    /// Negative dimensions reject; checked products and plan limits apply before allocation.
    [[nodiscard]] static core::Result<ArnoldiWorkspace>
    create(core::Index order, core::Index maximum_returned, core::ResourcePlan& plan) noexcept;
    /// Workspace identity and active leases cannot be copied.
    ArnoldiWorkspace(const ArnoldiWorkspace&) = delete;
    /// Workspace identity cannot be reassigned.
    ArnoldiWorkspace& operator=(const ArnoldiWorkspace&) = delete;
    /// Transfer inactive storage; moving while extraction is active is forbidden.
    ArnoldiWorkspace(ArnoldiWorkspace&&) noexcept;
    /// Transfer inactive storage; existing borrows expire. Neither workspace may be active.
    ArnoldiWorkspace& operator=(ArnoldiWorkspace&&) noexcept;
    /// Destroy after every synchronous extraction using it has returned.
    ~ArnoldiWorkspace();
    /// Prepared original order, independent of provider metadata calls.
    [[nodiscard]] core::Index order() const noexcept;
    /// Maximum logical result count for this workspace's real/complex planes.
    [[nodiscard]] core::Index capacity() const noexcept;

private:
    friend class detail::ArnoldiWorkspaceAccess;
    explicit ArnoldiWorkspace(std::unique_ptr<detail::ArnoldiWorkspaceData> data) noexcept;
    std::unique_ptr<detail::ArnoldiWorkspaceData> data_;
};

/// Owned qualified requested subset or subspace of one frozen algebraic operator.
/// This object cannot qualify physical stability, flutter coverage or a minimum critical speed.
/// Moved-from owners yield absent optional values/empty spans, never invalid references.
class ArnoldiSpectrum {
public:
    /// Transfer all owned arrays and proofs, leaving the source inactive.
    ArnoldiSpectrum(ArnoldiSpectrum&&) noexcept;
    /// Replace this packet by ownership transfer; existing borrows expire.
    ArnoldiSpectrum& operator=(ArnoldiSpectrum&&) noexcept;
    /// Implicit copies of owned arrays are forbidden.
    ArnoldiSpectrum(const ArnoldiSpectrum&) = delete;
    /// Copy assignment is forbidden.
    ArnoldiSpectrum& operator=(const ArnoldiSpectrum&) = delete;
    /// Release all arrays independently of providers and workspace.
    ~ArnoldiSpectrum();
    /// Whether this object retains ownership after moving.
    [[nodiscard]] bool active() const noexcept;
    /// Borrow actual consumed values; the owner must outlive this view.
    [[nodiscard]] std::optional<ArnoldiCandidateView> values() const& noexcept;
    std::optional<ArnoldiCandidateView> values() const&& = delete;
    /// Borrow actual retained pre-projection values and their original generation.
    [[nodiscard]] std::optional<ArnoldiCandidateView> retained_values() const& noexcept;
    std::optional<ArnoldiCandidateView> retained_values() const&& = delete;
    /// Borrow owned scales, seed, units and actual selection/backend metadata.
    [[nodiscard]] std::optional<ArnoldiMetadataView> metadata() const& noexcept;
    std::optional<ArnoldiMetadataView> metadata() const&& = delete;
    /// Explicit qualification, or no value after moving from this owner.
    [[nodiscard]] std::optional<ArnoldiQualificationState> qualification() const noexcept;
    /// Copy one actual returned pair proof; absent for the other output kind or invalid index/side.
    [[nodiscard]] std::optional<PairMetrics> pair_metric(ArnoldiSide side,
                                                         core::Index index) const noexcept;
    /// Copy an actual invariant-space proof, absent when that output/side was not requested.
    [[nodiscard]] std::optional<BasisMetrics> subspace_metric(ArnoldiSide side) const noexcept;
    /// Borrow per-pair sensitivity diagnostics. Empty for a subspace or moved-from owner.
    [[nodiscard]] std::span<const ArnoldiSensitivity> sensitivities() const& noexcept;
    std::span<const ArnoldiSensitivity> sensitivities() const&& = delete;
    /// Original problem identity, retained independently of callbacks and workspace.
    [[nodiscard]] std::optional<Context> context() const noexcept;
    /// Explicit optional parent/last accepted identifiers, not an ancestry proof.
    [[nodiscard]] std::optional<StateDiagnosticContext> state() const noexcept;
    /// Requested retention/projection policy with both generation identifiers.
    [[nodiscard]] std::optional<linear::ReturnedPrecision> returned_precision() const noexcept;
    /// Copy all extraction counts and the exact initial-vector fingerprint.
    [[nodiscard]] std::optional<ArnoldiDiagnostics> diagnostics() const noexcept;
    /// Algebraic extraction never constructs qualified physical coverage.
    [[nodiscard]] ArnoldiPhysicalCoverageUnresolved physical_coverage() const noexcept;
    /// Copy measured phase samples, or explicit disabled status.
    [[nodiscard]] std::optional<ArnoldiMeasurements> measurements() const noexcept;
    /// Declared backend admission; this is not allocator interception or process RSS.
    [[nodiscard]] std::optional<linear::BackendAdmission> backend_admission() const noexcept;
    /// Number of successful shifted-inverse linear proofs, including both requested sides.
    [[nodiscard]] std::size_t inverse_certificate_count() const noexcept;
    /// Copy the scaled forward-output bound and independent linear proof for one inverse call.
    [[nodiscard]] std::optional<ArnoldiInverseAssessment>
    inverse_assessment(core::Index index) const noexcept;
    /// Copy one proof in actual inverse-callback order, absent for an invalid index.
    [[nodiscard]] std::optional<core::ErrorCertificate>
    inverse_certificate(core::Index index) const noexcept;
    /// Copy an owned necessary consumer request, absent for an invalid index.
    [[nodiscard]] std::optional<linear::ConsumerCheck>
    required_check(core::Index index) const noexcept;
    /// Number of independently evaluated necessary consumer requirements.
    [[nodiscard]] std::size_t consumer_certificate_count() const noexcept;
    /// Copy one actual returned consumer proof, absent for an invalid index.
    [[nodiscard]] std::optional<core::ErrorCertificate>
    consumer_certificate(core::Index index) const noexcept;

private:
    template <StorageScalar Scalar>
    friend Result<ArnoldiSpectrum>
    arnoldi_spectrum(const ActionProvider<Scalar>&, const ArnoldiScales&,
                     const ArnoldiRequest<Scalar>&, const ArnoldiLimits&, ArnoldiWorkspace&,
                     core::ResourcePlan&) noexcept;
    explicit ArnoldiSpectrum(std::unique_ptr<detail::ArnoldiData> data) noexcept;
    std::unique_ptr<detail::ArnoldiData> data_;
};

/// Extract and qualify actual original eigenpairs or an invariant space using pinned Spectra.
/// No original dense n-by-n matrix is assembled. Callbacks and all input arrays remain borrowed
/// only until return. Prepare workspace through its explicit plan; owning results use result_plan.
/// A busy workspace rejects before accessing provider metadata or allocating result arrays.
/// Binary64 candidate generation is explicit. Unsupported input precision rejects before work;
/// unachieved original/consumer/inverse proofs never become successful partial results.
/// All exceptions become owned EigenFailure, with original inner causes preserved.
template <StorageScalar Scalar>
[[nodiscard]] Result<ArnoldiSpectrum>
arnoldi_spectrum(const ActionProvider<Scalar>& provider, const ArnoldiScales& scales,
                 const ArnoldiRequest<Scalar>& request, const ArnoldiLimits& limits,
                 ArnoldiWorkspace& workspace, core::ResourcePlan& result_plan) noexcept;

} // namespace uvlas::numerics::eigen
