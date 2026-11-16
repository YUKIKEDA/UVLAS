#pragma once

#include <span>
#include <uvlas/numerics/eigen_common.hpp>

namespace uvlas::numerics::eigen {

/// Borrow an original action's numerical words and componentwise absolute error bounds.
/// Both spans remain alive and unchanged for the synchronous metric evaluation.
template <StorageScalar Scalar>
struct ActionValues {
    WordView<Scalar> values; ///< Signed action words, without dropping their trailing plane.
    std::span<const double> absolute_error_upper; ///< Finite nonnegative bound for every entry.
};

/// Fixed component and eigenvalue scales of the original physical problem.
struct MetricScales {
    std::span<const double> components; ///< Positive finite diagonal entries of D.
    double eigenvalue_scale; ///< Positive finite eigenvalue scale, in the original units.
};

/// Identity assigned to the actual candidate whose metrics are being evaluated.
struct MetricContext {
    Context system;               ///< Original problem, layouts, policies and evaluation time.
    StateDiagnosticContext state; ///< Explicit reporting identities; never inferred.
    core::NumericalGeneration values_generation; ///< Generation of the exact candidate words.
};

/// Which original eigen equation the supplied action evaluates.
enum class PairSide : std::uint8_t {
    Right, ///< A*phi-lambda*phi, scaled by D^-1 and the eigenvalue scale.
    Left   ///< A^*psi-conj(lambda)*psi, scaled by D and the eigenvalue scale.
};

/// Independent budgets for one eigenpair's residual, vector norm and optional backward error.
struct PairMetricRequest {
    MetricRequirement
        residual; ///< Original dimensionless residual divided by vector infinity norm.
    MetricRequirement normalization; ///< Absolute error from unit dimensionless Euclidean norm.
    BackwardRequirement
        backward; ///< Additional componentwise criterion, or explicitly unrequested.
};

/// Independent budgets for a Schur/invariant or mass-orthonormal basis.
struct BasisMetricRequest {
    MetricRequirement residual; ///< Original matrix residual or maximum generalized-column metric.
    MetricRequirement orthogonality; ///< Infinity norm of the complete Gram error matrix.
    BackwardRequirement
        backward; ///< Additional componentwise criterion, or explicitly unrequested.
};

/// Independent budgets for a mass-normalized generalized full basis.
struct GeneralizedMetricRequest {
    MetricRequirement residual;      ///< Maximum original residual over all generalized columns.
    MetricRequirement normalization; ///< Maximum diagonal mass-Gram error from one.
    MetricRequirement orthogonality; ///< Infinity norm of the complete mass-Gram error.
    BackwardRequirement
        backward; ///< Additional componentwise criterion, or explicitly unrequested.
};

/// An evaluated optional criterion; NotRequested is never an achieved zero.
using BackwardCertificate = std::variant<NotRequested, core::ErrorCertificate>;

/// Fixed-size owned proofs for one pair; every present certificate is achieved.
struct PairMetrics {
    core::ErrorCertificate residual;      ///< Proof for the actual original-system residual.
    core::ErrorCertificate normalization; ///< Proof for the actual dimensionless vector norm.
    BackwardCertificate backward;         ///< Proof of the additional requirement, when requested.
    double
        vector_infinity_lower; ///< Positive certified lower bound used in residual normalization.
};

/// Fixed-size owned proofs for a complete evaluated basis, with no scratch borrow.
struct BasisMetrics {
    core::ErrorCertificate residual; ///< Proof covering all matrix entries or generalized columns.
    core::ErrorCertificate
        orthogonality;            ///< Complete Gram-error proof; nonempty bases also need <1.
    BackwardCertificate backward; ///< Optional componentwise proof, never substituted for residual.
};

/// Separate diagonal normalization and full mass-orthogonality proofs.
struct GeneralizedMetrics {
    core::ErrorCertificate residual;      ///< Proof covering every original generalized residual.
    core::ErrorCertificate normalization; ///< Worst diagonal mass-Gram error.
    core::ErrorCertificate orthogonality; ///< Full mass-Gram error, additionally <1 when nonempty.
    BackwardCertificate backward;         ///< Optional additional componentwise proof.
};

/// Qualify one nonzero actual physical eigenvector and a one-entry eigenvalue view.
/// action is A*vector for Right or A^*vector for Left, in physical coordinates.
/// Real and complex inputs retain both words. No allocation or mutation occurs.
/// Invalid shapes/negative inputs return InvalidInput; nonfinite supplied values return Nonfinite.
/// Unrepresentable arithmetic or unmet budgets return PrecisionUnresolved.
template <StorageScalar Scalar>
[[nodiscard]] Result<PairMetrics>
assess_pair_metrics(WordView<Scalar> vector, WordView<Scalar> eigenvalue,
                    ActionValues<Scalar> action, MetricScales scales, PairSide side,
                    const PairMetricRequest& request, const MetricContext& context) noexcept;

/// Qualify packed column-major dimensionless Q(n,k), T(k,k), and original Ahat*Q(n,k).
/// The triangular component errors include any physical-eigenvalue/scale conversion enclosure.
/// Require 0<=k<=n, exact upper-triangular T, and full Gram/residual checks.
/// A nonempty basis must have Gram error <1 as well as meeting the requested budget.
/// Empty bases have zero metrics and need no positive norm. No allocation or mutation occurs.
template <StorageScalar Scalar>
[[nodiscard]] Result<BasisMetrics>
assess_schur_metrics(core::Index rows, core::Index columns, WordView<Scalar> basis,
                     ActionValues<Scalar> triangular, ActionValues<Scalar> action,
                     const BasisMetricRequest& request, const MetricContext& context) noexcept;

/// Qualify a packed physical real full basis Phi(n,n), n eigenvalues and original K*Phi/M*Phi.
/// mass_scale is the positive finite scale of phi^T M phi, not a mass-only unit.
/// K/M symmetry and positive-definiteness qualification belong to the producer.
/// Evaluate original residuals and Phi^T M Phi/mass_scale-I without forming rounded scaled K/M.
/// Empty bases have zero metrics. Both planes are consumed; no allocation or mutation occurs.
[[nodiscard]] Result<GeneralizedMetrics>
assess_generalized_metrics(core::Index rows, WordView<double> basis, WordView<double> eigenvalues,
                           ActionValues<double> stiffness_action, ActionValues<double> mass_action,
                           MetricScales scales, double mass_scale,
                           const GeneralizedMetricRequest& request,
                           const MetricContext& context) noexcept;

} // namespace uvlas::numerics::eigen
