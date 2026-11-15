#pragma once

#include <span>
#include <utility>
#include <uvlas/core/error_certificate.hpp>
#include <uvlas/numerics/storage.hpp>
#include <variant>

namespace uvlas::numerics {

/// Finite high and low parts whose sum denotes one real value; no implicit scalar conversion.
struct TwoWordValue {
    double hi; ///< Rounded leading part.
    double lo; ///< Retained trailing part.
};

/// Return an error-free sum, or Nonfinite if inputs or the resulting sum are nonfinite.
/// Requires round-to-nearest and gradual underflow (FTZ/DAZ disabled); an unsupported
/// floating-point environment returns PrecisionUnsupported.
[[nodiscard]] core::Result<TwoWordValue> two_sum(double a, double b) noexcept;
/// Return an error-free product using an explicit fused multiply-add.
/// Nonfinite rejects nonfinite input or product overflow. PrecisionUnresolved rejects
/// a nonzero exact product requiring bits below the smallest binary64 subnormal.
/// Requires round-to-nearest and gradual underflow (FTZ/DAZ disabled); no unsupported value is
/// silently rounded into success.
[[nodiscard]] core::Result<TwoWordValue> two_product(double a, double b) noexcept;

/// Owning binary64 representation. Borrowed arrays expire with their matrix owner.
template <StorageScalar Scalar>
struct Binary64 {
    DenseMatrix<Scalar> values; ///< Column-major finite values, including finite padding.
};
/// Owning high/low representation; both matrices have exactly the same storage shape.
template <StorageScalar Scalar>
struct TwoWord {
    DenseMatrix<Scalar> hi; ///< Leading values.
    DenseMatrix<Scalar> lo; ///< Trailing values; never implicitly discarded.
};

template <StorageScalar Scalar>
class QualifiedProjection;
template <StorageScalar Scalar>
class ProjectionVerifier;

/// Move-only numerical values with an explicit arithmetic representation and generation.
/// Access is read-only; producing changed values requires a new owning buffer.
template <StorageScalar Scalar>
class PrecisionBuffer {
public:
    /// Preserve binary64 values by moving their ownership; validate before accepting them.
    [[nodiscard]] static core::Result<PrecisionBuffer>
    binary64(DenseMatrix<Scalar> values, core::NumericalGeneration generation) noexcept;
    /// Preserve both words by moving ownership; reject unequal shapes or nonfinite planes.
    [[nodiscard]] static core::Result<PrecisionBuffer>
    two_word(DenseMatrix<Scalar> hi, DenseMatrix<Scalar> lo,
             core::NumericalGeneration generation) noexcept;
    /// Transfer ownership; views into the source are invalidated.
    PrecisionBuffer(PrecisionBuffer&&) noexcept = default;
    /// Replace owned values, invalidating previous views.
    PrecisionBuffer& operator=(PrecisionBuffer&&) noexcept = default;
    /// Copying is forbidden; explicit factory copies belong to the storage layer.
    PrecisionBuffer(const PrecisionBuffer&) = delete;
    /// Copy assignment is forbidden.
    PrecisionBuffer& operator=(const PrecisionBuffer&) = delete;
    /// Release owned values.
    ~PrecisionBuffer() = default;
    /// Actual retained representation, not a requested arithmetic policy.
    [[nodiscard]] core::Precision precision() const noexcept;
    /// Generation of the exact retained values.
    [[nodiscard]] core::NumericalGeneration generation() const noexcept;
    /// Dimensions and padding shared by both planes.
    [[nodiscard]] DenseShape shape() const noexcept;
    /// Read the tagged owning representation while this buffer is alive and unmoved.
    [[nodiscard]] const std::variant<Binary64<Scalar>, TwoWord<Scalar>>& storage() const& noexcept;
    /// Temporaries cannot supply borrowed representations.
    [[nodiscard]] const std::variant<Binary64<Scalar>, TwoWord<Scalar>>& storage() const&& = delete;

private:
    template <StorageScalar T>
    friend core::Result<QualifiedProjection<T>>
    project_binary64(const PrecisionBuffer<T>&, const DenseMatrix<T>&, const DenseMatrix<T>&,
                     std::span<const double>, double, core::NumericalGeneration, core::ScopeId,
                     const ProjectionVerifier<T>&, core::ResourcePlan&) noexcept;
    PrecisionBuffer(std::variant<Binary64<Scalar>, TwoWord<Scalar>> storage,
                    core::NumericalGeneration generation) noexcept;
    std::variant<Binary64<Scalar>, TwoWord<Scalar>> storage_;
    core::NumericalGeneration generation_;
};

/// Explicit capability of a numerical consumer; this does not implement a kernel.
enum class SupportedPrecision : std::uint8_t {
    Binary64Only,      ///< Consumer has only a binary64 path.
    Binary64AndTwoWord ///< Consumer evaluates both planes without dropping the low word.
};
/// Reject unsupported representation or unknown capability with PrecisionUnsupported.
[[nodiscard]] core::Result<void> require_precision(core::Precision actual,
                                                   SupportedPrecision supported) noexcept;

/// Add a finite binary64 correction to both retained words and return new owning values.
/// The source stays unchanged; generation must differ from the source generation.
/// Allocation is budgeted by plan. Shape mismatch returns InvalidInput.
template <StorageScalar Scalar>
[[nodiscard]] core::Result<PrecisionBuffer<Scalar>>
add_correction(const PrecisionBuffer<Scalar>& source, const DenseMatrix<Scalar>& correction,
               core::NumericalGeneration generation, core::ResourcePlan& plan) noexcept;

/// Evaluate b - A*x with error-free real products and compensated accumulation.
/// Real and imaginary cross terms and both retained planes are consumed in fixed order.
/// Logical dimensions are A(m,n), x(n,k), b(m,k); padding is validated but never summed.
/// Return owning binary64 residuals, or a shape, finite-value, precision, or resource failure.
/// This evaluates the stored operator and does not certify its physical discretization.
template <StorageScalar Scalar>
[[nodiscard]] core::Result<DenseMatrix<Scalar>>
original_residual(const DenseMatrix<Scalar>& a, const PrecisionBuffer<Scalar>& x,
                  const DenseMatrix<Scalar>& b, core::ResourcePlan& plan) noexcept;

/// Identities of one source and its candidate binary64 projection for one consumer.
struct ProjectionContext {
    core::NumericalGeneration source; ///< Original retained-value generation.
    core::NumericalGeneration target; ///< Distinct generation assigned to the projected values.
    core::ScopeId consumer;           ///< Consumer whose complete necessary checks are evaluated.
};
/// Whether a consumer evaluated all necessary observation and compatibility checks.
enum class ProjectionAssessment : std::uint8_t {
    Verified,   ///< Every necessary consumer and observation check was evaluated and achieved.
    Unevaluated ///< One or more necessary checks have no evaluation.
};
/// Synchronous consumer-specific projection checker; no precomputed certificate is accepted.
/// The implementation must evaluate these exact source and candidate values, including
/// every necessary observation and the target consumer's binary64 compatibility.
/// The caller owns this checker and keeps it alive for the projection call.
template <StorageScalar Scalar>
class ProjectionVerifier {
public:
    /// Release through the interface.
    virtual ~ProjectionVerifier() = default;
    /// Inspect the exact candidate; borrowed references expire when this call returns.
    /// Return Verified only after all necessary checks meet their caller-defined budgets.
    /// Return Unevaluated or PrecisionUnresolved when any required check is unavailable;
    /// return PrecisionUnsupported for a consumer without a binary64 path. Throwing is
    /// translated into Failure by the projection entry; the source is never modified.
    [[nodiscard]] virtual core::Result<ProjectionAssessment>
    verify(const PrecisionBuffer<Scalar>& source, const DenseMatrix<Scalar>& candidate,
           const ProjectionContext& context) const = 0;
};

/// Owns the exact immutable candidate evaluated in the projection call.
/// Generation metadata labels the evaluation; it is not a global content hash.
template <StorageScalar Scalar>
class QualifiedProjection {
public:
    /// Transfer ownership and invalidate source views.
    QualifiedProjection(QualifiedProjection&&) noexcept = default;
    /// Replace owned result and invalidate destination views.
    QualifiedProjection& operator=(QualifiedProjection&&) noexcept = default;
    /// Copying is forbidden.
    QualifiedProjection(const QualifiedProjection&) = delete;
    /// Copy assignment is forbidden.
    QualifiedProjection& operator=(const QualifiedProjection&) = delete;
    /// Release owned result.
    ~QualifiedProjection() = default;
    /// Borrow projected values until move, assignment, or destruction.
    [[nodiscard]] const DenseMatrix<Scalar>& values() const& noexcept {
        return values_;
    }
    /// Temporaries cannot supply borrowed values.
    [[nodiscard]] const DenseMatrix<Scalar>& values() const&& = delete;
    /// Exact generations and consumer evaluated for this result.
    [[nodiscard]] ProjectionContext context() const noexcept {
        return context_;
    }

private:
    template <StorageScalar T>
    friend core::Result<QualifiedProjection<T>>
    project_binary64(const PrecisionBuffer<T>&, const DenseMatrix<T>&, const DenseMatrix<T>&,
                     std::span<const double>, double, core::NumericalGeneration, core::ScopeId,
                     const ProjectionVerifier<T>&, core::ResourcePlan&) noexcept;
    QualifiedProjection(DenseMatrix<Scalar> values, ProjectionContext context) noexcept
        : values_(std::move(values)), context_(context) {}
    DenseMatrix<Scalar> values_;
    ProjectionContext context_;
};

/// Explicitly round hi+lo, then reevaluate the original residual and consumer checks.
/// Scales are positive finite physical residual scales, one per logical b element in
/// packed column order, fixed independently of the computed solution. tolerance is
/// a nonnegative finite bound on max(abs(residual)/scale). Qualification encloses all
/// discarded compensation and final-rounding errors. A real-plus-imaginary magnitude
/// bound is a conservative sufficient condition for complex values. An insufficient
/// bound returns PrecisionUnresolved even when diagnostic residuals round to zero.
/// Zero-dimensional systems
/// are valid and still require the consumer checker. The target generation must differ.
/// A residual or observation miss returns PrecisionUnresolved; unsupported consumption
/// returns PrecisionUnsupported. No qualified result exists on failure, and source is unchanged.
template <StorageScalar Scalar>
[[nodiscard]] core::Result<QualifiedProjection<Scalar>>
project_binary64(const PrecisionBuffer<Scalar>& source, const DenseMatrix<Scalar>& a,
                 const DenseMatrix<Scalar>& b, std::span<const double> scales, double tolerance,
                 core::NumericalGeneration target_generation, core::ScopeId consumer,
                 const ProjectionVerifier<Scalar>& verifier, core::ResourcePlan& plan) noexcept;

} // namespace uvlas::numerics
