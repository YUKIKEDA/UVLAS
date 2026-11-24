#pragma once

#include <uvlas/structure/pg_reference_basis.hpp>
#include <variant>

namespace uvlas::structure {

/// Bounds for one independently enclosed entry of the dimensionless, unit-density
/// discrete geometric Gram matrix.
struct PgGramEntry {
    ExactRational lower, upper;
};

/// Symmetric upper triangle in column order (0,0), (0,1), (1,1), ... .
/// Each endpoint must be exact on the declared lattice 2^(-fractional_bits).
struct PgDiscreteGramEnclosure {
    std::span<const PgGramEntry> upper_triangle;
    std::uint32_t fractional_bits = 0;
};

/// An exact upper-triangular candidate T in the same packed column order.
struct PgConditioningCandidate {
    std::span<const ExactRational> upper_triangle;
    std::uint32_t fractional_bits = 0;
};

/// Explicit dyadic grids for the approximate factor and the certified result.
struct PgConditioningGeneration {
    std::uint32_t factor_fractional_bits = 0;
    std::uint32_t candidate_fractional_bits = 0;
};

struct PgConditioningLimits {
    std::uint64_t maximum_payload_bytes;
    std::uint64_t maximum_steps;
    std::stop_token cancellation;
};

struct PgConditioningInfo {
    std::uint64_t rank;
    std::uint64_t protected_rank;
    std::uint64_t payload_bytes;
    std::uint64_t steps;
};

/// One conditioned column at a point. Components use the corresponding point's g.
struct PgConditionedReferenceColumn {
    PgReferenceScalar value;                       ///< Dimensionless N.
    PgReferenceVector physical_gradient_per_metre; ///< Global Cartesian G.
};

/// Owned point order and all conditioned columns in point-major order.
/// Column j at point p is columns[p * column_count + j], zero-based.
struct PgConditionedReferenceBatch {
    std::uint64_t column_count;
    core::PlannedBuffer<PgReferenceValue> points;
    core::PlannedBuffer<PgConditionedReferenceColumn> columns;
    PgReferenceCallInfo call_info;
};

namespace detail {
struct PgConditioningState;
}

/// Certifies the complete discrete congruence against the caller's independently
/// qualified Gram enclosure. The caller supplies the enclosure for this actual basis.
/// The factory owns the basis and exact transform after successful certification.
class PgConditioningCertificate {
public:
    PgConditioningCertificate(PgConditioningCertificate&&) noexcept;
    PgConditioningCertificate& operator=(PgConditioningCertificate&&) noexcept;
    ~PgConditioningCertificate();

    [[nodiscard]] static core::Result<PgConditioningCertificate>
    certify(PgReferenceBasis&&, PgDiscreteGramEnclosure, PgConditioningCandidate,
            const PgConditioningLimits&, core::Workspace&, core::ResourcePlan&) noexcept;

    /// Generate a midpoint factor and certify the actual candidate against the
    /// original enclosure. The caller binds the enclosure to this actual basis.
    [[nodiscard]] static core::Result<PgConditioningCertificate>
    generate_and_certify(PgReferenceBasis&&, PgDiscreteGramEnclosure, PgConditioningGeneration,
                         const PgConditioningLimits&, core::Workspace&,
                         core::ResourcePlan&) noexcept;

    [[nodiscard]] core::Result<PgConditioningInfo> info() const noexcept;
    /// Exact infinity-norm upper bound.
    [[nodiscard]] core::Result<ExactRational> congruence_bound(const ExactLimits&,
                                                               core::ResourcePlan&) const noexcept;

    /// The same admitted physical reference and its unconditioned rank.
    [[nodiscard]] core::Result<const PgReferenceBasis*> reference_basis() const& noexcept;
    [[nodiscard]] core::Result<const PgReferenceBasis*> reference_basis() const&& = delete;

    /// The exact T entry, with a zero below its diagonal.
    [[nodiscard]] core::Result<ExactRational>
    transform_coefficient(PgRankBasisColumn row, PgRankBasisColumn column, const ExactLimits&,
                          core::ResourcePlan&) const noexcept;

    /// Evaluate one conditioned column of B T on the requested reference side.
    /// The reference, T, and the selected order all belong to this certificate.
    [[nodiscard]] core::Result<PgReferenceBasisValue> evaluate(const PgReferencePoint&,
                                                               PgRankBasisColumn,
                                                               const PgReferenceLimits&,
                                                               core::ResourcePlan&) const noexcept;

    /// Evaluate every conditioned column for every point in input order.
    /// Points retain the requested layer/face side, and the whole batch shares one budget.
    [[nodiscard]] core::Result<PgConditionedReferenceBatch>
    evaluate_reference(std::span<const PgReferencePoint>, const PgReferenceLimits&,
                       core::ResourcePlan&) const noexcept;

    /// Exact coefficient of B T, indexed by the original canonical coefficient row.
    [[nodiscard]] core::Result<ExactRational> coefficient(PgRankBasisColumn, core::Index row,
                                                          const PgReferenceLimits&,
                                                          core::ResourcePlan&) const noexcept;

    /// Exact coordinates in B T of a mandatory field, including all trailing zeros.
    [[nodiscard]] core::Result<core::PlannedBuffer<ExactRational>>
    required_coordinates(PgRankRequiredColumn, const PgReferenceLimits&, core::Workspace&,
                         core::ResourcePlan&) const noexcept;

private:
    [[nodiscard]] static core::Result<PgConditioningCertificate>
    certify_impl(PgReferenceBasis&&, PgDiscreteGramEnclosure,
                 std::variant<PgConditioningCandidate, PgConditioningGeneration>,
                 const PgConditioningLimits&, core::Workspace&, core::ResourcePlan&) noexcept;

    PgConditioningCertificate(PgReferenceBasis&&,
                              core::PlannedBuffer<detail::PgConditioningState>) noexcept;
    PgReferenceBasis basis_;
    core::PlannedBuffer<detail::PgConditioningState> state_;
};

} // namespace uvlas::structure
