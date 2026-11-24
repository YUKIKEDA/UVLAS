#pragma once

#include <uvlas/structure/pg_rank_preparation.hpp>
#include <uvlas/structure/pg_reference_geometry.hpp>

namespace uvlas::structure {

namespace detail {
struct BasisReferencePoint;
struct BasisReferenceSample;
} // namespace detail

/// An owned reference point and one scalar basis column.
/// All scalar components use the radicand in the same reference point.
struct PgReferenceBasisValue {
    PgReferenceValue reference;
    PgRankBasisColumn column;
    PgReferenceScalar value;                       ///< Dimensionless N.
    PgReferenceVector chart_derivatives;           ///< Partial derivatives in s,t,zeta order.
    PgReferenceVector physical_gradient_per_metre; ///< Global Cartesian order.
};

/// Couples a selected local basis to one immutable, physically qualified reference.
class PgReferenceBasis {
public:
    PgReferenceBasis(PgReferenceBasis&&) noexcept = default;
    PgReferenceBasis& operator=(PgReferenceBasis&&) noexcept = default;
    PgReferenceBasis(const PgReferenceBasis&) = delete;
    PgReferenceBasis& operator=(const PgReferenceBasis&) = delete;

    /// Revalidate the local rank against the actual reference input owners.
    /// Consume both candidates only after successful revalidation.
    [[nodiscard]] static core::Result<PgReferenceBasis>
    prepare(PgReferenceGeometry&&, PgRankPreparation&&, const PgRankPreparationLimits&,
            core::Workspace&, core::ResourcePlan&) noexcept;

    /// Return exact values for a selected column on the requested partition side.
    [[nodiscard]] core::Result<PgReferenceBasisValue> evaluate(const PgReferencePoint&,
                                                               PgRankBasisColumn,
                                                               const PgReferenceLimits&,
                                                               core::ResourcePlan&) const noexcept;

    [[nodiscard]] core::Result<std::uint64_t> size() const noexcept;
    [[nodiscard]] core::Result<PgReferenceIdentity> reference_identity() const noexcept;
    [[nodiscard]] core::Result<PgReferenceFrame> frame() const& noexcept;
    [[nodiscard]] core::Result<PgReferenceFrame> frame() const&& = delete;
    /// Borrow source maps and selected order for subsequent basis transforms.
    [[nodiscard]] core::Result<const PgRankPreparation*> rank() const& noexcept;
    [[nodiscard]] core::Result<const PgRankPreparation*> rank() const&& = delete;
    [[nodiscard]] core::Result<PgReferenceSource> source() const& noexcept;
    [[nodiscard]] core::Result<PgReferenceSource> source() const&& = delete;

private:
    friend class PgConditioningCertificate;
    friend class PgFlatGram;
    friend class PgDiscreteGram;
    [[nodiscard]] const ExactRational& length_metres() const& noexcept;
    [[nodiscard]] detail::BasisReferencePoint locate(detail::ReferenceArithmetic&,
                                                     const PgReferencePoint&) const;
    [[nodiscard]] std::array<std::array<ExactRational, 3>, 3>
    numerator_jet(detail::ReferenceArithmetic&, PgRankBasisColumn,
                  const detail::BasisReferencePoint&) const;
    [[nodiscard]] detail::BasisReferenceSample
    sample(detail::ReferenceArithmetic&, const detail::BasisReferencePoint&,
           const std::array<std::array<ExactRational, 3>, 3>&) const;
    [[nodiscard]] PgReferenceBasisValue
    assemble(detail::ReferenceArithmetic&, detail::BasisReferencePoint&&, PgRankBasisColumn,
             const std::array<std::array<ExactRational, 3>, 3>&) const;

    PgReferenceBasis(PgReferenceGeometry&&, PgRankPreparation&&) noexcept;

    PgReferenceGeometry geometry_;
    PgRankPreparation rank_;
};
} // namespace uvlas::structure
