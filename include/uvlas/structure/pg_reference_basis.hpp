#pragma once

#include <uvlas/structure/pg_rank_preparation.hpp>
#include <uvlas/structure/pg_reference_geometry.hpp>

namespace uvlas::structure {

/// An owned reference point and one selected, unconditioned scalar basis column.
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
    /// Borrow source maps and selected order for subsequent basis transforms.
    [[nodiscard]] core::Result<const PgRankPreparation*> rank() const& noexcept;
    [[nodiscard]] core::Result<const PgRankPreparation*> rank() const&& = delete;
    [[nodiscard]] core::Result<PgReferenceSource> source() const& noexcept;
    [[nodiscard]] core::Result<PgReferenceSource> source() const&& = delete;

private:
    PgReferenceBasis(PgReferenceGeometry&&, PgRankPreparation&&) noexcept;

    PgReferenceGeometry geometry_;
    PgRankPreparation rank_;
};
} // namespace uvlas::structure
