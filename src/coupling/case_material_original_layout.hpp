#pragma once

#include "case_material_rows.hpp"

namespace uvlas::coupling::detail {

class CaseStructuralSupport;
struct CaseMaterialOriginalAliasPoint;

struct CaseMaterialOriginalCorner {
    aero::detail::PhysicalFaceId face{0};
    aero::detail::FacetId facet{0};
    std::uint8_t local_vertex{};
    aero::detail::PointId canonical{0};
    CaseMaterialVertexSite site;
};

/// Maps the original checked facets of one ordered physical panel grid to Case
/// roles. Physical panel IDs come from the assembled topology, in chord-major
/// order. Child facet IDs, local corners and layer sides are derived internally.
class CaseMaterialOriginalLayout {
public:
    CaseMaterialOriginalLayout(CaseMaterialOriginalLayout&&) noexcept = default;
    CaseMaterialOriginalLayout& operator=(CaseMaterialOriginalLayout&&) noexcept = default;
    CaseMaterialOriginalLayout(const CaseMaterialOriginalLayout&) = delete;
    CaseMaterialOriginalLayout& operator=(const CaseMaterialOriginalLayout&) = delete;

    [[nodiscard]] static core::Result<CaseMaterialOriginalLayout>
    prepare(const CaseMaterialChart&, const CaseReferenceSurface&,
            const aero::detail::GlobalContactState&, std::span<const aero::detail::PhysicalFaceId>,
            const structure::PgReferenceLimits&, core::ResourcePlan&) noexcept;

    [[nodiscard]] core::Result<CaseMaterialOriginalCorner>
    corner(aero::detail::FacetId, std::uint8_t local_vertex) const noexcept;

    /// Route an actual original aerodynamic contribution to the selected PG row.
    [[nodiscard]] core::Result<CaseMaterialOriginalEnclosedForce>
    pullback(const DiscreteCirculationLoads&, std::size_t, const CaseMaterialChart&,
             const CaseReferenceSurface&, const structure::PgPreparedPatch&,
             const structure::PgPreparedReferenceProjection&, const structure::PgReferenceLimits&,
             const structure::PgVolumePositionLimits&, core::ResourcePlan&) const noexcept;

private:
    friend core::Result<void> qualify_case_material_original_alias(
        const CaseMaterialOriginalAliasPoint&, const CaseMaterialOriginalAliasPoint&,
        const CaseStructuralSupport&, const structure::PgReferenceLimits&, core::Workspace&,
        core::ResourcePlan&) noexcept;
    CaseMaterialOriginalLayout(CaseMaterialChartIdentity, aero::detail::GeometryContext,
                               std::array<std::byte, 32>,
                               core::PlannedBuffer<CaseMaterialOriginalCorner>) noexcept;
    CaseMaterialChartIdentity surface_;
    aero::detail::GeometryContext context_;
    std::array<std::byte, 32> geometry_;
    core::PlannedBuffer<CaseMaterialOriginalCorner> corners_;
};

/// Two original vertices selected from the same checked canonical point.
struct CaseMaterialOriginalAliasPoint {
    const CaseMaterialOriginalLayout& layout;
    aero::detail::FacetId facet;
    std::uint8_t local_vertex;
    const CaseMaterialChart& chart;
    const CaseReferenceSurface& thin;
    const structure::PgPreparedPatch& patch;
};

/// Verify an existing topological alias against one actual Case support owner.
[[nodiscard]] core::Result<void> qualify_case_material_original_alias(
    const CaseMaterialOriginalAliasPoint&, const CaseMaterialOriginalAliasPoint&,
    const CaseStructuralSupport&, const structure::PgReferenceLimits&, core::Workspace&,
    core::ResourcePlan&) noexcept;

} // namespace uvlas::coupling::detail
