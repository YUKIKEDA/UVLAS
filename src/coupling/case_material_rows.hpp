#pragma once

#include "case_material_reference.hpp"
#include "discrete_nonpenetration.hpp"

#include <optional>
#include <uvlas/structure/pg_material_tangents.hpp>

namespace uvlas::coupling::detail {

struct CaseMaterialCollocationSite {
    core::Index chord_cell{}, span_cell{};
};
struct CaseMaterialVertexSite {
    core::Index chord_column{}, span_knot{}, original_span_cell{};
};
using CaseMaterialSite = std::variant<CaseMaterialCollocationSite, CaseMaterialVertexSite>;

struct CaseMaterialPointKinematics {
    structure::PgVolumePositionResult position, velocity;
    std::optional<structure::PgVolumePositionResult> tangents;
};

struct CaseMaterialOriginalForce {
    CaseMaterialChartIdentity surface;
    CaseMaterialSite original_site;
    structure::PgPreparedPatchIdentity patch;
    structure::PgVolumePositionResult force;
};

struct CaseMaterialOriginalEnclosedForce {
    CaseMaterialChartIdentity surface;
    CaseMaterialVertexSite original_site;
    structure::PgPreparedPatchIdentity patch;
    core::PlannedBuffer<std::optional<DiscreteLoadVector>> force;
    std::uint64_t used_steps{};
};

struct CaseMaterialActionLimits {
    std::uint64_t maximum_steps{};
    double position_error{}, velocity_error{}, tangent_error{};
    std::stop_token cancellation;
};

/// Conditioned rows of one actual material point. A vertex records the original
/// incident span cell, even when the flow position uses a canonical vertex.
class CaseMaterialRows {
public:
    CaseMaterialRows(CaseMaterialRows&&) noexcept = default;
    CaseMaterialRows& operator=(CaseMaterialRows&&) noexcept = default;
    CaseMaterialRows(const CaseMaterialRows&) = delete;
    CaseMaterialRows& operator=(const CaseMaterialRows&) = delete;

    [[nodiscard]] static core::Result<CaseMaterialRows>
    prepare(const CaseMaterialChart&, const CaseReferenceSurface&,
            const structure::PgPreparedPatch&, CaseMaterialSite,
            const structure::PgPreparedReferenceProjection&, double maximum_tangent_row_error,
            const structure::PgReferenceLimits&, core::ResourcePlan&) noexcept;

    [[nodiscard]] core::Result<CaseMaterialPointKinematics>
    evaluate(structure::PgPreparedPatchIdentity, structure::PgCartesianWords q,
             structure::PgCartesianWords v, const CaseMaterialActionLimits&, core::Workspace&,
             core::ResourcePlan&) const noexcept;

    /// Send this contribution to its own conditioned structural row. The returned
    /// force remains in the original patch's layout, prior to global constraints.
    [[nodiscard]] core::Result<CaseMaterialOriginalForce>
    pullback(structure::PgPreparedPatchIdentity, structure::PgCartesianWords,
             const structure::PgVolumePositionLimits&, core::Workspace&,
             core::ResourcePlan&) const noexcept;

    /// Carry the force enclosure from the original aerodynamic contribution through
    /// the selected vertex row. The caller supplies the resolved original row.
    [[nodiscard]] core::Result<CaseMaterialOriginalEnclosedForce>
    pullback_enclosed(structure::PgPreparedPatchIdentity, const DiscreteLoadVector&,
                      const structure::PgVolumePositionLimits&, core::ResourcePlan&) const noexcept;

private:
    CaseMaterialRows(CaseMaterialChartIdentity, CaseMaterialSite,
                     structure::PgProjectedReferenceBatch&&,
                     std::optional<structure::PgMaterialTangentReference>&&) noexcept;
    CaseMaterialChartIdentity surface_;
    CaseMaterialSite site_;
    structure::PgProjectedReferenceBatch position_;
    std::optional<structure::PgMaterialTangentReference> tangents_;
};

} // namespace uvlas::coupling::detail
