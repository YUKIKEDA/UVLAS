#pragma once

#include "case_material_original_layout.hpp"
#include "case_pg_support_motion.hpp"

namespace uvlas::coupling::detail {

/// Physical panel IDs issued by the discrete geometry assembler, in Case
/// chord-major order. Every Case surface contributes one row to the table.
struct CaseAeroSurfacePanels {
    std::uint64_t surface{};
    std::span<const aero::detail::PhysicalFaceId> physical_panels;
};

struct CaseMaterialSurfaceModelLimits {
    CaseMaterialChartLimits material;
    structure::PgReferenceLimits alias;
    std::uint64_t maximum_surfaces{}, maximum_bound_panels{}, maximum_alias_checks{};
    std::uint64_t maximum_assembly_steps{};
};

struct CaseMaterialSurfaceCorner {
    CaseMaterialChartIdentity surface;
    CaseMaterialOriginalCorner original;
};

/// Owns all structural surfaces' original load routes for one checked contact.
/// Registered canonical vertices are qualified across each incident material site
/// before the model is published.
class CaseMaterialSurfaceModel {
public:
    CaseMaterialSurfaceModel(CaseMaterialSurfaceModel&&) noexcept;
    CaseMaterialSurfaceModel& operator=(CaseMaterialSurfaceModel&&) noexcept;
    CaseMaterialSurfaceModel(const CaseMaterialSurfaceModel&) = delete;
    CaseMaterialSurfaceModel& operator=(const CaseMaterialSurfaceModel&) = delete;
    ~CaseMaterialSurfaceModel();

    [[nodiscard]] static core::Result<CaseMaterialSurfaceModel>
    prepare(const io::SiCase&, std::span<const CaseAeroSurfacePanels>,
            const aero::detail::GlobalContactState&, const CaseStructuralSupport&,
            std::span<const structure::PgConstraintPatchBinding>,
            const CaseMaterialSurfaceModelLimits&, core::Workspace&, core::ResourcePlan&) noexcept;

    [[nodiscard]] std::size_t surface_count() const noexcept;

    [[nodiscard]] core::Result<CaseMaterialSurfaceCorner>
    corner(aero::detail::FacetId, std::uint8_t local_vertex) const noexcept;

    [[nodiscard]] core::Result<CaseMaterialOriginalEnclosedForce>
    pullback(const DiscreteCirculationLoads&, std::size_t, const structure::PgPreparedPatch&,
             const structure::PgPreparedReferenceProjection&, const structure::PgReferenceLimits&,
             const structure::PgVolumePositionLimits&, core::ResourcePlan&) const noexcept;

private:
    struct Surface;
    struct Route;
    CaseMaterialSurfaceModel(core::PlannedBuffer<std::optional<Surface>>,
                             core::PlannedBuffer<Route>) noexcept;
    [[nodiscard]] const Route* find(aero::detail::FacetId, std::uint8_t) const noexcept;
    core::PlannedBuffer<std::optional<Surface>> surfaces_;
    core::PlannedBuffer<Route> routes_;
};

} // namespace uvlas::coupling::detail
