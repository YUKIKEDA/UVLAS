#pragma once

#include "case_pg_reference.hpp"
#include "external_mesh_quantity.hpp"
#include "external_mesh_sidecar_binding.hpp"
#include "external_planar_pg.hpp"
#include "nastran_pcomp_physical.hpp"

namespace uvlas::io::detail {

struct ExternalSidecarPgLimits {
    ExternalSidecarPgLimits(const ExternalPgLimits& pg, GmshTopologyLimits topology_limits,
                            std::uint64_t thickness_vertices, std::uint64_t work_limit) noexcept
        : geometry(pg),
          topology(topology_limits),
          case_thickness_vertices(thickness_vertices),
          work(work_limit) {}

    ExternalPgLimits geometry;
    GmshTopologyLimits topology;
    std::uint64_t case_thickness_vertices;
    std::uint64_t work;
};

struct ExternalSidecarPgDiagnostic {
    std::optional<MeshSidecarDiagnostic> binding{};
    std::optional<ExternalPgDiagnostic> geometry{};
    std::optional<GmshTopologyDiagnostic> topology{};
    std::optional<ImportedQuantityMismatch> quantity{};
    std::optional<NastranPcompGeometryDiagnostic> laminate{};
};

struct QualifiedExternalSidecarPg {
    ExternalPlanarPg geometry;
    MeshShellId shell{0};
    MeshGeometryId pg{0};
    RootSnapshotLocation sidecar;
    SourceDigest sidecar_digest;
    bool thickness_compared{};
    std::uint64_t materials_compared{};
    std::uint64_t work{};
    // Original normal agrees with the PG normal in slot 0 and is reversed in slot 1.
    std::array<std::optional<NastranPcompGeometryMatch>, 2> laminate{};
};

struct ExternalPgConversionObservation {
    core::ObservationId id{0};
};

struct CertifiedExternalSidecarPgBounds {
    std::variant<GmshPgBounds, NastranPlanarPgBounds> geometry;
    MeshShellId shell{0};
    MeshGeometryId pg{0};
    RootSnapshotLocation sidecar;
    SourceDigest sidecar_digest;
    structure::ExactRational distance_budget_metres;
    structure::ExactRational normal_budget;
    core::PlannedBuffer<ExternalPgConversionObservation> observations;
};

struct QualifiedCaseSidecarPg {
    structure::PgReferenceGeometry reference;
    QualifiedExternalSidecarPg match;
};

// Derive the PG reference and units from the same Case and Shell as the retained match.
// The Nastran quantity owner must come from the same authenticated deck. Its material
// table is required for Nastran; Gmsh uses the Case material IDs already bound by the sidecar.
// This call certifies affine chart correspondence and supported explicit material and
// uniform-thickness matches, including ordered PCOMP plies. It retains no borrows.
[[nodiscard]] core::Result<QualifiedCaseSidecarPg> qualify_external_sidecar_planar_pg(
    const SiCase&, const ExternalMeshSidecarBinding&, MeshShellId, const NastranLaminateCards*,
    const CasePgReferenceLimits&, const ExternalSidecarPgLimits&, core::Workspace&,
    core::ResourcePlan&, const std::stop_token&, ExternalSidecarPgDiagnostic* = nullptr) noexcept;

struct CertifiedCaseSidecarPgBounds {
    structure::PgReferenceGeometry reference;
    CertifiedExternalSidecarPgBounds match;
};

// Verify the whole-face positional and oriented unit-normal bounds against
// the one declared Case or sidecar budget. Prepare the comparison reference from
// the selected shell of this SI Case. The returned observation identifiers still
// require independent physical evaluation before conversion acceptance.
[[nodiscard]] core::Result<CertifiedCaseSidecarPgBounds> certify_external_sidecar_pg_bounds(
    const SiCase&, const ExternalMeshSidecarBinding&, MeshShellId, const CasePgReferenceLimits&,
    const ExternalSidecarPgLimits&, std::uint64_t maximum_directed_work, core::Workspace&,
    core::ResourcePlan&, const std::stop_token&, ExternalSidecarPgDiagnostic* = nullptr) noexcept;

} // namespace uvlas::io::detail
