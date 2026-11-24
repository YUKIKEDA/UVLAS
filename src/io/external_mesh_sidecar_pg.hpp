#pragma once

#include "external_mesh_quantity.hpp"
#include "external_mesh_sidecar_binding.hpp"
#include "external_planar_pg.hpp"

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

// The caller derives the PG reference and units from the specified Case and Shell.
// The Nastran quantity owner must come from the same authenticated deck. Its material
// table is required for Nastran; Gmsh uses the Case material IDs already bound by the sidecar.
// This call certifies affine chart correspondence and supported explicit material and
// uniform-thickness matches. It retains no borrows and leaves both original owners intact.
[[nodiscard]] core::Result<QualifiedExternalSidecarPg> qualify_external_sidecar_planar_pg(
    const SchemaCase&, const ExternalMeshSidecarBinding&, MeshShellId,
    const structure::PgReferenceGeometry&, const SiUnits&, const NastranLaminateCards*,
    const ExternalSidecarPgLimits&, core::ResourcePlan&, const std::stop_token&,
    ExternalSidecarPgDiagnostic* = nullptr) noexcept;

// Verify the whole-face positional and oriented unit-normal bounds against
// the one declared Case or sidecar budget. The caller prepares pg from the
// selected shell of this SI Case. The returned observation identifiers still
// require independent physical evaluation before conversion acceptance.
[[nodiscard]] core::Result<CertifiedExternalSidecarPgBounds> certify_external_sidecar_pg_bounds(
    const SiCase&, const ExternalMeshSidecarBinding&, MeshShellId,
    const structure::PgReferenceGeometry&, const ExternalSidecarPgLimits&,
    std::uint64_t maximum_directed_work, core::ResourcePlan&, const std::stop_token&,
    ExternalSidecarPgDiagnostic* = nullptr) noexcept;

} // namespace uvlas::io::detail
