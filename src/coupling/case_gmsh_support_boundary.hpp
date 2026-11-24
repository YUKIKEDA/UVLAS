#pragma once

#include "case_external_pg.hpp"

namespace uvlas::coupling::detail {

struct CaseGmshSupportBoundaryLimits {
    structure::PgReferenceLimits arithmetic;
    std::uint64_t maximum_lines{};
};

enum class CaseGmshSupportBoundaryFailure : std::uint8_t {
    Binding,
    Support,
    Geometry,
    Incidence,
    DuplicateLine,
    OutsideBoundary,
    IncompleteBoundary,
};

struct CaseGmshSupportBoundaryDiagnostic {
    CaseGmshSupportBoundaryFailure reason = CaseGmshSupportBoundaryFailure::Binding;
    io::detail::MeshSupportId support{0};
    std::optional<io::detail::GmshKey> physical{};
    std::optional<core::Index> declaration{};
    std::optional<std::uint64_t> line{};
    std::optional<std::uint64_t> face{};
    std::optional<core::Index> boundary{};
};

struct CaseGmshSupportLine {
    std::uint64_t element{};
    std::array<std::uint64_t, 2> nodes{}; // Original Gmsh line order.
    std::uint64_t incident_face{};
};

struct CaseGmshSupportBoundary {
    io::detail::RootSnapshotLocation mesh;
    io::SourceDigest mesh_digest;
    io::detail::RootSnapshotLocation sidecar;
    io::SourceDigest sidecar_digest;
    io::detail::MeshSupportId support{0};
    io::detail::MeshShellId shell{0};
    io::detail::GmshKey physical{};
    core::Index declaration{};
    core::Index boundary{};
    structure::PgReferenceIdentity patch;
    core::PlannedBuffer<CaseGmshSupportLine> lines;
};

// Certify the original Gmsh line set against a named Case Shell outer boundary.
// A fit certificate is admissible here only if its relevant face has zero distance
// bound. The result attests to center-surface boundary provenance and coverage.
[[nodiscard]] core::Result<CaseGmshSupportBoundary> certify_case_gmsh_support_boundary(
    const io::SiCase&, const io::detail::ExternalMeshSidecarBinding&, const CaseGmshPgSet&,
    io::detail::MeshSupportId, io::detail::GmshKey, const CaseGmshSupportBoundaryLimits&,
    core::ResourcePlan&, std::optional<CaseGmshSupportBoundaryDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::coupling::detail
