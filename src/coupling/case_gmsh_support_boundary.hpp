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
    std::optional<std::uint8_t> world_component{}; // Zero-based Cartesian mismatch.
};

enum class CaseGmshSupportPositionEvidence : std::uint8_t {
    FaceEquality,
    EdgePolynomial,
};

struct CaseGmshSupportLine {
    std::uint64_t element{};
    std::array<std::uint64_t, 2> nodes{}; // Original Gmsh line order.
    std::uint64_t incident_face{};
    CaseGmshSupportPositionEvidence position = CaseGmshSupportPositionEvidence::FaceEquality;
};

struct CaseGmshEdgeEvidence {
    io::detail::MeshShellId shell{0};
    io::detail::GmshKey physical{};
    core::Index boundary{};
    structure::PgReferenceIdentity patch;
    core::PlannedBuffer<CaseGmshSupportLine> lines;
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
// A zero face distance bound suffices. Otherwise the original affine line and
// the restricted PG position must have identical exact polynomial coefficients.
// The result attests to center-surface boundary provenance and coverage.
[[nodiscard]] core::Result<CaseGmshSupportBoundary> certify_case_gmsh_support_boundary(
    const io::SiCase&, const io::detail::ExternalMeshSidecarBinding&, const CaseGmshPgSet&,
    io::detail::MeshSupportId, io::detail::GmshKey, const CaseGmshSupportBoundaryLimits&,
    core::ResourcePlan&, std::optional<CaseGmshSupportBoundaryDiagnostic>* = nullptr) noexcept;

struct CaseGmshJointBoundaryDiagnostic {
    io::detail::MeshJointId joint{0};
    std::optional<std::size_t> endpoint{}; // Case adhesive order, zero or one.
    CaseGmshSupportBoundaryDiagnostic edge;
};

struct CaseGmshJointBoundary {
    io::detail::RootSnapshotLocation mesh;
    io::SourceDigest mesh_digest;
    io::detail::RootSnapshotLocation sidecar;
    io::SourceDigest sidecar_digest;
    io::detail::MeshJointId joint{0};
    core::Index declaration{};
    std::array<CaseGmshEdgeEvidence, 2> endpoints;
};

/// Certify both original line groups against their respective Case adhesive regions.
[[nodiscard]] core::Result<CaseGmshJointBoundary> certify_case_gmsh_joint_boundary(
    const io::SiCase&, const io::detail::ExternalMeshSidecarBinding&, const CaseGmshPgSet&,
    io::detail::MeshJointId, const std::array<io::detail::GmshKey, 2>&,
    const CaseGmshSupportBoundaryLimits&, core::ResourcePlan&,
    std::optional<CaseGmshJointBoundaryDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::coupling::detail
