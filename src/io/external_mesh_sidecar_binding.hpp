#pragma once

#include "reference_resolution.hpp"

namespace uvlas::io::detail {

struct MeshSidecarBindingLimits {
    MeshSidecarLimits sidecar;
    GmshLimits gmsh;
    NastranBulkLimits nastran;
    NastranCoordinateLimits coordinates;
    NastranShellLimits shells;
    NastranSupportLimits support;
    std::uint64_t faces;
    std::uint64_t work;
};

struct BoundExternalFace {
    MeshShellId shell{0};
    MeshGeometryId pg{0};
    std::uint64_t element{};
    std::array<std::uint64_t, 4> nodes{};
    std::uint8_t count{};
    core::Index sidecar_node{};
    std::size_t source_ordinal{}; // MSH element index or original BDF card index.
    std::size_t sidecar_row{};    // Row of shells in the typed sidecar.
    std::size_t sidecar_face{};   // Index within that row for Gmsh, zero for Nastran.
};

// Output faces are in sidecar row order, then element-tag order within a Gmsh row.
// Chart and original field data remain owned and carry no references to the input graph.
struct ExternalMeshSidecarBinding {
    ParsedMeshSidecar sidecar;
    std::variant<GmshMesh, NastranShellCards> source;
    core::PlannedBuffer<BoundExternalFace> faces;
    std::uint64_t work;
};

// Resolve identity and tag correspondence. Chart coverage and qualification of physical
// quantities for the same material region use this owned input in the preparation phase.
[[nodiscard]] core::Result<ExternalMeshSidecarBinding>
bind_external_mesh_sidecar(const SchemaCase&, const ReferencedCase&, MeshGeometryId,
                           const MeshSidecarBindingLimits&, core::ResourcePlan&,
                           const std::stop_token&,
                           std::optional<MeshSidecarDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::io::detail
