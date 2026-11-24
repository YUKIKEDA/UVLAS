#pragma once

#include "reference_resolution.hpp"

namespace uvlas::io::detail {

struct MeshSidecarBindingLimits {
    MeshSidecarBindingLimits(MeshSidecarLimits mapping, GmshLimits mesh, NastranBulkLimits bulk,
                             NastranCoordinateLimits systems, NastranShellLimits elements,
                             NastranSupportLimits supports, std::uint64_t face_limit,
                             std::uint64_t work_limit) noexcept
        : sidecar(mapping),
          gmsh(mesh),
          nastran(bulk),
          coordinates(systems),
          shells(elements),
          support(supports),
          faces(face_limit),
          work(work_limit) {}

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
class ExternalMeshSidecarBinding {
public:
    using Source = std::variant<GmshMesh, NastranShellCards>;

    ExternalMeshSidecarBinding(ExternalMeshSidecarBinding&&) noexcept = default;
    ExternalMeshSidecarBinding& operator=(ExternalMeshSidecarBinding&&) = delete;

    [[nodiscard]] const ParsedMeshSidecar& sidecar() const& noexcept {
        return sidecar_;
    }
    [[nodiscard]] const ParsedMeshSidecar& sidecar() const&& = delete;
    [[nodiscard]] const Source& source() const& noexcept {
        return source_;
    }
    [[nodiscard]] const Source& source() const&& = delete;
    [[nodiscard]] std::span<const BoundExternalFace> faces() const& noexcept {
        return faces_.values();
    }
    [[nodiscard]] std::span<const BoundExternalFace> faces() const&& = delete;
    [[nodiscard]] std::uint64_t work() const noexcept {
        return work_;
    }

private:
    friend struct ExternalMeshSidecarBindingFactory;

    ExternalMeshSidecarBinding(ParsedMeshSidecar&& sidecar, Source&& source,
                               core::PlannedBuffer<BoundExternalFace>&& faces,
                               std::uint64_t work) noexcept
        : sidecar_(std::move(sidecar)),
          source_(std::move(source)),
          faces_(std::move(faces)),
          work_(work) {}

    ParsedMeshSidecar sidecar_;
    Source source_;
    core::PlannedBuffer<BoundExternalFace> faces_;
    std::uint64_t work_;
};

// Resolve identity and tag correspondence. Chart coverage and qualification of physical
// quantities for the same material region use this owned input in the preparation phase.
[[nodiscard]] core::Result<ExternalMeshSidecarBinding>
bind_external_mesh_sidecar(const SchemaCase&, const ReferencedCase&, MeshGeometryId,
                           const MeshSidecarBindingLimits&, core::ResourcePlan&,
                           const std::stop_token&,
                           std::optional<MeshSidecarDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::io::detail
