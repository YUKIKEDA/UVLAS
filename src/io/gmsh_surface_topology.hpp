#pragma once

#include "gmsh_mesh.hpp"

namespace uvlas::io::detail {

struct GmshTopologySource {
    RootSnapshotLocation location;
    SourceDigest digest;
};

struct GmshTopologyFace {
    GmshElement element;
    std::size_t first_edge = 0;
    std::size_t component = 0;
};

struct GmshTopologyEdge {
    std::uint64_t from = 0, to = 0;
    std::size_t face = 0, next = 0;
    std::optional<std::size_t> opposite;
};

struct GmshTopologyComponent {
    std::size_t first_face = 0; // First member ordinal; other members need not be contiguous.
    std::size_t face_count = 0;
    std::size_t boundary_count = 0;
};

enum class GmshTopologyFailure : std::uint8_t {
    Input,
    UndefinedFace,
    NonSurfaceElement,
    RepeatedSelection,
    RepeatedVertex,
    DuplicateFace,
    NonManifoldEdge,
    InconsistentOrientation,
    DisconnectedVertexFan,
    Limit
};

struct GmshTopologyDiagnostic {
    GmshTopologyFailure reason = GmshTopologyFailure::Input;
    std::uint64_t face = 0, other_face = 0;
    std::uint64_t vertex = 0, other_vertex = 0;
};

struct GmshTopologyLimits {
    std::uint64_t faces;
    std::uint64_t directed_edges;
    std::uint64_t boundary_edges;
    std::uint64_t components;
    std::uint64_t work;
};

struct GmshTopologyView {
    GmshTopologySource source;
    // Faces retain explicit selection order and their original element rings.
    std::span<const GmshTopologyFace> faces;
    std::span<const GmshTopologyEdge> edges;
    std::span<const std::size_t> boundary;             // Edge ordinals in face and ring order.
    std::span<const GmshTopologyComponent> components; // Ordered by first selected face.
    std::uint64_t work;
};

// Owns the integer incidence of an explicit face selection from one decoded mesh.
class GmshSurfaceTopology {
public:
    GmshSurfaceTopology(GmshSurfaceTopology&&) noexcept;
    GmshSurfaceTopology& operator=(GmshSurfaceTopology&&) = delete;

    [[nodiscard]] static core::Result<GmshSurfaceTopology>
    prepare(const GmshMesh&, std::span<const std::uint64_t> selected_tags, GmshTopologyLimits,
            core::ResourcePlan&, const std::stop_token&,
            GmshTopologyDiagnostic* = nullptr) noexcept;

    [[nodiscard]] std::optional<GmshTopologyView> view() const& noexcept;
    [[nodiscard]] std::optional<GmshTopologyView> view() const&& = delete;

private:
    GmshSurfaceTopology(GmshTopologySource, core::PlannedBuffer<GmshTopologyFace>&&,
                        core::PlannedBuffer<GmshTopologyEdge>&&, core::PlannedBuffer<std::size_t>&&,
                        core::PlannedBuffer<GmshTopologyComponent>&&, std::uint64_t) noexcept;

    std::optional<GmshTopologySource> source_;
    core::PlannedBuffer<GmshTopologyFace> faces_;
    core::PlannedBuffer<GmshTopologyEdge> edges_;
    core::PlannedBuffer<std::size_t> boundary_;
    core::PlannedBuffer<GmshTopologyComponent> components_;
    std::uint64_t work_;
};

} // namespace uvlas::io::detail
