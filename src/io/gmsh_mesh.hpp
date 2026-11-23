#pragma once

#include "reference_root_snapshot.hpp"

#include <compare>
#include <string_view>

namespace uvlas::io::detail {

struct GmshLimits {
    std::uint64_t file_bytes;
    std::uint64_t physical_names;
    std::uint64_t entities;
    std::uint64_t entity_links;
    std::uint64_t nodes;
    std::uint64_t elements;
    std::uint64_t work;
};

struct GmshKey {
    std::uint8_t dimension = 0;
    std::uint64_t tag = 0;
    auto operator<=>(const GmshKey&) const = default;
};

// Offsets into the raw bytes owned by the result, not rounded coordinates.
struct GmshLexeme {
    std::size_t offset = 0;
    std::size_t length = 0;
};

struct GmshNames {
    GmshKey group;
    GmshLexeme name;
};

struct GmshLink {
    GmshKey key;
    bool reversed = false;
};

struct GmshEntity {
    GmshKey key;
    std::array<GmshLexeme, 6> bounds;
    std::size_t first_physical = 0, physical_count = 0;
    std::size_t first_boundary = 0, boundary_count = 0;
};

struct GmshNode {
    std::uint64_t tag = 0;
    GmshKey entity;
    std::array<GmshLexeme, 3> coordinates;
    std::array<GmshLexeme, 3> parameters;
    bool parametric = false;
};

enum class GmshElementKind : std::uint8_t { Point = 15, Line = 1, Triangle = 2, Quadrangle = 3 };

struct GmshElement {
    std::uint64_t tag = 0;
    GmshKey entity;
    GmshElementKind kind = GmshElementKind::Point;
    std::array<std::uint64_t, 4> nodes{};

    [[nodiscard]] std::span<const std::uint64_t> connectivity() const& noexcept;
    [[nodiscard]] std::span<const std::uint64_t> connectivity() const&& = delete;
};

struct GmshMeshView {
    RootSnapshotLocation location;
    SourceDigest digest;
    std::string_view bytes;
    std::span<const GmshNames> physical_names;
    std::span<const GmshEntity> entities;
    std::span<const GmshLink> entity_links;
    std::span<const GmshNode> nodes;
    std::span<const GmshElement> elements;
    bool entities_declared;
    std::uint64_t work;

    [[nodiscard]] std::string_view text(GmshLexeme value) const noexcept;
};

// Owns a decoded MSH 4.1 ASCII payload. Missing optional metadata stays missing.
// This layer resolves tags and does not decide material, support, or reference geometry.
class GmshMesh {
public:
    GmshMesh(GmshMesh&&) noexcept = default;
    GmshMesh& operator=(GmshMesh&&) noexcept = default;

    [[nodiscard]] static core::Result<GmshMesh>
    prepare(const AcquiredRootSnapshot&, GmshLimits, core::ResourcePlan&, const std::stop_token&,
            std::size_t* failure_offset = nullptr) noexcept;

    [[nodiscard]] std::optional<GmshMeshView> view() const& noexcept;
    [[nodiscard]] std::optional<GmshMeshView> view() const&& = delete;

private:
    GmshMesh() = default;

    std::optional<RootSnapshotLocation> location_;
    SourceDigest digest_{};
    core::PlannedBuffer<char> bytes_;
    core::PlannedBuffer<GmshNames> names_;
    core::PlannedBuffer<GmshEntity> entities_;
    core::PlannedBuffer<GmshLink> links_;
    core::PlannedBuffer<GmshNode> nodes_;
    core::PlannedBuffer<GmshElement> elements_;
    bool entities_declared_ = false;
    std::uint64_t work_ = 0;
};

} // namespace uvlas::io::detail
