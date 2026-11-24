#pragma once

#include "gmsh_surface_topology.hpp"
#include "nastran_shell_cards.hpp"

#include <uvlas/io/si_case.hpp>
#include <uvlas/structure/pg_reference_geometry.hpp>

namespace uvlas::io::detail {

enum class ExternalPgWinding : std::uint8_t { Forward, Reversed };

struct ExternalPgCorner {
    std::uint64_t node;
    const structure::PgPartitionPoint& point;
};

struct ExternalPgFace {
    ExternalPgFace(std::uint64_t element_id, ExternalPgWinding declared_winding,
                   std::span<const ExternalPgCorner> declared_corners) noexcept
        : element(element_id), winding(declared_winding), corners(declared_corners) {}

    std::uint64_t element;
    ExternalPgWinding winding;
    std::span<const ExternalPgCorner> corners; // Original element ring, not sorted by tag.
};

struct NastranGridPlacement {
    NastranCoordinateId position_system{0};
    std::optional<std::uint64_t> cp_byte;
    std::optional<std::array<structure::ExactRational, 3>> unshifted_metres;
};

struct NastranReferenceOffset {
    structure::ExactRational signed_metres;
    std::uint64_t byte{};
};

struct ExternalPgOwnedCorner {
    std::uint64_t node{};
    structure::PgPartitionPoint chart;
    std::array<std::optional<std::uint64_t>, 3> source_bytes;
    std::array<structure::ExactRational, 3> position_metres;
    std::optional<NastranGridPlacement> nastran;
};

struct ExternalPgOwnedFace {
    std::uint64_t element{};
    ExternalPgWinding winding = ExternalPgWinding::Forward;
    std::uint8_t count{};
    std::array<ExternalPgOwnedCorner, 4> corners;
    std::optional<NastranReferenceOffset> nastran_offset;

    [[nodiscard]] std::span<const ExternalPgOwnedCorner> ring() const& noexcept {
        return std::span{corners}.first(count);
    }
    [[nodiscard]] std::span<const ExternalPgOwnedCorner> ring() const&& = delete;
};

enum class ExternalPgFailure : std::uint8_t {
    Source,
    CurvedPatch,
    Selection,
    PositionSystem,
    CoordinateFrame,
    Offset,
    SourceOrder,
    Position,
    Winding,
    Outside,
    Overlap,
    Coverage,
    Planarity,
    SourceWinding,
    EdgeIncidence,
};

struct ExternalPgDiagnostic {
    ExternalPgFailure reason = ExternalPgFailure::Source;
    std::optional<std::uint64_t> element, other, node, byte;
    std::optional<NastranCoordinateId> coordinate_system;
};

struct ExternalPgLimits {
    structure::PgReferenceLimits arithmetic;
    std::uint64_t maximum_faces;
};

struct ExternalPlanarPgView {
    RootSnapshotLocation source;
    SourceDigest digest;
    SiLengthUnit units;
    structure::PgReferenceIdentity patch;
    std::span<const ExternalPgOwnedFace> faces;
    structure::PgReferenceCallInfo work;
};

// Certifies an external face selection against one complete, affine PG mid-surface.
// The calling importer binds the already prepared PG owner to its Case member.
// Result owns source-ordered SI coordinates and exact material chart coordinates.
class ExternalPlanarPg {
public:
    ExternalPlanarPg(ExternalPlanarPg&&) noexcept = default;
    ExternalPlanarPg& operator=(ExternalPlanarPg&&) = delete;

    [[nodiscard]] static core::Result<ExternalPlanarPg>
    prepare(const structure::PgReferenceGeometry&, const GmshMesh&, const GmshSurfaceTopology&,
            SiLengthUnit, std::span<const ExternalPgFace>, const ExternalPgLimits&,
            core::ResourcePlan&, std::optional<ExternalPgDiagnostic>* = nullptr) noexcept;

    [[nodiscard]] static core::Result<ExternalPlanarPg>
    prepare(const structure::PgReferenceGeometry&, const NastranShellCards&, SiLengthUnit,
            std::span<const ExternalPgFace>, const ExternalPgLimits&, core::ResourcePlan&,
            std::optional<ExternalPgDiagnostic>* = nullptr) noexcept;

    [[nodiscard]] std::optional<ExternalPlanarPgView> view() const& noexcept;
    [[nodiscard]] ExternalPlanarPgView view() const&& = delete;

private:
    ExternalPlanarPg(RootSnapshotLocation, SourceDigest, SiLengthUnit,
                     structure::PgReferenceIdentity, core::PlannedBuffer<ExternalPgOwnedFace>&&,
                     structure::PgReferenceCallInfo) noexcept;
    RootSnapshotLocation source_;
    SourceDigest digest_;
    SiLengthUnit units_;
    structure::PgReferenceIdentity patch_;
    core::PlannedBuffer<ExternalPgOwnedFace> faces_;
    structure::PgReferenceCallInfo work_;
};

struct NastranPgBoundCorner {
    std::uint64_t node{};
    NastranCoordinateId cp{0};
    structure::PgPartitionPoint chart;
    std::optional<NastranBulkPosition> cp_field;
    std::array<std::optional<NastranBulkPosition>, 3> tuple_fields;
};

struct NastranPgFaceBound {
    std::uint64_t element{};
    NastranBulkPosition element_field;
    ExternalPgWinding winding = ExternalPgWinding::Forward;
    std::uint8_t count{};
    std::array<NastranPgBoundCorner, 4> corners;
    std::optional<NastranReferenceOffset> nastran_offset;
    structure::ExactRational distance_upper_metres;
    structure::ExactRational normal_chord_upper;

    [[nodiscard]] std::span<const NastranPgBoundCorner> ring() const& noexcept {
        return std::span{corners}.first(count);
    }
    [[nodiscard]] std::span<const NastranPgBoundCorner> ring() const&& = delete;
};

struct ExternalPgBoundLimits {
    ExternalPgLimits region;
    std::uint64_t maximum_directed_work = 0;
};

struct NastranPlanarPgBoundView {
    RootSnapshotLocation source;
    SourceDigest digest;
    SiLengthUnit units;
    structure::PgReferenceIdentity patch;
    std::span<const NastranPgFaceBound> faces;
    structure::PgReferenceCallInfo chart_work;
    std::uint64_t directed_work;
};

// Bounds the same-chart distance and oriented normal chord between complete
// planar Nastran faces and a prepared PG reference graph of degree at most four.
// Geometry bounds are an input to physical conversion acceptance.
class NastranPlanarPgBounds {
public:
    NastranPlanarPgBounds(NastranPlanarPgBounds&&) noexcept = default;
    NastranPlanarPgBounds& operator=(NastranPlanarPgBounds&&) = delete;

    [[nodiscard]] static core::Result<NastranPlanarPgBounds>
    prepare(const structure::PgReferenceGeometry&, const NastranShellCards&, SiLengthUnit,
            std::span<const ExternalPgFace>, const ExternalPgBoundLimits&, core::ResourcePlan&,
            std::optional<ExternalPgDiagnostic>* = nullptr) noexcept;

    [[nodiscard]] std::optional<NastranPlanarPgBoundView> view() const& noexcept;
    [[nodiscard]] std::optional<NastranPlanarPgBoundView> view() const&& = delete;

private:
    NastranPlanarPgBounds(RootSnapshotLocation, SourceDigest, SiLengthUnit,
                          structure::PgReferenceIdentity, core::PlannedBuffer<NastranPgFaceBound>&&,
                          structure::PgReferenceCallInfo, std::uint64_t) noexcept;
    RootSnapshotLocation source_;
    SourceDigest digest_;
    SiLengthUnit units_;
    structure::PgReferenceIdentity patch_;
    core::PlannedBuffer<NastranPgFaceBound> faces_;
    structure::PgReferenceCallInfo chart_work_;
    std::uint64_t directed_work_;
};

struct GmshPgBoundCorner {
    std::uint64_t node{};
    structure::PgPartitionPoint chart;
    std::array<std::uint64_t, 3> source_bytes{};
    std::array<structure::ExactRational, 3> position_metres;
};

struct GmshPgFaceBound {
    std::uint64_t element{};
    GmshKey entity;
    ExternalPgWinding winding = ExternalPgWinding::Forward;
    std::uint8_t count{};
    std::array<GmshPgBoundCorner, 4> corners;
    structure::ExactRational distance_upper_metres;
    structure::ExactRational normal_chord_upper;

    [[nodiscard]] std::span<const GmshPgBoundCorner> ring() const& noexcept {
        return std::span{corners}.first(count);
    }
    [[nodiscard]] std::span<const GmshPgBoundCorner> ring() const&& = delete;
};

struct GmshPgBoundView {
    RootSnapshotLocation source;
    SourceDigest digest;
    SiLengthUnit units;
    structure::PgReferenceIdentity patch;
    std::span<const GmshPgFaceBound> faces;
    structure::PgReferenceCallInfo chart_work;
    std::uint64_t directed_work;
};

// Bounds distance in metres and oriented unit-normal chord on the complete
// original triangle or bilinear quad with the declared PG chart interpolation.
class GmshPgBounds {
public:
    GmshPgBounds(GmshPgBounds&&) noexcept = default;
    GmshPgBounds& operator=(GmshPgBounds&&) = delete;

    [[nodiscard]] static core::Result<GmshPgBounds>
    prepare(const structure::PgReferenceGeometry&, const GmshMesh&, const GmshSurfaceTopology&,
            SiLengthUnit, std::span<const ExternalPgFace>, const ExternalPgBoundLimits&,
            core::ResourcePlan&, std::optional<ExternalPgDiagnostic>* = nullptr) noexcept;

    [[nodiscard]] std::optional<GmshPgBoundView> view() const& noexcept;
    [[nodiscard]] std::optional<GmshPgBoundView> view() const&& = delete;

private:
    GmshPgBounds(RootSnapshotLocation, SourceDigest, SiLengthUnit, structure::PgReferenceIdentity,
                 core::PlannedBuffer<GmshPgFaceBound>&&, structure::PgReferenceCallInfo,
                 std::uint64_t) noexcept;
    RootSnapshotLocation source_;
    SourceDigest digest_;
    SiLengthUnit units_;
    structure::PgReferenceIdentity patch_;
    core::PlannedBuffer<GmshPgFaceBound> faces_;
    structure::PgReferenceCallInfo chart_work_;
    std::uint64_t directed_work_;
};

} // namespace uvlas::io::detail
