#pragma once

#include "external_mesh_sidecar.hpp"
#include "nastran_laminate_cards.hpp"

#include <uvlas/structure/pg_reference_geometry.hpp>

namespace uvlas::io::detail {

enum class NastranPcompGeometryFailure : std::uint8_t {
    Source,
    Mode,
    LayerCount,
    ReferencePlane,
    FullThickness,
    Interface,
    Resource
};

struct NastranPcompGeometryDiagnostic {
    NastranPcompGeometryFailure reason = NastranPcompGeometryFailure::Source;
    std::optional<NastranBulkPosition> field{};
    std::optional<core::LayerId> case_layer{};
    std::optional<core::Index> thickness_vertex{};
};

struct NastranPcompGeometryLimits {
    NastranPcompGeometryLimits(const structure::ExactLimits& configured, std::uint64_t vertices,
                               std::uint64_t operations) noexcept
        : exact(configured), thickness_vertices(vertices), work(operations) {}

    structure::ExactLimits exact;
    std::uint64_t thickness_vertices;
    std::uint64_t work;
};

struct NastranPcompLayerMatch {
    core::LayerId layer{0};
    std::size_t case_ordinal{};
    std::size_t declared_ply{};
    NastranMaterialReference material;
    NastranBulkPosition thickness;
    std::optional<NastranShellReal> source_angle{};
    structure::ExactRational ply_thickness_metres;
};

struct NastranPcompLayerSection {
    RootSnapshotLocation source;
    SourceDigest digest;
    NastranPropertyId property{0};
    NastranBulkPosition property_field;
    std::optional<NastranBulkPosition> bottom_offset_field{};
    structure::PgReferenceIdentity pg;
    MeshSidecarOrientation winding = MeshSidecarOrientation::Same;
    structure::ExactRational full_thickness_metres;
    core::PlannedBuffer<NastranPcompLayerMatch> layers; // Case order.
    std::uint64_t work{};
};

// The layer-section result verifies the declared ply coordinates against the Case
// layer intervals. Full thickness must still be checked over the selected element regions.
[[nodiscard]] core::Result<NastranPcompLayerSection> compare_nastran_pcomp_layer_section(
    const NastranLaminateCards&, std::size_t property_ordinal, SiLengthUnit, MeshSidecarOrientation,
    const structure::PgReferenceGeometry&, const NastranPcompGeometryLimits&, core::ResourcePlan&,
    const std::stop_token&, std::optional<NastranPcompGeometryDiagnostic>* = nullptr) noexcept;

struct NastranPcompGeometryMatch : NastranPcompLayerSection {};

// Compare an adopted PCOMP stack after its signed element reference plane and whole
// region correspondence to the specified prepared PG have been certified. The original
// winding determines the Case order. This certifies thickness and interfaces, and hands
// off source-ply material identities for comparison to the Case material declarations.
[[nodiscard]] core::Result<NastranPcompGeometryMatch> compare_nastran_pcomp_geometry(
    const NastranLaminateCards&, std::size_t property_ordinal, SiLengthUnit, MeshSidecarOrientation,
    const structure::PgReferenceGeometry&, const NastranPcompGeometryLimits&, core::ResourcePlan&,
    const std::stop_token&, std::optional<NastranPcompGeometryDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::io::detail
