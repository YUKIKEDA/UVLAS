#pragma once

#include "external_planar_pg.hpp"
#include "nastran_laminate_cards.hpp"

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
    structure::ExactRational ply_thickness_metres;
};

struct NastranPcompGeometryMatch {
    RootSnapshotLocation source;
    SourceDigest digest;
    structure::PgReferenceIdentity pg;
    ExternalPgWinding winding = ExternalPgWinding::Forward;
    structure::ExactRational full_thickness_metres;
    core::PlannedBuffer<NastranPcompLayerMatch> layers; // Case order.
    std::uint64_t work{};
};

// Compare an adopted PCOMP stack after its signed element reference plane and whole
// region correspondence to the specified prepared PG have been certified. The original
// winding determines the Case order. This certifies thickness and interfaces, and hands
// off source-ply material identities for comparison to the Case material declarations.
[[nodiscard]] core::Result<NastranPcompGeometryMatch> compare_nastran_pcomp_geometry(
    const NastranLaminateCards&, std::size_t property_ordinal, SiLengthUnit, ExternalPgWinding,
    const structure::PgReferenceGeometry&, const NastranPcompGeometryLimits&, core::ResourcePlan&,
    const std::stop_token&, std::optional<NastranPcompGeometryDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::io::detail
