#pragma once
#include "exact_polynomial_detail.hpp"
#include "pg_graph_detail.hpp"

#include <uvlas/structure/pg_integration_partition.hpp>

namespace uvlas::structure::detail {
struct PgPartitionLine {
    std::array<ExactRational, 3> coefficients;
};
struct PgPartitionClipPiece {
    std::uint64_t first = 0, count = 0;
};
struct PgPartitionClipEdge {
    core::Index vertex = 0;
    core::Index line = 0;
};
struct PgPartitionEdgeOccurrence {
    core::Index first = 0, second = 0, piece = 0;
    bool forward = true;
};
struct PgPartitionData {
    std::optional<PgGraph> graph;
    core::PlannedBuffer<PgPartitionPoint> source_points;
    core::PlannedBuffer<PgPartitionSource> sources;
    core::PlannedBuffer<std::array<core::Index, 3>> thickness_triangles;
    core::PlannedBuffer<ExactRational> thickness_values;
    core::PlannedBuffer<PgPartitionOwnedLayer> layers;
    core::PlannedBuffer<PgPartitionBoundaryRecord> boundaries;
    core::PlannedBuffer<char> text;
    std::uint64_t domain_count = 0, thickness_offset = 0;
    core::PlannedBuffer<PgPartitionLine> lines;
    core::PlannedBuffer<core::Index> source_lines;
    core::PlannedBuffer<PgPartitionPoint> vertices;
    core::PlannedBuffer<core::Index> pairs;
    core::PlannedBuffer<PgPartitionCell> cells;
    core::PlannedBuffer<PgPartitionPiece> pieces;
    core::PlannedBuffer<core::Index> piece_vertices;
    core::PlannedBuffer<PgPartitionEdge> edges;
    core::PlannedBuffer<core::Index> edge_sources;
    core::PlannedBuffer<PgPartitionTriangle> triangles;
    core::PlannedBuffer<PgPartitionVolume> volumes;
    core::PlannedBuffer<PgPartitionHorizontalFace> horizontal;
    core::PlannedBuffer<PgPartitionLateralFace> lateral;
    core::PlannedBuffer<PgPartitionEdgeTrace> traces;
    core::PlannedBuffer<PgPartitionContact> contacts;
    PgPartitionIdentity identity{};
    std::uint64_t line_count = 0, vertex_count = 0, piece_count = 0, piece_vertex_count = 0;
    std::uint64_t edge_count = 0, edge_source_count = 0, triangle_count = 0;
    std::uint64_t volume_count = 0, horizontal_count = 0, lateral_count = 0;
    std::uint64_t trace_count = 0, contact_count = 0;
};
struct PgPartitionFactory {
    static const PgPartitionData* data(const PgIntegrationPartition&) noexcept;
    static void info(PgIntegrationPartition&, PgPartitionPlan) noexcept;
    static PgIntegrationPartition make(core::PlannedBuffer<PgPartitionData>,
                                       PgPartitionPlan) noexcept;
};
struct PgPartitionScratch {
    core::PlannedBuffer<PgPartitionLine> raw_lines;
    core::PlannedBuffer<PgPartitionPoint> raw_vertices;
    std::array<core::PlannedBuffer<PgPartitionClipPiece>, 2> pieces;
    std::array<core::PlannedBuffer<PgPartitionClipEdge>, 2> clip;
    core::PlannedBuffer<PgPartitionEdgeOccurrence> occurrences;
    core::PlannedBuffer<core::Index> order, merge, remap, queue;
};
} // namespace uvlas::structure::detail
