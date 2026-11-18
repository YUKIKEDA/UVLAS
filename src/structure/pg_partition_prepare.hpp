#pragma once
#include "pg_partition_hash.hpp"

namespace uvlas::structure::detail {
template <class T>
core::PlannedBuffer<T> partition_allocate(core::ResourcePlan& plan, std::uint64_t count) {
    return partition_take(plan.allocate<T>(partition_index(count), partition_operation));
}
inline PgIntegrationPartition partition_prepare(const PgPartitionRequest& request,
                                                const PgPartitionLimits& limits,
                                                core::Workspace& workspace,
                                                core::ResourcePlan& plan) {
    auto lease = partition_take(workspace.acquire());
    const auto before = plan.allocated_bytes();
    auto budget = partition_budget(request, limits);
    if (budget.report.owned_payload_bytes > plan.remaining_bytes())
        partition_fail(core::FailureCode::ResourceExceeded);
    auto storage = partition_allocate<PgPartitionData>(plan, 1);
    partition_cancel(limits);
    auto& d = storage.values()[0];
    const auto& p = budget.report;
    d.source_points = partition_allocate<PgPartitionPoint>(plan, p.input_points);
    d.sources = partition_allocate<PgPartitionSource>(plan, p.segments);
    d.thickness_triangles = partition_allocate<std::array<core::Index, 3>>(plan, budget.triangles);
    d.thickness_values = partition_allocate<ExactRational>(plan, request.thickness.values.size());
    d.layers = partition_allocate<PgPartitionOwnedLayer>(plan, budget.layers);
    d.boundaries = partition_allocate<PgPartitionBoundaryRecord>(plan, request.boundaries.size());
    d.text = partition_allocate<char>(plan, p.text_bytes);
    d.lines = partition_allocate<PgPartitionLine>(plan, p.line_capacity);
    d.source_lines = partition_allocate<core::Index>(plan, p.segments);
    d.vertices = partition_allocate<PgPartitionPoint>(plan, p.intersection_capacity);
    d.pairs = partition_allocate<core::Index>(plan, p.intersection_capacity);
    d.cells = partition_allocate<PgPartitionCell>(plan, budget.cells);
    d.pieces = partition_allocate<PgPartitionPiece>(plan, p.piece_capacity);
    d.piece_vertices = partition_allocate<core::Index>(plan, p.edge_capacity);
    d.edges = partition_allocate<PgPartitionEdge>(plan, p.edge_capacity);
    d.edge_sources = partition_allocate<core::Index>(plan, budget.provenance_capacity);
    d.triangles = partition_allocate<PgPartitionTriangle>(plan, p.edge_capacity);
    d.volumes =
        partition_allocate<PgPartitionVolume>(plan, partition_mul(p.edge_capacity, budget.layers));
    d.horizontal = partition_allocate<PgPartitionHorizontalFace>(
        plan, partition_mul(partition_mul(2, p.edge_capacity), budget.layers));
    d.lateral = partition_allocate<PgPartitionLateralFace>(
        plan, partition_mul(p.edge_capacity, budget.layers));
    d.traces = partition_allocate<PgPartitionEdgeTrace>(plan, budget.trace_capacity);
    d.contacts = partition_allocate<PgPartitionContact>(plan, budget.contact_capacity);
    PgPartitionScratch scratch;
    scratch.raw_lines = partition_allocate<PgPartitionLine>(plan, p.line_capacity);
    scratch.raw_vertices = partition_allocate<PgPartitionPoint>(plan, p.intersection_capacity);
    for (auto& array : scratch.pieces)
        array = partition_allocate<PgPartitionClipPiece>(plan, p.piece_capacity);
    for (auto& array : scratch.clip)
        array = partition_allocate<PgPartitionClipEdge>(plan, p.clip_vertex_capacity);
    scratch.occurrences = partition_allocate<PgPartitionEdgeOccurrence>(plan, p.edge_capacity);
    scratch.order = partition_allocate<core::Index>(plan, budget.scratch_capacity);
    scratch.merge = partition_allocate<core::Index>(plan, budget.scratch_capacity);
    scratch.remap = partition_allocate<core::Index>(plan, budget.scratch_capacity);
    scratch.queue = partition_allocate<core::Index>(plan, budget.scratch_capacity);
    PgPartitionAlgebra algebra{limits, budget, plan};
    PgPartitionGeometry geometry{d, scratch, algebra, request};
    const auto measured = limits.algebra.timing == ExactTimingMode::Enabled;
    core::StageRecorder input_timing{measured, partition_operation};
    core::StageRecorder intersection_timing{measured, partition_operation};
    core::StageRecorder planar_timing{measured, partition_operation};
    core::StageRecorder trace_timing{measured, partition_operation};
    {
        const core::StageScope phase{input_timing, core::MeasurementStage::Preparation};
        geometry.snapshot();
        geometry.validate_geometry();
    }
    {
        const core::StageScope phase{intersection_timing, core::MeasurementStage::Preparation};
        geometry.make_lines();
        geometry.make_vertices();
    }
    {
        const core::StageScope phase{planar_timing, core::MeasurementStage::Preparation};
        geometry.arrange();
        geometry.retain_pieces();
    }
    {
        const core::StageScope phase{trace_timing, core::MeasurementStage::Preparation};
        PgPartitionTopology topology{geometry};
        topology.make_edges();
        topology.components();
        topology.make_domains();
        topology.contacts();
        PgPartitionHash content{algebra, {}};
        content(d);
    }
    partition_cancel(limits);
    if (plan.allocated_bytes() - before > p.owned_payload_bytes ||
        algebra.observed_bits > p.value_bits)
        partition_fail(core::FailureCode::ResourceExceeded);
    auto report = p;
    report.allocated_owned_payload_bytes = plan.allocated_bytes() - before;
    report.observed_numeric_operations = algebra.numeric_used;
    report.observed_topology_visits = algebra.topology_used;
    report.observed_value_bits = algebra.observed_bits;
    if (measured) {
        const auto sample = [](const core::StageRecorder& recorder) {
            return partition_take(
                recorder.snapshot())[static_cast<std::size_t>(core::MeasurementStage::Preparation)];
        };
        report.timings = PgPartitionTimings{{},
                                            sample(input_timing),
                                            sample(intersection_timing),
                                            sample(planar_timing),
                                            sample(trace_timing)};
    }
    return PgPartitionFactory::make(std::move(storage), report);
}
} // namespace uvlas::structure::detail
