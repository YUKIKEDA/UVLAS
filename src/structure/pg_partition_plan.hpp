#pragma once
#include "pg_partition_common.hpp"

namespace uvlas::structure::detail {
struct PgPartitionBudget {
    PgPartitionPlan report;
    std::uint64_t original_scalars = 0, input_bits = 1, borrowed = 0;
    std::uint64_t triangles = 0, named_segments = 0, layers = 0, cells = 0;
    std::uint64_t contact_capacity = 0, provenance_capacity = 0, trace_capacity = 0,
                  scratch_capacity = 0;
    std::uint64_t graph_payload = 0;
};
inline void partition_cap(std::uint64_t count, std::uint64_t limit) {
    (void)partition_index(count);
    if (count > limit)
        partition_fail(core::FailureCode::ResourceExceeded);
}
inline PgPartitionBudget partition_budget(const PgPartitionRequest& request,
                                          const PgPartitionLimits& limits) {
    const auto add = partition_add, mul = partition_mul;
    partition_cancel(limits);
    const auto graph = partition_take(request.graph.view());
    if (limits.algebra.external_memory.valueless_by_exception() ||
        (limits.algebra.timing != ExactTimingMode::Disabled &&
         limits.algebra.timing != ExactTimingMode::Enabled))
        partition_fail(core::FailureCode::InvalidInput);
    if (std::holds_alternative<ExactRequireProvenBound>(limits.algebra.external_memory))
        partition_fail(core::FailureCode::Unsupported);
    for (const auto count :
         {limits.maximum_input_points, limits.maximum_segments, limits.maximum_lines,
          limits.maximum_pieces, limits.maximum_vertices, limits.maximum_edges,
          limits.maximum_layers, limits.maximum_numeric_operations, limits.maximum_topology_visits,
          limits.maximum_text_bytes, limits.algebra.maximum_input_characters,
          limits.algebra.maximum_terms, limits.algebra.maximum_raw_terms,
          limits.algebra.maximum_axis_degree, limits.algebra.maximum_total_degree,
          limits.algebra.maximum_value_bits, limits.algebra.maximum_work_bits})
        (void)partition_index(count);
    if (request.domain.size() < 3 || request.thickness.vertices.size() < 3 ||
        request.layers.empty() ||
        request.thickness.vertices.size() != request.thickness.values.size() ||
        request.thickness.triangles.empty())
        partition_fail(core::FailureCode::InvalidInput);
    if (limits.algebra.maximum_input_characters < 1)
        partition_fail(core::FailureCode::ResourceExceeded);
    PgPartitionBudget p;
    auto& o = p.report;
    p.layers = request.layers.size();
    p.triangles = request.thickness.triangles.size();
    o.input_points = add(request.domain.size(), request.thickness.vertices.size());
    o.segments = add(request.domain.size(), mul(3, p.triangles));
    for (const auto& boundary : request.boundaries) {
        if (boundary.vertices.size() < 2 || boundary.name.empty() ||
            boundary.name.find('\0') != std::string_view::npos)
            partition_fail(core::FailureCode::InvalidInput);
        o.input_points = add(o.input_points, boundary.vertices.size());
        p.named_segments = add(p.named_segments, boundary.vertices.size() - 1);
        o.text_bytes = add(o.text_bytes, boundary.name.size());
    }
    o.input_points = add(o.input_points, mul(2, request.refinements.size()));
    o.segments = add(o.segments, add(p.named_segments, request.refinements.size()));
    p.original_scalars =
        add(mul(2, o.input_points), add(request.thickness.values.size(), mul(2, p.layers)));
    partition_cap(o.input_points, limits.maximum_input_points);
    partition_cap(o.segments, limits.maximum_segments);
    partition_cap(p.layers, limits.maximum_layers);
    partition_cap(o.text_bytes, limits.maximum_text_bytes);
    for (const auto& tri : request.thickness.triangles) {
        for (const auto j : tri)
            if (j < 0 || static_cast<std::uint64_t>(j) >= request.thickness.vertices.size())
                partition_fail(core::FailureCode::InvalidInput);
        if (tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2])
            partition_fail(core::FailureCode::InvalidInput);
    }
    const auto scalar = [&](const ExactRational& x) {
        const auto b = partition_bits(x);
        if (b > limits.algebra.maximum_value_bits)
            partition_fail(core::FailureCode::ResourceExceeded);
        p.input_bits = std::max(p.input_bits, b);
        p.borrowed = add(p.borrowed, partition_retained(x));
    };
    const auto point = [&](const PgPartitionPoint& x) {
        scalar(x.s);
        scalar(x.t);
    };
    for (const auto& x : request.domain)
        point(x);
    for (const auto& x : request.thickness.vertices)
        point(x);
    for (const auto& x : request.thickness.values) {
        scalar(x);
        if (x.sign() <= 0)
            partition_fail(core::FailureCode::InvalidInput);
    }
    for (const auto& x : request.layers) {
        scalar(x.lower);
        scalar(x.upper);
    }
    for (const auto& x : request.boundaries)
        for (const auto& v : x.vertices)
            point(v);
    for (const auto& x : request.refinements) {
        point(x.first);
        point(x.second);
    }
    if (graph.height.terms().size() > limits.algebra.maximum_terms ||
        graph.height.terms().size() > limits.algebra.maximum_raw_terms)
        partition_fail(core::FailureCode::ResourceExceeded);
    for (const auto& term : graph.height.terms())
        if (term.powers.first < 0 || term.powers.second < 0 ||
            static_cast<std::uint64_t>(term.powers.first) > limits.algebra.maximum_axis_degree ||
            static_cast<std::uint64_t>(term.powers.second) > limits.algebra.maximum_axis_degree ||
            partition_add(static_cast<std::uint64_t>(term.powers.first),
                          static_cast<std::uint64_t>(term.powers.second)) >
                limits.algebra.maximum_total_degree)
            partition_fail(core::FailureCode::ResourceExceeded);
    scalar(graph.x0);
    scalar(graph.width);
    scalar(graph.skew);
    for (const auto& term : graph.height.terms())
        scalar(term.coefficient);
    const auto axis = static_cast<std::uint64_t>(2 * graph.n);
    p.cells = mul(axis, axis);
    const auto r =
        add(add(o.segments, mul(2, o.input_points)), 4 * static_cast<std::uint64_t>(graph.n) - 2);
    const auto m = o.line_capacity = add(r, 4);
    const auto c = o.piece_capacity = add(1, mul(r, add(r, 1)) / 2);
    const auto v = o.clip_vertex_capacity = mul(c, m);
    const auto pairs = o.intersection_capacity = mul(m, m - 1) / 2;
    const auto u = o.edge_capacity = mul(2, mul(m, m));
    partition_cap(m, limits.maximum_lines);
    partition_cap(c, limits.maximum_pieces);
    partition_cap(pairs, limits.maximum_vertices);
    partition_cap(u, limits.maximum_edges);
    p.contact_capacity = mul(p.cells, add(pairs, mul(m, m)));
    p.provenance_capacity = mul(u, add(o.segments, p.named_segments));
    p.trace_capacity = mul(u, p.named_segments);
    p.scratch_capacity = std::max({u, pairs, m, c});
    const auto b0 = std::max(p.input_bits, add(partition_log(axis), 2));
    const auto h = mul(2, b0), d = add(mul(2, h), 1), k = add(mul(2, d), 1);
    const auto centroid = add(mul(3, k), 3);
    const auto q = add(mul(u, add(mul(4, k), 2)), 2);
    const auto height = add(add(mul(6, k), mul(45, p.input_bits)), 32);
    o.value_bits = std::max({q, height, centroid});
    o.algebraic_work_bits =
        std::max({add(mul(2, q), 1), add(mul(2, height), 1), add(add(d, mul(2, centroid)), 2),
                  mul(16, add(p.input_bits, 1))});
    const auto w = static_cast<std::uint64_t>(IntegerBackend::limb_bits);
    const auto words =
        add(add(o.algebraic_work_bits / w, o.algebraic_work_bits % w != 0 ? 1 : 0), 2);
    o.planned_work_bits = mul(words, w);
    if (o.value_bits > limits.algebra.maximum_value_bits ||
        o.planned_work_bits > limits.algebra.maximum_work_bits)
        partition_fail(core::FailureCode::ResourceExceeded);
    const auto sort = [&](std::uint64_t x) { return mul(x, partition_log(x)); };
    const auto sort_ops =
        add(add(mul(3, sort(m)), mul(2, sort(pairs))), mul(mul(2, pairs), partition_log(pairs)));
    // The source points are homogenized per source edge instead of cached: 6
    // products supplement the contract's 17 source-line opcodes. Pair admission
    // includes four closed-square comparisons. Both clipped half-polygons
    // recompute endpoint signs; budget 20 rather than 5 calls per original edge.
    // Contact tests consume both endpoints, hence eight coordinate comparisons.
    auto numeric = add(add(p.original_scalars, 18), add(mul(3, o.input_points), mul(3, m)));
    for (const auto term :
         {mul(23, o.segments), mul(27, pairs), mul(52, mul(o.segments, o.segments)),
          mul(20, mul(r, v)), mul(7, mul(c, m)), mul(6, c),
          mul(c, add(add(mul(10, o.segments), mul(24, p.triangles)), mul(4, p.cells))), mul(4, u),
          mul(4, o.input_points), mul(4, c), mul(3, p.cells), mul(3, p.layers), std::uint64_t{8},
          mul(11, mul(v, pairs)), mul(146, u), mul(22, mul(u, add(o.segments, p.named_segments))),
          mul(8, mul(p.cells, add(pairs, u))), sort_ops})
        numeric = add(numeric, term);
    o.numeric_operations = numeric;
    partition_cap(numeric, limits.maximum_numeric_operations);
    // Three name passes (UTF-8 validation, copy, content hash), name comparison,
    // and each length-prefixed integer byte in the canonical source stream.
    const auto scalar_stream = add(17, mul(2, add(p.input_bits, 7) / 8));
    auto topology = add(mul(3, o.text_bytes), mul(request.boundaries.size(), o.text_bytes));
    topology =
        add(topology, add(133, add(mul(p.original_scalars, scalar_stream),
                                   add(mul(24, p.triangles),
                                       add(mul(8, p.layers), add(mul(24, request.boundaries.size()),
                                                                 mul(48, o.segments)))))));
    topology = add(topology, add(mul(p.layers, p.layers),
                                 mul(request.boundaries.size(), request.boundaries.size())));
    topology = add(topology, mul(p.cells, p.provenance_capacity));
    for (const auto term : {mul(o.segments, o.segments),
                            pairs,
                            mul(r, add(mul(3, v), mul(3, c))),
                            mul(c, m),
                            mul(c, add(add(o.segments, p.triangles), p.cells)),
                            mul(v, pairs),
                            mul(2, mul(u, u)),
                            mul(p.cells, add(u, pairs)),
                            mul(u, add(add(o.segments, p.named_segments), p.triangles)),
                            mul(u, add(p.layers, 1)),
                            sort(m),
                            sort(pairs),
                            sort(u),
                            sort(c),
                            o.input_points,
                            o.segments,
                            p.contact_capacity,
                            p.provenance_capacity,
                            mul(4, mul(u, p.layers)),
                            mul(2, p.scratch_capacity),
                            mul(2, v),
                            mul(2, c),
                            p.cells,
                            m})
        topology = add(topology, term);
    o.topology_visits = topology;
    partition_cap(topology, limits.maximum_topology_visits);
    const auto charge = [&](std::uint64_t count, std::size_t bytes) {
        (void)partition_index(count);
        o.owned_payload_bytes = add(o.owned_payload_bytes, mul(count, bytes));
    };
    charge(1, sizeof(PgPartitionData));
    charge(o.input_points, sizeof(PgPartitionPoint));
    charge(o.segments, sizeof(PgPartitionSource));
    charge(p.triangles, sizeof(std::array<core::Index, 3>));
    charge(request.thickness.values.size(), sizeof(ExactRational));
    charge(p.layers, sizeof(PgPartitionOwnedLayer));
    charge(request.boundaries.size(), sizeof(PgPartitionBoundaryRecord));
    charge(o.text_bytes, sizeof(char));
    charge(m, sizeof(PgPartitionLine));
    charge(o.segments, sizeof(core::Index));
    charge(pairs, sizeof(PgPartitionPoint));
    charge(m, sizeof(PgPartitionLine));
    charge(pairs, sizeof(PgPartitionPoint));
    charge(pairs, sizeof(core::Index));
    charge(p.cells, sizeof(PgPartitionCell));
    charge(c, sizeof(PgPartitionPiece));
    charge(u, sizeof(core::Index));
    charge(u, sizeof(PgPartitionEdge));
    charge(p.provenance_capacity, sizeof(core::Index));
    charge(u, sizeof(PgPartitionTriangle));
    charge(mul(u, p.layers), sizeof(PgPartitionVolume));
    charge(mul(mul(2, u), p.layers), sizeof(PgPartitionHorizontalFace));
    charge(mul(u, p.layers), sizeof(PgPartitionLateralFace));
    charge(p.trace_capacity, sizeof(PgPartitionEdgeTrace));
    charge(p.contact_capacity, sizeof(PgPartitionContact));
    charge(mul(2, c), sizeof(PgPartitionClipPiece));
    charge(mul(2, v), sizeof(PgPartitionClipEdge));
    charge(u, sizeof(PgPartitionEdgeOccurrence));
    charge(mul(4, p.scratch_capacity), sizeof(core::Index));
    // Each fresh scalar opcode creates at most one accounted RationalStorage.
    charge(numeric, sizeof(RationalStorage));
    const auto terms = graph.height.terms().size();
    p.graph_payload =
        add(add(sizeof(PgGraphData), mul(p.cells, sizeof(PgCellIndex))),
            add(mul(3, sizeof(RationalStorage)),
                add(mul(terms, add(add(exact_polynomial_raw_term_bytes(), sizeof(ExactTerm)),
                                   sizeof(RationalStorage))),
                    sizeof(PolynomialStorage))));
    o.owned_payload_bytes = add(o.owned_payload_bytes, p.graph_payload);
    const auto j = add(add(add(add(add(p.original_scalars, 18), mul(3, m)), mul(2, pairs)), c),
                       add(p.cells, 6));
    const auto capacity = words <= IntegerBackend::internal_limb_count ? 0 : mul(4, words);
    const auto resize = capacity == 0 ? 0 : words;
    const auto gcd = add(mul(6, add(words, 1)), mul(5, capacity));
    const auto division = mul(2, capacity), karatsuba = words >= 60 ? mul(5, words) : 0;
    o.estimated_external_bytes =
        add(p.borrowed,
            mul(add(add(mul(add(mul(2, j), 32), capacity), std::max({gcd, division, karatsuba})),
                    resize),
                sizeof(boost::multiprecision::limb_type)));
    const auto* admission = std::get_if<ExactEstimateAdmission>(&limits.algebra.external_memory);
    if (admission == nullptr)
        partition_fail(core::FailureCode::InvalidInput);
    if (o.estimated_external_bytes > admission->bytes)
        partition_fail(core::FailureCode::ResourceExceeded);
    partition_cancel(limits);
    return p;
}
} // namespace uvlas::structure::detail
