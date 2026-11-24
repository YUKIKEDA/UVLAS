#include "pg_conditioning_test_support.hpp"

#include <algorithm>
#include <string>
#include <uvlas/structure/pg_flat_gram.hpp>
#include <uvlas/structure/pg_full_trace.hpp>
#include <uvlas/structure/pg_global_constraints.hpp>
#include <uvlas/structure/pg_sparse_chart.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::structure;
using namespace partition_test;
using namespace pg_conditioning_test;

PgPreparedPatch patch(Fixture& f, const PgReferenceFrame& frame, core::Workspace& ws,
                      core::ResourcePlan& plan) {
    const auto limits = query_limits();
    auto basis = make_basis(f, frame, ws, plan);
    auto gram = take(PgFlatGram::prepare(std::move(basis), 256, {limits, steps}, ws, plan));
    auto accepted =
        take(std::move(gram).generate_and_certify({320, 246}, {profile, steps, {}}, ws, plan));
    return take(PgPreparedPatch::prepare(std::move(accepted), limits, ws, plan));
}

const PgIntegrationPartition& partition(const PgPreparedPatch& patch) {
    return take(take(patch.conditioning().reference_basis())->source()).partition;
}

struct EdgeRange {
    core::Index edge{};
    ExactRational first, last;
    bool increasing = true;
};

std::vector<EdgeRange> vertical_edges(Fixture& fixture, const PgPreparedPatch& patch, const Q& s,
                                      const Q& scale) {
    const auto& part = partition(patch);
    const auto vertices = take(part.vertices());
    const auto edges = take(part.edges());
    std::vector<EdgeRange> result;
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const auto& edge = edges[i];
        if (!edge.external)
            continue;
        const auto& a = vertices[static_cast<std::size_t>(edge.vertices[0])];
        const auto& b = vertices[static_cast<std::size_t>(edge.vertices[1])];
        if (q(a.s) != s || q(b.s) != s)
            continue;
        const Q x = q(a.t) * scale, y = q(b.t) * scale;
        const auto lo = std::min(x, y), hi = std::max(x, y);
        result.push_back({static_cast<core::Index>(i), rat(fixture, numerator(lo), denominator(lo)),
                          rat(fixture, numerator(hi), denominator(hi)), x < y});
    }
    std::ranges::sort(result, [](const auto& a, const auto& b) { return q(a.first) < q(b.first); });
    return result;
}

std::vector<PgTraceEdgeSpan> edge_spans(std::span<const EdgeRange> ranges,
                                        const ExactRational& zero, const ExactRational& one,
                                        std::optional<core::Index> layer = {}) {
    std::vector<PgTraceEdgeSpan> result;
    result.reserve(ranges.size());
    for (const auto& edge : ranges)
        result.push_back({edge.edge, edge.first, edge.last, edge.increasing ? zero : one,
                          edge.increasing ? one : zero, layer});
    return result;
}

std::vector<PgCartesianExactPosition> states(Fixture& f, std::size_t count, std::size_t offset) {
    std::vector<PgCartesianExactPosition> output(count);
    for (std::size_t i = 0; i < output.size(); ++i)
        for (std::size_t k = 0; k < output[i].size(); ++k)
            output[i][k] =
                f.r(std::to_string(static_cast<int>((4 * i + 2 * k + offset) % 11) - 5), "13");
    return output;
}

std::array<Q, 3> direct_trace(Fixture& fixture, const PgPreparedPatch& p,
                              std::span<const EdgeRange> edges, core::Index ordinal, const Q& eta,
                              const Q& z, const Q& half_thickness,
                              std::span<const PgCartesianExactPosition> state,
                              core::ResourcePlan& plan) {
    const auto e =
        std::ranges::find_if(edges, [&](const auto& edge) { return edge.edge == ordinal; });
    REQUIRE(e != edges.end());
    const Q a = (eta - q(e->first)) / (q(e->last) - q(e->first));
    const Q zeta = z / half_thickness;
    const auto& part = partition(p);
    const auto faces = take(part.lateral_faces());
    const auto layers = take(part.layers());
    const auto f = std::ranges::find_if(faces, [&](const auto& face) {
        const auto& layer = layers[static_cast<std::size_t>(face.layer)];
        return face.edge == ordinal && q(layer.lower) < zeta && zeta < q(layer.upper);
    });
    REQUIRE(f != faces.end());
    const auto& layer = layers[static_cast<std::size_t>(f->layer)];
    const Q fraction = f->canonical_direction == e->increasing ? a : 1 - a;
    const Q depth = (zeta - q(layer.lower)) / (q(layer.upper) - q(layer.lower));
    const auto along = rat(fixture, numerator(fraction), denominator(fraction));
    const auto thickness = rat(fixture, numerator(depth), denominator(depth));
    const PgReferencePoint request =
        PgReferenceLateralPoint{static_cast<core::Index>(f - faces.begin()), along, thickness};
    const auto values = take(p.evaluate_reference(std::span{&request, 1}, query_limits(), plan));
    REQUIRE(q(values.points.values().front().g) == 1);
    CHECK(q(values.points.values().front().position_metres[1].rational) == eta);
    std::array<Q, 3> result{};
    for (std::size_t j = 0; j < state.size(); ++j) {
        const auto& field = values.columns.values()[j].value;
        const Q n = q(field.rational) + q(field.radical);
        for (std::size_t k = 0; k < result.size(); ++k)
            result[k] += n * q(state[j][k]);
    }
    return result;
}

Q work(std::span<const PgCartesianExactPosition> a, std::span<const PgCartesianExactPosition> b) {
    REQUIRE(a.size() == b.size());
    Q result;
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a[i].size(); ++j)
            result += q(a[i][j]) * q(b[i][j]);
    return result;
}

} // namespace

TEST_CASE("PG nonconforming full trace owns coefficient rows and the original dual",
          "[SH-D05][SH-D06][N4-03][N4-08][verification]") {
    Fixture left, right{2};
    const auto zero = left.r("0"), one = left.r("1"), right_length = left.r("3", "2");
    const auto right_width = left.r("2", "3");
    const ExactChartPolynomial height;
    right.graph =
        take(PgGraph::prepare({zero, right_width, zero, height, 2, PgGraphConvention::Version1},
                              {exact_limits(), 64, {}}, right.input_plan));
    right.domain[2].t = right.r("2", "3");
    right.domain[3].t = right.r("2", "3");
    for (auto& h : left.heights)
        h = left.r("2");
    std::array<ExactRational, 9> rotation;
    for (const auto i : {0U, 4U, 8U})
        rotation[i] = left.r("1");
    const std::array<ExactRational, 3> left_origin{};
    std::array<ExactRational, 3> right_origin;
    right_origin[0] = left.r("1");
    core::Workspace ws{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan plan{profile};
    auto a = patch(left, {rotation, left_origin, one}, ws, plan);
    auto b = patch(right, {rotation, right_origin, right_length}, ws, plan);
    const auto a_edges = vertical_edges(left, a, Q{1}, Q{1});
    const auto b_edges = vertical_edges(left, b, Q{0}, Q("3/2"));
    REQUIRE(a_edges.size() == 2);
    REQUIRE(b_edges.size() == 3);
    const auto as = edge_spans(a_edges, zero, one), bs = edge_spans(b_edges, zero, one);
    PgTraceLimits limits{query_limits(), 100'000};
    limits.arithmetic.maximum_scalar_operations = 200'000'000;
    limits.arithmetic.maximum_partition_visits = 100'000'000;
    const auto verify = [&](std::span<const PgTraceEdgeSpan> le,
                            std::span<const PgTraceEdgeSpan> re, const PgTraceLimits& bound,
                            core::ResourcePlan& space) {
        return PgFullTrace::prepare(a, le, b, re, bound, ws, space);
    };
    core::ResourcePlan none{0};
    auto cancelled = limits;
    cancelled.arithmetic.cancellation = requested_stop_token();
    check_failure(verify(as, bs, cancelled, none), core::FailureCode::Cancelled);
    check_failure(verify(as, bs, limits, none), core::FailureCode::ResourceExceeded);
    CHECK(none.allocated_bytes() == 0);
    const auto none_left = edge_spans(a_edges, zero, one, core::Index{0});
    const auto none_right = edge_spans(b_edges, zero, one, core::Index{1});
    check_failure(verify(none_left, none_right, limits, plan), core::FailureCode::GeometryInvalid);
    REQUIRE(a.info());
    REQUIRE(b.info());

    const auto partial_first = left.r("1", "8"), partial_last = left.r("1", "4");
    // The physical endpoints are inside the first atomic edge of each unlike partition.
    const auto along_at = [&](const PgTraceEdgeSpan& span, const ExactRational& eta) {
        const Q fraction = (q(eta) - q(span.first)) / (q(span.last) - q(span.first));
        const Q value = q(span.along_first) + fraction * (q(span.along_last) - q(span.along_first));
        return rat(left, numerator(value), denominator(value));
    };
    const auto a_first = along_at(as.front(), partial_first);
    const auto a_last = along_at(as.front(), partial_last);
    const auto b_first = along_at(bs.front(), partial_first);
    const auto b_last = along_at(bs.front(), partial_last);
    const std::array partial_a{PgTraceEdgeSpan{as.front().edge, partial_first, partial_last,
                                               a_first, a_last, core::Index{0}}};
    const std::array partial_b{PgTraceEdgeSpan{bs.front().edge, partial_first, partial_last,
                                               b_first, b_last, core::Index{0}}};
    auto partial = take(verify(partial_a, partial_b, limits, plan));
    REQUIRE(partial.overlaps().size() == 1);
    const auto& clip = partial.overlaps().front();
    CHECK(q(clip.first) == Q("1/8"));
    CHECK(q(clip.last) == Q("1/4"));
    CHECK(clip.left_layer == 0);
    CHECK(clip.right_layer == 0);
    CHECK(q(clip.left_depth.lower_first) == -1);
    CHECK(q(clip.right_depth.lower_first) == Q("-1/2"));
    CHECK(q(clip.left_depth.upper_last) == Q("-1/3"));
    CHECK(q(clip.right_depth.upper_last) == Q("-1/6"));
    check_failure(PgFullTrace::prepare(a, partial_a, a, partial_a, limits, ws, plan),
                  core::FailureCode::GeometryInvalid);

    auto connection = take(verify(as, bs, limits, plan));
    const auto info = take(connection.info());
    REQUIRE(info.overlaps == 4);
    CHECK(info.left != info.right);
    CHECK(info.unique_rows > 0);
    CHECK(info.unique_rows < info.original_rows);
    CHECK(info.nonzeros > info.unique_rows);
    const std::array<Q, 5> breakpoints{Q{0}, Q("3/8"), Q("1/2"), Q("3/4"), Q{1}};
    for (std::size_t i = 0; i < info.overlaps; ++i) {
        const auto& overlap = connection.overlaps()[i];
        CHECK(q(overlap.first) == breakpoints[i]);
        CHECK(q(overlap.last) == breakpoints[i + 1]);
        CHECK(q(overlap.active_first) == 0);
        CHECK(q(overlap.active_last) == 1);
        CHECK(q(overlap.left_depth.lower_first) == -1);
        CHECK(q(overlap.right_depth.lower_first) == Q("-1/2"));
        CHECK(q(overlap.left_depth.upper_last) == 1);
        CHECK(q(overlap.right_depth.upper_last) == Q("1/2"));
    }
    const auto ca = states(left, a.reference_coordinates().size(), 1);
    const auto cb = states(left, b.reference_coordinates().size(), 6);
    const PgTraceCoordinates input{info.left, info.right, ca, cb};
    auto stale = input;
    stale.left = info.right;
    check_failure(connection.apply(stale, limits.arithmetic, ws, none),
                  core::FailureCode::StaleGeneration);
    const PgTraceCoordinates q0{info.left, info.right, a.reference_coordinates(),
                                b.reference_coordinates()};
    const auto zero_rows = take(connection.apply(q0, limits.arithmetic, ws, plan));
    for (const auto& row : zero_rows.values())
        for (const auto& x : row)
            CHECK(x.sign() == 0);
    const auto residual = take(connection.apply(input, limits.arithmetic, ws, plan));
    const Q eta_local("1/3");
    for (std::size_t i = 0; i < connection.overlaps().size(); ++i) {
        const auto& overlap = connection.overlaps()[i];
        const Q eta = q(overlap.first) + eta_local * (q(overlap.last) - q(overlap.first));
        for (const Q& z : std::array<Q, 3>{Q("-1/4"), Q{0}, Q("1/4")}) {
            auto actual = direct_trace(left, a, a_edges, overlap.left_edge, eta, z, Q{1}, ca, plan);
            const auto right_value =
                direct_trace(left, b, b_edges, overlap.right_edge, eta, z, Q("1/2"), cb, plan);
            std::array<Q, 3> from_rows{};
            for (std::size_t row = 0; row < connection.original_rows().size(); ++row) {
                const auto& key = connection.original_rows()[row];
                if (static_cast<std::size_t>(key.overlap) != i)
                    continue;
                const auto coefficient =
                    power(eta_local, key.power) * power(z, static_cast<core::Index>(key.field));
                for (std::size_t k = 0; k < from_rows.size(); ++k)
                    from_rows[k] += coefficient * q(residual.values()[row][k]);
            }
            for (std::size_t k = 0; k < actual.size(); ++k)
                actual[k] -= right_value[k];
            CHECK(from_rows == actual);
        }
    }

    const auto lambda = states(left, info.original_rows, 2);
    const auto force = take(connection.transpose(lambda, limits.arithmetic, ws, plan));
    CHECK(work(residual.values(), lambda) ==
          work(force.values().first(ca.size()), ca) + work(force.values().subspan(ca.size()), cb));
    const auto unique = states(left, info.unique_rows, 3);
    const auto restored = take(connection.restore_dual(unique, limits.arithmetic, ws, plan));
    const auto equivalent =
        take(connection.transpose(restored.values(), limits.arithmetic, ws, plan));
    const Q original_work = work(equivalent.values().first(ca.size()), ca) +
                            work(equivalent.values().subspan(ca.size()), cb);
    Q normalized_work;
    for (std::size_t row = 0; row < unique.size(); ++row)
        for (const auto& term : take(connection.unique_row(static_cast<core::Index>(row)))) {
            const auto j = static_cast<std::size_t>(term.column);
            const auto& values = j < ca.size() ? ca[j] : cb[j - ca.size()];
            for (std::size_t k = 0; k < values.size(); ++k)
                normalized_work += q(term.coefficient) * q(values[k]) * q(unique[row][k]);
        }
    CHECK(original_work == normalized_work);

    const auto only_first = std::span{as}.first(1);
    check_failure(PgFullTrace::prepare_support(a, only_first, info.right, ca, limits, ws, none),
                  core::FailureCode::StaleGeneration);
    const auto split_eta = left.r("1", "4");
    const auto split_along = along_at(as.front(), split_eta);
    const std::array support_spans{
        PgTraceEdgeSpan{
            as.front().edge, as.front().first, split_eta, as.front().along_first, split_along, {}},
        PgTraceEdgeSpan{
            as.front().edge, split_eta, as.front().last, split_along, as.front().along_last, {}}};
    const auto support =
        take(PgFullTrace::prepare_support(a, support_spans, info.left, ca, limits, ws, plan));
    const auto fixed = take(support.info());
    REQUIRE(fixed.overlaps == 2);
    const auto bending_slope = [](const PgTraceOriginalRow& row) {
        return row.field == PgCoefficientField::C && row.power == 1;
    };
    const auto raw = support.original_rows();
    const auto first_slope = std::ranges::find_if(raw, bending_slope);
    REQUIRE(first_slope != raw.end());
    const auto second_slope = std::ranges::find_if(first_slope + 1, raw.end(), bending_slope);
    REQUIRE(second_slope != raw.end());
    CHECK(first_slope->overlap == 0);
    CHECK(second_slope->overlap == 1);
    REQUIRE(first_slope->unique_row >= 0);
    CHECK(first_slope->unique_row == second_slope->unique_row);
    CHECK(q(first_slope->scale) == q(second_slope->scale));
    const auto normalized_slope = take(support.unique_row(first_slope->unique_row));
    REQUIRE_FALSE(normalized_slope.empty());
    CHECK(q(normalized_slope.front().coefficient) == 1);
    CHECK(fixed.kind == PgTraceKind::Support);
    CHECK(fixed.right_columns == 0);
    CHECK(std::ranges::any_of(support.original_rows(), [](const auto& row) {
        return std::ranges::any_of(row.prescribed, [](const auto& x) { return x.sign() != 0; });
    }));
    auto exactly_prescribed =
        take(support.residual({fixed.left, fixed.right, ca, {}}, limits.arithmetic, ws, plan));
    for (const auto& row : exactly_prescribed.values())
        for (const auto& x : row)
            CHECK(x.sign() == 0);
    const auto mismatch = take(support.residual(
        {fixed.left, fixed.right, a.reference_coordinates(), {}}, limits.arithmetic, ws, plan));
    CHECK(std::ranges::any_of(mismatch.values(), [](const auto& row) {
        return std::ranges::any_of(row, [](const auto& x) { return x.sign() != 0; });
    }));
}

TEST_CASE("PG global constraints retain original rows, right-hand sides and exact work",
          "[SH-D05][SH-D06][N4-03][N4-04][N4-08][verification]") {
    Fixture f, fb, fc;
    const auto zero = f.r("0"), one = f.r("1");
    std::array<ExactRational, 9> rotation;
    for (const auto i : {0U, 4U, 8U})
        rotation[i] = f.r("1");
    const std::array<ExactRational, 3> oa{};
    std::array<ExactRational, 3> ob, oc;
    ob[0] = f.r("1");
    oc[0] = f.r("2");
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan plan{profile};
    auto a = patch(f, {rotation, oa, one}, workspace, plan);
    auto b = patch(fb, {rotation, ob, one}, workspace, plan);
    auto c = patch(fc, {rotation, oc, one}, workspace, plan);
    const auto ai = take(a.info()).identity, bi = take(b.info()).identity,
               ci = take(c.info()).identity;

    const auto a_right = vertical_edges(f, a, Q{1}, Q{1});
    const auto b_left = vertical_edges(f, b, Q{0}, Q{1});
    const auto b_right = vertical_edges(f, b, Q{1}, Q{1});
    const auto c_left = vertical_edges(f, c, Q{0}, Q{1});
    const auto ar = edge_spans(a_right, zero, one), bl = edge_spans(b_left, zero, one),
               br = edge_spans(b_right, zero, one), cl = edge_spans(c_left, zero, one);
    PgGlobalConstraintLimits limits{query_limits(), 3, 4, 8192, 100'000};
    limits.arithmetic.maximum_scalar_operations = 300'000'000;
    limits.arithmetic.maximum_partition_visits = 100'000'000;
    const PgTraceLimits trace_limits{limits.arithmetic, limits.maximum_nonzeros};
    auto ab = take(PgFullTrace::prepare(a, ar, b, bl, trace_limits, workspace, plan));
    auto bc = take(PgFullTrace::prepare(b, br, c, cl, trace_limits, workspace, plan));
    const auto prescribed = states(f, c.reference_coordinates().size(), 8);
    auto fixed = take(PgFullTrace::prepare_support(c, std::span{cl}.first(1), ci, prescribed,
                                                   trace_limits, workspace, plan));
    const std::array patches{PgConstraintPatchBinding{numerics::PatchId{72}, b},
                             PgConstraintPatchBinding{numerics::PatchId{45}, c},
                             PgConstraintPatchBinding{numerics::PatchId{17}, a}};
    const std::array bindings{PgConstraintTraceBinding{PgConstraintId{240}, fixed, 1, {}},
                              PgConstraintTraceBinding{PgConstraintId{118}, ab, 2, core::Index{0}},
                              PgConstraintTraceBinding{PgConstraintId{315}, bc, 0, core::Index{1}},
                              PgConstraintTraceBinding{PgConstraintId{660}, ab, 2, core::Index{0}}};
    const auto prepare = [&](const auto& p, const auto& t, const PgGlobalConstraintLimits& bound,
                             core::ResourcePlan& budget) {
        return PgGlobalConstraints::prepare(p, t, bound, workspace, budget);
    };
    core::ResourcePlan none{0};
    auto cancelled = limits;
    cancelled.arithmetic.cancellation = requested_stop_token();
    check_failure(prepare(patches, bindings, cancelled, none), core::FailureCode::Cancelled);
    check_failure(prepare(patches, bindings, limits, none), core::FailureCode::ResourceExceeded);
    CHECK(none.allocated_bytes() == 0);
    const std::array duplicates{patches[0], PgConstraintPatchBinding{patches[0].id, c}, patches[2]};
    check_failure(prepare(duplicates, bindings, limits, plan), core::FailureCode::InvalidInput);
    const std::array reversed{PgConstraintTraceBinding{PgConstraintId{118}, ab, 0, core::Index{2}}};
    check_failure(prepare(patches, reversed, limits, plan), core::FailureCode::StaleGeneration);
    const std::array bad_support{
        PgConstraintTraceBinding{PgConstraintId{240}, fixed, 1, core::Index{0}}};
    check_failure(prepare(patches, bad_support, limits, plan), core::FailureCode::InvalidInput);
    const std::array repeated_id{bindings[1], bindings[1]};
    check_failure(prepare(patches, repeated_id, limits, plan), core::FailureCode::InvalidInput);

    auto global = take(prepare(patches, bindings, limits, plan));
    const auto layout = global.layout();
    const auto sources = global.sources();
    const auto info = take(global.info());
    REQUIRE(layout.size() == 3);
    REQUIRE(sources.size() == bindings.size());
    CHECK(layout[0].first == 0);
    CHECK(layout[1].first == layout[0].count);
    CHECK(layout[2].first == layout[0].count + layout[1].count);
    CHECK(info.columns == layout[2].first + layout[2].count);
    CHECK(layout[0].id == patches[0].id);
    CHECK(layout[1].basis == ci);
    CHECK(layout[2].basis == ai);
    CHECK(sources[0].left == patches[1].id);
    CHECK_FALSE(sources[0].right.has_value());
    CHECK(sources[1].left == patches[2].id);
    CHECK(sources[1].right == patches[0].id);
    CHECK(sources[2].left == patches[0].id);
    CHECK(sources[2].right == patches[1].id);
    CHECK(sources[3].first_row == sources[2].first_row + sources[2].rows);
    CHECK(info.rows == sources.back().first_row + sources.back().rows);

    const std::array reordered_patches{patches[2], patches[0], patches[1]};
    const std::array reordered_bindings{
        PgConstraintTraceBinding{bindings[2].id, bc, 1, core::Index{2}},
        PgConstraintTraceBinding{bindings[1].id, ab, 0, core::Index{1}},
        PgConstraintTraceBinding{bindings[0].id, fixed, 2, {}}};
    auto reordered = take(prepare(reordered_patches, reordered_bindings, limits, plan));
    const auto earlier_duplicate =
        std::array{PgConstraintTraceBinding{PgConstraintId{4}, ab, 2, core::Index{0}}, bindings[0],
                   bindings[1], bindings[2]};
    auto priority_changed = take(prepare(patches, earlier_duplicate, limits, plan));
    const std::vector<PgCartesianExactPosition> zero_support(c.reference_coordinates().size());
    const auto conflicting = take(PgFullTrace::prepare_support(
        c, std::span{cl}.first(1), ci, zero_support, trace_limits, workspace, plan));
    const std::array conflicting_bindings{
        bindings[0], PgConstraintTraceBinding{PgConstraintId{990}, conflicting, 1, {}}};
    auto incompatible = take(prepare(patches, conflicting_bindings, limits, plan));

    auto exact_count = limits;
    exact_count.maximum_rows = static_cast<std::uint64_t>(info.rows);
    exact_count.maximum_nonzeros = static_cast<std::uint64_t>(info.nonzeros);
    auto just_enough = take(prepare(patches, bindings, exact_count, plan));
    CHECK(take(just_enough.info()).nonzeros == info.nonzeros);
    const std::array<PgConstraintTraceBinding, 0> no_traces{};
    auto unconstrained = take(prepare(patches, no_traces, limits, plan));
    CHECK(take(unconstrained.info()).rows == 0);
    CHECK(take(unconstrained.info()).columns == info.columns);
    CHECK(unconstrained.original_rows().empty());

    auto short_rows = exact_count;
    --short_rows.maximum_rows;
    check_failure(prepare(patches, bindings, short_rows, plan),
                  core::FailureCode::ResourceExceeded);
    auto short_terms = exact_count;
    --short_terms.maximum_nonzeros;
    check_failure(prepare(patches, bindings, short_terms, plan),
                  core::FailureCode::ResourceExceeded);

    const auto column_count = static_cast<std::size_t>(info.columns);
    const auto values = states(f, column_count, 5);
    const auto lambda = states(f, static_cast<std::size_t>(info.rows), 3);
    const auto blocks = [&](std::size_t i) {
        return std::span{values}.subspan(static_cast<std::size_t>(layout[i].first),
                                         static_cast<std::size_t>(layout[i].count));
    };
    const std::array local_inputs{PgTraceCoordinates{ci, {}, blocks(1), {}},
                                  PgTraceCoordinates{ai, bi, blocks(2), blocks(0)},
                                  PgTraceCoordinates{bi, ci, blocks(0), blocks(1)},
                                  PgTraceCoordinates{ai, bi, blocks(2), blocks(0)}};
    const PgGlobalConstraintCoordinates state{layout, values};
    const auto forward = take(global.apply(state, limits.arithmetic, workspace, plan));
    const auto residual = take(global.residual(state, limits.arithmetic, workspace, plan));
    const auto force = take(global.transpose(lambda, limits.arithmetic, workspace, plan));
    std::vector<std::array<Q, 3>> expected_force(column_count);
    bool have_zero_row = false, have_nonzero_right_hand_side = false;
    for (std::size_t i = 0; i < bindings.size(); ++i) {
        const auto& local = bindings[i].trace;
        const auto first = static_cast<std::size_t>(sources[i].first_row);
        const auto count = static_cast<std::size_t>(sources[i].rows);
        const auto baseline =
            take(local.apply(local_inputs[i], limits.arithmetic, workspace, plan));
        const auto local_residual =
            take(local.residual(local_inputs[i], limits.arithmetic, workspace, plan));
        const auto local_force = take(local.transpose(std::span{lambda}.subspan(first, count),
                                                      limits.arithmetic, workspace, plan));
        const auto local_left = static_cast<std::size_t>(bindings[i].left_patch);
        const auto left_count = static_cast<std::size_t>(layout[local_left].count);
        const auto global_col = [&](core::Index column) {
            if (static_cast<std::size_t>(column) < left_count)
                return layout[local_left].first + column;
            return layout[static_cast<std::size_t>(bindings[i].right_patch.value())].first +
                   column - static_cast<core::Index>(left_count);
        };
        CHECK(sources[i].id == bindings[i].id);
        REQUIRE(sources[i].overlaps == static_cast<core::Index>(local.overlaps().size()));
        const auto overlap_start = static_cast<std::size_t>(sources[i].first_overlap);
        for (std::size_t j = 0; j < local.overlaps().size(); ++j) {
            const auto& reference = local.overlaps()[j];
            const auto& stored = global.overlaps()[overlap_start + j];
            CHECK(stored.left_edge == reference.left_edge);
            CHECK(stored.right_edge == reference.right_edge);
            CHECK(stored.first == reference.first);
            CHECK(stored.last == reference.last);
            CHECK(stored.left_depth.lower_first == reference.left_depth.lower_first);
            CHECK(stored.right_depth.upper_last == reference.right_depth.upper_last);
        }
        for (std::size_t j = 0; j < count; ++j) {
            const auto raw = take(global.row(static_cast<core::Index>(first + j)));
            const auto& original = local.original_rows()[j];
            const auto& gathered = global.original_rows()[first + j];
            CHECK(gathered.source == static_cast<core::Index>(i));
            CHECK(gathered.source_row == static_cast<core::Index>(j));
            CHECK(gathered.source_overlap == original.overlap);
            CHECK(gathered.field == original.field);
            CHECK(gathered.power == original.power);
            CHECK(gathered.local_representative == original.unique_row);
            CHECK(gathered.local_scale == original.scale);
            CHECK(gathered.prescribed == original.prescribed);
            CHECK(forward.values()[first + j] == baseline.values()[j]);
            CHECK(residual.values()[first + j] == local_residual.values()[j]);
            have_nonzero_right_hand_side |= std::ranges::any_of(
                gathered.prescribed, [](const auto& x) { return x.sign() != 0; });
            if (original.unique_row < 0) {
                CHECK(raw.empty());
                have_zero_row = true;
                continue;
            }
            const auto local_row = take(local.unique_row(original.unique_row));
            REQUIRE(raw.size() == local_row.size());
            CHECK(std::ranges::is_sorted(raw, std::ranges::less{}, &PgTraceSparseEntry::column));
            for (const auto& term : local_row) {
                const auto column = global_col(term.column);
                const auto found = std::ranges::lower_bound(raw, column, std::ranges::less{},
                                                            &PgTraceSparseEntry::column);
                REQUIRE(found != raw.end());
                CHECK(found->column == column);
                CHECK(q(found->coefficient) == q(original.scale) * q(term.coefficient));
            }
        }
        for (std::size_t j = 0; j < local_force.size(); ++j)
            for (std::size_t axis = 0; axis < 3; ++axis)
                expected_force[static_cast<std::size_t>(global_col(static_cast<core::Index>(j)))]
                              [axis] += q(local_force.values()[j][axis]);
    }
    CHECK(have_zero_row);
    CHECK(have_nonzero_right_hand_side);
    for (std::size_t j = 0; j < column_count; ++j)
        for (std::size_t axis = 0; axis < 3; ++axis)
            CHECK(q(force.values()[j][axis]) == expected_force[j][axis]);
    CHECK(work(forward.values(), lambda) == work(force.values(), values));

    const auto stale = std::array{layout[1], layout[0], layout[2]};
    check_failure(global.apply({stale, values}, limits.arithmetic, workspace, none),
                  core::FailureCode::StaleGeneration);
    check_failure(global.transpose(std::span{lambda}.first(1), limits.arithmetic, workspace, none),
                  core::FailureCode::InvalidInput);
    {
        const auto discard_ab = std::move(ab);
        const auto discard_bc = std::move(bc);
        const auto discard_support = std::move(fixed);
        const auto discard_a = std::move(a);
        const auto discard_b = std::move(b);
        const auto discard_c = std::move(c);
        CHECK(discard_ab.info().has_value());
        CHECK(discard_bc.info().has_value());
        CHECK(discard_support.info().has_value());
        CHECK(discard_a.info().has_value());
        CHECK(discard_b.info().has_value());
        CHECK(discard_c.info().has_value());
    }
    const auto independent = take(global.apply(state, limits.arithmetic, workspace, plan));
    CHECK(std::ranges::equal(independent.values(), forward.values()));

    auto chart_limits = PgSparseChartLimits{limits.arithmetic, 300, 400, 30'000, true};
    chart_limits.arithmetic.algebra.maximum_value_bits = 65536;
    chart_limits.arithmetic.algebra.maximum_work_bits = 262144;
    check_failure(PgSparseChart::prepare(std::move(incompatible), chart_limits, workspace, plan),
                  core::FailureCode::InvalidInput);
    REQUIRE(incompatible.info().has_value());
    auto no_pivots = chart_limits;
    no_pivots.maximum_pivots = 0;
    const auto no_factors = PgSparseChartLimits{chart_limits.arithmetic, 0, 0, 0};
    const auto empty_chart =
        take(PgSparseChart::prepare(std::move(unconstrained), no_factors, workspace, plan));
    CHECK(empty_chart.pivots().empty());
    CHECK(empty_chart.free_columns().size() == column_count);
    CHECK_FALSE(take(empty_chart.preparation()).timing.has_value());
    check_failure(PgSparseChart::prepare(std::move(global), no_pivots, workspace, plan),
                  core::FailureCode::ResourceExceeded);
    REQUIRE(global.info().has_value());
    auto chart = take(PgSparseChart::prepare(std::move(global), chart_limits, workspace, plan));
    auto permuted =
        take(PgSparseChart::prepare(std::move(reordered), chart_limits, workspace, plan));
    auto preferred =
        take(PgSparseChart::prepare(std::move(priority_changed), chart_limits, workspace, plan));
    const auto& original = chart.original();
    REQUIRE_FALSE(chart.pivots().empty());
    const auto report = take(chart.preparation());
    CHECK(report.rank == chart.pivots().size());
    CHECK(report.original_nonzeros == static_cast<std::uint64_t>(info.nonzeros));
    CHECK(report.lower_terms > 0);
    CHECK(report.upper_terms >= report.rank);
    CHECK(report.lower_terms + report.upper_terms <= report.factor_term_capacity);
    CHECK(report.pivot_capacity <= chart_limits.maximum_pivots);
    CHECK(report.work_row_capacity <= chart_limits.maximum_row_terms);
    CHECK(report.arithmetic.allocated_owned_payload_bytes > 0);
    CHECK(report.timing.value_or(core::StageTiming{}).calls == 1);
    CHECK(report.timing.value_or(core::StageTiming{}).elapsed_nanoseconds > 0);
    CHECK(
        chart.original()
            .sources()[static_cast<std::size_t>(
                original
                    .original_rows()[static_cast<std::size_t>(chart.pivots().front().original_row)]
                    .source)]
            .kind == PgTraceKind::Support);
    REQUIRE(chart.free_columns().size() + chart.pivots().size() == column_count);
    REQUIRE(chart.pivots().size() == permuted.pivots().size());
    REQUIRE(chart.pivots().size() == preferred.pivots().size());
    const auto physical = [](const PgGlobalConstraints& equations, core::Index column) {
        const auto physical_layout = equations.layout();
        const auto block = std::ranges::find_if(physical_layout, [column](const auto& candidate) {
            return column >= candidate.first && column < candidate.first + candidate.count;
        });
        REQUIRE(block != physical_layout.end());
        return std::pair{block->id.value(), column - block->first};
    };
    const auto raw_key = [](const PgGlobalConstraints& equations, core::Index row) {
        const auto& metadata = equations.original_rows()[static_cast<std::size_t>(row)];
        return std::pair{equations.sources()[static_cast<std::size_t>(metadata.source)].id.value(),
                         metadata.source_row};
    };
    bool changed_priority = false;
    for (std::size_t i = 0; i < chart.pivots().size(); ++i) {
        const auto x = chart.pivots()[i], y = permuted.pivots()[i], z = preferred.pivots()[i];
        const auto key = raw_key(original, x.original_row);
        CHECK(key == raw_key(permuted.original(), y.original_row));
        CHECK(physical(original, x.dependent_column) ==
              physical(permuted.original(), y.dependent_column));
        CHECK(physical(original, x.dependent_column) ==
              physical(preferred.original(), z.dependent_column));
        const auto changed = raw_key(preferred.original(), z.original_row);
        if (key.first == 118) {
            CHECK(changed == (std::pair{std::uint64_t{4}, key.second}));
            changed_priority = true;
        } else {
            CHECK(changed == key);
        }
        CHECK(key.first != 660);
    }
    CHECK(changed_priority);
    REQUIRE(chart.free_columns().size() == permuted.free_columns().size());
    for (std::size_t i = 0; i < chart.free_columns().size(); ++i)
        CHECK(physical(original, chart.free_columns()[i]) ==
              physical(permuted.original(), permuted.free_columns()[i]));

    const auto& limits_for_action = chart_limits.arithmetic;
    std::vector<PgCartesianExactPosition> independent_state(chart.free_columns().size());
    for (std::size_t i = 0; i < independent_state.size(); ++i)
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const auto column = static_cast<std::size_t>(chart.free_columns()[i]);
            const auto& exact = q(values[column][axis]);
            independent_state[i][axis] = rat(f, numerator(exact), denominator(exact));
        }
    const auto identity = take(chart.identity());
    const auto other_identity = take(permuted.identity());
    CHECK(identity != other_identity);
    const PgSparseChartInput rb{identity, independent_state};
    check_failure(
        chart.lift({other_identity, independent_state}, limits_for_action, workspace, none),
        core::FailureCode::StaleGeneration);
    {
        const auto lease = take(workspace.acquire());
        static_cast<void>(lease);
        check_failure(chart.lift(rb, limits_for_action, workspace, none),
                      core::FailureCode::InvalidInput);
    }
    const auto lifted = take(chart.lift(rb, limits_for_action, workspace, plan));
    const auto dispatched = take(original.residual({original.layout(), lifted.values()},
                                                   limits_for_action, workspace, plan));
    for (const auto& row : dispatched.values())
        CHECK(std::ranges::all_of(row, [](const auto& value) { return value.sign() == 0; }));

    const auto tangent = take(chart.tangent(rb, limits_for_action, workspace, plan));
    const auto full_force = states(f, column_count, 9);
    const auto dual =
        take(chart.transpose({identity, full_force}, limits_for_action, workspace, plan));
    CHECK(work(full_force, tangent.values()) ==
          work(dual.independent_force.values(), independent_state));
    const auto reaction = take(
        original.transpose(dual.original_multipliers.values(), limits_for_action, workspace, plan));
    std::vector<std::array<Q, 3>> recomposed(column_count);
    for (std::size_t i = 0; i < column_count; ++i)
        for (std::size_t axis = 0; axis < 3; ++axis)
            recomposed[i][axis] = q(reaction.values()[i][axis]);
    for (std::size_t i = 0; i < chart.free_columns().size(); ++i)
        for (std::size_t axis = 0; axis < 3; ++axis)
            recomposed[static_cast<std::size_t>(chart.free_columns()[i])][axis] +=
                q(dual.independent_force.values()[i][axis]);
    for (std::size_t i = 0; i < column_count; ++i)
        for (std::size_t axis = 0; axis < 3; ++axis)
            CHECK(recomposed[i][axis] == q(full_force[i][axis]));
    const auto combined = static_cast<std::size_t>(original.sources().back().first_row);
    for (const auto& multiplier : dual.original_multipliers.values().subspan(combined))
        CHECK(std::ranges::all_of(multiplier, [](const auto& value) { return value.sign() == 0; }));

    // For an exact time coefficient, the raw value must satisfy every equation.
    // C*values is a compatible coefficient. A ramp on only the repeated AB row
    // agrees with zero at t=0, but must be rejected for its nonzero coefficient.
    const auto compatible =
        take(chart.lift_rhs(rb, {identity, forward.values()}, limits_for_action, workspace, plan));
    CHECK(std::ranges::equal(compatible.values(), values));
    std::vector<PgCartesianExactPosition> inconsistent(static_cast<std::size_t>(info.rows));
    const auto original_ab = original.sources().back();
    bool have_duplicate = false;
    for (core::Index j = original_ab.first_row; j < original_ab.first_row + original_ab.rows; ++j) {
        if (!take(original.row(j)).empty()) {
            inconsistent[static_cast<std::size_t>(j)][0] = f.r("1");
            have_duplicate = true;
            break;
        }
    }
    REQUIRE(have_duplicate);
    check_failure(chart.lift_rhs(rb, {identity, inconsistent}, limits_for_action, workspace, plan),
                  core::FailureCode::InvalidInput);
}
