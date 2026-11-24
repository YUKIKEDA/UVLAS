#include "pg_conditioning_test_support.hpp"

#include <algorithm>
#include <string>
#include <uvlas/structure/pg_flat_gram.hpp>
#include <uvlas/structure/pg_full_trace.hpp>

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
    const auto support =
        take(PgFullTrace::prepare_support(a, only_first, info.left, ca, limits, ws, plan));
    const auto fixed = take(support.info());
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
