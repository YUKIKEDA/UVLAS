#include "pg_conditioning_test_support.hpp"

#include <algorithm>
#include <string>
#include <uvlas/structure/pg_discrete_gram.hpp>
#include <uvlas/structure/pg_full_trace.hpp>
#include <uvlas/structure/pg_prepared_patch.hpp>
#include <uvlas/structure/pg_prepared_volume.hpp>

using namespace uvlas;
using namespace uvlas::structure;
using namespace partition_test;
using namespace pg_conditioning_test;

namespace {

constexpr int curved_rule_order = 5;
constexpr int curved_rule_planar_count = (curved_rule_order + 1) * (curved_rule_order + 2) / 2;
constexpr int curved_rule_depth_count = 3;
constexpr auto curved_rule_points_per_volume = static_cast<std::size_t>(curved_rule_planar_count) *
                                               static_cast<std::size_t>(curved_rule_depth_count);
constexpr int curved_rule_weight_denominator =
    2 * curved_rule_planar_count * curved_rule_depth_count;

PgReferenceBasis parabolic_basis(Fixture& f, core::Workspace& scratch, core::ResourcePlan& plan,
                                 bool adjacent = false) {
    const auto zero = f.r("0"), one = f.r("1"), length = f.r("2");
    const auto coefficient = f.r("1", "16");
    const std::array<ExactInputTerm, 1> terms{{{{2, 0}, coefficient}}};
    const auto curve = take(ExactChartPolynomial::make(terms, exact_limits(), f.input_plan));
    f.graph = take(PgGraph::prepare({zero, one, zero, curve, 1, PgGraphConvention::Version1},
                                    {exact_limits(), 64, {}}, f.input_plan));
    std::array<ExactRational, 9> rotation;
    for (const auto i : {0U, 4U, 8U})
        rotation[i] = f.r("1");
    const std::array origin{ExactRational{}, adjacent ? f.r("2") : ExactRational{},
                            ExactRational{}};
    return make_basis(f, {rotation, origin, length}, scratch, plan);
}

struct CurvedRulePoint {
    core::Index volume;
    PgCellIndex cell;
    ExactRational u, v, depth;
    Q s, t, zeta, root, measure;
};

// Sample the parabola z=s²/16. Rational Pythagorean slopes let this
// finite-rule oracle use only exact rational arithmetic.
std::vector<CurvedRulePoint> pythagorean_rule(Fixture& f, const PgIntegrationPartition& partition) {
    const auto volumes = take(partition.volumes());
    const auto triangles = take(partition.triangles());
    const auto vertices = take(partition.vertices());
    const auto pieces = take(partition.pieces());
    const auto layers = take(partition.layers());
    const std::array<Q, curved_rule_depth_count> depths{Q("1/4"), Q("1/2"), Q("3/4")};
    constexpr int order = curved_rule_order;
    const Q half_thickness("1/4"), curvature("1/8");
    const Q weight = Q{1} / curved_rule_weight_denominator;
    std::vector<CurvedRulePoint> rule;
    for (std::size_t volume = 0; volume < volumes.size(); ++volume) {
        const auto& record = volumes[volume];
        const auto& tri = triangles[static_cast<std::size_t>(record.triangle)];
        const auto& layer = layers[static_cast<std::size_t>(record.layer)];
        const auto cell = pieces[static_cast<std::size_t>(tri.piece)].cell;
        const auto& a = vertices[static_cast<std::size_t>(tri.vertices[0])];
        const auto& b = vertices[static_cast<std::size_t>(tri.vertices[1])];
        const auto& c = vertices[static_cast<std::size_t>(tri.vertices[2])];
        const Q es = q(b.s) - q(a.s), et = q(b.t) - q(a.t);
        const Q fs = q(c.s) - q(a.s), ft = q(c.t) - q(a.t);
        const Q determinant = es * ft - et * fs;
        REQUIRE(determinant > 0);
        const Q width = q(layer.upper) - q(layer.lower);
        for (int i = 0; i <= order; ++i)
            for (int j = 0; j <= order - i; ++j) {
                const Q u0 = Q{128 * (i + 1) + (j + 1) * (j + 2)} / (128 * (order + 3));
                const Q v0 = Q{128 * (j + 1) + (i + 1) * (i + 3)} / (128 * (order + 3));
                const Q s0 = q(a.s) + u0 * es + v0 * fs;
                const Q t = q(a.t) + u0 * et + v0 * ft;
                const Q k = curvature * s0 / 2, k2 = k * k;
                const Q s = s0 / (1 - k2), root = (1 + k2) / (1 - k2);
                REQUIRE(root * root == 1 + curvature * curvature * s * s);
                const Q u = ((s - q(a.s)) * ft - (t - q(a.t)) * fs) / determinant;
                const Q v = (es * (t - q(a.t)) - et * (s - q(a.s))) / determinant;
                REQUIRE(u > 0);
                REQUIRE(v > 0);
                REQUIRE(u + v < 1);
                for (const auto& depth : depths) {
                    const Q zeta = q(layer.lower) + width * depth;
                    const Q offset = zeta * half_thickness;
                    // dX/ds scales by (1-offset*kappa) in the planar
                    // tangent with curvature f''/(1+f'²)^(3/2).
                    const Q measure = weight * determinant * width * half_thickness * root *
                                      (1 - offset * curvature / (root * root * root));
                    REQUIRE(measure > 0);
                    rule.push_back(
                        {static_cast<core::Index>(volume), cell,
                         rat(f, numerator(u), denominator(u)), rat(f, numerator(v), denominator(v)),
                         rat(f, numerator(depth), denominator(depth)), s, t, zeta, root, measure});
                }
            }
    }
    return rule;
}

Q independent_curve_value(const CurvedRulePoint& p, std::span<const Term> terms) {
    const Q u = 4 * p.s - (2 * p.cell.first + 1);
    const Q v = 4 * p.t - (2 * p.cell.second + 1);
    const Q z = p.zeta / 4;
    const std::array<Q, 3> through{Q{1}, z / (p.root * p.root * p.root), z * z / 2};
    Q value;
    for (const auto& term : terms) {
        if (term.row.cell == p.cell) {
            value += term.value * power(u, term.row.powers.first) *
                     power(v, term.row.powers.second) *
                     through[static_cast<std::size_t>(term.row.field)];
        }
    }
    return value;
}

const PgIntegrationPartition& trace_partition(const PgPreparedPatch& patch) {
    const auto* basis = take(patch.conditioning().reference_basis());
    return take(basis->source()).partition;
}

std::vector<PgTraceEdgeSpan> curved_edge(const PgPreparedPatch& patch, const ExactRational& fixed_t,
                                         const ExactRational& zero, const ExactRational& one) {
    const auto& partition = trace_partition(patch);
    const auto vertices = take(partition.vertices());
    const auto edges = take(partition.edges());
    std::vector<PgTraceEdgeSpan> spans;
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const auto& edge = edges[i];
        const auto& a = vertices[static_cast<std::size_t>(edge.vertices[0])];
        const auto& b = vertices[static_cast<std::size_t>(edge.vertices[1])];
        if (!edge.external || q(a.t) != q(fixed_t) || q(b.t) != q(fixed_t))
            continue;
        const bool forward = q(a.s) < q(b.s);
        spans.push_back({static_cast<core::Index>(i),
                         forward ? a.s : b.s,
                         forward ? b.s : a.s,
                         forward ? zero : one,
                         forward ? one : zero,
                         {}});
    }
    // Sort indices because the borrowed record itself is immutable.
    std::vector<const PgTraceEdgeSpan*> ordered;
    ordered.reserve(spans.size());
    for (const auto& span : spans)
        ordered.push_back(&span);
    std::ranges::sort(ordered,
                      [](const auto* a, const auto* b) { return q(a->first) < q(b->first); });
    std::vector<PgTraceEdgeSpan> sorted;
    sorted.reserve(ordered.size());
    for (const auto* span : ordered)
        sorted.push_back(*span);
    return sorted;
}

Q curved_point_value(Fixture& f, const PgPreparedPatch& patch, const PgTraceEdgeSpan& edge,
                     const Q& s, const Q& root, std::span<const PgCartesianExactPosition> state,
                     core::ResourcePlan& plan) {
    const auto faces = take(trace_partition(patch).lateral_faces());
    const auto selected = std::ranges::find_if(
        faces, [&](const auto& face) { return face.edge == edge.edge && face.layer == 1; });
    REQUIRE(selected != faces.end());
    const Q fraction = (s - q(edge.first)) / (q(edge.last) - q(edge.first));
    const bool forward = q(edge.along_first) < q(edge.along_last);
    const Q along = forward == selected->canonical_direction ? fraction : 1 - fraction;
    const auto a = rat(f, numerator(along), denominator(along));
    const auto depth = f.r("1", "2"); // zeta = 1/3, physical z = 1/6.
    const PgReferencePoint point =
        PgReferenceLateralPoint{static_cast<core::Index>(selected - faces.begin()), a, depth};
    const auto value = take(patch.evaluate_reference(std::span{&point, 1}, query_limits(), plan));
    CHECK(q(value.points.values().front().position_metres[1].rational) == 2);
    CHECK(q(value.points.values().front().g) == root * root);
    Q result;
    for (std::size_t i = 0; i < state.size(); ++i) {
        const auto& n = value.columns.values()[i].value;
        result += q(state[i][0]) * (q(n.rational) + root * q(n.radical));
    }
    return result;
}

void check_curved_trace(Fixture& f, const PgPreparedPatch& left, const PgPreparedPatch& right,
                        const PgReferenceLimits& limits, core::Workspace& workspace,
                        core::ResourcePlan& plan) {
    const auto zero = f.r("0"), one = f.r("1");
    const auto a = curved_edge(left, one, zero, one), b = curved_edge(right, zero, zero, one);
    REQUIRE(a.size() == 2);
    REQUIRE(b.size() == 2);
    const auto trace =
        take(PgFullTrace::prepare(left, a, right, b, {limits, 100'000}, workspace, plan));
    const auto info = take(trace.info());
    REQUIRE(info.overlaps == 2);
    std::vector<PgCartesianExactPosition> q_left(left.reference_coordinates().size()),
        q_right(right.reference_coordinates().size());
    for (std::size_t j = 0; j < q_left.size(); ++j)
        q_left[j][0] = f.r(std::to_string(static_cast<int>(j % 7) - 3), "17");
    for (std::size_t j = 0; j < q_right.size(); ++j)
        q_right[j][0] = f.r(std::to_string(static_cast<int>(j % 5) - 2), "11");
    const auto residual =
        take(trace.apply({info.left, info.right, q_left, q_right}, limits, workspace, plan));
    // f(s)=s²/16. The rational Pythagorean slope makes sqrt(1+(s/8)²) exact.
    const Q s("1024/4095"), root("4097/4095"), z("1/6");
    const auto& segment = trace.overlaps().front();
    const Q local = (s - q(segment.first)) / (q(segment.last) - q(segment.first));
    const Q normal_x = -s / 8;
    std::array<Q, 3> field{};
    for (std::size_t i = 0; i < trace.original_rows().size(); ++i) {
        const auto& row = trace.original_rows()[i];
        if (row.overlap == 0)
            field[static_cast<std::size_t>(row.field)] +=
                power(local, row.power) * q(residual.values()[i][0]);
    }
    const Q restored =
        field[0] + z * field[1] / (power(normal_x, 3) * power(root, 3)) + z * z * field[2];
    REQUIRE(field[1] != 0);
    CHECK(restored == curved_point_value(f, left, a.front(), s, root, q_left, plan) -
                          curved_point_value(f, right, b.front(), s, root, q_right, plan));
}

void check_curved_area_projection(Fixture& f, const PgPreparedPatch& patch,
                                  const PgReferenceLimits& limits, core::ResourcePlan& plan) {
    const auto& partition = trace_partition(patch);
    const auto faces = take(partition.horizontal_faces());
    const auto layers = take(partition.layers());
    const auto selected = std::ranges::find_if(faces, [&](const auto& face) {
        const auto& layer = layers[static_cast<std::size_t>(face.layer)];
        return q(face.side == PgPartitionFaceSide::Lower ? layer.lower : layer.upper) != 0;
    });
    REQUIRE(selected != faces.end());
    const auto& face = *selected;
    const auto& tri = take(partition.triangles())[static_cast<std::size_t>(face.triangle)];
    const auto& layer = layers[static_cast<std::size_t>(face.layer)];
    const auto vertices = take(partition.vertices());
    Q s;
    for (const auto index : tri.vertices)
        s += q(vertices[static_cast<std::size_t>(index)].s) / 3;
    const Q z = q(face.side == PgPartitionFaceSide::Lower ? layer.lower : layer.upper) / 2;
    REQUIRE(z != 0);
    const Q g = 1 + s * s / 64;
    const I grid = I{1} << 256;
    const Dec approximate = sqrt(numerator(g).convert_to<Dec>() / denominator(g).convert_to<Dec>());
    const auto floor_root = floor(approximate * grid.convert_to<Dec>()).convert_to<I>();
    const Q lower_root = Q{floor_root} / grid, upper_root = Q{floor_root + 1} / grid;
    REQUIRE(lower_root * lower_root < g);
    REQUIRE(upper_root * upper_root > g);
    // For X(s,t,z) = 2(s,t,s²/16) + z n, the tangent lengths on a
    // fixed-z face are 2 sqrt(g) - z/(8g) and 2, in orthogonal directions.
    const auto area = [&](const Q& root) -> Q {
        return q(tri.twice_area) * (4 * root - z / (4 * g));
    };
    const auto third = f.r("1", "3");
    const PgReferencePoint query = PgReferenceHorizontalPoint{
        static_cast<core::Index>(selected - faces.begin()), third, third};
    constexpr double tolerance = 0x1p-72;
    const PgPreparedReferenceProjection projection{
        core::Precision::TwoWord, tolerance, tolerance, tolerance, tolerance, tolerance};
    const auto projected =
        take(patch.project_reference(std::span{&query, 1}, projection, limits, plan));
    const Q measure = independent_double(projected.high_words()[3]) +
                      independent_double(projected.low_words()[3]);
    const Q budget = independent_double(tolerance);
    CHECK(measure - budget <= area(lower_root));
    CHECK(measure + budget >= area(upper_root));
}

void check_fixed_projection(const PgPreparedPatch& patch, std::span<const CurvedRulePoint> rule,
                            const PgReferenceLimits& limits, core::ResourcePlan& plan) {
    REQUIRE(rule.size() > 2);
    const auto source = take(patch.info());
    for (const auto precision : {core::Precision::Binary64, core::Precision::TwoWord}) {
        const bool paired = precision == core::Precision::TwoWord;
        const PgVolumeMassProjection projection{precision, paired ? 0x1p-82 : 0x1p-42,
                                                paired ? 0x1p-88 : 0x1p-45};
        core::ResourcePlan none{0};
        check_failure(
            PgPreparedVolumeReference::prepare(patch, rule.size(), 1, projection, limits, none),
            core::FailureCode::InvalidInput);
        check_failure(PgPreparedVolumeReference::prepare(patch, 1, 2, projection, limits, none),
                      core::FailureCode::ResourceExceeded);
        auto projected =
            take(PgPreparedVolumeReference::prepare(patch, 1, 2, projection, limits, plan));
        CHECK(projected.source() == source.identity);
        const auto rank = projected.column_count();
        REQUIRE(rank == source.basis.rank);
        const auto retained = projected.words();
        REQUIRE(retained.high.size() == 2 * (rank + 1));
        REQUIRE(retained.low.size() == (paired ? retained.high.size() : 0));
        const auto rational_words = [](double high, double low) -> Q {
            return independent_double(high) + independent_double(low);
        };
        const auto lo = [&](std::span<const double> low, std::size_t i) {
            return paired ? low[i] : 0.0;
        };
        std::vector<Q> exact_shapes(2 * rank);
        for (std::size_t row = 0; row < projected.points().size(); ++row) {
            const auto& actual = projected.points()[row];
            REQUIRE(actual.rule_ordinal == row + 1);
            CHECK(q(actual.weight) == Q{1} / curved_rule_weight_denominator);
            const auto& source_point = rule[row + 1];
            const PgReferencePoint query = PgReferenceVolumePoint{
                source_point.volume, source_point.u, source_point.v, source_point.depth};
            const auto exact = take(patch.evaluate_reference(std::span{&query, 1}, limits, plan));
            const auto& point = exact.points.values().front();
            CHECK(actual.reference.identity == point.identity);
            CHECK(actual.reference.partition_record == point.partition_record);
            CHECK(actual.reference.layer == point.layer);
            const auto slot0 = row * (rank + 1);
            const auto& jacobian = std::get<PgReferenceVolumeMeasure>(point.measure).jacobian;
            const auto exact_j = q(jacobian.rational) + source_point.root * q(jacobian.radical);
            CHECK(abs(rational_words(retained.high[slot0], lo(retained.low, slot0)) - exact_j) <=
                  independent_double(projection.maximum_volume_jacobian_error_metres_cubed));
            for (std::size_t j = 0; j < rank; ++j) {
                const auto& n = exact.columns.values()[j].value;
                auto& scalar = exact_shapes[row * rank + j];
                scalar = q(n.rational) + source_point.root * q(n.radical);
                const auto i = slot0 + j + 1;
                CHECK(abs(rational_words(retained.high[i], lo(retained.low, i)) - scalar) <=
                      independent_double(projection.maximum_shape_error));
            }
        }
        std::vector<double> coords(3 * rank), inputs_low(paired ? coords.size() : 0);
        for (std::size_t j = 0; j < coords.size(); ++j)
            coords[j] = static_cast<double>(static_cast<int>(j % 9) - 4) / 16;
        const std::array<double, 6> lambda_high{0.5, -0.25, 0.75, -0.125, 0.375, -0.5};
        const std::array<double, 6> lambda_low{0x1p-60, 0, 0, 0, 0, 0};
        const PgCartesianWords q_words{precision, coords, inputs_low};
        const PgCartesianWords lambda{precision, lambda_high,
                                      paired ? std::span<const double>{lambda_low}
                                             : std::span<const double>{}};
        auto output = take(PgVolumePositionResult::allocate(static_cast<core::Index>(coords.size()),
                                                            precision, plan));
        const auto workspace = take(core::Workspace::create(output.workspace_bytes(), plan));
        const double bound = paired ? 0x1p-68 : 0x1p-29;
        const PgVolumePositionLimits application{100'000, bound, {}};
        const auto map = PgVolumePositionDirection::ToPoints;
        const auto adjoint = PgVolumePositionDirection::TransposeToColumns;
        check_failure(apply_pg_prepared_volume_position(projected, {}, map, q_words, application,
                                                        *workspace, output),
                      core::FailureCode::StaleGeneration);
        REQUIRE(apply_pg_prepared_volume_position(projected, source.identity, map, q_words,
                                                  application, *workspace, output));
        auto action = output.words();
        REQUIRE(action.high.size() == lambda_high.size());
        for (std::size_t row = 0; row < 2; ++row)
            for (std::size_t axis = 0; axis < 3; ++axis) {
                Q expected;
                for (std::size_t j = 0; j < rank; ++j)
                    expected +=
                        exact_shapes[row * rank + j] * independent_double(coords[3 * j + axis]);
                const auto i = 3 * row + axis;
                CHECK(abs(rational_words(action.high[i], lo(action.low, i)) - expected) <=
                      independent_double(bound));
            }
        REQUIRE(apply_pg_prepared_volume_position(projected, source.identity, adjoint, lambda,
                                                  application, *workspace, output));
        action = output.words();
        REQUIRE(action.high.size() == coords.size());
        for (std::size_t j = 0; j < rank; ++j)
            for (std::size_t axis = 0; axis < 3; ++axis) {
                Q expected;
                for (std::size_t row = 0; row < 2; ++row) {
                    const auto i = 3 * row + axis;
                    expected += exact_shapes[row * rank + j] *
                                rational_words(lambda.high[i], lo(lambda.low, i));
                }
                const auto i = 3 * j + axis;
                CHECK(abs(rational_words(action.high[i], lo(action.low, i)) - expected) <=
                      independent_double(bound));
            }
    }
}

void check_projected_weight_admission(const PgPreparedPatch& patch, const PgReferenceLimits& limits,
                                      core::ResourcePlan& plan) {
    const auto source = take(patch.info());
    REQUIRE(source.basis.rank == 130);
    constexpr std::uint64_t limb_bits = detail::IntegerBackend::limb_bits;
    const auto limbs = limits.algebra.maximum_work_bits / limb_bits;
    const auto integer_payload = 4 * limbs * sizeof(boost::multiprecision::limb_type);
    // An isolated exact point has 336 work slots, 74 point slots and 8 per basis column.
    const auto point_rationals = 336 + 74 + 8 * source.basis.rank;
    const auto point_payload = (2 * point_rationals + 1) * integer_payload;
    const auto isolated = take(
        detail::exact_rational_admission(limits.algebra, (limbs - 2) * limb_bits, point_payload));

    const PgVolumeMassProjection projection{core::Precision::Binary64, 0x1p-42, 0x1p-45};
    auto restricted = limits;
    auto& budget = std::get<ExactEstimateAdmission>(restricted.algebra.external_memory).bytes;
    budget = isolated.estimated_external_payload_bytes;
    REQUIRE(PgPreparedVolumeReference::prepare(patch, 1, 1, projection, restricted, plan));
    REQUIRE(PgPreparedVolumeReference::prepare(patch, 2, 1, projection, restricted, plan));
    check_failure(PgPreparedVolumeReference::prepare(patch, 1, 2, projection, restricted, plan),
                  core::FailureCode::ResourceExceeded);
    CHECK(take(patch.info()).identity == source.identity);

    budget += 2 * integer_payload;
    const auto pair =
        take(PgPreparedVolumeReference::prepare(patch, 1, 2, projection, restricted, plan));
    CHECK(pair.points().size() == 2);
}

} // namespace

TEST_CASE("PG curved discrete Gram encloses the selected rule and certifies its transform",
          "[SH-D03][SH-D07][N4-02][N4-03][N4-08][verification]") {
    Fixture f;
    core::Workspace scratch{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan setup{profile};
    auto basis = parabolic_basis(f, scratch, setup);
    const auto columns = static_cast<std::size_t>(take(basis.size()));
    CAPTURE(columns);
    const auto source = take(basis.source());
    Fixture adjacent;
    auto adjacent_basis = parabolic_basis(adjacent, scratch, setup, true);
    const auto* selected_a = take(basis.rank());
    const auto* selected_b = take(adjacent_basis.rank());
    CHECK(std::ranges::equal(take(selected_a->selected()), take(selected_b->selected())));
    auto rule = pythagorean_rule(f, source.partition);
    auto raw = basis_terms(*take(basis.rank()), columns, setup);
    auto weight = f.r("1", std::to_string(curved_rule_weight_denominator));
    std::vector<PgDiscreteVolumePoint> points;
    points.reserve(rule.size());
    for (const auto& p : rule)
        points.push_back({{p.volume, p.u, p.v, p.depth}, weight});
    const PgDiscreteGramGrid grid{352, 272};
    PgDiscreteGramLimits caps{query_limits(), steps};
    caps.arithmetic.maximum_scalar_operations = 500'000'000;
    caps.arithmetic.maximum_partition_visits = 250'000'000;
    const auto prepare = [&](std::span<const PgDiscreteVolumePoint> queries,
                             PgDiscreteGramGrid grids, const PgDiscreteGramLimits& limits,
                             core::ResourcePlan& plan) {
        return PgDiscreteGram::prepare(std::move(basis), queries, grids, limits, scratch, plan);
    };
    core::ResourcePlan empty{0};
    auto stop = caps;
    stop.arithmetic.cancellation = requested_stop_token();
    check_failure(prepare(points, grid, stop, empty), core::FailureCode::Cancelled);
    check_failure(prepare(points, grid, caps, empty), core::FailureCode::ResourceExceeded);
    check_failure(prepare(points, {271, 272}, caps, empty), core::FailureCode::InvalidInput);
    check_failure(
        prepare(std::span{points}.subspan(curved_rule_points_per_volume), grid, caps, empty),
        core::FailureCode::InvalidInput);
    // Keep the full volume coverage and change only the weight or order.
    const std::array invalid_weights{f.r("0"), f.r("-1")};
    for (const auto& invalid_weight : invalid_weights) {
        std::vector<PgDiscreteVolumePoint> invalid;
        invalid.reserve(points.size());
        invalid.push_back({points.front().point, invalid_weight});
        for (const auto& point : std::span{points}.subspan(1))
            invalid.push_back(point);
        check_failure(prepare(invalid, grid, caps, empty), core::FailureCode::InvalidInput);
        REQUIRE(basis.size());
    }
    std::vector<PgDiscreteVolumePoint> reversed;
    reversed.reserve(points.size());
    constexpr auto block = curved_rule_points_per_volume;
    for (std::size_t k = 0; k < points.size(); ++k)
        reversed.push_back(points[k == block + 1 ? 0 : k]);
    check_failure(prepare(reversed, grid, caps, empty), core::FailureCode::InvalidInput);
    REQUIRE(basis.size());
    CHECK(empty.allocated_bytes() == 0);
    REQUIRE(basis.size());

    core::ResourcePlan bounded{profile};
    auto no_steps = caps;
    no_steps.maximum_integer_steps = 0;
    check_failure(prepare(points, grid, no_steps, bounded), core::FailureCode::ResourceExceeded);
    REQUIRE(basis.size());
    core::ResourcePlan plan{profile};
    auto enclosure_owner = take(prepare(points, grid, caps, plan));
    CHECK_FALSE(basis.size());
    const auto info = enclosure_owner.info();
    CHECK(info.points == rule.size());
    CHECK(info.integer_steps > 0);
    CHECK(info.arithmetic.allocated_owned_payload_bytes == plan.allocated_bytes());
    const auto enclosure = take(enclosure_owner.enclosure());
    REQUIRE(enclosure.fractional_bits == grid.output_fractional_bits);
    REQUIRE(enclosure.upper_triangle.size() == slot(0, columns));
    const Q unit = Q{1} / (I{1} << grid.output_fractional_bits);
    const auto verify_entry = [&](std::size_t i, std::size_t j) {
        CAPTURE(i, j);
        Q expected;
        for (const auto& p : rule)
            expected +=
                p.measure * independent_curve_value(p, raw[i]) * independent_curve_value(p, raw[j]);
        const auto& entry = enclosure.upper_triangle[slot(i, j)];
        CHECK(q(entry.lower) <= expected);
        CHECK(q(entry.upper) >= expected);
        CHECK(q(entry.upper) - q(entry.lower) <= 2 * unit);
    };
    for (const auto j : std::array<std::size_t, 5>{0, 2, columns / 3, columns / 2, columns - 1}) {
        verify_entry(0, j);
        verify_entry(j, j);
        if (j > 1)
            verify_entry(j - 1, j);
    }
    std::vector<Q> centres;
    for (const auto& entry : enclosure.upper_triangle)
        centres.push_back((q(entry.lower) + q(entry.upper)) / 2);
    using ReferenceDecimal =
        boost::multiprecision::number<boost::multiprecision::cpp_dec_float<160>>;
    const auto independent_transform = proposal<ReferenceDecimal>(f, centres, columns, 216);
    // This is the same chart, thickness, local basis and finite rule under
    // a rigid origin translation, which leaves N and det(J) unchanged.
    auto adjacent_certificate = take(PgConditioningCertificate::generate_and_certify(
        std::move(adjacent_basis), enclosure, {512, 216}, {profile, steps, {}}, scratch, plan));
    const auto generate = [&](core::ResourcePlan& budget) {
        return std::move(enclosure_owner)
            .generate_and_certify({512, 216}, {profile, steps, {}}, scratch, budget);
    };
    core::ResourcePlan refused{0};
    check_failure(generate(refused), core::FailureCode::ResourceExceeded);
    REQUIRE(enclosure_owner.enclosure());
    const Q first = q(rule.front().u);
    const Q original_weight = q(weight);
    points.clear();
    rule.front().u = f.r("0");
    weight = f.r("2");
    REQUIRE(enclosure_owner.rule().size() == rule.size());
    CHECK(q(enclosure_owner.rule().front().u) == first);
    auto generated = generate(plan);
    if (!generated) {
        CAPTURE(static_cast<unsigned>(generated.error().code));
        REQUIRE(generated);
    }
    auto certified = take(std::move(generated));
    CHECK_FALSE(enclosure_owner.enclosure());
    CHECK(enclosure_owner.rule().empty());
    const auto frozen = certified.rule();
    REQUIRE(frozen.size() == rule.size());
    for (std::size_t i = 0; i < frozen.size(); ++i) {
        CHECK(frozen[i].volume == rule[i].volume);
        CHECK(q(frozen[i].u) == (i == 0 ? first : q(rule[i].u)));
        CHECK(q(frozen[i].v) == q(rule[i].v));
        CHECK(q(frozen[i].depth) == q(rule[i].depth));
        CHECK(q(frozen[i].weight) == original_weight);
    }
    const auto& conditioning = certified.conditioning();
    CHECK(take(conditioning.info()).rank == columns);
    const auto residual = q(take(conditioning.congruence_bound(caps.arithmetic.algebra, plan)));
    CHECK(residual <= Q("1/100000000000000000000"));
    check_generated_required(conditioning, 216, independent_transform, scratch, plan);
    const auto source_basis = take(take(conditioning.reference_basis())->reference_identity());
    auto patch =
        take(PgPreparedPatch::prepare(std::move(certified), caps.arithmetic, scratch, plan));
    CHECK(take(patch.info()).geometry == source_basis);
    const auto available_rule = patch.fixed_rule();
    REQUIRE(available_rule.has_value());
    const auto owned_rule = available_rule.value_or(std::span<const PgOwnedDiscreteVolumePoint>{});
    REQUIRE(owned_rule.size() == rule.size());
    const auto& p = owned_rule.front();
    const auto& exact = rule.front();
    const PgReferencePoint query = PgReferenceVolumePoint{p.volume, p.u, p.v, p.depth};
    const auto evaluated =
        take(patch.evaluate_reference(std::span{&query, 1}, caps.arithmetic, plan));
    const auto root = exact.root;
    const std::array<Q, 3> expected{2 * exact.s - exact.zeta * exact.s / (16 * root), 2 * exact.t,
                                    exact.s * exact.s / 8 + exact.zeta / (2 * root)};
    std::array<Q, 3> rebuilt{};
    Q translated{};
    for (std::size_t j = 0; j < columns; ++j) {
        const auto& n = evaluated.columns.values()[j].value;
        const Q scalar = q(n.rational) + q(n.radical) * root;
        for (std::size_t d = 0; d < rebuilt.size(); ++d)
            rebuilt[d] += scalar * q(patch.reference_coordinates()[j][d]);
        translated += scalar * q(patch.translation()[j]);
    }
    CHECK(rebuilt == expected);
    CHECK(translated == 1);
    const auto adjacent_patch = take(
        PgPreparedPatch::prepare(std::move(adjacent_certificate), caps.arithmetic, scratch, plan));
    check_curved_trace(f, patch, adjacent_patch, caps.arithmetic, scratch, plan);
    check_fixed_projection(patch, rule, caps.arithmetic, plan);
    check_projected_weight_admission(patch, caps.arithmetic, plan);
    check_curved_area_projection(f, patch, caps.arithmetic, plan);
}

TEST_CASE("PG irrational curved samples enclose signed Gram entries",
          "[SH-D03][SH-D07][N4-02][N4-08][verification]") {
    Fixture f;
    core::Workspace scratch{core::PlannedBuffer<std::byte>{}};
    core::ResourcePlan setup{profile}, plan{profile};
    auto basis = parabolic_basis(f, scratch, setup);
    const auto& partition = take(basis.source()).partition;
    const auto volumes = take(partition.volumes());
    const auto triangles = take(partition.triangles());
    const auto pieces = take(partition.pieces());
    const auto vertices = take(partition.vertices());
    const auto layers = take(partition.layers());
    const auto third = f.r("1", "3"), half = f.r("1", "2");
    std::vector<PgDiscreteVolumePoint> points;
    points.reserve(volumes.size());
    std::vector<std::vector<Term>> raw = basis_terms(*take(basis.rank()), 3, setup);
    const std::array<std::pair<std::size_t, std::size_t>, 2> checks{{{0, 0}, {1, 2}}};
    std::array<std::array<Q, 2>, checks.size()> bounds{};
    std::size_t irrational = 0;
    const I oracle_denominator = I{1} << 256;
    const Dec scale = oracle_denominator.convert_to<Dec>();
    for (std::size_t id = 0; id < volumes.size(); ++id) {
        const auto& volume = volumes[id];
        const auto& triangle = triangles[static_cast<std::size_t>(volume.triangle)];
        const auto& layer = layers[static_cast<std::size_t>(volume.layer)];
        const auto cell = pieces[static_cast<std::size_t>(triangle.piece)].cell;
        std::array<Q, 2> sum{};
        for (const auto v : triangle.vertices) {
            const auto& vertex = vertices[static_cast<std::size_t>(v)];
            sum[0] += q(vertex.s);
            sum[1] += q(vertex.t);
        }
        const Q s = sum[0] / 3, t = sum[1] / 3;
        const Q z = (q(layer.lower) + q(layer.upper)) / 8;
        const Q g = 1 + s * s / 64;
        const I n = numerator(g), d = denominator(g);
        const I a = sqrt(n), b = sqrt(d);
        irrational += a * a != n || b * b != d ? 1U : 0U;
        const Dec root = sqrt(n.convert_to<Dec>() / d.convert_to<Dec>());
        const I k = floor(root * scale).convert_to<I>();
        const Q lo = Q{k} / oracle_denominator, hi = Q{k + 1} / oracle_denominator;
        REQUIRE(lo * lo <= g);
        REQUIRE(hi * hi > g);

        using Radical = std::array<Q, 2>;
        const auto product = [&](const Radical& x, const Radical& y) -> Radical {
            return {x[0] * y[0] + g * x[1] * y[1], x[0] * y[1] + x[1] * y[0]};
        };
        const Q u = 4 * s - (2 * cell.first + 1), v = 4 * t - (2 * cell.second + 1);
        const auto value = [&](std::size_t col) -> Radical {
            std::array<Q, 3> field{};
            for (const auto& term : raw[col])
                if (term.row.cell == cell)
                    field[static_cast<std::size_t>(term.row.field)] +=
                        term.value * power(u, term.row.powers.first) *
                        power(v, term.row.powers.second);
            return {field[0] + z * z * field[2] / 2, z * field[1] / (g * g)};
        };
        const Q factor = q(triangle.twice_area) * (q(layer.upper) - q(layer.lower)) / 8;
        const Radical measure{-factor * z / (8 * g), factor};
        for (std::size_t j = 0; j < checks.size(); ++j) {
            const auto entry =
                product(measure, product(value(checks[j].first), value(checks[j].second)));
            const Q x = entry[0] + entry[1] * lo, y = entry[0] + entry[1] * hi;
            bounds[j][0] += std::min(x, y);
            bounds[j][1] += std::max(x, y);
        }
        points.push_back({{static_cast<core::Index>(id), third, third, half}, half});
    }
    REQUIRE(irrational > 0);
    const auto rank = take(basis.size());
    REQUIRE(rank == 130);
    PgDiscreteGramLimits limited{query_limits(), steps};
    limited.arithmetic.algebra.external_memory =
        ExactEstimateAdmission{std::uint64_t{32} * 1024 * 1024};
    // Preparation transfers the basis only after all admission checks succeed.
    const auto prepare_gram = [&](const PgDiscreteGramLimits& limits, core::ResourcePlan& memory) {
        return PgDiscreteGram::prepare(std::move(basis), points, {224, 152}, limits, scratch,
                                       memory);
    };
    core::ResourcePlan rejected_plan{profile};
    check_failure(prepare_gram(limited, rejected_plan), core::FailureCode::ResourceExceeded);
    CHECK(rejected_plan.allocated_bytes() == 0);
    REQUIRE(take(basis.size()) == rank);

    auto gram = take(prepare_gram({query_limits(), steps}, plan));
    const auto entries = take(gram.enclosure()).upper_triangle;
    for (std::size_t i = 0; i < checks.size(); ++i) {
        const auto& entry = entries[slot(checks[i].first, checks[i].second)];
        CHECK(q(entry.lower) <= bounds[i][0]);
        CHECK(q(entry.upper) >= bounds[i][1]);
    }
    const auto generate = [&] {
        return std::move(gram).generate_and_certify({240, 160}, {profile, steps, {}}, scratch,
                                                    plan);
    };
    check_failure(generate(), core::FailureCode::PrecisionUnresolved);
    REQUIRE(gram.enclosure());
    CHECK(gram.rule().size() == points.size());
}
