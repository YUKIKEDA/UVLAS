#include "pg_conditioning_test_support.hpp"

#include <algorithm>
#include <string>
#include <uvlas/structure/pg_discrete_gram.hpp>

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

PgReferenceBasis parabolic_basis(Fixture& f, core::Workspace& scratch, core::ResourcePlan& plan) {
    const auto zero = f.r("0"), one = f.r("1"), length = f.r("2");
    const auto coefficient = f.r("1", "16");
    const std::array<ExactInputTerm, 1> terms{{{{2, 0}, coefficient}}};
    const auto curve = take(ExactChartPolynomial::make(terms, exact_limits(), f.input_plan));
    f.graph = take(PgGraph::prepare({zero, one, zero, curve, 1, PgGraphConvention::Version1},
                                    {exact_limits(), 64, {}}, f.input_plan));
    return make_basis(f, length, scratch, plan);
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
    core::ResourcePlan rejected_plan{profile};
    check_failure(PgDiscreteGram::prepare(std::move(basis), points, {224, 152}, limited, scratch,
                                         rejected_plan),
                  core::FailureCode::ResourceExceeded);
    CHECK(rejected_plan.allocated_bytes() == 0);
    REQUIRE(take(basis.size()) == rank);

    auto gram = take(PgDiscreteGram::prepare(std::move(basis), points, {224, 152},
                                             {query_limits(), steps}, scratch, plan));
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
