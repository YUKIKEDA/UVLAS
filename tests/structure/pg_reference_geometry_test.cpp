#include "exact_polynomial_detail.hpp"
#include "pg_partition_test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <uvlas/structure/pg_reference_basis.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::structure;
using partition_test::Fixture;
using partition_test::take;
using Q = detail::Rational;
using R = ExactRational;
constexpr auto cap = partition_test::cap;

[[nodiscard]] PgReferenceLimits reference_limits() {
    auto limits = partition_test::exact_limits();
    limits.maximum_value_bits = 8192;
    limits.maximum_work_bits = 32768;
    return {limits, 10000000, 100000, {}};
}

[[nodiscard]] PgReferenceInputs prepare_inputs(Fixture& fixture, core::ResourcePlan& plan) {
    auto workspace = core::Workspace{core::PlannedBuffer<std::byte>{}};
    const auto coefficient_limits =
        PgCoefficientLimits{partition_test::exact_limits(), 2000, 20000000, 100000, {}};
    auto dictionary = take(PgCoefficientDictionary::prepare(fixture.chart_graph(),
                                                            coefficient_limits, workspace, plan));
    auto partition = take(PgIntegrationPartition::prepare(
        fixture.request(), partition_test::limits(), workspace, plan));
    auto graph = std::move(fixture.chart_graph());
    return {std::move(graph), std::move(dictionary), std::move(partition)};
}

struct Frame {
    std::array<R, 9> q;
    std::array<R, 3> origin;
    R length;
    explicit Frame(Fixture& f)
        : q{f.r("0"), f.r("0"), f.r("1"), f.r("1"), f.r("0"),
            f.r("0"), f.r("0"), f.r("1"), f.r("0")},
          origin{f.r("1"), f.r("-2"), f.r("3")},
          length{f.r("2")} {}
    [[nodiscard]] PgReferenceFrame view() const {
        return {q, origin, length};
    }
};

[[nodiscard]] core::Result<PgReferenceGeometry> prepare_candidate(PgReferenceInputs& candidate,
                                                                  const PgReferenceFrame& frame,
                                                                  const PgReferenceLimits& limits,
                                                                  core::ResourcePlan& plan) {
    // The rvalue factory transfers only after all checks succeed.
    return PgReferenceGeometry::prepare(std::move(candidate), frame, limits, plan);
}

[[nodiscard]] PgReferenceGeometry transfer(PgReferenceGeometry& source) {
    return std::move(source);
}

template <class T>
[[nodiscard]] const T& required(const std::optional<T>& value) {
    if (!value)
        throw std::runtime_error("Missing reference result");
    return *value;
}

void set_heights(Fixture& f, const std::array<std::string_view, 4>& values) {
    for (std::size_t i = 0; i < values.size(); ++i) {
        const auto split = values[i].find('/');
        const auto n = values[i].substr(0, split);
        const auto d = split == std::string_view::npos ? "1" : values[i].substr(split + 1);
        f.heights[i] = f.r(n, d);
    }
}
[[nodiscard]] const Q& rational(const R& value) {
    return detail::ExactFactory::value(value);
}
[[nodiscard]] Q rational_value(const PgReferenceScalar& value) {
    // Used only for planar tests, for which the known g=1.
    return rational(value.rational) + rational(value.radical);
}

void bounds(Fixture& f, std::string_view maximum_n, std::string_view maximum_d) {
    const auto point = [&](bool x, bool y) {
        return PgPartitionPoint{x ? f.r(maximum_n, maximum_d) : f.r("0"),
                                y ? f.r(maximum_n, maximum_d) : f.r("0")};
    };
    f.domain.clear();
    f.thickness.clear();
    for (const auto side : std::array<std::array<bool, 2>, 4>{
             {{false, false}, {true, false}, {true, true}, {false, true}}}) {
        f.domain.push_back(point(side[0], side[1]));
        f.thickness.push_back(point(side[0], side[1]));
    }
}

void set_graph(Fixture& f, core::Index n, const R& w, const R& a,
               std::span<const ExactInputTerm> terms) {
    const auto zero = f.r("0");
    const auto limits = partition_test::exact_limits();
    const auto f_poly = take(ExactChartPolynomial::make(terms, limits, f.input_plan));
    f.graph = take(PgGraph::prepare({zero, w, a, f_poly, n, PgGraphConvention::Version1},
                                    {limits, 64, {}}, f.input_plan));
}

[[nodiscard]] R fraction(Fixture& fixture, const Q& value) {
    return fixture.r(boost::multiprecision::numerator(value).str(),
                     boost::multiprecision::denominator(value).str());
}

struct Barycentric {
    R u, v;
};

[[nodiscard]] std::optional<Barycentric> inside(Fixture& f, const PgIntegrationPartition& partition,
                                                const PgPartitionTriangle& tri, const Q& s,
                                                const Q& t) {
    // Invert the affine 3x3 homogeneous simplex matrix by signed subareas.
    const auto points = take(partition.vertices());
    const auto& a = points[static_cast<std::size_t>(tri.vertices[0])];
    const auto& b = points[static_cast<std::size_t>(tri.vertices[1])];
    const auto& c = points[static_cast<std::size_t>(tri.vertices[2])];
    const auto signed_area = [&](const PgPartitionPoint& p, const PgPartitionPoint& q) -> Q {
        return rational(p.s) * rational(q.t) - rational(p.t) * rational(q.s) +
               (rational(p.t) - rational(q.t)) * s + (rational(q.s) - rational(p.s)) * t;
    };
    auto u = signed_area(c, a), v = signed_area(a, b), remaining = signed_area(b, c);
    if (u < 0 || v < 0 || remaining < 0)
        return std::nullopt;
    const Q total = u + v + remaining;
    return Barycentric{fraction(f, u / total), fraction(f, v / total)};
}

void check_exact_scalar(const PgReferenceScalar& actual, std::string_view a, std::string_view b) {
    CHECK(rational(actual.rational) == Q(std::string(a)));
    CHECK(rational(actual.radical) == Q(std::string(b)));
}
void check_exact_vector(const PgReferenceVector& actual, const std::array<std::string_view, 3>& a,
                        const std::array<std::string_view, 3>& b) {
    for (std::size_t i = 0; i < actual.size(); ++i) {
        CAPTURE(i);
        check_exact_scalar(actual[i], a[i], b[i]);
    }
}

[[nodiscard]] PgCoefficientColumn find_column(const PgCoefficientDictionary& dictionary,
                                              PgCoefficientKind kind, std::uint8_t component) {
    const auto columns = take(dictionary.columns());
    for (std::size_t i = 0; i < columns.size(); ++i)
        if (columns[i].kind == kind && columns[i].component == component)
            return {static_cast<core::Index>(i)};
    FAIL("Reference test dictionary column is absent");
}

[[nodiscard]] PgRankPreparationLimits rank_limits() {
    return {std::uint64_t{48} << 30U, 10000000000000ULL, 128, {}};
}

[[nodiscard]] PgRankPreparation make_rank(const PgReferenceGeometry& geometry,
                                          core::Workspace& workspace, core::ResourcePlan& plan) {
    const auto source = take(geometry.source());
    return take(PgRankPreparation::prepare(source.dictionary, source.partition, rank_limits(),
                                           workspace, plan));
}

[[nodiscard]] core::Result<PgReferenceBasis> bind_basis(PgReferenceGeometry& geometry,
                                                        PgRankPreparation& rank,
                                                        const PgRankPreparationLimits& limits,
                                                        core::Workspace& workspace,
                                                        core::ResourcePlan& plan) {
    return PgReferenceBasis::prepare(std::move(geometry), std::move(rank), limits, workspace, plan);
}

[[nodiscard]] PgRankBasisColumn selected_required(const PgRankPreparation& rank,
                                                  std::uint32_t target) {
    const auto selected = take(rank.selected());
    const auto position = std::ranges::find(selected, target);
    REQUIRE(position != selected.end());
    return {static_cast<core::Index>(position - selected.begin())};
}

[[nodiscard]] PgRankBasisColumn centre_c(const PgRankPreparation& rank) {
    const auto selected = take(rank.selected());
    const auto sources = take(rank.provenance());
    for (std::size_t i = 0; i < selected.size(); ++i) {
        if (selected[i] < 33)
            continue;
        const auto& origin = sources[selected[i] - 33];
        if (origin.kind == PgCoefficientKind::BasicC && origin.first == 1 && origin.second == 1)
            return {static_cast<core::Index>(i)};
    }
    FAIL("Selected basis has no centre C column");
}

struct IndependentRadical {
    Q a, b;
};

[[nodiscard]] IndependentRadical oracle_dot(const PgReferenceVector& gradient,
                                            const std::array<IndependentRadical, 3>& column,
                                            const Q& g) {
    IndependentRadical product{};
    for (std::size_t i = 0; i < gradient.size(); ++i) {
        const auto& x = rational(gradient[i].rational);
        const auto& y = rational(gradient[i].radical);
        product.a += x * column[i].a + g * y * column[i].b;
        product.b += x * column[i].b + y * column[i].a;
    }
    return product;
}

double r5_quadrature_volume(Fixture& fixture, const PgReferenceGeometry& geometry,
                            std::size_t volume_count, core::ResourcePlan& plan) {
    constexpr std::array<std::array<std::string_view, 2>, 9> plane{{
        {"0.015919880246186955082212", "0.040637194180787205985946"},
        {"0.081984446336682102850285", "0.090324080347428702029236"},
        {"0.193314283649704801345649", "0.130305348201467731159371"},
        {"0.337873288298095535480731", "0.156173538520001420034315"},
        {"0.500000000000000000000000", "0.165119677500629881582264"},
        {"0.662126711701904464519269", "0.156173538520001420034315"},
        {"0.806685716350295198654351", "0.130305348201467731159371"},
        {"0.918015553663317897149715", "0.090324080347428702029236"},
        {"0.984080119753813044917788", "0.040637194180787205985946"},
    }};
    const auto exact_decimal = [&](std::string_view value) {
        return take(R::from_decimal(value, partition_test::exact_limits(), fixture.input_plan));
    };
    const std::array depth{exact_decimal("0.211324865405187117745426"),
                           exact_decimal("0.788675134594812882254574")};
    double result = 0.0;
    for (const auto& radial : plane) {
        const auto u = exact_decimal(radial[0]);
        const Q first_weight = rational(exact_decimal(radial[1])) * (1 - rational(u));
        for (const auto& transverse : plane) {
            const auto v =
                fraction(fixture, (1 - rational(u)) * rational(exact_decimal(transverse[0])));
            const double weight =
                Q(first_weight * rational(exact_decimal(transverse[1])) / 2).convert_to<double>();
            for (std::size_t i = 0; i < volume_count; ++i)
                for (const auto& d : depth) {
                    const auto sample = take(geometry.evaluate(
                        PgReferenceVolumePoint{static_cast<core::Index>(i), u, v, d}, std::nullopt,
                        reference_limits(), plan));
                    const auto& jacobian =
                        std::get<PgReferenceVolumeMeasure>(sample.measure).jacobian;
                    result += weight * (rational(jacobian.rational).convert_to<double>() +
                                        rational(jacobian.radical).convert_to<double>() *
                                            std::sqrt(rational(sample.g).convert_to<double>()));
                }
        }
    }
    return result;
}

} // namespace

TEST_CASE("PG reference geometry matches the frozen swept-taper R5 footprint",
          "[GEO-01][SH-D02][SH-D07][verification]") {
    for (const core::Index resolution : {1, 2}) {
        for (const auto denominator : {"100", "10000"}) {
            CAPTURE(resolution, denominator);
            Fixture fixture;
            // R5, polynomial_reference_surface.coefficients((.4,.2,0.,.1,-.06,.025)).
            const auto s2 = fixture.r("1", "5"), t2 = fixture.r("1", "10");
            const auto s3 = fixture.r("1", "30"), t3 = fixture.r("-1", "50");
            const auto mixed = fixture.r("1", "40");
            const std::array<ExactInputTerm, 5> height_terms{
                {{{2, 0}, s2}, {{0, 2}, t2}, {{3, 0}, s3}, {{0, 3}, t3}, {{2, 2}, mixed}}};
            const auto zero = fixture.r("0"), one = fixture.r("1");
            set_graph(fixture, resolution, one, zero, height_terms);
            fixture.domain.clear();
            fixture.domain.push_back(fixture.point("0", "0"));
            fixture.domain.push_back({fixture.r("1"), fixture.r("3", "10")});
            fixture.domain.push_back({fixture.r("1"), fixture.r("9", "10")});
            fixture.domain.push_back(fixture.point("0", "1"));
            const Q base_thickness{std::string("1/") + denominator};
            for (std::size_t i = 0; i < fixture.thickness.size(); ++i) {
                const auto& vertex = fixture.thickness[i];
                const Q thickness =
                    base_thickness * (1 + rational(vertex.s) / 10 + rational(vertex.t) / 20);
                fixture.heights[i] = fraction(fixture, thickness);
            }
            fixture.layers.clear();
            fixture.layers.push_back(
                {core::LayerId{11}, fixture.zeta.front(), fixture.zeta.back()});

            std::array<R, 9> rotation;
            for (const auto diagonal : {0U, 4U, 8U})
                rotation[diagonal] = fixture.r("1");
            const std::array<R, 3> origin;
            core::ResourcePlan plan{cap};
            auto inputs = prepare_inputs(fixture, plan);
            const auto geometry =
                take(prepare_candidate(inputs, {rotation, origin, one}, reference_limits(), plan));
            const auto& partition = take(geometry.source()).partition;
            const auto vertices = take(partition.vertices());
            const auto triangles = take(partition.triangles());
            const auto volumes = take(partition.volumes());
            REQUIRE(volumes.size() == triangles.size());
            const auto third = fixture.r("1", "3"), depth = fixture.r("3", "4");
            Q area{}, first_s{}, first_t{};
            for (std::size_t v = 0; v < volumes.size(); ++v) {
                const auto& triangle = triangles[static_cast<std::size_t>(volumes[v].triangle)];
                Q s{}, t{};
                for (const auto vertex : triangle.vertices) {
                    const auto& point = vertices[static_cast<std::size_t>(vertex)];
                    s += rational(point.s) / 3;
                    t += rational(point.t) / 3;
                }
                const Q triangle_area = rational(triangle.twice_area) / 2;
                area += triangle_area;
                first_s += triangle_area * s;
                first_t += triangle_area * t;

                const Q hs = 2 * s / 5 + s * s / 10 + s * t * t / 20;
                const Q ht = t / 5 - 3 * t * t / 50 + s * s * t / 20;
                const Q hss = Q("2/5") + s / 5 + t * t / 20;
                const Q htt = Q("1/5") - 3 * t / 25 + s * s / 20;
                const Q hst = s * t / 10;
                const Q g = 1 + hs * hs + ht * ht;
                const Q h = base_thickness * (1 + s / 10 + t / 20);
                const Q z = h / 4;
                const Q gaussian_numerator = hss * htt - hst * hst;
                const Q mean_numerator =
                    (1 + hs * hs) * htt - 2 * hs * ht * hst + (1 + ht * ht) * hss;
                // det[dX/ds,dX/dt,dX/dz] from the two shape-operator invariants.
                const Q determinant_rational = -z * mean_numerator / g;
                const Q determinant_radical = 1 + z * z * gaussian_numerator / (g * g);
                const auto sample = take(geometry.evaluate(
                    PgReferenceVolumePoint{static_cast<core::Index>(v), third, third, depth},
                    std::nullopt, reference_limits(), plan));
                CHECK(rational(sample.s) == s);
                CHECK(rational(sample.t) == t);
                CHECK(rational(sample.g) == g);
                CHECK(rational(sample.full_thickness_metres) == h);
                CHECK(rational(sample.position_metres[0].rational) == s);
                CHECK(rational(sample.position_metres[0].radical) == -z * hs / g);
                CHECK(rational(sample.position_metres[1].radical) == -z * ht / g);
                CHECK(rational(sample.position_metres[2].radical) == z / g);
                const auto& measure = std::get<PgReferenceVolumeMeasure>(sample.measure).jacobian;
                // The unit depth spans all thickness and the triangle map has area factor 2A.
                CHECK(rational(measure.rational) == determinant_rational * h * 2 * triangle_area);
                CHECK(rational(measure.radical) == determinant_radical * h * 2 * triangle_area);
            }
            CHECK(area == Q("4/5"));
            CHECK(first_s == Q("11/30"));
            CHECK(first_t == Q("131/300"));
            // Independently integrated over chart strips and exact depth moments.
            constexpr double full_volume = 0.008930137852542443516987;
            constexpr double thin_volume = 0.000089301323707466072234;
            const auto expected = base_thickness == Q("1/100") ? full_volume : thin_volume;
            const auto integrated = r5_quadrature_volume(fixture, geometry, volumes.size(), plan);
            CHECK(std::abs(integrated - expected) <= 1e-12 * expected);
        }
    }
}

TEST_CASE("PG physical map and fourth-order jet match an independent quartic reference",
          "[SH-D01][SH-D02][SH-D07][N4-06][V-PG][verification]") {
    Fixture f;
    bounds(f, "1", "2");
    const auto w = f.r("4", "5"), a = f.r("3", "20");
    const auto three_quarters = f.r("3", "4"), half = f.r("1", "2");
    const std::array<ExactInputTerm, 2> terms{{{{3, 1}, three_quarters}, {{2, 2}, half}}};
    set_graph(f, 1, w, a, terms);
    set_heights(f, {"1/2", "3/5", "13/20", "11/20"});
    f.zeta[1] = f.r("0");
    const Frame frame{f};
    core::ResourcePlan plan{cap};
    auto inputs = prepare_inputs(f, plan);
    auto geometry = take(prepare_candidate(inputs, frame.view(), reference_limits(), plan));
    REQUIRE_FALSE(inputs.graph.view());
    const auto source = take(geometry.source());
    const auto& partition = source.partition;
    const auto triangles = take(partition.triangles());
    const auto volumes = take(partition.volumes());
    const auto depth = f.r("1", "3");
    for (std::size_t i = 0; i < volumes.size(); ++i) {
        if (volumes[i].layer != 1)
            continue;
        const auto bary =
            inside(f, partition, triangles[static_cast<std::size_t>(volumes[i].triangle)], Q("1/3"),
                   Q("1/5"));
        if (!bary)
            continue;
        const auto sample = take(geometry.evaluate(
            PgReferenceVolumePoint{static_cast<core::Index>(i), bary->u, bary->v, depth},
            find_column(source.dictionary, PgCoefficientKind::NormalD, 2), reference_limits(),
            plan));
        CHECK(rational(sample.full_thickness_metres) == Q("44/75"));
        CHECK(rational(sample.thickness_gradient_metres[0]) == Q("1/5"));
        CHECK(rational(sample.thickness_gradient_metres[1]) == Q("1/10"));
        CHECK(rational(sample.g) == Q("23217889/23040000"));
        constexpr std::array<std::string_view, 15> height{
            "7/900", "19/300", "1/20", "17/50", "23/60", "1/9", "9/10", "19/10",
            "2/3",   "0",      "0",    "9/2",   "2",     "0",   "0"};
        for (std::size_t j = 0; j < height.size(); ++j)
            CHECK(rational(sample.height_chart_jet[j]) == Q(height[j]));
        constexpr std::array<std::string_view, 3> zero{"0", "0", "0"};
        check_exact_vector(sample.normal, zero,
                           {"23040000/23217889", "-1824000/23217889", "-878400/23217889"});
        check_exact_vector(sample.position_metres, {"457/450", "-211/150", "17/5"},
                           {"2252800/23217889", "-535040/69653667", "-85888/23217889"});
        check_exact_vector(sample.jacobian_columns_metres[0], {"19/150", "8/5", "0"},
                           {"15452557670400/539070369616321", "-23453063126400/539070369616321",
                            "-51915054633760/1617211108848963"});
        check_exact_vector(sample.jacobian_columns_metres[1], {"1/10", "3/10", "2"},
                           {"6869079424000/539070369616321", "-76820218100800/1617211108848963",
                            "-20827247934640/4851633326546889"});
        check_exact_vector(sample.jacobian_columns_metres[2], zero,
                           {"6758400/23217889", "-535040/23217889", "-257664/23217889"});

        // For this raw column the scalar director is n_z. Its independent derivative
        // is -(v dot v_a)/g^(3/2), calculated here from the independent height literals.
        REQUIRE(sample.column);
        const Q g{"23217889/23040000"}, vx{"-19/240"}, vy{"-61/1600"};
        const Q vxs{"-17/40"}, vyt{Q("3/20") * Q("23/48") - Q("1/9")};
        const Q vxt{"-23/48"}, vys{Q("3/20") * Q("17/40") - Q("23/60")};
        const Q z{"11/225"};
        const auto& basis = required(sample.column);
        CHECK(rational(basis.value.rational) == 0);
        CHECK(rational(basis.value.radical) == z / g);
        CHECK(rational(basis.chart_derivatives[0].radical) ==
              Q("1/60") / g - z * (vx * vxs + vy * vys) / (g * g));
        CHECK(rational(basis.chart_derivatives[1].radical) ==
              Q("1/120") / g - z * (vx * vxt + vy * vyt) / (g * g));
        CHECK(rational(basis.chart_derivatives[2].radical) == Q("11/75") / g);
        return;
    }
    FAIL("Independent chart point is outside the product partition");
}

TEST_CASE("PG evaluates both material faces and real lateral faces for every adopted grid",
          "[SH-D02][SH-D07][N4-08][V-PG][verification]") {
    for (const core::Index n : {1, 2, 4}) {
        for (const bool cut : {false, true}) {
            CAPTURE(n, cut);
            Fixture f{n};
            if (cut)
                partition_test::incidence_fixture(f);
            set_heights(f, {"1", "2", "1", "3"});
            f.zeta[1] = f.r("0");
            const Frame frame{f};
            core::ResourcePlan plan{cap};
            auto inputs = prepare_inputs(f, plan);
            auto geometry = take(prepare_candidate(inputs, frame.view(), reference_limits(), plan));
            const auto& partition = take(geometry.source()).partition;
            const auto zero = f.r("0");
            const auto one = f.r("1");
            const auto volumes = take(partition.volumes());
            const auto triangles = take(partition.triangles());
            const auto volume = take(geometry.evaluate(PgReferenceVolumePoint{0, zero, zero, one},
                                                       std::nullopt, reference_limits(), plan));
            const auto triangle = static_cast<std::size_t>(volumes[0].triangle);
            REQUIRE(rational(volume.g) == 1);
            const Q physical_det = 2 * rational(volume.full_thickness_metres);
            CHECK(rational_value(volume.determinant_metres_cubed) == physical_det);
            CHECK(rational_value(std::get<PgReferenceVolumeMeasure>(volume.measure).jacobian) ==
                  physical_det * rational(triangles[triangle].twice_area));
            CHECK(rational(**geometry.determinant_lower_bound()) <= physical_det);

            const auto horizontal = take(partition.horizontal_faces());
            std::optional<PgReferenceValue> upper, lower;
            for (std::size_t i = 0; i < horizontal.size(); ++i) {
                const auto& face = horizontal[i];
                if (face.triangle != volumes[0].triangle)
                    continue;
                if (face.layer == 0 && face.side == PgPartitionFaceSide::Upper)
                    upper.emplace(take(geometry.evaluate(
                        PgReferenceHorizontalPoint{static_cast<core::Index>(i), zero, zero},
                        std::nullopt, reference_limits(), plan)));
                if (face.layer == 1 && face.side == PgPartitionFaceSide::Lower)
                    lower.emplace(take(geometry.evaluate(
                        PgReferenceHorizontalPoint{static_cast<core::Index>(i), zero, zero},
                        std::nullopt, reference_limits(), plan)));
            }
            REQUIRE(upper);
            REQUIRE(lower);
            const auto& upper_result = required(upper);
            const auto& up = std::get<PgReferenceHorizontalMeasure>(upper_result.measure);
            CHECK(up.side == PgPartitionFaceSide::Upper);
            const auto& upper_measure = up.area;
            const auto& lower_result = required(lower);
            const auto& lo = std::get<PgReferenceHorizontalMeasure>(lower_result.measure);
            CHECK(lo.side == PgPartitionFaceSide::Lower);
            const auto& lower_measure = lo.area;
            const Q area = 4 * rational(triangles[triangle].twice_area);
            CHECK(rational_value(upper_measure.outward[0]) == area);
            CHECK(rational_value(lower_measure.outward[0]) == -area);
            CHECK(rational_value(upper_measure.squared_norm) == area * area);
            CHECK(rational_value(lower_measure.squared_norm) == area * area);
            CHECK(upper_result.layer != lower_result.layer);
            CHECK(rational_value(upper_result.position_metres[0]) ==
                  rational_value(lower_result.position_metres[0]));

            const auto faces = take(partition.lateral_faces());
            const auto edges = take(partition.edges());
            const auto vertices = take(partition.vertices());
            for (std::size_t i = 0; i < faces.size(); ++i) {
                const auto& face = faces[i];
                const auto& edge = edges[static_cast<std::size_t>(face.edge)];
                const auto& a = vertices[static_cast<std::size_t>(
                    edge.vertices[face.canonical_direction ? 0 : 1])];
                const auto& b = vertices[static_cast<std::size_t>(
                    edge.vertices[face.canonical_direction ? 1 : 0])];
                const auto side = take(geometry.evaluate(
                    PgReferenceLateralPoint{static_cast<core::Index>(i), zero, zero}, std::nullopt,
                    reference_limits(), plan));
                const auto& measure = std::get<PgReferenceLateralMeasure>(side.measure).area;
                const auto h = rational(side.full_thickness_metres);
                const Q y = h * (rational(b.t) - rational(a.t));
                const Q z = -h * (rational(b.s) - rational(a.s));
                CHECK(rational_value(measure.outward[0]) == 0);
                CHECK(rational_value(measure.outward[1]) == y);
                CHECK(rational_value(measure.outward[2]) == z);
                CHECK(rational_value(measure.squared_norm) == y * y + z * z);
            }
        }
    }
}

TEST_CASE("PG raw coefficient jets and thickness gradients retain separate one-sided cells",
          "[SH-D01][SH-D02][SH-D07][V-PG][verification]") {
    Fixture f;
    set_graph(f, 1, f.r("4", "5"), f.r("3", "20"), {});
    set_heights(f, {"1", "2", "1", "3"});
    f.zeta[1] = f.r("0");
    const Frame frame{f};
    core::ResourcePlan plan{cap};
    auto inputs = prepare_inputs(f, plan);
    auto geometry = take(prepare_candidate(inputs, frame.view(), reference_limits(), plan));
    const auto source = take(geometry.source());
    const auto columns = take(source.dictionary.columns());
    std::optional<PgCoefficientColumn> column;
    for (std::size_t i = 0; i < columns.size(); ++i) {
        const auto& item = columns[i];
        if (item.kind == PgCoefficientKind::Paired && item.axis == PgPairedAxis::S &&
            item.first == 0 && item.second == 0 && item.component == 2)
            column = PgCoefficientColumn{static_cast<core::Index>(i)};
    }
    REQUIRE(column);
    const auto volumes = take(source.partition.volumes());
    const auto triangles = take(source.partition.triangles());
    const auto pieces = take(source.partition.pieces());
    const auto half = f.r("1", "2");
    std::array<bool, 2> seen{};
    for (std::size_t i = 0; i < volumes.size(); ++i) {
        const auto& volume = volumes[i];
        if (volume.layer != 1)
            continue;
        const auto& tri = triangles[static_cast<std::size_t>(volume.triangle)];
        const auto& piece = pieces[static_cast<std::size_t>(tri.piece)];
        if (piece.cell.first != piece.thickness_triangle)
            continue;
        const auto bary = inside(f, source.partition, tri, Q("1/2"), Q("1/2"));
        if (!bary)
            continue;
        const auto sample = take(geometry.evaluate(
            PgReferenceVolumePoint{static_cast<core::Index>(i), bary->u, bary->v, half}, column,
            reference_limits(), plan));
        const bool after = piece.cell.first == 1;
        seen[after ? 1 : 0] = true;
        REQUIRE(sample.column);
        const auto& r = required(sample.column).numerators[0];
        CHECK(rational(r[0]) == Q("1/8"));
        CHECK(rational(r[1]) == Q("1/2"));
        CHECK(rational(r[3]) == (after ? 0 : 1));
        CHECK(rational(sample.full_thickness_metres) == 1);
        CHECK(rational(sample.thickness_gradient_metres[0]) == (after ? -2 : 1));
        CHECK(rational(sample.thickness_gradient_metres[1]) == (after ? 2 : -1));
        CHECK(rational_value(sample.jacobian_columns_metres[0][0]) ==
              (after ? Q("-1/2") : Q("1/4")));
    }
    CHECK(seen[0]);
    CHECK(seen[1]);
}

TEST_CASE("PG strict normal tube bounds distinguish uncertainty and observed inversion",
          "[SH-D02][SH-D07][N4-08][V-PG][negative]") {
    Fixture f;
    const auto one = f.r("1"), zero = f.r("0"), two = f.r("2");
    const std::array<ExactInputTerm, 1> terms{{{{2, 0}, two}}};
    set_graph(f, 1, one, zero, terms);
    set_heights(f, {"2", "2", "2", "2"});
    f.zeta[1] = f.r("0");
    Frame frame{f};
    core::ResourcePlan plan{cap};
    auto inputs = prepare_inputs(f, plan);
    const auto assert_code = [](const auto& actual, const auto code) {
        REQUIRE_FALSE(actual);
        CHECK(actual.error().code == code);
    };
    frame.length = f.r("4"); // H f_xx / (2L) == 1.
    assert_code(prepare_candidate(inputs, frame.view(), reference_limits(), plan),
                core::FailureCode::PrecisionUnresolved);
    REQUIRE(inputs.graph.view());
    REQUIRE(inputs.partition.identity());
    frame.length = f.r("2");
    assert_code(prepare_candidate(inputs, frame.view(), reference_limits(), plan),
                core::FailureCode::PrecisionUnresolved);
    const auto triangles = take(inputs.partition.triangles());
    const auto volumes = take(inputs.partition.volumes());
    for (std::size_t i = 0; i < volumes.size(); ++i) {
        if (volumes[i].layer != 1)
            continue;
        const auto bary =
            inside(f, inputs.partition, triangles[static_cast<std::size_t>(volumes[i].triangle)],
                   Q("0"), Q("0"));
        if (!bary)
            continue;
        const auto half = f.r("1", "2");
        const auto quarter = f.r("1", "4");
        assert_code(PgReferenceGeometry::diagnose_point(
                        inputs, frame.view(),
                        PgReferenceVolumePoint{static_cast<core::Index>(i), bary->u, bary->v, one},
                        reference_limits(), plan),
                    core::FailureCode::GeometryInvalid);
        assert_code(PgReferenceGeometry::diagnose_point(
                        inputs, frame.view(),
                        PgReferenceVolumePoint{static_cast<core::Index>(i), bary->u, bary->v, half},
                        reference_limits(), plan),
                    core::FailureCode::GeometryInvalid);
        CHECK(PgReferenceGeometry::diagnose_point(
            inputs, frame.view(),
            PgReferenceVolumePoint{static_cast<core::Index>(i), bary->u, bary->v, quarter},
            reference_limits(), plan));
        frame.length = f.r("401", "100");
        const auto allowed =
            take(prepare_candidate(inputs, frame.view(), reference_limits(), plan));
        const Q theta_squared{Q("400/401") * Q("400/401")};
        const Q gap{1 - theta_squared};
        CHECK(rational(**allowed.determinant_lower_bound()) ==
              Q("401/200") * Q("401/200") * gap * gap);
        return;
    }
    FAIL("The origin is outside the actual product partition");
}

TEST_CASE("PG exact frames and reference candidates survive rejected preparation or point requests",
          "[SH-D02][N4-02][N4-08][V-PG][negative][property]") {
    Fixture f;
    Frame frame{f};
    core::ResourcePlan plan{cap};
    auto inputs = prepare_inputs(f, plan);
    const auto initial_identity = take(inputs.partition.identity());
    frame.q[2] = f.r("-1");
    const auto reflected = prepare_candidate(inputs, frame.view(), reference_limits(), plan);
    REQUIRE_FALSE(reflected);
    CHECK(reflected.error().code == core::FailureCode::GeometryInvalid);
    frame.q[2] = f.r("2");
    const auto stretched = prepare_candidate(inputs, frame.view(), reference_limits(), plan);
    REQUIRE_FALSE(stretched);
    CHECK(stretched.error().code == core::FailureCode::GeometryInvalid);
    frame.q[2] = f.r("1");
    const auto source_liveness = take(inputs.graph.view()).identity;
    const auto w = f.r("1"), a = f.r("1", "10");
    set_graph(f, 1, w, a, {});
    std::swap(inputs.graph, f.chart_graph());
    const auto stale = prepare_candidate(inputs, frame.view(), reference_limits(), plan);
    REQUIRE_FALSE(stale);
    CHECK(stale.error().code == core::FailureCode::StaleGeneration);
    std::swap(inputs.graph, f.chart_graph());
    CHECK(take(inputs.graph.view()).identity == source_liveness);
    CHECK(take(inputs.partition.identity()) == initial_identity);

    auto limits = reference_limits();
    limits.cancellation = partition_test::requested_stop_token();
    const auto previous_plan = plan.allocated_bytes();
    const auto cancelled = prepare_candidate(inputs, frame.view(), limits, plan);
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == core::FailureCode::Cancelled);
    CHECK(plan.allocated_bytes() == previous_plan);
    CHECK(take(inputs.graph.view()).identity == source_liveness);
    const auto zero = f.r("0");
    const auto one = f.r("1");
    const auto minus = f.r("-1");
    const auto query = PgReferenceVolumePoint{0, zero, zero, one};

    limits = reference_limits();
    limits.algebra.external_memory = ExactRequireProvenBound{cap};
    const auto unsupported = prepare_candidate(inputs, frame.view(), limits, plan);
    REQUIRE_FALSE(unsupported);
    CHECK(unsupported.error().code == core::FailureCode::Unsupported);

    limits = reference_limits();
    const auto init = plan.allocated_bytes();
    auto geom = take(prepare_candidate(inputs, frame.view(), limits, plan));
    CHECK(geom.preparation().allocated_owned_payload_bytes == plan.allocated_bytes() - init);
    const auto invalid =
        geom.evaluate(PgReferenceVolumePoint{0, minus, zero, zero}, std::nullopt, limits, plan);
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == core::FailureCode::InvalidInput);
    const auto outside =
        geom.evaluate(PgReferenceVolumePoint{0, one, one, zero}, std::nullopt, limits, plan);
    REQUIRE_FALSE(outside);
    CHECK(outside.error().code == core::FailureCode::InvalidInput);
    const auto missing =
        geom.evaluate(PgReferenceLateralPoint{-1, zero, zero}, std::nullopt, limits, plan);
    REQUIRE_FALSE(missing);
    CHECK(missing.error().code == core::FailureCode::InvalidInput);

    core::ResourcePlan discovery{cap};
    const auto observed = take(geom.evaluate(query, std::nullopt, limits, discovery)).call_info;
    REQUIRE(observed.allocated_owned_payload_bytes > 0);
    REQUIRE(observed.scalar_operations > 0);
    REQUIRE(observed.partition_visits == 1);
    limits.maximum_scalar_operations = observed.scalar_operations;
    limits.maximum_partition_visits = observed.partition_visits;
    core::ResourcePlan exact{observed.allocated_owned_payload_bytes};
    const auto output = take(geom.evaluate(query, std::nullopt, limits, exact));
    CHECK(exact.remaining_bytes() == 0);
    CHECK(output.call_info.allocated_owned_payload_bytes == observed.allocated_owned_payload_bytes);
    core::ResourcePlan short_plan{observed.allocated_owned_payload_bytes - 1};
    const auto short_owned = geom.evaluate(query, std::nullopt, limits, short_plan);
    REQUIRE_FALSE(short_owned);
    CHECK(short_owned.error().code == core::FailureCode::ResourceExceeded);
    --limits.maximum_scalar_operations;
    const auto short_operation = geom.evaluate(query, std::nullopt, limits, plan);
    REQUIRE_FALSE(short_operation);
    CHECK(short_operation.error().code == core::FailureCode::ResourceExceeded);
    ++limits.maximum_scalar_operations;
    --limits.maximum_partition_visits;
    const auto short_visit = geom.evaluate(query, std::nullopt, limits, plan);
    REQUIRE_FALSE(short_visit);
    CHECK(short_visit.error().code == core::FailureCode::ResourceExceeded);
    CHECK(rational_value(output.position_metres[0]) == Q("5/6"));

    auto moved = transfer(geom);
    CHECK_FALSE(geom.source());
    CHECK_FALSE(geom.determinant_lower_bound());
    CHECK_FALSE(geom.evaluate(query, std::nullopt, reference_limits(), plan));
    CHECK(moved.evaluate(query, std::nullopt, reference_limits(), plan));
}

TEST_CASE("PG reference identity follows an exact non-axis frame and stamps owned output",
          "[SH-D02][N4-06][V-PG][verification]") {
    Fixture first, second;
    const Frame frame{first};
    Frame angled{second};
    angled.q[0] = second.r("3", "5");
    angled.q[1] = second.r("-4", "5");
    angled.q[2] = second.r("0");
    angled.q[3] = second.r("4", "5");
    angled.q[4] = second.r("3", "5");
    angled.q[5] = second.r("0");
    angled.q[6] = second.r("0");
    angled.q[7] = second.r("0");
    angled.q[8] = second.r("1");
    core::ResourcePlan plan{cap};
    auto source1 = prepare_inputs(first, plan), source2 = prepare_inputs(second, plan);
    CHECK(take(source1.graph.view()).identity == take(source2.graph.view()).identity);
    CHECK(take(source1.partition.identity()) == take(source2.partition.identity()));
    auto g1 = take(
        PgReferenceGeometry::prepare(std::move(source1), frame.view(), reference_limits(), plan));
    auto g2 = take(
        PgReferenceGeometry::prepare(std::move(source2), angled.view(), reference_limits(), plan));
    CHECK(take(g1.identity()) != take(g2.identity()));
    const auto zero = second.r("0"), one = second.r("1");
    const auto out = take(g2.evaluate(PgReferenceVolumePoint{0, zero, zero, one}, std::nullopt,
                                      reference_limits(), plan));
    CHECK(out.identity == take(g2.identity()));
    CHECK(rational_value(out.jacobian_columns_metres[0][0]) == Q("6/5"));
    CHECK(rational_value(out.jacobian_columns_metres[0][1]) == Q("8/5"));
    CHECK(rational_value(out.jacobian_columns_metres[1][0]) == Q("-8/5"));
    CHECK(rational_value(out.jacobian_columns_metres[1][1]) == Q("6/5"));
    CHECK(rational_value(out.normal[2]) == 1);
    CHECK(rational_value(out.determinant_metres_cubed) == 2);
}

TEST_CASE("PG reference caps and cancels distinct graph comparison before geometry work",
          "[SH-D02][N4-08][V-PG][negative][verification]") {
    Fixture f;
    const auto width = f.r("1"), skew = f.r("0");
    // The coefficient is 2^-130, so its canonical denominator takes five 32-bit words.
    const auto height = f.r("1", "1361129467683753853853498429727072845824");
    const std::array<ExactInputTerm, 1> terms{{{{2, 0}, height}}};
    set_graph(f, 1, width, skew, terms);
    core::ResourcePlan producer_plan{cap};
    auto inputs = prepare_inputs(f, producer_plan);
    Frame frame{f};
    frame.length = f.r("-1"); // Deliberate failure immediately after input comparison.
    const auto zero = f.r("0");
    const auto query = PgReferenceVolumePoint{0, zero, zero, zero};

    // One fixed digest; each graph has one header, four cells, one term, and
    // canonical equalities for x0=0, width=1, skew=0, and height=2^-130.
    constexpr std::uint64_t comparison = 1 + 2 * (1 + 4 + 1 + 3 + 5 + 3 + 13);
    const auto check_code = [](const auto& result, core::FailureCode code) {
        if (result)
            FAIL("The negative reference comparison control succeeded");
        else
            CHECK(result.error().code == code);
    };
    auto limits = reference_limits();
    core::ResourcePlan call_plan{cap};
    for (const auto visit_cap : {comparison - 1, comparison}) {
        limits.maximum_partition_visits = visit_cap;
        const auto expected = visit_cap == comparison ? core::FailureCode::GeometryInvalid
                                                      : core::FailureCode::ResourceExceeded;
        check_code(prepare_candidate(inputs, frame.view(), limits, call_plan), expected);
        check_code(
            PgReferenceGeometry::diagnose_point(inputs, frame.view(), query, limits, call_plan),
            expected);
        CHECK(call_plan.allocated_bytes() == 0);
    }

    // A cancelled request takes precedence over a stale graph, even when no
    // comparison budget is available.
    set_graph(f, 1, width, skew, {});
    std::swap(inputs.graph, f.chart_graph());
    limits.cancellation = partition_test::requested_stop_token();
    limits.maximum_partition_visits = 0;
    check_code(prepare_candidate(inputs, frame.view(), limits, call_plan),
               core::FailureCode::Cancelled);
    check_code(PgReferenceGeometry::diagnose_point(inputs, frame.view(), query, limits, call_plan),
               core::FailureCode::Cancelled);
    CHECK(call_plan.allocated_bytes() == 0);
    std::swap(inputs.graph, f.chart_graph());

    const auto pieces = take(inputs.partition.pieces()).size();
    const auto vertices = take(inputs.partition.piece_vertex_indices()).size();
    frame.length = f.r("2");
    limits = reference_limits();
    limits.maximum_partition_visits = comparison + pieces + vertices;
    const auto prepared = take(prepare_candidate(inputs, frame.view(), limits, call_plan));
    CHECK(prepared.preparation().partition_visits == limits.maximum_partition_visits);
}

TEST_CASE("Selected PG director basis satisfies the quartic physical gradient equations",
          "[SH-D02][SH-D03][N4-02][N4-06][V-PG][verification]") {
    Fixture f;
    bounds(f, "1", "2");
    const auto w = f.r("4", "5"), a = f.r("3", "20");
    const auto three_quarters = f.r("3", "4"), half = f.r("1", "2");
    const std::array<ExactInputTerm, 2> terms{{{{3, 1}, three_quarters}, {{2, 2}, half}}};
    set_graph(f, 1, w, a, terms);
    set_heights(f, {"1/2", "3/5", "13/20", "11/20"});
    f.zeta[1] = f.r("0");
    const Frame frame{f};
    core::ResourcePlan plan{cap};
    auto inputs = prepare_inputs(f, plan);
    auto geometry = take(prepare_candidate(inputs, frame.view(), reference_limits(), plan));
    const auto id = take(geometry.identity());
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    auto rank = make_rank(geometry, workspace, plan);
    const auto column = selected_required(rank, 29);
    const auto coordinate_x = selected_required(rank, 0);
    const auto rank_size = take(rank.preparation()).rank;
    const auto basis = take(bind_basis(geometry, rank, rank_limits(), workspace, plan));
    CHECK(take(basis.size()) == rank_size);
    CHECK(take(take(basis.rank())->preparation()).rank == rank_size);
    const auto& partition = take(basis.source()).partition;
    const auto triangles = take(partition.triangles());
    const auto volumes = take(partition.volumes());
    const auto depth = f.r("1", "3");
    for (std::size_t i = 0; i < volumes.size(); ++i) {
        if (volumes[i].layer != 1)
            continue;
        const auto bary =
            inside(f, partition, triangles[static_cast<std::size_t>(volumes[i].triangle)], Q("1/3"),
                   Q("1/5"));
        if (!bary)
            continue;
        const PgReferenceVolumePoint query{static_cast<core::Index>(i), bary->u, bary->v, depth};
        auto point_plan = core::ResourcePlan{cap};
        const auto sample = take(basis.evaluate(query, column, reference_limits(), point_plan));
        CHECK(sample.reference.identity == id);
        CHECK(sample.column.value == column.value);
        CHECK(sample.reference.call_info.allocated_owned_payload_bytes ==
              point_plan.allocated_bytes());

        // at (s,t,zeta)=(1/3,1/5,1/3), with Q(x,y,z)=(z,x,y) and L=2.
        const Q g{"23217889/23040000"}, vx{"-19/240"}, vy{"-61/1600"};
        const Q vxs{"-17/40"}, vyt{Q("3/20") * Q("23/48") - Q("1/9")};
        const Q vxt{"-23/48"}, vys{Q("3/20") * Q("17/40") - Q("23/60")};
        const Q z{"11/225"};
        const std::array<Q, 3> derivatives{
            Q("1/60") / g - z * (vx * vxs + vy * vys) / (g * g),
            Q("1/120") / g - z * (vx * vxt + vy * vyt) / (g * g),
            Q("11/75") / g,
        };
        CHECK(rational(sample.reference.g) == g);
        CHECK(rational(sample.value.rational) == 0);
        CHECK(rational(sample.value.radical) == z / g);
        // Independent Jacobian from analytic height/thickness derivatives.
        const std::array<std::array<IndependentRadical, 3>, 3> jacobian{
            std::array{IndependentRadical{Q("19/150"), Q("15452557670400/539070369616321")},
                       IndependentRadical{Q("8/5"), Q("-23453063126400/539070369616321")},
                       IndependentRadical{Q("0"), Q("-51915054633760/1617211108848963")}},
            std::array{IndependentRadical{Q("1/10"), Q("6869079424000/539070369616321")},
                       IndependentRadical{Q("3/10"), Q("-76820218100800/1617211108848963")},
                       IndependentRadical{Q("2"), Q("-20827247934640/4851633326546889")}},
            std::array{IndependentRadical{Q("0"), Q("6758400/23217889")},
                       IndependentRadical{Q("0"), Q("-535040/23217889")},
                       IndependentRadical{Q("0"), Q("-257664/23217889")}},
        };
        for (std::size_t d = 0; d < derivatives.size(); ++d) {
            CHECK(rational(sample.chart_derivatives[d].rational) == 0);
            CHECK(rational(sample.chart_derivatives[d].radical) == derivatives[d]);
            const auto pullback = oracle_dot(sample.physical_gradient_per_metre, jacobian[d], g);
            CHECK(pullback.a == 0);
            CHECK(pullback.b == derivatives[d]);
        }

        const auto position_x = take(basis.evaluate(query, coordinate_x, reference_limits(), plan));
        CHECK(rational(position_x.value.rational) == Q("89/300"));
        CHECK(rational(position_x.value.radical) == 0);
        const std::array<Q, 3> x_derivatives{Q("4/5"), Q("3/20"), Q("0")};
        for (std::size_t d = 0; d < x_derivatives.size(); ++d) {
            CHECK(rational(position_x.chart_derivatives[d].rational) == x_derivatives[d]);
            CHECK(rational(position_x.chart_derivatives[d].radical) == 0);
            const auto pullback =
                oracle_dot(position_x.physical_gradient_per_metre, jacobian[d], g);
            CHECK(pullback.a == x_derivatives[d]);
            CHECK(pullback.b == 0);
        }

        auto exact_caps = reference_limits();
        exact_caps.maximum_scalar_operations = sample.reference.call_info.scalar_operations;
        exact_caps.maximum_partition_visits = sample.reference.call_info.partition_visits;
        const auto bytes = sample.reference.call_info.allocated_owned_payload_bytes;
        core::ResourcePlan equal{bytes};
        REQUIRE(basis.evaluate(query, column, exact_caps, equal));
        core::ResourcePlan short_bytes{bytes - 1};
        const auto byte_failure = basis.evaluate(query, column, exact_caps, short_bytes);
        REQUIRE_FALSE(byte_failure);
        CHECK(byte_failure.error().code == core::FailureCode::ResourceExceeded);
        --exact_caps.maximum_scalar_operations;
        const auto scalar_failure = basis.evaluate(query, column, exact_caps, plan);
        REQUIRE_FALSE(scalar_failure);
        CHECK(scalar_failure.error().code == core::FailureCode::ResourceExceeded);
        ++exact_caps.maximum_scalar_operations;
        --exact_caps.maximum_partition_visits;
        const auto visit_failure = basis.evaluate(query, column, exact_caps, plan);
        REQUIRE_FALSE(visit_failure);
        CHECK(visit_failure.error().code == core::FailureCode::ResourceExceeded);
        exact_caps.cancellation = partition_test::requested_stop_token();
        core::ResourcePlan no_allocations{0};
        const auto stopped = basis.evaluate(query, column, exact_caps, no_allocations);
        REQUIRE_FALSE(stopped);
        CHECK(stopped.error().code == core::FailureCode::Cancelled);
        CHECK(no_allocations.allocated_bytes() == 0);
        return;
    }
    FAIL("Independent point is absent");
}

TEST_CASE("Selected PG quadratic-thickness basis preserves both cell sides and face kinds",
          "[SH-D02][SH-D03][SH-D07][N4-06][V-PG][verification]") {
    Fixture f{2};
    f.zeta[1] = f.r("0");
    const Frame frame{f};
    auto inputs = prepare_inputs(f, f.input_plan);
    auto geometry = take(prepare_candidate(inputs, frame.view(), reference_limits(), f.input_plan));
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    auto rank = make_rank(geometry, workspace, f.input_plan);
    const auto column = centre_c(rank);
    const auto basis = take(bind_basis(geometry, rank, rank_limits(), workspace, f.input_plan));
    const auto& partition = take(basis.source()).partition;
    const auto triangles = take(partition.triangles());
    const auto pieces = take(partition.pieces());
    const auto volumes = take(partition.volumes());
    const auto depth = f.r("1", "2");
    std::array<bool, 2> sides{};
    for (std::size_t i = 0; i < volumes.size(); ++i) {
        if (volumes[i].layer != 1)
            continue;
        const auto& tri = triangles[static_cast<std::size_t>(volumes[i].triangle)];
        const auto& piece = pieces[static_cast<std::size_t>(tri.piece)];
        const auto side = piece.cell.first == 1 ? std::size_t{0} : std::size_t{1};
        const auto bary = inside(f, partition, tri, Q("1/2"), Q("1/3"));
        if (!bary || sides[side])
            continue;
        const auto sample = take(basis.evaluate(
            PgReferenceVolumePoint{static_cast<core::Index>(i), bary->u, bary->v, depth}, column,
            reference_limits(), f.input_plan));
        const auto derivatives = side == 0 ? Q("1/96") : Q("-1/96");
        const Q z_hat{"1/8"};
        // C is the tensor of centred linear hats on s,t. H=1m and L=2m.
        REQUIRE(rational(sample.reference.g) == 1);
        CHECK(rational_value(sample.value) == z_hat * z_hat / 3);
        CHECK(rational_value(sample.chart_derivatives[0]) == derivatives);
        CHECK(rational_value(sample.chart_derivatives[1]) == Q("1/64"));
        CHECK(rational_value(sample.chart_derivatives[2]) == Q("1/48"));
        CHECK(rational_value(sample.physical_gradient_per_metre[0]) == Q("1/24"));
        CHECK(rational_value(sample.physical_gradient_per_metre[1]) == derivatives / 2);
        CHECK(rational_value(sample.physical_gradient_per_metre[2]) == Q("1/128"));
        sides[side] = true;
    }
    REQUIRE(sides[0]);
    REQUIRE(sides[1]);

    const auto zero = f.r("0");
    const auto horizontal = take(partition.horizontal_faces());
    std::array<bool, 2> faces{};
    for (std::size_t i = 0; i < horizontal.size(); ++i) {
        const auto& face = horizontal[i];
        if (face.layer < 0 || face.layer > 1 || face.triangle != volumes.front().triangle)
            continue;
        const auto layer = static_cast<std::size_t>(face.layer);
        const auto expected = layer == 0 ? PgPartitionFaceSide::Upper : PgPartitionFaceSide::Lower;
        if (face.side != expected)
            continue;
        const auto sample =
            take(basis.evaluate(PgReferenceHorizontalPoint{static_cast<core::Index>(i), zero, zero},
                                column, reference_limits(), f.input_plan));
        CHECK(rational(sample.reference.zeta) == 0);
        CHECK(std::get<PgReferenceHorizontalMeasure>(sample.reference.measure).side == expected);
        CHECK(rational_value(sample.value) == 0);
        faces[layer] = true;
    }
    REQUIRE(faces[0]);
    REQUIRE(faces[1]);
    const auto side = take(basis.evaluate(PgReferenceLateralPoint{0, zero, depth}, column,
                                          reference_limits(), f.input_plan));
    CHECK(std::holds_alternative<PgReferenceLateralMeasure>(side.reference.measure));
    CHECK(rational_value(side.value) == 0); // Boundary vanishes for the interior C hat.
}

TEST_CASE("PG basis binding rejects stale rank or work limits before moving the owners",
          "[SH-D03][N4-08][negative][verification]") {
    Fixture f;
    const Frame frame{f};
    core::ResourcePlan plan{cap};
    auto inputs = prepare_inputs(f, plan);
    auto geometry = take(prepare_candidate(inputs, frame.view(), reference_limits(), plan));
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    auto rank = make_rank(geometry, workspace, plan);

    Fixture reordered;
    std::rotate(reordered.domain.begin(), reordered.domain.begin() + 1, reordered.domain.end());
    auto other_inputs = prepare_inputs(reordered, plan);
    auto different = take(prepare_candidate(other_inputs, frame.view(), reference_limits(), plan));
    core::ResourcePlan no_allocations{0};
    const auto stale = bind_basis(different, rank, rank_limits(), workspace, no_allocations);
    REQUIRE_FALSE(stale);
    CHECK(stale.error().code == core::FailureCode::StaleGeneration);
    REQUIRE(different.source());
    REQUIRE(rank.preparation());
    CHECK(no_allocations.allocated_bytes() == 0);

    auto limit = rank_limits();
    limit.maximum_steps = 0;
    const auto capped = bind_basis(geometry, rank, limit, workspace, no_allocations);
    REQUIRE_FALSE(capped);
    CHECK(capped.error().code == core::FailureCode::ResourceExceeded);
    REQUIRE(geometry.source());
    REQUIRE(rank.preparation());
    limit.cancellation = partition_test::requested_stop_token();
    const auto stopped = bind_basis(geometry, rank, limit, workspace, no_allocations);
    REQUIRE_FALSE(stopped);
    CHECK(stopped.error().code == core::FailureCode::Cancelled);
    REQUIRE(geometry.source());
    REQUIRE(rank.preparation());
    CHECK(no_allocations.allocated_bytes() == 0);
    const auto accepted = take(bind_basis(geometry, rank, rank_limits(), workspace, plan));
    REQUIRE(accepted.source());
    CHECK_FALSE(geometry.source());
    CHECK_FALSE(rank.preparation());
    const auto zero = f.r("0");
    const auto invalid = accepted.evaluate(PgReferenceVolumePoint{0, zero, zero, zero}, {-1},
                                           reference_limits(), no_allocations);
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == core::FailureCode::InvalidInput);
    CHECK(no_allocations.allocated_bytes() == 0);
}
