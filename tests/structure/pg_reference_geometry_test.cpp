#include "exact_polynomial_detail.hpp"
#include "pg_partition_test_support.hpp"

#include <array>
#include <optional>
#include <span>
#include <stop_token>
#include <string>
#include <string_view>
#include <utility>
#include <uvlas/structure/pg_reference_geometry.hpp>

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
    auto dictionary =
        take(PgCoefficientDictionary::prepare(*fixture.graph, coefficient_limits, workspace, plan));
    auto partition = take(PgIntegrationPartition::prepare(
        fixture.request(), partition_test::limits(), workspace, plan));
    auto graph = std::move(*fixture.graph);
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

} // namespace

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
    std::swap(inputs.graph, *f.graph);
    const auto stale = prepare_candidate(inputs, frame.view(), reference_limits(), plan);
    REQUIRE_FALSE(stale);
    CHECK(stale.error().code == core::FailureCode::StaleGeneration);
    std::swap(inputs.graph, *f.graph);
    CHECK(take(inputs.graph.view()).identity == source_liveness);
    CHECK(take(inputs.partition.identity()) == initial_identity);

    std::stop_source stop;
    REQUIRE(stop.request_stop());
    auto limits = reference_limits();
    limits.cancellation = stop.get_token();
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
    std::swap(inputs.graph, *f.graph);
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    limits.cancellation = stop.get_token();
    limits.maximum_partition_visits = 0;
    check_code(prepare_candidate(inputs, frame.view(), limits, call_plan),
               core::FailureCode::Cancelled);
    check_code(PgReferenceGeometry::diagnose_point(inputs, frame.view(), query, limits, call_plan),
               core::FailureCode::Cancelled);
    CHECK(call_plan.allocated_bytes() == 0);
    std::swap(inputs.graph, *f.graph);

    const auto pieces = take(inputs.partition.pieces()).size();
    const auto vertices = take(inputs.partition.piece_vertex_indices()).size();
    frame.length = f.r("2");
    limits = reference_limits();
    limits.maximum_partition_visits = comparison + pieces + vertices;
    const auto prepared = take(prepare_candidate(inputs, frame.view(), limits, call_plan));
    CHECK(prepared.preparation().partition_visits == limits.maximum_partition_visits);
}
