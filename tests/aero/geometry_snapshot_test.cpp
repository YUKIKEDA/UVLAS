#include "aero/geometry_snapshot_storage.hpp"
#include "geometry_snapshot_fixture.hpp"

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <exception>
#include <limits>
#include <stop_token>
#include <thread>
#include <type_traits>

namespace {
using namespace uvlas;
using namespace uvlas::test_geometry;
using namespace subject;
using Rat = boost::multiprecision::cpp_rational;
using Int = boost::multiprecision::cpp_int;
Rat exact(ScalarWords words) {
    if (const auto* p = std::get_if<SingleWord>(&words))
        return Rat{p->value};
    const auto* p = std::get_if<DoubleWord>(&words);
    if (p == nullptr)
        throw std::runtime_error("Invalid independent scalar");
    const Rat value = Rat{p->high} + Rat{p->low};
    return boost::multiprecision::numerator(value) == 0 ? Rat{0} : value;
}
Rat distance(const ExactTriangleResult& result) {
    const auto* outcome = result.outcome();
    if (outcome == nullptr)
        throw std::runtime_error("Missing owned predicate");
    const auto* value = std::get_if<NondegenerateTriangle>(outcome);
    if (value == nullptr)
        throw std::runtime_error("Unexpected degeneracy");
    const Rat n{Int{value->distance_squared.numerator}};
    if (boost::multiprecision::numerator(n) == 0)
        return Rat{0};
    CHECK(value->distance_squared.unit_exponent == -4296);
    return n / (Rat{Int{value->distance_squared.denominator}} * Rat{Int{1} << 4296});
}
CanonicalSupport support(const OwnedGeometry& g, QueryDefinition d, GeometryWorkspace& w) {
    core::ResourcePlan plan{budget};
    auto result = take(resolve_query(g, {QueryId{0}, d}, g.context(), limits(), w, plan));
    CHECK(plan.allocated_bytes() == ResolvedQuery::payload_bytes());
    return take(result.support());
}
} // namespace

TEST_CASE("Owned geometry resolves every stable identity and keeps all P1 triangle basis",
          "[UI-D01][UI-D05][N4-02][N4-03][verification]") {
    core::ResourcePlan preparation{budget}, output{budget};
    auto w = workspace(preparation);
    TriangleFixture fixture;
    const auto input = fixture.input();
    auto geometry = take(prepare_owned_geometry(input, limits(), *w, output));
    CHECK(output.allocated_bytes() == take(geometry.payload_bytes()));
    CHECK(take(geometry.owner(PointId{11})) == PointId{11});
    CHECK(take(geometry.facet_points(FacetId{70})) == fixture.points);
    CHECK(take(geometry.facet_edges(FacetId{70})) ==
          std::array{EdgeId{10}, EdgeId{20}, EdgeId{30}});
    CHECK(take(geometry.incidence(FacetId{70})) == std::array<std::int8_t, 3>{-1, 1, 1});
    CHECK(take(geometry.facet_face(FacetId{70})) == PhysicalFaceId{50});
    CHECK(take(geometry.face_band(PhysicalFaceId{50})) == BandId{100});
    const auto map = take(geometry.scatter(FacetId{70}));
    const auto* independent = std::get_if<IndependentP1Scatter>(&map);
    REQUIRE(independent != nullptr);
    CHECK(independent->coefficients == std::array<std::int8_t, 9>{1, 0, 0, 0, 1, 0, 0, 0, 1});
    CHECK(exact(take(geometry.core_radius())) == Rat{1} / 4);
    const auto old_identity = take(geometry.identity());
    fixture.current.fill(std::numeric_limits<double>::quiet_NaN());
    fixture.points[0] = PointId{999};
    CHECK(take(geometry.identity()) == old_identity);
    CHECK(exact(take(geometry.point_words(PointId{11}, Configuration::Current))[0]) == 0);
    CHECK(take(geometry.facet_points(FacetId{70}))[0] == PointId{11});
    failure(geometry.owner(PointId{999}), core::FailureCode::InvalidInput);
}

TEST_CASE("Every query variant retains its origin while normalizing exact support",
          "[UI-D03][UI-D05][N4-03][N4-06][verification][negative]") {
    core::ResourcePlan preparation{budget}, output{budget};
    auto w = workspace(preparation);
    const TriangleFixture fixture;
    auto g = take(prepare_owned_geometry(fixture.input(), limits(), *w, output));
    CHECK(std::holds_alternative<NoMaterialSupport>(
        support(g, FreePoint{FreeSingle{{.75, 1, 0}}}, *w)));
    CHECK(std::get<PointId>(support(g, SourceVertex{PointId{11}}, *w)) == PointId{11});
    // Edge10's canonical direction is Point3 -> Point11, opposite the raw input.
    CHECK(std::get<PointId>(support(g, SourceEdgePoint{EdgeId{10}, SingleWord{0}}, *w)) ==
          PointId{3});
    CHECK(std::get<PointId>(support(g, SourceEdgePoint{EdgeId{10}, SingleWord{1}}, *w)) ==
          PointId{11});
    CHECK(std::get<EdgeId>(support(g, SourceEdgePoint{EdgeId{10}, SingleWord{.25}}, *w)) ==
          EdgeId{10});
    const std::array parameters{std::array{0., 0.},  std::array{1., 0.},   std::array{0., 1.},
                                std::array{.25, 0.}, std::array{.25, .75}, std::array{0., .25},
                                std::array{.25, .25}};
    const std::array<CanonicalSupport, 7> expected{PointId{11}, PointId{3}, PointId{7}, EdgeId{10},
                                                   EdgeId{20},  EdgeId{30}, FacetId{70}};
    for (std::size_t i = 0; i < parameters.size(); ++i) {
        core::ResourcePlan query_plan{budget}, result_plan{budget};
        const QueryRequest request{QueryId{std::numeric_limits<std::uint64_t>::max()},
                                   SourceFacePoint{FacetId{70}, SingleWord{parameters[i][0]},
                                                   SingleWord{parameters[i][1]}}};
        auto query = take(resolve_query(g, request, g.context(), limits(), *w, query_plan));
        CHECK(take(query.support()) == expected[i]);
        const auto retained = take(query.definition());
        CHECK(retained.id == request.id);
        CHECK(std::get<SourceFacePoint>(retained.definition).facet == FacetId{70});
        const auto mathematical =
            take(classify_resolved_query(g, query, FacetId{70}, limits(), *w, result_plan));
        CHECK(distance(mathematical) == 0);
    }
    CHECK(std::get<EdgeId>(support(g, SourceEdgePoint{EdgeId{20}, DoubleWord{-.25, .5}}, *w)) ==
          EdgeId{20});
    CHECK(std::get<PointId>(support(g, SourceEdgePoint{EdgeId{20}, DoubleWord{.25, -.25}}, *w)) ==
          PointId{3});
    core::ResourcePlan rejected_plan{budget};
    failure(resolve_query(g, {QueryId{1}, SourceEdgePoint{EdgeId{10}, DoubleWord{.25, -.5}}},
                          g.context(), limits(), *w, rejected_plan),
            core::FailureCode::InvalidInput);
    failure(resolve_query(
                g, {QueryId{1}, SourceFacePoint{FacetId{70}, SingleWord{.75}, SingleWord{.5}}},
                g.context(), limits(), *w, rejected_plan),
            core::FailureCode::InvalidInput);
    failure(resolve_query(g, {QueryId{1}, SourceVertex{PointId{999}}}, g.context(), limits(), *w,
                          rejected_plan),
            core::FailureCode::InvalidInput);
}

TEST_CASE("Exact alias equality is a physical word sum rather than raw representation",
          "[UI-D05][UW-D07][N4-03][N4-06][negative][verification]") {
    TriangleFixture fixture;
    fixture.two_word = true;
    const std::array ids{PointId{11}, PointId{3}, PointId{7}, PointId{0}};
    const std::array bindings{
        binding(), binding(), binding(),
        PointBinding{ReferenceMapId{1}, MappingGeneration{1}, AliasOf{PointId{11}}}};
    std::array<double, 12> high{0, 0, 0, 3, 0, 0, 0, 4, 0, 1, -0., .5};
    std::array<double, 12> low{0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0., -.5};
    auto input = fixture.input();
    input.points = ids;
    input.bindings = bindings;
    input.reference = DoublePositions{high, low};
    input.current = DoublePositions{high, low};
    core::ResourcePlan prep{budget}, out{budget};
    auto w = workspace(prep);
    auto geometry = take(prepare_owned_geometry(input, limits(), *w, out));
    CHECK(take(geometry.owner(PointId{0})) == PointId{11});
    const auto alias = take(geometry.point_words(PointId{0}, Configuration::Current));
    CHECK(exact(alias[0]) == 0);
    CHECK(std::get<DoubleWord>(alias[0]).high == 1);
    const auto old = take(geometry.identity());
    low[9] = -.5;
    core::ResourcePlan denied{budget};
    failure(prepare_owned_geometry(input, limits(), *w, denied),
            core::FailureCode::GeometryInvalid);
    CHECK(take(geometry.identity()) == old);
    low[9] = -1;
    auto bad_bindings = bindings;
    bad_bindings[3] = PointBinding{ReferenceMapId{1}, MappingGeneration{2}, AliasOf{PointId{11}}};
    input.bindings = bad_bindings;
    failure(prepare_owned_geometry(input, limits(), *w, denied),
            core::FailureCode::StaleGeneration);
    bad_bindings[3] = PointBinding{ReferenceMapId{1}, MappingGeneration{1}, AliasOf{PointId{0}}};
    failure(prepare_owned_geometry(input, limits(), *w, denied), core::FailureCode::InvalidInput);
    bad_bindings[3] = PointBinding{ReferenceMapId{1}, MappingGeneration{1}, AliasOf{PointId{999}}};
    failure(prepare_owned_geometry(input, limits(), *w, denied), core::FailureCode::InvalidInput);
}

TEST_CASE("Quad sources retain fixed children diagonal incidence and both LR maps",
          "[UI-D01][UI-D05][N4-02][verification][negative]") {
    const std::array p{PointId{101}, PointId{7}, PointId{55}, PointId{3}};
    const std::array b{binding(), binding(), binding(), binding()};
    const std::array<double, 12> xyz{0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0};
    const std::array edges{EdgeInput{EdgeId{1}, {p[0], p[1]}}, EdgeInput{EdgeId{2}, {p[1], p[2]}},
                           EdgeInput{EdgeId{3}, {p[2], p[3]}}, EdgeInput{EdgeId{4}, {p[3], p[0]}},
                           EdgeInput{EdgeId{5}, {p[0], p[2]}}};
    const std::array faces{FaceInput{PhysicalFaceId{8}, BandId{0}}};
    const std::array bands{BandId{0}};
    const std::array quads{
        QuadInput{faces[0].id,
                  p,
                  {EdgeId{1}, EdgeId{2}, EdgeId{3}, EdgeId{4}},
                  EdgeId{5},
                  {FacetId{0}, FacetId{std::numeric_limits<std::uint64_t>::max()}}}};
    const GeometryInput input{p,
                              b,
                              SinglePositions{xyz},
                              SinglePositions{xyz},
                              edges,
                              faces,
                              bands,
                              {},
                              quads,
                              SingleWord{0},
                              context()};
    core::ResourcePlan prep{budget}, out{budget};
    auto w = workspace(prep);
    auto geometry = take(prepare_owned_geometry(input, limits(), *w, out));
    const auto first = FacetId{0}, second = quads[0].facets[1];
    CHECK(take(geometry.facet_points(first)) == std::array{p[0], p[2], p[1]});
    CHECK(take(geometry.facet_points(second)) == std::array{p[0], p[3], p[2]});
    CHECK(take(geometry.facet_edges(first)) == std::array{EdgeId{5}, EdgeId{2}, EdgeId{1}});
    CHECK(take(geometry.facet_edges(second)) == std::array{EdgeId{4}, EdgeId{3}, EdgeId{5}});
    CHECK(take(geometry.incidence(first))[0] == -take(geometry.incidence(second))[2]);
    CHECK(std::get<QuadLRScatter>(take(geometry.scatter(first))).coefficients ==
          std::array<std::int8_t, 6>{1, 0, 0, 1, 1, 0});
    CHECK(std::get<QuadLRScatter>(take(geometry.scatter(second))).coefficients ==
          std::array<std::int8_t, 6>{1, 0, 0, 1, 0, 1});
    // A known signed incidence at every vertex independently has zero boundary.
    for (const auto id : {first, second}) {
        const auto vertices = take(geometry.facet_points(id));
        const auto es = take(geometry.facet_edges(id));
        const auto signs = take(geometry.incidence(id));
        std::array<int, 4> sum{};
        for (std::size_t edge = 0; edge < 3; ++edge) {
            const auto& original = edges[static_cast<std::size_t>(es[edge].value() - 1)].points;
            const auto lo = std::min(original[0], original[1]),
                       hi = std::max(original[0], original[1]);
            for (std::size_t point = 0; point < 4; ++point) {
                if (p[point] == lo)
                    sum[point] -= signs[edge];
                if (p[point] == hi)
                    sum[point] += signs[edge];
            }
        }
        CHECK(std::all_of(sum.begin(), sum.end(), [](int x) { return x == 0; }));
        CHECK(vertices[0] == p[0]);
    }
}
TEST_CASE("Material query products reach exact predicates without a two word projection",
          "[UI-D03][UI-D05][N4-06][verification][negative]") {
    const double tiny = std::numeric_limits<double>::denorm_min();
    const std::array p{PointId{1}, PointId{2}, PointId{3}, PointId{4}, PointId{5}, PointId{6}};
    const std::array bindings{binding(), binding(), binding(), binding(), binding(), binding()};
    const std::array<double, 18> xyz{1, 0, 0, 2, tiny, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1};
    const std::array edges{
        EdgeInput{EdgeId{10}, {p[0], p[1]}}, EdgeInput{EdgeId{20}, {p[1], p[2]}},
        EdgeInput{EdgeId{30}, {p[2], p[0]}}, EdgeInput{EdgeId{40}, {p[3], p[4]}},
        EdgeInput{EdgeId{50}, {p[4], p[5]}}, EdgeInput{EdgeId{60}, {p[5], p[3]}}};
    const std::array faces{FaceInput{PhysicalFaceId{50}, BandId{100}},
                           FaceInput{PhysicalFaceId{60}, BandId{100}}};
    const std::array bands{BandId{100}};
    const std::array triangles{
        TriangleInput{
            FacetId{70}, faces[0].id, {p[0], p[1], p[2]}, {edges[0].id, edges[1].id, edges[2].id}},
        TriangleInput{
            FacetId{80}, faces[1].id, {p[3], p[4], p[5]}, {edges[3].id, edges[4].id, edges[5].id}}};
    const GeometryInput input{p,
                              bindings,
                              SinglePositions{xyz},
                              SinglePositions{xyz},
                              edges,
                              faces,
                              bands,
                              triangles,
                              {},
                              SingleWord{0},
                              context()};
    core::ResourcePlan prep{budget}, out{budget};
    auto w = workspace(prep);
    auto g = take(prepare_owned_geometry(input, limits(), *w, out));
    const DoubleWord t{.5, tiny};
    core::ResourcePlan qp{budget}, rp{budget};
    auto q = take(resolve_query(g, {QueryId{1}, SourceEdgePoint{EdgeId{10}, t}}, g.context(),
                                limits(), *w, qp));
    const auto original = take(classify_resolved_query(g, q, FacetId{70}, limits(), *w, rp));
    CHECK(distance(original) == 0);
    const auto separate = take(classify_resolved_query(g, q, FacetId{80}, limits(), *w, rp));
    const Rat x = Rat{1} + exact(t);
    CHECK(distance(separate) == x * x);
    // The exact y contains 2^-2148, below the smallest possible binary64 word.
    const Rat y = exact(t) * Rat{tiny};
    CHECK(y == Rat{tiny} / 2 + Rat{tiny} * Rat{tiny});
    CHECK(y != Rat{tiny} / 2);
    const DoubleWord u{.25, 0x1p-1000}, v{.25, tiny};
    auto face = take(resolve_query(g, {QueryId{2}, SourceFacePoint{FacetId{70}, u, v}}, g.context(),
                                   limits(), *w, qp));
    const auto face_original =
        take(classify_resolved_query(g, face, FacetId{70}, limits(), *w, rp));
    CHECK(distance(face_original) == 0);
    const auto face_separate =
        take(classify_resolved_query(g, face, FacetId{80}, limits(), *w, rp));
    const Rat face_x = Rat{1} + exact(u);
    CHECK(distance(face_separate) == face_x * face_x);
}

TEST_CASE("Geometry state matching excludes reporting additions but retains true identity",
          "[UI-D05][N4-02][N4-03][verification][negative]") {
    TriangleFixture fixture;
    core::ResourcePlan prep{budget}, out{budget};
    auto w = workspace(prep);
    auto g = take(prepare_owned_geometry(fixture.input(), limits(), *w, out));
    const QueryRequest request{QueryId{1}, SourceVertex{PointId{11}}};
    core::ResourcePlan qp{budget}, rp{budget};
    auto q = take(resolve_query(g, request, g.context(), limits(), *w, qp));
    auto reporting = g.context();
    reporting.diagnostics.component = core::ComponentId{99};
    reporting.diagnostics.problem = core::ProblemId{core::PanelId{98}};
    auto allowed = resolve_query(g, request, reporting, limits(), *w, qp);
    REQUIRE(allowed.has_value());
    fixture.stamp = reporting;
    auto equivalent = take(prepare_owned_geometry(fixture.input(), limits(), *w, out));
    CHECK(take(g.identity()) == take(equivalent.identity()));
    REQUIRE(classify_resolved_query(equivalent, q, FacetId{70}, limits(), *w, rp).has_value());
    for (unsigned changed = 0; changed < 14; ++changed) {
        auto stale = g.context();
        switch (changed) {
        case 0:
            stale.model = core::ModelGeneration{100};
            break;
        case 1:
            stale.reference = ReferenceGeometryGeneration{100};
            break;
        case 2:
            stale.current = core::GeometryGeneration{100};
            break;
        case 3:
            stale.mapping = MappingGeneration{100};
            break;
        case 4:
            stale.connectivity = ConnectivityGeneration{100};
            break;
        case 5:
            stale.policy = KernelPolicyId{100};
            break;
        case 6:
            stale.kernel = KernelGeneration{100};
            break;
        case 7:
            stale.radius = RadiusGeneration{100};
            break;
        case 8:
            stale.precision = core::PrecisionGeneration{100};
            break;
        case 9:
            stale.numerical = core::NumericalGeneration{100};
            break;
        case 10:
            stale.time = 2;
            break;
        case 11:
            stale.diagnostics.parent = core::AcceptedStateId{100};
            break;
        case 12:
            stale.diagnostics.trial = core::TrialId{100};
            break;
        default:
            stale.diagnostics.last_accepted = core::AcceptedStateId{100};
            break;
        }
        failure(resolve_query(g, request, stale, limits(), *w, qp),
                core::FailureCode::StaleGeneration);
    }
    auto invalid = g.context();
    invalid.time = std::numeric_limits<double>::quiet_NaN();
    failure(resolve_query(g, request, invalid, limits(), *w, qp), core::FailureCode::Nonfinite);
    fixture.current[3] = 6;
    fixture.stamp = g.context();
    auto changed = take(prepare_owned_geometry(fixture.input(), limits(), *w, out));
    CHECK(take(changed.identity()) != take(g.identity()));
    failure(classify_resolved_query(changed, q, FacetId{70}, limits(), *w, rp),
            core::FailureCode::StaleGeneration);
    auto caps = limits();
    caps.counts.points = 2;
    failure(resolve_query(g, request, g.context(), caps, *w, qp),
            core::FailureCode::ResourceExceeded);
    caps = limits();
    caps.predicate.external_memory = numerics::linear::RequireProvenBound{28672};
    failure(resolve_query(g, request, g.context(), caps, *w, qp), core::FailureCode::Unsupported);
    CHECK(distance(take(classify_resolved_query(g, q, FacetId{70}, limits(), *w, rp))) == 0);
}
TEST_CASE("Geometry limits cancellation and workspace ownership are explicit",
          "[UI-D05][UI-D06][N4-02][N4-03][N4-09][negative][verification]") {
    const TriangleFixture fixture;
    core::ResourcePlan prep{budget}, out{budget};
    auto w = workspace(prep);
    auto g = take(prepare_owned_geometry(fixture.input(), limits(), *w, out));
    core::ResourcePlan qp{budget};
    const QueryRequest request{QueryId{1}, SourceVertex{PointId{11}}};
    const auto held = w->acquire();
    REQUIRE(held.has_value());
    failure(resolve_query(g, request, g.context(), limits(), *w, qp),
            core::FailureCode::InvalidInput);
    CHECK(qp.allocated_bytes() == 0);
    // A separate workspace owns independent scratch and can operate during that lease.
    core::ResourcePlan other_preparation{budget};
    auto other = workspace(other_preparation);
    auto q = take(resolve_query(g, request, g.context(), limits(), *other, qp));
    core::ResourcePlan rp{budget};
    CHECK(distance(take(classify_resolved_query(g, q, FacetId{70}, limits(), *other, rp))) == 0);
    // MSVC14.51 request_stop is nonconst; libstdc++ declares it const.
    // NOLINTNEXTLINE(misc-const-correctness)
    std::stop_source stop;
    stop.request_stop();
    failure(prepare_owned_geometry(fixture.input(), limits(), *other, out, stop.get_token()),
            core::FailureCode::Cancelled);
    auto limited = limits();
    limited.max_steps = 0;
    failure(resolve_query(g, request, g.context(), limited, *other, qp),
            core::FailureCode::ResourceExceeded);
    core::ResourcePlan exact_query{ResolvedQuery::payload_bytes()};
    REQUIRE(resolve_query(g, request, g.context(), limits(), *other, exact_query).has_value());
    core::ResourcePlan short_query{ResolvedQuery::payload_bytes() - 1};
    failure(resolve_query(g, request, g.context(), limits(), *other, short_query),
            core::FailureCode::ResourceExceeded);
    CHECK(short_query.allocated_bytes() == 0);
    auto moved = std::move(g);
    // The API defines an inactive moved-from owner.
    // NOLINTNEXTLINE(bugprone-use-after-move)
    CHECK_FALSE(g.active());
    failure(g.identity(), core::FailureCode::InvalidInput);
    CHECK(distance(take(classify_resolved_query(moved, q, FacetId{70}, limits(), *other, rp))) ==
          0);
    auto moved_query = std::move(q);
    // The API defines an inactive moved-from query.
    // NOLINTNEXTLINE(bugprone-use-after-move)
    CHECK_FALSE(q.active());
    failure(q.support(), core::FailureCode::InvalidInput);
    CHECK(distance(take(classify_resolved_query(moved, moved_query, FacetId{70}, limits(), *other,
                                                rp))) == 0);
}

TEST_CASE("Geometry rejects malformed words connectivity and both configuration degeneracies",
          "[UI-D01][UI-D05][UI-D06][N4-06][negative][verification]") {
    core::ResourcePlan prep{budget}, out{budget};
    auto w = workspace(prep);
    for (unsigned domain = 0; domain < 4; ++domain)
        for (std::size_t word = 0; word < 9; ++word) {
            TriangleFixture f;
            f.two_word = true;
            switch (domain) {
            case 0:
                f.reference[word] = std::numeric_limits<double>::infinity();
                break;
            case 1:
                f.reference_low[word] = std::numeric_limits<double>::quiet_NaN();
                break;
            case 2:
                f.current[word] = std::numeric_limits<double>::infinity();
                break;
            default:
                f.current_low[word] = std::numeric_limits<double>::quiet_NaN();
                break;
            }
            failure(prepare_owned_geometry(f.input(), limits(), *w, out),
                    core::FailureCode::Nonfinite);
        }
    TriangleFixture f;
    auto input = f.input();
    input.current = DoublePositions{f.current, {}};
    failure(prepare_owned_geometry(input, limits(), *w, out), core::FailureCode::InvalidInput);
    f.points[1] = f.points[0];
    failure(prepare_owned_geometry(f.input(), limits(), *w, out), core::FailureCode::InvalidInput);
    f = TriangleFixture{};
    f.edges[0] = EdgeInput{EdgeId{10}, {PointId{11}, PointId{999}}};
    failure(prepare_owned_geometry(f.input(), limits(), *w, out), core::FailureCode::InvalidInput);
    f = TriangleFixture{};
    std::swap(f.triangles[0].edges[0], f.triangles[0].edges[1]);
    failure(prepare_owned_geometry(f.input(), limits(), *w, out), core::FailureCode::InvalidInput);
    for (const bool reference : {false, true}) {
        f = TriangleFixture{};
        auto& xyz = reference ? f.reference : f.current;
        xyz[7] = 0;
        failure(prepare_owned_geometry(f.input(), limits(), *w, out),
                core::FailureCode::GeometryInvalid);
    }
    for (const ScalarWords radius :
         std::array<ScalarWords, 3>{SingleWord{-1}, DoubleWord{1, -2}, DoubleWord{0, -0x1p-1074}}) {
        f = TriangleFixture{};
        f.radius = radius;
        failure(prepare_owned_geometry(f.input(), limits(), *w, out),
                core::FailureCode::InvalidInput);
    }
    f = TriangleFixture{};
    f.radius = DoubleWord{-.25, .5};
    REQUIRE(prepare_owned_geometry(f.input(), limits(), *w, out).has_value());
}
TEST_CASE(
    "Immutable geometry and queries survive source plans and concurrent independent workspaces",
    "[UI-D05][N4-02][N4-03][verification]") {
    auto retained = [] {
        const TriangleFixture f;
        core::ResourcePlan prep{budget}, output{budget}, query_plan{budget};
        auto w = workspace(prep);
        auto g = take(prepare_owned_geometry(f.input(), limits(), *w, output));
        auto q = take(resolve_query(g, {QueryId{9}, FreePoint{FreeSingle{{.75, 1, 20}}}},
                                    g.context(), limits(), *w, query_plan));
        return std::pair{std::move(g), std::move(q)};
    }();
    auto& g = retained.first;
    auto& q = retained.second;
    const auto identity = take(g.identity());
    // Each thread owns its Plan/workspace/result; immutable owners alone are shared.
    std::array<bool, 2> passed{};
    std::array<std::exception_ptr, 2> errors{};
    std::array<std::thread, 2> workers;
    for (std::size_t i = 0; i < workers.size(); ++i)
        workers[i] = std::thread([&, i] {
            try {
                core::ResourcePlan prep{budget}, result_plan{budget};
                auto w = GeometryWorkspace::create(prep);
                if (!w)
                    return;
                auto r = classify_resolved_query(g, q, FacetId{70}, limits(), **w, result_plan);
                if (!r)
                    return;
                const auto* outcome = r->outcome();
                if (outcome == nullptr)
                    return;
                const auto* triangle = std::get_if<NondegenerateTriangle>(outcome);
                if (triangle == nullptr)
                    return;
                // Exact 20^2, compared without Catch machinery in the worker.
                const Int expected =
                    Int{triangle->distance_squared.denominator} * (Int{400} << 4296);
                passed[i] = triangle->distance_squared.numerator == expected;
            } catch (...) {
                errors[i] = std::current_exception();
            }
        });
    for (auto& worker : workers)
        worker.join();
    for (std::size_t i = 0; i < workers.size(); ++i) {
        CHECK(errors[i] == nullptr);
        CHECK(passed[i]);
    }
    CHECK(take(g.identity()) == identity);
}
TEST_CASE("All six triangle orderings preserve original material interpolation",
          "[UI-D01][UI-D03][UI-D05][N4-06][verification]") {
    std::array<std::size_t, 3> permutation{0, 1, 2};
    do {
        TriangleFixture f;
        const auto original = f.points;
        for (std::size_t i = 0; i < 3; ++i)
            f.triangles[0].points[i] = original[permutation[i]];
        for (std::size_t i = 0; i < 3; ++i) {
            const auto a = f.triangles[0].points[i], b = f.triangles[0].points[(i + 1) % 3];
            for (const auto& edge : f.edges)
                if ((edge.points[0] == a && edge.points[1] == b) ||
                    (edge.points[0] == b && edge.points[1] == a))
                    f.triangles[0].edges[i] = edge.id;
        }
        core::ResourcePlan prep{budget}, out{budget}, qp{budget}, rp{budget};
        auto w = workspace(prep);
        auto g = take(prepare_owned_geometry(f.input(), limits(), *w, out));
        const std::array<double, 3> weights{.25, .5, .25};
        auto q = take(resolve_query(
            g,
            {QueryId{0}, SourceFacePoint{FacetId{70}, SingleWord{weights[permutation[1]]},
                                         SingleWord{weights[permutation[2]]}}},
            g.context(), limits(), *w, qp));
        const auto result = take(classify_resolved_query(g, q, FacetId{70}, limits(), *w, rp));
        CHECK(distance(result) == 0);
        CHECK(std::get<FacetId>(take(q.support())) == FacetId{70});
    } while (std::next_permutation(permutation.begin(), permutation.end()));
}
TEST_CASE("Empty topology and coincident canonical point identities retain their exact domains",
          "[UI-D05][N4-02][N4-03][negative][verification]") {
    const auto stamp = context();
    core::ResourcePlan prep{budget}, out{budget};
    auto w = workspace(prep);
    const GeometryInput empty{
        {}, {}, SinglePositions{{}}, SinglePositions{{}}, {}, {}, {}, {}, {}, SingleWord{0}, stamp};
    auto zero = take(prepare_owned_geometry(empty, limits(), *w, out));
    CHECK(take(zero.counts()).points == 0);
    core::ResourcePlan qp{budget};
    failure(resolve_query(zero, {QueryId{0}, SourceVertex{PointId{0}}}, stamp, limits(), *w, qp),
            core::FailureCode::InvalidInput);
    const std::array p{PointId{0}, PointId{std::numeric_limits<std::uint64_t>::max()}};
    const std::array b{binding(), binding()};
    const std::array<double, 6> xyz{1, 0, 0, 1, 0, 0};
    const GeometryInput input{
        p, b, SinglePositions{xyz}, SinglePositions{xyz}, {}, {}, {}, {}, {}, SingleWord{0}, stamp};
    auto g = take(prepare_owned_geometry(input, limits(), *w, out));
    CHECK(take(g.owner(p[0])) == p[0]);
    CHECK(take(g.owner(p[1])) == p[1]);
    CHECK(take(g.owner(p[0])) != take(g.owner(p[1])));
    for (const auto point : p) {
        auto q = take(resolve_query(g, {QueryId{point.value()}, SourceVertex{point}}, stamp,
                                    limits(), *w, qp));
        CHECK(std::get<PointId>(take(q.support())) == point);
        CHECK(take(q.definition()).id.value() == point.value());
    }
}
TEST_CASE("Material word products retain extreme dyadic scales and either orientation",
          "[UI-D01][UI-D03][UI-D05][N4-06][verification]") {
    // Boundary exponents exercise subnormal decoding, the normal transition,
    // large exponent shifts while retaining both finite coordinate words.
    for (const int exponent : std::array{-1074, -1073, -1022, -1000, -1, 0, 500, 1000, 1023})
        for (const int sign : {-1, 1}) {
            const double h = static_cast<double>(sign) * std::ldexp(1., exponent);
            const double tail =
                exponent == -1074 ? 0. : static_cast<double>(sign) * std::ldexp(1., exponent - 1);
            const std::array p{PointId{1}, PointId{2}, PointId{3},
                               PointId{4}, PointId{5}, PointId{6}};
            const std::array bindings{binding(), binding(), binding(),
                                      binding(), binding(), binding()};
            const std::array<double, 18> high{0, 0, 0, h, 0, 0, 0, h, 0, 0, 0, 0, 0, h, 0, 0, 0, h};
            const std::array<double, 18> low{0, 0, 0, tail, 0,    0, 0, tail, 0,
                                             0, 0, 0, 0,    tail, 0, 0, 0,    tail};
            const std::array edges{
                EdgeInput{EdgeId{10}, {p[0], p[1]}}, EdgeInput{EdgeId{20}, {p[1], p[2]}},
                EdgeInput{EdgeId{30}, {p[2], p[0]}}, EdgeInput{EdgeId{40}, {p[3], p[4]}},
                EdgeInput{EdgeId{50}, {p[4], p[5]}}, EdgeInput{EdgeId{60}, {p[5], p[3]}}};
            const std::array faces{FaceInput{PhysicalFaceId{50}, BandId{100}},
                                   FaceInput{PhysicalFaceId{60}, BandId{100}}};
            const std::array bands{BandId{100}};
            const std::array triangles{TriangleInput{FacetId{70},
                                                     faces[0].id,
                                                     {p[0], p[1], p[2]},
                                                     {edges[0].id, edges[1].id, edges[2].id}},
                                       TriangleInput{FacetId{80},
                                                     faces[1].id,
                                                     {p[3], p[4], p[5]},
                                                     {edges[3].id, edges[4].id, edges[5].id}}};
            const GeometryInput input{p,
                                      bindings,
                                      DoublePositions{high, low},
                                      DoublePositions{high, low},
                                      edges,
                                      faces,
                                      bands,
                                      triangles,
                                      {},
                                      SingleWord{0},
                                      context()};
            core::ResourcePlan prep{budget}, out{budget}, qp{budget}, rp{budget};
            auto w = workspace(prep);
            auto g = take(prepare_owned_geometry(input, limits(), *w, out));
            const DoubleWord u{.5, std::numeric_limits<double>::denorm_min()};
            auto q = take(
                resolve_query(g, {QueryId{0}, SourceFacePoint{FacetId{70}, u, SingleWord{.25}}},
                              g.context(), limits(), *w, qp));
            CHECK(distance(take(classify_resolved_query(g, q, FacetId{70}, limits(), *w, rp))) ==
                  0);
            const Rat coordinate = exact(u) * (Rat{h} + Rat{tail});
            CHECK(distance(take(classify_resolved_query(g, q, FacetId{80}, limits(), *w, rp))) ==
                  coordinate * coordinate);
        }
}
TEST_CASE("Free queries construct the selected word alternative and preserve every original word",
          "[UI-D03][UI-D05][N4-03][N4-06][verification]") {
    const TriangleFixture fixture;
    core::ResourcePlan prep{budget}, out{budget};
    auto w = workspace(prep);
    const auto g = take(prepare_owned_geometry(fixture.input(), limits(), *w, out));
    const std::array<PointWords, 2> words{FreeSingle{{.75, 1, 20}},
                                          FreeDouble{{.75, 1, 20}, {0x1p-80, -0x1p-90, 0x1p-100}}};
    for (const auto& source : words) {
        core::ResourcePlan qp{budget}, rp{budget};
        const auto q =
            take(resolve_query(g, {QueryId{17}, FreePoint{source}}, g.context(), limits(), *w, qp));
        CHECK(std::holds_alternative<NoMaterialSupport>(take(q.support())));
        const auto retained = take(q.definition());
        const auto* point = std::get_if<FreePoint>(&retained.definition);
        REQUIRE(point != nullptr);
        CHECK(point->words.index() == source.index());
        Rat height{20};
        if (const auto* single = std::get_if<FreeSingle>(&source)) {
            const auto* saved = std::get_if<FreeSingle>(&point->words);
            REQUIRE(saved != nullptr);
            CHECK(saved->words == single->words);
        } else {
            const auto* pair = std::get_if<FreeDouble>(&source);
            REQUIRE(pair != nullptr);
            const auto* saved = std::get_if<FreeDouble>(&point->words);
            REQUIRE(saved != nullptr);
            CHECK(saved->high == pair->high);
            CHECK(saved->low == pair->low);
            height += Rat{pair->low[2]};
        }
        CHECK(distance(take(classify_resolved_query(g, q, FacetId{70}, limits(), *w, rp))) ==
              height * height);
    }
}
TEST_CASE("Connection diagnostics select the lowest facet ID independently of source row order",
          "[UI-D05][UI-D06][N4-03][negative][verification]") {
    const TriangleFixture fixture;
    std::array rows{
        TriangleInput{
            FacetId{90}, PhysicalFaceId{999}, fixture.points, {EdgeId{10}, EdgeId{20}, EdgeId{30}}},
        TriangleInput{FacetId{10},
                      PhysicalFaceId{999},
                      fixture.points,
                      {EdgeId{10}, EdgeId{20}, EdgeId{30}}}};
    core::ResourcePlan prep{budget};
    auto w = workspace(prep);
    for (unsigned order = 0; order < 2; ++order) {
        auto input = fixture.input();
        input.triangles = rows;
        core::ResourcePlan plan{budget};
        const auto rejected = prepare_owned_geometry(input, limits(), *w, plan);
        failure(rejected, core::FailureCode::InvalidInput);
        if (rejected || !rejected.error().problem) {
            FAIL_CHECK("Missing typed facet failure");
            return;
        }
        const auto* problem = std::get_if<FacetId>(&*rejected.error().problem);
        if (problem == nullptr) {
            FAIL_CHECK("Wrong typed geometry problem");
            return;
        }
        CHECK(*problem == FacetId{10});
        std::swap(rows[0], rows[1]);
    }
}
namespace {
void expected_problem(const AeroResult<OwnedGeometry>& result, core::FailureCode code,
                      GeometryProblem wanted) {
    failure(result, code);
    if (result || !result.error().problem) {
        FAIL_CHECK("Missing typed geometry failure");
        return;
    }
    CHECK(*result.error().problem == wanted);
}
} // namespace
TEST_CASE("Quad connection checks use child IDs and preserve original split ordinals",
          "[UI-D01][UI-D05][UI-D06][N4-03][negative][verification]") {
    const std::array p{PointId{1}, PointId{2}, PointId{3}, PointId{4}};
    const std::array b{binding(), binding(), binding(), binding()};
    const std::array<double, 12> xyz{0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0};
    const std::array edges{EdgeInput{EdgeId{1}, {p[0], p[1]}}, EdgeInput{EdgeId{2}, {p[1], p[2]}},
                           EdgeInput{EdgeId{3}, {p[2], p[3]}}, EdgeInput{EdgeId{4}, {p[3], p[0]}},
                           EdgeInput{EdgeId{5}, {p[0], p[2]}}};
    const std::array faces{FaceInput{PhysicalFaceId{50}, BandId{100}},
                           FaceInput{PhysicalFaceId{51}, BandId{100}}};
    const std::array bands{BandId{100}};
    std::array rows{QuadInput{PhysicalFaceId{999},
                              p,
                              {EdgeId{1}, EdgeId{2}, EdgeId{3}, EdgeId{4}},
                              EdgeId{5},
                              {FacetId{90}, FacetId{100}}},
                    QuadInput{PhysicalFaceId{999},
                              p,
                              {EdgeId{1}, EdgeId{2}, EdgeId{3}, EdgeId{4}},
                              EdgeId{5},
                              {FacetId{70}, FacetId{10}}}};
    core::ResourcePlan prep{budget};
    auto w = workspace(prep);
    for (unsigned order = 0; order < 2; ++order) {
        GeometryInput input{p,
                            b,
                            SinglePositions{xyz},
                            SinglePositions{xyz},
                            edges,
                            faces,
                            bands,
                            {},
                            rows,
                            SingleWord{0},
                            context()};
        core::ResourcePlan plan{budget};
        expected_problem(prepare_owned_geometry(input, limits(), *w, plan),
                         core::FailureCode::InvalidInput, FacetId{10});
        const std::array triangle{TriangleInput{FacetId{900},
                                                PhysicalFaceId{999},
                                                {p[0], p[2], p[1]},
                                                {EdgeId{5}, EdgeId{2}, EdgeId{1}}}};
        input.triangles = triangle;
        expected_problem(prepare_owned_geometry(input, limits(), *w, plan),
                         core::FailureCode::InvalidInput, FacetId{900});
        input.triangles = {};
        rows[0].face = PhysicalFaceId{50};
        rows[1].face = PhysicalFaceId{51};
        auto owned = take(prepare_owned_geometry(input, limits(), *w, plan));
        CHECK(take(owned.facet_points(FacetId{10})) == std::array{p[0], p[3], p[2]});
        CHECK(std::get<QuadLRScatter>(take(owned.scatter(FacetId{10}))).coefficients ==
              std::array<std::int8_t, 6>{1, 0, 0, 1, 0, 1});
        rows[0].face = PhysicalFaceId{999};
        rows[1].face = PhysicalFaceId{999};
        std::swap(rows[0], rows[1]);
    }
}
TEST_CASE("Finite point diagnostics complete the admitted pass before canonical selection",
          "[UI-D05][UI-D06][N4-02][N4-03][N4-06][negative][verification]") {
    TriangleFixture f;
    f.two_word = true;
    f.points = {PointId{90}, PointId{10}, PointId{30}};
    f.edges = {EdgeInput{EdgeId{10}, {f.points[0], f.points[1]}},
               EdgeInput{EdgeId{20}, {f.points[1], f.points[2]}},
               EdgeInput{EdgeId{30}, {f.points[2], f.points[0]}}};
    f.triangles[0] = TriangleInput{
        FacetId{70}, PhysicalFaceId{50}, f.points, {EdgeId{10}, EdgeId{20}, EdgeId{30}}};
    core::ResourcePlan prep{budget};
    auto w = workspace(prep);
    f.reference[0] = std::numeric_limits<double>::quiet_NaN();
    f.current_low[3] = std::numeric_limits<double>::infinity();
    f.stamp.time = std::numeric_limits<double>::quiet_NaN();
    f.radius = SingleWord{std::numeric_limits<double>::infinity()};
    for (unsigned order = 0; order < 2; ++order) {
        core::ResourcePlan plan{budget};
        expected_problem(prepare_owned_geometry(f.input(), limits(), *w, plan),
                         core::FailureCode::Nonfinite, PointId{10});
        auto low = limits();
        low.max_steps = 17; // Six scalar checks per point, including both words.
        const auto budget_failure = prepare_owned_geometry(f.input(), low, *w, plan);
        failure(budget_failure, core::FailureCode::ResourceExceeded);
        if (!budget_failure)
            CHECK_FALSE(budget_failure.error().problem);
        low.max_steps = 18; // Equality admits the complete scan and finds the canonical bad point.
        expected_problem(prepare_owned_geometry(f.input(), low, *w, plan),
                         core::FailureCode::Nonfinite, PointId{10});
        // MSVC request_stop is nonconst; preserve the shared cancellation state handle.
        // NOLINTNEXTLINE(misc-const-correctness)
        std::stop_source stop;
        CHECK(stop.request_stop());
        low.max_steps = 0;
        const auto cancelled = prepare_owned_geometry(f.input(), low, *w, plan, stop.get_token());
        failure(cancelled, core::FailureCode::Cancelled);
        if (!cancelled)
            CHECK_FALSE(cancelled.error().problem);
        CHECK(plan.allocated_bytes() == 0);
        std::swap(f.points[0], f.points[1]);
        for (std::size_t axis = 0; axis < 3; ++axis) {
            std::swap(f.reference[axis], f.reference[3 + axis]);
            std::swap(f.current[axis], f.current[3 + axis]);
            std::swap(f.reference_low[axis], f.reference_low[3 + axis]);
            std::swap(f.current_low[axis], f.current_low[3 + axis]);
        }
    }
}
TEST_CASE("Alias edge face and geometric failures retain their existing canonical table order",
          "[UI-D05][UI-D06][N4-03][N4-06][negative][verification]") {
    const TriangleFixture fixture;
    core::ResourcePlan prep{budget};
    auto w = workspace(prep);
    std::array points{PointId{11}, PointId{3}, PointId{7}, PointId{90}, PointId{10}};
    std::array bindings{
        binding(), binding(), binding(),
        PointBinding{ReferenceMapId{1}, MappingGeneration{1}, AliasOf{PointId{999}}},
        PointBinding{ReferenceMapId{1}, MappingGeneration{1}, AliasOf{PointId{999}}}};
    const std::array<double, 15> xyz{0, 0, 0, 3, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0};
    std::array edges{fixture.edges[0], fixture.edges[1], fixture.edges[2],
                     EdgeInput{EdgeId{90}, {PointId{999}, PointId{1000}}},
                     EdgeInput{EdgeId{40}, {PointId{999}, PointId{1000}}}};
    std::array faces{fixture.faces[0], FaceInput{PhysicalFaceId{90}, BandId{999}},
                     FaceInput{PhysicalFaceId{40}, BandId{999}}};
    std::array triangles{
        TriangleInput{
            FacetId{90}, PhysicalFaceId{50}, fixture.points, {EdgeId{10}, EdgeId{20}, EdgeId{30}}},
        TriangleInput{
            FacetId{10}, PhysicalFaceId{50}, fixture.points, {EdgeId{10}, EdgeId{20}, EdgeId{30}}}};
    auto degenerate = fixture.current;
    degenerate[7] = 0;
    for (unsigned order = 0; order < 2; ++order) {
        core::ResourcePlan plan{budget};
        auto input = fixture.input();
        input.points = points;
        input.bindings = bindings;
        input.reference = SinglePositions{xyz};
        input.current = SinglePositions{xyz};
        expected_problem(prepare_owned_geometry(input, limits(), *w, plan),
                         core::FailureCode::InvalidInput, PointId{10});
        input = fixture.input();
        input.edges = edges;
        expected_problem(prepare_owned_geometry(input, limits(), *w, plan),
                         core::FailureCode::InvalidInput, EdgeId{40});
        input = fixture.input();
        input.faces = faces;
        expected_problem(prepare_owned_geometry(input, limits(), *w, plan),
                         core::FailureCode::InvalidInput, PhysicalFaceId{40});
        input = fixture.input();
        input.triangles = triangles;
        input.current = SinglePositions{degenerate};
        expected_problem(prepare_owned_geometry(input, limits(), *w, plan),
                         core::FailureCode::GeometryInvalid, FacetId{10});
        std::swap(points[3], points[4]);
        std::swap(bindings[3], bindings[4]);
        std::swap(edges[3], edges[4]);
        std::swap(faces[1], faces[2]);
        std::swap(triangles[0], triangles[1]);
    }
}
