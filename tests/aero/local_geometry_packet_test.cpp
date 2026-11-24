#include "aero/global_geometry.hpp"
#include "aero/local_geometry_bounds_private.hpp"
#include "geometry_snapshot_fixture.hpp"

#include <boost/multiprecision/cpp_int.hpp>

#include <cmath>
#include <thread>

namespace {
using namespace uvlas;
using namespace uvlas::test_geometry;
using namespace subject;
using Rat = boost::multiprecision::cpp_rational;
using Int = boost::multiprecision::cpp_int;
Rat two(int exponent) {
    return exponent >= 0 ? Rat{Int{1} << exponent} : Rat{1} / Rat{Int{1} << (-exponent)};
}
Rat center(const ScaledEnclosure& value) {
    Rat sum = 0;
    const auto words = value.center();
    if (const auto* one = std::get_if<SingleWord>(&words))
        sum = Rat{one->value};
    else if (const auto* pair = std::get_if<DoubleWord>(&words))
        sum = Rat{pair->high} + Rat{pair->low};
    return sum * two(value.exponent());
}
const ScaledEnclosure* packet_value(const LocalGeometryPacket& packet, LocalQuantity quantity) {
    return &take(packet.value(quantity)).get();
}
void contains(const ScaledEnclosure* value, const Rat& exact) {
    REQUIRE(value != nullptr);
    const Rat radius = Rat{value->radius()} * two(value->exponent());
    CHECK(center(*value) - radius <= exact);
    CHECK(center(*value) + radius >= exact);
}
std::unique_ptr<LocalPacketWorkspace> packet_workspace(core::ResourcePlan& plan) {
    auto result = LocalPacketWorkspace::create(plan);
    if (!result) {
        FAIL_CHECK("Packet workspace failed");
        throw std::runtime_error("workspace");
    }
    return std::move(*result);
}
} // namespace
TEST_CASE("Local packet LG-F01 owns exact geometry and verified face incidence",
          "[UI-D01][UI-D03][UI-D04][N4-02][N4-06][verification]") {
    core::ResourcePlan setup{budget}, result_plan{budget};
    auto old = workspace(setup);
    auto current = packet_workspace(setup);
    const TriangleFixture fixture;
    auto geometry = take(prepare_owned_geometry(fixture.input(), limits(), *old, result_plan));
    auto query = take(resolve_query(
        geometry, {QueryId{5}, SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}},
        geometry.context(), limits(), *old, result_plan));
    auto packet =
        take(prepare_local_packet(geometry, query, FacetId{70}, {core::Precision::TwoWord, 0},
                                  limits(), *current, result_plan));
    contains(packet_value(packet, LocalQuantity::PositionX), Rat{3} / 4);
    contains(packet_value(packet, LocalQuantity::PositionY), Rat{1});
    contains(packet_value(packet, LocalQuantity::NormalZ), Rat{1});
    contains(packet_value(packet, LocalQuantity::AreaZ), Rat{6});
    contains(packet_value(packet, LocalQuantity::Length), Rat{5});
    contains(packet_value(packet, LocalQuantity::Distance), Rat{0});
    for (unsigned i = 0; i < 12; ++i) {
        const auto* value = packet_value(packet, static_cast<LocalQuantity>(i));
        REQUIRE(value);
        CHECK(value->radius() == 0);
    }
    CHECK(std::holds_alternative<VerifiedIncidence>(take(packet.qualification())));
    CHECK(take(packet.snapshot_identity()) == take(geometry.identity()));
    CHECK(take(packet.query()).id == QueryId{5});
}

TEST_CASE("Local alias packet verifies the registered canonical midpoint at its own entry",
          "[UI-D03][UW-D07][N4-02][negative][verification]") {
    core::ResourcePlan setup{budget}, packets{budget};
    auto query_workspace = workspace(setup);
    auto local = packet_workspace(setup);
    const TriangleFixture fixture;
    const EdgeId canonical{10}, raw{40};
    const std::array edges{fixture.edges[0], fixture.edges[1], fixture.edges[2],
                           EdgeInput{raw, {fixture.points[1], fixture.points[0]}}};
    auto facets = fixture.triangles;
    facets[0].edges[0] = raw;
    auto specification = fixture.input();
    specification.edges = edges;
    specification.triangles = facets;
    const auto geometry =
        take(prepare_owned_geometry(specification, limits(), *query_workspace, setup));

    const std::array roles{FaceRoleInput{fixture.faces[0].id, AeroFaceRole::Bound}};
    const std::array registered{EdgeBindingInput{canonical, CanonicalEdge{}},
                                EdgeBindingInput{fixture.edges[1].id, CanonicalEdge{}},
                                EdgeBindingInput{fixture.edges[2].id, CanonicalEdge{}},
                                EdgeBindingInput{raw, canonical}};
    auto unregistered = registered;
    unregistered.back().ownership = CanonicalEdge{};
    const GlobalLimits cap{UINT64_C(100000000), numerics::linear::EstimateAdmission{4096}};
    auto global = GlobalWorkspace::create(setup, cap.metadata);
    if (!global)
        FAIL("Global workspace creation failed");
    const auto checked = [&](std::span<const EdgeBindingInput> bindings) {
        auto topology =
            build_global_topology(geometry, {roles, bindings, {}}, cap, **global, setup);
        if (!topology)
            FAIL("External topology could not be qualified");
        auto contact = check_global_geometry(*topology, geometry, cap, **global, setup);
        if (!contact)
            FAIL("Global contact could not be qualified");
        return std::move(*contact);
    };
    const auto active = checked(registered);
    const auto absent = checked(unregistered);
    const auto query_at = [&](const OwnedGeometry& owner, QueryDefinition expression) {
        return take(resolve_query(owner, {QueryId{95}, expression}, owner.context(), limits(),
                                  *query_workspace, setup));
    };
    const auto at_midpoint = query_at(geometry, SourceEdgePoint{canonical, SingleWord{.5}});
    const auto qualify = [&](const GlobalContactState& contact, const ResolvedQuery& query,
                             const GeometryLimits& cap_for_packet = limits(),
                             const std::stop_token& stop = {}) {
        return prepare_aliased_midpoint_packet(contact, query, FacetId{70}, true,
                                               {core::Precision::TwoWord, 1e-10}, cap_for_packet,
                                               *local, packets, stop);
    };
    const auto accepted = take(qualify(active, at_midpoint));
    const auto verify = [&](const LocalGeometryPacket& packet) {
        const auto incidence = take(packet.qualification());
        const auto* member = std::get_if<VerifiedIncidence>(&incidence);
        REQUIRE(member != nullptr);
        CHECK(member->support == std::variant<PointId, EdgeId, FacetId>{canonical});
        contains(packet_value(packet, LocalQuantity::PositionX), Rat{3, 2});
        contains(packet_value(packet, LocalQuantity::PositionY), Rat{0});
    };
    verify(accepted);
    CHECK(take(accepted.steps()) <= limits().max_steps);
    const auto split =
        query_at(geometry, SourceEdgePoint{canonical, DoubleWord{.5 + 0x1p-53, -0x1p-53}});
    verify(take(qualify(active, split)));

    const auto bytes = packets.allocated_bytes();
    failure(qualify(absent, at_midpoint), core::FailureCode::GeometryInvalid);
    const auto from_alias = query_at(geometry, SourceEdgePoint{raw, SingleWord{.5}});
    failure(qualify(active, from_alias), core::FailureCode::GeometryInvalid);
    for (const auto coordinate : {0., .25, 1.}) {
        const auto on_edge = query_at(geometry, SourceEdgePoint{canonical, SingleWord{coordinate}});
        failure(qualify(active, on_edge), core::FailureCode::InvalidInput);
    }
    const auto near_half = query_at(geometry, SourceEdgePoint{canonical, DoubleWord{.5, 0x1p-75}});
    failure(qualify(active, near_half), core::FailureCode::InvalidInput);
    const auto face =
        query_at(geometry, SourceFacePoint{FacetId{70}, SingleWord{.5}, SingleWord{0}});
    failure(qualify(active, face), core::FailureCode::InvalidInput);

    auto changed_reference = fixture.reference;
    changed_reference[8] = .125;
    auto other_specification = specification;
    other_specification.reference = SinglePositions{changed_reference};
    const auto other =
        take(prepare_owned_geometry(other_specification, limits(), *query_workspace, setup));
    const auto foreign = query_at(other, SourceEdgePoint{canonical, SingleWord{.5}});
    failure(qualify(active, foreign), core::FailureCode::StaleGeneration);
    auto short_budget = limits();
    short_budget.max_steps = 1;
    failure(qualify(active, at_midpoint, short_budget), core::FailureCode::ResourceExceeded);
    // NOLINTNEXTLINE(misc-const-correctness) MSVC 14.51 request_stop is non-const.
    std::stop_source cancel;
    REQUIRE(cancel.request_stop());
    failure(qualify(active, at_midpoint, limits(), cancel.get_token()),
            core::FailureCode::Cancelled);
    CHECK(packets.allocated_bytes() == bytes);
    CHECK(active.active());
    CHECK(absent.active());
    CHECK(at_midpoint.active());
    verify(accepted);
    verify(take(qualify(active, at_midpoint)));
}

TEST_CASE("Local packet LG-F02 F03 and F05 distinguish self contact and detached",
          "[UI-D03][UI-D05][UI-D06][N4-03][negative][verification]") {
    for (const double radius : {0., .25}) {
        core::ResourcePlan plan{budget};
        auto old = workspace(plan);
        auto current = packet_workspace(plan);
        TriangleFixture fixture;
        fixture.radius = SingleWord{radius};
        auto geometry = take(prepare_owned_geometry(fixture.input(), limits(), *old, plan));
        auto material = take(resolve_query(
            geometry, {QueryId{5}, SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}},
            geometry.context(), limits(), *old, plan));
        if (radius == 0)
            failure(prepare_local_packet(geometry, material, FacetId{70},
                                         {core::Precision::TwoWord, 1e-10}, limits(), *current,
                                         plan),
                    core::FailureCode::GeometryInvalid);
        for (const double height : {0., 2.}) {
            auto query =
                take(resolve_query(geometry, {QueryId{6}, FreePoint{FreeSingle{{.75, 1., height}}}},
                                   geometry.context(), limits(), *old, plan));
            auto result =
                prepare_local_packet(geometry, query, FacetId{70},
                                     {core::Precision::TwoWord, 1e-10}, limits(), *current, plan);
            if (height == 0)
                failure(result, core::FailureCode::GeometryInvalid);
            else {
                auto packet = take(std::move(result));
                CHECK(std::holds_alternative<Detached>(take(packet.qualification())));
                contains(packet_value(packet, LocalQuantity::Distance), Rat{2});
            }
        }
    }
}
TEST_CASE("Local square root contains independent rational squares and rejects half-width budget",
          "[UI-D06][N4-06][N4-08][verification][negative]") {
    LocalRootScratch scratch;
    LocalWork work{0, 10000, {}};
    const ExactInteger n{3}, d{2};
    auto root = LocalBounds::root(n, d, 0, false, core::Precision::TwoWord, scratch, work);
    REQUIRE(root);
    const Rat h = center(*root), r = Rat{root->radius()} * two(root->exponent());
    CHECK(h - r >= 0);
    CHECK((h - r) * (h - r) <= Rat{3} / 2);
    CHECK((h + r) * (h + r) >= Rat{3} / 2);
    CHECK(r < Rat{1} / Rat{Int{1} << 100});
    const ExactInteger one{1};
    // A direct dyadic conversion has a nonzero discarded residual in Binary64.
    ExactInteger value{1};
    value <<= 106;
    value += 1;
    auto rounded = LocalBounds::dyadic(value, -106, core::Precision::Binary64, scratch);
    REQUIRE(rounded);
    REQUIRE(rounded->radius() > 0);
    const double radius = std::scalbn(rounded->radius(), rounded->exponent());
    CHECK_FALSE(LocalBounds::width(*rounded, radius, one, false, true, scratch));
    CHECK(LocalBounds::width(*rounded, 2 * radius, one, false, true, scratch));
}

TEST_CASE("Local packet LG-F06 F07 F14 use finite distance exact far and zero-width policy",
          "[UI-D03][UI-D05][UI-D06][N4-06][verification][negative]") {
    for (const bool right_triangle : {false, true}) {
        core::ResourcePlan plan{budget};
        auto old = workspace(plan);
        auto current = packet_workspace(plan);
        TriangleFixture fixture;
        if (right_triangle) {
            fixture.reference = {0, 0, 0, 2, 0, 0, 0, 2, 0};
            fixture.current = fixture.reference;
        }
        auto geometry = take(prepare_owned_geometry(fixture.input(), limits(), *old, plan));
        if (right_triangle) {
            auto query =
                take(resolve_query(geometry, {QueryId{8}, FreePoint{FreeSingle{{3, 3, 0}}}},
                                   geometry.context(), limits(), *old, plan));
            auto packet = take(prepare_local_packet(geometry, query, FacetId{70},
                                                    {core::Precision::TwoWord, 1e-10}, limits(),
                                                    *current, plan));
            const auto* value = packet_value(packet, LocalQuantity::Distance);
            const Rat c = center(*value), r = Rat{value->radius()} * two(value->exponent());
            CHECK((c - r) * (c - r) <= 8);
            CHECK((c + r) * (c + r) >= 8);
            auto material = take(resolve_query(
                geometry,
                {QueryId{9}, SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}},
                geometry.context(), limits(), *old, plan));
            failure(prepare_local_packet(geometry, material, FacetId{70},
                                         {core::Precision::TwoWord, 0}, limits(), *current, plan),
                    core::FailureCode::PrecisionUnresolved);
        } else {
            const std::array heights{std::nextafter(20., 0.), 20., std::nextafter(20., 21.)};
            const std::array expected{FarRelation::Near, FarRelation::Boundary, FarRelation::Far};
            for (std::size_t i = 0; i < heights.size(); ++i) {
                auto query = take(resolve_query(
                    geometry, {QueryId{8}, FreePoint{FreeSingle{{.75, 1., heights[i]}}}},
                    geometry.context(), limits(), *old, plan));
                auto packet = take(prepare_local_packet(geometry, query, FacetId{70},
                                                        {core::Precision::TwoWord, 1e-10}, limits(),
                                                        *current, plan));
                CHECK(take(packet.far()) == expected[i]);
                contains(packet_value(packet, LocalQuantity::Distance), Rat{heights[i]});
            }
        }
    }
}
TEST_CASE("Local packet LG-F08 F10 F11 F12 preserve thin scaled low-word and rotated geometry",
          "[UI-D01][UI-D03][UI-D06][N4-06][verification]") {
    for (const int scenario : {0, 1, 2, 3, 4}) {
        core::ResourcePlan plan{budget};
        auto old = workspace(plan);
        auto current = packet_workspace(plan);
        TriangleFixture fixture;
        int scale = 0;
        if (scenario == 0) {
            fixture.reference = {0, 0, 0, 1, 0, 0, 1, 0x1p-60, 0};
            fixture.current = fixture.reference;
        }
        if (scenario == 1 || scenario == 2) {
            scale = scenario == 1 ? -900 : 900;
            for (auto& x : fixture.reference)
                x = std::scalbn(x, scale);
            fixture.current = fixture.reference;
            fixture.radius = SingleWord{std::scalbn(.25, scale)};
        }
        if (scenario == 3) {
            fixture.two_word = true;
            fixture.reference_low = fixture.reference;
            fixture.current_low = fixture.current;
            fixture.reference.fill(0x1p500);
            fixture.current = fixture.reference;
        }
        if (scenario == 4) {
            for (std::size_t i = 0; i < 3; ++i) {
                fixture.current[3 * i + 1] = -fixture.current[3 * i + 1];
                fixture.current[3 * i + 2] = -fixture.current[3 * i + 2];
            }
        }
        auto geometry = take(prepare_owned_geometry(fixture.input(), limits(), *old, plan));
        auto query = take(resolve_query(
            geometry, {QueryId{5}, SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}},
            geometry.context(), limits(), *old, plan));
        auto packet =
            take(prepare_local_packet(geometry, query, FacetId{70},
                                      {core::Precision::TwoWord, 1e-10}, limits(), *current, plan));
        contains(packet_value(packet, LocalQuantity::NormalZ), Rat{scenario == 4 ? -1 : 1});
        const Rat area = scenario == 0 ? two(-61) : Rat{scenario == 4 ? -6 : 6} * two(2 * scale);
        contains(packet_value(packet, LocalQuantity::AreaZ), area);
        if (scenario == 3) {
            contains(packet_value(packet, LocalQuantity::PositionX), two(500) + Rat{3} / 4);
            contains(packet_value(packet, LocalQuantity::PositionY), two(500) + 1);
        }
        if (scenario == 1 || scenario == 2)
            contains(packet_value(packet, LocalQuantity::Length), Rat{5} * two(scale));
    }
}
TEST_CASE("Local packet keeps original words and ambient topology without inputs or workspace",
          "[UI-D01][UI-D04][N4-02][N4-03][verification][negative]") {
    core::ResourcePlan plan{budget};
    auto old = workspace(plan);
    auto current = packet_workspace(plan);
    TriangleFixture fixture;
    auto geometry = take(prepare_owned_geometry(fixture.input(), limits(), *old, plan));
    auto query = take(resolve_query(
        geometry, {QueryId{5}, SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}},
        geometry.context(), limits(), *old, plan));
    auto packet = take(prepare_local_packet(
        geometry, query, FacetId{70}, {core::Precision::TwoWord, 1e-10}, limits(), *current, plan));
    const auto layout = take(packet.ambient());
    CHECK(layout.facet == FacetId{70});
    CHECK(layout.face == PhysicalFaceId{50});
    CHECK(layout.band == BandId{100});
    CHECK(layout.triangle_vertices == fixture.points);
    CHECK(layout.triangle_edges == std::array{EdgeId{10}, EdgeId{20}, EdgeId{30}});
    CHECK(layout.positive_boundary_incidence == std::array<std::int8_t, 3>{-1, 1, 1});
    current.reset();
    old.reset();
    fixture.current.fill(999);
    const auto& snapshot = take(packet.geometry()).get();
    const auto words = take(snapshot.point_words(PointId{3}, Configuration::Current));
    CHECK(std::get<SingleWord>(words[0]).value == 3);
    contains(packet_value(packet, LocalQuantity::AreaZ), Rat{6});
    auto moved = std::move(packet);
    // NOLINTNEXTLINE(bugprone-use-after-move): verify documented moved-from rejection.
    failure(packet.value(LocalQuantity::AreaZ), core::FailureCode::InvalidInput);
    // NOLINTNEXTLINE(bugprone-use-after-move): verify documented moved-from rejection.
    failure(packet.geometry(), core::FailureCode::InvalidInput);
    failure(moved.value(static_cast<LocalQuantity>(255)), core::FailureCode::InvalidInput);
    CHECK(take(moved.query()).id == QueryId{5});
}

TEST_CASE("Local packet LG-F04 retains every support and all six source orientations",
          "[UI-D01][UI-D03][UI-D04][UI-D05][N4-06][verification]") {
    std::array<unsigned, 3> order{0, 1, 2};
    do {
        TriangleFixture fixture;
        for (std::size_t i = 0; i < 3; ++i) {
            fixture.triangles[0].points[i] = fixture.points[order[i]];
            const auto a = order[i], b = order[(i + 1) % 3];
            const auto edge = (a == 0 && b == 1) || (a == 1 && b == 0)
                                  ? 0U
                                  : ((a == 1 && b == 2) || (a == 2 && b == 1) ? 1U : 2U);
            fixture.triangles[0].edges[i] = fixture.edges[edge].id;
        }
        core::ResourcePlan plan{budget};
        auto old = workspace(plan);
        auto current = packet_workspace(plan);
        auto geometry = take(prepare_owned_geometry(fixture.input(), limits(), *old, plan));
        unsigned inversions = 0;
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = i + 1; j < 3; ++j)
                if (order[i] > order[j])
                    ++inversions;
        const int sign = inversions % 2 == 0 ? 1 : -1;
        std::array<QueryDefinition, 7> definitions{
            SourceVertex{fixture.points[0]},
            SourceVertex{fixture.points[1]},
            SourceVertex{fixture.points[2]},
            SourceEdgePoint{EdgeId{10}, SingleWord{.25}},
            SourceEdgePoint{EdgeId{20}, SingleWord{.25}},
            SourceEdgePoint{EdgeId{30}, SingleWord{.25}},
            SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}};
        for (std::size_t index = 0; index < definitions.size(); ++index) {
            auto query = take(resolve_query(geometry, {QueryId{index}, definitions[index]},
                                            geometry.context(), limits(), *old, plan));
            auto packet =
                take(prepare_local_packet(geometry, query, FacetId{70},
                                          {core::Precision::TwoWord, 0}, limits(), *current, plan));
            CHECK(std::holds_alternative<VerifiedIncidence>(take(packet.qualification())));
            contains(packet_value(packet, LocalQuantity::NormalZ), Rat{sign});
            contains(packet_value(packet, LocalQuantity::AreaZ), Rat{6 * sign});
            CHECK(take(packet.ambient()).triangle_vertices == fixture.triangles[0].points);
            if (index < 3) {
                CHECK(std::get<PointId>(take(packet.closest())) == fixture.points[index]);
            }
        }
    } while (std::next_permutation(order.begin(), order.end()));
}
TEST_CASE("Local packet owns both quad ambient orders and all LR basis scattering",
          "[UI-D01][UI-D04][UI-D05][N4-02][verification]") {
    const std::array p{PointId{101}, PointId{7}, PointId{55}, PointId{3}};
    const std::array bindings{binding(), binding(), binding(), binding()};
    const std::array<double, 12> xyz{0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0};
    const std::array edges{EdgeInput{EdgeId{1}, {p[0], p[1]}}, EdgeInput{EdgeId{2}, {p[1], p[2]}},
                           EdgeInput{EdgeId{3}, {p[2], p[3]}}, EdgeInput{EdgeId{4}, {p[3], p[0]}},
                           EdgeInput{EdgeId{5}, {p[0], p[2]}}, EdgeInput{EdgeId{6}, {p[2], p[0]}}};
    const std::array faces{FaceInput{PhysicalFaceId{8}, BandId{0}}};
    const std::array bands{BandId{0}};
    const std::array quads{QuadInput{faces[0].id,
                                     p,
                                     {EdgeId{1}, EdgeId{2}, EdgeId{3}, EdgeId{4}},
                                     EdgeId{5},
                                     {FacetId{0}, FacetId{99}}}};
    const GeometryInput input{p,
                              bindings,
                              SinglePositions{xyz},
                              SinglePositions{xyz},
                              edges,
                              faces,
                              bands,
                              {},
                              quads,
                              SingleWord{.25},
                              context()};
    core::ResourcePlan plan{budget};
    auto old = workspace(plan);
    auto current = packet_workspace(plan);
    auto geometry = take(prepare_owned_geometry(input, limits(), *old, plan));
    const std::array roles{FaceRoleInput{faces[0].id, AeroFaceRole::Bound}};
    const std::array edge_bindings{EdgeBindingInput{edges[0].id, CanonicalEdge{}},
                                   EdgeBindingInput{edges[1].id, CanonicalEdge{}},
                                   EdgeBindingInput{edges[2].id, CanonicalEdge{}},
                                   EdgeBindingInput{edges[3].id, CanonicalEdge{}},
                                   EdgeBindingInput{edges[4].id, edges[5].id},
                                   EdgeBindingInput{edges[5].id, CanonicalEdge{}}};
    const GlobalLimits cap{UINT64_C(100000000000), numerics::linear::EstimateAdmission{4096}};
    auto global = GlobalWorkspace::create(plan, cap.metadata);
    if (!global)
        FAIL("Global workspace creation failed");
    const auto topology =
        build_global_topology(geometry, {roles, edge_bindings, {}}, cap, **global, plan);
    if (!topology)
        FAIL("Quad topology could not be qualified");
    const auto contact = check_global_geometry(*topology, geometry, cap, **global, plan);
    if (!contact)
        FAIL("Quad contact could not be qualified");
    const auto midpoint =
        take(resolve_query(geometry, {QueryId{6}, SourceEdgePoint{edges[5].id, SingleWord{.5}}},
                           geometry.context(), limits(), *old, plan));
    for (const auto target : quads[0].facets) {
        const auto aliased = take(prepare_aliased_midpoint_packet(*contact, midpoint, target, true,
                                                                  {core::Precision::TwoWord, 1e-10},
                                                                  limits(), *current, plan, {}));
        CHECK(std::holds_alternative<BoundInterior>(take(aliased.qualification())));
        auto query = take(resolve_query(
            geometry, {QueryId{5}, SourceFacePoint{target, SingleWord{.25}, SingleWord{.25}}},
            geometry.context(), limits(), *old, plan));
        auto packet = take(prepare_local_packet(
            geometry, query, target, {core::Precision::TwoWord, 1e-10}, limits(), *current, plan));
        const auto descriptor = take(packet.ambient());
        CHECK(std::get<std::array<PointId, 4>>(descriptor.physical_vertices) == p);
        CHECK(descriptor.triangle_vertices == take(geometry.facet_points(target)));
        const auto expected = target == FacetId{0} ? std::array<std::int8_t, 6>{1, 0, 0, 1, 1, 0}
                                                   : std::array<std::int8_t, 6>{1, 0, 0, 1, 0, 1};
        CHECK(std::get<QuadLRScatter>(descriptor.scatter).coefficients == expected);
        contains(packet_value(packet, LocalQuantity::NormalZ), Rat{1});
    }
}
TEST_CASE("Local packet preserves zero-radius word sums and source-content failures",
          "[UI-D03][UI-D05][UI-D06][N4-03][N4-06][negative][verification]") {
    for (const auto radius : {DoubleWord{-.25, .5}, DoubleWord{.25, -.25}, DoubleWord{.25, -.5}}) {
        core::ResourcePlan plan{budget};
        auto old = workspace(plan);
        auto current = packet_workspace(plan);
        TriangleFixture fixture;
        fixture.radius = radius;
        auto prepared = prepare_owned_geometry(fixture.input(), limits(), *old, plan);
        if (radius.low == -.5) {
            failure(prepared, core::FailureCode::InvalidInput);
            continue;
        }
        auto geometry = take(std::move(prepared));
        auto query = take(resolve_query(geometry, {QueryId{5}, SourceVertex{PointId{11}}},
                                        geometry.context(), limits(), *old, plan));
        auto result =
            prepare_local_packet(geometry, query, FacetId{70}, {core::Precision::TwoWord, 1e-10},
                                 limits(), *current, plan);
        if (radius.low == -.25)
            failure(result, core::FailureCode::GeometryInvalid);
        else {
            auto packet = take(std::move(result));
            contains(packet_value(packet, LocalQuantity::Radius), Rat{1} / 4);
        }
    }
    core::ResourcePlan plan{budget};
    auto old = workspace(plan);
    auto current = packet_workspace(plan);
    TriangleFixture fixture;
    auto first = take(prepare_owned_geometry(fixture.input(), limits(), *old, plan));
    auto query = take(resolve_query(first, {QueryId{5}, SourceVertex{PointId{11}}}, first.context(),
                                    limits(), *old, plan));
    fixture.current[3] = 4;
    auto second = take(prepare_owned_geometry(fixture.input(), limits(), *old, plan));
    failure(prepare_local_packet(second, query, FacetId{70}, {core::Precision::TwoWord, 1e-10},
                                 limits(), *current, plan),
            core::FailureCode::StaleGeneration);
    auto changed = fixture.input();
    changed.context.time = 2;
    auto third = take(prepare_owned_geometry(changed, limits(), *old, plan));
    failure(prepare_local_packet(third, query, FacetId{70}, {core::Precision::TwoWord, 1e-10},
                                 limits(), *current, plan),
            core::FailureCode::StaleGeneration);
    auto collinear = fixture;
    collinear.reference = {0, 0, 0, 1, 0, 0, 2, 0, 0};
    collinear.current = collinear.reference;
    failure(prepare_owned_geometry(collinear.input(), limits(), *old, plan),
            core::FailureCode::GeometryInvalid);
}

namespace {
std::array<Rat, 2> reference_root(const Rat& square) {
    REQUIRE(square >= 0);
    const Int n = boost::multiprecision::numerator(square);
    const Int d = boost::multiprecision::denominator(square);
    if (n == 0)
        return {Rat{0}, Rat{0}};
    const Int scaled = (n * d) << 512;
    const Int root = boost::multiprecision::sqrt(scaled);
    CHECK(root * root <= scaled);
    CHECK((root + 1) * (root + 1) > scaled);
    const Rat divisor = Rat{d} * two(256);
    const Rat lower = Rat{root} / divisor;
    const Rat upper = root * root == scaled ? lower : Rat{root + 1} / divisor;
    CHECK(lower * lower <= square);
    CHECK(upper * upper >= square);
    CHECK(upper - lower <= two(-256));
    return {lower, upper};
}
void contains_root(const ScaledEnclosure& value, const Rat& square, bool negative = false) {
    auto bounds = reference_root(square);
    if (negative) {
        const Rat low = -bounds[1];
        bounds[1] = -bounds[0];
        bounds[0] = low;
    }
    const Rat radius = Rat{value.radius()} * two(value.exponent());
    CHECK(center(value) - radius <= bounds[0]);
    CHECK(center(value) + radius >= bounds[1]);
}
} // namespace
TEST_CASE("Local packet irrational normals and lengths enclose an independent 256-bit reference",
          "[UI-D06][N4-06][N4-08][verification]") {
    for (const int scale : {-900, 0, 900}) {
        TriangleFixture fixture;
        fixture.reference = {0, 0, 0, 1, .5, .25, -.25, 2, .5};
        for (auto& x : fixture.reference)
            x = std::scalbn(x, scale);
        fixture.current = fixture.reference;
        fixture.radius = SingleWord{std::scalbn(.25, scale)};
        core::ResourcePlan plan{budget};
        auto old = workspace(plan);
        auto current = packet_workspace(plan);
        auto geometry = take(prepare_owned_geometry(fixture.input(), limits(), *old, plan));
        auto query = take(resolve_query(
            geometry, {QueryId{5}, SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.375}}},
            geometry.context(), limits(), *old, plan));
        auto packet =
            take(prepare_local_packet(geometry, query, FacetId{70},
                                      {core::Precision::TwoWord, 1e-10}, limits(), *current, plan));
        const std::array<Rat, 3> e{Rat{1} * two(scale), Rat{1} / 2 * two(scale),
                                   Rat{1} / 4 * two(scale)};
        const std::array<Rat, 3> f{-Rat{1} / 4 * two(scale), Rat{2} * two(scale),
                                   Rat{1} / 2 * two(scale)};
        std::array<Rat, 3> c;
        Rat j2 = 0, longest = 0;
        for (std::size_t i = 0; i < 3; ++i) {
            c[i] = e[(i + 1) % 3] * f[(i + 2) % 3] - e[(i + 2) % 3] * f[(i + 1) % 3];
            j2 += c[i] * c[i];
        }
        for (unsigned edge = 0; edge < 3; ++edge) {
            Rat length = 0;
            for (std::size_t i = 0; i < 3; ++i) {
                const Rat delta = edge == 0 ? e[i] : (edge == 1 ? f[i] - e[i] : f[i]);
                length += delta * delta;
            }
            if (length > longest)
                longest = length;
        }
        for (std::size_t i = 0; i < 3; ++i) {
            contains(packet_value(packet, static_cast<LocalQuantity>(i)),
                     e[i] / 4 + Rat{3} * f[i] / 8);
            contains(packet_value(packet, static_cast<LocalQuantity>(i + 6)), c[i] / 2);
            const auto& normal = *packet_value(packet, static_cast<LocalQuantity>(i + 3));
            contains_root(normal, c[i] * c[i] / j2, c[i] < 0);
            CHECK(Rat{2} * Rat{normal.radius()} * two(normal.exponent()) <= Rat{1e-10});
        }
        const auto& length = *packet_value(packet, LocalQuantity::Length);
        contains_root(length, longest);
        const Rat width = Rat{2} * Rat{length.radius()} * two(length.exponent());
        CHECK(width * width <= Rat{1e-10} * Rat{1e-10} * longest);
    }
}
TEST_CASE("Local conversion preserves discarded subgrid residuals and nonzero word gaps",
          "[UI-D06][N4-06][N4-08][negative][verification]") {
    LocalRootScratch scratch;
    for (const int gap : {60, 500, 1100, 2100}) {
        ExactInteger original{1};
        original <<= static_cast<unsigned>(gap);
        original += 1;
        auto value = LocalBounds::dyadic(original, -gap, core::Precision::TwoWord, scratch);
        REQUIRE(value);
        contains(&*value, Rat{1} + two(-gap));
        const Rat c = center(*value), radius = Rat{value->radius()} * two(value->exponent());
        CHECK(radius <= two(-100));
        if (gap <= 1074)
            CHECK(radius == 0);
        else
            CHECK(radius > 0);
        CHECK(c + radius >= Rat{1} + two(-gap));
    }
}
TEST_CASE("Local packet workspaces run independently and reject a live lease without mutation",
          "[UI-D06][N4-02][N4-03][verification][negative]") {
    core::ResourcePlan plan{budget};
    auto old = workspace(plan);
    const TriangleFixture fixture;
    auto geometry = take(prepare_owned_geometry(fixture.input(), limits(), *old, plan));
    auto query = take(resolve_query(geometry, {QueryId{5}, SourceVertex{PointId{11}}},
                                    geometry.context(), limits(), *old, plan));
    std::array<std::optional<LocalGeometryPacket>, 2> packets;
    std::array<bool, 2> failed{false, false};
    const auto worker = [&](std::size_t i) {
        try {
            core::ResourcePlan own{budget};
            auto w = packet_workspace(own);
            auto result = prepare_local_packet(geometry, query, FacetId{70},
                                               {core::Precision::TwoWord, 0}, limits(), *w, own);
            if (result)
                packets[i].emplace(std::move(*result));
            else
                failed[i] = true;
        } catch (...) {
            failed[i] = true;
        }
    };
    std::thread first{worker, 0}, second{worker, 1};
    first.join();
    second.join();
    for (const bool worker_failed : failed)
        REQUIRE_FALSE(worker_failed);
    const auto& left = packets[0];
    const auto& right = packets[1];
    if (!left || !right) {
        FAIL_CHECK("A completed worker did not publish its owned packet");
        return;
    }
    contains(packet_value(*left, LocalQuantity::AreaZ), Rat{6});
    contains(packet_value(*right, LocalQuantity::AreaZ), Rat{6});
    CHECK(center(*packet_value(*left, LocalQuantity::Length)) ==
          center(*packet_value(*right, LocalQuantity::Length)));
}
TEST_CASE(
    "Local packet covariance uses the fixed geometric tolerance without widening qualification",
    "[UI-D01][UI-D06][N4-06][verification]") {
    for (const double scale : {1e-3, 1e3}) {
        TriangleFixture fixture;
        for (std::size_t i = 0; i < 3; ++i) {
            const double y = fixture.reference[3 * i + 1];
            fixture.current[3 * i] = scale * fixture.reference[3 * i] + .125;
            fixture.current[3 * i + 1] = scale * .6 * y - 3;
            fixture.current[3 * i + 2] = scale * .8 * y + 7;
        }
        fixture.reference = fixture.current;
        fixture.radius = SingleWord{scale * .25};
        core::ResourcePlan plan{budget};
        auto old = workspace(plan);
        auto current = packet_workspace(plan);
        auto geometry = take(prepare_owned_geometry(fixture.input(), limits(), *old, plan));
        auto query = take(resolve_query(
            geometry, {QueryId{5}, SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}},
            geometry.context(), limits(), *old, plan));
        auto packet =
            take(prepare_local_packet(geometry, query, FacetId{70},
                                      {core::Precision::TwoWord, 1e-10}, limits(), *current, plan));
        const Rat tolerance{1e-8};
        CHECK(abs(center(*packet_value(packet, LocalQuantity::NormalY)) + Rat{.8}) <= tolerance);
        CHECK(abs(center(*packet_value(packet, LocalQuantity::NormalZ)) - Rat{.6}) <= tolerance);
        CHECK(abs(center(*packet_value(packet, LocalQuantity::Length)) / Rat{scale} - 5) <=
              tolerance);
    }
}
TEST_CASE("Local packet rejects malformed accuracy and preserves existing values",
          "[UI-D06][N4-03][N4-06][negative][verification]") {
    core::ResourcePlan plan{budget};
    auto old = workspace(plan);
    auto current = packet_workspace(plan);
    const TriangleFixture fixture;
    auto geometry = take(prepare_owned_geometry(fixture.input(), limits(), *old, plan));
    auto query = take(resolve_query(geometry, {QueryId{5}, SourceVertex{PointId{11}}},
                                    geometry.context(), limits(), *old, plan));
    const auto make = [&](LocalAccuracy accuracy) {
        return prepare_local_packet(geometry, query, FacetId{70}, accuracy, limits(), *current,
                                    plan);
    };
    failure(make({core::Precision::TwoWord, -1}), core::FailureCode::InvalidInput);
    failure(make({core::Precision::TwoWord, std::numeric_limits<double>::quiet_NaN()}),
            core::FailureCode::Nonfinite);
    failure(make({static_cast<core::Precision>(255), 1e-10}), core::FailureCode::InvalidInput);
    auto normal = take(make({core::Precision::Binary64, 1e-10}));
    CHECK(
        std::holds_alternative<SingleWord>(packet_value(normal, LocalQuantity::Length)->center()));
    failure(prepare_local_packet(geometry, query, FacetId{999}, {core::Precision::TwoWord, 1e-10},
                                 limits(), *current, plan),
            core::FailureCode::InvalidInput);
    contains(packet_value(normal, LocalQuantity::Length), Rat{5});
}
TEST_CASE("Local packet instrumentation leaves all qualified numerical words unchanged",
          "[UI-D06][N4-06][N4-09][verification]") {
    core::ResourcePlan plan{budget};
    auto old = workspace(plan);
    auto current = packet_workspace(plan);
    const TriangleFixture fixture;
    auto geometry = take(prepare_owned_geometry(fixture.input(), limits(), *old, plan));
    auto query = take(resolve_query(geometry, {QueryId{5}, FreePoint{FreeSingle{{.75, 1, 2}}}},
                                    geometry.context(), limits(), *old, plan));
    auto plain = take(prepare_local_packet(geometry, query, FacetId{70},
                                           {core::Precision::TwoWord, 1e-10, false}, limits(),
                                           *current, plan));
    auto measured = take(prepare_local_packet(geometry, query, FacetId{70},
                                              {core::Precision::TwoWord, 1e-10, true}, limits(),
                                              *current, plan));
    for (unsigned i = 0; i < 12; ++i) {
        const auto& a = *packet_value(plain, static_cast<LocalQuantity>(i));
        const auto& b = *packet_value(measured, static_cast<LocalQuantity>(i));
        CHECK(a.exponent() == b.exponent());
        CHECK(std::bit_cast<std::uint64_t>(a.radius()) == std::bit_cast<std::uint64_t>(b.radius()));
        const auto aw = std::get<DoubleWord>(a.center()), bw = std::get<DoubleWord>(b.center());
        CHECK(std::bit_cast<std::uint64_t>(aw.high) == std::bit_cast<std::uint64_t>(bw.high));
        CHECK(std::bit_cast<std::uint64_t>(aw.low) == std::bit_cast<std::uint64_t>(bw.low));
    }
    CHECK(take(plain.steps()) == take(measured.steps()));
    CHECK(take(plain.closest()) == take(measured.closest()));
    CHECK(take(plain.far()) == take(measured.far()));
    CHECK(take(plain.query()).id == take(measured.query()).id);
    CHECK_FALSE(take(plain.timings()).has_value());
    const auto timing = take(measured.timings());
    if (!timing) {
        FAIL_CHECK("Missing enabled timing report");
        return;
    }
    CHECK(timing->classification.calls == 1);
    CHECK(timing->qualification.calls == 1);
    CHECK(timing->enclosure.calls == 1);
    CHECK(timing->root.calls == 5);
    CHECK(timing->width.calls == 12);
}
TEST_CASE("Local root enclosure retains uncertainty below the output word grid",
          "[UI-D06][N4-06][N4-08][negative][verification]") {
    LocalRootScratch scratch;
    LocalWork work{0, 1000, {}};
    ExactInteger n{1};
    n <<= 300;
    const ExactInteger d = n;
    n += 1;
    auto value = LocalBounds::root(n, d, 0, false, core::Precision::TwoWord, scratch, work);
    REQUIRE(value);
    contains_root(*value, Rat{Int{n}} / Rat{Int{d}});
    CHECK(value->radius() > 0);
    const ExactInteger unit_length = ExactInteger{1} << 4296;
    ExactInteger input{1};
    input <<= 106;
    input += 1;
    auto rounded = LocalBounds::dyadic(input, -106, core::Precision::Binary64, scratch);
    REQUIRE(rounded);
    const double radius = std::scalbn(rounded->radius(), rounded->exponent());
    for (const bool area : {false, true}) {
        CHECK_FALSE(LocalBounds::width(*rounded, radius, unit_length, area, false, scratch));
        CHECK(LocalBounds::width(*rounded, 2 * radius, unit_length, area, false, scratch));
    }
}
TEST_CASE("Removing only common powers of two keeps every certified output word and work count",
          "[UI-D06][N4-06][N4-09][verification]") {
    const ExactInteger n{3}, d{2};
    ExactInteger lifted_n = n, lifted_d = d;
    lifted_n <<= 12000;
    lifted_d <<= 12000;
    LocalRootScratch scratch;
    LocalWork first{0, 1000, {}}, second{0, 1000, {}};
    auto a = LocalBounds::root(n, d, 0, false, core::Precision::TwoWord, scratch, first);
    auto b =
        LocalBounds::root(lifted_n, lifted_d, 0, false, core::Precision::TwoWord, scratch, second);
    REQUIRE(a);
    REQUIRE(b);
    contains_root(*a, Rat{3} / 2);
    contains_root(*b, Rat{3} / 2);
    const auto aw = std::get<DoubleWord>(a->center()), bw = std::get<DoubleWord>(b->center());
    CHECK(std::bit_cast<std::uint64_t>(aw.high) == std::bit_cast<std::uint64_t>(bw.high));
    CHECK(std::bit_cast<std::uint64_t>(aw.low) == std::bit_cast<std::uint64_t>(bw.low));
    CHECK(std::bit_cast<std::uint64_t>(a->radius()) == std::bit_cast<std::uint64_t>(b->radius()));
    CHECK(a->exponent() == b->exponent());
    CHECK(first.used == second.used);
}
