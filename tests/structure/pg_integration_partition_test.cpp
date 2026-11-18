#include "exact_polynomial_detail.hpp"
#include "pg_partition_test_support.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <fstream>
#include <limits>
#include <set>
#include <thread>
#include <tuple>
#include <type_traits>

namespace {
using namespace partition_test;
using Oracle = boost::multiprecision::cpp_rational;
using Json = nlohmann::json;
[[nodiscard]] Oracle actual(const ExactRational& value) {
    return Oracle{detail::ExactFactory::value(value)};
}
[[nodiscard]] Oracle fraction(const Json& value) {
    return Oracle{value.at(0).get<std::string>()} / Oracle{value.at(1).get<std::string>()};
}
[[nodiscard]] Oracle determinant(const PgPartitionPoint& a, const PgPartitionPoint& b,
                                 const PgPartitionPoint& c) {
    return (actual(b.s) - actual(a.s)) * (actual(c.t) - actual(a.t)) -
           (actual(b.t) - actual(a.t)) * (actual(c.s) - actual(a.s));
}
[[nodiscard]] PgIntegrationPartition prepare(Fixture& f) {
    core::ResourcePlan plan{cap};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    return take(PgIntegrationPartition::prepare(f.request(), limits(), workspace, plan));
}
void require_code(Fixture& f, core::FailureCode code) {
    core::ResourcePlan plan{cap};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    const auto result = PgIntegrationPartition::prepare(f.request(), limits(), workspace, plan);
    REQUIRE_FALSE(result);
    CHECK(result.error().code == code);
}
[[nodiscard]] PgPartitionPoint json_point(Fixture& f, const Json& value) {
    return {f.r(value.at(0).at(0).get_ref<const std::string&>(),
                value.at(0).at(1).get_ref<const std::string&>()),
            f.r(value.at(1).at(0).get_ref<const std::string&>(),
                value.at(1).at(1).get_ref<const std::string&>())};
}
void check_topology(const PgIntegrationPartition& result, const std::array<Oracle, 3>& total) {
    const auto vertices = take(result.vertices());
    const auto pieces = take(result.pieces());
    const auto triangles = take(result.triangles());
    const auto edges = take(result.edges());
    const auto sources = take(result.sources());
    const auto indices = take(result.edge_source_indices());
    const auto layers = take(result.layers());
    const auto horizontal = take(result.horizontal_faces());
    const auto volumes = take(result.volumes());
    REQUIRE(horizontal.size() == triangles.size() * layers.size() * 2);
    REQUIRE(volumes.size() == triangles.size() * layers.size());
    for (std::size_t i = 0; i < volumes.size(); ++i) {
        CHECK(volumes[i].triangle == static_cast<core::Index>(i / layers.size()));
        CHECK(volumes[i].layer == static_cast<core::Index>(i % layers.size()));
        CHECK(horizontal[2 * i].triangle == volumes[i].triangle);
        CHECK(horizontal[2 * i + 1].triangle == volumes[i].triangle);
        CHECK(horizontal[2 * i].layer == volumes[i].layer);
        CHECK(horizontal[2 * i + 1].layer == volumes[i].layer);
        CHECK(horizontal[2 * i].side == PgPartitionFaceSide::Lower);
        CHECK(horizontal[2 * i + 1].side == PgPartitionFaceSide::Upper);
    }
    Oracle boundary_area2 = 0;
    std::size_t external = 0;
    for (const auto& edge : edges) {
        REQUIRE(edge.vertices[0] >= 0);
        REQUIRE(edge.vertices[1] > edge.vertices[0]);
        REQUIRE(static_cast<std::size_t>(edge.vertices[1]) < vertices.size());
        bool domain = false;
        for (core::Index j = 0; j < edge.source_count; ++j) {
            const auto source = indices[static_cast<std::size_t>(edge.first_source + j)];
            REQUIRE(source >= 0);
            REQUIRE(static_cast<std::size_t>(source) < sources.size());
            domain = domain || std::holds_alternative<PgPartitionDomainSource>(
                                   sources[static_cast<std::size_t>(source)].descriptor);
        }
        CHECK(edge.external == (edge.left.has_value() != edge.right.has_value()));
        CHECK(edge.external == domain);
        if (edge.left) {
            const auto left_piece = *edge.left;
            CHECK(static_cast<std::size_t>(left_piece) < pieces.size());
        }
        if (edge.right) {
            const auto right_piece = *edge.right;
            CHECK(static_cast<std::size_t>(right_piece) < pieces.size());
        }
        if (edge.external) {
            ++external;
            const auto& a = vertices[static_cast<std::size_t>(edge.vertices[0])];
            const auto& b = vertices[static_cast<std::size_t>(edge.vertices[1])];
            const Oracle moment = actual(a.s) * actual(b.t) - actual(a.t) * actual(b.s);
            boundary_area2 += edge.left ? moment : -moment;
        }
    }
    CHECK(boundary_area2 == 2 * total[0]);
    const auto lateral = take(result.lateral_faces());
    CHECK(lateral.size() == external * layers.size());
    std::set<std::pair<core::Index, core::Index>> unique;
    for (const auto& face : lateral) {
        REQUIRE(face.edge >= 0);
        REQUIRE(static_cast<std::size_t>(face.edge) < edges.size());
        const auto& edge = edges[static_cast<std::size_t>(face.edge)];
        CHECK(edge.external);
        CHECK(face.canonical_direction == edge.left.has_value());
        CHECK(unique.emplace(face.edge, face.layer).second);
    }
}
/// Serialize the actual mathematical result, excluding padding, clocks and resource forecasts.
[[nodiscard]] std::string mathematical_payload(const PgIntegrationPartition& result) {
    Json out = Json::array();
    const auto value = [](const ExactRational& x) {
        const auto& raw = detail::ExactFactory::value(x).backend();
        return Json::array(
            {raw.num().str(0, std::ios_base::dec), raw.denom().str(0, std::ios_base::dec)});
    };
    for (const auto& x : take(result.vertices()))
        out.push_back(Json::array({"vertex", value(x.s), value(x.t)}));
    for (const auto& x : take(result.activity_map()))
        out.push_back(Json::array(
            {"cell", x.cell.first, x.cell.second, value(x.area), x.components, x.active}));
    for (const auto& x : take(result.pieces()))
        out.push_back(
            Json::array({"piece", x.cell.first, x.cell.second, x.component, x.thickness_triangle,
                         x.first_vertex, x.vertex_count, value(x.area)}));
    for (const auto x : take(result.piece_vertex_indices()))
        out.push_back(Json::array({"index", x}));
    for (const auto& x : take(result.edges()))
        out.push_back(Json::array(
            {"edge", x.vertices[0], x.vertices[1], x.left ? Json{*x.left} : Json{},
             x.right ? Json{*x.right} : Json{}, x.first_source, x.source_count, x.external}));
    for (const auto x : take(result.edge_source_indices()))
        out.push_back(Json::array({"source-index", x}));
    for (const auto& x : take(result.triangles()))
        out.push_back(Json::array({"triangle", x.piece, x.vertices[0], x.vertices[1], x.vertices[2],
                                   value(x.twice_area)}));
    for (const auto& x : take(result.volumes()))
        out.push_back(Json::array({"volume", x.triangle, x.layer}));
    for (const auto& x : take(result.horizontal_faces()))
        out.push_back(
            Json::array({"horizontal", x.triangle, x.layer, static_cast<unsigned>(x.side)}));
    for (const auto& x : take(result.lateral_faces()))
        out.push_back(Json::array({"lateral", x.edge, x.layer, x.canonical_direction}));
    for (const auto& x : take(result.edge_traces()))
        out.push_back(Json::array({"trace", x.edge, x.source, x.piece, x.canonical_direction}));
    for (const auto& x : take(result.contacts()))
        out.push_back(Json::array(
            {"contact", x.cell.first, x.cell.second, x.source, x.vertices[0], x.vertices[1]}));
    return out.dump();
}
template <class T>
concept RvalueVertices = requires(const T&& value) { std::move(value).vertices(); };
static_assert(!RvalueVertices<PgIntegrationPartition>);
static_assert(!std::is_default_constructible_v<PgIntegrationPartition>);
static_assert(!std::is_copy_constructible_v<PgIntegrationPartition>);
} // namespace

TEST_CASE("PG integration partition matches independent all-cell moments and full topology",
          "[SH-D07][SH-D03][N4-02][N4-08][verification]") {
    std::ifstream stream{UVLAS_PG_PARTITION_FIXTURE};
    REQUIRE(stream.good());
    const Json fixture = Json::parse(stream);
    REQUIRE(fixture.at("cases").size() == 12);
    std::size_t comparisons = 0;
    for (const auto& entry : fixture.at("cases")) {
        const auto name = entry.at("name").get<std::string>();
        const auto n = entry.at("n").get<core::Index>();
        CAPTURE(name, n);
        Fixture f{n};
        f.domain.clear();
        for (const auto& p : entry.at("domain"))
            f.domain.push_back(json_point(f, p));
        const auto result = prepare(f);
        const auto cells = take(result.activity_map());
        const auto pieces = take(result.pieces());
        const auto vertices = take(result.vertices());
        const auto triangles = take(result.triangles());
        REQUIRE(cells.size() == static_cast<std::size_t>(4 * n * n));
        std::vector<std::array<Oracle, 3>> observed(cells.size());
        for (const auto& tri : triangles) {
            REQUIRE(tri.piece >= 0);
            REQUIRE(static_cast<std::size_t>(tri.piece) < pieces.size());
            const auto& piece = pieces[static_cast<std::size_t>(tri.piece)];
            const auto& a = vertices[static_cast<std::size_t>(tri.vertices[0])];
            const auto& b = vertices[static_cast<std::size_t>(tri.vertices[1])];
            const auto& c = vertices[static_cast<std::size_t>(tri.vertices[2])];
            const Oracle twice = determinant(a, b, c);
            CHECK(twice > 0);
            CHECK(actual(tri.twice_area) == twice);
            auto& value =
                observed[static_cast<std::size_t>(piece.cell.first * 2 * n + piece.cell.second)];
            value[0] += twice / 2;
            value[1] += twice * (actual(a.s) + actual(b.s) + actual(c.s)) / 6;
            value[2] += twice * (actual(a.t) + actual(b.t) + actual(c.t)) / 6;
        }
        std::array<Oracle, 3> total{};
        for (std::size_t i = 0; i < cells.size(); ++i) {
            const auto& expected = entry.at("cells").at(i);
            CHECK(cells[i].cell.first == expected.at("cell").at(0).get<core::Index>());
            CHECK(cells[i].cell.second == expected.at("cell").at(1).get<core::Index>());
            CHECK(actual(cells[i].area) == fraction(expected.at("moments").at(0)));
            CHECK(cells[i].active == (fraction(expected.at("moments").at(0)) > 0));
            CHECK((cells[i].components != 0) == cells[i].active);
            for (std::size_t j = 0; j < 3; ++j) {
                CHECK(observed[i][j] == fraction(expected.at("moments").at(j)));
                total[j] += observed[i][j];
                ++comparisons;
            }
        }
        for (std::size_t j = 0; j < 3; ++j)
            CHECK(total[j] == fraction(entry.at("total_moments").at(j)));
        if (name == "two-components" && n == 1)
            CHECK(cells[2].components == 2);
        check_topology(result, total);
        CHECK(result.validate_for(f.request()));
        CHECK(result.preparation().observed_numeric_operations <=
              result.preparation().numeric_operations);
        CHECK(result.preparation().observed_topology_visits <=
              result.preparation().topology_visits);
        CHECK(result.preparation().observed_value_bits <= result.preparation().value_bits);
    }
    CHECK(comparisons == 1008);
}

TEST_CASE("PG partition preserves continuous different-gradient thickness and partial named traces",
          "[SH-D07][SH-D03][N4-02][N4-08][verification][property]") {
    for (const auto n : {core::Index{1}, core::Index{2}, core::Index{4}}) {
        Fixture f{n};
        f.heights[1] = f.r("2");
        f.heights[2] = f.r("2");
        f.heights[3] = f.r("2");
        f.named.push_back({f.r("1", "4"), f.r("0")});
        f.named.push_back({f.r("3", "4"), f.r("0")});
        f.names.push_back({"partial", f.named});
        const auto cut_a = PgPartitionPoint{f.r("1", "10"), f.r("4", "5")};
        const auto cut_b = PgPartitionPoint{f.r("9", "10"), f.r("1", "5")};
        f.refinements.push_back({std::numeric_limits<std::uint64_t>::max(), 7, cut_a, cut_b});
        const auto result = prepare(f);
        const auto thickness = take(result.thickness());
        REQUIRE(thickness.values.size() == 4);
        for (std::size_t i = 0; i < 4; ++i)
            CHECK(actual(thickness.values[i]) == actual(f.heights[i]));
        const auto sources = take(result.sources());
        const auto edges = take(result.edges());
        const auto vertices = take(result.vertices());
        const auto traces = take(result.edge_traces());
        const auto pieces = take(result.pieces());
        bool saw_cut = false, saw_shared_thickness = false;
        for (const auto& source : sources)
            if (const auto* cut = std::get_if<PgPartitionRefinementSource>(&source.descriptor)) {
                saw_cut = true;
                CHECK(cut->id == std::numeric_limits<std::uint64_t>::max());
                CHECK(cut->version == 7);
            }
        CHECK(saw_cut);
        if (n == 1) {
            const auto before = take(result.identity());
            f.refinements[0].version = 8;
            const auto changed = prepare(f);
            CHECK(take(changed.identity()) != before);
            const auto stale = result.validate_for(f.request());
            REQUIRE_FALSE(stale);
            CHECK(stale.error().code == core::FailureCode::StaleGeneration);
            CHECK(mathematical_payload(changed) == mathematical_payload(result));
            f.refinements[0].version = 7;
        }
        for (const auto& edge : edges)
            if (edge.left && edge.right) {
                const auto a = pieces[static_cast<std::size_t>(*edge.left)].thickness_triangle;
                const auto b = pieces[static_cast<std::size_t>(*edge.right)].thickness_triangle;
                if (a != b) {
                    saw_shared_thickness = true;
                    for (const auto i : edge.vertices)
                        CHECK(actual(vertices[static_cast<std::size_t>(i)].s) ==
                              actual(vertices[static_cast<std::size_t>(i)].t));
                }
            }
        CHECK(saw_shared_thickness);
        Oracle length = 0;
        for (const auto& trace : traces) {
            const auto& edge = edges[static_cast<std::size_t>(trace.edge)];
            const auto& a = vertices[static_cast<std::size_t>(edge.vertices[0])];
            const auto& b = vertices[static_cast<std::size_t>(edge.vertices[1])];
            CHECK(actual(a.t) == 0);
            CHECK(actual(b.t) == 0);
            CHECK(trace.canonical_direction);
            length += actual(b.s) - actual(a.s);
        }
        CHECK(length == Oracle{1} / 2);
        CHECK_FALSE(traces.empty());
        const auto names = take(result.boundaries());
        const auto text = take(result.boundary_text());
        REQUIRE(names.size() == 1);
        CHECK(std::string_view{text.data(), text.size()} == "partial");
        check_topology(result, {Oracle{1}, Oracle{1} / 2, Oracle{1} / 2});
    }
}

TEST_CASE("PG partition uses one common split at thickness T-junctions and material interfaces",
          "[SH-D07][SH-D03][verification]") {
    Fixture f{2};
    f.thickness.push_back({f.r("1", "2"), f.r("1", "2")});
    f.triangles = {{0, 1, 2}, {0, 4, 3}, {4, 2, 3}};
    f.heights[0] = f.r("1");
    f.heights[1] = f.r("2");
    f.heights[2] = f.r("4");
    f.heights[3] = f.r("3");
    f.heights.push_back(f.r("5", "2"));
    const auto result = prepare(f);
    check_topology(result, {Oracle{1}, Oracle{1} / 2, Oracle{1} / 2});
    const auto pieces = take(result.pieces());
    std::array<bool, 3> used{};
    for (const auto& piece : pieces) {
        REQUIRE(piece.thickness_triangle >= 0);
        REQUIRE(piece.thickness_triangle < 3);
        used[static_cast<std::size_t>(piece.thickness_triangle)] = true;
    }
    CHECK(std::all_of(used.begin(), used.end(), [](bool x) { return x; }));
}

TEST_CASE("PG partition returns exact closed-cell zero contacts without activating them",
          "[SH-D07][SH-D03][verification][negative]") {
    Fixture f;
    f.domain[2] = {f.r("1"), f.r("1", "2")};
    f.domain[3] = {f.r("0"), f.r("1", "2")};
    const auto result = prepare(f);
    const auto cells = take(result.activity_map());
    const auto contacts = take(result.contacts());
    const auto vertices = take(result.vertices());
    REQUIRE(cells.size() == 4);
    CHECK(cells[0].active);
    CHECK_FALSE(cells[1].active);
    CHECK(cells[2].active);
    CHECK_FALSE(cells[3].active);
    CHECK_FALSE(contacts.empty());
    for (const auto& contact : contacts) {
        CHECK(contact.cell.second == 1);
        CHECK(actual(vertices[static_cast<std::size_t>(contact.vertices[0])].t) == Oracle{1} / 2);
        CHECK(actual(vertices[static_cast<std::size_t>(contact.vertices[1])].t) == Oracle{1} / 2);
    }
}

TEST_CASE("PG partition rejects every malformed geometry family without repairing input",
          "[SH-D07][N4-08][negative]") {
    SECTION("self crossing") {
        Fixture f;
        std::swap(f.domain[1], f.domain[2]);
        require_code(f, core::FailureCode::GeometryInvalid);
    }
    SECTION("zero edge") {
        Fixture f;
        f.domain[1] = f.point("0", "0");
        require_code(f, core::FailureCode::GeometryInvalid);
    }
    SECTION("outside square") {
        Fixture f;
        f.domain[1].s = f.r("2");
        require_code(f, core::FailureCode::InvalidInput);
    }
    SECTION("folded collinear boundary") {
        Fixture f;
        f.domain.insert(f.domain.begin() + 1, {f.r("1", "2"), f.r("0")});
        f.domain.insert(f.domain.begin() + 2, {f.r("1", "4"), f.r("0")});
        require_code(f, core::FailureCode::GeometryInvalid);
    }
    SECTION("thickness overlap") {
        Fixture f;
        f.triangles.push_back(f.triangles[0]);
        require_code(f, core::FailureCode::GeometryInvalid);
    }
    SECTION("thickness missing") {
        Fixture f;
        f.triangles.pop_back();
        require_code(f, core::FailureCode::GeometryInvalid);
    }
    SECTION("thickness zero") {
        Fixture f;
        f.heights[0] = f.r("0");
        require_code(f, core::FailureCode::InvalidInput);
    }
    SECTION("degenerate thickness") {
        Fixture f;
        f.thickness[2] = f.point("2", "0");
        require_code(f, core::FailureCode::GeometryInvalid);
    }
    SECTION("bad thickness index") {
        Fixture f;
        f.triangles[0][0] = -1;
        require_code(f, core::FailureCode::InvalidInput);
    }
    SECTION("layer gap") {
        Fixture f;
        const auto middle = f.r("0");
        f.layers.pop_back();
        f.layers.push_back({core::LayerId{29}, middle, f.zeta[2]});
        require_code(f, core::FailureCode::InvalidInput);
    }
    SECTION("duplicate layer identity") {
        Fixture f;
        f.layers[1].id = core::LayerId{11};
        require_code(f, core::FailureCode::InvalidInput);
    }
    SECTION("named boundary through interior") {
        Fixture f;
        f.named.push_back(f.point("0", "0"));
        f.named.push_back(f.point("1", "1"));
        f.names.push_back({"invalid", f.named});
        require_code(f, core::FailureCode::GeometryInvalid);
    }
    SECTION("invalid UTF8 name") {
        Fixture f;
        f.named.push_back(f.point("0", "0"));
        f.named.push_back(f.point("1", "0"));
        const std::string bad{"bad\xc0"};
        f.names.push_back({bad, f.named});
        require_code(f, core::FailureCode::InvalidInput);
    }
    SECTION("named foldback") {
        Fixture f;
        f.named.push_back(f.point("0", "0"));
        f.named.push_back(f.point("1", "0"));
        f.named.push_back({f.r("1", "2"), f.r("0")});
        f.names.push_back({"fold", f.named});
        require_code(f, core::FailureCode::GeometryInvalid);
    }
    SECTION("zero refinement") {
        Fixture f;
        const auto p = f.point("0", "0");
        f.refinements.push_back({0, 0, p, p});
        require_code(f, core::FailureCode::InvalidInput);
    }
    SECTION("jump on duplicated seam vertices") {
        Fixture f;
        f.thickness.clear();
        f.heights.clear();
        for (const auto& p : std::array<std::array<std::string_view, 2>, 6>{
                 {{"0", "0"}, {"1", "0"}, {"1", "1"}, {"0", "0"}, {"1", "1"}, {"0", "1"}}}) {
            f.thickness.push_back(f.point(p[0], p[1]));
            f.heights.push_back(f.r("1"));
        }
        f.triangles = {{0, 1, 2}, {3, 4, 5}};
        const auto epsilon =
            take(ExactRational::from_decimal("1e-400", exact_limits(), f.input_plan));
        f.heights[3] = take(f.heights[3].add(epsilon, exact_limits(), f.input_plan));
        require_code(f, core::FailureCode::GeometryInvalid);
    }
}

TEST_CASE(
    "PG partition ownership content identity and workspace parallelism survive caller destruction",
    "[SH-D07][N4-02][N4-08][property][negative]") {
    Fixture f;
    const auto original = prepare(f);
    const auto original_id = take(original.identity());
    const auto original_payload = mathematical_payload(original);
    auto second = prepare(f);
    CHECK(take(second.identity()) == original_id);
    auto moved = std::move(second);
    // NOLINTBEGIN(bugprone-use-after-move): verify documented moved-from query failures.
    CHECK_FALSE(second.identity());
    CHECK_FALSE(second.vertices());
    CHECK_FALSE(second.validate_for(f.request()));
    // NOLINTEND(bugprone-use-after-move)
    CHECK(take(moved.identity()) == original_id);
    f.heights[0] = f.r("2");
    CHECK_FALSE(original.validate_for(f.request()));
    CHECK(take(original.identity()) == original_id);
    f.heights[0] = f.r("1");
    CHECK(original.validate_for(f.request()));
    core::Workspace busy{core::PlannedBuffer<std::byte>{}};
    auto lease = take(busy.acquire());
    core::ResourcePlan rejected{cap};
    const auto reentrant = PgIntegrationPartition::prepare(f.request(), limits(), busy, rejected);
    REQUIRE_FALSE(reentrant);
    CHECK(reentrant.error().code == core::FailureCode::InvalidInput);
    CHECK(rejected.allocated_bytes() == 0);
    std::array<std::optional<PgPartitionIdentity>, 2> parallel;
    std::array<std::string, 2> payloads;
    const auto worker = [&](std::size_t i) {
        core::ResourcePlan plan{cap};
        core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
        const auto result = PgIntegrationPartition::prepare(f.request(), limits(), workspace, plan);
        if (result) {
            parallel[i] = take(result->identity());
            payloads[i] = mathematical_payload(*result);
        }
    };
    std::thread left{worker, 0}, right{worker, 1};
    left.join();
    right.join();
    const auto& a = parallel[0];
    const auto& b = parallel[1];
    if (!a || !b) {
        FAIL("Both independent preparations must return an identity");
        return;
    }
    const auto left_id = *a;
    const auto right_id = *b;
    CHECK(left_id == original_id);
    CHECK(right_id == original_id);
    CHECK(payloads[0] == original_payload);
    CHECK(payloads[1] == original_payload);
    CHECK(mathematical_payload(original) == original_payload);
    const auto owned = []() {
        Fixture temporary;
        return prepare(temporary);
    }();
    CHECK(take(owned.activity_map()).size() == 4);
    CHECK(take(owned.thickness()).values.size() == 4);
    CHECK(take(owned.layers())[1].id == core::LayerId{29});
    CHECK(owned.graph());
}

TEST_CASE("PG preparation timing preserves the full mathematical payload",
          "[SH-D07][N4-02][N4-08][property]") {
    Fixture f;
    const auto unmeasured = prepare(f);
    auto enabled = limits();
    enabled.algebra.timing = ExactTimingMode::Enabled;
    core::ResourcePlan plan{cap};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    const auto measured =
        take(PgIntegrationPartition::prepare(f.request(), enabled, workspace, plan));
    CHECK_FALSE(unmeasured.preparation().timings);
    const auto& timing = measured.preparation().timings;
    REQUIRE(timing.has_value());
    if (!timing)
        return;
    CHECK(timing->whole.calls == 1);
    CHECK(timing->input.calls == 1);
    CHECK(timing->intersections.calls == 1);
    CHECK(timing->planar.calls == 1);
    CHECK(timing->traces.calls == 1);
    CHECK(mathematical_payload(measured) == mathematical_payload(unmeasured));
    CHECK(take(measured.identity()) == take(unmeasured.identity()));
    CHECK(measured.preparation().observed_numeric_operations ==
          unmeasured.preparation().observed_numeric_operations);
}
