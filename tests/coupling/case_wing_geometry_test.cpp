#include "aero/geometry_snapshot_fixture.hpp"
#include "coupling/case_wing_geometry.hpp"
#include "io/case_schema_test_adapter.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <numbers>
#include <ranges>

namespace {
using namespace uvlas;
namespace a = aero::detail;
namespace c = coupling::detail;
namespace p = a::pchip;
namespace tokens = test::schema_adapter;
using enum core::FailureCode;
using tokens::JsonToken;

constexpr std::uint64_t capacity = 1ULL << 42U;
const numerics::linear::BackendMemoryLimit metadata = numerics::linear::EstimateAdmission{capacity};
constexpr p::PchipInputStamp stamp{p::PchipIssuerId{41}, p::PchipInputGeneration{6}};
constexpr std::array<std::uint32_t, 1> unit_word{1};
constexpr std::array<std::uint32_t, 2> tight_denominator{0, 1U << 12U};
constexpr std::array<std::uint32_t, 3> two_word_denominator{0, 0, 1U << 16U};
constexpr p::RationalView one{1, unit_word, unit_word}, zero{0, {}, unit_word};

structure::ExactLimits exact() {
    return {100000,
            100000,
            100000,
            64,
            64,
            100000,
            200000,
            structure::ExactEstimateAdmission{capacity},
            structure::ExactTimingMode::Disabled};
}

c::CaseWingLimits limits(std::stop_token stop = {}) {
    const p::CurveLimits curve{{256, 100000, capacity, 30000}, metadata};
    return {{exact(), capacity, capacity, stop},
            {{exact(), 100000, c::case_motion_frame_interval_scratch_bytes},
             {exact(), curve, 10000, 100000},
             2048,
             10000000},
            {{exact(), curve, 10000}, 100, 10000},
            capacity,
            {{10000, 10000, 10000, 10000, 10000, 10000}, capacity, {metadata, capacity}},
            {std::numeric_limits<std::uint64_t>::max() / 4, metadata}};
}

c::CaseWingProjection projection(core::Precision format) {
    const auto denominator = format == core::Precision::Binary64
                                 ? std::span<const std::uint32_t>{tight_denominator}
                                 : std::span<const std::uint32_t>{two_word_denominator};
    const c::CaseFlowComponentBudget component{one, {1, unit_word, denominator}, zero};
    const std::array vector{component, component, component};
    c::CaseWingProjection selection{format, {2048, capacity}, vector, vector,    vector,
                                    vector, vector,           {},     component, component};
    selection.motion_pose.fill(component);
    return selection;
}

template <class T>
T take(core::Result<T> result) {
    if (!result) {
        if (const auto& cap = result.error().details.resources)
            INFO("resource " << static_cast<int>(cap->kind) << " needs " << cap->required
                             << " observed " << cap->observed);
        FAIL("Wing geometry failed with " << static_cast<int>(result.error().code));
    }
    return std::move(*result);
}

template <class T>
void reject(const core::Result<T>& result, core::FailureCode expected) {
    REQUIRE_FALSE(result);
    CHECK(result.error().code == expected);
}

io::SiCase read(std::string_view wire) {
    core::ResourcePlan plan{capacity};
    const io::InputLimits syntax{
        1000000, 256,     100000,  100000,
        1000000, 1000000, 1000000, {io::ParserAllocationMode::EstimateAdmission, capacity}};
    std::optional<io::InputDiagnostic> input_diag;
    const auto doc =
        take(io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                     io::SourceId{19}, syntax, plan, input_diag));
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> schema_caps{};
    schema_caps.fill(capacity);
    io::SchemaDiagnostic diagnostic;
    const auto schema =
        io::validate_case_schema(doc, io::SchemaLimits{schema_caps}, plan, &diagnostic);
    if (!schema)
        FAIL("Case rejected at " << diagnostic.path().value_or(""));
    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> si_caps{};
    si_caps.fill(capacity);
    const io::SiLimits requirements{si_caps, exact()};
    const auto catalog = take(io::make_si_catalog(*schema, requirements, plan, nullptr));
    std::vector<io::SiRequirementRow> allocations;
    for (std::uint64_t i = 0; i < take(catalog.size()); ++i) {
        const auto entry = take(catalog.entry(i));
        allocations.push_back({entry.key,
                               {io::SiDecimalInput{"1"}, entry.recipe.dimensions},
                               io::SiDecimalInput{"1e-25"},
                               io::SiDecimalInput{"1e-25"}});
    }
    const auto request =
        take(io::make_si_request(catalog, allocations, requirements, plan, nullptr));
    return take(io::convert_case_to_si(*schema, catalog, request, requirements, plan, nullptr));
}

std::string frozen(int angle, int chord, int span) {
    const auto name = "d" + std::to_string(angle) + "-m" + std::to_string(chord) + "-n" +
                      std::to_string(span) + ".json";
    const auto path = std::filesystem::path{UVLAS_WING_EDGE_CASES} / name;
    std::ifstream file{path, std::ios::binary};
    REQUIRE(file.is_open());
    return {std::istreambuf_iterator<char>{file}, std::istreambuf_iterator<char>{}};
}

struct Run {
    core::ResourcePlan memory{capacity};
    std::unique_ptr<a::GeometryWorkspace> geometry = test_geometry::workspace(memory);
    std::unique_ptr<a::GlobalWorkspace> global = take(a::GlobalWorkspace::create(memory, metadata));
    p::CurveWorkspace motion;
    a::GeometryContext context = test_geometry::context();

    Run() {
        context.time = 0;
    }

    core::Result<c::CaseWingGeometry> prepare(const io::SiCase& input,
                                              std::optional<c::CaseFlowClockQuery> address = {},
                                              const c::CaseWingLimits& bounds = limits(),
                                              core::Precision format = core::Precision::Binary64,
                                              c::CaseWingDiagnostic* diagnostic = nullptr) {
        return c::CaseWingGeometry::prepare(input, context, address, projection(format), bounds,
                                            *geometry, *global, motion, memory, diagnostic);
    }
};

double word(io::SiEnclosure input) {
    return std::visit(
        [](const auto& x) {
            if constexpr (std::is_same_v<std::decay_t<decltype(x)>, io::ScaledBinary64>)
                return std::ldexp(x.center, static_cast<int>(x.exponent2));
            else
                return std::ldexp(x.high + x.low, static_cast<int>(x.exponent2));
        },
        input);
}

void compare(std::span<const io::SiEnclosure, 3> actual, std::array<double, 3> expected,
             double tolerance = 5e-14) {
    for (std::size_t k = 0; k < expected.size(); ++k)
        CHECK(std::abs(word(actual[k]) - expected[k]) <= tolerance);
}

double cosine_coordinate(std::uint64_t i, int count) {
    return (1.0 - std::cos(std::numbers::pi * static_cast<double>(i) / count)) / 2.0;
}

std::array<double, 3> independent_position(std::uint64_t id, int angle, double xi, double eta) {
    const bool main = id <= 2;
    const double half_span = main ? 4.0 : 2.5;
    const double root = main ? 1.5 : 1.0;
    const double taper = main ? 0.9 : 0.6;
    const double sweep = main ? 0x1.66819a3a0bf7ap-1 : 0x1.74b49cf3902d4p-2;
    const double bank = main && angle != 0 ? 0x1.665a8349d55e1p-4 : 0.0;
    const double y = half_span * (eta - (id % 2 != 0 ? 1.0 : 0.0));
    const double abs_y = std::abs(y);
    return {(main ? 0.0 : 5.0) + sweep * abs_y + xi * (root - taper * abs_y / half_span), y,
            (main ? 0.0 : 0.6) + bank * abs_y};
}

JsonToken rational(std::string numerator, std::string denominator = "1") {
    auto token = tokens::parse(R"({"numerator":"0","denominator":"1"})");
    token.set("numerator", tokens::string(numerator));
    token.set("denominator", tokens::string(denominator));
    return token;
}

JsonToken coordinate(JsonToken x, JsonToken y) {
    auto point = tokens::parse("[]");
    point.elements = {std::move(x), std::move(y)};
    return point;
}

JsonToken pairings(const std::vector<std::array<std::string, 3>>& values) {
    auto knots = tokens::parse("[]");
    for (const auto& row : values) {
        auto knot = tokens::parse("{}");
        const auto parse = [](const std::string& x) {
            const auto i = x.find('/');
            return i == std::string::npos ? rational(x) : rational(x.substr(0, i), x.substr(i + 1));
        };
        knot.set("eta", parse(row[0]));
        knot.set("leading", parse(row[1]));
        knot.set("trailing", parse(row[2]));
        knots.elements.push_back(std::move(knot));
    }
    auto pairing = tokens::parse("{}");
    pairing.set("knots", std::move(knots));
    return pairing;
}

JsonToken& graph(JsonToken& root, std::size_t index = 0) {
    return root.field("model").field("geometry").elements.at(index);
}

JsonToken& surface(JsonToken& root, std::size_t index = 0) {
    return root.field("model").field("aero_surfaces").elements.at(index);
}

void perimeter(JsonToken& geometry, std::vector<JsonToken> leading,
               std::vector<JsonToken> trailing) {
    auto& boundaries = geometry.field("boundaries").elements;
    REQUIRE(boundaries.size() >= 2);
    const auto assign = [](JsonToken& bound, const std::vector<JsonToken>& vertices) {
        auto list = tokens::parse("[]");
        list.elements = vertices;
        bound.set("vertices", std::move(list));
    };
    assign(boundaries[0], leading);
    assign(boundaries[1], trailing);
    auto domain = tokens::parse("[]");
    domain.elements = std::move(leading);
    std::ranges::move(trailing | std::views::reverse, std::back_inserter(domain.elements));
    geometry.set("domain", std::move(domain));
}
} // namespace

TEST_CASE("Typed V-MULTI Cases own all frozen vertices, collocations, seams and TE rows",
          "[GEO-01][PHY-02][AE-02][UW-D01][UW-D03][UW-D06][V-MULTI][verification]") {
    for (const int degrees : {0, 5})
        for (const auto& [m, n] : {std::pair{2, 8}, std::pair{4, 16}, std::pair{8, 32}}) {
            CAPTURE(degrees, m, n);
            const auto source = read(frozen(degrees, m, n));
            Run run;
            const auto prepared = take(run.prepare(source));
            REQUIRE(prepared.validate_for(source, run.context, {}));
            CHECK(prepared.surface_count() == 4);
            CHECK(prepared.sharing().size() == 4);
            CHECK_FALSE(prepared.clock_address());
            for (std::size_t patch = 0; patch < prepared.surface_count(); ++patch) {
                const auto chart = take(prepared.surface(patch));
                CHECK(chart.surface == patch + 1);
                CHECK(chart.geometry == patch + 1);
                CHECK(chart.source_leading == 0);
                CHECK(chart.source_trailing == 1);
                CHECK(chart.leading.size() == 2);
                CHECK(chart.trailing.size() == 2);
                CHECK(chart.knots.size() == 2);
            }
            REQUIRE(prepared.points().size() ==
                    static_cast<std::size_t>(4 * (m + 1) * (n / 2 + 1)));
            REQUIRE(prepared.cells().size() == static_cast<std::size_t>(2 * m * n));
            REQUIRE(prepared.trailing_bands().size() == static_cast<std::size_t>(2 * n));
            const auto& global = take(prepared.global()).get();
            const auto& kinematics = take(prepared.kinematics()).get();
            const auto retained = global.geometry();
            REQUIRE(retained);
            const auto& geometry = retained->get();
            REQUIRE(take(global.reference_topology()).trailing_point_count() ==
                    static_cast<std::size_t>(2 * (n + 1)));
            const auto half = static_cast<std::size_t>(n / 2);
            const auto width = half + 1;
            const auto patch_points = static_cast<std::size_t>(m + 1) * width;
            const auto patch_cells = static_cast<std::size_t>(m) * half;
            for (std::size_t k = 0; k < prepared.points().size(); ++k) {
                const auto& point = prepared.points()[k];
                const auto patch = k / patch_points;
                const auto i = k % patch_points / width, j = k % width;
                const auto low = cosine_coordinate(i, m);
                const auto xi = i == static_cast<std::size_t>(m)
                                    ? 1.0
                                    : (3.0 * low + cosine_coordinate(i + 1, m)) / 4.0;
                const auto eta = cosine_coordinate(j, n / 2);
                const auto expected = independent_position(patch + 1, degrees, xi, eta);
                CHECK(point.original == a::PointId{k + 1});
                CHECK(point.surface == patch + 1);
                CHECK(point.chord == i);
                CHECK(point.span == j);
                CHECK(point.owner_chord == std::min(i, static_cast<std::size_t>(m - 1)));
                CHECK(point.owner_span == (j == 0 ? 0 : j - 1));
                CHECK(point.last_incident_span == std::min(j, half - 1));
                compare(point.reference, expected);
                compare(point.position, expected);
                compare(point.velocity, {0, 0, 0}, 0);
                const auto checked_velocity = kinematics.velocity(point.original, 0);
                REQUIRE(checked_velocity);
                CHECK(checked_velocity->sign == 0);
                CHECK(checked_velocity->limbs.empty());
                const auto actual = test_geometry::take(
                    geometry.point_words(point.original, a::Configuration::Current));
                for (std::size_t axis = 0; axis < actual.size(); ++axis)
                    CHECK(std::get<a::SingleWord>(actual[axis]).value ==
                          word(point.position[axis]));
                const bool aliased = patch % 2 == 1 && j == 0;
                CHECK(point.canonical == a::PointId{aliased ? k + 1 - patch_points + half : k + 1});
            }
            for (std::size_t k = 0; k < prepared.cells().size(); ++k) {
                const auto& cell = prepared.cells()[k];
                const auto patch = k / patch_cells;
                const auto i = k % patch_cells / half, j = k % half;
                const auto xi = (cosine_coordinate(i, m) + 3.0 * cosine_coordinate(i + 1, m)) / 4.0;
                const auto eta =
                    (cosine_coordinate(j, n / 2) + cosine_coordinate(j + 1, n / 2)) / 2.0;
                CHECK(cell.surface == patch + 1);
                CHECK(cell.chord == i);
                CHECK(cell.span == j);
                CHECK(cell.face == a::PhysicalFaceId{k + 1});
                REQUIRE(cell.collocation);
                CHECK(cell.collocation->context() == run.context);
                const auto query = test_geometry::take(cell.collocation->definition());
                CHECK(query.id == a::QueryId{k + 1});
                const auto point =
                    std::get<a::FreeSingle>(std::get<a::FreePoint>(query.definition).words);
                for (std::size_t axis = 0; axis < point.words.size(); ++axis)
                    CHECK(point.words[axis] == word(cell.position[axis]));
                CHECK(cell.facets == std::array{a::FacetId{2 * k + 1}, a::FacetId{2 * k + 2}});
                compare(cell.position, independent_position(patch + 1, degrees, xi, eta));
                compare(cell.velocity, {0, 0, 0}, 0);
                const double slope = degrees != 0 && patch < 2
                                         ? (patch == 0 ? -1.0 : 1.0) * 0x1.665a8349d55e1p-4
                                         : 0.0;
                const double length = std::hypot(slope, 1.0);
                compare(cell.normal, {0, -slope / length, 1 / length});
            }
            for (std::size_t k = 0; k < prepared.trailing_bands().size(); ++k) {
                const auto& band = prepared.trailing_bands()[k];
                const auto patch = k / half, j = k % half;
                const auto bound = a::PhysicalFaceId{
                    patch * patch_cells + static_cast<std::size_t>(m - 1) * half + j + 1};
                CHECK(band.surface == patch + 1);
                CHECK(band.source_boundary == 1);
                CHECK(band.span == j);
                CHECK(band.transfer.bound_face == bound);
                CHECK(band.transfer.coefficient == 1);
                CHECK(band.transfer.wake_sign == 1);
                const auto edge_result = global.uvlm_load_edge(band.transfer.edge, capacity);
                REQUIRE(edge_result);
                CHECK(edge_result->trailing);
                const auto edge = edge_result->endpoints;
                CHECK(edge == std::array{prepared
                                             .points()[patch * patch_points +
                                                       static_cast<std::size_t>(m) * width + j]
                                             .canonical,
                                         prepared
                                             .points()[patch * patch_points +
                                                       static_cast<std::size_t>(m) * width + j + 1]
                                             .canonical});
            }
        }
}

TEST_CASE("Explicit pairing repairs unequal span ranges and still qualifies the exact domain",
          "[GEO-01][AE-02][UW-D01][verification][negative]") {
    auto input = tokens::root("aerodynamic_steady");
    const auto qp = [](std::string x, std::string y) {
        return coordinate(rational(x, "4"), rational(y, "4"));
    };
    auto& geometry = graph(input);
    perimeter(geometry, {qp("0", "0"), qp("0", "3")}, {qp("4", "1"), qp("4", "4")});
    Run run;
    reject(run.prepare(read(input.wire())), InvalidInput);
    auto& wing = surface(input);
    const auto affine = pairings({{"0", "0", "0"}, {"1", "1", "1"}});
    wing.set("edge_pairing", affine);
    const auto valid = take(run.prepare(read(input.wire())));
    REQUIRE(valid.cells().size() == 1);
    compare(valid.cells()[0].position, {0.75, 0.5625, 0}, 0);
    compare(valid.cells()[0].chord_tangent, {1, 0.25, 0}, 0);
    compare(valid.cells()[0].span_tangent, {0, 0.75, 0}, 0);
    CHECK(valid.trailing_bands()[0].transfer.coefficient == 1);

    SECTION("reversed surface reverses the area and wake traversal") {
        wing.set("orientation", tokens::string("reversed"));
        const auto reversed =
            take(run.prepare(read(input.wire()), {}, limits(), core::Precision::TwoWord));
        CHECK(reversed.cells()[0].reversed);
        compare(reversed.cells()[0].span_tangent, {0, 0.75, 0}, 0);
        compare(reversed.cells()[0].normal, {0, 0, -1}, 0);
        CHECK(reversed.trailing_bands()[0].transfer.coefficient == 1);
        CHECK(reversed.trailing_bands()[0].transfer.wake_sign == -1);
    }
    SECTION("crossing correspondence") {
        wing.set("edge_pairing", pairings({{"0", "0", "1"}, {"1", "1", "0"}}));
        reject(run.prepare(read(input.wire())), GeometryInvalid);
    }
    SECTION("uncovered polygon area") {
        auto& polygon = geometry.field("domain").elements;
        polygon.insert(polygon.begin() + 2, qp("2", "2"));
        reject(run.prepare(read(input.wire())), GeometryInvalid);
    }
    SECTION("an interior collocation knot requires agreement of both derivatives") {
        wing.set("edge_pairing",
                 pairings({{"0", "0", "0"}, {"1/2", "1/4", "1/4"}, {"1", "1", "1"}}));
        c::CaseWingDiagnostic failure;
        reject(run.prepare(read(input.wire()), {}, limits(), core::Precision::Binary64, &failure),
               GeometryInvalid);
        CHECK(failure.phase == c::CaseWingPhase::Collocation);
        CHECK(failure.surface == 1);
        CHECK(failure.chord == 0);
        CHECK(failure.span == 0);
        CHECK(failure.strip == 1);
        CHECK(failure.other_strip == 0);
    }
    SECTION("a straight landing through a redundant knot preserves the result") {
        wing.set("edge_pairing",
                 pairings({{"0", "0", "0"}, {"1/2", "1/2", "1/2"}, {"1", "1", "1"}}));
        const auto redundant = take(run.prepare(read(input.wire())));
        compare(redundant.cells()[0].position, {0.75, 0.5625, 0}, 0);
        compare(redundant.cells()[0].span_tangent, {0, 0.75, 0}, 0);
    }
}

TEST_CASE("Exact sharing rejects false numerical matches and foreign provenance",
          "[GEO-01][UW-D01][UW-D03][UW-D06][V-MULTI][verification][negative]") {
    const auto original = frozen(0, 2, 8);
    auto input = tokens::parse(original);
    Run run;
    SECTION("wrong sharing orientation") {
        for (unsigned index : {0U, 1U})
            surface(input, index)
                .field("sharing")
                .elements.at(0)
                .set("orientation", tokens::string("reversed"));
        reject(run.prepare(read(input.wire())), GeometryInvalid);
    }
    SECTION("different numbers of named seam edges") {
        surface(input, 1).set("chord", tokens::parse("[0,1]"));
        reject(run.prepare(read(input.wire())), GeometryInvalid);
    }
    SECTION("same endpoints and seam vertices do not prove a curved seam") {
        for (unsigned index : {0U, 1U})
            surface(input, index).set("chord", tokens::parse("[0,1]"));
        // z(s) = s(s - 3/32)(s - 3/8) on the positive-x wing root.
        // The LE, quarter-chord vertex and TE all have z == 0.
        graph(input).set("degree", tokens::parse("[3,0]"));
        graph(input).set(
            "coefficients",
            tokens::parse(
                R"([{"powers":[1,0],"coefficient":{"numerator":"9","denominator":"256"}},{"powers":[2,0],"coefficient":{"numerator":"-15","denominator":"32"}},{"powers":[3,0],"coefficient":{"numerator":"1","denominator":"1"}}])"));
        reject(run.prepare(read(input.wire())), GeometryInvalid);
    }
    SECTION("primary and tail graphs may not be substituted for a retained result") {
        const auto live = read(original);
        const auto prepared = take(run.prepare(live));
        surface(input, 0).set("geometry", tokens::string("3"));
        surface(input, 1).set("geometry", tokens::string("4"));
        surface(input, 2).set("geometry", tokens::string("1"));
        surface(input, 3).set("geometry", tokens::string("2"));
        const auto swapped = read(input.wire());
        reject(prepared.validate_for(swapped, run.context, {}), StaleGeneration);
        const auto changed = take(run.prepare(swapped));
        CHECK(word(changed.points()[0].position[2]) == 0.6);
        auto stale = run.context;
        stale.connectivity = a::ConnectivityGeneration{8};
        reject(prepared.validate_for(live, stale, {}), StaleGeneration);
    }
}

TEST_CASE("Convex strips that turn back into earlier strips cannot publish a lattice",
          "[GEO-01][AE-02][UW-D01][verification][negative]") {
    auto input = tokens::root("aerodynamic_steady");
    auto& geometry = graph(input);
    const auto polygon = geometry.field("domain");
    const auto point = [](int x, int y) {
        return coordinate(rational(std::to_string(x), "8"), rational(std::to_string(y), "8"));
    };
    const std::array path{std::pair{-1, 2},  std::pair{2, 2},  std::pair{2, -2},
                          std::pair{-2, -2}, std::pair{-2, 2}, std::pair{1, 2}};
    std::vector<JsonToken> leading, trailing;
    for (const auto& [x, y] : path) {
        leading.push_back(point(4 + x, 4 + y));
        trailing.push_back(point(4 + 2 * x, 4 + 2 * y));
    }
    perimeter(geometry, std::move(leading), std::move(trailing));
    geometry.set("domain", polygon);
    surface(input).set("edge_pairing", pairings({{"0", "0", "0"}, {"1", "5", "5"}}));
    Run run;
    reject(run.prepare(read(input.wire())), GeometryInvalid);
}

TEST_CASE("Prepared rigid Case geometry binds exact stage and prescribed point velocities",
          "[GEO-01][AE-02][AE-07][UW-D01][UW-D03][verification][negative]") {
    auto input = tokens::root("aerodynamic_prescribed_transient");
    auto& motion = surface(input).field("attachment").field("motion");
    input.field("units").set("angle", tokens::string("deg"));
    motion.set("rotation_vector", tokens::parse(R"({"kind":"constant","value":[0,0,90]})"));
    const auto source = read(input.wire());
    Run run;
    run.context.time = 1;
    const c::CaseFlowClockQuery end{stamp, 0, 1, p::RequestedSide::Left};
    const auto result = take(run.prepare(source, end, limits(), core::Precision::TwoWord));
    REQUIRE(result.validate_for(source, run.context, end));
    REQUIRE(result.cells().size() == 1);
    compare(result.points()[0].reference, {0.25, 0, 0}, 0);
    compare(result.points()[0].position, {1, 0.25, 0});
    compare(result.points()[0].velocity, {3, 0, 0}, 0);
    compare(result.cells()[0].position, {0.5, 0.75, 0});
    compare(result.cells()[0].velocity, {3, 0, 0}, 0);
    compare(result.cells()[0].chord_tangent, {0, 1, 0});
    compare(result.cells()[0].span_tangent, {-1, 0, 0});
    compare(result.cells()[0].normal, {0, 0, 1});
    CHECK(result.trailing_bands()[0].transfer.coefficient == 1);
    auto right = end;
    right.side = p::RequestedSide::Right;
    reject(result.validate_for(source, run.context, right), StaleGeneration);
    reject(run.prepare(source, {}, limits()), StaleGeneration);
    run.context.time = 0;
    reject(run.prepare(source, end, limits()), StaleGeneration);
}

TEST_CASE("A shared seam requires velocity agreement as well as position agreement",
          "[GEO-01][AE-07][UW-D03][verification][negative]") {
    auto input = tokens::root("aerodynamic_prescribed_transient");
    const auto prototype = surface(input).field("attachment");
    const auto frozen_input = tokens::parse(frozen(0, 2, 8));
    input.set("model", frozen_input.field("model"));
    for (unsigned k : {0U, 1U}) {
        surface(input, k).set("attachment", prototype);
        surface(input, k)
            .field("attachment")
            .field("motion")
            .set(
                "translation",
                tokens::parse(
                    R"({"kind":"polynomial","origin":0,"scale":1,"domain":[0,1],"coefficients":[[0,0,0],[0,0,1]]})"));
    }
    const c::CaseFlowClockQuery start{stamp, 0, 0, p::RequestedSide::Right};
    Run run;
    const auto valid = take(run.prepare(read(input.wire()), start));
    const auto& left = valid.points()[4];
    const auto& right = valid.points()[15];
    CHECK(left.canonical == right.canonical);
    compare(left.position, {0.1875, 0, 0});
    compare(right.velocity, {0, 0, 1}, 0);
    const auto& kinematics = take(valid.kinematics()).get();
    const auto left_exact = kinematics.velocity(left.original, 2);
    const auto right_exact = kinematics.velocity(right.original, 2);
    REQUIRE(left_exact);
    REQUIRE(right_exact);
    CHECK(left_exact->sign == 1);
    CHECK(right_exact->sign == 1);
    CHECK(left_exact->first_limb == right_exact->first_limb);
    CHECK(std::ranges::equal(left_exact->limbs, right_exact->limbs));
    surface(input, 1)
        .field("attachment")
        .field("motion")
        .set(
            "translation",
            tokens::parse(
                R"({"kind":"polynomial","origin":0,"scale":1,"domain":[0,1],"coefficients":[[0,0,0],[0,0,2]]})"));
    reject(run.prepare(read(input.wire()), start), GeometryInvalid);
}

TEST_CASE("Cancelled and over-budget preparation leaves a prior owner active",
          "[GEO-01][UW-D06][verification][negative]") {
    const auto source = read(frozen(0, 2, 8));
    Run run;
    const auto prepared = take(run.prepare(source));
    std::stop_source cancelled;
    cancelled.request_stop();
    reject(run.prepare(source, {}, limits(cancelled.get_token())), Cancelled);
    auto tiny = limits();
    tiny.chart.maximum_partition_visits = 1;
    reject(run.prepare(source, {}, tiny), ResourceExceeded);
    tiny = limits();
    tiny.geometry.counts.points = 59;
    reject(run.prepare(source, {}, tiny), ResourceExceeded);
    core::ResourcePlan no_memory{0};
    reject(c::CaseWingGeometry::prepare(source, run.context, {},
                                        projection(core::Precision::Binary64), limits(),
                                        *run.geometry, *run.global, run.motion, no_memory),
           ResourceExceeded);
    REQUIRE(prepared.validate_for(source, run.context, {}));
    REQUIRE(prepared.trailing_bands().size() == 16);
}

TEST_CASE("A matched internal trailing boundary is canonical but releases no wake",
          "[GEO-01][UW-D03][UW-D06][verification]") {
    auto input = tokens::root("aerodynamic_steady");
    auto& model = input.field("model");
    auto second_graph = graph(input);
    second_graph.set("id", tokens::string("2"));
    second_graph.field("frame").set("origin", tokens::parse("[1,0,0]"));
    auto& boundaries = second_graph.field("boundaries").elements;
    std::swap(boundaries[0], boundaries[1]);
    model.field("geometry").elements.push_back(std::move(second_graph));
    auto second_surface = surface(input);
    second_surface.set("id", tokens::string("2"));
    second_surface.set("geometry", tokens::string("2"));
    second_surface.set("orientation", tokens::string("reversed"));
    second_surface.set(
        "sharing",
        tokens::parse(
            R"([{"other":"1","own_boundary":1,"other_boundary":1,"orientation":"same"}])"));
    model.field("aero_surfaces").elements.push_back(std::move(second_surface));
    surface(input).set(
        "sharing",
        tokens::parse(
            R"([{"other":"2","own_boundary":1,"other_boundary":1,"orientation":"same"}])"));
    Run run;
    const auto result = take(run.prepare(read(input.wire())));
    REQUIRE(result.points().size() == 8);
    CHECK(result.points()[6].canonical == result.points()[2].original);
    CHECK(result.points()[7].canonical == result.points()[3].original);
    CHECK(result.trailing_bands().empty());
    CHECK(take(result.global()).get().trailing_edge_count() == 0);
}
