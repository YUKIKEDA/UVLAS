#include "aero/geometry_snapshot_fixture.hpp"
#include "coupling/case_wing_geometry.hpp"
#include "coupling/case_wing_harmonic.hpp"
#include "io/case_schema_test_adapter.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <numbers>
#include <ranges>
#include <set>
#include <tuple>

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

TEST_CASE("Four Case wing clock stages share a checked reference topology and trailing rows",
          "[GEO-01][AE-02][AE-07][UW-D03][verification][negative]") {
    auto input = tokens::root("aerodynamic_prescribed_transient");
    input.field("analysis")
        .set("intervals",
             tokens::parse(R"([{"start":0,"end":0.75,"step":0.25,"restart":"initial"}])"));
    auto& model = input.field("model");
    auto fixed_graph = model.field("geometry").elements.front();
    fixed_graph.set("id", tokens::string("2"));
    fixed_graph.field("frame").set("origin", tokens::parse("[0,3,0]"));
    model.field("geometry").elements.push_back(std::move(fixed_graph));
    auto fixed_surface = model.field("aero_surfaces").elements.front();
    fixed_surface.set("id", tokens::string("2"));
    fixed_surface.set("geometry", tokens::string("2"));
    fixed_surface.set("attachment", tokens::parse(R"({"kind":"fixed_reference"})"));
    model.field("aero_surfaces").elements.push_back(std::move(fixed_surface));
    const auto source = read(input.wire());
    Run run;
    const auto initial_context = run.context;
    const c::CaseFlowClockQuery start{stamp, 0, 0, p::RequestedSide::Right};
    const auto anchor = take(run.prepare(source, start, limits(), core::Precision::TwoWord));
    const auto anchor_topology = take(take(anchor.global()).get().reference_topology());
    REQUIRE(anchor.trailing_bands().size() == 2);
    const auto prepare_at = [&](const c::CaseWingGeometry& wing, const io::SiCase& typed,
                                c::CaseFlowClockQuery address, const c::CaseWingLimits& caps) {
        return wing.evaluate_stage(typed, run.context, address,
                                   projection(core::Precision::TwoWord), caps, *run.geometry,
                                   *run.global, run.motion, run.memory);
    };
    const auto check = [&](const c::CaseWingGeometry& wing, c::CaseFlowClockQuery address) {
        const auto& contact = take(wing.global()).get();
        CHECK(take(contact.reference_topology()).same_owner(anchor_topology));
        REQUIRE(wing.validate_for(source, run.context, address));
        REQUIRE(wing.points().size() == 8);
        REQUIRE(wing.cells().size() == 2);
        REQUIRE(wing.trailing_bands().size() == 2);
        const double time = static_cast<double>(address.step) / 4;
        for (const auto& vertex : wing.points()) {
            const auto moving = vertex.surface == 1;
            const double base_x = vertex.chord == 0 ? .25 : 1.0;
            const double base_y = static_cast<double>(vertex.span) + (moving ? 0. : 3.);
            compare(vertex.reference, {base_x, base_y, 0}, 0);
            compare(vertex.position, {base_x + (moving ? time * time * time : 0), base_y, 0}, 0);
            compare(vertex.velocity, {moving ? 3 * time * time : 0, 0, 0}, 0);
        }
        for (std::size_t i = 0; i < wing.cells().size(); ++i) {
            const auto& cell = wing.cells()[i];
            const bool moving = cell.surface == 1;
            compare(cell.position, {.75 + (moving ? time * time * time : 0), moving ? .5 : 3.5, 0},
                    0);
            compare(cell.velocity, {moving ? 3 * time * time : 0, 0, 0}, 0);
            compare(cell.normal, {0, 0, 1}, 0);
            REQUIRE(cell.collocation);
            CHECK(cell.collocation->context() == run.context);
            CHECK(cell.facets == anchor.cells()[i].facets);
            const auto row = take(wing.collocation(i));
            CHECK(row.facet == cell.facets[0]);
            CHECK(&row.query == &*cell.collocation);
            CHECK(cell.face == anchor.cells()[i].face);
            const auto& current = wing.trailing_bands()[i].transfer;
            const auto& retained = anchor.trailing_bands()[i].transfer;
            CHECK(current.edge == retained.edge);
            CHECK(current.bound_face == cell.face);
            CHECK(current.bound_face == retained.bound_face);
            CHECK(current.coefficient == retained.coefficient);
            CHECK(current.wake_sign == retained.wake_sign);
        }
    };
    check(anchor, start);
    std::optional<c::CaseWingGeometry> previous;
    for (std::uint64_t i = 1; i < 4; ++i) {
        run.context.time = static_cast<double>(i) / 4;
        run.context.current = core::GeometryGeneration{initial_context.current.value() + i};
        const c::CaseFlowClockQuery address{
            stamp, 0, i, i == 3 ? p::RequestedSide::Left : p::RequestedSide::Automatic};
        auto current = take(prepare_at(previous ? *previous : anchor, source, address, limits()));
        check(current, address);
        previous.emplace(std::move(current));
    }
    const c::CaseFlowClockQuery end{stamp, 0, 3, p::RequestedSide::Left};
    const auto independent = take(run.prepare(source, end, limits(), core::Precision::TwoWord));
    CHECK_FALSE(
        take(take(independent.global()).get().reference_topology()).same_owner(anchor_topology));
    auto changed_stamp = end;
    changed_stamp.stamp = {p::PchipIssuerId{99}, p::PchipInputGeneration{6}};
    reject(prepare_at(anchor, source, changed_stamp, limits()), StaleGeneration);
    auto changed_case = input;
    changed_case.field("flow").field("uniform").set("value", tokens::parse("[2,0,0]"));
    reject(prepare_at(anchor, read(changed_case.wire()), end, limits()), StaleGeneration);
    const c::CaseWingStraightWake closure{{128, 0, 0}};
    const auto with_closure = take(c::CaseWingGeometry::prepare(
        source, initial_context, start, projection(core::Precision::TwoWord), limits(),
        *run.geometry, *run.global, run.motion, run.memory, nullptr, closure));
    REQUIRE(with_closure.wake_bands().size() == 2);
    c::CaseWingDiagnostic diagnostic;
    const auto staged_closure = with_closure.evaluate_stage(
        source, run.context, end, projection(core::Precision::TwoWord), limits(), *run.geometry,
        *run.global, run.motion, run.memory, &diagnostic);
    reject(staged_closure, Unsupported);
    CHECK(diagnostic.phase == c::CaseWingPhase::TrailingEdge);
    REQUIRE(with_closure.validate_for(source, initial_context, start));
    const auto harmonic_source = read(frozen(0, 2, 8));
    const auto harmonic = take(c::CaseWingGeometry::prepare(
        harmonic_source, initial_context, {}, projection(core::Precision::Binary64), limits(),
        *run.geometry, *run.global, run.motion, run.memory, nullptr, closure));
    REQUIRE_FALSE(harmonic.wake_bands().empty());
    reject(prepare_at(harmonic, harmonic_source, end, limits()), Unsupported);
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    reject(prepare_at(anchor, source, end, limits(stop.get_token())), Cancelled);
    run.context.current = initial_context.current;
    reject(prepare_at(anchor, source, end, limits()), StaleGeneration);
    run.context.current = core::GeometryGeneration{initial_context.current.value() + 3};
    run.context.reference = a::ReferenceGeometryGeneration{initial_context.reference.value() + 1};
    reject(prepare_at(anchor, source, end, limits()), StaleGeneration);
    REQUIRE(anchor.validate_for(source, initial_context, start));
    CHECK(take(take(anchor.global()).get().reference_topology()).same_owner(anchor_topology));

    // The translation crosses the other wing during staging, although the
    // two reference lattices at x=[1/4,1] and [11/8,17/8] are disjoint.
    auto crossing = input;
    crossing.field("model")
        .field("geometry")
        .elements[1]
        .field("frame")
        .set("origin", tokens::parse("[1.125,0,0]"));
    const auto crossing_source = read(crossing.wire());
    run.context = initial_context;
    const auto before_crossing =
        take(run.prepare(crossing_source, start, limits(), core::Precision::TwoWord));
    run.context.time = .75;
    run.context.current = core::GeometryGeneration{initial_context.current.value() + 3};
    reject(prepare_at(before_crossing, crossing_source, end, limits()), GeometryInvalid);
    REQUIRE(before_crossing.validate_for(crossing_source, initial_context, start));
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

namespace {
using Json = nlohmann::json;
using Complex = std::complex<long double>;
using Position = std::array<double, 3>;
using ComplexVector = std::array<Complex, 3>;
namespace ls = numerics::linear;

ls::Context harmonic_context(const a::GeometryContext& geometry) {
    REQUIRE(geometry.diagnostics.trial);
    return {{geometry.fingerprint, geometry.model, geometry.current, core::ChartGeneration{1},
             core::SparsityGeneration{2}, core::NumericalGeneration{3}, core::ScaleGeneration{4},
             geometry.precision, core::Precision::Binary64, core::WakePolicyGeneration{5},
             core::TimeIntervalGeneration{6}, core::BackendGeneration{7}},
            *geometry.diagnostics.trial,
            core::ScopeId{8},
            core::ScopeId{9},
            geometry.time};
}

Complex complex_value(const Json& value) {
    REQUIRE(value.size() == 2);
    return {value[0].get<long double>(), value[1].get<long double>()};
}

long double center(const a::ScaledEnclosure& value) {
    const auto words = value.center();
    if (const auto* single = std::get_if<a::SingleWord>(&words))
        return std::ldexp(static_cast<long double>(single->value), value.exponent());
    const auto pair = std::get<a::DoubleWord>(words);
    return std::ldexp(static_cast<long double>(pair.high) + pair.low, value.exponent());
}

long double radius(const a::ScaledEnclosure& value) {
    return std::ldexp(static_cast<long double>(value.radius()), value.exponent());
}

ComplexVector loads(const c::DiscreteLoadVector& real, const c::DiscreteLoadVector& imag) {
    ComplexVector result;
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = {center(real[i]), center(imag[i])};
    return result;
}

template <class A, class B>
auto cross(const std::array<A, 3>& a, const std::array<B, 3>& b) {
    return std::array{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                      a[0] * b[1] - a[1] * b[0]};
}

Position subtract(const Position& a, const Position& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

void add(ComplexVector& to, const ComplexVector& from) {
    for (std::size_t i = 0; i < to.size(); ++i)
        to[i] += from[i];
}

void compare_load(const ComplexVector& actual, const ComplexVector& expected) {
    long double difference = 0, scale = 0;
    for (std::size_t i = 0; i < actual.size(); ++i) {
        difference += std::norm(actual[i] - expected[i]);
        scale += std::norm(expected[i]);
    }
    REQUIRE(scale > 0);
    CHECK(std::sqrt(difference / scale) < 1e-9L);
}

enum class HarmonicCore { Baseline, Small };

Json harmonic_reference(const char* name, HarmonicCore core = HarmonicCore::Baseline) {
    std::ifstream input{name, std::ios::binary};
    REQUIRE(input);
    const auto result = Json::parse(input);
    REQUIRE(result.at("schema") == "uvlas.vmulti-harmonic-positive-core-companion/1");
    REQUIRE(result.at("core_radius_m") == (core == HarmonicCore::Small ? 1e-7 : 1e-5));
    REQUIRE(result.at("wake_length_m") == 128.);
    return result;
}

std::size_t reference_index(const c::CaseWingCell& cell, std::size_t m, std::size_t n) {
    const auto side = static_cast<std::size_t>(cell.surface - 1);
    return side / 2 * m * n + static_cast<std::size_t>(cell.chord) * n + side % 2 * (n / 2) +
           static_cast<std::size_t>(cell.span);
}

std::vector<Position> independent_points(const c::CaseWingGeometry& wing, const Json& record,
                                         std::size_t n) {
    std::vector<Position> result;
    result.reserve(wing.points().size());
    for (const auto& vertex : wing.points()) {
        const auto side = static_cast<std::size_t>(vertex.surface - 1);
        const auto span = side % 2 * (n / 2) + static_cast<std::size_t>(vertex.span);
        result.push_back(
            record.at("vertices_m").at(side / 2).at(vertex.chord).at(span).get<Position>());
    }
    return result;
}

void verify_original_loads(const c::CaseWingGeometry& wing, const c::HarmonicLoadBreakdown& loads_,
                           const Json& record, const Json& shared, std::size_t drive, std::size_t m,
                           std::size_t n) {
    const auto points = independent_points(wing, record, n);
    const auto density = shared.at("density_kg_m3").get<long double>();
    const auto omega = shared.at("omega_rad_s").get<long double>();
    const auto flow = shared.at("velocity_m_s").get<std::array<long double, 3>>();
    const auto origin = shared.at("moment_origin_m").get<Position>();
    const auto& contact = take(wing.global()).get();
    const auto retained = contact.geometry();
    REQUIRE(retained);
    const auto& geometry = retained->get();
    std::array<ComplexVector, 2> force{}, moment{};
    std::vector<std::size_t> area_count(wing.cells().size()), edge_count(wing.cells().size());
    std::set<std::tuple<a::PhysicalFaceId, a::FacetId, a::PointId, std::optional<a::EdgeId>>> seen;
    const auto point = [&](a::PointId id) -> const Position& {
        return points.at(static_cast<std::size_t>(id.value() - 1));
    };
    const auto vector = [](const Position& value) {
        return std::array<long double, 3>{value[0], value[1], value[2]};
    };
    constexpr ComplexVector shift{Complex{.33L, .08L}, Complex{-.18L, .06L}, Complex{.25L, -.03L}};
    constexpr ComplexVector rotation{Complex{.02L, .1L}, Complex{.05L, -.03L},
                                     Complex{-.12L, .03L}};
    std::array<Complex, 2> work{};
    for (const auto& maybe : loads_.original.values().first(loads_.original_count)) {
        REQUIRE(maybe);
        const auto& term = *maybe;
        const auto row = static_cast<std::size_t>(term.face.value() - 1);
        const auto& cell = wing.cells()[row];
        const auto w = static_cast<std::size_t>((cell.surface - 1) / 2);
        REQUIRE(term.original_surface == cell.surface);
        REQUIRE((term.facet == cell.facets[0] || term.facet == cell.facets[1]));
        const auto canonical_facet = test_geometry::take(geometry.facet_points(term.facet));
        REQUIRE(seen.emplace(term.face, term.facet, term.point, term.canonical_edge).second);
        const auto width = n / 2 + 1;
        const auto patch = static_cast<std::size_t>(cell.surface - 1);
        const auto corner_ordinal = patch * (m + 1) * width + cell.chord * width + cell.span + 1;
        const std::array corners{a::PointId{corner_ordinal}, a::PointId{corner_ordinal + 1},
                                 a::PointId{corner_ordinal + width + 1},
                                 a::PointId{corner_ordinal + width}};
        const auto facet = term.facet == cell.facets[0]
                               ? std::array{corners[0], corners[2], corners[1]}
                               : std::array{corners[0], corners[3], corners[2]};
        const auto owner = [&](a::PointId id) {
            return wing.points()[static_cast<std::size_t>(id.value() - 1)].canonical;
        };
        CHECK(canonical_facet == std::array{owner(facet[0]), owner(facet[1]), owner(facet[2])});
        CHECK(term.canonical == owner(term.point));
        CHECK(std::ranges::find(facet, term.point) != facet.end());
        const auto gamma =
            complex_value(record.at("gamma").at(reference_index(cell, m, n)).at(drive));
        ComplexVector expected{};
        if (term.canonical_edge) {
            ++edge_count[row];
            const auto edge = contact.uvlm_load_edge(*term.canonical_edge, capacity);
            REQUIRE(edge);
            CHECK_FALSE(edge->trailing);
            const auto vertex =
                wing.points()[static_cast<std::size_t>(term.point.value() - 1)].canonical;
            CHECK((edge->endpoints[0] == vertex || edge->endpoints[1] == vertex));
            std::optional<std::size_t> side;
            for (std::size_t i = 0; i < facet.size(); ++i) {
                const auto raw = test_geometry::take(geometry.facet_edges(term.facet))[i];
                const auto bound_edge = contact.uvlm_load_edge(raw, capacity);
                REQUIRE(bound_edge);
                if (bound_edge->edge == *term.canonical_edge)
                    side = i;
            }
            REQUIRE(side);
            const auto first = facet[*side], last = facet[(*side + 1) % 3];
            const auto first_owner =
                wing.points()[static_cast<std::size_t>(first.value() - 1)].canonical;
            CHECK(term.edge_sign == (first_owner == edge->endpoints[0] ? -1 : 1));
            const auto direction = cross(flow, vector(subtract(point(first), point(last))));
            for (std::size_t i = 0; i < 3; ++i)
                expected[i] = density * gamma * direction[i] / 2.L;
        } else {
            ++area_count[row];
            CHECK(term.edge_sign == 0);
            const auto direction = cross(vector(subtract(point(facet[1]), point(facet[0]))),
                                         vector(subtract(point(facet[2]), point(facet[0]))));
            for (std::size_t i = 0; i < 3; ++i)
                expected[i] = Complex{0, density * omega} * gamma * direction[i] / 6.L;
        }
        const auto actual = loads(term.force_real, term.force_imaginary);
        long double scale = 0;
        for (std::size_t i = 0; i < actual.size(); ++i)
            scale += std::norm(complex_value(record.at("dual").at(w).at(i).at(drive)));
        REQUIRE(scale > 0);
        for (std::size_t i = 0; i < actual.size(); ++i)
            CHECK(std::abs(actual[i] - expected[i]) < 1e-9L * std::sqrt(scale));
        const auto arm = vector(subtract(point(term.point), origin));
        add(force[w], actual);
        add(moment[w], cross(arm, actual));
        const auto angular = cross(rotation, arm);
        for (std::size_t i = 0; i < actual.size(); ++i)
            work[w] += actual[i] * std::conj(shift[i] + angular[i]);
    }
    for (std::size_t i = 0; i < wing.cells().size(); ++i) {
        CHECK(area_count[i] == 6);
        CHECK(edge_count[i] == (wing.cells()[i].chord + 1 == m ? 6 : 8));
    }
    for (std::size_t w = 0; w < force.size(); ++w) {
        const auto& group = *loads_.groups.values()[w];
        compare_load(force[w], loads(group.totals.force_real, group.totals.force_imaginary));
        compare_load(moment[w], loads(group.totals.moment_real, group.totals.moment_imaginary));
        Complex virtual_work{};
        long double scale = 0;
        for (std::size_t i = 0; i < 3; ++i) {
            const auto f = complex_value(record.at("dual").at(w).at(i).at(drive));
            const auto mm = complex_value(record.at("dual").at(w).at(i + 3).at(drive));
            virtual_work += f * std::conj(shift[i]) + mm * std::conj(rotation[i]);
            scale += std::norm(f * shift[i]) + std::norm(mm * rotation[i]);
        }
        CHECK(std::abs(work[w] - virtual_work) < 1e-9L * std::sqrt(scale));
    }
    add(force[0], force[1]);
    add(moment[0], moment[1]);
    compare_load(force[0], loads(loads_.total.force_real, loads_.total.force_imaginary));
    compare_load(moment[0], loads(loads_.total.moment_real, loads_.total.moment_imaginary));
}

void verify_rotated_shared_edge(const c::CaseWingGeometry& prepared,
                                const c::DiscreteHarmonicCirculation& solved) {
    const auto& result = *solved.drive_loads.values().front();
    const auto& global = take(prepared.global()).get();
    using Key = std::pair<a::EdgeId, a::PointId>;
    struct Pair {
        std::array<ComplexVector, 2> origin;
        std::array<int, 2> sign{};
        std::array<unsigned, 2> count{};
    };
    std::map<Key, Pair> matching;
    for (const auto& maybe : result.original.values().first(result.original_count)) {
        REQUIRE(maybe);
        const auto& term = *maybe;
        if (!term.canonical_edge || term.original_surface > 2)
            continue;
        const auto point =
            prepared.points()[static_cast<std::size_t>(term.point.value() - 1)].canonical;
        auto& pair = matching[{*term.canonical_edge, point}];
        const auto side = static_cast<std::size_t>(term.original_surface - 1);
        add(pair.origin[side], loads(term.force_real, term.force_imaginary));
        pair.sign[side] += term.edge_sign;
        ++pair.count[side];
    }
    unsigned cancelling = 0;
    for (const auto& [key, pair] : matching) {
        if (pair.count != std::array{1U, 1U})
            continue;
        const auto edge = global.uvlm_load_edge(key.first, capacity);
        REQUIRE(edge);
        CHECK_FALSE(edge->trailing);
        CHECK(pair.sign[0] == -pair.sign[1]);
        const auto& first = pair.origin[0];
        const auto& second = pair.origin[1];
        const auto scale = std::abs(first[1]);
        REQUIRE(scale > 1e-3L);
        CHECK(std::abs(first[1] + second[1]) < 1e-9L * scale);
        ++cancelling;
    }
    CHECK(cancelling == 4);
}

void run_harmonic(const Json& reference, std::array<int, 3> selected, bool controls,
                  bool rotated_shared = false, HarmonicCore core = HarmonicCore::Baseline) {
    CAPTURE(selected);
    const auto found = std::ranges::find_if(
        reference.at("records"), [&](const auto& record) { return record.at("case") == selected; });
    REQUIRE(found != reference.at("records").end());
    const auto& ref = *found;
    const auto [dihedral, chord, span] = selected;
    const auto m = static_cast<std::size_t>(chord), n = static_cast<std::size_t>(span);
    auto source_bytes = frozen(dihedral, chord, span);
    if (core == HarmonicCore::Small) {
        constexpr std::string_view baseline = R"("core_radius":0.00001)";
        const auto at = source_bytes.find(baseline);
        REQUIRE(at != std::string::npos);
        REQUIRE(source_bytes.find(baseline, at + baseline.size()) == std::string::npos);
        source_bytes.replace(at, baseline.size(), R"("core_radius":1e-7)");
    }
    if (rotated_shared) {
        auto tilted = tokens::parse(source_bytes);
        for (std::size_t i = 0; i < 2; ++i) {
            auto& frame = graph(tilted, i).field("frame");
            frame.set("rotation", tokens::parse("[[0.6,0,0.8],[0,1,0],[-0.8,0,0.6]]"));
        }
        source_bytes = tilted.wire();
    }
    const auto source = read(source_bytes);
    Run run;
    auto caps = limits();
    caps.geometry.max_steps = 4'000'000;
    const c::CaseWingStraightWake terminal{{128, 0, 0}};
    const auto prepared = take(c::CaseWingGeometry::prepare(
        source, run.context, {}, projection(core::Precision::Binary64), caps, *run.geometry,
        *run.global, run.motion, run.memory, nullptr, terminal));
    REQUIRE(prepared.wake_bands().size() == 2 * n);
    const auto count = prepared.cells().size();
    REQUIRE(count == 2 * m * n);
    auto ages = ref.at("panel_starts_seconds").get<std::vector<double>>();
    // Use two Gauss16 panels per base interval. The recorded Gauss48 reference has four.
    for (std::size_t i = 1; i < ages.size() / 2; ++i)
        ages[i] = ages[2 * i];
    ages.resize(ages.size() / 2);
    const auto omega = reference.at("omega_rad_s").get<double>();
    const auto speed = reference.at("velocity_m_s").at(0).get<double>();
    std::vector<std::size_t> order;
    std::array<std::vector<c::HarmonicCirculationRhs>, 4> drive;
    std::array<long double, 4> rhs_squared{};
    for (const auto& cell : prepared.cells()) {
        order.push_back(reference_index(cell, m, n));
        const auto w = static_cast<std::size_t>((cell.surface - 1) / 2);
        const auto nz = word(cell.normal[2]);
        const auto x = word(cell.position[0]);
        const double pitch_axis = w == 0 ? .375 : 5.25;
        for (std::size_t i = 0; i < drive.size(); ++i) {
            const auto value =
                i / 2 != w ? Complex{}
                           : (i % 2 == 0 ? Complex{0., omega * nz}
                                         : Complex{-speed * nz, -omega * (x - pitch_axis) * nz});
            if (!rotated_shared) {
                const auto expected = complex_value(ref.at("rhs").at(order.back()).at(i));
                CHECK(std::abs(value - expected) < 1e-12L);
            }
            rhs_squared[i] += std::norm(value);
            drive[i].push_back({cell.face, a::SingleWord{static_cast<double>(value.real())},
                                a::SingleWord{static_cast<double>(value.imag())}});
        }
    }
    const auto& contact = take(prepared.global()).get();
    for (const auto& band : prepared.wake_bands()) {
        const auto patch = contact.uvlm_quad(band.facet, capacity);
        REQUIRE(patch);
        REQUIRE(patch->attachment);
        CHECK(patch->attachment->edge == band.trailing_edge);
        CHECK(patch->attachment->bound_face == band.bound_face);
        CHECK(patch->corners ==
              std::array{band.front[0], band.front[1], band.terminal[1], band.terminal[0]});
    }
    const std::vector<double> scales(count, 20.);
    const std::vector<ls::UnitDimensions> gamma_units(count, {2, 0, -1, 0, 0, 0, 0});
    const std::vector<ls::UnitDimensions> velocity_units(count, {1, 0, -1, 0, 0, 0, 0});
    const double min_rhs = static_cast<double>(std::sqrt(*std::ranges::min_element(rhs_squared)));
    c::DiscreteCirculationRequest solving{
        {scales, scales, gamma_units, velocity_units, harmonic_context(run.context)},
        {static_cast<core::Index>(count), 4, static_cast<core::Index>(2 * n), 5, metadata},
        {ls::Norm::L2, 1e-12, 0., 1e-12, core::ObservationId{20}, core::ObservationId{21}},
        core::RhsGeneration{3},
        core::NumericalGeneration{4},
        6e-9,
        1e-13,
        .99e-11 * min_rhs / 20.0,
        0.};
    solving.certified_closed_contour = true;
    solving.contour_workers = 2;
    const c::DiscreteNonpenetrationLimits action{
        {core::Precision::TwoWord, 1e-12},
        {core::Precision::TwoWord, 1e-12, 1e-12},
        {caps.geometry,
         {numerics::linear::EstimateAdmission{65536}, 1000000000},
         UINT64_C(150000000000)},
        UINT64_C(150000000000),
        UINT64_C(70000000000000)};
    const c::CaseWingHarmonicRequest input{{{c::CaseWingHarmonicGroup{1, 2}, {3, 4}}},
                                           {drive[0], drive[1], drive[2], drive[3]},
                                           omega,
                                           ages,
                                           16 * ages.size()};
    const auto budget = projection(core::Precision::TwoWord);
    const c::CaseWingHarmonicProjection loading{
        budget.representation, budget.exact, budget.position[0], budget.velocity, budget.position};
    const auto solve = [&](const io::SiCase& c, const c::CaseWingGeometry& geometry,
                           const c::CaseWingHarmonicRequest& request) {
        return c::solve_case_wing_harmonic(c, geometry, run.context, request, loading, caps.chart,
                                           action, solving, run.memory);
    };
    if (controls) {
        auto wrong = tokens::parse(source_bytes);
        auto& wake = wrong.field("analysis").field("wake");
        wake.set("core_radius", tokens::parse("0.00002"));
        auto newer = read(wrong.wire());
        reject(solve(newer, prepared, input), StaleGeneration);
        auto observed = input;
        observed.observations[0].moment_observation = 4;
        reject(solve(source, prepared, observed), InvalidInput);
        auto wrong_rows = drive[0];
        std::swap(wrong_rows[0].face, wrong_rows[n / 2].face);
        auto reordered = input;
        reordered.drives[0] = wrong_rows;
        reject(solve(source, prepared, reordered), InvalidInput);
        const c::CaseWingStraightWake wrong_terminal{{127, 0, 0}};
        const auto wrong_te = take(c::CaseWingGeometry::prepare(
            source, run.context, {}, projection(core::Precision::Binary64), caps, *run.geometry,
            *run.global, run.motion, run.memory, nullptr, wrong_terminal));
        reject(solve(source, wrong_te, input), StaleGeneration);
    }
    const auto result = take(solve(source, prepared, input));
    const auto& solved = result.solved;
    REQUIRE(solved.drive_loads.size() == 4);
    if (rotated_shared) {
        verify_rotated_shared_edge(prepared, solved);
        return;
    }
    REQUIRE(solved.residual.size() == count * 4);
    REQUIRE(solved.original_coefficients.size() == count * count);
    REQUIRE(solved.physical_values().shape().rows == static_cast<core::Index>(count));
    REQUIRE(solved.physical_values().shape().columns == 4);
    const auto& circulation =
        std::get<numerics::TwoWord<std::complex<double>>>(solved.physical_values().storage());
    long double dropped = 0, coefficient_scale = 0, max_row = 0;
    Json metrics{{"case", selected}, {"drives", Json::array()}};
    for (std::size_t i = 0; i < count; ++i) {
        long double diff = 0, scale = 0;
        for (std::size_t j = 0; j < count; ++j) {
            const auto expected = complex_value(ref.at("A").at(order[i]).at(order[j]));
            const auto& actual = solved.original_coefficients.values()[i * count + j];
            const auto error = [](const a::ContourInterval& bounds, long double expectation) {
                return std::max(std::abs(bounds.lower - expectation),
                                std::abs(bounds.upper - expectation));
            };
            diff += std::norm(Complex{error(actual.real, expected.real()),
                                      error(actual.imaginary, expected.imag())});
            scale += std::norm(expected);
            if (order[i] / (m * n) != order[j] / (m * n))
                dropped += std::norm(expected);
        }
        REQUIRE(scale > 0);
        const auto row_error = std::sqrt(diff / scale);
        CHECK(row_error < 1e-9L);
        max_row = std::max(max_row, row_error);
        coefficient_scale += scale;
    }
    // A solver omitting the other wing's action would fail the same matrix criterion.
    const auto omitted = std::sqrt(dropped / coefficient_scale);
    CHECK(omitted > 1e-9L);
    metrics["max_coefficient_row_relative"] = max_row;
    metrics["omitted_other_wing_relative"] = omitted;
    for (std::size_t drive_index = 0; drive_index < 4; ++drive_index) {
        CAPTURE(drive_index);
        long double diff = 0, scale = 0, residual = 0;
        for (std::size_t i = 0; i < count; ++i) {
            const auto row = drive_index * count + i;
            const auto expected = complex_value(ref.at("gamma").at(order[i]).at(drive_index));
            const Complex actual =
                Complex{circulation.hi.values()[row]} + Complex{circulation.lo.values()[row]};
            diff += std::norm(actual - expected);
            scale += std::norm(expected);
            const auto& value = solved.residual.values()[row];
            REQUIRE(value);
            CHECK(value->face == prepared.cells()[i].face);
            residual +=
                std::norm(Complex{std::abs(center(value->real)) + radius(value->real),
                                  std::abs(center(value->imaginary)) + radius(value->imaginary)});
        }
        const auto circulation_relative = std::sqrt(diff / scale);
        const auto original_relative = std::sqrt(residual / rhs_squared[drive_index]);
        CHECK(circulation_relative < 1e-9L);
        CHECK(original_relative < 1e-11L);
        Json drive_record{{"circulation_relative", circulation_relative},
                          {"all_rows_residual_relative", original_relative},
                          {"wing_loads", Json::array()}};
        const auto& load = solved.drive_loads.values()[drive_index];
        REQUIRE(load);
        REQUIRE(load->groups.size() == 2);
        for (std::size_t wing = 0; wing < 2; ++wing) {
            const auto& got = load->groups.values()[wing];
            REQUIRE(got);
            CHECK(got->id == wing);
            ComplexVector force{}, moment{};
            for (std::size_t i = 0; i < 3; ++i) {
                force[i] = complex_value(ref.at("dual").at(wing).at(i).at(drive_index));
                moment[i] = complex_value(ref.at("dual").at(wing).at(i + 3).at(drive_index));
            }
            const auto compare = [&](const auto& real, const auto& imaginary,
                                     const ComplexVector& expectation) {
                long double numerator = 0, denominator = 0;
                const auto value = loads(real, imaginary);
                for (std::size_t i = 0; i < value.size(); ++i) {
                    numerator += std::norm(Complex{
                        std::abs(value[i].real() - expectation[i].real()) + radius(real[i]),
                        std::abs(value[i].imag() - expectation[i].imag()) + radius(imaginary[i])});
                    denominator += std::norm(expectation[i]);
                }
                REQUIRE(denominator > 0);
                const auto relative = std::sqrt(numerator / denominator);
                CHECK(relative < 1e-9L);
                return relative;
            };
            drive_record["wing_loads"].push_back(
                {{"force_relative",
                  compare(got->totals.force_real, got->totals.force_imaginary, force)},
                 {"moment_relative",
                  compare(got->totals.moment_real, got->totals.moment_imaginary, moment)}});
        }
        verify_original_loads(prepared, *load, ref, reference, drive_index, m, n);
        metrics["drives"].push_back(std::move(drive_record));
    }
    std::cout << "V-MULTI " << metrics.dump() << '\n';
}
} // namespace

TEST_CASE("Typed two-wing Case solves all four finite-core harmonic drives",
          "[PHY-02][UW-D05][UW-D06][V-MULTI][verification]") {
    run_harmonic(harmonic_reference(UVLAS_VMULTI_SMALL_FIXTURE), {0, 2, 8}, true);
}

TEST_CASE("A tilted shared edge carries opposite nonzero original harmonic forces",
          "[PHY-02][UW-D05][UW-D06][V-MULTI][verification]") {
    run_harmonic(harmonic_reference(UVLAS_VMULTI_SMALL_FIXTURE), {0, 2, 8}, false, true);
}

#ifdef UVLAS_VMULTI_FORMAL_FIXTURE
#define UVLAS_VMULTI_TEST(D, M, N)                                                                 \
    TEST_CASE("Typed V-MULTI d" #D "/m" #M "/n" #N,                                                \
              "[.vmulti-harmonic][d" #D "m" #M "n" #N "][V-MULTI][verification]") {                \
        run_harmonic(harmonic_reference(UVLAS_VMULTI_FORMAL_FIXTURE), {D, M, N}, false);           \
    }
UVLAS_VMULTI_TEST(0, 4, 16)
UVLAS_VMULTI_TEST(0, 8, 32)
UVLAS_VMULTI_TEST(5, 2, 8)
UVLAS_VMULTI_TEST(5, 4, 16)
UVLAS_VMULTI_TEST(5, 8, 32)
#undef UVLAS_VMULTI_TEST
#endif

#ifdef UVLAS_VMULTI_SMALL_CORE_FIXTURES
#define UVLAS_VMULTI_SMALL_CORE_TEST(D, M, N)                                                      \
    TEST_CASE("Typed V-MULTI small core d" #D "/m" #M "/n" #N,                                     \
              "[.vmulti-small-core][d" #D "-m" #M "-n" #N "][V-MULTI][verification]") {            \
        const auto path = std::filesystem::path{UVLAS_VMULTI_SMALL_CORE_FIXTURES} /                \
                          "d" #D "-m" #M "-n" #N "-reference.json";                                \
        run_harmonic(harmonic_reference(path.string().c_str(), HarmonicCore::Small), {D, M, N},    \
                     false, false, HarmonicCore::Small);                                           \
    }
UVLAS_VMULTI_SMALL_CORE_TEST(0, 2, 8)
UVLAS_VMULTI_SMALL_CORE_TEST(0, 4, 16)
UVLAS_VMULTI_SMALL_CORE_TEST(0, 8, 32)
UVLAS_VMULTI_SMALL_CORE_TEST(5, 2, 8)
UVLAS_VMULTI_SMALL_CORE_TEST(5, 4, 16)
UVLAS_VMULTI_SMALL_CORE_TEST(5, 8, 32)
#undef UVLAS_VMULTI_SMALL_CORE_TEST
#endif
