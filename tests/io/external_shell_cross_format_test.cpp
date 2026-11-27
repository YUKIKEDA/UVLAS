#include "../structure/pg_partition_test_support.hpp"
#include "coupling/case_external_pg.hpp"
#include "coupling/case_external_pg_chart.hpp"
#include "coupling/case_nastran_static_point_support.hpp"
#include "coupling/case_pg_shell_material.hpp"
#include "reference_file_test_fixture.hpp"
#include "reference_resolution.hpp"
#include "source_digest.hpp"

#include <catch2/generators/catch_generators.hpp>

#include <source_location>
#include <stdexcept>

namespace {
namespace io = uvlas::io;
namespace core = uvlas::core;
namespace rd = uvlas::io::detail;
namespace coupling = uvlas::coupling::detail;
namespace fixture = uvlas::test::declaration_fixture;
namespace files = uvlas::test::reference_file_fixture;
using namespace uvlas::structure;

template <class T>
[[nodiscard]] const T&
required_evidence(const std::optional<T>& value,
                  std::source_location location = std::source_location::current()) {
    if (!value)
        throw std::runtime_error("Missing external chart evidence at line " +
                                 std::to_string(location.line()));
    return *value;
}

constexpr std::string_view mesh = R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 1 0
17 1000 2000 3000 1000 4000 5000 1 7 0
$EndEntities
$Nodes
1 4 10 99
2 17 0 4
10 20 40 99
1000 2000 3000
1000 4000 3000
1000 4000 5000
1000 2000 5000
$EndNodes
$Elements
1 1 701 701
2 17 3 1
701 10 20 40 99
$EndElements
)MSH";

std::string deck(bool bad) {
    return std::string{"GRID,10,,1.0,2.0,3.0\nGRID,20,,1.0,4.0,3.0\n"
                       "GRID,40,,1.0,4.0,5.0\nGRID,99,,1.0,2.0,5.0\n"
                       "CQUAD4,701,90,10,20,40,99,0.0\nMAT1,6,1.0,,0.0,1.0\n"} +
           (bad ? "PSHELL,90,6,0.5\n" : "PSHELL,90,6,0.25\n");
}

std::string hex_digest(std::string_view body) {
    constexpr std::string_view digits = "0123456789abcdef";
    std::string text;
    for (const auto byte : rd::source_digest(files::bytes_of(body)).bytes) {
        const auto value = std::to_integer<unsigned>(byte);
        text.push_back(digits[value >> 4U]);
        text.push_back(digits[value & 15U]);
    }
    return text;
}

nlohmann::json native_sidecar(bool gmsh, bool bad, std::string_view geometry) {
    const nlohmann::json zero{{"numerator", "0"}, {"denominator", "1"}};
    const nlohmann::json one{{"numerator", "1"}, {"denominator", "1"}};
    const nlohmann::json chart{
        {"orientation", "same"},
        {"vertices", nlohmann::json::array({{zero, zero}, {one, zero}, {one, one}, {zero, one}})}};
    auto result =
        nlohmann::json{{"schema", gmsh ? "uvlas.gmsh-sidecar/1" : "uvlas.nastran-sidecar/1"},
                       {"geometry", geometry},
                       {"supports", nlohmann::json::array()},
                       {"joints", nlohmann::json::array()},
                       {"shells", nlohmann::json::array()}};
    if (gmsh) {
        result["length"] = "mm";
        result["shells"].push_back(
            {{"physical", {{"dimension", 2}, {"tag", 7}}},
             {"shell", "1"},
             {"full_thickness", bad ? 500 : 250},
             {"faces", nlohmann::json::array({{{"element", "701"}, {"chart", chart}}})},
             {"material_layers", nlohmann::json::array({{{"layer", "1"}, {"material", "1"}}})}});
    } else {
        result["generator"] = "cross-format-oracle";
        result["beams"] = nlohmann::json::array();
        result["masses"] = nlohmann::json::array();
        result["units"] = {{"length", "m"}, {"mass", "kg"}, {"time", "s"}, {"angle", "deg"}};
        result["shells"].push_back({{"element", "701"}, {"shell", "1"}, {"chart", chart}});
        result["materials"] = nlohmann::json::array({{{"mid", "6"}, {"material", "1"}}});
    }
    return result;
}

void write_external(fixture::JsonToken& model, const files::Tree& tree, bool gmsh, unsigned id,
                    unsigned first_source, std::string_view native, std::string_view sidecar) {
    const auto geometry = std::to_string(id);
    const auto native_file = geometry + (gmsh ? ".msh" : ".bdf");
    const auto sidecar_file = geometry + ".json";
    auto external = fixture::fragment(gmsh ? "Geometry.gmsh" : "Geometry.nastran");
    external.set("id", fixture::string(geometry));
    if (gmsh) {
        external.set("reference", fixture::fragment("ReferenceGeometry.fit"));
    }
    auto& input = external.field(gmsh ? "mesh" : "deck");
    input.set("source_id", fixture::string(std::to_string(first_source)));
    files::declare_path(input, native_file, hex_digest(native));
    auto& map = external.field("sidecar");
    map.set("source_id", fixture::string(std::to_string(first_source + 1)));
    files::declare_path(map, sidecar_file, hex_digest(sidecar));
    files::write_file(tree.root / native_file, native);
    files::write_file(tree.root / sidecar_file, sidecar);
    model.field("geometry").elements.push_back(std::move(external));
}

void add_external(fixture::JsonToken& model, const files::Tree& tree, bool gmsh, bool bad,
                  unsigned id, unsigned first_source) {
    write_external(model, tree, gmsh, id, first_source, gmsh ? std::string{mesh} : deck(bad),
                   native_sidecar(gmsh, bad, std::to_string(id)).dump());
}

rd::MeshSidecarBindingLimits make_binding_limits() {
    return {{files::limits().document,
             {256, 256, 256, 8, 8, 4096, 8192, ExactEstimateAdmission{fixture::budget},
              ExactTimingMode::Disabled},
             32,
             8,
             32'000},
            {4096, 16, 16, 64, 16, 16, 32'000},
            {4096, 80, 32, 32, 256, 32'000},
            {32, 32, 8, 32'000},
            {32, 32, 32'000},
            {32, 32, 32'000},
            8,
            32'000};
}

rd::ReferencedCase read_graph(const files::Tree& tree, const fixture::JsonToken& config,
                              core::ResourcePlan& plan) {
    const auto binding_limits = make_binding_limits();
    files::write_file(tree.root / "case.json", config.wire());
    rd::ReferenceInputDiagnostic input_problem;
    const auto input = rd::ReferenceInput::from_file(tree.root, "case.json", io::SourceId{12},
                                                     files::manifest_limits, files::limits(), plan,
                                                     {}, &input_problem);
    CAPTURE(input_problem.stage, input_problem.schema.path());
    REQUIRE(input);
    auto declarations = files::prepare(*input);
    rd::LeafReferenceDiagnostic reference_problem;
    auto result = rd::resolve_leaf_references(
        *input, std::move(declarations),
        {{{4096, 4097, 2, 4096}, 64, 256}, 8, 32'768, 32'776, 16, 32'768, 64}, {8, 8, 1, 512, 512},
        {binding_limits.gmsh, binding_limits.nastran, binding_limits.sidecar}, plan, {}, nullptr,
        &reference_problem);
    CAPTURE(reference_problem.reason, reference_problem.pointer());
    if (const auto& sidecar_problem = reference_problem.sidecar) {
        INFO(static_cast<unsigned>(sidecar_problem->reason));
        REQUIRE(result);
    }
    return fixture::take(std::move(result));
}

io::SiCase convert(const io::SchemaCase& schema, const ExactLimits& exact,
                   core::ResourcePlan& plan) {
    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> caps{};
    caps.fill(fixture::budget);
    const io::SiLimits limits{caps, exact};
    const auto catalog = fixture::take(io::make_si_catalog(schema, limits, plan, nullptr));
    std::vector<io::SiRequirementRow> requirements;
    for (std::uint64_t i = 0; i < fixture::take(catalog.size()); ++i) {
        const auto recipe = fixture::take(catalog.entry(i));
        requirements.push_back({recipe.key,
                                {io::SiDecimalInput{"1"}, recipe.recipe.dimensions},
                                io::SiDecimalInput{"1e-25"},
                                io::SiDecimalInput{"1e-25"}});
    }
    const auto request =
        fixture::take(io::make_si_request(catalog, requirements, limits, plan, nullptr));
    return fixture::take(io::convert_case_to_si(schema, catalog, request, limits, plan, nullptr));
}

rd::NastranLaminateCards quantities(const rd::ReferencedCase& graph,
                                    const io::SchemaIdentity& parent, io::SourceId source,
                                    const rd::MeshSidecarBindingLimits& limits,
                                    core::ResourcePlan& plan) {
    auto input = fixture::take(
        rd::prepare_referenced_nastran_deck(graph, parent, source, limits.nastran, plan, {}));
    auto coordinates = fixture::take(
        rd::NastranCoordinateCards::prepare(std::move(input), limits.coordinates, plan, {}));
    auto shells = fixture::take(
        rd::NastranShellCards::prepare(std::move(coordinates), limits.shells, plan, {}));
    auto mats = fixture::take(
        rd::NastranIsotropicShellCards::prepare(std::move(shells), {32, 32'000}, plan, {}));
    return fixture::take(
        rd::NastranLaminateCards::prepare(std::move(mats), {32, 32'000}, plan, {}));
}

coupling::CasePgChartLimits chart_preparation(const coupling::CasePgReferenceLimits& reference,
                                              std::uint64_t patches, std::uint64_t sources) {
    auto exact = reference.geometry;
    exact.algebra.maximum_value_bits = 65536;
    exact.algebra.maximum_work_bits = 262144;
    exact.maximum_scalar_operations = 300'000'000;
    exact.maximum_partition_visits = 100'000'000;
    constexpr auto steps = std::uint64_t{10'000'000'000'000};
    constexpr auto profile = std::uint64_t{48} << 30U;
    return {
        reference,
        {profile, steps, 128, {}},
        {768, 704},
        {exact, steps},
        {768, 192},
        {profile, steps, {}},
        exact,
        {reference.geometry.algebra, 10'000, coupling::case_motion_frame_interval_scratch_bytes},
        {exact, 120'000},
        {exact, patches, sources, 500, 120'000},
        {exact, 512, 960, 140'000},
        100'000};
}

class RectangleOracle {
public:
    RectangleOracle(const PgReferenceLimits& limits, core::ResourcePlan& plan)
        : limits_(limits), plan_(plan) {}

    [[nodiscard]] ExactRational r(std::string_view n, std::string_view d = "1") const {
        return fixture::take(ExactRational::from_integer_pair(n, d, limits_.algebra, plan_));
    }

    template <class Face>
    void corners(const Face& face) const {
        constexpr std::array<std::uint64_t, 4> ids{10, 20, 40, 99};
        constexpr std::array<std::string_view, 4> chart_s{"0", "1", "1", "0"};
        constexpr std::array<std::string_view, 4> chart_t{"0", "0", "1", "1"};
        constexpr std::array<std::string_view, 4> world_y{"2", "4", "4", "2"};
        constexpr std::array<std::string_view, 4> world_z{"3", "3", "5", "5"};
        CHECK(face.element == 701);
        CHECK(face.winding == rd::ExternalPgWinding::Forward);
        const auto ring = face.ring();
        REQUIRE(ring.size() == ids.size());
        for (std::size_t i = 0; i < ids.size(); ++i) {
            CHECK(ring[i].node == ids[i]);
            CHECK(ring[i].chart.s == r(chart_s[i]));
            CHECK(ring[i].chart.t == r(chart_t[i]));
            CHECK(ring[i].position_metres[0] == r("1"));
            CHECK(ring[i].position_metres[1] == r(world_y[i]));
            CHECK(ring[i].position_metres[2] == r(world_z[i]));
        }
    }

    void thickness(const rd::ExternalPgThicknessReport& report,
                   const PgReferenceIdentity& identity) const {
        CHECK(report.pg == identity);
        REQUIRE(report.matches.size() == 1);
        const auto& face = report.matches.values().front();
        CHECK(face.element == 701);
        CHECK(std::get<ExactRational>(face.value_metres) == r("1", "4"));
        CHECK(face.twice_chart_area == r("2"));
        CHECK(face.pieces > 0);
    }

    void physical(const PgReferenceGeometry& reference) const {
        const auto& partition = fixture::take(reference.source()).partition;
        const auto volumes = fixture::take(partition.volumes());
        const auto horizontal = fixture::take(partition.horizontal_faces());
        const auto triangles = fixture::take(partition.triangles());
        const auto vertices = fixture::take(partition.vertices());
        REQUIRE_FALSE(volumes.empty());
        REQUIRE(horizontal.size() == 2 * volumes.size());
        const auto quarter = r("1", "4");
        const auto depth = r("3", "4");
        const auto two = r("2");
        const auto three = r("3");
        auto triangle = [&](core::Index index) {
            const auto& row = triangles[static_cast<std::size_t>(index)];
            const auto& a = vertices[static_cast<std::size_t>(row.vertices[0])];
            const auto& b = vertices[static_cast<std::size_t>(row.vertices[1])];
            const auto& c = vertices[static_cast<std::size_t>(row.vertices[2])];
            return std::array{quarter_point(a.s, b.s, c.s), quarter_point(a.t, b.t, c.t),
                              minus(times(minus(b.s, a.s), minus(c.t, a.t)),
                                    times(minus(b.t, a.t), minus(c.s, a.s)))};
        };
        ExactRational total_volume;
        ExactRational total_upper_area;
        for (std::size_t i = 0; i < volumes.size(); ++i) {
            CAPTURE(i);
            const auto [s, t, area] = triangle(volumes[i].triangle);
            const auto point = fixture::take(reference.evaluate(
                PgReferenceVolumePoint{static_cast<core::Index>(i), quarter, quarter, depth},
                std::nullopt, limits_, plan_));
            CHECK(point.layer == core::LayerId{1});
            CHECK(point.full_thickness_metres == quarter);
            CHECK(point.g == r("1"));
            CHECK(point.s == s);
            CHECK(point.t == t);
            CHECK(point.height_chart_jet[0] == r("0"));
            expect(point.normal[0], r("1"));
            expect(point.normal[1], r("0"));
            expect(point.normal[2], r("0"));
            expect(point.position_metres[0], r("17", "16"));
            expect(point.position_metres[1], sum(two, times(s, two)));
            expect(point.position_metres[2], sum(three, times(t, two)));
            expect(point.determinant_metres_cubed, r("1", "2"));
            const auto& measure = std::get<PgReferenceVolumeMeasure>(point.measure).jacobian;
            expect(measure, area);
            total_volume = sum(total_volume, half(sum(measure.rational, measure.radical)));
        }
        for (std::size_t i = 0; i < horizontal.size(); ++i) {
            CAPTURE(i);
            const auto [s, t, area] = triangle(horizontal[i].triangle);
            const auto point = fixture::take(reference.evaluate(
                PgReferenceHorizontalPoint{static_cast<core::Index>(i), quarter, quarter},
                std::nullopt, limits_, plan_));
            const auto& measure = std::get<PgReferenceHorizontalMeasure>(point.measure);
            const bool upper = horizontal[i].side == PgPartitionFaceSide::Upper;
            const auto area_scale = times(area, r("4"));
            CHECK(measure.side == horizontal[i].side);
            CHECK(point.g == r("1"));
            CHECK(point.s == s);
            CHECK(point.t == t);
            expect(point.position_metres[0], upper ? r("9", "8") : r("7", "8"));
            expect(point.position_metres[1], sum(two, times(s, two)));
            expect(point.position_metres[2], sum(three, times(t, two)));
            expect(measure.area.squared_norm, times(area_scale, area_scale));
            expect(measure.area.outward[0],
                   upper ? times(area_scale, r("1")) : times(area_scale, r("-1")));
            expect(measure.area.outward[1], r("0"));
            expect(measure.area.outward[2], r("0"));
            if (upper) {
                total_upper_area = sum(
                    total_upper_area,
                    half(sum(measure.area.outward[0].rational, measure.area.outward[0].radical)));
            }
        }
        CHECK(total_volume == r("1"));
        CHECK(total_upper_area == r("4"));
    }

private:
    void expect(const PgReferenceScalar& actual, const ExactRational& expected,
                const std::source_location where = std::source_location::current()) const {
        INFO("Formula at line " << where.line());
        // The flat rectangle has g=1; the representation may distribute a rational value
        // between the two coefficients.
        CHECK(sum(actual.rational, actual.radical) == expected);
    }
    [[nodiscard]] ExactRational sum(const ExactRational& a, const ExactRational& b) const {
        return fixture::take(a.add(b, limits_.algebra, plan_));
    }
    [[nodiscard]] ExactRational half(const ExactRational& x) const {
        return fixture::take(x.divide(r("2"), limits_.algebra, plan_));
    }
    [[nodiscard]] ExactRational minus(const ExactRational& a, const ExactRational& b) const {
        return fixture::take(a.subtract(b, limits_.algebra, plan_));
    }
    [[nodiscard]] ExactRational times(const ExactRational& a, const ExactRational& b) const {
        return fixture::take(a.multiply(b, limits_.algebra, plan_));
    }
    [[nodiscard]] ExactRational quarter_point(const ExactRational& a, const ExactRational& b,
                                              const ExactRational& c) const {
        return times(sum(sum(times(a, r("2")), b), c), r("1", "4"));
    }
    const PgReferenceLimits& limits_;
    core::ResourcePlan& plan_;
};
} // namespace

TEST_CASE("One planar Case has the same measured Shell through Nastran and Gmsh",
          "[FD-03][IO4-01][IO4-03][IO4-04][R-05][V-USE][integration][verification][negative]") {
    core::ResourcePlan plan{fixture::budget};
    auto config = fixture::root("coupled_transient");
    auto& model = config.field("model");
    auto& pg = model.field("geometry").elements.front();
    pg.field("frame").set("origin", fixture::json("[1,2,3]"));
    pg.field("frame").set("rotation", fixture::json("[[0,0,1],[1,0,0],[0,1,0]]"));
    pg.set("reference_length", fixture::json("2"));
    pg.field("thickness").set("values", fixture::json("[0.25,0.25,0.25,0.25]"));
    auto& layer = model.field("shells").elements.front().field("layers").elements.front();
    layer.field("axis").field("direction").set("value", fixture::json("[0,1,0]"));
    const auto binding = make_binding_limits();
    const auto graph = [&] {
        const files::Tree tree;
        add_external(model, tree, true, false, 2, 1);
        add_external(model, tree, false, false, 3, 3);
        add_external(model, tree, true, true, 4, 5);
        add_external(model, tree, false, true, 5, 7);
        return read_graph(tree, config, plan);
    }();
    const auto schema = fixture::checked(config);
    const auto parent = fixture::take(schema.identity());
    const auto graph_identity = fixture::take(graph.identity());
    CHECK(graph_identity.source == parent.source);
    CHECK(graph_identity.contract_revision == parent.contract_revision);
    CHECK(graph_identity.representation_revision == parent.representation_revision);
    const PgReferenceLimits pg_limits{partition_test::exact_limits(), 20'000'000, 10'000'000, {}};
    const auto si = convert(schema, pg_limits.algebra, plan);
    const auto bind = [&](unsigned id) {
        return fixture::take(rd::bind_external_mesh_sidecar(schema, graph, rd::MeshGeometryId{id},
                                                            binding, plan, {}, nullptr));
    };
    const auto gmsh = bind(2);
    const auto nastran = bind(3);
    const auto bad_gmsh = bind(4);
    const auto bad_nastran = bind(5);
    const auto native = quantities(graph, parent, io::SourceId{3}, binding, plan);
    const auto bad_native = quantities(graph, parent, io::SourceId{7}, binding, plan);

    const coupling::CasePgReferenceLimits preparation{
        {pg_limits.algebra, 4, {}},
        {pg_limits.algebra, 2000, 20'000'000, 100'000, {}},
        partition_test::limits(),
        pg_limits,
        100'000,
        {}};
    const rd::ExternalSidecarPgLimits qualification{
        {pg_limits, 8}, {8, 32, 32, 8, 32'000}, 8, 32'000};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    const auto prepare_gmsh = [&](const auto& corresponding, auto* problem) {
        return coupling::CaseGmshPgSet::prepare(si, corresponding, preparation, qualification,
                                                {1, 1, 200'000}, workspace, plan, problem);
    };
    const auto prepare_nastran = [&](const auto& corresponding, const auto& cards, auto* problem) {
        return coupling::CaseNastranPlanarPgSet::prepare(si, corresponding, cards, preparation,
                                                         qualification, {1, 1, 1}, workspace, plan,
                                                         problem);
    };

    coupling::CaseGmshPgSetDiagnostic gmsh_diag;
    coupling::CaseNastranPgSetDiagnostic nastran_diag;
    const auto g_set = fixture::take(prepare_gmsh(gmsh, &gmsh_diag));
    const auto n_set = fixture::take(prepare_nastran(nastran, native, &nastran_diag));
    REQUIRE(g_set.size() == 1);
    REQUIRE(n_set.size() == 1);
    const auto& g = fixture::take(g_set.at(0)).get();
    const auto& n = fixture::take(n_set.at(0)).get();
    const auto identity = fixture::take(g.reference().identity());
    CHECK(fixture::take(n.reference().identity()) == identity);
    const RectangleOracle oracle{pg_limits, plan};

    const auto& g_match = g.match();
    const auto& n_match = n.match().planar;
    CHECK(g_match.shell == rd::MeshShellId{1});
    CHECK(n_match.shell == rd::MeshShellId{1});
    CHECK(g_match.pg == rd::MeshGeometryId{1});
    CHECK(n_match.pg == rd::MeshGeometryId{1});
    REQUIRE(g_match.membership.size() == 1);
    CHECK(g_match.membership.values().front().physical == rd::GmshKey{2, 7});
    CHECK(g_match.membership.values().front().element == 701);
    const auto& fit = std::get<rd::GmshPgConversionEvidence>(g_match.geometry);
    if (const auto gv = fit.geometry.view()) {
        CHECK(gv->source.parent.source == parent.source);
        CHECK(gv->source.source == io::SourceId{1});
        REQUIRE(gv->faces.size() == 1);
        oracle.corners(gv->faces.front());
        CHECK(gv->faces.front().distance_upper_metres.sign() == 0);
        CHECK(gv->faces.front().normal_chord_upper.sign() == 0);
    } else {
        FAIL("Missing Gmsh bounded face");
    }
    if (const auto nv = n_match.geometry.view()) {
        CHECK(nv->source.parent.source == parent.source);
        CHECK(nv->source.source == io::SourceId{3});
        REQUIRE(nv->faces.size() == 1);
        oracle.corners(nv->faces.front());
    } else {
        FAIL("Missing Nastran exact face");
    }

    if (const auto& gt = g_match.thickness) {
        oracle.thickness(*gt, identity);
    } else {
        FAIL("Missing Gmsh thickness comparison");
    }
    if (const auto& nt = n_match.thickness_regions) {
        oracle.thickness(*nt, identity);
        CHECK(nt->matches.values().front().property == rd::NastranPropertyId{90});
        CHECK(nt->matches.values().front().field);
    } else {
        FAIL("Missing Nastran thickness comparison");
    }
    oracle.physical(g.reference());
    oracle.physical(n.reference());

    const auto g_failure = prepare_gmsh(bad_gmsh, &gmsh_diag);
    const auto n_failure = prepare_nastran(bad_nastran, bad_native, &nastran_diag);
    REQUIRE_FALSE(g_failure);
    REQUIRE_FALSE(n_failure);
    CHECK(g_failure.error().code == core::FailureCode::InvalidInput);
    CHECK(n_failure.error().code == core::FailureCode::InvalidInput);
    CHECK(gmsh_diag.shell == rd::MeshShellId{1});
    CHECK(nastran_diag.shell == rd::MeshShellId{1});
    if (const auto& problem = gmsh_diag.member.regional_thickness) {
        CHECK(problem->reason == rd::ExternalPgThicknessFailure::Difference);
        CHECK(problem->element == 701);
        const auto& source = std::get<rd::GmshSidecar>(bad_gmsh.sidecar().payload());
        REQUIRE(source.shells.size() == 1);
        const auto& group = source.shells.values().front();
        CHECK(group.physical == rd::GmshKey{2, 7});
        REQUIRE(group.full_thickness_node);
        CHECK(problem->sidecar_node == group.full_thickness_node);
        CHECK(problem->case_triangle);
    } else {
        FAIL("Missing originating physical-group thickness diagnostic");
    }
    if (const auto& problem = nastran_diag.member.planar.regional_thickness) {
        CHECK(problem->reason == rd::ExternalPgThicknessFailure::Difference);
        CHECK(problem->element == 701);
        CHECK(problem->property == rd::NastranPropertyId{90});
        if (const auto& position = problem->field) {
            CHECK(position->byte == deck(true).find("0.5"));
            CHECK(position->line == 7);
            CHECK(position->column == 13);
        } else {
            FAIL("Missing PSHELL T position");
        }
        CHECK(problem->case_triangle);
    } else {
        FAIL("Missing originating deck thickness diagnostic");
    }
    const auto g_again = fixture::take(prepare_gmsh(gmsh, &gmsh_diag));
    const auto n_again = fixture::take(prepare_nastran(nastran, native, &nastran_diag));
    CHECK(fixture::take(fixture::take(g_again.at(0)).get().reference().identity()) == identity);
    CHECK(fixture::take(fixture::take(n_again.at(0)).get().reference().identity()) == identity);
}

TEST_CASE(
    "Authenticated external reference owners enter the same full-depth Case chart",
    "[.external-case-chart][FD-03][IO4-03][IO4-04][R-05][V-USE][V-PG][integration][negative]") {
    using Code = core::FailureCode;
    core::ResourcePlan plan{fixture::budget};
    auto config = fixture::root("coupled_transient");
    auto& model = config.field("model");
    model.field("supports").elements[0] = fixture::fragment("Support.trace_constraint");
    auto& pg = model.field("geometry").elements.front();
    pg.field("frame").set("origin", fixture::json("[1,2,3]"));
    pg.field("frame").set("rotation", fixture::json("[[0,0,1],[1,0,0],[0,1,0]]"));
    pg.set("reference_length", fixture::json("2"));
    pg.field("thickness").set("values", fixture::json("[0.25,0.25,0.25,0.25]"));
    model.field("shells")
        .elements.front()
        .field("layers")
        .elements.front()
        .field("axis")
        .field("direction")
        .set("value", fixture::json("[0,1,0]"));
    const auto with_line = [](std::string_view line) {
        std::string msh{mesh};
        const auto entities = msh.find("0 0 1 0\n");
        REQUIRE(entities != std::string::npos);
        msh.replace(entities, std::string_view{"0 0 1 0\n"}.size(),
                    "0 1 1 0\n11 1000 2000 3000 1000 4000 5000 1 7 0\n");
        constexpr std::string_view old = "1 1 701 701\n";
        const auto elements = msh.find(old);
        REQUIRE(elements != std::string::npos);
        msh.replace(elements, old.size(), "2 2 501 701\n1 11 1 1\n501 " + std::string{line} + "\n");
        return msh;
    };
    const auto map = [](unsigned id) {
        auto sidecar = native_sidecar(true, false, std::to_string(id));
        sidecar["supports"].push_back(
            {{"physical", {{"dimension", 1}, {"tag", 7}}}, {"support", "1"}});
        return sidecar.dump();
    };
    const auto graph = [&] {
        const files::Tree tree;
        write_external(model, tree, true, 2, 1, with_line("99 10"), map(2));
        add_external(model, tree, false, false, 3, 3);
        write_external(model, tree, true, 4, 5, with_line("20 40"), map(4));
        add_external(model, tree, false, true, 5, 7);
        return read_graph(tree, config, plan);
    }();
    const auto schema = fixture::checked(config);
    const auto parent = fixture::take(schema.identity());
    const auto binding_limits = make_binding_limits();
    const auto bind = [&](unsigned id) {
        return fixture::take(rd::bind_external_mesh_sidecar(schema, graph, rd::MeshGeometryId{id},
                                                            binding_limits, plan, {}, nullptr));
    };
    const auto g = bind(2), n = bind(3), wrong_edge = bind(4), wrong_thickness = bind(5);
    const auto cards = quantities(graph, parent, io::SourceId{3}, binding_limits, plan);
    const auto bad_cards = quantities(graph, parent, io::SourceId{7}, binding_limits, plan);

    const PgReferenceLimits arithmetic{partition_test::exact_limits(), 20'000'000, 10'000'000, {}};
    const auto si = convert(schema, arithmetic.algebra, plan);
    const coupling::CasePgReferenceLimits preparation{
        {arithmetic.algebra, 4, {}},
        {arithmetic.algebra, 2000, 20'000'000, 100'000, {}},
        partition_test::limits(),
        arithmetic,
        100'000,
        {}};
    const rd::ExternalSidecarPgLimits correspondence{
        {arithmetic, 8}, {8, 32, 32, 8, 32'000}, 8, 32'000};
    const auto chart_limits = chart_preparation(preparation, 1, 1);
    const coupling::CaseGmshPgChartLimits g_limits{
        chart_limits, correspondence, {1, 1, 200'000}, {arithmetic, 1}, 1};
    const coupling::CaseNastranPgChartLimits n_limits{chart_limits, correspondence, {1, 1, 1}};
    const RectangleOracle oracle{arithmetic, plan};
    const auto plane_weight = oracle.r("1", "6");
    const auto depth_weight = oracle.r("1", "3");
    std::array<ExactRational, 6> positions;
    std::vector<PgLineQuadratureNode> plane;
    plane.reserve(positions.size());
    for (std::size_t i = 0; i < positions.size(); ++i) {
        positions[i] = oracle.r(std::to_string(2 * i + 1), "12");
        plane.push_back({positions[i], plane_weight});
    }
    const auto quarter = oracle.r("1", "4"), half = oracle.r("1", "2"), last = oracle.r("3", "4");
    const std::array depth{PgLineQuadratureNode{quarter, depth_weight},
                           PgLineQuadratureNode{half, depth_weight},
                           PgLineQuadratureNode{last, depth_weight}};
    const std::array rules{
        coupling::CasePgRuleBinding{1, PgVolumeProductRule{plane, plane, depth}}};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    coupling::CaseGmshPgChartDiagnostic g_problem;
    coupling::CaseNastranPgSetDiagnostic n_problem;
    core::ResourcePlan none{0};
    const auto stopped = coupling::prepare_case_gmsh_pg_chart(
        si, g, rules, g_limits, workspace, none, partition_test::requested_stop_token());
    REQUIRE_FALSE(stopped);
    CHECK(stopped.error().code == Code::Cancelled);
    auto short_limits = g_limits;
    short_limits.maximum_support_mappings = 0;
    const auto denied =
        coupling::prepare_case_gmsh_pg_chart(si, g, rules, short_limits, workspace, none);
    REQUIRE_FALSE(denied);
    CHECK(denied.error().code == Code::ResourceExceeded);
    CHECK(none.allocated_bytes() == 0);

    const auto rejected_edge = coupling::prepare_case_gmsh_pg_chart(
        si, wrong_edge, rules, g_limits, workspace, plan, {}, &g_problem);
    REQUIRE_FALSE(rejected_edge);
    const auto& edge_problem = required_evidence(g_problem.support);
    CHECK(rejected_edge.error().code == Code::InvalidInput);
    CHECK(edge_problem.reason == coupling::CaseGmshSupportBoundaryFailure::OutsideBoundary);
    CHECK(edge_problem.line == 501);
    const auto rejected_thickness = coupling::prepare_case_nastran_pg_chart(
        si, wrong_thickness, bad_cards, rules, n_limits, workspace, plan, {}, &n_problem);
    REQUIRE_FALSE(rejected_thickness);
    CHECK(rejected_thickness.error().code == Code::InvalidInput);
    CHECK(required_evidence(n_problem.member.planar.regional_thickness).element == 701);

    core::ResourcePlan g_plan{fixture::budget}, n_plan{fixture::budget},
        direct_plan{fixture::budget};
    const auto g_chart = fixture::take(coupling::prepare_case_gmsh_pg_chart(
        si, g, rules, g_limits, workspace, g_plan, {}, &g_problem));
    const auto n_chart = fixture::take(coupling::prepare_case_nastran_pg_chart(
        si, n, cards, rules, n_limits, workspace, n_plan, {}, &n_problem));
    const auto reference = fixture::take(
        coupling::prepare_case_pg_chart(si, rules, chart_limits, workspace, direct_plan));
    core::ResourcePlan g_separate{fixture::budget}, n_separate{fixture::budget};
    const auto g_set = fixture::take(coupling::CaseGmshPgSet::prepare(
        si, g, preparation, correspondence, g_limits.shells, workspace, g_separate));
    REQUIRE(coupling::certify_case_gmsh_support_boundary(
        si, g, g_set, rd::MeshSupportId{1}, rd::GmshKey{1, 7}, g_limits.support, g_separate));
    REQUIRE(coupling::CaseNastranPlanarPgSet::prepare(si, n, cards, preparation, correspondence,
                                                      n_limits.shells, workspace, n_separate));
    CHECK(g_plan.allocated_bytes() < g_separate.allocated_bytes() + direct_plan.allocated_bytes());
    CHECK(n_plan.allocated_bytes() < n_separate.allocated_bytes() + direct_plan.allocated_bytes());
    const auto expected = fixture::take(reference.shell_reference(uvlas::numerics::PatchId{1}, 1));
    const auto g_shell =
        fixture::take(g_chart.chart.shell_reference(uvlas::numerics::PatchId{1}, 1));
    const auto n_shell =
        fixture::take(n_chart.chart.shell_reference(uvlas::numerics::PatchId{1}, 1));
    REQUIRE(g_chart.shells.size() == 1);
    const auto& g_evidence = required_evidence(g_chart.shells.values().front());
    REQUIRE(n_chart.shells.size() == 1);
    const auto& n_evidence = required_evidence(n_chart.shells.values().front());
    REQUIRE(g_chart.supports.size() == 1);
    const auto& support = required_evidence(g_chart.supports.values().front());
    REQUIRE(support.lines.size() == 1);
    CHECK(support.lines.values().front().nodes == std::array<std::uint64_t, 2>{99, 10});
    CHECK(support.lines.values().front().incident_face == 701);
    const auto identity = fixture::take(
        fixture::take(expected.patch().conditioning().reference_basis())->reference_identity());
    CHECK(support.patch == identity);
    CHECK(required_evidence(g_evidence.thickness).pg == identity);
    CHECK(required_evidence(n_evidence.planar.thickness_regions).pg == identity);
    for (const auto* observed : {&g_shell, &n_shell}) {
        CHECK(observed->chart_identity() == expected.chart_identity());
        CHECK(observed->columns().basis == expected.columns().basis);
        CHECK(observed->fixed_rule().size() == expected.fixed_rule().size());
        const auto& original = observed->chart().original();
        REQUIRE(original.sources().size() == 1);
        CHECK(original.sources()[0].kind == PgTraceKind::Support);
        REQUIRE_FALSE(original.original_rows().empty());
        const auto x0 = observed->patch().reference_coordinates();
        const auto residual = fixture::take(original.residual(
            {original.layout(), x0}, chart_limits.chart.arithmetic, workspace, plan));
        REQUIRE(residual.values().size() == original.original_rows().size());
        for (const auto& value : residual.values())
            for (const auto& component : value)
                CHECK(component.sign() == 0);
    }
}

TEST_CASE("Both Gmsh joint sides retain provenance in the assembled Case connection",
          "[.external-case-chart][FD-03][FD-05][IO4-04][R-05][V-PG][integration][negative]") {
    constexpr std::string_view joined_mesh = R"msh($MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 2 2 0
11 1000 0 0 1000 1000 0 1 11 0
12 1000 0 0 1000 1000 0 1 12 0
17 0 0 0 1000 1000 0 1 7 0
18 1000 0 0 2000 1000 0 1 7 0
$EndEntities
$Nodes
2 8 10 80
2 17 0 4
10 20 30 40
0 0 0
1000 0 0
1000 1000 0
0 1000 0
2 18 0 4
50 60 70 80
1000 0 0
2000 0 0
2000 1000 0
1000 1000 0
$EndNodes
$Elements
4 4 501 702
1 11 1 1
501 20 30
1 12 1 1
502 80 50
2 17 3 1
701 10 20 30 40
2 18 3 1
702 50 60 70 80
$EndElements
)msh";
    constexpr std::string_view adhesive = R"json(
{"kind":"adhesive","id":"1",
 "first":{"target":{"kind":"shell","id":"1"},"kind":"edge","geometry":"1","boundary":1,"side":"interior"},
 "second":{"target":{"kind":"shell","id":"2"},"kind":"edge","geometry":"6","boundary":0,"side":"interior"},
 "mapping":[
   {"first_cell":[1,0],"second_cell":[0,0],
    "domain":[[{"numerator":"0","denominator":"1"},{"numerator":"-1","denominator":"1"}],
              [{"numerator":"1","denominator":"2"},{"numerator":"-1","denominator":"1"}],
              [{"numerator":"1","denominator":"2"},{"numerator":"1","denominator":"1"}],
              [{"numerator":"0","denominator":"1"},{"numerator":"1","denominator":"1"}]],
    "first":[[{"powers":[0,0],"coefficient":{"numerator":"1","denominator":"1"}}],
             [{"powers":[1,0],"coefficient":{"numerator":"1","denominator":"1"}}],
             [{"powers":[0,1],"coefficient":{"numerator":"1","denominator":"1"}}]],
    "second":[[],[{"powers":[1,0],"coefficient":{"numerator":"1","denominator":"1"}}],
                 [{"powers":[0,1],"coefficient":{"numerator":"1","denominator":"1"}}]],
    "orientation":"same"},
   {"first_cell":[1,1],"second_cell":[0,1],
    "domain":[[{"numerator":"1","denominator":"2"},{"numerator":"-1","denominator":"1"}],
              [{"numerator":"1","denominator":"1"},{"numerator":"-1","denominator":"1"}],
              [{"numerator":"1","denominator":"1"},{"numerator":"1","denominator":"1"}],
              [{"numerator":"1","denominator":"2"},{"numerator":"1","denominator":"1"}]],
    "first":[[{"powers":[0,0],"coefficient":{"numerator":"1","denominator":"1"}}],
             [{"powers":[1,0],"coefficient":{"numerator":"1","denominator":"1"}}],
             [{"powers":[0,1],"coefficient":{"numerator":"1","denominator":"1"}}]],
    "second":[[],[{"powers":[1,0],"coefficient":{"numerator":"1","denominator":"1"}}],
                 [{"powers":[0,1],"coefficient":{"numerator":"1","denominator":"1"}}]],
    "orientation":"same"}],"split_at_attachment":false}
)json";
    core::ResourcePlan plan{fixture::budget};
    auto config = fixture::root("coupled_transient");
    auto& model = config.field("model");
    model.field("supports").elements[0] = fixture::fragment("Support.trace_constraint");
    model.field("joints").elements.push_back(fixture::json(adhesive));
    auto& pg = model.field("geometry").elements.front();
    pg.field("thickness").set("values", fixture::json("[0.25,0.25,0.25,0.25]"));
    auto other_geometry = pg;
    other_geometry.set("id", fixture::string("6"));
    other_geometry.field("frame").set("origin", fixture::json("[1,0,0]"));
    model.field("geometry").elements.push_back(std::move(other_geometry));
    auto other_shell = model.field("shells").elements.front();
    other_shell.set("id", fixture::string("2"));
    other_shell.set("geometry", fixture::string("6"));
    other_shell.field("layers").elements.front().set("id", fixture::string("2"));
    model.field("shells").elements.push_back(std::move(other_shell));

    const auto mapping = [](unsigned geometry, bool swapped) {
        auto result = native_sidecar(true, false, std::to_string(geometry));
        auto right = result["shells"][0];
        right["shell"] = "2";
        right["faces"][0]["element"] = "702";
        right["material_layers"][0]["layer"] = "2";
        result["shells"].push_back(std::move(right));
        result["joints"].push_back({{"first", {{"dimension", 1}, {"tag", swapped ? 12 : 11}}},
                                    {"second", {{"dimension", 1}, {"tag", swapped ? 11 : 12}}},
                                    {"joint", "1"}});
        return result.dump();
    };
    const auto graph = [&] {
        const files::Tree tree;
        write_external(model, tree, true, 2, 1, joined_mesh, mapping(2, false));
        write_external(model, tree, true, 3, 3, joined_mesh, mapping(3, true));
        return read_graph(tree, config, plan);
    }();
    const auto schema = fixture::checked(config);
    const auto bind = [&](unsigned id) {
        return fixture::take(rd::bind_external_mesh_sidecar(
            schema, graph, rd::MeshGeometryId{id}, make_binding_limits(), plan, {}, nullptr));
    };
    const auto bound = bind(2), swapped = bind(3);
    const PgReferenceLimits arithmetic{partition_test::exact_limits(), 20'000'000, 10'000'000, {}};
    const auto si = convert(schema, arithmetic.algebra, plan);
    const coupling::CasePgReferenceLimits pg_limits{
        {arithmetic.algebra, 4, {}},
        {arithmetic.algebra, 2000, 20'000'000, 100'000, {}},
        partition_test::limits(),
        arithmetic,
        100'000,
        {}};
    const auto chart_limits = chart_preparation(pg_limits, 2, 2);
    const coupling::CaseGmshPgChartLimits limits{
        chart_limits,
        {{arithmetic, 8}, {8, 32, 32, 8, 32'000}, 8, 32'000},
        {2, 8, 200'000},
        {arithmetic, 1},
        0,
        1};
    const RectangleOracle oracle{arithmetic, plan};
    const auto weight = oracle.r("1", "6"), third = oracle.r("1", "3");
    std::array<ExactRational, 6> positions;
    std::vector<PgLineQuadratureNode> plane;
    plane.reserve(positions.size());
    for (std::size_t i = 0; i < positions.size(); ++i) {
        positions[i] = oracle.r(std::to_string(2 * i + 1), "12");
        plane.push_back({positions[i], weight});
    }
    const auto low = oracle.r("1", "4"), half = oracle.r("1", "2"), high = oracle.r("3", "4");
    const std::array depth{PgLineQuadratureNode{low, third}, PgLineQuadratureNode{half, third},
                           PgLineQuadratureNode{high, third}};
    const PgVolumeProductRule product{plane, plane, depth};
    const std::array rules{coupling::CasePgRuleBinding{1, product},
                           coupling::CasePgRuleBinding{2, product}};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    coupling::CaseGmshPgChartDiagnostic diagnostic;
    const auto rejected = coupling::prepare_case_gmsh_pg_chart(si, swapped, rules, limits,
                                                               workspace, plan, {}, &diagnostic);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::InvalidInput);
    const auto& joint_problem = required_evidence(diagnostic.joint);
    CHECK(joint_problem.joint == rd::MeshJointId{1});
    CHECK(joint_problem.endpoint == 0);
    CHECK(joint_problem.edge.reason == coupling::CaseGmshSupportBoundaryFailure::Incidence);
    CHECK(joint_problem.edge.line == 502);

    const auto prepared = fixture::take(coupling::prepare_case_gmsh_pg_chart(
        si, bound, rules, limits, workspace, plan, {}, &diagnostic));
    CHECK_FALSE(diagnostic.joint);
    REQUIRE(prepared.joints.size() == 1);
    const auto& joint = required_evidence(prepared.joints.values().front());
    CHECK(joint.joint == rd::MeshJointId{1});
    for (std::size_t i = 0; i < joint.endpoints.size(); ++i) {
        const auto& edge = joint.endpoints[i];
        CHECK(edge.shell == rd::MeshShellId{i + 1});
        REQUIRE(edge.lines.size() == 1);
        CHECK(edge.lines.values()[0].element == 501 + i);
        CHECK(edge.lines.values()[0].incident_face == 701 + i);
        const auto shell =
            fixture::take(prepared.chart.shell_reference(uvlas::numerics::PatchId{i + 1}, 2));
        const auto* basis = fixture::take(shell.patch().conditioning().reference_basis());
        CHECK(fixture::take(basis->reference_identity()) == edge.patch);
    }
    const auto& original = prepared.chart.chart().original();
    REQUIRE(original.sources().size() == 2);
    CHECK(std::ranges::count_if(original.sources(), [](const auto& entry) {
              return entry.kind == PgTraceKind::Connection;
          }) == 1);
    REQUIRE_FALSE(original.original_rows().empty());
    CHECK_FALSE(prepared.chart.chart().pivots().empty());

    std::vector<PgCartesianExactPosition> coordinates;
    REQUIRE(original.layout().size() == 2);
    for (const auto& entry : original.layout()) {
        REQUIRE(static_cast<std::size_t>(entry.first) == coordinates.size());
        const auto shell = fixture::take(prepared.chart.shell_reference(entry.id, 2));
        const auto reference = shell.patch().reference_coordinates();
        REQUIRE(static_cast<std::size_t>(entry.count) == reference.size());
        for (const auto& point : reference) {
            PgCartesianExactPosition copy;
            for (std::size_t d = 0; d < copy.size(); ++d)
                copy[d] = fixture::take(point[d].clone(arithmetic.algebra, plan));
            coordinates.push_back(std::move(copy));
        }
    }
    const auto residual = fixture::take(original.residual(
        {original.layout(), coordinates}, chart_limits.chart.arithmetic, workspace, plan));
    REQUIRE(residual.values().size() == original.original_rows().size());
    for (const auto& value : residual.values())
        for (const auto& component : value)
            CHECK(component.sign() == 0);
}

TEST_CASE("An original SPCADD point support agrees with the same Case Shell",
          "[FD-03][FD-12][IO4-03][MAT-D08][V-MAT][R-05][integration][verification][negative]") {
    using Code = core::FailureCode;
    using Reason = coupling::CaseNastranPointSupportFailure;
    enum class Trial : std::uint8_t {
        Match,
        DifferentComponents,
        PrescribedMotion,
        DifferentPoint,
        Offset
    };
    const auto trial = GENERATE(Trial::Match, Trial::DifferentComponents, Trial::PrescribedMotion,
                                Trial::DifferentPoint, Trial::Offset);
    CAPTURE(static_cast<unsigned>(trial));
    core::ResourcePlan plan{fixture::budget};
    auto config = fixture::root("coupled_transient");
    auto& model = config.field("model");
    auto& pg = model.field("geometry").elements.front();
    pg.field("thickness").set("values", fixture::json("[0.25,0.25,0.25,0.25]"));
    pg.field("frame").set("origin",
                          fixture::json(trial == Trial::Offset ? "[2,3,5.125]" : "[2,3,5]"));
    auto support = fixture::fragment("Support.constraint");
    support.set("components",
                fixture::json(trial == Trial::DifferentComponents ? "[0,2]" : "[1,2]"));
    if (trial == Trial::DifferentPoint)
        support.field("target").set("coordinates", fixture::json("[0.5,0,0]"));
    auto& supports = model.field("supports").elements;
    supports.clear();
    supports.push_back(std::move(support));
    std::string native{R"BDF(GRID,10,,2.0,3.0,5.0,12
GRID,20,,3.0,3.0,5.0
GRID,40,,3.0,4.0,5.0
GRID,99,,2.0,4.0,5.0
CORD2R,12,,0.0,0.0,0.0,0.0,0.0,1.0
,0.0,1.0,0.0
CQUAD4,701,90,10,20,40,99,0.0)BDF"};
    if (trial == Trial::Offset)
        native += ",0.125";
    native += "\nMAT1,6,1.0,,0.0,1.0\nPSHELL,90,6,0.25\nSPC1,31,1,10\nSPC,31,10,3,";
    native += trial == Trial::PrescribedMotion ? "0.001" : "-0.0";
    native += "\nSPC1,61,2,10\nSPCADD,61,31\n";
    auto mapping = native_sidecar(false, false, "2");
    mapping["supports"].push_back(
        {{"set", "61"}, {"grid", "10"}, {"components", {1, 3}}, {"support", "1"}});
    const auto graph = [&] {
        const files::Tree tree;
        write_external(model, tree, false, 2, 1, native, mapping.dump());
        return read_graph(tree, config, plan);
    }();
    const auto schema = fixture::checked(config);
    const auto parent = fixture::take(schema.identity());
    const auto binding_limits = make_binding_limits();
    const auto bound = fixture::take(rd::bind_external_mesh_sidecar(
        schema, graph, rd::MeshGeometryId{2}, binding_limits, plan, {}));
    const auto cards = quantities(graph, parent, io::SourceId{1}, binding_limits, plan);
    const PgReferenceLimits arithmetic{partition_test::exact_limits(), 20'000'000, 10'000'000, {}};
    const auto si = convert(schema, arithmetic.algebra, plan);
    const coupling::CasePgReferenceLimits preparation{
        {arithmetic.algebra, 4, {}},
        {arithmetic.algebra, 2000, 20'000'000, 100'000, {}},
        partition_test::limits(),
        arithmetic,
        100'000,
        {}};
    const rd::ExternalSidecarPgLimits qualification{
        {arithmetic, 8}, {8, 32, 32, 8, 32'000}, 8, 32'000};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    coupling::CaseNastranPgSetDiagnostic geometry_failure;
    auto prepared =
        coupling::CaseNastranPlanarPgSet::prepare(si, bound, cards, preparation, qualification,
                                                  {1, 1, 1}, workspace, plan, &geometry_failure);
    CAPTURE(geometry_failure.member.planar.geometry, geometry_failure.member.material_axes);
    REQUIRE(prepared);
    const auto& shell = fixture::take(prepared->at(0)).get();
    const auto patch = fixture::take(shell.reference().identity());
    std::optional<coupling::CaseNastranPointSupportDiagnostic> diagnostic;
    const auto certify = [&](const coupling::CaseNastranPointSupportLimits& limits) {
        return coupling::certify_case_nastran_static_point_support(
            si, bound, *prepared, rd::NastranSupportSetId{61}, rd::NastranGridId{10}, limits, plan,
            &diagnostic);
    };
    const RectangleOracle oracle{arithmetic, plan};
    const auto plane_weight = oracle.r("1", "6");
    const auto depth_weight = oracle.r("1", "3");
    std::array<ExactRational, 6> positions;
    std::vector<PgLineQuadratureNode> plane;
    plane.reserve(positions.size());
    for (std::size_t i = 0; i < positions.size(); ++i) {
        positions[i] = oracle.r(std::to_string(2 * i + 1), "12");
        plane.push_back({positions[i], plane_weight});
    }
    const auto quarter = oracle.r("1", "4"), half = oracle.r("1", "2"), last = oracle.r("3", "4");
    const std::array depth{PgLineQuadratureNode{quarter, depth_weight},
                           PgLineQuadratureNode{half, depth_weight},
                           PgLineQuadratureNode{last, depth_weight}};
    const std::array rules{
        coupling::CasePgRuleBinding{1, PgVolumeProductRule{plane, plane, depth}}};
    const auto chart_limits = chart_preparation(preparation, 1, 1);
    const coupling::CaseNastranPointPgChartLimits entry_limits{
        {chart_limits, qualification, {1, 1, 1}}, {arithmetic, 2}, 1};
    coupling::CaseNastranPointPgChartDiagnostic entry_diagnostic;
    auto entry = coupling::prepare_case_nastran_pg_chart(si, bound, cards, rules, entry_limits,
                                                         workspace, plan, {}, &entry_diagnostic);
    if (trial != Trial::Match) {
        const auto rejected = certify({arithmetic, 2});
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == Code::InvalidInput);
        REQUIRE_FALSE(entry);
        CHECK(entry.error().code == Code::InvalidInput);
        REQUIRE(entry_diagnostic.support);
        CHECK(entry_diagnostic.support.transform(
                  [](const auto& problem) { return problem.grid; }) == rd::NastranGridId{10});
        CHECK(entry_diagnostic.support.transform([](const auto& problem) {
            return problem.reason;
        }) == diagnostic.transform([](const auto& problem) { return problem.reason; }));
        if (const auto& problem = diagnostic) {
            CHECK(problem->set == rd::NastranSupportSetId{61});
            CHECK(problem->grid == rd::NastranGridId{10});
            REQUIRE(problem->field);
            const auto byte = [](const auto& field) { return field.byte; };
            const auto field_byte = problem->field.transform(byte);
            if (trial == Trial::DifferentComponents) {
                CHECK(problem->reason == Reason::Components);
                CHECK(field_byte == native.find(",1,10") + 1);
                REQUIRE(problem->cd_field);
                CHECK(problem->cd_field.transform(byte) == native.find(",12") + 1);
            } else if (trial == Trial::PrescribedMotion) {
                CHECK(problem->reason == Reason::Motion);
                CHECK(field_byte == native.find("0.001"));
            } else {
                CHECK(problem->reason == Reason::Position);
                CHECK(field_byte == native.find("10"));
                if (trial == Trial::Offset)
                    CHECK(problem->offset_byte == native.find("0.125"));
            }
        } else {
            FAIL("Missing support diagnostic");
        }
        CHECK(fixture::take(shell.reference().identity()) == patch);
        return;
    }
    REQUIRE(entry);
    CHECK_FALSE(entry_diagnostic.support);
    REQUIRE(entry->supports.size() == 1);
    REQUIRE(entry->shells.size() == 1);
    REQUIRE(entry->structural.patch_count() == 1);
    const auto& entry_support = required_evidence(entry->supports.values()[0]);
    const auto* entry_patch = fixture::take(entry->structural.patch(0));
    CHECK(entry_support.mapping.grid == rd::NastranGridId{10});
    CHECK(entry_support.patch == fixture::take(entry_patch->info()).geometry);
    CHECK(entry_support.world_components == 6U);
    auto manual_input = config;
    manual_input.field("model").field("geometry").elements.pop_back();
    const auto manual_case = convert(fixture::checked(manual_input), arithmetic.algebra, plan);
    const auto manual = fixture::take(coupling::prepare_case_component_pg_chart(
        manual_case, rules, chart_limits, workspace, plan));
    CHECK(fixture::take(manual.chart().identity()) ==
          fixture::take(entry->structural.chart().identity()));
    CHECK(std::ranges::equal(manual.chart().free_columns(),
                             entry->structural.chart().free_columns()));

    auto pair_input = manual_input;
    auto& pair_model = pair_input.field("model");
    auto high_geometry = pair_model.field("geometry").elements.front();
    high_geometry.set("id", fixture::string("6"));
    high_geometry.field("frame").set("origin", fixture::json("[2,3,9]"));
    pair_model.field("geometry").elements.push_back(std::move(high_geometry));
    auto high_shell = pair_model.field("shells").elements.front();
    high_shell.set("id", fixture::string("2"));
    high_shell.set("geometry", fixture::string("6"));
    high_shell.field("layers").elements.front().set("id", fixture::string("2"));
    pair_model.field("shells").elements.push_back(std::move(high_shell));
    pair_model.field("shells").elements.front().field("layers").elements.front().set(
        "initial_green", fixture::json(R"({
          "kind":"expression",
          "expression":{
            "domain":{"kind":"reference","coordinates":"chart","lower":[0,0,-1],"upper":[1,1,1]},
            "nodes":[{"op":"constant","value":0.125,"dimensions":[0,0,0,0]},
                     {"op":"constant","value":0,"dimensions":[0,0,0,0]}],
            "outputs":[0,1,1,1,1,1]
          },
          "error":[0,0,0,0,0,0]
        })"));
    auto high_support = pair_model.field("supports").elements.front();
    high_support.set("id", fixture::string("2"));
    high_support.field("target").field("target").set("id", fixture::string("2"));
    pair_model.field("supports").elements.push_back(std::move(high_support));
    const auto pair_case = convert(fixture::checked(pair_input), arithmetic.algebra, plan);
    const std::array pair_rules{rules.front(),
                                coupling::CasePgRuleBinding{2, rules.front().quadrature}};
    const auto pair_limits = chart_preparation(preparation, 2, 2);
    const auto pair_chart = fixture::take(coupling::prepare_case_component_pg_chart(
        pair_case, pair_rules, pair_limits, workspace, plan));
    const auto pair_source = fixture::take(pair_case.source_case());
    const auto pair_model_record = fixture::take(pair_source.get().root()).model.record();
    const auto condition = [&](core::Index index) {
        return std::get<io::schema::SupportConstraint>(
            fixture::take(pair_model_record.supports.at(index)).record());
    };
    const auto lower = condition(0), upper = condition(1);
    const auto lower_patch = fixture::take(pair_chart.point_patch(0));
    const auto upper_patch = fixture::take(pair_chart.point_patch(1));
    const auto assemble = [&](const auto& point, const auto& bound_patch) {
        return coupling::assemble_case_static_pg_point_row(pair_case, point, bound_patch,
                                                           pair_limits.patch, plan);
    };
    CHECK(fixture::take(assemble(lower, lower_patch)).prescribed_metres[2] == oracle.r("5"));
    CHECK(fixture::take(assemble(upper, upper_patch)).prescribed_metres[2] == oracle.r("9"));
    const auto& volume_patch = lower_patch.patch();
    const auto volume_id = fixture::take(volume_patch.info()).identity;
    const auto* volume_basis = fixture::take(volume_patch.conditioning().reference_basis());
    const auto volume_partition = fixture::take(volume_basis->source());
    const auto volume_count = fixture::take(volume_partition.partition.volumes()).size();
    std::vector<PgOwnedDiscreteVolumePoint> volume_rule;
    for (std::size_t j = 0; j < volume_count; ++j)
        volume_rule.push_back({static_cast<core::Index>(j), oracle.r("1", "3"), oracle.r("1", "3"),
                               oracle.r("1", "2"), oracle.r("1", "2")});
    auto shell_volume = fixture::take(coupling::CasePgShellVolume::prepare(
        pair_case, lower_patch, volume_rule, {}, pair_limits.patch, plan));
    // The prepared volume owns its rule. Assemble after releasing the input.
    volume_rule.clear();
    const auto q0 = volume_patch.reference_coordinates();
    std::vector<PgCartesianExactPosition> stretched(q0.size()), stretch(q0.size());
    const auto two = oracle.r("2");
    for (std::size_t j = 0; j < q0.size(); ++j) {
        stretch[j][0] = fixture::take(q0[j][0].clone(arithmetic.algebra, plan));
        stretched[j][0] = fixture::take(q0[j][0].multiply(two, arithmetic.algebra, plan));
        for (std::size_t a = 1; a < 3; ++a)
            stretched[j][a] = fixture::take(q0[j][a].clone(arithmetic.algebra, plan));
    }
    const auto width = oracle.r("1", "100000000000000000000");
    const PgShellVolumeResolution resolution{96, width, width, width};
    const auto action = fixture::take(shell_volume.elastic(
        lower_patch, {volume_id, stretched}, PgShellPatchDirection{volume_id, stretch}, resolution,
        pair_limits.patch, plan));
    const auto energy = oracle.r("121", "512");
    CHECK(fixture::take(action.energy_joules.lower.compare(energy, arithmetic.algebra, plan)) <= 0);
    CHECK(fixture::take(action.energy_joules.upper.compare(energy, arithmetic.algebra, plan)) >= 0);
    CHECK(action.tangent_action_newtons.size() == q0.size());
    const auto work_encloses = [&](std::span<const PgShellCartesianEnclosure> force,
                                   std::span<const PgCartesianExactPosition> direction,
                                   const ExactRational& expected) {
        const auto& exact = arithmetic.algebra;
        ExactRational lo, hi;
        REQUIRE(force.size() == direction.size());
        for (std::size_t j = 0; j < force.size(); ++j) {
            for (std::size_t i = 0; i < 3; ++i) {
                const auto& d = direction[j][i];
                const auto& bounds = force[j][i];
                const bool negative = d.sign() < 0;
                const auto a =
                    fixture::take(d.multiply(negative ? bounds.upper : bounds.lower, exact, plan));
                const auto b =
                    fixture::take(d.multiply(negative ? bounds.lower : bounds.upper, exact, plan));
                lo = fixture::take(lo.add(a, exact, plan));
                hi = fixture::take(hi.add(b, exact, plan));
            }
        }
        CHECK(fixture::take(lo.compare(expected, exact, plan)) <= 0);
        CHECK(fixture::take(hi.compare(expected, exact, plan)) >= 0);
    };
    // One unit square with thickness 1/4, F11=2, E=1, nu=0 and E0,11=1/8:
    // Elastic E11 is 11/8; the work is (1/4)*P11 and (1/4)*dP11/dF11.
    work_encloses(action.internal_force_newtons.values(), stretch, oracle.r("11", "16"));
    work_encloses(action.tangent_action_newtons.values(), stretch, oracle.r("43", "32"));

    const auto& exact = arithmetic.algebra;
    const auto one = oracle.r("1"), five = oracle.r("5");
    // Augment the in-plane displacement by the through-thickness field Z-5.
    for (std::size_t j = 0; j < q0.size(); ++j) {
        const auto middle =
            fixture::take(volume_patch.translation()[j].multiply(five, exact, plan));
        stretch[j][2] = fixture::take(q0[j][2].subtract(middle, exact, plan));
    }
    constexpr std::string_view scale = "1000000000000000000000000000000000000000000";
    const auto outer = oracle.r("69431844202973712388026755553595247452137", scale);
    const auto inner = oracle.r("330009478207571867598667120448377656399712", scale);
    using PointWeight = std::pair<ExactRational, ExactRational>;
    // Four-point Gauss nodes and weights rounded at 42 decimal places.
    const std::array<PointWeight, 4> gauss{
        {{fixture::take(outer.clone(exact, plan)),
          oracle.r("173927422568726928686531974610999703617674", scale)},
         {fixture::take(inner.clone(exact, plan)),
          oracle.r("326072577431273071313468025389000296382326", scale)},
         {fixture::take(one.subtract(inner, exact, plan)),
          oracle.r("326072577431273071313468025389000296382326", scale)},
         {fixture::take(one.subtract(outer, exact, plan)),
          oracle.r("173927422568726928686531974610999703617674", scale)}}};
    const std::array<PointWeight, 7> seven{{{oracle.r("0"), oracle.r("41", "840")},
                                            {oracle.r("1", "6"), oracle.r("216", "840")},
                                            {oracle.r("1", "3"), oracle.r("27", "840")},
                                            {oracle.r("1", "2"), oracle.r("272", "840")},
                                            {oracle.r("2", "3"), oracle.r("27", "840")},
                                            {oracle.r("5", "6"), oracle.r("216", "840")},
                                            {oracle.r("1"), oracle.r("41", "840")}}};
    const auto low = oracle.r("1", "6"), high = oracle.r("2", "3");
    const std::array triangle{std::pair{&low, &low}, std::pair{&high, &low},
                              std::pair{&low, &high}};
    const auto mass_for_rule = [&](std::span<const PointWeight> depth_rule) {
        std::vector<PgOwnedDiscreteVolumePoint> mass_rule;
        for (std::size_t j = 0; j < volume_count; ++j) {
            for (const auto& [u, v] : triangle) {
                for (const auto& [z, w] : depth_rule) {
                    mass_rule.push_back(
                        {static_cast<core::Index>(j), fixture::take(u->clone(exact, plan)),
                         fixture::take(v->clone(exact, plan)), fixture::take(z.clone(exact, plan)),
                         fixture::take(w.multiply(low, exact, plan))});
                }
            }
        }
        const PgShellMassResolution mass_resolution{96, width, width};
        const auto mass =
            fixture::take(shell_volume.mass(lower_patch, mass_rule, {volume_id, stretch},
                                            mass_resolution, pair_limits.patch, plan));
        const auto quarter_mass = oracle.r("1", "4");
        CHECK(fixture::take(mass.mass_kg.lower.compare(quarter_mass, exact, plan)) <= 0);
        CHECK(fixture::take(mass.mass_kg.upper.compare(quarter_mass, exact, plan)) >= 0);
        // Integral of X²+(Z-5)² over [2,3] × [3,4] × [5-1/8,5+1/8].
        work_encloses(mass.action_kg_metres.values(), stretch, oracle.r("1217", "768"));
    };
    mass_for_rule(gauss);
    mass_for_rule(seven);
    core::ResourcePlan no_allocation{0};
    const auto swapped_volume =
        shell_volume.elastic(upper_patch, {volume_id, stretched}, std::nullopt, resolution,
                             pair_limits.patch, no_allocation);
    REQUIRE_FALSE(swapped_volume);
    CHECK(swapped_volume.error().code == Code::StaleGeneration);
    CHECK(no_allocation.allocated_bytes() == 0);
    const auto before_swap = plan.allocated_bytes();
    const auto reject_swap = [&](const auto& point, const auto& wrong_patch) {
        const auto swapped = assemble(point, wrong_patch);
        REQUIRE_FALSE(swapped);
        CHECK(swapped.error().code == Code::InvalidInput);
    };
    reject_swap(lower, upper_patch);
    reject_swap(upper, lower_patch);
    CHECK(plan.allocated_bytes() == before_swap);

    auto mixed_input = manual_input;
    auto edge = fixture::fragment("Support.trace_constraint");
    edge.set("id", fixture::json("\"2\""));
    edge.field("region").set("boundary", fixture::json("1"));
    mixed_input.field("model").field("supports").elements.push_back(std::move(edge));
    auto mixed_limits = entry_limits;
    mixed_limits.base.chart.constraints.maximum_sources = 2;
    const auto mixed_manual = fixture::take(coupling::prepare_case_component_pg_chart(
        convert(fixture::checked(mixed_input), arithmetic.algebra, plan), rules,
        mixed_limits.base.chart, workspace, plan));
    const auto mixed_graph = [&] {
        const files::Tree tree;
        write_external(mixed_input.field("model"), tree, false, 2, 1, native, mapping.dump());
        return read_graph(tree, mixed_input, plan);
    }();
    const auto mixed_schema = fixture::checked(mixed_input);
    const auto mixed_bound = fixture::take(rd::bind_external_mesh_sidecar(
        mixed_schema, mixed_graph, rd::MeshGeometryId{2}, binding_limits, plan, {}));
    const auto mixed_cards = quantities(mixed_graph, fixture::take(mixed_schema.identity()),
                                        io::SourceId{1}, binding_limits, plan);
    const auto mixed = fixture::take(coupling::prepare_case_nastran_pg_chart(
        convert(mixed_schema, arithmetic.algebra, plan), mixed_bound, mixed_cards, rules,
        mixed_limits, workspace, plan));
    const auto& mixed_chart = mixed.structural.chart();
    CHECK(fixture::take(mixed_chart.identity()) == fixture::take(mixed_manual.chart().identity()));
    const auto mixed_sources = mixed_chart.original().sources();
    REQUIRE(mixed_sources.size() == 2);
    CHECK(mixed_sources[0].id == PgConstraintId{2});
    CHECK(mixed_sources[0].world_components == 7U);
    CHECK(mixed_sources[1].id == PgConstraintId{1});
    CHECK(mixed_sources[1].world_components == 6U);
    auto no_support_capacity = entry_limits;
    no_support_capacity.maximum_support_mappings = 0;
    core::ResourcePlan empty{0};
    const auto denied = coupling::prepare_case_nastran_pg_chart(
        si, bound, cards, rules, no_support_capacity, workspace, empty, {}, &entry_diagnostic);
    REQUIRE_FALSE(denied);
    CHECK(denied.error().code == Code::ResourceExceeded);
    CHECK(empty.allocated_bytes() == 0);
    auto cancelled_limits = arithmetic;
    // MSVC implements request_stop() as a non-const member.
    // NOLINTNEXTLINE(misc-const-correctness)
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    cancelled_limits.cancellation = stop.get_token();
    const auto cancelled = certify({cancelled_limits, 2});
    REQUIRE_FALSE(cancelled);
    CHECK(cancelled.error().code == Code::Cancelled);

    config.field("model")
        .field("supports")
        .elements.front()
        .set("components", fixture::json("[0,2]"));
    const auto other_schema = fixture::checked(config);
    const auto other_si = convert(other_schema, arithmetic.algebra, plan);
    const auto stale = coupling::certify_case_nastran_static_point_support(
        other_si, bound, *prepared, rd::NastranSupportSetId{61}, rd::NastranGridId{10},
        {arithmetic, 2}, plan, &diagnostic);
    REQUIRE_FALSE(stale);
    CHECK(stale.error().code == Code::StaleGeneration);
    REQUIRE(diagnostic);
    CHECK(diagnostic.transform([](const auto& problem) { return problem.reason; }) ==
          Reason::Binding);

    const auto limited = certify({arithmetic, 1});
    REQUIRE_FALSE(limited);
    CHECK(limited.error().code == Code::ResourceExceeded);
    CHECK(fixture::take(shell.reference().identity()) == patch);
    const auto result = certify({arithmetic, 2});
    REQUIRE(result);
    CHECK_FALSE(diagnostic);
    CHECK(result->patch == patch);
    CHECK(result->deck.location.source == io::SourceId{1});
    CHECK(result->deck.digest == rd::source_digest(files::bytes_of(native)));
    CHECK(result->sidecar_digest == bound.sidecar().digest());
    CHECK(result->mapping.support == rd::MeshSupportId{1});
    CHECK(result->mapping.component_bits == 5U); // Native x,z.
    CHECK(result->world_components == 6U);       // World y,z from the declared CORD2R.
    CHECK(result->element == 701);
    CHECK(result->chart.s.sign() == 0);
    CHECK(result->chart.t.sign() == 0);
    const auto coordinate = [&](std::string_view word) {
        return fixture::take(ExactRational::from_decimal(word, arithmetic.algebra, plan));
    };
    CHECK(result->world_metres[0] == coordinate("2"));
    CHECK(result->world_metres[1] == coordinate("3"));
    CHECK(result->world_metres[2] == coordinate("5"));

    // Prepare the reference basis without the point constraint. The exact
    // point row must attach to this same physical reference and basis order.
    auto reference_input = config;
    reference_input.field("model").field("supports").elements[0] =
        fixture::fragment("Support.trace_constraint");
    const auto reference_case =
        convert(fixture::checked(reference_input), arithmetic.algebra, plan);
    const auto reference_chart = fixture::take(
        coupling::prepare_case_pg_chart(reference_case, rules, chart_limits, workspace, plan));
    const auto* conditioned = fixture::take(reference_chart.patch(0));
    REQUIRE(conditioned != nullptr);
    const auto point_row = fixture::take(coupling::assemble_case_nastran_static_point_row(
        *result, *conditioned, chart_limits.patch, plan));
    CHECK(point_row.basis == fixture::take(conditioned->info()).identity);
    CHECK(point_row.world_components == 6U);
    REQUIRE(point_row.terms.size() != 0);
    const auto translation = conditioned->translation();
    ExactRational response;
    for (const auto& term : point_row.terms.values()) {
        const auto j = static_cast<std::size_t>(term.column);
        REQUIRE(j < translation.size());
        const auto contribution = fixture::take(
            term.coefficient.multiply(translation[j], chart_limits.patch.algebra, plan));
        response = fixture::take(response.add(contribution, chart_limits.patch.algebra, plan));
    }
    CHECK(response == coordinate("1"));
    CHECK(point_row.prescribed_metres[1] == coordinate("3"));
    CHECK(point_row.prescribed_metres[2] == coordinate("5"));
    const std::array patch_bindings{
        PgConstraintPatchBinding{uvlas::numerics::PatchId{1}, *conditioned}};
    const std::array point_bindings{PgConstraintPointBinding{
        PgConstraintId{point_row.support.value()}, point_row.basis, 0, point_row.world_components,
        point_row.terms.values(), point_row.prescribed_metres}};
    auto global_limits = chart_limits.constraints;
    global_limits.maximum_sources = 2;
    auto original = fixture::take(PgGlobalConstraints::prepare(patch_bindings, {}, point_bindings,
                                                               global_limits, workspace, plan));
    CHECK_FALSE(fixture::take(original.info()).common_cartesian);
    REQUIRE(original.sources().size() == 1);
    CHECK(original.sources()[0].world_components == 6U);
    CHECK(original.original_rows()[0].source_overlap == -1);
    const auto residual =
        fixture::take(original.residual({original.layout(), conditioned->reference_coordinates()},
                                        chart_limits.patch, workspace, plan));
    REQUIRE(residual.size() == 1);
    for (const auto& axis : residual.values()[0])
        CHECK(axis.sign() == 0);

    std::vector<PgCartesianExactPosition> horizontal_motion(translation.size());
    for (std::size_t j = 0; j < translation.size(); ++j)
        horizontal_motion[j][0] =
            fixture::take(translation[j].clone(chart_limits.patch.algebra, plan));
    const auto applied = fixture::take(original.apply({original.layout(), horizontal_motion},
                                                      chart_limits.patch, workspace, plan));
    for (const auto& axis : applied.values()[0])
        CHECK(axis.sign() == 0);
    const std::array<PgCartesianExactPosition, 1> dual{
        PgCartesianExactPosition{coordinate("17"), coordinate("2"), coordinate("-3")}};
    const auto reaction =
        fixture::take(original.transpose(dual, chart_limits.patch, workspace, plan));
    for (const auto& coefficient : reaction.values())
        CHECK(coefficient[0].sign() == 0);
    for (const auto& term : point_row.terms.values()) {
        const auto& coefficient = reaction.values()[static_cast<std::size_t>(term.column)];
        CHECK(coefficient[1] == fixture::take(term.coefficient.multiply(
                                    dual[0][1], chart_limits.patch.algebra, plan)));
        CHECK(coefficient[2] == fixture::take(term.coefficient.multiply(
                                    dual[0][2], chart_limits.patch.algebra, plan)));
    }
    const std::array duplicate{point_bindings[0], point_bindings[0]};
    const auto duplicate_support =
        PgGlobalConstraints::prepare(patch_bindings, {}, duplicate, global_limits, workspace, plan);
    REQUIRE_FALSE(duplicate_support);
    CHECK(duplicate_support.error().code == Code::InvalidInput);
    const auto unusable_legacy =
        PgSparseChart::prepare(std::move(original), chart_limits.chart, workspace, plan);
    REQUIRE_FALSE(unusable_legacy);
    CHECK(unusable_legacy.error().code == Code::Unsupported);
    // NOLINTNEXTLINE(bugprone-use-after-move): rejected preparation preserves the original owner.
    REQUIRE(original.info());
    const auto direct_chart = fixture::take(
        PgComponentSparseChart::prepare(std::move(original), chart_limits.chart, workspace, plan));
    const auto& component_chart = entry->structural.chart();
    CHECK(fixture::take(direct_chart.identity()) == fixture::take(component_chart.identity()));
    const auto component_report = fixture::take(component_chart.preparation());
    CHECK(component_report.factors == 1);
    CHECK(component_report.scalar_rank == 2);
    const auto chart_id = fixture::take(component_chart.identity());
    const auto free_columns = component_chart.free_columns();
    const auto all_coordinates = conditioned->reference_coordinates();
    REQUIRE(free_columns.size() == 3 * all_coordinates.size() - 2);
    std::vector<ExactRational> reference_free, motion_free;
    reference_free.reserve(free_columns.size());
    motion_free.reserve(free_columns.size());
    for (const auto index : free_columns) {
        const auto flat = static_cast<std::size_t>(index);
        reference_free.push_back(fixture::take(
            all_coordinates[flat / 3][flat % 3].clone(chart_limits.patch.algebra, plan)));
        motion_free.push_back(fixture::take(
            horizontal_motion[flat / 3][flat % 3].clone(chart_limits.patch.algebra, plan)));
    }
    const auto lifted = fixture::take(
        component_chart.lift({chart_id, reference_free}, chart_limits.patch, workspace, plan));
    const auto manual_lift = fixture::take(
        manual.chart().lift({chart_id, reference_free}, chart_limits.patch, workspace, plan));
    const auto tangent = fixture::take(
        component_chart.tangent({chart_id, motion_free}, chart_limits.patch, workspace, plan));
    for (std::size_t j = 0; j < all_coordinates.size(); ++j) {
        CHECK(lifted.values()[j] == all_coordinates[j]);
        CHECK(manual_lift.values()[j] == lifted.values()[j]);
        CHECK(tangent.values()[j] == horizontal_motion[j]);
    }
    std::vector<PgCartesianExactPosition> force(all_coordinates.size());
    for (auto& cartesian_force : force)
        cartesian_force = {coordinate("1"), coordinate("2"), coordinate("-3")};
    const auto pulled = fixture::take(
        component_chart.transpose({chart_id, force}, chart_limits.patch, workspace, plan));
    REQUIRE(pulled.original_multipliers.size() == 1);
    CHECK(pulled.original_multipliers.values()[0][0].sign() == 0);
    const auto constrained = fixture::take(component_chart.original().transpose(
        pulled.original_multipliers.values(), chart_limits.patch, workspace, plan));
    std::size_t free_at = 0;
    ExactRational original_work, reduced_work;
    for (std::size_t flat = 0; flat < 3 * force.size(); ++flat) {
        const auto column = flat / 3;
        const auto axis = flat % 3;
        const auto add_product = [&](const ExactRational& sum, const ExactRational& a,
                                     const ExactRational& b) {
            const auto product = fixture::take(a.multiply(b, chart_limits.patch.algebra, plan));
            return fixture::take(sum.add(product, chart_limits.patch.algebra, plan));
        };
        original_work =
            add_product(original_work, force[column][axis], tangent.values()[column][axis]);
        auto recovered = fixture::take(
            constrained.values()[column][axis].clone(chart_limits.patch.algebra, plan));
        if (free_at < free_columns.size() &&
            free_columns[free_at] == static_cast<core::Index>(flat)) {
            const auto& reduced_force = pulled.independent_force.values()[free_at];
            recovered =
                fixture::take(recovered.add(reduced_force, chart_limits.patch.algebra, plan));
            reduced_work = add_product(reduced_work, reduced_force, motion_free[free_at]);
            ++free_at;
        }
        CHECK(recovered == force[column][axis]);
    }
    CHECK(free_at == free_columns.size());
    CHECK(original_work == reduced_work);
    const auto conflicting_position =
        PgCartesianExactPosition{coordinate("2"), coordinate("4"), coordinate("5")};
    const std::array conflict{
        point_bindings[0],
        PgConstraintPointBinding{PgConstraintId{2}, point_row.basis, 0, point_row.world_components,
                                 point_row.terms.values(), conflicting_position}};
    auto overconstrained = fixture::take(
        PgGlobalConstraints::prepare(patch_bindings, {}, conflict, global_limits, workspace, plan));
    const auto refused = PgComponentSparseChart::prepare(std::move(overconstrained),
                                                         chart_limits.chart, workspace, plan);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().code == Code::InvalidInput);
    // NOLINTNEXTLINE(bugprone-use-after-move): rejected preparation preserves the original owner.
    REQUIRE(overconstrained.info());
    auto limited_row = chart_limits.patch;
    limited_row.maximum_partition_visits = 0;
    const auto exhausted =
        coupling::assemble_case_nastran_static_point_row(*result, *conditioned, limited_row, plan);
    REQUIRE_FALSE(exhausted);
    CHECK(exhausted.error().code == Code::ResourceExceeded);
    CHECK(result->cd == rd::NastranCoordinateId{12});
    CHECK(result->grid_field.byte == native.find("10"));
    const auto terms = result->terms.values();
    REQUIRE(terms.size() == 2);
    CHECK(terms[0].source_set == rd::NastranSupportSetId{31});
    CHECK(terms[1].source_set == rd::NastranSupportSetId{31});
    CHECK(terms[0].local_components == 1U);
    CHECK(terms[1].local_components == 4U);
    CHECK(terms[0].card.byte == native.find("SPC1,31"));
    CHECK(terms[1].card.byte == native.find("SPC,31"));
    CHECK_FALSE(terms[0].motion);
    REQUIRE(terms[1].motion);
    CHECK(terms[1].motion.transform([](const auto& field) { return field.byte; }) ==
          native.find("-0.0"));
}
