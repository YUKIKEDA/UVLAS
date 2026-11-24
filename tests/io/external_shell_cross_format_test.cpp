#include "../structure/pg_partition_test_support.hpp"
#include "coupling/case_external_pg.hpp"
#include "coupling/case_nastran_static_point_support.hpp"
#include "reference_file_test_fixture.hpp"
#include "reference_resolution.hpp"
#include "source_digest.hpp"

#include <catch2/generators/catch_generators.hpp>

#include <source_location>

namespace {
namespace io = uvlas::io;
namespace core = uvlas::core;
namespace rd = uvlas::io::detail;
namespace coupling = uvlas::coupling::detail;
namespace fixture = uvlas::test::declaration_fixture;
namespace files = uvlas::test::reference_file_fixture;
using namespace uvlas::structure;

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

const rd::MeshSidecarBindingLimits binding_limits{
    {files::limits().document,
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

rd::ReferencedCase read_graph(const files::Tree& tree, const fixture::JsonToken& config,
                              core::ResourcePlan& plan) {
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
    const auto& binding = binding_limits;
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

TEST_CASE("An original SPCADD point support agrees with the same Case Shell",
          "[FD-03][IO4-03][R-05][integration][verification][negative]") {
    using Code = core::FailureCode;
    using Reason = coupling::CaseNastranPointSupportFailure;
    enum class Trial { Match, DifferentComponents, PrescribedMotion, DifferentPoint, Offset };
    const auto trial = GENERATE(Trial::Match, Trial::DifferentComponents, Trial::PrescribedMotion,
                                Trial::DifferentPoint, Trial::Offset);
    CAPTURE(static_cast<unsigned>(trial));
    core::ResourcePlan plan{fixture::budget};
    auto config = fixture::root("coupled_transient");
    auto& model = config.field("model");
    auto& pg = model.field("geometry").elements.front();
    pg.field("thickness").set("values", fixture::json("[0.25,0.25,0.25,0.25]"));
    if (trial == Trial::Offset)
        pg.field("frame").set("origin", fixture::json("[0,0,0.125]"));
    auto support = fixture::fragment("Support.constraint");
    support.set("components",
                fixture::json(trial == Trial::DifferentComponents ? "[0,2]" : "[1,2]"));
    if (trial == Trial::DifferentPoint)
        support.field("target").set("coordinates", fixture::json("[0.5,0,0]"));
    auto& supports = model.field("supports").elements;
    supports.clear();
    supports.push_back(std::move(support));
    std::string native{R"BDF(GRID,10,,0.0,0.0,0.0,12
GRID,20,,1.0,0.0,0.0
GRID,40,,1.0,1.0,0.0
GRID,99,,0.0,1.0,0.0
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
    if (trial != Trial::Match) {
        const auto rejected = certify({arithmetic, 2});
        REQUIRE_FALSE(rejected);
        CHECK(rejected.error().code == Code::InvalidInput);
        REQUIRE(diagnostic);
        CHECK(diagnostic->set == rd::NastranSupportSetId{61});
        CHECK(diagnostic->grid == rd::NastranGridId{10});
        REQUIRE(diagnostic->field);
        if (trial == Trial::DifferentComponents) {
            CHECK(diagnostic->reason == Reason::Components);
            CHECK(diagnostic->field->byte == native.find(",1,10") + 1);
            REQUIRE(diagnostic->cd_field);
            CHECK(diagnostic->cd_field->byte == native.find(",12") + 1);
        } else if (trial == Trial::PrescribedMotion) {
            CHECK(diagnostic->reason == Reason::Motion);
            CHECK(diagnostic->field->byte == native.find("0.001"));
        } else {
            CHECK(diagnostic->reason == Reason::Position);
            CHECK(diagnostic->field->byte == native.find("10"));
            if (trial == Trial::Offset)
                CHECK(diagnostic->offset_byte == native.find("0.125"));
        }
        CHECK(fixture::take(shell.reference().identity()) == patch);
        return;
    }
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
    CHECK(diagnostic->reason == Reason::Binding);

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
    CHECK(terms[1].motion->byte == native.find("-0.0"));
}
