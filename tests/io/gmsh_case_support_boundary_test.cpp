#include "../structure/pg_partition_test_support.hpp"
#include "coupling/case_gmsh_support_boundary.hpp"
#include "reference_file_test_fixture.hpp"
#include "reference_resolution.hpp"
#include "source_digest.hpp"

namespace {
namespace io = uvlas::io;
namespace core = uvlas::core;
namespace rd = uvlas::io::detail;
namespace coupling = uvlas::coupling::detail;
namespace fixture = uvlas::test::declaration_fixture;
namespace files = uvlas::test::reference_file_fixture;
using namespace uvlas::structure;
using Reason = coupling::CaseGmshSupportBoundaryFailure;
using Code = core::FailureCode;

enum class Variant : std::uint8_t {
    Complete,
    Missing,
    WrongEdge,
    Duplicate,
    Interior,
    CurvedExact,
    CurvedGap
};

std::string mesh(Variant variant) {
    constexpr std::string_view prefix = R"msh($MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 1 1 0
11 1000 2000 3000 1000 4000 5000 1 7 0
17 1000 2000 3000 1000 4000 5000 1 7 0
$EndEntities
$Nodes
1 6 10 60
2 17 0 6
10 20 30 40 50 60
1000 2000 3000
1000 4000 3000
1000 2000 4000
1000 4000 4000
1000 2000 5000
1000 4000 5000
$EndNodes
$Elements
)msh";
    std::string wire{prefix};
    if (variant == Variant::Missing) {
        wire += "2 3 501 702\n1 11 1 1\n501 30 10\n";
    } else {
        wire += "2 4 501 702\n1 11 1 2\n";
        wire += variant == Variant::WrongEdge  ? "501 10 20\n"
                : variant == Variant::Interior ? "501 30 40\n"
                                               : "501 30 10\n";
        wire += variant == Variant::Duplicate ? "502 10 30\n" : "502 30 50\n";
    }
    wire += "2 17 3 2\n701 10 20 40 30\n702 30 40 60 50\n$EndElements\n";
    return wire;
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

nlohmann::json sidecar() {
    const nlohmann::json zero{{"numerator", "0"}, {"denominator", "1"}};
    const nlohmann::json half{{"numerator", "1"}, {"denominator", "2"}};
    const nlohmann::json one{{"numerator", "1"}, {"denominator", "1"}};
    const auto chart = [&](const auto& t0, const auto& t1) {
        return nlohmann::json{
            {"orientation", "same"},
            {"vertices", nlohmann::json::array({{zero, t0}, {one, t0}, {one, t1}, {zero, t1}})}};
    };
    return {
        {"schema", "uvlas.gmsh-sidecar/1"},
        {"geometry", "2"},
        {"length", "mm"},
        {"supports",
         nlohmann::json::array({{{"physical", {{"dimension", 1}, {"tag", 7}}}, {"support", "1"}}})},
        {"joints", nlohmann::json::array()},
        {"shells", nlohmann::json::array(
                       {{{"physical", {{"dimension", 2}, {"tag", 7}}},
                         {"shell", "1"},
                         {"faces", nlohmann::json::array(
                                       {{{"element", "701"}, {"chart", chart(zero, half)}},
                                        {{"element", "702"}, {"chart", chart(half, one)}}})}}})}};
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

TEST_CASE("A Gmsh support group covers the exact named Case boundary",
          "[FD-03][FD-05][IO4-04][R-05][N4-02][integration][verification][negative]") {
    Variant variant = Variant::Complete;
    SECTION("Both original lines cover the whole Case boundary") {}
    SECTION("A proper subset touches the shell but leaves a gap") {
        variant = Variant::Missing;
    }
    SECTION("An original free edge lies on another named boundary") {
        variant = Variant::WrongEdge;
    }
    SECTION("Distinct element IDs duplicate the same original edge") {
        variant = Variant::Duplicate;
    }
    SECTION("A line touching two selected faces is interior") {
        variant = Variant::Interior;
    }
    SECTION("A nonzero face gap vanishes along the whole support line") {
        variant = Variant::CurvedExact;
    }
    SECTION("Matching endpoints conceal a curved gap along the support line") {
        variant = Variant::CurvedGap;
    }

    auto config = fixture::root("coupled_transient");
    auto& model = config.field("model");
    model.field("supports").elements[0] = fixture::fragment("Support.trace_constraint");
    auto external = fixture::fragment("Geometry.gmsh");
    auto fit = fixture::fragment("ReferenceGeometry.fit");
    auto& pg = model.field("geometry").elements.front();
    pg.field("frame").set("origin", fixture::json("[1,2,3]"));
    pg.field("frame").set("rotation", fixture::json("[[0,0,1],[1,0,0],[0,1,0]]"));
    pg.set("reference_length", fixture::json("2"));
    if (variant == Variant::CurvedExact) {
        pg.set("degree", fixture::json("[2,0]"));
        pg.set(
            "coefficients",
            fixture::json(
                R"([{"powers":[1,0],"coefficient":{"numerator":"1","denominator":"10"}},{"powers":[2,0],"coefficient":{"numerator":"-1","denominator":"10"}}])"));
        fit.set("degree", fixture::json("2"));
    }
    if (variant == Variant::CurvedGap) {
        // t(t-1/2)(t-1)/10. Every native endpoint is on the PG,
        // but the physical +X gap at t=1/4 is 2*(3/640) metres.
        pg.set("degree", fixture::json("[0,3]"));
        pg.set(
            "coefficients",
            fixture::json(
                R"([{"powers":[0,1],"coefficient":{"numerator":"1","denominator":"20"}},{"powers":[0,2],"coefficient":{"numerator":"-3","denominator":"20"}},{"powers":[0,3],"coefficient":{"numerator":"1","denominator":"10"}}])"));
        fit.set("degree", fixture::json("3"));
    }
    const auto msh_text = mesh(variant);
    const auto map_text = sidecar().dump();
    external.set("reference", std::move(fit));
    files::declare_path(external.field("mesh"), "external.msh", hex_digest(msh_text));
    files::declare_path(external.field("sidecar"), "sidecar.json", hex_digest(map_text));
    model.field("geometry").elements.push_back(std::move(external));

    core::ResourcePlan plan{fixture::budget};
    const rd::MeshSidecarBindingLimits binding{
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
    const auto graph = [&] {
        const files::Tree tree;
        files::write_file(tree.root / "case.json", config.wire());
        files::write_file(tree.root / "external.msh", msh_text);
        files::write_file(tree.root / "sidecar.json", map_text);
        rd::ReferenceInputDiagnostic problem;
        const auto input = rd::ReferenceInput::from_file(tree.root, "case.json", io::SourceId{12},
                                                         files::manifest_limits, files::limits(),
                                                         plan, {}, &problem);
        CAPTURE(problem.schema.path());
        REQUIRE(input);
        auto declarations = files::prepare(*input);
        return fixture::take(rd::resolve_leaf_references(
            *input, std::move(declarations),
            {{{4096, 4097, 2, 4096}, 8, 256}, 2, 8192, 8194, 4, 8192, 8}, {2, 2, 1, 128, 128},
            {binding.gmsh, binding.nastran, binding.sidecar}, plan, {}, nullptr, nullptr));
    }();

    const auto schema = fixture::checked(config);
    const PgReferenceLimits arithmetic{partition_test::exact_limits(), 20'000'000, 10'000'000, {}};
    const auto si = convert(schema, arithmetic.algebra, plan);
    const auto bound = fixture::take(rd::bind_external_mesh_sidecar(
        schema, graph, rd::MeshGeometryId{2}, binding, plan, {}, nullptr));
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
    const auto set = fixture::take(coupling::CaseGmshPgSet::prepare(
        si, bound, preparation, qualification, {1, 3, 200'000}, workspace, plan));
    std::optional<coupling::CaseGmshSupportBoundaryDiagnostic> diagnostic;
    const rd::MeshSupportId support{1};
    const auto choose = [&](const coupling::CaseGmshSupportBoundaryLimits& limits) {
        return coupling::certify_case_gmsh_support_boundary(
            si, bound, set, support, rd::GmshKey{1, 7}, limits, plan, &diagnostic);
    };
    const bool curved = variant == Variant::CurvedExact || variant == Variant::CurvedGap;
    if (curved) {
        const auto& shell = fixture::take(set.at(0)).get();
        const auto& conversion = std::get<rd::GmshPgConversionEvidence>(shell.match().geometry);
        if (const auto view = conversion.geometry.view()) {
            REQUIRE(view->faces.size() == 2);
            const auto gap = fixture::take(ExactRational::from_integer_pair(
                variant == Variant::CurvedExact ? "1" : "3",
                variant == Variant::CurvedExact ? "20" : "320", arithmetic.algebra, plan));
            for (const auto& face : view->faces) {
                CHECK(std::is_gteq(fixture::take(
                    face.distance_upper_metres.compare(gap, arithmetic.algebra, plan))));
            }
        } else {
            FAIL("Missing authenticated fit geometry");
        }
    }
    if (variant != Variant::Complete && variant != Variant::CurvedExact) {
        const auto result = choose({arithmetic, 2});
        REQUIRE_FALSE(result);
        CHECK(result.error().code == Code::InvalidInput);
        if (const auto& problem = diagnostic) {
            CHECK(problem->support == support);
            CHECK(problem->physical == rd::GmshKey{1, 7});
            CHECK(problem->boundary == 0);
            CHECK(problem->declaration);
            switch (variant) {
            case Variant::Missing:
                CHECK(problem->reason == Reason::IncompleteBoundary);
                CHECK_FALSE(problem->line);
                break;
            case Variant::WrongEdge:
                CHECK(problem->reason == Reason::OutsideBoundary);
                CHECK(problem->line == 501);
                CHECK(problem->face == 701);
                break;
            case Variant::Duplicate:
                CHECK(problem->reason == Reason::DuplicateLine);
                CHECK(problem->line == 502);
                break;
            case Variant::Interior:
                CHECK(problem->reason == Reason::Incidence);
                CHECK(problem->line == 501);
                break;
            case Variant::CurvedGap:
                CHECK(problem->reason == Reason::Geometry);
                CHECK(problem->line == 501);
                CHECK(problem->face == 701);
                CHECK(problem->world_component == 0);
                break;
            case Variant::CurvedExact:
            case Variant::Complete:
                FAIL("Unexpected complete case in failure check");
            }
        } else {
            FAIL("Missing support boundary provenance");
        }
        return;
    }

    const auto& g = std::get<rd::GmshSidecar>(bound.sidecar().payload());
    REQUIRE(g.supports.size() == 1);
    const auto before = plan.allocated_bytes();
    const auto short_result = choose({arithmetic, 1});
    REQUIRE_FALSE(short_result);
    CHECK(short_result.error().code == Code::ResourceExceeded);
    CHECK(plan.allocated_bytes() == before);
    auto cancelled = arithmetic;
    cancelled.cancellation = partition_test::requested_stop_token();
    const auto stopped = choose({cancelled, 2});
    REQUIRE_FALSE(stopped);
    CHECK(stopped.error().code == Code::Cancelled);
    CHECK_FALSE(diagnostic);
    CHECK(plan.allocated_bytes() == before);

    const auto result = choose({arithmetic, 2});
    REQUIRE(result);
    CHECK_FALSE(diagnostic);
    CHECK(result->mesh.source == io::SourceId{1});
    CHECK(result->mesh_digest == rd::source_digest(files::bytes_of(msh_text)));
    CHECK(result->sidecar.source == io::SourceId{2});
    CHECK(result->sidecar_digest == rd::source_digest(files::bytes_of(map_text)));
    CHECK(result->support == support);
    CHECK(result->shell == rd::MeshShellId{1});
    CHECK(result->boundary == 0);
    CHECK(result->physical == rd::GmshKey{1, 7});
    CHECK(result->declaration == g.supports.values().front().node);
    const auto& patch = fixture::take(set.at(0)).get().reference();
    CHECK(result->patch == fixture::take(patch.identity()));
    const auto lines = result->lines.values();
    REQUIRE(lines.size() == 2);
    CHECK(lines[0].element == 501);
    CHECK(lines[0].nodes == std::array<std::uint64_t, 2>{30, 10});
    CHECK(lines[0].incident_face == 701);
    CHECK(lines[1].element == 502);
    CHECK(lines[1].nodes == std::array<std::uint64_t, 2>{30, 50});
    CHECK(lines[1].incident_face == 702);
    const auto expected = curved ? coupling::CaseGmshSupportPositionEvidence::EdgePolynomial
                                 : coupling::CaseGmshSupportPositionEvidence::FaceEquality;
    CHECK(lines[0].position == expected);
    CHECK(lines[1].position == expected);
}
} // namespace
