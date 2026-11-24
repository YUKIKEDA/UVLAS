#include "../structure/pg_partition_test_support.hpp"
#include "coupling/case_external_pg.hpp"
#include "reference_file_test_fixture.hpp"
#include "reference_resolution.hpp"
#include "source_digest.hpp"

#include <compare>

namespace {
namespace io = uvlas::io;
namespace core = uvlas::core;
namespace rd = uvlas::io::detail;
namespace fixture = uvlas::test::declaration_fixture;
using namespace uvlas::test::reference_file_fixture;

} // namespace

TEST_CASE("One Gmsh fit prepares every mapped Case Shell",
          "[FD-03][IO4-01][IO4-04][R-05][N4-02][integration][verification][negative]") {
    bool wrong_thickness = false;
    bool short_distance = false;
    SECTION("Matching group thickness and distance") {}
    SECTION("The later curved Shell still requires exact thickness") {
        wrong_thickness = true;
    }
    SECTION("A later curved Shell can fail the shared distance budget") {
        short_distance = true;
    }

    constexpr std::string_view two_planes = R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 1 2 0
11 0 0 0 1000 0 0 1 7 0
17 0 0 0 1000 1000 0 1 7 0
18 0 0 1000 1000 1000 1000 1 7 0
$EndEntities
$Nodes
2 8 10 199
2 17 0 4
10 20 40 99
0 0 0
1000 0 0
1000 1000 0
0 1000 0
2 18 0 4
110 120 140 199
0 0 1000
1000 0 1000
1000 1000 1000
0 1000 1000
$EndNodes
$Elements
3 3 501 702
1 11 1 1
501 10 20
2 17 3 1
701 10 20 40 99
2 18 3 1
702 110 120 140 199
$EndElements
)MSH";
    const nlohmann::json zero{{"numerator", "0"}, {"denominator", "1"}};
    const nlohmann::json one{{"numerator", "1"}, {"denominator", "1"}};
    const nlohmann::json chart{
        {"orientation", "same"},
        {"vertices", nlohmann::json::array({{zero, zero}, {one, zero}, {one, one}, {zero, one}})}};
    const nlohmann::json group{{"dimension", 2}, {"tag", 7}};
    const nlohmann::json sidecar{
        {"schema", "uvlas.gmsh-sidecar/1"},
        {"geometry", "2"},
        {"length", "mm"},
        {"supports",
         nlohmann::json::array({{{"physical", {{"dimension", 1}, {"tag", 7}}}, {"support", "1"}}})},
        {"joints", nlohmann::json::array()},
        {"shells",
         nlohmann::json::array(
             {{{"physical", group},
               {"shell", "2"},
               {"full_thickness", 250},
               {"faces", nlohmann::json::array({{{"element", "702"}, {"chart", chart}}})},
               {"material_layers", nlohmann::json::array({{{"layer", "2"}, {"material", "1"}}})}},
              {{"physical", group},
               {"shell", "1"},
               {"faces", nlohmann::json::array({{{"element", "701"}, {"chart", chart}}})},
               {"material_layers",
                nlohmann::json::array({{{"layer", "1"}, {"material", "1"}}})}}})}};
    const auto hex = [](std::string_view body) {
        constexpr std::string_view digits = "0123456789abcdef";
        const auto digest = rd::source_digest(bytes_of(body));
        std::string text;
        for (const auto byte : digest.bytes) {
            const auto value = std::to_integer<unsigned>(byte);
            text.push_back(digits[value >> 4U]);
            text.push_back(digits[value & 15U]);
        }
        return text;
    };
    auto sidecar_text = sidecar.dump();
    auto config = fixture::root("coupled_transient");
    config.field("model")
        .field("geometry")
        .elements[0]
        .field("thickness")
        .set("values", fixture::json("[0.25,0.25,0.25,0.25]"));
    auto external = fixture::fragment("Geometry.gmsh");
    external.set("reference", fixture::fragment("ReferenceGeometry.fit"));
    auto& fit = external.field("reference");
    fit.set("degree", fixture::json("2"));
    fit.set("distance_budget", fixture::json(short_distance ? "0.024" : "0.051"));
    declare_path(external.field("mesh"), "external.dat", hex(two_planes));
    declare_path(external.field("sidecar"), "sidecar.json", hex(sidecar_text));
    config.field("model").field("geometry").elements.push_back(std::move(external));

    auto& model = config.field("model");
    auto second_pg = model.field("geometry").elements.front();
    second_pg.set("id", fixture::json(R"("3")"));
    second_pg.field("frame").set("origin", fixture::json("[0,0,1]"));
    model.field("geometry").elements.push_back(std::move(second_pg));

    auto& shells = model.field("shells").elements;
    auto& first_layer = shells.front().field("layers").elements.front();
    first_layer.set(
        "axis",
        fixture::json(
            R"({"kind":"projected","direction":{"kind":"constant","value":[1,0,0]},"angle":{"kind":"constant","value":[0]}})"));
    auto second_shell = shells.front();
    second_shell.set("id", fixture::json(R"("2")"));
    second_shell.set("geometry", fixture::json(R"("3")"));
    second_shell.field("layers").elements.front().set("id", fixture::json(R"("2")"));
    shells.push_back(std::move(second_shell));
    {
        auto& curved = model.field("geometry").elements.front();
        curved.set("degree", fixture::json("[2,0]"));
        curved.set(
            "coefficients",
            fixture::json(
                R"([{"powers":[1,0],"coefficient":{"numerator":"1","denominator":"10"}},{"powers":[2,0],"coefficient":{"numerator":"-1","denominator":"10"}}])"));
        if (wrong_thickness) {
            curved.field("thickness").set("values", fixture::json("[0.5,0.5,0.5,0.5]"));
        }
    }
    const rd::MeshSidecarBindingLimits binding{
        {limits().document,
         {256, 256, 256, 8, 8, 4096, 8192,
          uvlas::structure::ExactEstimateAdmission{fixture::budget},
          uvlas::structure::ExactTimingMode::Disabled},
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
    auto graph = [&] {
        const Tree tree;
        write_file(tree.root / "case.json", config.wire());
        write_file(tree.root / "external.dat", two_planes);
        write_file(tree.root / "sidecar.json", sidecar_text);
        core::ResourcePlan input_plan{fixture::budget};
        const auto input = fixture::take(
            rd::ReferenceInput::from_file(tree.root, "case.json", io::SourceId{12}, manifest_limits,
                                          limits(), input_plan, {}, nullptr));
        auto declarations = prepare(input);
        core::ResourcePlan graph_plan{fixture::budget};
        rd::LeafReferenceDiagnostic diagnostic;
        auto result = rd::resolve_leaf_references(
            input, std::move(declarations),
            {{{4096, 4097, 2, 4096}, 8, 256}, 2, 8192, 8194, 4, 8192, 8}, {2, 2, 1, 128, 128},
            {binding.gmsh, binding.nastran, binding.sidecar}, graph_plan, {}, nullptr, &diagnostic);
        return result;
    }(); // The graph now owns the exact source bytes. All input files have been removed.

    REQUIRE(graph);
    const auto schema = fixture::checked(config);
    core::ResourcePlan plan{fixture::budget};
    const auto result = rd::bind_external_mesh_sidecar(schema, *graph, rd::MeshGeometryId{2},
                                                       binding, plan, {}, nullptr);
    REQUIRE(result);
    const auto faces = result->faces();
    REQUIRE(faces.size() == 2);
    using namespace uvlas::structure;
    const PgReferenceLimits pg_limits{partition_test::exact_limits(), 20'000'000, 10'000'000, {}};
    const uvlas::coupling::detail::CasePgReferenceLimits preparation{
        {pg_limits.algebra, 4, {}},
        {pg_limits.algebra, 2000, 20'000'000, 100'000, {}},
        partition_test::limits(),
        pg_limits,
        100'000,
        {}};
    core::Workspace workspace{core::PlannedBuffer<std::byte>{}};
    const rd::ExternalSidecarPgLimits qualification{
        {pg_limits, 8}, {8, 32, 32, 8, 32'000}, 8, 32'000};
    const auto to_si = [&](const io::SchemaCase& typed_case) {
        std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> si_caps{};
        si_caps.fill(fixture::budget);
        const io::SiLimits si_limits{si_caps, pg_limits.algebra};
        const auto catalog =
            fixture::take(io::make_si_catalog(typed_case, si_limits, plan, nullptr));
        std::vector<io::SiRequirementRow> requirements;
        for (std::uint64_t i = 0; i < fixture::take(catalog.size()); ++i) {
            const auto entry = fixture::take(catalog.entry(i));
            requirements.push_back({entry.key,
                                    {io::SiDecimalInput{"1"}, entry.recipe.dimensions},
                                    io::SiDecimalInput{"1e-25"},
                                    io::SiDecimalInput{"1e-25"}});
        }
        const auto request =
            fixture::take(io::make_si_request(catalog, requirements, si_limits, plan, nullptr));
        return fixture::take(
            io::convert_case_to_si(typed_case, catalog, request, si_limits, plan, nullptr));
    };
    const auto si = to_si(schema);
    using uvlas::coupling::detail::CaseGmshPgSet;
    using uvlas::coupling::detail::CaseGmshPgSetDiagnostic;
    using uvlas::coupling::detail::CaseGmshPgSetLimits;
    const CaseGmshPgSetLimits batch_limits{2, 3, 200'000};
    CaseGmshPgSetDiagnostic batch_diagnostic;
    const auto prepared = CaseGmshPgSet::prepare(si, *result, preparation, qualification,
                                                 batch_limits, workspace, plan, &batch_diagnostic);
    CHECK(faces[0].shell == rd::MeshShellId{2});
    CHECK(faces[1].shell == rd::MeshShellId{1});
    if (wrong_thickness || short_distance) {
        REQUIRE_FALSE(prepared);
        CHECK(prepared.error().code == (short_distance ? core::FailureCode::PrecisionUnresolved
                                                       : core::FailureCode::InvalidInput));
        CHECK(batch_diagnostic.shell == rd::MeshShellId{1});
        CHECK(batch_diagnostic.first_declaration == faces[1].sidecar_node);
        if (short_distance) {
            if (const auto& detail = batch_diagnostic.member.binding) {
                CHECK(detail->reason == rd::MeshSidecarFailure::Chart);
                CHECK(detail->element == 701);
            } else {
                FAIL("Missing later curved face diagnostic");
            }
        } else if (const auto& detail = batch_diagnostic.member.regional_thickness) {
            CHECK(detail->reason == rd::ExternalPgThicknessFailure::Difference);
            CHECK(detail->element == 701);
            CHECK(detail->sidecar_node);
            CHECK(detail->case_triangle);
        } else {
            FAIL("Missing exact thickness diagnostic on bounded map");
        }
        const auto rebuilt = fixture::take(uvlas::coupling::detail::prepare_case_pg_reference(
            si, 2, preparation, workspace, plan));
        const auto first_again = rd::prepare_gmsh_sidecar_pg(
            si, *result, rd::MeshShellId{2}, rebuilt, qualification, 200'000, plan, {});
        REQUIRE(first_again);
        CHECK(first_again->membership.size() == 1);
        return;
    }
    REQUIRE(prepared);
    REQUIRE(prepared->size() == 2);
    CHECK_FALSE(batch_diagnostic.shell);
    const auto one_quarter =
        fixture::take(ExactRational::from_integer_pair("1", "4", pg_limits.algebra, plan));
    const auto actual_gap =
        fixture::take(ExactRational::from_integer_pair("1", "40", pg_limits.algebra, plan));
    const auto distance_budget =
        fixture::take(ExactRational::from_integer_pair("51", "1000", pg_limits.algebra, plan));
    const auto normal_budget =
        fixture::take(ExactRational::from_integer_pair("1", "1", pg_limits.algebra, plan));
    for (core::Index i = 0; i < 2; ++i) {
        const auto& member = fixture::take(prepared->at(i)).get();
        const auto& match = member.match();
        const auto element = i == 0 ? 702U : 701U;
        CHECK(match.shell == rd::MeshShellId{i == 0 ? 2U : 1U});
        CHECK(match.pg == rd::MeshGeometryId{i == 0 ? 3U : 1U});
        CHECK(match.sidecar_digest == result->sidecar().digest());
        const auto membership = match.membership.values();
        REQUIRE(membership.size() == 1);
        CHECK(membership.front().physical == rd::GmshKey{2, 7});
        CHECK(membership.front().element == element);
        CHECK(membership.front().declaration == faces[static_cast<std::size_t>(i)].sidecar_node);
        const auto* conversion_result = std::get_if<rd::GmshPgConversionEvidence>(&match.geometry);
        REQUIRE(conversion_result != nullptr);
        CHECK(conversion_result->distance_budget_metres == distance_budget);
        CHECK(conversion_result->normal_budget == normal_budget);
        REQUIRE(conversion_result->observations.size() == 1);
        CHECK(conversion_result->observations.values().front().id == core::ObservationId{1});
        if (const auto view = conversion_result->geometry.view()) {
            CHECK(view->patch == fixture::take(member.reference().identity()));
            CHECK(view->source.source == io::SourceId{1});
            REQUIRE(view->faces.size() == 1);
            CHECK(view->faces.front().element == element);
            const auto& distance = view->faces.front().distance_upper_metres;
            if (i == 0) {
                CHECK(std::is_lt(
                    fixture::take(distance.compare(actual_gap, pg_limits.algebra, plan))));
            } else {
                CHECK(std::is_gteq(
                    fixture::take(distance.compare(actual_gap, pg_limits.algebra, plan))));
            }
        } else {
            FAIL("Missing owned Gmsh bounded view");
        }
        if (const auto& thickness = match.thickness) {
            REQUIRE(thickness->matches.size() == 1);
            CHECK(thickness->matches.values().front().element == element);
            CHECK(std::get<ExactRational>(thickness->matches.values().front().value_metres) ==
                  one_quarter);
        } else {
            FAIL("Missing thickness inherited from the same physical group");
        }
    }
    auto short_limits = batch_limits;
    --short_limits.maximum_shells;
    const auto short_result = CaseGmshPgSet::prepare(si, *result, preparation, qualification,
                                                     short_limits, workspace, plan);
    REQUIRE_FALSE(short_result);
    CHECK(short_result.error().code == core::FailureCode::ResourceExceeded);
}
