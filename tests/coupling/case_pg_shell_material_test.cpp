#include "../io/reference_declarations_fixture.hpp"
#include "coupling/case_pg_shell_material.hpp"

#include <catch2/catch_test_macros.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::coupling::detail;
namespace fixture = uvlas::test::declaration_fixture;
namespace pg = uvlas::structure;
constexpr std::uint64_t budget = 64ULL << 20U;

pg::ExactLimits algebra() {
    return {256,
            256,
            256,
            8,
            8,
            8192,
            32768,
            pg::ExactEstimateAdmission{budget},
            pg::ExactTimingMode::Disabled};
}

io::SiCase materialize(const fixture::JsonToken& input) {
    auto checked = fixture::checked(input);
    core::ResourcePlan memory{budget};
    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> quotas;
    quotas.fill(budget);
    const io::SiLimits caps{quotas, algebra()};
    const auto catalog = fixture::take(io::make_si_catalog(checked, caps, memory, nullptr));
    std::vector<io::SiRequirementRow> requirements;
    for (std::uint64_t i = 0; i < fixture::take(catalog.size()); ++i) {
        const auto entry = fixture::take(catalog.entry(i));
        requirements.push_back({entry.key,
                                {io::SiDecimalInput{"1"}, entry.recipe.dimensions},
                                io::SiDecimalInput{"1e-25"},
                                io::SiDecimalInput{"1e-25"}});
    }
    const auto request =
        fixture::take(io::make_si_request(catalog, requirements, caps, memory, nullptr));
    return fixture::take(io::convert_case_to_si(checked, catalog, request, caps, memory, nullptr));
}

TEST_CASE("The Case shell binds ordered layer fields to reusable SI elastic laws",
          "[MAT-D06][MAT-D08][MAT-D09][MAT-D10][V-MAT][negative][verification]") {
    auto document = fixture::root("structural_static");
    document.field("units").set("length", fixture::json(R"("mm")"));
    auto& model = document.field("model");
    auto& materials = model.field("materials").elements;
    materials.front().set("young", fixture::json("100"));
    materials.front().set("poisson", fixture::json("0.25"));
    auto orthotropic = fixture::fragment("Material.orthotropic3d");
    orthotropic.set("id", fixture::json(R"("2")"));
    orthotropic.set("young", fixture::json("[6,8,10]"));
    orthotropic.set("shear", fixture::json("[2,3,4]"));
    materials.push_back(std::move(orthotropic));
    auto effective = fixture::fragment("Material.effective3d");
    effective.set("id", fixture::json(R"("3")"));
    effective.set("components", fixture::json(R"("engineering")"));
    materials.push_back(std::move(effective));

    auto& layers = model.field("shells").elements.front().field("layers").elements;
    using LayerInterval = std::string_view;
    const std::array<LayerInterval, 4> intervals{
        R"([{"numerator":"-1","denominator":"1"},{"numerator":"-1","denominator":"2"}])",
        R"([{"numerator":"-1","denominator":"2"},{"numerator":"0","denominator":"1"}])",
        R"([{"numerator":"0","denominator":"1"},{"numerator":"1","denominator":"2"}])",
        R"([{"numerator":"1","denominator":"2"},{"numerator":"1","denominator":"1"}])"};
    const std::array layer_ids{R"("4")", R"("2")", R"("3")", R"("1")"};
    const std::array material_ids{R"("1")", R"("2")", R"("3")", R"("1")"};
    const auto prototype = layers.front();
    layers.clear();
    for (std::size_t i = 0; i < intervals.size(); ++i) {
        auto layer = prototype;
        layer.set("id", fixture::json(layer_ids[i]));
        layer.set("interval", fixture::json(intervals[i]));
        layer.set("material", fixture::json(material_ids[i]));
        layers.push_back(std::move(layer));
    }
    layers.front().set("initial_green",
                       fixture::json(R"({"kind":"constant","value":[0.125,0,0,0,0,0]})"));
    layers.back().set(
        "initial_green",
        fixture::json(
            R"({"kind":"polynomial","coordinates":"chart","origin":[0,0,0],"lengths":[1,1,1],"components":[[{"powers":[1,0,0],"coefficient":0.125}],[],[],[],[],[]],"error":[0.001,0,0,0,0,0]})"));
    const auto si = materialize(document);
    core::ResourcePlan memory{budget};
    const pg::PgReferenceLimits limits{algebra(), 100'000, 100'000, {}};
    const auto q = [&](std::string_view n, std::string_view d = "1") {
        return fixture::take(pg::ExactRational::from_integer_pair(n, d, limits.algebra, memory));
    };
    const auto owned = CasePgShellMaterials::prepare(si, 1, limits, memory);
    REQUIRE(owned);
    CHECK(owned->law_count() == 3);
    const auto prepared = owned->layers();
    REQUIRE(prepared.size() == 4);
    CHECK(prepared[0].id == core::LayerId{4});
    CHECK(prepared[0].upper == q("-1", "2"));
    CHECK(prepared[3].id == core::LayerId{1});
    CHECK(prepared[3].lower == q("1", "2"));
    const auto& shared = fixture::take(owned->law_for(0)).get();
    CHECK(&shared == &fixture::take(owned->law_for(3)).get());
    CHECK(shared.stiffness()[0].rational == q("120000"));
    CHECK(shared.density_kg_m3() == q("1000000000"));
    CHECK(fixture::take(owned->law_for(1)).get().stiffness()[21].rational == q("4000"));
    CHECK(fixture::take(owned->law_for(2)).get().stiffness()[21].rational == q("2000"));

    const auto first = fixture::take(owned->source_layer(si, 0)).initial_green.record();
    const auto last = fixture::take(owned->source_layer(si, 3)).initial_green.record();
    REQUIRE(std::holds_alternative<io::schema::ReferenceFieldConstant>(first));
    REQUIRE(std::holds_alternative<io::schema::ReferenceFieldPolynomial>(last));
    const std::array chart{q("1", "4"), q("0"), q("-3", "4")};
    pg::PgReferenceVector world, local;
    world[0] = {q("2", "1000"), q("1", "1000")};
    local[0].rational = q("5", "1000");
    const auto g = q("2");
    const CasePgMaterialCoordinates point{chart, local, world, g};
    auto at_first = owned->sample_initial_green(si, 0, point, limits, memory);
    REQUIRE(at_first);
    CHECK(at_first->layer == core::LayerId{4});
    CHECK(at_first->mandel[0].rational == q("1", "8"));
    CHECK(at_first->declared_error[0].sign() == 0);
    std::array<pg::ExactRational, 6> sampled;
    for (std::size_t i = 0; i < sampled.size(); ++i) {
        CHECK(at_first->mandel[i].radical.sign() == 0);
        sampled[i] = std::move(at_first->mandel[i].rational);
    }
    const std::array<pg::PgReferenceScalar, 6> zero;
    const auto result = pg::evaluate_pg_shell_material(
        shared, {sampled, pg::PgMaterialComponents::Mandel}, zero, limits, memory);
    REQUIRE(result);
    CHECK(result->second_pk_mandel_pascal[0].rational == q("-15000"));
    CHECK(result->energy_joules_per_cubic_metre.rational == q("1875", "2"));
    CHECK_FALSE(owned->law_for(4));

    const auto outside = owned->sample_initial_green(si, 3, point, limits, memory);
    REQUIRE_FALSE(outside);
    CHECK(outside.error().code == core::FailureCode::InvalidInput);
    const std::array chart_upper{q("1", "4"), q("0"), q("3", "4")};
    const CasePgMaterialCoordinates upper{chart_upper, local, world, g};
    const auto polynomial = owned->sample_initial_green(si, 3, upper, limits, memory);
    REQUIRE(polynomial);
    CHECK(polynomial->mandel[0].rational == q("1", "32"));
    CHECK(polynomial->mandel[0].radical.sign() == 0);
    CHECK(polynomial->declared_error[0] == q("1", "1000"));
    const std::array<pg::PgMaterialEnclosure, 6> unstrained;
    const auto bounded =
        pg::enclose_pg_shell_material(shared, polynomial->enclosure(), unstrained, limits, memory);
    REQUIRE(bounded);
    // For E=0, the two endpoints come from e0xx=1/32 +/- 1/1000.
    CHECK(bounded->second_pk_mandel_pascal[0].lower == q("-3870"));
    CHECK(bounded->second_pk_mandel_pascal[0].upper == q("-3630"));
    CHECK(bounded->energy_joules_per_cubic_metre.lower == q("43923", "800"));
    CHECK(bounded->energy_joules_per_cubic_metre.upper == q("49923", "800"));
    const std::array chart_right{q("3", "4"), q("0"), q("3", "4")};
    const auto second =
        owned->sample_initial_green(si, 3, {chart_right, local, world, g}, limits, memory);
    REQUIRE(second);
    CHECK(second->mandel[0].rational == q("3", "32"));

    for (const bool use_world : {false, true}) {
        auto positioned = document;
        auto& field = positioned.field("model")
                          .field("shells")
                          .elements.front()
                          .field("layers")
                          .elements.back()
                          .field("initial_green");
        field.set("coordinates",
                  fixture::json(use_world ? R"("world_reference")" : R"("material")"));
        const auto source = materialize(positioned);
        const auto model_at = CasePgShellMaterials::prepare(source, 1, limits, memory);
        REQUIRE(model_at);
        const auto read = model_at->sample_initial_green(source, 3, upper, limits, memory);
        REQUIRE(read);
        CHECK(read->mandel[0].rational == (use_world ? q("1", "4") : q("5", "8")));
        CHECK(read->mandel[0].radical == (use_world ? q("1", "8") : q("0")));
        CHECK(read->g == q("2"));
        CHECK(read->declared_error[0] == q("1", "1000"));
    }

    materials.front().set("density", fixture::json("2"));
    const auto other = materialize(document);
    const auto stale = owned->source_layer(other, 0);
    REQUIRE_FALSE(stale);
    CHECK(stale.error().code == core::FailureCode::StaleGeneration);
}
} // namespace
