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
    document.field("units").set("angle", fixture::json(R"("deg")"));
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
    layers.front().set(
        "axis",
        fixture::json(
            R"({"kind":"explicit","columns":{"kind":"constant","value":[1,0,0,0,1,0,0,0,1]}})"));
    layers.back().set(
        "axis",
        fixture::json(
            R"({"kind":"projected","direction":{"kind":"constant","value":[1,0,0]},"angle":{"kind":"polynomial","coordinates":"chart","origin":[0,0,0],"lengths":[1,1,1],"components":[[{"powers":[1,0,0],"coefficient":360}]],"error":[0]}})"));
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
    const auto explicit_axes = owned->sample_axes(si, 0, point, limits, memory);
    REQUIRE(explicit_axes);
    const auto& explicit_columns = std::get<CasePgExplicitAxisValue>(explicit_axes->value).columns;
    CHECK(explicit_columns[4].rational == q("1"));
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
    // An undeformed point retains the stress from the Case's initial field.
    pg::PgShellTensorEnclosure identity;
    for (std::size_t i = 0; i < 3; ++i) {
        identity[4 * i] = {q("1"), q("1")};
    }
    const auto case_axes = owned->sample_axes(si, 3, upper, limits, memory);
    REQUIRE(case_axes);
    const auto& axis = std::get<CasePgProjectedAxisValue>(case_axes->value);
    CHECK(axis.cosine.lower.sign() == 0);
    CHECK(axis.cosine.upper.sign() == 0);
    CHECK(axis.sine.lower == q("1"));
    CHECK(axis.sine.upper == q("1"));
    pg::PgReferenceValue reference{};
    reference.g = q("2");
    reference.layer = core::LayerId{1};
    for (std::size_t i = 0; i < 3; ++i)
        reference.jacobian_columns_metres[i][i].rational = q("1");
    const pg::PgShellPatchPoint trial{{}, std::move(reference), {}, {}, std::move(identity), {}};
    const auto material_point = pg::enclose_pg_shell_material_point(
        shared, polynomial->enclosure(), case_axes->axes(), trial, limits, memory);
    REQUIRE(material_point);
    CHECK(material_point->reference_axes[3].lower == q("1"));
    CHECK(material_point->second_pk_world_pascal[4].lower == q("-3870"));
    CHECK(material_point->first_pk_world_pascal[4].lower == q("-3870"));
    CHECK(material_point->cauchy_world_pascal[4].upper == q("-3630"));
    const std::array chart_right{q("3", "4"), q("0"), q("3", "4")};
    const auto second =
        owned->sample_initial_green(si, 3, {chart_right, local, world, g}, limits, memory);
    REQUIRE(second);
    CHECK(second->mandel[0].rational == q("3", "32"));
    const auto turned = owned->sample_axes(si, 3, {chart_right, local, world, g}, limits, memory);
    REQUIRE(turned);
    CHECK(std::get<CasePgProjectedAxisValue>(turned->value).sine.lower == q("-1"));
    const std::array chart_diagonal{q("1", "8"), q("0"), q("3", "4")};
    const auto diagonal =
        owned->sample_axes(si, 3, {chart_diagonal, local, world, g}, limits, memory);
    REQUIRE(diagonal);
    const auto& diagonal_axis = std::get<CasePgProjectedAxisValue>(diagonal->value);
    const auto half = q("1", "2");
    for (const auto* component : {&diagonal_axis.cosine, &diagonal_axis.sine}) {
        REQUIRE(component->lower.sign() > 0);
        const auto low_squared =
            fixture::take(component->lower.multiply(component->lower, limits.algebra, memory));
        const auto high_squared =
            fixture::take(component->upper.multiply(component->upper, limits.algebra, memory));
        CHECK(fixture::take(low_squared.compare(half, limits.algebra, memory)) <= 0);
        CHECK(fixture::take(high_squared.compare(half, limits.algebra, memory)) >= 0);
    }

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
TEST_CASE("Fixed shell expressions certify their domain and preserve analytical field errors",
          "[FD-12][MAT-D06][V-MAT][negative][verification]") {
    auto input = fixture::root("structural_static");
    input.field("units").set("angle", fixture::json(R"("deg")"));
    auto& layer =
        input.field("model").field("shells").elements.front().field("layers").elements.front();
    layer.set("initial_green", fixture::json(R"({
      "kind":"expression",
      "expression":{
        "domain":{"kind":"reference","coordinates":"chart","lower":[0,0,-1],"upper":[1,1,1]},
        "nodes":[
          {"op":"variable","axis":0},
          {"op":"exp","operand":0},
          {"op":"log","operand":1},
          {"op":"sqrt","operand":0},
          {"op":"add","first":2,"second":3},
          {"op":"constant","value":2,"dimensions":[0,0,0,0]},
          {"op":"integer_power","operand":5,"exponent":-1},
          {"op":"multiply","first":4,"second":6},
          {"op":"constant","value":90,"dimensions":[0,0,0,1]},
          {"op":"sin","operand":8},
          {"op":"cos","operand":8},
          {"op":"constant","value":1,"dimensions":[0,0,0,0]},
          {"op":"negate","operand":11},
          {"op":"add","first":9,"second":12},
          {"op":"subtract","first":10,"second":10}
        ],
        "outputs":[7,10,13,14,14,14]
      },
      "error":[0.001,0,0,0,0,0]
    })"));
    layer.set("axis", fixture::json(R"({
      "kind":"projected",
      "direction":{"kind":"constant","value":[1,0,0]},
      "angle":{"kind":"expression",
        "expression":{
          "domain":{"kind":"reference","coordinates":"chart","lower":[0,0,-1],"upper":[1,1,1]},
          "nodes":[{"op":"constant","value":90,"dimensions":[0,0,0,1]}],
          "outputs":[0]
        },
        "error":[0]
      }
    })"));
    core::ResourcePlan memory{budget};
    const pg::PgReferenceLimits limits{algebra(), 200'000, 100'000, {}};
    const auto q = [&](std::string_view n, std::string_view d = "1") {
        return fixture::take(pg::ExactRational::from_integer_pair(n, d, limits.algebra, memory));
    };
    const auto source = materialize(input);
    const auto shell = fixture::take(CasePgShellMaterials::prepare(source, 1, limits, memory));
    const auto one = q("1");
    const pg::PgReferenceVector origin;
    const std::array chart{q("1", "4"), q("0"), q("0")};
    const CasePgMaterialCoordinates point{chart, origin, origin, one};
    const auto sampled =
        fixture::take(shell.sample_initial_green(source, 0, point, limits, memory));
    const auto contains = [&](std::size_t i, const pg::ExactRational& expected) {
        const auto& value = sampled.mandel[i];
        CHECK(value.radical.sign() == 0);
        const auto lo = fixture::take(
            value.rational.subtract(sampled.declared_error[i], limits.algebra, memory));
        const auto hi =
            fixture::take(value.rational.add(sampled.declared_error[i], limits.algebra, memory));
        CHECK(fixture::take(lo.compare(expected, limits.algebra, memory)) <= 0);
        CHECK(fixture::take(hi.compare(expected, limits.algebra, memory)) >= 0);
    };
    contains(0, q("3", "8"));
    for (std::size_t i = 1; i < 6; ++i)
        contains(i, q("0"));
    REQUIRE(fixture::take(
                sampled.declared_error[0].compare(q("1", "1000"), limits.algebra, memory)) >= 0);
    CHECK(fixture::take(
              sampled.declared_error[0].compare(q("11", "10000"), limits.algebra, memory)) < 0);
    const auto limit = q("1", "100000000000000000000");
    CHECK(fixture::take(sampled.declared_error[1].compare(limit, limits.algebra, memory)) < 0);
    const auto orientation = fixture::take(shell.sample_axes(source, 0, point, limits, memory));
    const auto& rotation = std::get<CasePgProjectedAxisValue>(orientation.value);
    CHECK(rotation.cosine.lower.sign() <= 0);
    CHECK(rotation.cosine.upper.sign() >= 0);
    const auto cosine_width = fixture::take(
        rotation.cosine.upper.subtract(rotation.cosine.lower, limits.algebra, memory));
    CHECK(fixture::take(cosine_width.compare(limit, limits.algebra, memory)) < 0);
    CHECK(fixture::take(rotation.sine.lower.compare(one, limits.algebra, memory)) <= 0);
    CHECK(fixture::take(rotation.sine.upper.compare(one, limits.algebra, memory)) >= 0);
    const std::array outside{q("5", "4"), q("0"), q("0")};
    const auto rejected =
        shell.sample_initial_green(source, 0, {outside, origin, origin, one}, limits, memory);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::InvalidInput);

    layer.set("initial_green", fixture::json(R"({
      "kind":"expression",
      "expression":{
        "domain":{"kind":"reference","coordinates":"chart","lower":[0,0,-1],"upper":[1,1,1]},
        "nodes":[{"op":"variable","axis":0},
                 {"op":"constant","value":1,"dimensions":[0,0,0,0]},
                 {"op":"divide","first":1,"second":0}],
        "outputs":[2,2,2,2,2,2]
      },
      "error":[0,0,0,0,0,0]
    })"));
    const auto singular_source = materialize(input);
    const auto singular_shell =
        fixture::take(CasePgShellMaterials::prepare(singular_source, 1, limits, memory));
    // The requested point is above zero, but the declared box includes zero.
    const auto singular =
        singular_shell.sample_initial_green(singular_source, 0, point, limits, memory);
    REQUIRE_FALSE(singular);
    CHECK(singular.error().code == core::FailureCode::PrecisionUnresolved);
}
} // namespace
