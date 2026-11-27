#include "pg_conditioning_test_support.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <uvlas/structure/pg_shell_material.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::structure;
using namespace pg_conditioning_test;
using Json = nlohmann::json;

ExactRational word(const Json& value, const ExactLimits& limits, core::ResourcePlan& plan) {
    return take(ExactRational::from_binary64(value.get<double>(), limits, plan));
}

Q decimal(const Json& value, const ExactLimits& limits, core::ResourcePlan& plan) {
    const auto text = value.is_string() ? value.get<std::string>() : value.dump();
    return q(take(ExactRational::from_decimal(text, limits, plan)));
}

PgShellElasticLaw law(const Json& input, const PgReferenceLimits& limits,
                      core::ResourcePlan& plan) {
    const auto rho = word(input.at("rho"), limits.algebra, plan);
    const auto read = [&](const char* key) { return word(input.at(key), limits.algebra, plan); };
    if (input.at("kind") == "isotropic") {
        const auto young = read("E"), poisson = read("nu");
        return take(
            PgShellElasticLaw::prepare(PgIsotropicElasticity{young, poisson}, rho, limits, plan));
    }
    if (input.at("kind") == "orthotropic") {
        const std::array young{read("E1"), read("E2"), read("E3")};
        const std::array poisson{read("nu12"), read("nu13"), read("nu23")};
        const std::array shear{read("G23"), read("G13"), read("G12")};
        return take(PgShellElasticLaw::prepare(PgOrthotropicElasticity{young, poisson, shear}, rho,
                                               limits, plan));
    }
    REQUIRE(input.at("kind") == "effective_3d");
    std::array<ExactRational, 36> stiffness;
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            stiffness[6 * i + j] = word(input.at("C").at(i).at(j), limits.algebra, plan);
    return take(PgShellElasticLaw::prepare(
        PgEffectiveElasticity{stiffness, PgMaterialComponents::Mandel}, rho, limits, plan));
}

PgReferenceGeometry geometry(Fixture& fixture, const Json& coefficients, const Json& materials,
                             const PgReferenceLimits& limits, core::ResourcePlan& plan) {
    std::vector<ExactRational> values;
    std::vector<ExactInputTerm> terms;
    values.reserve(coefficients.size());
    for (const auto& term : coefficients) {
        values.push_back(fixture.r(term.at(2).get<std::string>(), term.at(3).get<std::string>()));
        terms.push_back(
            {{term.at(0).get<core::Index>(), term.at(1).get<core::Index>()}, values.back()});
    }
    const auto prep_limits = partition_test::exact_limits();
    const auto stage = []<class Value>(core::Result<Value> result, std::string_view name) {
        INFO(name);
        const auto code = result ? 0 : static_cast<int>(result.error().code);
        CAPTURE(code);
        REQUIRE(result);
        return std::move(*result);
    };
    const auto height = stage(ExactChartPolynomial::make(terms, prep_limits, plan), "height");
    const auto zero = fixture.r("0"), one = fixture.r("1");
    fixture.graph =
        stage(PgGraph::prepare({zero, one, zero, height, 1, PgGraphConvention::Version1},
                               {prep_limits, 64, {}}, plan),
              "graph");
    // H=(1+s/10+t/20)/20 m from the frozen source.
    fixture.heights[0] = fixture.r("1", "20");
    fixture.heights[1] = fixture.r("11", "200");
    fixture.heights[2] = fixture.r("23", "400");
    fixture.heights[3] = fixture.r("21", "400");
    std::vector<std::array<ExactRational, 2>> endpoints;
    std::vector<PgPartitionLayer> laminate;
    endpoints.reserve(materials.size());
    laminate.reserve(materials.size());
    for (std::size_t i = 0; i < materials.size(); ++i) {
        const auto& item = materials.at(i);
        auto lo = take(
            ExactRational::from_decimal(item.at("lower").get<std::string>(), limits.algebra, plan));
        auto hi = take(
            ExactRational::from_decimal(item.at("upper").get<std::string>(), limits.algebra, plan));
        endpoints.push_back({std::move(lo), std::move(hi)});
        laminate.emplace_back(core::LayerId{i + 1}, endpoints.back()[0], endpoints.back()[1]);
    }
    core::Workspace scratch{core::PlannedBuffer<std::byte>{}};
    auto dictionary = stage(
        PgCoefficientDictionary::prepare(
            fixture.chart_graph(), {prep_limits, 2000, 20'000'000, 100'000, {}}, scratch, plan),
        "dictionary");
    auto request = fixture.request();
    request.layers = laminate;
    auto partition =
        stage(PgIntegrationPartition::prepare(request, partition_test::limits(), scratch, plan),
              "partition");
    const std::array rotation{fixture.r("1"), fixture.r("0"), fixture.r("0"),
                              fixture.r("0"), fixture.r("1"), fixture.r("0"),
                              fixture.r("0"), fixture.r("0"), fixture.r("1")};
    const std::array<ExactRational, 3> origin;
    return stage(PgReferenceGeometry::prepare({std::move(fixture.chart_graph()),
                                               std::move(dictionary), std::move(partition)},
                                              {rotation, origin, one}, limits, plan),
                 "physical reference");
}

PgReferenceValue evaluate(Fixture& fixture, const PgReferenceGeometry& geometry, const Json& input,
                          const PgReferenceLimits& limits, core::ResourcePlan& plan) {
    const auto& coordinates = input.at("coordinates");
    const auto s = decimal(coordinates.at(0), limits.algebra, plan);
    const auto t = decimal(coordinates.at(1), limits.algebra, plan);
    const auto zeta = decimal(coordinates.at(2), limits.algebra, plan);
    const auto layer = input.at("layer").get<core::Index>();
    const auto& partition = take(geometry.source()).partition;
    const auto vertices = take(partition.vertices());
    const auto triangles = take(partition.triangles());
    const auto volumes = take(partition.volumes());
    const auto layers = take(partition.layers());
    const auto material = static_cast<std::size_t>(layer);
    REQUIRE(material < layers.size());
    const Q depth = (zeta - q(layers[material].lower)) /
                    (q(layers[material].upper) - q(layers[material].lower));
    REQUIRE(depth >= 0);
    REQUIRE(depth <= 1);
    for (std::size_t i = 0; i < volumes.size(); ++i) {
        if (volumes[i].layer != layer)
            continue;
        const auto& corners = triangles[static_cast<std::size_t>(volumes[i].triangle)].vertices;
        const auto& a = vertices[static_cast<std::size_t>(corners[0])];
        const auto& b = vertices[static_cast<std::size_t>(corners[1])];
        const auto& c = vertices[static_cast<std::size_t>(corners[2])];
        const Q ds = s - q(a.s), dt = t - q(a.t);
        const Q bs = q(b.s) - q(a.s), bt = q(b.t) - q(a.t);
        const Q cs = q(c.s) - q(a.s), ct = q(c.t) - q(a.t);
        const Q det = bs * ct - bt * cs;
        const Q u = (ds * ct - dt * cs) / det, v = (bs * dt - bt * ds) / det;
        if (u < 0 || v < 0 || u + v > 1)
            continue;
        const auto convert = [&](const Q& value) {
            return rat(fixture, numerator(value), denominator(value));
        };
        const auto ur = convert(u), vr = convert(v), dr = convert(depth);
        return take(
            geometry.evaluate({PgReferenceVolumePoint{static_cast<core::Index>(i), ur, vr, dr}},
                              std::nullopt, limits, plan));
    }
    FAIL("Frozen V-MAT point is outside the physical partition");
    throw std::runtime_error("Missing V-MAT point");
}

void check(const PgMaterialEnclosure& actual, const Json& expected) {
    const auto value = expected.get<double>();
    const auto budget = 1e-10 * std::max(1.0, std::abs(value));
    CHECK(q(actual.lower).convert_to<double>() >= value - budget);
    CHECK(q(actual.upper).convert_to<double>() <= value + budget);
}

TEST_CASE("Export selected PG coefficients for the frozen V-MAT transfer", "[.fixture][V-MAT]") {
    std::ifstream input{UVLAS_PG_FROZEN_MATERIAL_PROBES};
    REQUIRE(input);
    const auto frozen = Json::parse(input);
    const auto limits = query_limits();
    core::ResourcePlan plan{profile};
    Fixture fixture;
    const auto patch = geometry(fixture, frozen.at("height_coefficients"),
                                frozen.at("material_input"), limits, plan);
    const auto source = take(patch.source());
    const PgRankPreparationLimits rank_limits{profile, steps, 128, {}};
    core::Workspace scratch{core::PlannedBuffer<std::byte>{}};
    const auto rank = take(PgRankPreparation::prepare(source.dictionary, source.partition,
                                                      rank_limits, scratch, plan));
    const auto rows = take(rank.rows());
    const auto count = static_cast<core::Index>(take(rank.preparation()).rank);
    for (core::Index column = 0; column < count; ++column) {
        for (const auto index : take(rank.basis_indices({column}))) {
            const auto& row = rows[index];
            const auto value = take(rank.basis({column}, index, limits.algebra, plan));
            std::cout << "PG_COEFFICIENT " << column << ' ' << static_cast<int>(row.field) << ' '
                      << row.cell.first << ' ' << row.cell.second << ' ' << row.powers.first << ' '
                      << row.powers.second << ' ' << q(value) << '\n';
        }
    }
    REQUIRE(std::cout.good());
}

TEST_CASE("Frozen V-MAT curved probes reproduce each full 3D material and stress measure",
          "[MAT-D08][MAT-D10][V-MAT][verification]") {
    std::ifstream input{UVLAS_PG_FROZEN_MATERIAL_PROBES};
    REQUIRE(input);
    const auto frozen = Json::parse(input);
    REQUIRE(frozen.at("schema") == "uvlas.pg-frozen-material-probes/1");
    auto limits = query_limits();
    core::ResourcePlan plan{profile};
    Fixture fixture;
    const auto& materials = frozen.at("material_input");
    const auto patch = geometry(fixture, frozen.at("height_coefficients"), materials, limits, plan);
    const auto one = fixture.r("1"), zero = fixture.r("0");
    REQUIRE(frozen.at("probes").size() == 7);
    for (const auto& probe : frozen.at("probes")) {
        CAPTURE(probe.at("coordinates").dump(), probe.at("layer").dump());
        const auto layer = probe.at("layer").get<std::size_t>();
        const auto& material = materials.at(layer);
        const auto& law_input = material.at("material");
        const auto elastic = law(law_input, limits, plan);
        const auto reference = evaluate(fixture, patch, probe, limits, plan);
        CHECK(reference.layer == core::LayerId{layer + 1});
        const auto root = std::sqrt(q(reference.g).convert_to<double>());
        std::array<ExactRational, 3> fibre;
        for (std::size_t i = 0; i < 3; ++i) {
            fibre[i] = word(material.at("fiber_global").at(i), limits.algebra, plan);
            const auto& actual = reference.position_metres[i];
            const auto position = q(actual.rational).convert_to<double>() +
                                  root * q(actual.radical).convert_to<double>();
            CHECK(std::abs(position - probe.at("reference_position_metres").at(i).get<double>()) <
                  1e-12);
        }
        std::array<PgReferenceScalar, 6> initial;
        const std::array<ExactRational, 6> error;
        for (std::size_t i = 0; i < 6; ++i)
            initial[i].rational = word(law_input.at("e0").at(i), limits.algebra, plan);
        PgShellTensorEnclosure gradient;
        const auto& frozen_gradient = probe.at("deformation_gradient");
        for (std::size_t i = 0; i < 9; ++i) {
            const auto& value = frozen_gradient.at(i / 3).at(i % 3);
            gradient[i] = {word(value, limits.algebra, plan), word(value, limits.algebra, plan)};
        }
        const auto result = take(enclose_pg_shell_material_point(
            elastic, {one, initial, error}, PgShellProjectedAxes{reference, {fibre, one, zero}},
            gradient, limits, plan));
        check(result.material.energy_joules_per_cubic_metre,
              probe.at("energy_joules_per_cubic_metre"));
        for (std::size_t i = 0; i < 9; ++i) {
            const auto expected = [&](const char* field) -> const Json& {
                return probe.at(field).at(i / 3).at(i % 3);
            };
            check(result.reference_axes[i], expected("reference_axes"));
            check(result.second_pk_world_pascal[i], expected("second_pk_world_pascal"));
            check(result.first_pk_world_pascal[i], expected("first_pk_world_pascal"));
            check(result.cauchy_world_pascal[i], expected("cauchy_world_pascal"));
        }
    }
}
} // namespace
