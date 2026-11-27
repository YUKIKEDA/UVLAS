#include "pg_conditioned_volume_row.hpp"
#include "pg_conditioning_test_support.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <uvlas/structure/pg_discrete_gram.hpp>
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

struct ProbeLocation {
    core::Index volume;
    ExactRational u, v, depth;

    [[nodiscard]] PgReferenceVolumePoint point() const& noexcept {
        return {volume, u, v, depth};
    }
    [[nodiscard]] PgReferenceVolumePoint point() const&& = delete;
};

ProbeLocation locate_chart(Fixture& fixture, const PgIntegrationPartition& partition,
                           const std::array<Q, 3>& coordinates, core::Index layer) {
    const auto& [s, t, zeta] = coordinates;
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
        return {static_cast<core::Index>(i), convert(u), convert(v), convert(depth)};
    }
    FAIL("Frozen V-MAT point is outside the physical partition");
    throw std::runtime_error("Missing V-MAT point");
}

ProbeLocation locate(Fixture& fixture, const PgIntegrationPartition& partition, const Json& input,
                     const PgReferenceLimits& limits, core::ResourcePlan& plan) {
    const auto& coordinates = input.at("coordinates");
    const std::array converted{decimal(coordinates.at(0), limits.algebra, plan),
                               decimal(coordinates.at(1), limits.algebra, plan),
                               decimal(coordinates.at(2), limits.algebra, plan)};
    return locate_chart(fixture, partition, converted, input.at("layer").get<core::Index>());
}

template <class Value>
Value checked(core::Result<Value> result, std::string_view operation) {
    INFO(operation);
    CAPTURE(result ? 0 : static_cast<int>(result.error().code));
    REQUIRE(result);
    return std::move(*result);
}

PgPreparedPatch prepare_patch(PgReferenceGeometry&& geometry, Fixture& fixture,
                              const PgReferenceLimits& limits, core::Workspace& scratch,
                              core::ResourcePlan& plan) {
    const auto source = take(geometry.source());
    const PgRankPreparationLimits rank_limits{profile, steps, 128, {}};
    auto rank = checked(
        PgRankPreparation::prepare(source.dictionary, source.partition, rank_limits, scratch, plan),
        "selected rank");
    auto basis = checked(
        PgReferenceBasis::prepare(std::move(geometry), std::move(rank), rank_limits, scratch, plan),
        "reference basis");
    // An overdetermined positive sampling rule qualifies the basis independently
    // of the rule used later for elastic or mass actions.
    const std::array planar{fixture.r("1", "12"), fixture.r("3", "12"), fixture.r("5", "12"),
                            fixture.r("7", "12"), fixture.r("9", "12"), fixture.r("11", "12")};
    const auto sixth = fixture.r("1", "6"), third = fixture.r("1", "3");
    const std::array xy{
        PgLineQuadratureNode{planar[0], sixth}, PgLineQuadratureNode{planar[1], sixth},
        PgLineQuadratureNode{planar[2], sixth}, PgLineQuadratureNode{planar[3], sixth},
        PgLineQuadratureNode{planar[4], sixth}, PgLineQuadratureNode{planar[5], sixth}};
    const std::array depths{fixture.r("1", "6"), fixture.r("1", "2"), fixture.r("5", "6")};
    const std::array z{PgLineQuadratureNode{depths[0], third},
                       PgLineQuadratureNode{depths[1], third},
                       PgLineQuadratureNode{depths[2], third}};
    const PgDiscreteGramLimits gram_limits{limits, steps};
    auto gram = checked(PgDiscreteGram::prepare_product(std::move(basis), {xy, xy, z}, {352, 272},
                                                        gram_limits, scratch, plan),
                        "curved geometric Gram");
    auto conditioning = checked(
        std::move(gram).generate_and_certify({512, 216}, {profile, steps, {}}, scratch, plan),
        "Gram conditioning");
    return checked(PgPreparedPatch::prepare(std::move(conditioning), limits, scratch, plan),
                   "prepared patch");
}

struct TrialState {
    std::vector<PgCartesianExactPosition> current;
    std::vector<PgCartesianExactPosition> displacement;
};

TrialState transferred_state(const Json& coefficients, const PgPreparedPatch& patch,
                             Fixture& fixture, const PgReferenceLimits& limits,
                             core::ResourcePlan& plan) {
    const auto reference = patch.reference_coordinates();
    const auto count = reference.size();
    REQUIRE(coefficients.size() == count);
    TrialState state{std::vector<PgCartesianExactPosition>(count),
                     std::vector<PgCartesianExactPosition>(count)};
    std::vector<std::array<Q, 3>> solution(count);
    const I dyadic_scale = I{1} << 256U;
    const auto convert = [&](const Q& value) {
        return checked(ExactRational::from_integer_pair(
                           numerator(value).str(), denominator(value).str(), limits.algebra, plan),
                       "Cartesian coefficient");
    };
    for (std::size_t row = count; row-- > 0;) {
        std::array<Q, 3> right;
        for (std::size_t axis = 0; axis < 3; ++axis)
            right[axis] = Q{coefficients.at(row).at(axis).get<std::string>()};
        const PgRankBasisColumn i{static_cast<core::Index>(row)};
        const auto transform = [&](std::size_t j) -> ExactRational {
            return take(patch.conditioning().transform_coefficient(i, {static_cast<core::Index>(j)},
                                                                   limits.algebra, plan));
        };
        for (std::size_t j = row + 1; j < count; ++j) {
            const auto coefficient = transform(j);
            for (std::size_t axis = 0; axis < 3; ++axis)
                right[axis] -= q(coefficient) * solution[j][axis];
        }
        const auto diagonal = transform(row);
        REQUIRE(diagonal.sign() > 0);
        for (std::size_t axis = 0; axis < 3; ++axis) {
            // Quantizing each unknown avoids the unbounded denominators of an
            // exact triangular solve. The residual test uses the actual T entry.
            const Q scaled = right[axis] * dyadic_scale / q(diagonal);
            const I units = numerator(scaled) / denominator(scaled);
            solution[row][axis] = Q{units} / dyadic_scale;
            CHECK(abs(right[axis] - q(diagonal) * solution[row][axis]) <
                  Q{"1/10000000000000000000000000000000000000000"});
            state.displacement[row][axis] = rat(fixture, units, dyadic_scale);
            state.current[row][axis] = convert(q(reference[row][axis]) + solution[row][axis]);
        }
    }
    return state;
}

void check(const PgMaterialEnclosure& actual, const Json& expected, double tolerance = 1e-10) {
    const auto value = expected.get<double>();
    const auto budget = tolerance * std::max(1.0, std::abs(value));
    CHECK(q(actual.lower).convert_to<double>() >= value - budget);
    CHECK(q(actual.upper).convert_to<double>() <= value + budget);
}

struct MaterialInput {
    PgShellElasticLaw elastic;
    std::array<ExactRational, 3> fibre;
    std::array<PgReferenceScalar, 6> initial;
};

MaterialInput read_material(const Json& input, const PgReferenceLimits& limits,
                            core::ResourcePlan& plan) {
    const auto& law_input = input.at("material");
    MaterialInput material{law(law_input, limits, plan), {}, {}};
    for (std::size_t axis = 0; axis < material.fibre.size(); ++axis)
        material.fibre[axis] = word(input.at("fiber_global").at(axis), limits.algebra, plan);
    for (std::size_t i = 0; i < material.initial.size(); ++i)
        material.initial[i].rational = word(law_input.at("e0").at(i), limits.algebra, plan);
    return material;
}

void check_material_response(const PgShellMaterialPointResponse& actual, const Json& probe) {
    check(actual.material.energy_joules_per_cubic_metre, probe.at("energy_joules_per_cubic_metre"));
    for (std::size_t i = 0; i < 9; ++i) {
        const auto expected = [&](const char* field) -> const Json& {
            return probe.at(field).at(i / 3).at(i % 3);
        };
        check(actual.reference_axes[i], expected("reference_axes"));
        check(actual.second_pk_world_pascal[i], expected("second_pk_world_pascal"));
        check(actual.first_pk_world_pascal[i], expected("first_pk_world_pascal"));
        check(actual.cauchy_world_pascal[i], expected("cauchy_world_pascal"));
    }
}

void export_volume_midpoints(const PgIntegrationPartition& partition, std::ostream& output) {
    const auto vertices = take(partition.vertices());
    const auto triangles = take(partition.triangles());
    const auto layers = take(partition.layers());
    for (const auto& volume : take(partition.volumes())) {
        const auto& triangle = triangles[static_cast<std::size_t>(volume.triangle)];
        const auto& a = vertices[static_cast<std::size_t>(triangle.vertices[0])];
        const auto& b = vertices[static_cast<std::size_t>(triangle.vertices[1])];
        const auto& c = vertices[static_cast<std::size_t>(triangle.vertices[2])];
        const auto& layer = layers[static_cast<std::size_t>(volume.layer)];
        const Q area =
            ((q(b.s) - q(a.s)) * (q(c.t) - q(a.t)) - (q(b.t) - q(a.t)) * (q(c.s) - q(a.s))) / 2;
        const Q weight = area * (q(layer.upper) - q(layer.lower));
        REQUIRE(weight > 0);
        output << "PG_VOLUME " << volume.layer << ' ' << (q(a.s) + q(b.s) + q(c.s)) / 3 << ' '
               << (q(a.t) + q(b.t) + q(c.t)) / 3 << ' ' << (q(layer.lower) + q(layer.upper)) / 2
               << ' ' << weight << '\n';
    }
}

void check_relative_enclosure(const Q& lower, const Q& upper, const Json& expected) {
    const Q reference{expected.get<double>()};
    REQUIRE(reference != 0);
    const Q tolerance = abs(reference) / 100'000'000;
    CHECK(lower <= upper);
    CHECK(lower >= reference - tolerance);
    CHECK(upper <= reference + tolerance);
}

void check_action(std::span<const PgShellCartesianEnclosure> action,
                  std::span<const PgCartesianExactPosition> displacement, const Json& expected) {
    REQUIRE(action.size() == displacement.size());
    Q lo{}, hi{};
    for (std::size_t i = 0; i < action.size(); ++i) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const Q a = q(action[i][axis].lower) * q(displacement[i][axis]);
            const Q b = q(action[i][axis].upper) * q(displacement[i][axis]);
            lo += std::min(a, b);
            hi += std::max(a, b);
        }
    }
    check_relative_enclosure(lo, hi, expected);
}

void compare_mass_row(const PgPreparedPatch& patch, const PgReferenceVolumePoint& point,
                      const PgReferenceLimits& limits, core::ResourcePlan& plan) {
    const auto slim =
        checked(detail::PgConditionedVolumeRow::evaluate(patch.conditioning(), point, limits, plan),
                "conditioned volume row");
    const PgReferencePoint query{point};
    const auto full = checked(patch.evaluate_reference(std::span{&query, 1}, limits, plan),
                              "complete conditioned reference row");
    const auto& reference = full.points.values().front();
    REQUIRE(slim.reference.layer == reference.layer);
    const auto& measure = std::get<PgReferenceVolumeMeasure>(slim.reference.measure).jacobian;
    const auto& full_measure = std::get<PgReferenceVolumeMeasure>(reference.measure).jacobian;
    CHECK(q(slim.reference.g) == q(reference.g));
    CHECK(q(measure.rational) == q(full_measure.rational));
    CHECK(q(measure.radical) == q(full_measure.radical));
    const auto row = full.columns.values();
    REQUIRE(slim.values.size() == row.size());
    for (std::size_t j = 0; j < row.size(); ++j) {
        CHECK(q(slim.values.values()[j].rational) == q(row[j].value.rational));
        CHECK(q(slim.values.values()[j].radical) == q(row[j].value.radical));
    }
    CHECK(slim.work.allocated_owned_payload_bytes < full.call_info.allocated_owned_payload_bytes);
}

void check_volume(const PgPreparedPatch& patch, PgPreparedPatchIdentity patch_id,
                  const TrialState& state, const Json& material_input, const Json& expected,
                  Fixture& fixture, const PgReferenceLimits& limits, core::ResourcePlan& plan) {
    const auto& partition = take(take(patch.conditioning().reference_basis())->source()).partition;
    const auto layers = take(partition.layers());
    const auto volumes = take(partition.volumes());
    REQUIRE(volumes.size() == 24);
    REQUIRE(volumes.size() == expected.at("points").get<std::size_t>());
    REQUIRE(material_input.size() == layers.size());
    std::vector<MaterialInput> materials;
    std::vector<ExactRational> densities;
    materials.reserve(layers.size());
    densities.reserve(layers.size());
    for (const auto& entry : material_input) {
        materials.push_back(read_material(entry, limits, plan));
        densities.push_back(word(entry.at("material").at("rho"), limits.algebra, plan));
    }
    const auto one = fixture.r("1"), zero = fixture.r("0");
    const auto third = fixture.r("1", "3"), half = fixture.r("1", "2");
    const std::array<ExactRational, 6> error;
    std::vector<PgShellLayerDensity> mass_layers;
    mass_layers.reserve(layers.size());
    for (std::size_t i = 0; i < layers.size(); ++i)
        mass_layers.push_back({layers[i].id, densities[i]});
    std::vector<PgShellVolumeSample> samples;
    std::vector<PgDiscreteVolumePoint> mass_rule;
    samples.reserve(volumes.size());
    mass_rule.reserve(volumes.size());
    std::array<std::size_t, 3> layer_counts{};
    for (std::size_t i = 0; i < volumes.size(); ++i) {
        const auto layer = static_cast<std::size_t>(volumes[i].layer);
        ++layer_counts.at(layer);
        const auto& material = materials.at(layer);
        const PgDiscreteVolumePoint midpoint{{static_cast<core::Index>(i), third, third, half},
                                             half};
        samples.push_back({midpoint,
                           layers[layer].id,
                           material.elastic,
                           {one, material.initial, error},
                           PgShellProjectedFibre{material.fibre, one, zero}});
        mass_rule.push_back(midpoint);
    }
    REQUIRE(layer_counts == std::array<std::size_t, 3>{8, 8, 8});
    for (const auto i : {std::size_t{0}, volumes.size() / 2, volumes.size() - 1})
        compare_mass_row(patch, mass_rule[i].point, limits, plan);
    const auto width = fixture.r("1", "1000000000000000000000");
    const PgShellPatchDirection displacement{patch_id, state.displacement};
    auto assembly_limits = limits;
    // The two outward rational endpoints per row carry limb-product charges.
    // Budget those charges on this full rule separately from point evaluation.
    assembly_limits.maximum_scalar_operations = 2'000'000'000;
    assembly_limits.algebra.external_memory = ExactEstimateAdmission{profile};
    assembly_limits.algebra.maximum_value_bits = 16'384;
    assembly_limits.algebra.maximum_work_bits = 65'536;
    core::ResourcePlan elastic_plan{profile};
    auto elastic_action =
        assemble_pg_shell_elastic_volume(patch, {patch_id, state.current}, samples, displacement,
                                         {112, width, width, width}, assembly_limits, elastic_plan);
    CAPTURE(elastic_plan.allocated_bytes());
    const auto elastic = checked(std::move(elastic_action), "frozen trial elastic assembly");
    core::ResourcePlan mass_plan{profile};
    auto mass_action =
        apply_pg_shell_consistent_mass(patch, mass_rule, mass_layers, displacement,
                                       {112, width, width}, assembly_limits, mass_plan);
    CAPTURE(mass_plan.allocated_bytes());
    const auto mass = checked(std::move(mass_action), "frozen trial consistent mass");
    CAPTURE(elastic.work.scalar_operations, elastic.work.partition_visits,
            mass.work.scalar_operations, mass.work.partition_visits);
    CHECK(elastic.work.allocated_owned_payload_bytes == elastic_plan.allocated_bytes());
    CHECK(mass.work.allocated_owned_payload_bytes == mass_plan.allocated_bytes());
    const auto check_total = [&](const PgMaterialEnclosure& actual, std::string_view field) {
        INFO(field);
        check_relative_enclosure(q(actual.lower), q(actual.upper), expected.at(field));
    };
    check_total(elastic.energy_joules, "energy_joules");
    check_total(mass.mass_kg, "mass_kg");
    check_action(elastic.internal_force_newtons.values(), state.displacement,
                 expected.at("displacement_internal_work_joules"));
    check_action(elastic.tangent_action_newtons.values(), state.displacement,
                 expected.at("displacement_tangent_quadratic_joules"));
    check_action(mass.action_kg_metres.values(), state.displacement,
                 expected.at("displacement_mass_quadratic_kg_metres_squared"));
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
    export_volume_midpoints(source.partition, std::cout);
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
        const auto material = read_material(materials.at(layer), limits, plan);
        const auto location = locate(fixture, take(patch.source()).partition, probe, limits, plan);
        const auto reference = take(patch.evaluate({location.point()}, std::nullopt, limits, plan));
        CHECK(reference.layer == core::LayerId{layer + 1});
        const auto root = std::sqrt(q(reference.g).convert_to<double>());
        for (std::size_t i = 0; i < 3; ++i) {
            const auto& actual = reference.position_metres[i];
            const auto position = q(actual.rational).convert_to<double>() +
                                  root * q(actual.radical).convert_to<double>();
            CHECK(std::abs(position - probe.at("reference_position_metres").at(i).get<double>()) <
                  1e-12);
        }
        const std::array<ExactRational, 6> error;
        PgShellTensorEnclosure gradient;
        const auto& frozen_gradient = probe.at("deformation_gradient");
        for (std::size_t i = 0; i < 9; ++i) {
            const auto& value = frozen_gradient.at(i / 3).at(i % 3);
            gradient[i] = {word(value, limits.algebra, plan), word(value, limits.algebra, plan)};
        }
        const auto result = take(enclose_pg_shell_material_point(
            material.elastic, {one, material.initial, error},
            PgShellProjectedAxes{reference, {material.fibre, one, zero}}, gradient, limits, plan));
        check_material_response(result, probe);
    }
}

TEST_CASE("Frozen V-MAT selected coefficients produce the same trial through the prepared patch",
          "[SH-D04][MAT-D08][MAT-D10][V-MAT][verification]") {
    std::ifstream input{UVLAS_PG_FROZEN_MATERIAL_PROBES};
    std::ifstream transfer{UVLAS_PG_FROZEN_MATERIAL_TRANSFER};
    REQUIRE(input);
    REQUIRE(transfer);
    const auto frozen = Json::parse(input);
    const auto selected = Json::parse(transfer);
    REQUIRE(selected.at("schema") == "uvlas.pg-frozen-material-transfer/1");
    REQUIRE(selected.at("archive_sha256") == frozen.at("archive").at("sha256"));
    auto limits = query_limits();
    limits.maximum_scalar_operations = 500'000'000;
    limits.maximum_partition_visits = 250'000'000;
    core::ResourcePlan plan{profile};
    core::Workspace scratch{core::PlannedBuffer<std::byte>{}};
    Fixture fixture;
    const auto& materials = frozen.at("material_input");
    const auto patch =
        prepare_patch(geometry(fixture, frozen.at("height_coefficients"), materials, limits, plan),
                      fixture, limits, scratch, plan);
    const auto state = transferred_state(selected.at("transfer").at("displacement_metres"), patch,
                                         fixture, limits, plan);
    const auto info = take(patch.info());
    const auto& partition = take(take(patch.conditioning().reference_basis())->source()).partition;
    const auto one = fixture.r("1"), zero = fixture.r("0");
    const std::array<ExactRational, 6> error;
    for (const auto& probe : frozen.at("probes")) {
        CAPTURE(probe.at("coordinates").dump(), probe.at("layer").dump());
        const auto location = locate(fixture, partition, probe, limits, plan);
        const auto actual =
            checked(enclose_pg_shell_patch_point(patch, {info.identity, state.current},
                                                 {location.point()}, limits, plan),
                    "frozen trial kinematics");
        CHECK(actual.reference.layer == core::LayerId{probe.at("layer").get<std::size_t>() + 1});
        for (std::size_t i = 0; i < 3; ++i)
            check(actual.current_position_metres[i], probe.at("current_position_metres").at(i),
                  1e-12);
        for (std::size_t i = 0; i < 9; ++i)
            check(actual.deformation_gradient[i],
                  probe.at("deformation_gradient").at(i / 3).at(i % 3));
        const auto material =
            read_material(materials.at(probe.at("layer").get<std::size_t>()), limits, plan);
        const auto result =
            checked(enclose_pg_shell_material_point(
                        material.elastic, {one, material.initial, error},
                        PgShellProjectedFibre{material.fibre, one, zero}, actual, limits, plan),
                    "frozen trial material response");
        check_material_response(result, probe);
    }
    check_volume(patch, info.identity, state, materials, selected.at("midpoint_actions"), fixture,
                 limits, plan);
}
} // namespace
