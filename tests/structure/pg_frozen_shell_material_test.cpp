#include "pg_composed_basis.hpp"
#include "pg_conditioned_volume_row.hpp"
#include "pg_conditioning_test_support.hpp"
#include "pg_material_grid.hpp"

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

struct VolumeRuleSite {
    ProbeLocation location;
    ExactRational weight;
};

[[maybe_unused]] std::vector<VolumeRuleSite>
import_volume_rule(const Json& points, const PgIntegrationPartition& partition, Fixture& fixture) {
    const auto volumes = take(partition.volumes());
    const auto triangles = take(partition.triangles());
    const auto vertices = take(partition.vertices());
    const auto layers = take(partition.layers());
    std::vector<VolumeRuleSite> rule;
    rule.reserve(points.size());
    for (const auto& point : points) {
        const auto& input = point.at("coordinates");
        const std::array chart{Q{input.at(0).get<std::string>()}, Q{input.at(1).get<std::string>()},
                               Q{input.at(2).get<std::string>()}};
        auto site = locate_chart(fixture, partition, chart, point.at("layer").get<core::Index>());
        const auto& volume = volumes[static_cast<std::size_t>(site.volume)];
        const auto& tri = triangles[static_cast<std::size_t>(volume.triangle)];
        const auto& a = vertices[static_cast<std::size_t>(tri.vertices[0])];
        const auto& b = vertices[static_cast<std::size_t>(tri.vertices[1])];
        const auto& c = vertices[static_cast<std::size_t>(tri.vertices[2])];
        const Q determinant =
            (q(b.s) - q(a.s)) * (q(c.t) - q(a.t)) - (q(b.t) - q(a.t)) * (q(c.s) - q(a.s));
        const auto& layer = layers[static_cast<std::size_t>(volume.layer)];
        const Q weight = Q{point.at("weight").get<std::string>()} /
                         (determinant * (q(layer.upper) - q(layer.lower)));
        REQUIRE(weight > 0);
        rule.push_back({std::move(site), rat(fixture, numerator(weight), denominator(weight))});
    }
    std::ranges::stable_sort(rule, {}, [](const auto& site) { return site.location.volume; });
    return rule;
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

void check_same_enclosure(const PgMaterialEnclosure& a, const PgMaterialEnclosure& b) {
    CHECK(q(a.lower) == q(b.lower));
    CHECK(q(a.upper) == q(b.upper));
}

void compare_composed_rows(const PgPreparedPatch& patch, const PgPreparedPatch& other,
                           std::span<const PgDiscreteVolumePoint> rule, Fixture& fixture,
                           const PgReferenceLimits& limits) {
    core::ResourcePlan comparison_plan{profile};
    detail::ReferenceArithmetic preparation{limits, comparison_plan};
    const auto cache = detail::PgComposedBasis::prepare(patch.conditioning(), preparation);
    REQUIRE(cache);
    for (const auto& point : rule) {
        CAPTURE(point.point.volume);
        detail::ReferenceArithmetic evaluation{limits, comparison_plan};
        const auto actual = cache->evaluate(patch.conditioning(), point.point, evaluation);
        REQUIRE(actual);
        const PgReferencePoint query{point.point};
        const auto direct =
            checked(patch.evaluate_reference(std::span{&query, 1}, limits, comparison_plan),
                    "direct conditioned reference row");
        const auto input = actual->columns.values();
        const auto reference = direct.columns.values();
        REQUIRE(input.size() == reference.size());
        const auto compare = [](const PgReferenceScalar& a, const PgReferenceScalar& b) {
            CHECK(q(a.rational) == q(b.rational));
            CHECK(q(a.radical) == q(b.radical));
        };
        CHECK(q(actual->points.values().front().g) == q(direct.points.values().front().g));
        for (std::size_t j = 0; j < input.size(); ++j) {
            compare(input[j].value, reference[j].value);
            for (std::size_t a = 0; a < 3; ++a)
                compare(input[j].physical_gradient_per_metre[a],
                        reference[j].physical_gradient_per_metre[a]);
        }
        CHECK(actual->call_info.allocated_owned_payload_bytes <
              direct.call_info.allocated_owned_payload_bytes);
    }
    const auto small = fixture.r("1", "35184372088833");
    const auto& original = rule.front().point;
    detail::ReferenceArithmetic wide{limits, comparison_plan};
    CHECK_FALSE(cache->evaluate(patch.conditioning(),
                                {original.volume, small, original.v, original.depth}, wide));
    core::ResourcePlan no_storage{0};
    detail::ReferenceArithmetic rejected{limits, no_storage};
    CHECK_FALSE(detail::PgComposedBasis::prepare(patch.conditioning(), rejected));
    CHECK(no_storage.allocated_bytes() == 0);

    detail::ReferenceArithmetic different{limits, comparison_plan};
    const auto stale = detail::partition_entry([&]() -> core::Result<void> {
        static_cast<void>(cache->evaluate(other.conditioning(), original, different));
        return {};
    });
    REQUIRE_FALSE(stale);
    CHECK(stale.error().code == core::FailureCode::StaleGeneration);
    CHECK(different.info().allocated_owned_payload_bytes == 0);

    auto stopped = limits;
    stopped.cancellation = requested_stop_token();
    detail::ReferenceArithmetic cancelled{stopped, comparison_plan};
    const auto failure = detail::partition_entry([&]() -> core::Result<void> {
        static_cast<void>(cache->evaluate(patch.conditioning(), original, cancelled));
        return {};
    });
    REQUIRE_FALSE(failure);
    CHECK(failure.error().code == core::FailureCode::Cancelled);
}

void check_volume(const PgPreparedPatch& patch, PgPreparedPatchIdentity patch_id,
                  const TrialState& state, const Json& material_input, const Json& expected,
                  Fixture& fixture, const PgReferenceLimits& limits, core::ResourcePlan& plan,
                  const PgPreparedPatch* other_owner = nullptr) {
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
    auto elastic_action = assemble_pg_shell_elastic_volume(
        patch, {patch_id, state.current}, samples, displacement,
        {112, width, width, width, PgShellReferencePath::Composed}, assembly_limits, elastic_plan);
    CAPTURE(elastic_plan.allocated_bytes());
    const auto elastic = checked(std::move(elastic_action), "frozen trial elastic assembly");
    if (other_owner) {
        compare_composed_rows(patch, *other_owner, mass_rule, fixture, assembly_limits);
        core::ResourcePlan direct_plan{profile};
        auto direct_limits = assembly_limits;
        direct_limits.maximum_scalar_operations = 550'000'000;
        const auto direct = checked(assemble_pg_shell_elastic_volume(
                                        patch, {patch_id, state.current}, samples, displacement,
                                        {112, width, width, width}, direct_limits, direct_plan),
                                    "direct elastic assembly");
        CAPTURE(elastic.work.scalar_operations, direct.work.scalar_operations,
                elastic_plan.allocated_bytes(), direct_plan.allocated_bytes());
        CHECK(elastic_plan.allocated_bytes() < direct_plan.allocated_bytes());
        check_same_enclosure(elastic.energy_joules, direct.energy_joules);
        for (std::size_t j = 0; j < state.displacement.size(); ++j)
            for (std::size_t a = 0; a < 3; ++a) {
                check_same_enclosure(elastic.internal_force_newtons.values()[j][a],
                                     direct.internal_force_newtons.values()[j][a]);
                check_same_enclosure(elastic.tangent_action_newtons.values()[j][a],
                                     direct.tangent_action_newtons.values()[j][a]);
            }
        std::clog << "VMAT_COMPOSED_24," << elastic_plan.allocated_bytes() << ','
                  << elastic.work.scalar_operations << ',' << elastic.work.partition_visits << ','
                  << direct_plan.allocated_bytes() << ',' << direct.work.scalar_operations << ','
                  << direct.work.partition_visits << '\n';
    }

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
TEST_CASE("The frozen composed jet and elastic action match the original exact path",
          "[.composed-vmat][V-MAT][SH-D04]") {
    const auto frozen = Json::parse(std::ifstream{UVLAS_PG_FROZEN_MATERIAL_PROBES});
    const auto selected = Json::parse(std::ifstream{UVLAS_PG_FROZEN_MATERIAL_TRANSFER});
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
    const auto other =
        prepare_patch(geometry(fixture, frozen.at("height_coefficients"), materials, limits, plan),
                      fixture, limits, scratch, plan);
    check_volume(patch, take(patch.info()).identity, state, materials,
                 selected.at("midpoint_actions"), fixture, limits, plan, &other);
}

TEST_CASE("Elastic volume grid rounds signed endpoints in their outward directions",
          "[V-MAT][SH-D04][negative]") {
    Fixture fixture;
    for (const auto source_bits : {256U, 9000U}) {
        CAPTURE(source_bits);
        const I denominator = (I{1} << source_bits) + 1;
        const I numerator = denominator + 1;
        const auto positive = rat(fixture, numerator, denominator);
        const auto negative = rat(fixture, -numerator, denominator);
        constexpr std::uint32_t grid = 112;
        const Q increment = Q{1} / (I{1} << grid);
        auto limits = query_limits();
        limits.algebra.maximum_work_bits = positive.numerator_bits() + grid;
        core::ResourcePlan plan{profile};
        detail::ReferenceArithmetic math{limits, plan};
        const auto project = [&](const ExactRational& value, detail::GridEnd side) {
            return q(detail::outward_material_grid(value, grid, side, math));
        };
        CHECK(project(positive, detail::GridEnd::Lower) == Q{1});
        CHECK(project(positive, detail::GridEnd::Upper) == Q{1} + increment);
        CHECK(project(negative, detail::GridEnd::Lower) == Q{-1} - increment);
        CHECK(project(negative, detail::GridEnd::Upper) == Q{-1});

        --limits.algebra.maximum_work_bits;
        core::ResourcePlan rejected{profile};
        detail::ReferenceArithmetic short_work{limits, rejected};
        const auto failure = detail::partition_entry([&]() -> core::Result<ExactRational> {
            return detail::outward_material_grid(positive, grid, detail::GridEnd::Upper,
                                                 short_work);
        });
        REQUIRE_FALSE(failure);
        CHECK(failure.error().code == core::FailureCode::ResourceExceeded);
        CHECK(rejected.allocated_bytes() == 0);
    }
}
#ifdef UVLAS_FORMAL_VMAT_ELASTIC_FIXTURE
TEST_CASE("Frozen V-MAT elastic action at 3672 points", "[.high-elastic][V-MAT][SH-D04]") {
    const auto input = Json::parse(std::ifstream{UVLAS_PG_FROZEN_MATERIAL_PROBES});
    const auto transfer = Json::parse(std::ifstream{UVLAS_PG_FROZEN_MATERIAL_TRANSFER});
    const auto cubature = Json::parse(std::ifstream{UVLAS_FORMAL_VMAT_ELASTIC_FIXTURE});
    const auto& archive = input.at("archive").at("sha256");
    REQUIRE(transfer.at("archive_sha256") == archive);
    REQUIRE(cubature.at("schema") == "uvlas.pg-frozen-material-elastic-rule/1");
    REQUIRE(cubature.at("archive_sha256") == archive);
    auto limits = query_limits();
    limits.maximum_scalar_operations = 500'000'000;
    limits.maximum_partition_visits = 250'000'000;
    core::ResourcePlan plan{profile};
    core::Workspace scratch{core::PlannedBuffer<std::byte>{}};
    Fixture fixture;
    const auto& materials = input.at("material_input");
    const auto patch =
        prepare_patch(geometry(fixture, input.at("height_coefficients"), materials, limits, plan),
                      fixture, limits, scratch, plan);
    const auto state = transferred_state(transfer.at("transfer").at("displacement_metres"), patch,
                                         fixture, limits, plan);
    const auto& partition = take(take(patch.conditioning().reference_basis())->source()).partition;
    const auto volumes = take(partition.volumes());
    const auto layers = take(partition.layers());
    const auto rule = import_volume_rule(cubature.at("points"), partition, fixture);
    REQUIRE(rule.size() == 3'672);

    std::vector<MaterialInput> material_inputs;
    material_inputs.reserve(layers.size());
    for (const auto& material : materials)
        material_inputs.push_back(read_material(material, limits, plan));
    REQUIRE(material_inputs.size() == layers.size());

    const auto one = fixture.r("1"), zero = fixture.r("0");
    const std::array<ExactRational, 6> error;
    std::vector<PgShellVolumeSample> samples;
    samples.reserve(rule.size());
    std::array<std::size_t, 3> layer_points{};
    REQUIRE(layers.size() == layer_points.size());
    for (const auto& site : rule) {
        const auto& volume = volumes[static_cast<std::size_t>(site.location.volume)];
        const auto layer = static_cast<std::size_t>(volume.layer);
        ++layer_points.at(layer);
        const auto& material = material_inputs.at(layer);
        samples.push_back({{site.location.point(), site.weight},
                           layers[layer].id,
                           material.elastic,
                           {one, material.initial, error},
                           PgShellProjectedFibre{material.fibre, one, zero}});
    }
    REQUIRE(layer_points == std::array<std::size_t, 3>{1'224, 1'224, 1'224});
    auto assembly_limits = limits;
    // The frozen rule's sampled per-point charges project to about 576 billion.
    assembly_limits.maximum_scalar_operations = 1'000'000'000'000ULL;
    assembly_limits.maximum_partition_visits = 2'000'000'000;
    assembly_limits.algebra.maximum_value_bits = 16'384;
    assembly_limits.algebra.maximum_work_bits = 65'536;
    assembly_limits.algebra.external_memory = ExactEstimateAdmission{profile};

    core::ResourcePlan elastic_plan{profile};
    const auto identity = take(patch.info()).identity;
    const auto width = fixture.r("1", "1000000000000000000000");
    const PgShellPatchDirection displacement{identity, state.displacement};
    auto outcome = assemble_pg_shell_elastic_volume(
        patch, {identity, state.current}, samples, displacement,
        {112, width, width, width, PgShellReferencePath::Composed}, assembly_limits, elastic_plan);
    CAPTURE(elastic_plan.allocated_bytes());
    const auto elastic = checked(std::move(outcome), "high order elastic assembly");
    CHECK(elastic.work.allocated_owned_payload_bytes == elastic_plan.allocated_bytes());
    std::clog << "VMAT_ELASTIC_3672," << elastic_plan.allocated_bytes() << ','
              << elastic.work.scalar_operations << ',' << elastic.work.partition_visits << '\n';

    const auto& expected = cubature.at("reference");
    check_relative_enclosure(q(elastic.energy_joules.lower), q(elastic.energy_joules.upper),
                             expected.at("energy_joules"));
    check_action(elastic.internal_force_newtons.values(), state.displacement,
                 expected.at("displacement_internal_work_joules"));
    check_action(elastic.tangent_action_newtons.values(), state.displacement,
                 expected.at("displacement_tangent_quadratic_joules"));
}
#endif
} // namespace
