#include "aero/geometry_snapshot_fixture.hpp"
#include "aero/global_geometry.hpp"
#include "aero/p1_shared_action.hpp"
#include "coupling/case_uniform_flow.hpp"
#include "coupling/case_uniform_flow_projection.hpp"
#include "coupling/discrete_nonpenetration.hpp"
#include "io/case_schema_test_adapter.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <optional>
#include <uvlas/numerics/dense_solve.hpp>
#include <vector>

namespace {
using namespace uvlas;
namespace a = aero::detail;
namespace ls = numerics::linear;
using Rat = boost::multiprecision::cpp_rational;

constexpr std::size_t chord_count = 2;
constexpr std::size_t band_count = 12;
constexpr std::size_t wing_count = 2;
constexpr std::size_t wing_order = chord_count * band_count;
constexpr std::size_t order = wing_order * wing_count;
constexpr std::size_t points_per_wing = (chord_count + 2) * (band_count + 1);
constexpr double wake_length = 1024.;
constexpr double core_radius = 1e-5;
constexpr auto action_steps = UINT64_C(20000000);
using Point = std::array<double, 3>;
using Matrix = std::array<double, order * order>;

double decimal(const nlohmann::json& number) {
    const auto token = number.get<std::string>();
    std::size_t parsed{};
    const double value = std::stod(token, &parsed);
    REQUIRE(parsed == token.size());
    REQUIRE(std::isfinite(value));
    return value;
}

Point point(const nlohmann::json& value) {
    REQUIRE(value.size() == 3);
    return {decimal(value.at(0)), decimal(value.at(1)), decimal(value.at(2))};
}

template <class T>
T take(core::Result<T> result) {
    if (!result)
        FAIL("Product operation failed with " << static_cast<int>(result.error().code));
    return std::move(*result);
}

template <class T>
T take_global(a::GlobalResult<T> result) {
    if (!result)
        FAIL("Global geometry failed with " << static_cast<int>(result.error().origin.cause.code));
    return std::move(*result);
}

a::P1ComplexWords scalar(double value) {
    return {a::SingleWord{value}, a::SingleWord{0}};
}

struct Enclosure {
    Rat center, radius;
};

Rat word(const a::ScalarWords& value) {
    return std::visit(
        [](const auto& v) -> Rat {
            if constexpr (std::is_same_v<std::decay_t<decltype(v)>, a::SingleWord>)
                return Rat{v.value};
            else
                return Rat{v.high} + Rat{v.low};
        },
        value);
}

Enclosure enclosure(const a::ScaledEnclosure& value) {
    const Rat scale = std::ldexp(1., value.exponent());
    return {word(value.center()) * scale, Rat{value.radius()} * scale};
}

void check(const Enclosure& got, const Rat& expected, const Rat& tolerance) {
    CHECK(abs(got.center - expected) + got.radius <= tolerance);
}

struct Mesh {
    explicit Mesh(const nlohmann::json& coordinates) {
        REQUIRE(coordinates.size() == wing_count);
        const auto point_id = [](std::size_t wing, std::size_t chord, std::size_t span) {
            return a::PointId{1 + wing * points_per_wing + chord * (band_count + 1) + span};
        };
        std::map<std::array<a::PointId, 2>, a::EdgeId> edge_ids;
        const auto edge = [&](a::PointId from, a::PointId to) {
            const std::array ends{std::min(from, to), std::max(from, to)};
            const auto proposed = a::EdgeId{1000 + edges.size()};
            const auto [found, inserted] = edge_ids.try_emplace(ends, proposed);
            if (inserted) {
                edges.push_back({proposed, {from, to}});
                edge_bindings.emplace_back(proposed, a::CanonicalEdge{});
            }
            return found->second;
        };
        for (std::size_t wing = 0; wing < wing_count; ++wing) {
            const auto& sections = coordinates.at(wing);
            REQUIRE(sections.size() == chord_count + 1);
            for (std::size_t chord = 0; chord < chord_count + 2; ++chord) {
                const auto& section = sections.at(std::min(chord, chord_count));
                REQUIRE(section.size() == band_count + 1);
                for (std::size_t span = 0; span <= band_count; ++span) {
                    auto xyz = point(section.at(span));
                    if (chord > chord_count)
                        xyz[0] += wake_length;
                    points.push_back(point_id(wing, chord, span));
                    bindings.push_back(test_geometry::binding());
                    values.insert(values.end(), xyz.begin(), xyz.end());
                }
            }
            for (std::size_t chord = 0; chord <= chord_count; ++chord) {
                for (std::size_t span = 0; span < band_count; ++span) {
                    if (chord == 0)
                        bands.emplace_back(300 + wing * band_count + span);
                    const auto band = bands.at(wing * band_count + span);
                    const auto index = faces.size();
                    const a::PhysicalFaceId face{100 + index};
                    const std::array corners{
                        point_id(wing, chord, span), point_id(wing, chord, span + 1),
                        point_id(wing, chord + 1, span + 1), point_id(wing, chord + 1, span)};
                    const std::array contour{
                        edge(corners[0], corners[1]), edge(corners[1], corners[2]),
                        edge(corners[2], corners[3]), edge(corners[3], corners[0])};
                    const std::array facets{a::FacetId{500 + 2 * index},
                                            a::FacetId{501 + 2 * index}};
                    faces.push_back({face, band});
                    quads.emplace_back(face, corners, contour, edge(corners[0], corners[2]),
                                       facets);
                    const bool wake = chord == chord_count;
                    roles.push_back({face, wake ? a::AeroFaceRole::Wake : a::AeroFaceRole::Bound});
                    if (wake) {
                        trailing.push_back({contour[0], {corners[0], corners[1]}});
                        expected_owners.push_back(faces.at(index - band_count).id);
                        wake_facets.push_back(facets[0]);
                    } else {
                        bound.push_back(face);
                        bound_facets.push_back(facets[0]);
                    }
                }
            }
        }
    }

    [[nodiscard]] a::GeometryInput input(double radius) const noexcept {
        const a::SinglePositions xyz{values};
        return {points,
                bindings,
                xyz,
                xyz,
                edges,
                faces,
                bands,
                {},
                quads,
                a::SingleWord{radius},
                test_geometry::context()};
    }

    std::vector<a::PointId> points;
    std::vector<a::PointBinding> bindings;
    std::vector<double> values;
    std::vector<a::EdgeInput> edges;
    std::vector<a::EdgeBindingInput> edge_bindings;
    std::vector<a::FaceInput> faces;
    std::vector<a::BandId> bands;
    std::vector<a::QuadInput> quads;
    std::vector<a::FaceRoleInput> roles;
    std::vector<a::TrailingEdgeInput> trailing;
    std::vector<a::PhysicalFaceId> bound, expected_owners;
    std::vector<a::FacetId> bound_facets, wake_facets;
};

struct Product {
    explicit Product(const Mesh& mesh, double radius = core_radius)
        : geometry_work(test_geometry::workspace(memory)),
          geometry(test_geometry::take(a::prepare_owned_geometry(
              mesh.input(radius), geometry_limits, *geometry_work, memory))),
          global(take(a::GlobalWorkspace::create(memory, global_limits.metadata))),
          topology(take_global(
              a::build_global_topology(geometry, {mesh.roles, mesh.edge_bindings, mesh.trailing},
                                       global_limits, *global, memory))),
          contact(take_global(
              a::check_global_geometry(topology, geometry, global_limits, *global, memory))),
          local(take(a::LocalPacketWorkspace::create(memory))),
          kernel(take(a::P1Workspace::create(memory))) {
        REQUIRE(topology.trailing_edges().size() == band_count * wing_count);
        for (const auto face : mesh.bound) {
            maps.emplace_back(face, a::P0BoundColumn{});
            strengths.push_back(a::P1StrengthId::bound(face));
        }
        for (std::size_t i = 0; i < mesh.wake_facets.size(); ++i) {
            CAPTURE(i);
            const auto row =
                take_global(contact.uvlm_quad(mesh.wake_facets[i], geometry_limits.max_steps));
            if (const auto& attachment = row.attachment) {
                CHECK(attachment->edge == mesh.trailing.at(i).edge);
                CHECK(attachment->bound_face == mesh.expected_owners.at(i));
                REQUIRE(attachment->coefficient == 1);
                const auto key = a::P1StrengthId::bound(attachment->bound_face);
                maps.emplace_back(row.face, std::array{key, key});
            } else {
                FAIL("The wake front has no qualified trailing-edge owner");
            }
        }
        std::ranges::sort(maps, {}, &a::P1FaceStrengths::face);
    }

    [[nodiscard]] a::ResolvedQuery query(std::size_t index, const Point& where) {
        return test_geometry::take(
            a::resolve_query(geometry, {a::QueryId{1 + index}, a::FreePoint{a::FreeSingle{where}}},
                             geometry.context(), geometry_limits, *geometry_work, memory));
    }

    [[nodiscard]] a::SharedP1Action prepare(const a::P1ActionSource& source) {
        return test_geometry::take(a::prepare_bound_wake_induction(
            source, maps, strengths, {core::Precision::TwoWord, 1e-12},
            {core::Precision::TwoWord, 1e-12, 1e-12}, action_limits, shared, *local, *kernel,
            memory));
    }

    const a::GeometryLimits geometry_limits{
        {128, 512, 128, 32, 0, 128}, 20000000, {ls::EstimateAdmission{65536}, 1000000}};
    const a::GlobalLimits global_limits{std::numeric_limits<std::uint64_t>::max() / 4,
                                        ls::EstimateAdmission{65536}};
    const a::P1AssemblyLimits action_limits{geometry_limits,
                                            {ls::EstimateAdmission{65536}, UINT64_C(1000000000)},
                                            UINT64_C(150000000000)};
    core::ResourcePlan memory{UINT64_C(2000000000)};
    std::unique_ptr<a::GeometryWorkspace> geometry_work;
    a::OwnedGeometry geometry;
    std::unique_ptr<a::GlobalWorkspace> global;
    a::GlobalTopology topology;
    a::GlobalContactState contact;
    a::SharedP1Workspace shared;
    std::unique_ptr<a::LocalPacketWorkspace> local;
    std::unique_ptr<a::P1Workspace> kernel;
    std::vector<a::P1FaceStrengths> maps;
    std::vector<a::P1StrengthId> strengths;
};

using Coefficients = std::array<a::P1ComplexWords, order>;
using WorkColumns = std::array<Enclosure, order>;

struct RetainedRow {
    a::ResolvedQuery query;
    a::SharedP1Action action;
};

Enclosure add_scaled(const Enclosure& initial, const Enclosure& increment, const Rat& factor) {
    return {initial.center + increment.center * factor,
            initial.radius + increment.radius * abs(factor)};
}

RetainedRow compare_row(Product& product, const nlohmann::json& fixture, std::size_t index,
                        const Coefficients& direction, Matrix& matrix, Enclosure& work_forward,
                        WorkColumns& work_reverse) {
    CAPTURE(index);
    auto query = product.query(index, point(fixture.at("queries_m").at(index)));
    const a::P1ActionSource source{product.contact, query, {700, 9}};
    auto action = product.prepare(source);
    const Coefficients zero{};
    const auto transpose = test_geometry::take(
        a::adjoint_shared_p1(action, source, zero, {scalar(0), scalar(0), scalar(1)}, action_steps,
                             product.shared, product.memory));
    const auto forward = test_geometry::take(
        a::apply_shared_p1(action, source, direction, {}, action_steps, product.shared));
    REQUIRE(action.strength_count() == order);
    const auto& expected = fixture.at("matrix_per_m").at(index);
    REQUIRE(expected.size() == order);
    const Rat receiver_weight = Rat{static_cast<int>(index % 7) - 3} / 8;
    Rat independent_velocity{};
    for (std::size_t j = 0; j < order; ++j) {
        CAPTURE(j);
        CHECK(test_geometry::take(action.strength(j)) == product.strengths.at(j));
        const auto coefficient = enclosure(test_geometry::take(transpose.strength(j)).real);
        const Rat independent = decimal(expected.at(j));
        check(coefficient, independent, Rat{2e-10});
        matrix[j * order + index] = coefficient.center.convert_to<double>();
        independent_velocity += independent * word(direction[j].real);
        work_reverse[j] = add_scaled(work_reverse[j], coefficient, receiver_weight);
    }
    const auto velocity = enclosure(forward.velocity()[2].real);
    check(velocity, independent_velocity, Rat{2e-10});
    work_forward = add_scaled(work_forward, velocity, receiver_weight);
    return {std::move(query), std::move(action)};
}

class AlgebraicOnly final : public ls::ConsumerEvaluator<double> {
public:
    [[nodiscard]] numerics::SupportedPrecision supported_precision() const noexcept override {
        return numerics::SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const numerics::PrecisionBuffer<double>&, const numerics::PrecisionBuffer<double>&,
           const ls::Context&, const ls::ConsumerCheck& check,
           core::CertificateContext expected) const override {
        FAIL_CHECK("Physical residual is evaluated by the retained mixed velocity action");
        return core::ErrorCertificate::unevaluated(expected, check.budget);
    }
};

ls::Context trial_context(const a::GeometryContext& geometry) {
    const auto trial = geometry.diagnostics.trial;
    if (!trial)
        FAIL("Expected the original circulation trial");
    return {{geometry.fingerprint, geometry.model, geometry.current, core::ChartGeneration{1},
             core::SparsityGeneration{2}, core::NumericalGeneration{3}, core::ScaleGeneration{4},
             geometry.precision, core::Precision::Binary64, core::WakePolicyGeneration{5},
             core::TimeIntervalGeneration{6}, core::BackendGeneration{7}},
            *trial,
            core::ScopeId{8},
            core::ScopeId{9},
            geometry.time};
}

struct SolveSetup {
    SolveSetup(core::Index size, const a::GeometryContext& geometry)
        : context(trial_context(geometry)),
          limits{size, 1, 0, 5, ls::EstimateAdmission{1U << 20U}},
          scale(static_cast<std::size_t>(size), 20.),
          gamma_units(scale.size(), {2, 0, -1, 0, 0, 0, 0}),
          speed_units(scale.size(), {1, 0, -1, 0, 0, 0, 0}) {}

    [[nodiscard]] ls::Scales scales() const {
        return {scale, scale, gamma_units, speed_units, context};
    }

    ls::Context context;
    ls::Limits limits;
    ls::Accuracy accuracy{ls::Norm::Infinity,     1e-25, 0., 1e-22, core::ObservationId{20},
                          core::ObservationId{21}};
    std::vector<double> scale;
    std::vector<ls::UnitDimensions> gamma_units, speed_units;
};

std::vector<a::ScalarWords> solve(std::span<const double> coefficients, core::Index size,
                                  const a::GeometryContext& geometry, core::ResourcePlan& memory) {
    constexpr auto operation = core::Operation::SolveLinear;
    const SolveSetup setup{size, geometry};
    const auto scales = setup.scales();
    const auto count = static_cast<std::size_t>(size);
    const auto matrix = take(
        numerics::DenseMatrix<double>::copy({size, size, size}, coefficients, memory, operation));
    // Both planes have upward normal +z. The physical ambient direction is (0,0,20).
    const std::vector<double> rhs(count, -20.);
    auto source =
        take(numerics::DenseMatrix<double>::copy({size, 1, size}, rhs, memory, operation));
    const auto right = take(numerics::PrecisionBuffer<double>::binary64(
        std::move(source), core::NumericalGeneration{30}));
    const ls::DenseOperator op{matrix, setup.context};
    auto prepared = take(ls::DensePrepared::prepare(op, scales, setup.limits, memory));
    const AlgebraicOnly algebraic;
    const ls::Precision<double> precision{
        ls::RetainedPrecision{core::NumericalGeneration{31}}, {}, algebraic};
    const auto solved =
        take(ls::solve_linear(op, ls::Rhs{right, setup.context, core::RhsGeneration{40}}, scales,
                              precision, setup.accuracy, setup.limits, prepared, memory));

    std::vector<a::ScalarWords> result;
    result.reserve(count);
    std::visit(
        [&](const auto& stored) {
            using T = std::decay_t<decltype(stored)>;
            for (std::size_t i = 0; i < count; ++i) {
                if constexpr (std::is_same_v<T, numerics::Binary64<double>>)
                    result.emplace_back(a::SingleWord{stored.values.values()[i]});
                else
                    result.emplace_back(
                        a::DoubleWord{stored.hi.values()[i], stored.lo.values()[i]});
            }
        },
        solved.values().storage());
    return result;
}

std::array<Rat, 2> first_order_wrench(std::span<const a::ScalarWords> gamma,
                                      const nlohmann::json& sections) {
    REQUIRE(gamma.size() == wing_order);
    std::array<Rat, 2> result{};
    // At the zero-incidence equilibrium the old sheet and bound circulation are zero.
    // The spanwise force derivative is rho U_x times the strength jump and edge span.
    for (std::size_t c = 0; c < chord_count; ++c) {
        for (std::size_t s = 0; s < band_count; ++s) {
            Rat strength = word(gamma[c * band_count + s]);
            if (c != 0)
                strength -= word(gamma[(c - 1) * band_count + s]);
            const auto left = point(sections.at(c).at(s));
            const auto right = point(sections.at(c).at(s + 1));
            const Rat force = 24 * strength * (Rat{right[1]} - Rat{left[1]});
            result[0] += force;
            result[1] -= (Rat{left[0]} + Rat{right[0]}) * force / 2;
        }
    }
    return result;
}

void check_wrenches(const nlohmann::json& fixture, const Matrix& matrix,
                    std::span<const a::ScalarWords> gamma, Product& product) {
    for (std::size_t wing = 0; wing < wing_count; ++wing) {
        CAPTURE(wing);
        const auto first = wing * wing_order;
        std::array<double, wing_order * wing_order> isolated{};
        for (std::size_t row = 0; row < wing_order; ++row)
            for (std::size_t col = 0; col < wing_order; ++col)
                isolated[col * wing_order + row] = matrix[(first + col) * order + first + row];
        const auto alone = solve(isolated, static_cast<core::Index>(wing_order),
                                 product.geometry.context(), product.memory);
        const auto& sections = fixture.at("vertices_m").at(wing);
        const auto coupled = first_order_wrench(gamma.subspan(first, wing_order), sections);
        const auto single = first_order_wrench(alone, sections);
        for (std::size_t component = 0; component < coupled.size(); ++component) {
            CAPTURE(component);
            const auto& expected = fixture.at("wrenches").at(wing);
            const Rat expected_coupled = decimal(expected.at(0).at(component));
            const Rat expected_single = decimal(expected.at(1).at(component));
            CHECK(abs(coupled[component] - expected_coupled) <= Rat{2e-8});
            CHECK(abs(single[component] - expected_single) <= Rat{2e-8});
            CHECK(abs((coupled[component] - single[component]) -
                      (expected_coupled - expected_single)) <= Rat{4e-8});
            CHECK(abs(coupled[component] - single[component]) > Rat{1});
        }
    }
}

void check_mismatch(const Product& other, const a::ResolvedQuery& other_query,
                    const RetainedRow& row, const Coefficients& values, Product& original) {
    const a::P1ActionSource wrong{other.contact, other_query, {700, 9}};
    const auto rejected =
        a::apply_shared_p1(row.action, wrong, values, {}, action_steps, original.shared);
    if (rejected)
        FAIL("A prepared action accepted a different qualified geometry or core");
    CHECK(rejected.error().cause.code == core::FailureCode::StaleGeneration);
}

namespace cf = coupling::detail;
namespace p = a::pchip;

cf::ProjectedCaseUniformFlow project_ambient(core::ResourcePlan& memory, p::PchipInputStamp stamp) {
    using namespace test::schema_adapter;
    constexpr auto admission = std::uint64_t{1} << 48U;
    const structure::ExactLimits exact{100000,
                                       100000,
                                       100000,
                                       64,
                                       100000,
                                       200000,
                                       100000,
                                       structure::ExactEstimateAdmission{admission},
                                       structure::ExactTimingMode::Disabled};
    const cf::CaseUniformFlowLimits flow_limits{
        exact, {{128, 100000, admission, 10000}, ls::EstimateAdmission{admission}}, 100};
    auto recipe = root("aerodynamic_steady");
    recipe.field("flow").field("uniform").set("value", parse("[20,0,0]"));
    recipe.field("flow").set("density", parse("1.2"));
    recipe.field("analysis").field("speed").set("value", parse("20"));
    const auto wire = recipe.wire();
    const io::InputLimits syntax{
        1000000, 256,     100000,  100000,
        1000000, 1000000, 1000000, {io::ParserAllocationMode::EstimateAdmission, admission}};
    std::optional<io::InputDiagnostic> parsing;
    const auto doc =
        take(io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                     io::SourceId{620}, syntax, memory, parsing));
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> schema_caps{};
    schema_caps.fill(admission);
    const auto checked =
        take(io::validate_case_schema(doc, io::SchemaLimits{schema_caps}, memory, nullptr));
    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> si_caps{};
    si_caps.fill(admission);
    const io::SiLimits si_limits{si_caps, exact};
    const auto catalog = take(io::make_si_catalog(checked, si_limits, memory, nullptr));
    std::vector<io::SiRequirementRow> rows;
    for (std::uint64_t i = 0; i < take(catalog.size()); ++i) {
        const auto quantity = take(catalog.entry(i));
        rows.push_back({quantity.key,
                        {io::SiDecimalInput{"1"}, quantity.recipe.dimensions},
                        io::SiDecimalInput{"1e-25"},
                        io::SiDecimalInput{"1e-25"}});
    }
    const auto request = take(io::make_si_request(catalog, rows, si_limits, memory, nullptr));
    const auto si =
        take(io::convert_case_to_si(checked, catalog, request, si_limits, memory, nullptr));
    a::UniformFlowWorkspace workspace;
    const auto field =
        take(cf::prepare_case_uniform_flow(si, stamp, {}, workspace, flow_limits, memory));
    constexpr std::array<std::uint32_t, 1> five{5}, four{4}, one{1}, billion{1000000000};
    const p::RationalView time{1, five, four}, scale{1, one, one}, absolute{1, one, billion},
        zero{0, {}, one};
    auto sampled = take(cf::evaluate_case_uniform_flow(field, field, {stamp, time}, workspace,
                                                       flow_limits.flow, memory));
    const cf::CaseFlowComponentBudget component{scale, absolute, zero};
    const std::array bounds{component, component, component};
    return take(cf::project_case_uniform_flow(
        std::move(sampled), {core::Precision::TwoWord, bounds, bounds, {4096, admission}}));
}

cf::NonpenetrationVector vector_value(const Point& pnt) {
    return {a::SingleWord{pnt[0]}, a::SingleWord{pnt[1]}, a::SingleWord{pnt[2]}};
}

Enclosure multiply(const Enclosure& left, const Enclosure& right) {
    return {left.center * right.center, abs(left.center) * right.radius +
                                            abs(right.center) * left.radius +
                                            left.radius * right.radius};
}

std::size_t wing_of(a::PointId id) {
    REQUIRE(id.value() > 0);
    const auto wing = (id.value() - 1) / points_per_wing;
    REQUIRE(wing < wing_count);
    return static_cast<std::size_t>(wing);
}

void compare_native_loads(Product& product, const Mesh& mesh, const nlohmann::json& fixture) {
    constexpr a::P1ActionStamp trial{700, 9};
    const p::PchipInputStamp stamp{p::PchipIssuerId{72}, p::PchipInputGeneration{1}};
    constexpr core::NumericalGeneration generation{31};
    const auto zero = vector_value({});
    const cf::ExplicitCollocationTangents planar{vector_value({1, 0, 0}), vector_value({0, 1, 0})};
    std::vector<a::ResolvedQuery> queries;
    queries.reserve(order);
    for (std::size_t i = 0; i < order; ++i)
        queries.push_back(product.query(i, point(fixture.at("queries_m").at(i))));
    std::vector<cf::DiscreteCollocation> collocations;
    std::vector<cf::BoundCirculationRate> rates;
    std::vector<cf::BoundLoadDirection> gamma_direction;
    for (std::size_t i = 0; i < order; ++i) {
        collocations.emplace_back(mesh.bound_facets.at(i), queries.at(i), planar, zero,
                                  a::SingleWord{0}, a::SingleWord{0});
        rates.emplace_back(mesh.bound.at(i), a::SingleWord{0});
        const auto derivative = decimal(fixture.at("circulation_slope_m2_s").at(i));
        gamma_direction.emplace_back(mesh.bound.at(i), a::SingleWord{derivative}, a::SingleWord{0});
    }
    std::vector<cf::DiscreteWakePanel> wake;
    std::vector<cf::FrozenWakeTrace> history;
    std::vector<a::EdgeId> te_edges, force_edges;
    for (std::size_t i = 0; i < mesh.wake_facets.size(); ++i) {
        const a::P1StrengthId id{800 + i};
        wake.emplace_back(mesh.wake_facets.at(i), cf::AttachedWakeFront{}, id);
        history.push_back({id, a::SingleWord{0}});
        te_edges.push_back(mesh.trailing.at(i).edge);
    }
    for (const auto& quad : mesh.quads) {
        if (std::ranges::binary_search(mesh.bound, quad.face))
            force_edges.insert(force_edges.end(), quad.edges.begin(), quad.edges.end());
    }
    std::ranges::sort(force_edges);
    auto duplicates = std::ranges::unique(force_edges);
    force_edges.erase(duplicates.begin(), duplicates.end());
    std::erase_if(force_edges,
                  [&](a::EdgeId id) { return std::ranges::find(te_edges, id) != te_edges.end(); });
    std::vector<cf::EdgeMidpointVelocity> edges;
    edges.reserve(force_edges.size());
    for (const auto edge : force_edges)
        edges.emplace_back(edge, zero);
    const std::vector<cf::NonpenetrationVector> displacements(mesh.points.size(), zero);
    const cf::DiscreteNonpenetrationInput input{product.contact,
                                                collocations,
                                                wake,
                                                {core::AcceptedStateId{2}, history},
                                                trial,
                                                stamp,
                                                p::RequestedSide::Automatic};
    const cf::DiscreteVectorLoadInput loads{trial, generation,         edges,
                                            rates, a::SingleWord{1.2}, zero};
    const cf::DiscreteVectorLoadDirection variation{trial,
                                                    generation,
                                                    mesh.points,
                                                    displacements,
                                                    edges,
                                                    gamma_direction,
                                                    vector_value({0, 0, 20})};
    const cf::DiscreteNonpenetrationLimits limits{{core::Precision::TwoWord, 1e-12},
                                                  {core::Precision::TwoWord, 1e-12, 1e-12},
                                                  product.action_limits,
                                                  action_steps,
                                                  UINT64_C(70000000000000)};
    const SolveSetup setup{static_cast<core::Index>(order), product.geometry.context()};
    const cf::DiscreteCirculationRequest request{setup.scales(),
                                                 setup.limits,
                                                 setup.accuracy,
                                                 core::RhsGeneration{3},
                                                 generation,
                                                 1e-13,
                                                 1e-13,
                                                 1e-11,
                                                 0.};
    const auto result = take(cf::solve_discrete_circulation_load_direction(
        input, loads, variation, project_ambient(product.memory, stamp), limits, request,
        *product.geometry_work, product.shared, *product.local, *product.kernel, product.memory));
    using Wrench = std::array<Enclosure, 2>;
    std::array<Wrench, wing_count> per_wing{};
    REQUIRE(result.loads().edge_count() == edges.size());
    for (std::size_t i = 0; i < result.loads().edge_count(); ++i) {
        CAPTURE(i);
        const auto location = take(result.loads().edge(i));
        const auto derivative = take(result.edge(i));
        const auto wing = wing_of(location.endpoints[0]);
        REQUIRE(wing_of(location.endpoints[1]) == wing);
        const auto force_z = enclosure(derivative.force[2]);
        const auto moment_y = multiply(enclosure(location.position[0]), force_z);
        check(enclosure(derivative.force[0]), Rat{0}, Rat{2e-8});
        check(enclosure(derivative.force[1]), Rat{0}, Rat{2e-8});
        per_wing[wing][0] = add_scaled(per_wing[wing][0], force_z, Rat{1});
        per_wing[wing][1] = add_scaled(per_wing[wing][1], moment_y, Rat{-1});
    }
    Wrench total{};
    for (std::size_t wing = 0; wing < wing_count; ++wing) {
        CAPTURE(wing);
        const auto& expected = fixture.at("wrenches").at(wing).at(0);
        for (std::size_t component = 0; component < total.size(); ++component) {
            check(per_wing[wing][component], Rat{decimal(expected.at(component))}, Rat{2e-8});
            total[component] = add_scaled(total[component], per_wing[wing][component], Rat{1});
        }
    }
    const auto force = result.total_force_direction();
    const auto moment = result.total_moment_direction();
    check(add_scaled(enclosure(force[2]), total[0], Rat{-1}), Rat{0}, Rat{2e-8});
    check(add_scaled(enclosure(moment[1]), total[1], Rat{-1}), Rat{0}, Rat{2e-8});

    Enclosure nodal_work{};
    const std::array<Rat, wing_count> translation{Rat{1} / 4, -Rat{1} / 8};
    const std::array<Rat, wing_count> rotation{-Rat{1} / 16, Rat{1} / 32};
    for (std::size_t i = 0; i < result.loads().point_count(); ++i) {
        const auto node = take(result.point(i));
        const auto wing = wing_of(node.point);
        const auto offset = 3 * static_cast<std::size_t>(node.point.value() - 1);
        const Rat x = mesh.values.at(offset);
        const Rat z = mesh.values.at(offset + 2);
        nodal_work = add_scaled(nodal_work, enclosure(node.force[0]), rotation[wing] * z);
        nodal_work = add_scaled(nodal_work, enclosure(node.force[2]),
                                translation[wing] - rotation[wing] * x);
    }
    Enclosure wrench_work{};
    for (std::size_t wing = 0; wing < wing_count; ++wing) {
        wrench_work = add_scaled(wrench_work, per_wing[wing][0], translation[wing]);
        wrench_work = add_scaled(wrench_work, per_wing[wing][1], rotation[wing]);
    }
    check(add_scaled(nodal_work, wrench_work, Rat{-1}), Rat{0}, Rat{2e-8});
}
} // namespace

TEST_CASE("Frozen main and tail wings share the same physical wake action",
          "[V-MULTI][UW-D01][UW-D02][UW-D05][verification]") {
    std::ifstream input{UVLAS_MULTIPLE_SURFACE_FIXTURE};
    REQUIRE(input);
    const auto fixture = nlohmann::json::parse(input);
    REQUIRE(fixture.at("m") == chord_count);
    REQUIRE(fixture.at("n") == band_count);
    REQUIRE(decimal(fixture.at("radius_m")) == core_radius);
    const Mesh mesh{fixture.at("vertices_m")};
    auto product = std::make_unique<Product>(mesh);

    if (UVLAS_MULTIPLE_SURFACE_NATIVE_LOADS != 0) {
        compare_native_loads(*product, mesh, fixture);
        return;
    }
    Matrix matrix{};
    Coefficients direction{};
    for (std::size_t j = 0; j < order; ++j)
        direction[j] = scalar((static_cast<double>(j % 5) - 2.) / 8. + (j < wing_order ? .5 : -.5));
    Enclosure work_forward{};
    WorkColumns work_reverse{};
    std::vector<RetainedRow> rows;
    rows.reserve(order);
    for (std::size_t i = 0; i < order; ++i)
        rows.push_back(
            compare_row(*product, fixture, i, direction, matrix, work_forward, work_reverse));
    Enclosure paired{};
    for (std::size_t j = 0; j < order; ++j)
        paired = add_scaled(paired, work_reverse[j], word(direction[j].real));
    CHECK(abs(paired.center - work_forward.center) + paired.radius + work_forward.radius <=
          Rat{2e-10});

    const auto gamma = solve(matrix, static_cast<core::Index>(order), product->geometry.context(),
                             product->memory);
    Coefficients solved{};
    for (std::size_t i = 0; i < order; ++i) {
        CAPTURE(i);
        CHECK(abs(word(gamma.at(i)) - Rat{decimal(fixture.at("circulation_slope_m2_s").at(i))}) <=
              Rat{2e-10});
        solved[i] = {gamma.at(i), a::SingleWord{0}};
    }
    for (std::size_t i = 0; i < order; ++i) {
        CAPTURE(i);
        const auto& row = rows.at(i);
        const a::P1ActionSource source{product->contact, row.query, {700, 9}};
        const auto forward = test_geometry::take(
            a::apply_shared_p1(row.action, source, solved, {}, action_steps, product->shared));
        check(enclosure(forward.velocity()[2].real), Rat{-20}, Rat{2e-10});
    }
    check_wrenches(fixture, matrix, gamma, *product);

    const auto where = point(fixture.at("queries_m").at(0));
    {
        auto changed_core = std::make_unique<Product>(mesh, 2 * core_radius);
        const auto alternate = changed_core->query(0, where);
        check_mismatch(*changed_core, alternate, rows.front(), solved, *product);
    }
    {
        auto shifted = mesh;
        for (std::size_t i = 3 * points_per_wing; i < shifted.values.size(); i += 3)
            shifted.values[i] += .125;
        auto changed_geometry = std::make_unique<Product>(shifted);
        const auto alternate = changed_geometry->query(0, where);
        check_mismatch(*changed_geometry, alternate, rows.front(), solved, *product);
    }
}
