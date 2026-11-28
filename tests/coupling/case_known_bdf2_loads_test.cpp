#include "aero/geometry_snapshot_fixture.hpp"
#include "coupling/case_flow_clock.hpp"
#include "coupling/discrete_nonpenetration.hpp"
#include "coupling/fixed_bdf2_circulation.hpp"
#include "io/case_schema_test_adapter.hpp"

#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <limits>
#include <source_location>
#include <vector>

namespace {
using namespace uvlas;
using namespace uvlas::test::schema_adapter;
namespace a = aero::detail;
namespace p = a::pchip;
namespace cf = coupling::detail;
namespace ls = numerics::linear;

constexpr std::uint64_t budget = std::uint64_t{1} << 48U;
const p::PchipInputStamp stamp{p::PchipIssuerId{71}, p::PchipInputGeneration{8}};

structure::ExactLimits exact_limits() noexcept {
    return {100000,
            100000,
            100000,
            64,
            64,
            100000,
            200000,
            structure::ExactEstimateAdmission{budget},
            structure::ExactTimingMode::Disabled};
}
const cf::CaseUniformFlowLimits limits{
    exact_limits(),
    {{128, 100000, budget, 10000}, numerics::linear::EstimateAdmission{budget}},
    100};

template <class T>
T take(core::Result<T> value, const std::source_location at = std::source_location::current()) {
    if (!value) {
        INFO(at.file_name() << ':' << at.line());
        const auto& resources = value.error().details.resources;
        if (resources) {
            const auto& usage = *resources;
            INFO("Resource kind " << static_cast<int>(usage.kind) << ", required " << usage.required
                                  << ", observed " << usage.observed);
            FAIL("Unexpected case-flow failure " << static_cast<int>(value.error().code));
        }
        FAIL("Unexpected case-flow failure " << static_cast<int>(value.error().code));
    }
    return std::move(*value);
}

template <class T, class Diagnostic>
T take_io(core::Result<T> value, const Diagnostic& diag) {
    if (!value) {
        const auto path = diag.path();
        if (path)
            FAIL("Unexpected case-flow input at " << *path << " code "
                                                  << static_cast<int>(value.error().code));
        FAIL("Unexpected case-flow preflight " << static_cast<int>(value.error().code));
    }
    return std::move(*value);
}

template <class T>
void reject(const core::Result<T>& value, core::FailureCode code,
            const std::source_location at = std::source_location::current()) {
    INFO(at.line());
    if (value)
        FAIL("Unexpected case-flow success");
    CHECK(value.error().code == code);
}

io::SiCase materialize(const JsonToken& source, io::SourceId id) {
    core::ResourcePlan memory{budget};
    const auto wire = source.wire();
    const io::InputLimits syntax{
        1000000, 256,     100000,  100000,
        1000000, 1000000, 1000000, {io::ParserAllocationMode::EstimateAdmission, budget}};
    std::optional<io::InputDiagnostic> parsing;
    const auto doc = take(io::read_input_document(
        std::as_bytes(std::span{wire.data(), wire.size()}), id, syntax, memory, parsing));
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> schema_caps{};
    schema_caps.fill(budget);
    io::SchemaDiagnostic schema_diagnostic;
    const auto checked = take_io(
        io::validate_case_schema(doc, io::SchemaLimits{schema_caps}, memory, &schema_diagnostic),
        schema_diagnostic);

    std::array<std::uint64_t, static_cast<std::size_t>(io::SiQuota::Count)> si_caps{};
    si_caps.fill(budget);
    const io::SiLimits si_limits{si_caps, exact_limits()};
    io::SiDiagnostic si_diagnostic;
    const auto catalog =
        take_io(io::make_si_catalog(checked, si_limits, memory, &si_diagnostic), si_diagnostic);
    std::vector<io::SiRequirementRow> requirements;
    for (std::uint64_t i = 0; i < take(catalog.size()); ++i) {
        const auto quantity = take(catalog.entry(i));
        requirements.push_back({quantity.key,
                                {io::SiDecimalInput{"1"}, quantity.recipe.dimensions},
                                io::SiDecimalInput{"1e-25"},
                                io::SiDecimalInput{"1e-25"}});
    }
    const auto request =
        take_io(io::make_si_request(catalog, requirements, si_limits, memory, &si_diagnostic),
                si_diagnostic);
    return take_io(
        io::convert_case_to_si(checked, catalog, request, si_limits, memory, &si_diagnostic),
        si_diagnostic);
}

struct Number {
    Number(std::uint32_t numerator, std::uint32_t denominator = 1, int sign = 1) noexcept
        : n{numerator}, d{denominator}, sign_(numerator == 0 ? 0 : sign) {}
    [[nodiscard]] p::RationalView view() const& noexcept {
        return {sign_, sign_ == 0 ? std::span<const std::uint32_t>{} : n, d};
    }
    [[nodiscard]] p::RationalView view() const&& = delete;
    std::array<std::uint32_t, 1> n, d;
    int sign_;
};
JsonToken model_case() {
    auto input = root("coupled_transient");
    input.field("units").set("length", string("m"));
    input.field("units").set("time", string("s"));
    input.field("units").set("angle", string("deg"));
    input.field("flow").set("uniform", parse(R"({"kind":"constant","value":[0,0,1]})"));
    input.field("flow").set(
        "gust",
        parse(
            R"({"kind":"sine","amplitude":{"kind":"constant","value":[1]},"direction":[0,0,1],"wavevector":[180,0,0],"angular_frequency":180,"phase":0,"origin":[0,0,0],"time_origin":0})"));
    input.field("analysis")
        .set("intervals", parse(R"([{"start":0,"end":0.75,"step":0.25,"restart":"initial"}])"));
    input.field("analysis")
        .field("wake")
        .set("convection",
             parse(R"({"kind":"prescribed","velocity":{"kind":"constant","value":[1,0,0]}})"));
    return input;
}

const cf::CaseKnownFlowLimits gust_limits{limits, 2048, 10000};

cf::CaseKnownFlowProjection gust_projection(core::Precision precision, const Number& absolute) {
    const static Number unit{1}, zero{0};
    const cf::CaseFlowComponentBudget per_component{unit.view(), absolute.view(), zero.view()};
    const std::array<cf::CaseFlowComponentBudget, 3> vector{per_component, per_component,
                                                            per_component};
    const std::array<cf::CaseFlowComponentBudget, 9> matrix{
        per_component, per_component, per_component, per_component, per_component,
        per_component, per_component, per_component, per_component};
    return {precision, vector, vector, matrix};
}

cf::NonpenetrationVector coordinates(double x, double y, double z) noexcept {
    return {a::SingleWord{x}, a::SingleWord{y}, a::SingleWord{z}};
}

template <class T>
T take_global(a::GlobalResult<T> result,
              const std::source_location at = std::source_location::current()) {
    if (!result) {
        INFO(at.line());
        FAIL("Unexpected global qualification failure "
             << static_cast<int>(result.error().origin.cause.code));
    }
    return std::move(*result);
}
ls::Context circulation_context(const a::GeometryContext& geometry) {
    const auto trial = geometry.diagnostics.trial;
    if (!trial)
        FAIL("The circulation fixture requires a trial");
    return {{geometry.fingerprint, geometry.model, geometry.current, core::ChartGeneration{1},
             core::SparsityGeneration{2}, core::NumericalGeneration{3}, core::ScaleGeneration{4},
             geometry.precision, core::Precision::Binary64, core::WakePolicyGeneration{5},
             core::TimeIntervalGeneration{6}, core::BackendGeneration{7}},
            *trial,
            core::ScopeId{8},
            core::ScopeId{9},
            geometry.time};
}

cf::DiscreteCirculationRequest circulation_request(const ls::Scales& scales, core::Index order,
                                                   ls::Norm norm) {
    return {scales,
            {order, 1, 0, 5, ls::EstimateAdmission{1U << 20U}},
            {norm, 1e-28, 0., 1e-26, core::ObservationId{20}, core::ObservationId{21}},
            core::RhsGeneration{3},
            core::NumericalGeneration{4},
            1e-14,
            1e-14,
            1e-13,
            0.};
}

struct DisjointCirculationPanels {
    std::array<a::PointId, 8> points{a::PointId{1}, a::PointId{2}, a::PointId{3}, a::PointId{4},
                                     a::PointId{5}, a::PointId{6}, a::PointId{7}, a::PointId{8}};
    std::array<a::PointBinding, 8> bindings{test_geometry::binding(), test_geometry::binding(),
                                            test_geometry::binding(), test_geometry::binding(),
                                            test_geometry::binding(), test_geometry::binding(),
                                            test_geometry::binding(), test_geometry::binding()};
    std::array<double, 24> xyz{0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0,
                               0, 2, 0, 0, 3, 0, 1, 3, 0, 1, 2, 0};
    std::array<a::BandId, 2> bands{a::BandId{50}, a::BandId{51}};
    std::array<a::FaceInput, 2> faces{a::FaceInput{a::PhysicalFaceId{100}, bands[0]},
                                      a::FaceInput{a::PhysicalFaceId{101}, bands[1]}};
    std::array<a::FaceRoleInput, 2> roles{a::FaceRoleInput{faces[0].id, a::AeroFaceRole::Bound},
                                          a::FaceRoleInput{faces[1].id, a::AeroFaceRole::Bound}};
    std::array<a::EdgeInput, 10> edges{a::EdgeInput{a::EdgeId{10}, {points[0], points[1]}},
                                       a::EdgeInput{a::EdgeId{11}, {points[1], points[2]}},
                                       a::EdgeInput{a::EdgeId{12}, {points[2], points[3]}},
                                       a::EdgeInput{a::EdgeId{13}, {points[3], points[0]}},
                                       a::EdgeInput{a::EdgeId{14}, {points[0], points[2]}},
                                       a::EdgeInput{a::EdgeId{15}, {points[4], points[5]}},
                                       a::EdgeInput{a::EdgeId{16}, {points[5], points[6]}},
                                       a::EdgeInput{a::EdgeId{17}, {points[6], points[7]}},
                                       a::EdgeInput{a::EdgeId{18}, {points[7], points[4]}},
                                       a::EdgeInput{a::EdgeId{19}, {points[4], points[6]}}};
    std::array<a::EdgeBindingInput, 10> owners{
        a::EdgeBindingInput{edges[0].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[1].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[2].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[3].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[4].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[5].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[6].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[7].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[8].id, a::CanonicalEdge{}},
        a::EdgeBindingInput{edges[9].id, a::CanonicalEdge{}}};
    std::array<a::QuadInput, 2> quads{
        a::QuadInput{faces[0].id,
                     {points[0], points[1], points[2], points[3]},
                     {edges[0].id, edges[1].id, edges[2].id, edges[3].id},
                     edges[4].id,
                     {a::FacetId{200}, a::FacetId{201}}},
        a::QuadInput{faces[1].id,
                     {points[4], points[5], points[6], points[7]},
                     {edges[5].id, edges[6].id, edges[7].id, edges[8].id},
                     edges[9].id,
                     {a::FacetId{202}, a::FacetId{203}}}};

    [[nodiscard]] a::GeometryInput input() const noexcept {
        const a::SinglePositions positions{xyz};
        return {points,
                bindings,
                positions,
                positions,
                edges,
                faces,
                bands,
                {},
                quads,
                a::SingleWord{.125},
                test_geometry::context()};
    }
};

using High = boost::multiprecision::cpp_dec_float_50;
using Pair = std::array<High, 2>;
using PlanePoint = std::array<double, 2>;

const High pi = boost::math::constants::pi<High>();
constexpr std::array<PlanePoint, 2> collocation{{{.75, .5}, {.625, 2.375}}};
constexpr std::array<std::array<double, 2>, 3> wall{{{.125, .25}, {.375, -.125}, {.75, .5}}};
constexpr std::array wall_direction{.125, -.25};
constexpr double time_direction = .125;
constexpr double density = 1.25;

High wind(double x, const High& t) {
    return 1 + sin(pi * (High{x} - t));
}

// Integrate the regularized straight-filament kernel along the unit length segment.
High influence(PlanePoint query, PlanePoint start, PlanePoint end) {
    const High dx = end[0] - start[0], dy = end[1] - start[1];
    const High x = query[0] - start[0], y = query[1] - start[1];
    const High s = x * dx + y * dy, n = dx * y - dy * x;
    const High smooth = n * n + High{1} / 64;
    const auto integral = [&](const High& endpoint) -> High {
        return endpoint / sqrt(endpoint * endpoint + smooth);
    };
    return n * (integral(s) - integral(s - 1)) / (4 * pi * smooth);
}

High influence(PlanePoint query, std::size_t panel) {
    const double y = panel == 0 ? 0. : 2.;
    const std::array<PlanePoint, 4> vertices{{{0., y}, {0., y + 1.}, {1., y + 1.}, {1., y}}};
    High sum = 0;
    for (std::size_t i = 0; i < vertices.size(); ++i)
        sum += influence(query, vertices[i], vertices[(i + 1) % vertices.size()]);
    return sum;
}

Pair independent_gamma(const Pair& rhs) {
    const High a00 = influence(collocation[0], 0), a01 = influence(collocation[0], 1);
    const High a10 = influence(collocation[1], 0), a11 = influence(collocation[1], 1);
    const High determinant = a00 * a11 - a01 * a10;
    return {(rhs[0] * a11 - rhs[1] * a01) / determinant,
            (a00 * rhs[1] - a10 * rhs[0]) / determinant};
}

High center(const a::ScaledEnclosure& enclosure) {
    const auto sig = std::visit(
        [](const auto& words) -> High {
            if constexpr (std::is_same_v<std::decay_t<decltype(words)>, a::SingleWord>)
                return High{words.value};
            else
                return High{words.high} + High{words.low};
        },
        enclosure.center());
    return ldexp(sig, enclosure.exponent());
}

// N4-05 fixes the scale for small original residuals and force balance.
void expect(const a::ScaledEnclosure& result, const High& answer,
            const High& tolerance = High{"1e-11"}) {
    CHECK(abs(center(result) - answer) + ldexp(High{result.radius()}, result.exponent()) <=
          tolerance);
}
} // namespace

TEST_CASE("Case gust and retained BDF2 circulation drive endpoint loads and implicit actions",
          "[AE-07][UW-D06][J1-D01][FD-07][N4-06][verification][negative]") {
    constexpr auto format = core::Precision::TwoWord;
    const auto source_case = model_case();
    constexpr io::SourceId source_id{1190};
    constexpr core::AcceptedStateId parent{2};
    core::ResourcePlan memory{256ULL * 1024 * 1024};
    a::UniformFlowWorkspace flow_work;
    const auto si = materialize(source_case, source_id);
    const cf::CaseClockedKnownFlowLimits clock_limits{{limits, 1, 4}, gust_limits};
    const auto clock =
        take(cf::prepare_case_clocked_known_flows(si, stamp, flow_work, clock_limits, memory));
    const auto provider =
        take(cf::prepare_case_known_flow(si, stamp, {}, flow_work, gust_limits, memory));
    const Number one{1}, zero{0}, absolute{1, 1000000000};
    const cf::CaseFlowComponentBudget allowed{one.view(), absolute.view(), zero.view()};
    const std::array components{allowed, allowed, allowed};
    const cf::CaseFlowProjectionRequest wake{format, components, components, {4096, budget}};
    const auto ambient = gust_projection(format, absolute);
    const auto side = [](std::size_t i) {
        return i == 0 ? p::RequestedSide::Right : p::RequestedSide::Automatic;
    };
    const std::array times{Number{0}, Number{1, 4}, Number{1, 2}};
    const auto node = [&](std::size_t i) {
        const std::array position{zero.view(), zero.view(), zero.view()};
        return take(cf::evaluate_case_clocked_known_flows(clock, clock, {stamp, 0, i, side(i)},
                                                          position, ambient, wake, flow_work,
                                                          clock_limits, memory));
    };

    const DisjointCirculationPanels fixture;
    const auto base = fixture.input();
    const a::GlobalLimits cap{std::numeric_limits<std::uint64_t>::max() / 4,
                              ls::EstimateAdmission{65536}};
    auto geometry_work = test_geometry::workspace(memory);
    auto global = take(a::GlobalWorkspace::create(memory, cap.metadata));
    const auto reference = test_geometry::take(
        a::prepare_owned_geometry(base, test_geometry::limits(), *geometry_work, memory));
    const a::GlobalTopologyInput roles{fixture.roles, fixture.owners, {}};
    const auto topology =
        take_global(a::build_global_topology(reference, roles, cap, *global, memory));
    const auto contact_at = [&](std::size_t i, const a::GlobalTopology& owner) {
        auto positions = base;
        positions.context.time = .25 * static_cast<double>(i);
        positions.context.current = core::GeometryGeneration{3 + i};
        positions.exact_time = times[i].view();
        const auto geometry = test_geometry::take(
            a::prepare_owned_geometry(positions, test_geometry::limits(), *geometry_work, memory));
        return take_global(a::check_global_geometry(owner, geometry, cap, *global, memory));
    };
    const auto geometry_at = [&](const a::GlobalContactState& contact) -> const a::OwnedGeometry& {
        return take_global(contact.geometry()).get();
    };
    const auto queries_at = [&](const a::GlobalContactState& contact) {
        const auto& g = geometry_at(contact);
        const auto query = [&](std::size_t panel) {
            const auto x = collocation[panel];
            return test_geometry::take(a::resolve_query(
                g, {a::QueryId{880 + panel}, a::FreePoint{a::FreeSingle{{x[0], x[1], 0.}}}},
                g.context(), test_geometry::limits(), *geometry_work, memory));
        };
        return std::array{query(0), query(1)};
    };
    const auto rows_at = [&](std::size_t i, const std::array<a::ResolvedQuery, 2>& queries) {
        const auto row = [&](std::size_t panel) {
            const std::array frame{coordinates(1, 0, 0), coordinates(0, 1, 0)};
            return cf::DiscreteCollocation{
                fixture.quads[panel].facets[0],    queries[panel],   frame,
                coordinates(0, 0, wall[i][panel]), a::SingleWord{8}, a::SingleWord{0}};
        };
        // Exercise retained physical IDs across a change of caller row order.
        return i == 1 ? std::array{row(1), row(0)} : std::array{row(0), row(1)};
    };
    const auto flow_at = [&](std::size_t i) {
        return cf::DiscreteKnownFlowInput{
            provider,       provider,    {stamp, times[i].view(), side(i)},
            ambient,        gust_limits, test_geometry::limits(),
            *geometry_work, flow_work};
    };
    const cf::DiscreteNonpenetrationLimits load_limits{
        {format, 1e-10},
        {format, 1e-10, 1e-10},
        {test_geometry::limits(), {ls::EstimateAdmission{65536}, 1000000000}, UINT64_C(5000000000)},
        100000,
        UINT64_C(24000000000)};
    a::SharedP1Workspace shared;
    auto local = take(a::LocalPacketWorkspace::create(memory));
    auto kernel = take(a::P1Workspace::create(memory));
    const std::array scales{1., 1.};
    const ls::UnitDimensions gamma_units{2, 0, -1, 0, 0, 0, 0};
    const ls::UnitDimensions speed_units{1, 0, -1, 0, 0, 0, 0};
    const std::array gamma{gamma_units, gamma_units}, speed{speed_units, speed_units};
    const auto request_at = [&](std::size_t i, const a::GlobalContactState& contact) {
        auto setup = circulation_request(
            {scales, scales, gamma, speed, circulation_context(geometry_at(contact).context())}, 2,
            ls::Norm::Infinity);
        setup.rhs_generation = core::RhsGeneration{31 + i};
        setup.retained_generation = core::NumericalGeneration{41 + i};
        return setup;
    };
    const auto history = [&](std::size_t i, const a::GlobalTopology& owner) {
        const auto contact = contact_at(i, owner);
        const auto queries = queries_at(contact);
        const auto rows = rows_at(i, queries);
        const cf::DiscreteNonpenetrationInput input{contact, rows,  {},     {parent, {}},
                                                    {50, i}, stamp, side(i)};
        return take(cf::solve_discrete_known_circulation(input, flow_at(i), std::monostate{},
                                                         load_limits, request_at(i, contact),
                                                         shared, *local, *kernel, memory));
    };
    const auto earlier = history(0, topology);
    const auto previous = history(1, topology);
    const auto t0 = node(0);
    const auto t1 = node(1);
    const auto prepare = [&](const cf::CaseClockedKnownFlowSample& one_node,
                             const cf::CaseClockedKnownFlowSample& two_node,
                             const cf::DiscreteKnownCirculation& prev,
                             const cf::DiscreteKnownCirculation& old) {
        return cf::FixedBdf2Circulation::prepare(node(2), one_node, two_node, prev, old, parent,
                                                 4096, 10000, memory);
    };
    reject(cf::FixedBdf2Circulation::prepare(node(1), t0, t0, earlier, earlier, parent, 4096, 10000,
                                             memory),
           core::FailureCode::StaleGeneration);
    reject(prepare(t0, t1, previous, earlier), core::FailureCode::StaleGeneration);
    reject(prepare(t1, t0, earlier, previous), core::FailureCode::StaleGeneration);
    const auto other_topology =
        take_global(a::build_global_topology(reference, roles, cap, *global, memory));
    const auto foreign = history(0, other_topology);
    reject(prepare(t1, t0, previous, foreign), core::FailureCode::StaleGeneration);
    const auto stencil = take(prepare(t1, t0, previous, earlier));

    const auto contact = contact_at(2, topology);
    const auto queries = queries_at(contact);
    const auto rows = rows_at(2, queries);
    const cf::DiscreteNonpenetrationInput input{contact, rows,  {},     {parent, {}},
                                                {50, 2}, stamp, side(2)};
    const auto request = request_at(2, contact);
    const auto no_motion = coordinates(0, 0, 0);
    const auto z_axis = coordinates(0, 0, 1);
    const std::array midpoints{cf::EdgeMidpointVelocity{a::EdgeId{10}, no_motion},
                               cf::EdgeMidpointVelocity{a::EdgeId{11}, no_motion},
                               cf::EdgeMidpointVelocity{a::EdgeId{12}, no_motion},
                               cf::EdgeMidpointVelocity{a::EdgeId{13}, no_motion},
                               cf::EdgeMidpointVelocity{a::EdgeId{15}, no_motion},
                               cf::EdgeMidpointVelocity{a::EdgeId{16}, no_motion},
                               cf::EdgeMidpointVelocity{a::EdgeId{17}, no_motion},
                               cf::EdgeMidpointVelocity{a::EdgeId{18}, no_motion}};
    cf::DiscreteVectorLoadInput load{
        input.trial_stamp, request.retained_generation, midpoints, {}, a::SingleWord{density},
        no_motion};
    load.fixed_bdf2 = &stencil;
    const auto solve = [&](const cf::DiscreteNonpenetrationInput& in,
                           const cf::DiscreteVectorLoadInput& loads,
                           const cf::DiscreteKnownFlowInput& flow,
                           const cf::DiscreteCirculationRequest& solve_request,
                           const cf::KnownLoadDerivative& derivative) {
        return cf::solve_discrete_known_loads(in, loads, flow, derivative, load_limits,
                                              solve_request, shared, *local, *kernel, memory);
    };
    const auto flow = flow_at(2);
    const auto run = [&](const cf::KnownLoadDerivative& derivative) {
        return take(solve(input, load, flow, request, derivative));
    };
    const auto primal = run(std::monostate{});

    std::array<Pair, 3> expected_gamma{};
    for (std::size_t i = 0; i < times.size(); ++i) {
        Pair rhs{};
        for (std::size_t panel = 0; panel < collocation.size(); ++panel)
            rhs[panel] = High{wall[i][panel]} - wind(collocation[panel][0], High{i} / 4);
        expected_gamma[i] = independent_gamma(rhs);
    }
    Pair rhs_direction{};
    for (std::size_t panel = 0; panel < collocation.size(); ++panel) {
        rhs_direction[panel] =
            High{wall_direction[panel]} +
            pi * cos(pi * (High{collocation[panel][0]} - High{.5})) * time_direction;
    }
    const auto expected_dgamma = independent_gamma(rhs_direction);
    Pair rate{};
    for (std::size_t panel = 0; panel < collocation.size(); ++panel) {
        rate[panel] = 2 * (3 * expected_gamma[2][panel] - 4 * expected_gamma[1][panel] +
                           expected_gamma[0][panel]);
        const High bdf1 = 4 * (expected_gamma[2][panel] - expected_gamma[1][panel]);
        REQUIRE(abs(rate[panel] - bdf1) > High{"1e-4"});
    }

    std::array<cf::CirculationCollocationDirection, 2> drows{
        {{fixture.faces[0].id,
          no_motion,
          {no_motion, no_motion},
          coordinates(0, 0, wall_direction[0])},
         {fixture.faces[1].id,
          no_motion,
          {no_motion, no_motion},
          coordinates(0, 0, wall_direction[1])}}};
    const std::array displacements{no_motion, no_motion, no_motion, no_motion,
                                   no_motion, no_motion, no_motion, no_motion};
    const cf::DiscreteCirculationDirection state{base.points, displacements, drows, no_motion};
    const cf::DiscreteImplicitVectorLoadDirection variation{
        load.stamp, load.circulation, state, midpoints, {}};
    auto tangent_request = request;
    tangent_request.rhs_generation = core::RhsGeneration{52};
    tangent_request.retained_generation = core::NumericalGeneration{62};
    const auto tangent =
        run(cf::KnownLoadTangent{variation, a::SingleWord{time_direction}, tangent_request});
    const cf::DiscreteLoadSampleWeights none{no_motion, no_motion};
    const std::array edge_seeds{none, none, none, none, none, none, none, none};
    const std::array facet_seeds{none, none, none, none};
    const cf::DiscreteVectorLoadWeights weights{
        load.stamp, load.circulation, edge_seeds, facet_seeds, displacements, z_axis, no_motion};
    const ls::UnitDimensions dual_gamma_units{-2, 0, 1, 0, 0, 0, 0};
    const ls::UnitDimensions dual_speed_units{-1, 0, 1, 0, 0, 0, 0};
    const std::array dual_gamma{dual_gamma_units, dual_gamma_units};
    const std::array dual_speed{dual_speed_units, dual_speed_units};
    auto dual_context = request.scales.context;
    dual_context.cache.matrix = core::NumericalGeneration{70};
    dual_context.cache.scales = core::ScaleGeneration{71};
    dual_context.unknown_layout = core::ScopeId{72};
    dual_context.residual_layout = core::ScopeId{73};
    auto transpose = circulation_request({scales, scales, dual_speed, dual_gamma, dual_context}, 2,
                                         ls::Norm::Infinity);
    transpose.rhs_generation = core::RhsGeneration{74};
    transpose.retained_generation = core::NumericalGeneration{75};
    const auto dual = run(cf::KnownLoadAdjoint{weights, transpose});

    REQUIRE(primal.bdf2_rate_count() == 2);
    REQUIRE(tangent.bdf2_rate_count() == 2);
    REQUIRE(dual.bdf2_rate_count() == 2);
    for (std::size_t panel = 0; panel < collocation.size(); ++panel) {
        const auto& circulation = take(primal.circulation()).get();
        expect(take(circulation.row(panel)).original_residual, High{0}, High{"1e-13"});
        const auto observed = take(primal.bdf2_rate(panel));
        CHECK(observed.face == fixture.faces[panel].id);
        expect(observed.rate, rate[panel]);
        expect(take(tangent.bdf2_rate(panel)).rate, rate[panel]);
        expect(take(dual.bdf2_rate(panel)).rate, rate[panel]);
        for (std::size_t facet = 2 * panel; facet < 2 * panel + 2; ++facet) {
            expect(take(primal.facet(facet)).force[2], density * rate[panel] / 2);
            expect(take(tangent.facet_direction(facet)).force[2],
                   3 * density * expected_dgamma[panel]);
        }
    }
    expect(take(primal.total_force())[2], density * (rate[0] + rate[1]));
    expect(take(primal.total_moment())[0], density * (rate[0] / 2 + 5 * rate[1] / 2));
    expect(take(primal.total_moment())[1], -density * (rate[0] + rate[1]) / 2);
    const High expected_vertical_direction =
        6 * density * (expected_dgamma[0] + expected_dgamma[1]);
    expect(take(tangent.force_direction())[2], expected_vertical_direction);
    const auto pullback_time = take(dual.flow_time_dual());
    REQUIRE(abs(center(pullback_time)) > High{".01"});
    High pairing = center(pullback_time) * time_direction;
    High expected_time_dual = 0;
    for (std::size_t panel = 0; panel < collocation.size(); ++panel) {
        Pair unit_rhs{};
        unit_rhs[panel] = 1;
        const auto action = independent_gamma(unit_rhs);
        const High expected_wall_dual = 6 * density * (action[0] + action[1]);
        const auto observed_wall_dual = take(dual.collocation_dual(panel)).velocity[2];
        expect(observed_wall_dual, expected_wall_dual);
        pairing += center(observed_wall_dual) * wall_direction[panel];
        expected_time_dual +=
            expected_wall_dual * pi * cos(pi * (High{collocation[panel][0]} - High{.5}));
    }
    expect(pullback_time, expected_time_dual);
    CHECK(abs(pairing - expected_vertical_direction) < High{"1e-11"});

    const PlanePoint midpoint{0., .5};
    const High edge_speed = wind(0., High{.5}) + influence(midpoint, 0) * expected_gamma[2][0] +
                            influence(midpoint, 1) * expected_gamma[2][1];
    expect(take(primal.edge(0)).force[0], -density * expected_gamma[2][0] * edge_speed);

    auto wrong_parent = input;
    wrong_parent.history.parent = core::AcceptedStateId{9};
    reject(solve(wrong_parent, load, flow, request, std::monostate{}),
           core::FailureCode::StaleGeneration);
    auto wrong_side = flow;
    wrong_side.time.side = p::RequestedSide::Left;
    auto wrong_side_input = input;
    wrong_side_input.flow_side = wrong_side.time.side;
    reject(solve(wrong_side_input, load, wrong_side, request, std::monostate{}),
           core::FailureCode::StaleGeneration);
    const auto wrong_time_contact = contact_at(1, topology);
    const auto wrong_time_queries = queries_at(wrong_time_contact);
    const auto wrong_time_rows = rows_at(1, wrong_time_queries);
    const cf::DiscreteNonpenetrationInput wrong_time{
        wrong_time_contact, wrong_time_rows, {}, {parent, {}}, input.trial_stamp, stamp, side(1)};
    auto wrong_time_request = request_at(1, wrong_time_contact);
    wrong_time_request.retained_generation = request.retained_generation;
    reject(solve(wrong_time, load, flow_at(1), wrong_time_request, std::monostate{}),
           core::FailureCode::StaleGeneration);
    auto missing_panel = input;
    missing_panel.collocations = missing_panel.collocations.first(1);
    auto partial_load = load;
    partial_load.edges = partial_load.edges.first(4);
    reject(solve(missing_panel, partial_load, flow, request, std::monostate{}),
           core::FailureCode::StaleGeneration);
    auto stale_request = request;
    stale_request.retained_generation = core::NumericalGeneration{42};
    auto stale_load = load;
    stale_load.circulation = stale_request.retained_generation;
    reject(solve(input, stale_load, flow, stale_request, std::monostate{}),
           core::FailureCode::StaleGeneration);

    // The stencil also rejects a reconstructed source that shares the user's SourceId.
    auto altered = source_case;
    altered.field("flow").field("gust").set("angular_frequency", parse("90"));
    const auto edited_si = materialize(altered, source_id);
    const auto edited =
        take(cf::prepare_case_known_flow(edited_si, stamp, {}, flow_work, gust_limits, memory));
    const cf::DiscreteKnownFlowInput wrong_source{
        provider,       edited,   flow.time, ambient, gust_limits, test_geometry::limits(),
        *geometry_work, flow_work};
    reject(solve(input, load, wrong_source, request, std::monostate{}),
           core::FailureCode::StaleGeneration);
}
