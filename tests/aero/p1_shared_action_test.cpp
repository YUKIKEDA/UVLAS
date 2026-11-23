#include "aero/p1_shared_action.hpp"
#include "geometry_snapshot_fixture.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <fstream>

namespace {
using namespace uvlas;
namespace p = aero::detail;
using Big = boost::multiprecision::cpp_bin_float_100;
using Words = std::array<p::ScalarWords, 3>;

template <class T>
T global(p::GlobalResult<T> result) {
    if (!result)
        FAIL("Unexpected global failure " << static_cast<int>(result.error().origin.cause.code));
    return std::move(*result);
}
template <class T>
void rejected(const p::AeroResult<T>& result, core::FailureCode code) {
    if (result)
        FAIL("Unexpected shared P1 success");
    CHECK(result.error().cause.code == code);
    CHECK(result.error().cause.operation == core::Operation::EvaluateAero);
}

p::ScalarWords word(double a) noexcept {
    return p::SingleWord{a};
}
Words vector(double x, double y, double z) noexcept {
    return {word(x), word(y), word(z)};
}
p::P1ComplexWords scalar(double re, double im = 0) noexcept {
    return {word(re), word(im)};
}

struct GeometryFixture {
    std::array<p::PointId, 5> point{p::PointId{1}, p::PointId{2}, p::PointId{3}, p::PointId{4},
                                    p::PointId{5}};
    std::array<p::PointBinding, 5> bindings{test_geometry::binding(), test_geometry::binding(),
                                            test_geometry::binding(), test_geometry::binding(),
                                            test_geometry::binding()};
    std::array<double, 15> coordinates{0, 0, 0, 0, 1, 0, 1, 1, .25, 1, 0, 0, -.75, .5, -.125};
    std::array<double, 15> offset{};
    std::array<p::EdgeInput, 8> edges{p::EdgeInput{p::EdgeId{10}, {point[0], point[1]}},
                                      p::EdgeInput{p::EdgeId{11}, {point[1], point[2]}},
                                      p::EdgeInput{p::EdgeId{12}, {point[2], point[3]}},
                                      p::EdgeInput{p::EdgeId{13}, {point[3], point[0]}},
                                      p::EdgeInput{p::EdgeId{14}, {point[0], point[2]}},
                                      p::EdgeInput{p::EdgeId{15}, {point[1], point[4]}},
                                      p::EdgeInput{p::EdgeId{16}, {point[4], point[0]}},
                                      p::EdgeInput{p::EdgeId{17}, {point[1], point[0]}}};
    std::array<p::EdgeBindingInput, 8> edge_bindings{
        p::EdgeBindingInput{edges[0].id, p::CanonicalEdge{}},
        p::EdgeBindingInput{edges[1].id, p::CanonicalEdge{}},
        p::EdgeBindingInput{edges[2].id, p::CanonicalEdge{}},
        p::EdgeBindingInput{edges[3].id, p::CanonicalEdge{}},
        p::EdgeBindingInput{edges[4].id, p::CanonicalEdge{}},
        p::EdgeBindingInput{edges[5].id, p::CanonicalEdge{}},
        p::EdgeBindingInput{edges[6].id, p::CanonicalEdge{}},
        p::EdgeBindingInput{edges[7].id, edges[0].id}};
    std::array<p::BandId, 2> bands{p::BandId{50}, p::BandId{51}};
    std::array<p::FaceInput, 2> faces{p::FaceInput{p::PhysicalFaceId{100}, bands[0]},
                                      p::FaceInput{p::PhysicalFaceId{101}, bands[1]}};
    std::array<p::FaceRoleInput, 2> roles{p::FaceRoleInput{faces[0].id, p::AeroFaceRole::Wake},
                                          p::FaceRoleInput{faces[1].id, p::AeroFaceRole::Wake}};
    std::array<p::QuadInput, 1> quad{
        p::QuadInput{faces[0].id,
                     {point[0], point[1], point[2], point[3]},
                     {edges[0].id, edges[1].id, edges[2].id, edges[3].id},
                     edges[4].id,
                     {p::FacetId{200}, p::FacetId{201}}}};
    std::array<p::TriangleInput, 1> triangle{
        p::TriangleInput{p::FacetId{202},
                         faces[1].id,
                         {point[0], point[1], point[4]},
                         {edges[0].id, edges[5].id, edges[6].id}}};
    bool shifted = false, aliased = false;
    double radius = .125;

    void share_by_alias() noexcept {
        aliased = true;
        triangle[0].edges[0] = edges[7].id;
    }

    [[nodiscard]] p::GeometryInput input() const {
        const p::PositionView words = shifted
                                          ? p::PositionView{p::DoublePositions{offset, coordinates}}
                                          : p::PositionView{p::SinglePositions{coordinates}};
        return {point,
                bindings,
                words,
                words,
                std::span<const p::EdgeInput>{edges}.first(aliased ? 8 : 7),
                faces,
                bands,
                triangle,
                quad,
                word(radius),
                test_geometry::context()};
    }
    [[nodiscard]] p::GlobalTopologyInput topology() const {
        return {
            roles, std::span<const p::EdgeBindingInput>{edge_bindings}.first(aliased ? 8 : 7), {}};
    }
    [[nodiscard]] p::QueryRequest query(double z = .75) const {
        if (shifted)
            return {p::QueryId{29}, p::FreePoint{p::FreeDouble{{offset[0], offset[1], offset[2]},
                                                               {.375, .4375, z}}}};
        return {p::QueryId{29}, p::FreePoint{p::FreeSingle{{.375, .4375, z}}}};
    }
};

struct Bound {
    Big center, radius;
};
Bound decode(const p::ScaledEnclosure& value) {
    const Big center = std::visit(
        [](const auto& data) -> Big {
            using T = std::decay_t<decltype(data)>;
            if constexpr (std::is_same_v<T, p::SingleWord>)
                return Big{data.value};
            else
                return Big{data.high} + data.low;
        },
        value.center());
    return {ldexp(center, value.exponent()), ldexp(Big{value.radius()}, value.exponent())};
}
void compare(const p::P1ComplexEnclosure& value, const nlohmann::json& independent) {
    // GL48 reference, successive GL36→48 maximum change 1.91e-25.
    const Big tolerance{"1e-22"};
    const auto real = decode(value.real);
    const auto imag = decode(value.imaginary);
    CHECK(abs(real.center - Big{independent.at(0).get<std::string>()}) <= tolerance + real.radius);
    CHECK(abs(imag.center - Big{independent.at(1).get<std::string>()}) <= tolerance + imag.radius);
    CHECK(real.radius < Big{"1e-12"});
    CHECK(imag.radius < Big{"1e-12"});
}
void include_dot(Bound& sum, const Bound& bound, double coefficient) {
    sum.center += bound.center * coefficient;
    sum.radius += bound.radius * abs(Big{coefficient});
}

p::P1AssemblyLimits limits() {
    return {test_geometry::limits(),
            {numerics::linear::EstimateAdmission{65536}, UINT64_C(1000000000)},
            UINT64_C(4000000000)};
}
constexpr std::uint64_t apply_steps = 100000;
const std::array<p::P1StrengthId, 3> strength{p::P1StrengthId{5}, p::P1StrengthId{8},
                                              p::P1StrengthId{13}};
const std::array direction_points{vector(.0625, 0, -.03125), vector(0, .0625, .03125),
                                  vector(0, -.03125, .0625), vector(.03125, -.0625, 0),
                                  vector(-.03125, 0, .0625)};
const Words direction_query = vector(-.0625, .03125, .0625);
const std::array direction_mu{scalar(.125, -.125), scalar(.25, .0625), scalar(-.125, .125)};
const std::array mu{scalar(.75, .25), scalar(-.375, .125), scalar(.5, -.25)};
const std::array dual{scalar(.25, -.0625), scalar(-.5, .125), scalar(.0625, .375)};

struct Prepared {
    core::ResourcePlan memory{128ULL * 1024 * 1024};
    std::unique_ptr<p::GeometryWorkspace> geometry_workspace = test_geometry::workspace(memory);
    p::OwnedGeometry geometry;
    std::unique_ptr<p::GlobalWorkspace> global_workspace;
    p::GlobalTopology topology;
    p::GlobalContactState contact;
    std::unique_ptr<p::LocalPacketWorkspace> local_workspace;
    std::unique_ptr<p::P1Workspace> kernel_workspace;
    p::SharedP1Workspace shared_workspace;

    explicit Prepared(const GeometryFixture& input)
        : geometry(test_geometry::take(p::prepare_owned_geometry(input.input(), limits().geometry,
                                                                 *geometry_workspace, memory))),
          global_workspace(global_workspace_create(memory)),
          topology(global(p::build_global_topology(geometry, input.topology(), global_limits(),
                                                   *global_workspace, memory))),
          contact(global(p::check_global_geometry(topology, geometry, global_limits(),
                                                  *global_workspace, memory))),
          local_workspace(required(p::LocalPacketWorkspace::create(memory))),
          kernel_workspace(required(p::P1Workspace::create(memory))) {}

    template <class T>
    static T required(core::Result<T> result) {
        if (!result)
            FAIL("Unexpected workspace failure " << static_cast<int>(result.error().code));
        return std::move(*result);
    }
    static std::unique_ptr<p::GlobalWorkspace> global_workspace_create(core::ResourcePlan& memory) {
        return required(
            p::GlobalWorkspace::create(memory, numerics::linear::EstimateAdmission{4096}));
    }
    static p::GlobalLimits global_limits() {
        return {std::numeric_limits<std::uint64_t>::max() / 4,
                numerics::linear::EstimateAdmission{4096}};
    }
    p::ResolvedQuery query(const p::QueryRequest& query) {
        return test_geometry::take(p::resolve_query(
            geometry, query, geometry.context(), limits().geometry, *geometry_workspace, memory));
    }
    p::AeroResult<p::SharedP1Action> prepare(const p::P1ActionSource& source,
                                             std::span<const p::P1FaceStrengths> faces,
                                             core::ResourcePlan& resources,
                                             p::P1AssemblyLimits cap = limits(),
                                             const std::stop_token& stop = {}) {
        return p::prepare_shared_p1(source, faces, strength, {core::Precision::TwoWord, 1e-10},
                                    {core::Precision::TwoWord, 1e-10, 1e-10}, cap, shared_workspace,
                                    *local_workspace, *kernel_workspace, resources, stop);
    }
};

std::array<p::P1FaceStrengths, 2> maps(const GeometryFixture& fixture, bool equal) {
    const auto [left, right, third] = strength;
    return {p::P1FaceStrengths{fixture.faces[0].id, std::array{left, right}},
            p::P1FaceStrengths{fixture.faces[1].id, equal ? std::array{left, left, third}
                                                          : std::array{right, left, third}}};
}

double high(const p::ScalarWords& value) {
    return std::get<p::SingleWord>(value).value;
}

void duality(const p::P1ActionValue& value, const p::P1ActionAdjoint& transpose) {
    const auto& direction = value.jvp();
    if (!direction)
        FAIL("Missing JVP");
    Bound primal{}, adjoint{};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        include_dot(primal, decode((*direction)[axis].real), high(dual[axis].real));
        include_dot(primal, decode((*direction)[axis].imaginary), high(dual[axis].imaginary));
        const auto query = test_geometry::take(transpose.free_query(axis));
        if (!query)
            FAIL("Missing free-query cotangent");
        include_dot(adjoint, decode(*query), high(direction_query[axis]));
        for (std::size_t point = 0; point < direction_points.size(); ++point)
            include_dot(adjoint, decode(test_geometry::take(transpose.point(point, axis))),
                        high(direction_points[point][axis]));
    }
    for (std::size_t index = 0; index < strength.size(); ++index) {
        const auto bound = test_geometry::take(transpose.strength(index));
        include_dot(adjoint, decode(bound.real), high(direction_mu[index].real));
        include_dot(adjoint, decode(bound.imaginary), high(direction_mu[index].imaginary));
    }
    include_dot(adjoint, decode(test_geometry::take(transpose.radius())), .015625);
    CHECK(abs(primal.center - adjoint.center) <= primal.radius + adjoint.radius);
}
} // namespace

TEST_CASE("Shared P1 quad action and all derivatives match independent area quadrature",
          "[UI-D01][UI-D04][UI-D05][UI-D06][N4-06][verification]") {
    std::ifstream stream{UVLAS_P1_SHARED_FIXTURE};
    REQUIRE(stream.good());
    const auto fixture_oracle = nlohmann::json::parse(stream).at("cases");
    for (const unsigned mode : {0U, 1U, 2U}) {
        GeometryFixture geometry;
        geometry.shifted = mode == 1;
        if (geometry.shifted)
            geometry.offset.fill(0x1p900);
        if (mode == 2)
            geometry.share_by_alias();
        Prepared prepared{geometry};
        const auto query = prepared.query(geometry.query());
        const p::P1ActionSource source{prepared.contact, query, {7, 11}};
        for (const bool equal : {true, false}) {
            CAPTURE(mode, equal);
            const auto map = maps(geometry, equal);
            const auto forecast =
                test_geometry::take(p::forecast_shared_p1(source, map, strength, limits()));
            core::ResourcePlan memory{forecast.payload_bytes};
            const auto action = test_geometry::take(prepared.prepare(source, map, memory));
            const auto report = test_geometry::take(action.report());
            CHECK(memory.allocated_bytes() == forecast.payload_bytes);
            CHECK(report.facets == 3);
            CHECK(report.far_facets == 0);
            CHECK(report.cancelled_analytic_pairs == (equal ? 2 : 1));
            CHECK(report.used_steps <= report.plan.reserved_steps);
            CHECK(action.point_count() == direction_points.size());
            CHECK(action.strength_count() == strength.size());
            for (std::size_t index = 0; index < strength.size(); ++index)
                CHECK(test_geometry::take(action.strength(index)) == strength[index]);
            const p::P1ActionDirection direction{direction_points, direction_mu, direction_query,
                                                 word(.015625)};
            const auto value = test_geometry::take(p::apply_shared_p1(
                action, source, mu, direction, apply_steps, prepared.shared_workspace));
            const auto transpose = test_geometry::take(p::adjoint_shared_p1(
                action, source, mu, dual, apply_steps, prepared.shared_workspace, prepared.memory));
            const auto& reference = fixture_oracle.at(equal ? "equal_trace" : "unequal_trace");
            const auto& jvp = value.jvp();
            if (!jvp)
                FAIL("Missing JVP");
            for (std::size_t axis = 0; axis < 3; ++axis) {
                compare(value.velocity()[axis], reference.at("velocity").at(axis));
                compare((*jvp)[axis], reference.at("jvp").at(axis));
            }
            duality(value, transpose);
        }
    }
}

TEST_CASE("P1 applies real and low-word strengths while retaining distinct equal-valued columns",
          "[UI-D04][UI-D05][UI-D06][N4-06][property]") {
    const GeometryFixture geometry;
    Prepared prepared{geometry};
    const auto query = prepared.query(geometry.query());
    const p::P1ActionSource source{prepared.contact, query, {7, 11}};
    const auto equal_action =
        test_geometry::take(prepared.prepare(source, maps(geometry, true), prepared.memory));
    const auto unequal_action =
        test_geometry::take(prepared.prepare(source, maps(geometry, false), prepared.memory));
    const std::array real{scalar(1), scalar(1), scalar(.5)}, seed{scalar(1), scalar(-1), scalar(0)};
    const auto zero = vector(0, 0, 0);
    const std::array zero_points{zero, zero, zero, zero, zero};
    const p::P1ActionDirection direction{zero_points, seed, vector(0, 0, 0), word(0)};
    const auto evaluate = [&](const p::SharedP1Action& action,
                              const std::array<p::P1ComplexWords, 3>& strengths,
                              std::optional<p::P1ActionDirection> d = std::nullopt) {
        return test_geometry::take(p::apply_shared_p1(action, source, strengths, d, apply_steps,
                                                      prepared.shared_workspace));
    };
    const auto with_shared = evaluate(equal_action, real, direction);
    const auto with_distinct = evaluate(unequal_action, real, direction);
    const auto& shared_jvp = with_shared.jvp();
    const auto& distinct_jvp = with_distinct.jvp();
    if (!shared_jvp || !distinct_jvp)
        FAIL("Missing coefficient JVP");
    bool derivative_differs = false;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const auto a = decode(with_shared.velocity()[axis].real);
        const auto b = decode(with_distinct.velocity()[axis].real);
        CHECK(abs(a.center - b.center) <= a.radius + b.radius);
        CHECK(decode(with_distinct.velocity()[axis].imaginary).center == 0);
        CHECK(decode(with_distinct.velocity()[axis].imaginary).radius == 0);
        const auto first = decode((*shared_jvp)[axis].real);
        const auto second = decode((*distinct_jvp)[axis].real);
        derivative_differs |= abs(first.center - second.center) > first.radius + second.radius;
    }
    CHECK(derivative_differs); // Equal primal words must not erase an independent coefficient.

    auto low = real;
    low[0].real = p::DoubleWord{1, 0x1p-80};
    const auto retained = evaluate(equal_action, low);
    const auto unit = evaluate(equal_action, {scalar(1), scalar(0), scalar(0)});
    bool resolved_low_word = false;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const auto actual = decode(retained.velocity()[axis].real);
        const auto base = decode(with_shared.velocity()[axis].real);
        const auto basis = decode(unit.velocity()[axis].real);
        const auto change = actual.center - base.center;
        const auto combined_radius = actual.radius + base.radius;
        CHECK(abs(change - ldexp(basis.center, -80)) <= combined_radius + ldexp(basis.radius, -80));
        resolved_low_word |= abs(change) > combined_radius;
    }
    CHECK(resolved_low_word);
}

TEST_CASE("P1 material query scatter includes each shared point occurrence in JVP and VJP",
          "[UI-D04][UI-D05][N4-03][N4-06][property]") {
    const GeometryFixture geometry;
    Prepared prepared{geometry};
    const auto query =
        prepared.query({p::QueryId{31}, p::SourceEdgePoint{geometry.edges[0].id, word(.375)}});
    const p::P1ActionSource source{prepared.contact, query, {7, 11}};
    const auto action =
        test_geometry::take(prepared.prepare(source, maps(geometry, true), prepared.memory));
    const std::array shift{vector(.125, -.25, .5), vector(.125, -.25, .5), vector(.125, -.25, .5),
                           vector(.125, -.25, .5), vector(.125, -.25, .5)};
    const std::array zero_mu{scalar(0), scalar(0), scalar(0)};
    const p::P1ActionDirection direction{shift, zero_mu, std::nullopt, word(0)};
    const auto value = test_geometry::take(
        p::apply_shared_p1(action, source, mu, direction, apply_steps, prepared.shared_workspace));
    const auto& jvp = value.jvp();
    if (!jvp)
        FAIL("Missing translation JVP");
    const auto dual_value = test_geometry::take(p::adjoint_shared_p1(
        action, source, mu, dual, apply_steps, prepared.shared_workspace, prepared.memory));
    Bound total{};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        CHECK_FALSE(test_geometry::take(dual_value.free_query(axis)));
        const auto re = decode((*jvp)[axis].real);
        const auto im = decode((*jvp)[axis].imaginary);
        CHECK(abs(re.center) <= re.radius);
        CHECK(abs(im.center) <= im.radius);
        for (std::size_t row = 0; row < shift.size(); ++row)
            include_dot(total, decode(test_geometry::take(dual_value.point(row, axis))),
                        high(shift[row][axis]));
    }
    CHECK(abs(total.center) <= total.radius);
}

TEST_CASE("P1 retains the selected far formula and refuses incompatible replay and resources",
          "[UI-D02][UI-D05][N4-03][N4-08][N4-09][negative]") {
    const GeometryFixture geometry;
    Prepared prepared{geometry};
    const auto query = prepared.query(geometry.query(5));
    const p::P1ActionSource source{prepared.contact, query, {7, 11}};
    const auto map = maps(geometry, true);
    const auto plan = test_geometry::take(p::forecast_shared_p1(source, map, strength, limits()));

    auto work_short = limits();
    work_short.steps = plan.reserved_steps - 1;
    auto opaque_short = limits();
    opaque_short.kernel.external_memory =
        numerics::linear::EstimateAdmission{plan.estimated_metadata_bytes - 1};
    auto proven = limits();
    proven.kernel.external_memory = numerics::linear::RequireProvenBound{UINT64_MAX};
    for (const auto& cap : {work_short, opaque_short, proven}) {
        core::ResourcePlan memory{plan.payload_bytes};
        rejected(prepared.prepare(source, map, memory, cap),
                 std::holds_alternative<numerics::linear::RequireProvenBound>(
                     cap.kernel.external_memory.selection())
                     ? core::FailureCode::Unsupported
                     : core::FailureCode::ResourceExceeded);
        CHECK(memory.allocated_bytes() == 0);
    }
    core::ResourcePlan byte_short{plan.payload_bytes - 1};
    rejected(prepared.prepare(source, map, byte_short), core::FailureCode::ResourceExceeded);
    CHECK(byte_short.allocated_bytes() == 0);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 declares request_stop non-const.
    std::stop_source stopped;
    stopped.request_stop();
    core::ResourcePlan cancelled{plan.payload_bytes};
    rejected(prepared.prepare(source, map, cancelled, limits(), stopped.get_token()),
             core::FailureCode::Cancelled);
    CHECK(cancelled.allocated_bytes() == 0);

    auto equal_cap = limits();
    equal_cap.steps = plan.reserved_steps;
    equal_cap.kernel.external_memory =
        numerics::linear::EstimateAdmission{plan.estimated_metadata_bytes};
    core::ResourcePlan memory{plan.payload_bytes};
    const auto action = test_geometry::take(prepared.prepare(source, map, memory, equal_cap));
    const auto report = test_geometry::take(action.report());
    CHECK(report.far_facets == 1);
    CHECK(report.cancelled_analytic_pairs == 1);
    for (std::size_t facet = 0; facet < report.facets; ++facet) {
        const auto diagnostic = test_geometry::take(action.facet_diagnostic(facet));
        if (std::holds_alternative<p::P1QuadratureDifference>(diagnostic))
            CHECK(test_geometry::take(action.report()).far_facets == 1);
    }
    const auto wrong_query = prepared.query(geometry.query(5.125));
    const p::P1ActionSource wrong{prepared.contact, wrong_query, source.stamp};
    const p::P1ActionSource old{prepared.contact, query, {7, 10}};
    auto changed_core = geometry;
    changed_core.radius = .25;
    Prepared other{changed_core};
    const auto other_query = other.query(changed_core.query(5));
    const p::P1ActionSource other_core{other.contact, other_query, source.stamp};
    for (const auto& switched : {wrong, old, other_core})
        rejected(
            p::apply_shared_p1(action, switched, mu, {}, apply_steps, prepared.shared_workspace),
            core::FailureCode::StaleGeneration);
    const auto valid = test_geometry::take(
        p::apply_shared_p1(action, source, mu, {}, apply_steps, prepared.shared_workspace));
    CHECK_FALSE(valid.jvp());
    auto nonfinite = mu;
    nonfinite[1].imaginary = word(std::numeric_limits<double>::infinity());
    rejected(
        p::apply_shared_p1(action, source, nonfinite, {}, apply_steps, prepared.shared_workspace),
        core::FailureCode::Nonfinite);
    rejected(p::apply_shared_p1(action, source, mu, {}, 0, prepared.shared_workspace),
             core::FailureCode::ResourceExceeded);
    CHECK(test_geometry::take(action.report()).used_steps == report.used_steps);
    CHECK(decode(valid.velocity()[2].real).center != 0);
}
