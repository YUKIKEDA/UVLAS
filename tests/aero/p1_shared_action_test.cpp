#include "aero/closed_contour.hpp"
#include "aero/p1_interval_private.hpp"
#include "aero/p1_shared_action.hpp"
#include "aero/p1_shared_kernel.hpp"
#include "geometry_snapshot_fixture.hpp"

#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <bit>
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
    bool shifted = false, aliased = false, empty = false;
    double radius = .125;

    void share_by_alias() noexcept {
        aliased = true;
        triangle[0].edges[0] = edges[7].id;
    }

    [[nodiscard]] p::GeometryInput input() const {
        if (empty)
            return {{}, {}, p::SinglePositions{{}}, p::SinglePositions{{}},  {}, {}, {},
                    {}, {}, word(radius),           test_geometry::context()};
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
        if (empty)
            return {{}, {}, {}};
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
    core::ResourcePlan memory;
    std::unique_ptr<p::GeometryWorkspace> geometry_workspace;
    p::OwnedGeometry geometry;
    std::unique_ptr<p::GlobalWorkspace> global_workspace;
    p::GlobalTopology topology;
    p::GlobalContactState contact;
    std::unique_ptr<p::LocalPacketWorkspace> local_workspace;
    std::unique_ptr<p::P1Workspace> kernel_workspace;
    p::SharedP1Workspace shared_workspace;

    template <class Fixture>
    explicit Prepared(const Fixture& input)
        : memory(128ULL * 1024 * 1024),
          geometry_workspace(test_geometry::workspace(memory)),
          geometry(test_geometry::take(p::prepare_owned_geometry(input.input(), limits().geometry,
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
            const auto* derivative = jvp ? &*jvp : nullptr;
            if (derivative == nullptr)
                FAIL("Missing JVP");
            for (std::size_t axis = 0; axis < 3; ++axis) {
                compare(value.velocity()[axis], reference.at("velocity").at(axis));
                compare((*derivative)[axis], reference.at("jvp").at(axis));
            }
            duality(value, transpose);
        }
    }
}

TEST_CASE("Gauss16 wake samples enclose the fixed age rule and complex phase",
          "[UW-D01][UW-D02][verification][negative]") {
    const auto context = test_geometry::context();
    const std::array boundaries{0., .125, .5};
    const auto make = [&](double frequency, std::uint64_t sample_limit, std::uint64_t work,
                          core::ResourcePlan& plan) {
        return p::prepare_straight_wake_age_quadrature(context, std::span{boundaries}.first<2>(),
                                                       {.5, .5}, frequency, sample_limit, work,
                                                       plan);
    };
    using Gauss = boost::math::quadrature::gauss<Big, 16>;
    const auto& nodes = Gauss::abscissa();
    const auto& weights = Gauss::weights();
    constexpr std::uint64_t sample_count = 32;
    constexpr auto bytes = sample_count * sizeof(p::StraightWakeAgeSample);
    const auto encloses = [](p::ContourInterval interval, const Big& target) {
        CHECK(Big(interval.lower) <= target);
        CHECK(target <= Big(interval.upper));
        CHECK(interval.upper - interval.lower < 1e-14);
    };
    for (const double frequency : {32., -32., 0., 1e-20}) {
        CAPTURE(frequency);
        core::ResourcePlan plan{bytes};
        const auto result = make(frequency, sample_count, 1'000'000, plan);
        REQUIRE(result);
        CHECK(plan.allocated_bytes() == bytes);
        REQUIRE(result->samples.size() == sample_count);
        for (std::size_t panel = 0; panel < boundaries.size() - 1; ++panel) {
            const Big left{boundaries[panel]}, right{boundaries[panel + 1]};
            const Big scale = (right - left) / 2;
            const Big midpoint = (left + right) / 2;
            for (std::size_t i = 0; i < sample_count / 2; ++i) {
                const auto positive = i >= nodes.size();
                const auto j = positive ? i - nodes.size() : nodes.size() - 1 - i;
                const Big age = midpoint + (positive ? nodes[j] : -nodes[j]) * scale;
                const Big quadrature = scale * weights[j];
                const Big phase = age * frequency;
                const auto& sample = result->samples.values()[panel * 16 + i];
                encloses(sample.age_seconds, age);
                encloses(sample.real_weight_seconds, quadrature * (cos(phase) - 1));
                encloses(sample.imaginary_weight_seconds, -quadrature * sin(phase));
                if (frequency == 0.) {
                    CHECK(sample.real_weight_seconds.lower == 0);
                    CHECK(sample.real_weight_seconds.upper == 0);
                    CHECK(sample.imaginary_weight_seconds.lower == 0);
                    CHECK(sample.imaginary_weight_seconds.upper == 0);
                }
            }
        }
        core::ResourcePlan exact{bytes};
        const auto replay = make(frequency, sample_count, result->interval_steps, exact);
        REQUIRE(replay);
        CHECK(replay->interval_steps == result->interval_steps);
        core::ResourcePlan short_steps{bytes};
        rejected(make(frequency, sample_count, result->interval_steps - 1, short_steps),
                 core::FailureCode::ResourceExceeded);
    }
    core::ResourcePlan short_samples{bytes};
    rejected(make(32., sample_count - 1, 1'000'000, short_samples),
             core::FailureCode::ResourceExceeded);
    CHECK(short_samples.allocated_bytes() == 0);
    core::ResourcePlan short_storage{bytes - 1};
    rejected(make(32., sample_count, 1'000'000, short_storage),
             core::FailureCode::ResourceExceeded);
    CHECK(short_storage.allocated_bytes() == 0);
    const std::array invalid_ages{0., 0.};
    core::ResourcePlan invalid{bytes};
    rejected(p::prepare_straight_wake_age_quadrature(context, invalid_ages, {.5, .5}, 32.,
                                                     sample_count, 1'000'000, invalid),
             core::FailureCode::InvalidInput);
    CHECK(invalid.allocated_bytes() == 0);
}

TEST_CASE("Straight wake age update agrees with the finite-core strip integral",
          "[UW-D01][UW-D02][UW-D07][UI-D05][verification][negative]") {
    const auto vector_words = [](const std::array<double, 3>& point) {
        return p::ContourVector{{{point[0], point[0]}, {point[1], point[1]}, {point[2], point[2]}}};
    };
    const auto origin = vector_words({0, -.5, 0});
    const auto tip = vector_words({.25, .5, .125});
    const auto velocity = vector_words({2, 0, 0});
    const auto normal = vector_words({.25, -.5, 1.});
    constexpr p::ContourInterval radius{.125, .125};
    const std::array<p::StraightWakeAgeSample, 3> samples{{
        {{0, 0}, {-.125, -.125}, {-.25, -.25}},
        {{.125, .125}, {-.5, -.5}, {.125, .125}},
        {{.75, .75}, {-.25, -.25}, {-.125, -.125}},
    }};
    const auto context = test_geometry::context();
    const auto projection = [&](const auto& first, const auto& second, const auto& flow, auto core,
                                const auto& point, auto nodes, std::uint64_t work,
                                std::stop_token stop = {}) {
        return p::project_straight_wake_column_update(context, first, second, flow, core, point,
                                                      normal, nodes, work, stop);
    };

    for (const std::array point : {std::array{-.125, 0., .125}, std::array{-.25, -.75, .0625},
                                   std::array{8., 32., -16.}, std::array{0., -.5, 0.}}) {
        CAPTURE(point);
        const auto projected =
            projection(origin, tip, velocity, radius, vector_words(point), samples, samples.size());
        REQUIRE(projected);
        const auto reversed =
            projection(tip, origin, velocity, radius, vector_words(point), samples, samples.size());
        REQUIRE(reversed);
        std::array<Big, 2> reference{};
        for (const auto& sample : samples) {
            const Big age{sample.age_seconds.lower};
            // Apply the Green-function Hessian over the infinitesimal strip.
            // dX/ds cross dX/d(age) = (0, 1/4, -2) for this oriented band.
            const auto density = [&](const Big& s) -> Big {
                const Big x = Big(point[0]) - s / 4 - 2 * age;
                const Big y = Big(point[1]) + Big(.5) - s;
                const Big z = Big(point[2]) - s / 8;
                const Big r2 = x * x + y * y + z * z;
                const Big epsilon2 = Big(1) / 64;
                const Big a2 = r2 + epsilon2;
                const Big projected_normal = Big(-17) / 8;
                const Big dot_r = x / 4 - y / 2 + z;
                const Big area_dot_r = y / 4 - 2 * z;
                return (projected_normal * (r2 - 2 * epsilon2) - 3 * dot_r * area_dot_r) /
                       (4 * boost::math::constants::pi<Big>() * a2 * a2 * sqrt(a2));
            };
            const auto integrate = [&](unsigned count) {
                Big total{};
                for (unsigned i = 0; i < count; ++i)
                    total += boost::math::quadrature::gauss<Big, 30>::integrate(
                        density, Big(i) / count, Big(i + 1) / count);
                return total;
            };
            const Big area = integrate(32);
            REQUIRE(abs(area - integrate(64)) < Big("1e-35"));
            reference[0] += area * sample.real_weight_seconds.lower;
            reference[1] += area * sample.imaginary_weight_seconds.lower;
        }
        const auto encloses = [](p::ContourInterval interval, const Big& target) {
            CHECK(Big(interval.lower) <= target);
            CHECK(target <= Big(interval.upper));
            CHECK(interval.upper - interval.lower < 1e-10);
        };
        encloses(projected->real, reference[0]);
        encloses(projected->imaginary, reference[1]);
        encloses(reversed->real, -reference[0]);
        encloses(reversed->imaginary, -reference[1]);
    }

    const auto point = vector_words({0, 0, 0});
    rejected(projection(origin, tip, velocity, radius, point, samples, samples.size() - 1),
             core::FailureCode::ResourceExceeded);
    rejected(projection(origin, origin, velocity, radius, point, samples, samples.size()),
             core::FailureCode::InvalidInput);
    rejected(
        projection(origin, tip, velocity, p::ContourInterval{}, point, samples, samples.size()),
        core::FailureCode::InvalidInput);
    auto bad = samples;
    bad[1].age_seconds = {-1, -.5};
    rejected(projection(origin, tip, velocity, radius, point, bad, bad.size()),
             core::FailureCode::InvalidInput);
    const std::array<p::StraightWakeAgeSample, 1> zero{{{{.25, .25}, {}, {}}}};
    const auto zero_update = projection(origin, tip, velocity, radius, point, zero, zero.size());
    REQUIRE(zero_update);
    CHECK(zero_update->real.lower == 0);
    CHECK(zero_update->real.upper == 0);
    CHECK(zero_update->imaginary.lower == 0);
    CHECK(zero_update->imaginary.upper == 0);
    std::stop_source stop;
    REQUIRE(stop.request_stop());
    rejected(
        projection(origin, tip, velocity, radius, point, samples, samples.size(), stop.get_token()),
        core::FailureCode::Cancelled);
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

TEST_CASE("P1 shared preparation rejects invalid accuracy independently of facet count",
          "[UI-D01][UI-D04][N4-02][negative]") {
    constexpr auto two_word = core::Precision::TwoWord;
    constexpr auto invalid_precision = static_cast<core::Precision>(255);
    constexpr auto nan = std::numeric_limits<double>::quiet_NaN();
    const p::LocalAccuracy local{two_word, 1e-10};
    const p::P1Accuracy kernel{two_word, 1e-10, 1e-10};
    struct AccuracyCase {
        const char* name;
        p::LocalAccuracy packet;
        p::P1Accuracy p1;
        core::FailureCode expected;
    };
    using enum core::FailureCode;
    const std::array cases{
        AccuracyCase{"local nonfinite", {two_word, nan}, kernel, Nonfinite},
        AccuracyCase{"local negative", {two_word, -1}, kernel, InvalidInput},
        AccuracyCase{"local representation", {invalid_precision, 0}, kernel, InvalidInput},
        AccuracyCase{"basis nonfinite", local, {two_word, nan, 0}, Nonfinite},
        AccuracyCase{"Jacobian nonfinite", local, {two_word, 0, nan}, Nonfinite},
        AccuracyCase{"basis negative", local, {two_word, -1, 0}, InvalidInput},
        AccuracyCase{"Jacobian negative", local, {two_word, 0, -1}, InvalidInput},
        AccuracyCase{
            "basis representation", local, {invalid_precision, 0, 0}, PrecisionUnsupported},
    };
    for (const bool empty : {false, true}) {
        GeometryFixture geometry;
        geometry.empty = empty;
        Prepared prepared{geometry};
        const auto query = prepared.query(geometry.query());
        const p::P1ActionSource source{prepared.contact, query, {12, 0}};
        const auto map = maps(geometry, true);
        const auto faces = empty ? std::span<const p::P1FaceStrengths>{}
                                 : std::span<const p::P1FaceStrengths>{map};
        const auto strengths =
            empty ? std::span<const p::P1StrengthId>{} : std::span<const p::P1StrengthId>{strength};
        const auto forecast =
            test_geometry::take(p::forecast_shared_p1(source, faces, strengths, limits()));
        for (const auto& entry : cases) {
            CAPTURE(empty, entry.name);
            core::ResourcePlan memory{forecast.payload_bytes};
            rejected(p::prepare_shared_p1(source, faces, strengths, entry.packet, entry.p1,
                                          limits(), prepared.shared_workspace,
                                          *prepared.local_workspace, *prepared.kernel_workspace,
                                          memory),
                     entry.expected);
            CHECK(memory.allocated_bytes() == 0);
        }
    }
}

TEST_CASE("P1 empty contact topology is an exact zero action",
          "[UI-D01][UI-D04][N4-02][property]") {
    GeometryFixture geometry;
    geometry.empty = true;
    Prepared prepared{geometry};
    const auto query = prepared.query(geometry.query());
    const p::P1ActionSource source{prepared.contact, query, {12, 0}};
    const auto forecast = test_geometry::take(p::forecast_shared_p1(source, {}, {}, limits()));
    core::ResourcePlan memory{forecast.payload_bytes};
    const auto action = test_geometry::take(p::prepare_shared_p1(
        source, {}, {}, {core::Precision::TwoWord, 0}, {core::Precision::TwoWord, 0, 0}, limits(),
        prepared.shared_workspace, *prepared.local_workspace, *prepared.kernel_workspace, memory));
    CHECK(memory.allocated_bytes() == forecast.payload_bytes);
    CHECK(action.point_count() == 0);
    CHECK(action.strength_count() == 0);
    CHECK(test_geometry::take(action.report()).facets == 0);
    const p::P1ActionDirection direction{{}, {}, vector(3, 4, 5), word(6)};
    const auto result = test_geometry::take(
        p::apply_shared_p1(action, source, {}, direction, apply_steps, prepared.shared_workspace));
    const auto& jvp = result.jvp();
    if (!jvp)
        FAIL("Missing zero JVP");
    const auto vjp = test_geometry::take(p::adjoint_shared_p1(
        action, source, {}, dual, apply_steps, prepared.shared_workspace, prepared.memory));
    for (std::size_t axis = 0; axis < 3; ++axis) {
        CHECK(decode(result.velocity()[axis].real).center == 0);
        CHECK(decode(result.velocity()[axis].real).radius == 0);
        CHECK(decode((*jvp)[axis].real).center == 0);
        CHECK(decode((*jvp)[axis].imaginary).radius == 0);
        const auto covector = test_geometry::take(vjp.free_query(axis));
        if (!covector)
            FAIL("Missing zero query cotangent");
        CHECK(decode(*covector).center == 0);
    }
    CHECK(decode(test_geometry::take(vjp.radius())).center == 0);
}

namespace {
using BigVector = std::array<Big, 3>;
BigVector minus(const BigVector& a, const BigVector& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
Big dot(const BigVector& a, const BigVector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
BigVector cross(const BigVector& a, const BigVector& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// Independent dimensional segment formula from analytic Biot–Savart integration.
// It uses the unnormalized endpoints directly, and has no sheet, facet, or scatter inputs.
BigVector ring(const GeometryFixture& geometry, double z, const Big& fraction, double x = .375,
               double y = .4375) {
    const BigVector query{Big{x} + fraction * high(direction_query[0]),
                          Big{y} + fraction * high(direction_query[1]),
                          Big{z} + fraction * high(direction_query[2])};
    const Big radius = Big{geometry.radius} + fraction * Big{.015625};
    const Big pi = acos(Big{-1});
    BigVector result{};
    for (std::size_t i = 0; i < 4; ++i) {
        BigVector a{}, b{};
        const auto next = (i + 1) % 4;
        for (std::size_t j = 0; j < 3; ++j) {
            a[j] = Big{geometry.coordinates[3 * i + j]} + fraction * high(direction_points[i][j]);
            b[j] = Big{geometry.coordinates[3 * next + j]} +
                   fraction * high(direction_points[next][j]);
        }
        const auto chord = minus(b, a);
        const auto start = minus(query, a), end = minus(query, b);
        const auto numerator = cross(chord, start);
        const Big denominator = dot(numerator, numerator) + radius * radius * dot(chord, chord);
        const Big weight = (dot(chord, start) / sqrt(dot(start, start) + radius * radius) -
                            dot(chord, end) / sqrt(dot(end, end) + radius * radius)) /
                           (4 * pi * denominator);
        for (std::size_t j = 0; j < 3; ++j)
            result[j] += numerator[j] * weight;
    }
    return result;
}
} // namespace

TEST_CASE("Mixed induction follows the bound circulation contour and the wake scatter",
          "[UW-D01][UW-D02][UW-D07][UI-D05][N4-06][verification]") {
    GeometryFixture geometry;
    geometry.roles[0].role = p::AeroFaceRole::Bound;
    geometry.share_by_alias();
    Prepared prepared{geometry};
    const auto face = geometry.faces[0].id;
    const auto wake_same_number = p::P1StrengthId{face.value()};
    const auto circulation = p::P1StrengthId::bound(face);
    CHECK(wake_same_number != circulation);
    CHECK(wake_same_number.wake_value() == face.value());
    CHECK(circulation.bound_face() == face);
    const std::array ids{strength[0], wake_same_number, circulation};
    const std::array wake_ids{ids[0], ids[1], p::P1StrengthId{101}};
    const std::array bound_wake{
        p::P1FaceStrengths{geometry.faces[0].id, p::P0BoundColumn{}},
        p::P1FaceStrengths{geometry.faces[1].id, std::array{ids[0], ids[1], ids[0]}}};
    const std::array wake_only{
        p::P1FaceStrengths{geometry.faces[0].id, std::array{wake_ids[2], wake_ids[2]}},
        bound_wake[1]};
    const auto prepare_mixed = [&](const p::P1ActionSource& source, auto maps, auto keys,
                                   core::ResourcePlan& memory,
                                   const p::P1AssemblyLimits& caps = limits()) {
        return p::prepare_bound_wake_induction(
            source, maps, keys, {core::Precision::TwoWord, 1e-10},
            {core::Precision::TwoWord, 1e-10, 1e-10}, caps, prepared.shared_workspace,
            *prepared.local_workspace, *prepared.kernel_workspace, memory);
    };
    for (const double z : {.75, 12.}) {
        CAPTURE(z);
        const auto query = prepared.query(geometry.query(z));
        const p::P1ActionSource source{prepared.contact, query, {3, 5}};
        const auto forecast = test_geometry::take(
            p::forecast_bound_wake_induction(source, bound_wake, ids, limits()));
        auto tight = limits();
        const auto kernel_steps = p::P1Workspace::maximum_steps();
        REQUIRE(kernel_steps > 0);
        REQUIRE(kernel_steps < tight.kernel.max_steps);
        tight.kernel.max_steps = kernel_steps;
        tight.geometry.max_steps = UINT64_C(500000000);
        tight.steps = forecast.reserved_steps;
        const auto at_limit =
            test_geometry::take(p::forecast_bound_wake_induction(source, bound_wake, ids, tight));
        CHECK(at_limit.reserved_steps == forecast.reserved_steps);
        core::ResourcePlan memory{forecast.payload_bytes};
        const auto mixed =
            test_geometry::take(prepare_mixed(source, bound_wake, ids, memory, tight));
        const auto report = test_geometry::take(mixed.report());
        CHECK(report.used_steps <= forecast.reserved_steps);
        --tight.kernel.max_steps;
        rejected(p::forecast_bound_wake_induction(source, bound_wake, ids, tight),
                 core::FailureCode::ResourceExceeded);
        CHECK(memory.allocated_bytes() == forecast.payload_bytes);
        CHECK(report.far_facets == (z > 1 ? 1 : 0));
        CHECK(report.cancelled_analytic_pairs == 1);
        CHECK(test_geometry::take(mixed.strength(1)) == wake_same_number);
        CHECK(test_geometry::take(mixed.strength(2)) == circulation);
        if (z < 1) {
            auto explicit_tie = bound_wake;
            explicit_tie[1].columns = std::array{circulation, circulation, ids[0]};
            const auto tied =
                test_geometry::take(prepare_mixed(source, explicit_tie, ids, prepared.memory));
            CHECK(test_geometry::take(tied.report()).cancelled_analytic_pairs == 2);
        }
        rejected(p::forecast_bound_wake_induction(source, wake_only, wake_ids, limits()),
                 core::FailureCode::InvalidInput);
        rejected(p::forecast_bound_wake_induction(source, bound_wake, std::span{ids}.first<2>(),
                                                  limits()),
                 core::FailureCode::InvalidInput);
        rejected(p::forecast_shared_p1(source, bound_wake, ids, limits()),
                 core::FailureCode::InvalidInput);
        core::ResourcePlan short_memory{forecast.payload_bytes - 1};
        rejected(prepare_mixed(source, bound_wake, ids, short_memory),
                 core::FailureCode::ResourceExceeded);
        CHECK(short_memory.allocated_bytes() == 0);

        const auto wake = test_geometry::take(p::prepare_shared_p1(
            source, wake_only, wake_ids, {core::Precision::TwoWord, 1e-10},
            {core::Precision::TwoWord, 1e-10, 1e-10}, limits(), prepared.shared_workspace,
            *prepared.local_workspace, *prepared.kernel_workspace, prepared.memory));
        auto wake_values = mu, wake_deltas = direction_mu;
        wake_values[2] = wake_deltas[2] = scalar(0);
        const p::P1ActionDirection direction{direction_points, direction_mu, direction_query,
                                             word(.015625)};
        const p::P1ActionDirection wake_direction{direction_points, wake_deltas, direction_query,
                                                  word(.015625)};
        const auto value = test_geometry::take(p::apply_shared_p1(
            mixed, source, mu, direction, apply_steps, prepared.shared_workspace));
        const auto background = test_geometry::take(p::apply_shared_p1(
            wake, source, wake_values, wake_direction, apply_steps, prepared.shared_workspace));
        const auto transpose = test_geometry::take(p::adjoint_shared_p1(
            mixed, source, mu, dual, apply_steps, prepared.shared_workspace, prepared.memory));
        duality(value, transpose);
        const Big h{"1e-25"};
        const auto center = ring(geometry, z, Big{0});
        const auto plus = ring(geometry, z, h), neg = ring(geometry, z, -h);
        const auto& jvp = value.jvp();
        const auto& wake_jvp = background.jvp();
        if (!jvp || !wake_jvp)
            FAIL("Missing mixed derivative");
        for (std::size_t axis = 0; axis < center.size(); ++axis) {
            const Big derivative = (plus[axis] - neg[axis]) / (2 * h);
            const auto check = [](const p::ScaledEnclosure& total, const p::ScaledEnclosure& base,
                                  const Big& expected) {
                const auto total_bound = decode(total), wake_bound = decode(base);
                CHECK(abs(total_bound.center - wake_bound.center - expected) <=
                      total_bound.radius + wake_bound.radius + Big{"1e-45"});
                CHECK(total_bound.radius < Big{"1e-15"});
            };
            check(value.velocity()[axis].real, background.velocity()[axis].real,
                  center[axis] * high(mu[2].real));
            check(value.velocity()[axis].imaginary, background.velocity()[axis].imaginary,
                  center[axis] * high(mu[2].imaginary));
            check((*jvp)[axis].real, (*wake_jvp)[axis].real,
                  derivative * high(mu[2].real) + center[axis] * high(direction_mu[2].real));
            check((*jvp)[axis].imaginary, (*wake_jvp)[axis].imaginary,
                  derivative * high(mu[2].imaginary) +
                      center[axis] * high(direction_mu[2].imaginary));
        }
    }
}

TEST_CASE("A constant wake quad uses an analytic contour at near and far queries",
          "[UW-D01][UW-D02][UW-D07][UI-D05][N4-06][verification]") {
    GeometryFixture geometry;
    geometry.share_by_alias();
    Prepared prepared{geometry};
    const std::array maps{
        p::P1FaceStrengths{geometry.faces[0].id, p::ConstantWakeColumn{strength[2]}},
        p::P1FaceStrengths{geometry.faces[1].id,
                           std::array{strength[0], strength[1], strength[0]}}};
    auto values = mu, deltas = direction_mu;
    values[0] = values[1] = deltas[0] = deltas[1] = scalar(0);
    const p::P1ActionDirection direction{direction_points, deltas, direction_query, word(.015625)};
    for (const double z : {.75, 12.}) {
        CAPTURE(z);
        const auto query = prepared.query(geometry.query(z));
        const p::P1ActionSource source{prepared.contact, query, {3, 5}};
        const auto prepare = [&](const auto& columns) {
            return test_geometry::take(p::prepare_bound_wake_induction(
                source, columns, strength, {core::Precision::TwoWord, 1e-10},
                {core::Precision::TwoWord, 1e-10, 1e-10}, limits(), prepared.shared_workspace,
                *prepared.local_workspace, *prepared.kernel_workspace, prepared.memory));
        };
        const auto closed = prepare(maps);
        const auto result = test_geometry::take(p::apply_shared_p1(
            closed, source, values, direction, apply_steps, prepared.shared_workspace));
        const auto transpose = test_geometry::take(p::adjoint_shared_p1(
            closed, source, values, dual, apply_steps, prepared.shared_workspace, prepared.memory));
        auto variable = maps;
        variable[0].columns = std::array{strength[2], strength[2]};
        if (z < 1) {
            const auto baseline = prepare(variable);
            const auto original = test_geometry::take(p::apply_shared_p1(
                baseline, source, values, direction, apply_steps, prepared.shared_workspace));
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const auto a = decode(result.velocity()[axis].real);
                const auto b = decode(original.velocity()[axis].real);
                CHECK(abs(a.center - b.center) <= a.radius + b.radius);
            }
        }
        const auto value = ring(geometry, z, Big{0});
        const Big h{"1e-25"};
        const auto above = ring(geometry, z, h), below = ring(geometry, z, -h);
        const auto& tangent = result.jvp();
        const auto* tangent_values = tangent ? &*tangent : nullptr;
        if (tangent_values == nullptr)
            FAIL("Missing constant contour derivative");
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const auto actual = decode(result.velocity()[axis].real);
            const auto derivative = decode((*tangent_values)[axis].real);
            const Big expected = value[axis] * high(values[2].real);
            const Big expected_derivative =
                (above[axis] - below[axis]) * high(values[2].real) / (2 * h) +
                value[axis] * high(deltas[2].real);
            CHECK(abs(actual.center - expected) <= actual.radius + Big{"1e-45"});
            CHECK(abs(derivative.center - expected_derivative) <= derivative.radius + Big{"1e-45"});
        }
        const auto bound = test_geometry::take(transpose.strength(2));
        const auto d = decode(bound.real);
        Big expected = 0;
        for (std::size_t axis = 0; axis < 3; ++axis)
            expected += value[axis] * high(dual[axis].real);
        CHECK(abs(d.center - expected) <= d.radius + Big{"1e-45"});
    }
}

TEST_CASE("Mixed induction reuses unowned facet packets and bases",
          "[UW-D01][UW-D07][N4-06][verification]") {
    GeometryFixture geometry;
    geometry.roles[0].role = p::AeroFaceRole::Bound;
    geometry.share_by_alias();
    Prepared prepared{geometry};
    const std::array ids{strength[0], strength[1], p::P1StrengthId::bound(geometry.faces[0].id)};
    const std::array maps{
        p::P1FaceStrengths{geometry.faces[0].id, p::P0BoundColumn{}},
        p::P1FaceStrengths{geometry.faces[1].id, std::array{ids[0], ids[1], ids[0]}}};
    const auto near = prepared.query(geometry.query(.75));
    const auto far = prepared.query(geometry.query(12.));
    const p::P1ActionSource near_source{prepared.contact, near, {3, 5}};
    const p::P1ActionSource far_source{prepared.contact, far, {3, 5}};
    const p::P1ActionDirection direction{direction_points, direction_mu, direction_query,
                                         word(.015625)};
    const auto near_plan =
        test_geometry::take(p::forecast_bound_wake_induction(near_source, maps, ids, limits()));
    const auto far_plan =
        test_geometry::take(p::forecast_bound_wake_induction(far_source, maps, ids, limits()));
    const auto serial = test_geometry::take(
        p::forecast_bound_wake_serial_induction(far_source, maps, ids, limits()));
    REQUIRE(serial.reusable);
    CHECK(serial.ordinary.payload_bytes == far_plan.payload_bytes);

    const auto recycled = [&](const p::P1ActionSource& src, p::BoundWakeFacetStorage& storage,
                              const p::P1AssemblyLimits& maximum,
                              const std::stop_token& token = {}) {
        return p::prepare_bound_wake_induction_reusing_facets(
            src, maps, ids, {core::Precision::TwoWord, 1e-10},
            {core::Precision::TwoWord, 1e-10, 1e-10}, maximum, prepared.shared_workspace,
            *prepared.local_workspace, *prepared.kernel_workspace, storage, token);
    };
    const auto bits = [](double x) { return std::bit_cast<std::uint64_t>(x); };
    const auto same = [&](const p::ScaledEnclosure& a, const p::ScaledEnclosure& b) {
        CHECK(a.exponent() == b.exponent());
        CHECK(bits(a.radius()) == bits(b.radius()));
        const auto a_words = a.center(), b_words = b.center();
        const auto* a_pair = std::get_if<p::DoubleWord>(&a_words);
        const auto* b_pair = std::get_if<p::DoubleWord>(&b_words);
        if (a_pair == nullptr || b_pair == nullptr)
            FAIL("Expected a two-word enclosure");
        CHECK(bits(a_pair->high) == bits(b_pair->high));
        CHECK(bits(a_pair->low) == bits(b_pair->low));
    };
    const auto same_complex = [&](const p::P1ComplexEnclosure& a, const p::P1ComplexEnclosure& b) {
        same(a.real, b.real);
        same(a.imaginary, b.imaginary);
    };

    for (const auto* source : {&near_source, &far_source}) {
        const auto budget =
            source == &near_source ? near_plan.payload_bytes : far_plan.payload_bytes;
        core::ResourcePlan coefficient_memory{3 * budget};
        const auto full = test_geometry::take(p::prepare_bound_wake_induction(
            *source, maps, ids, {core::Precision::TwoWord, 1e-10},
            {core::Precision::TwoWord, 1e-10, 1e-10}, limits(), prepared.shared_workspace,
            *prepared.local_workspace, *prepared.kernel_workspace, coefficient_memory));
        p::BoundWakeFacetStorage coefficient_storage{coefficient_memory};
        {
            const auto values = test_geometry::take(p::prepare_bound_wake_values(
                *source, maps, ids, {core::Precision::TwoWord, 1e-10},
                {core::Precision::TwoWord, 1e-10, 1e-10}, limits(), prepared.shared_workspace,
                *prepared.local_workspace, *prepared.kernel_workspace, coefficient_memory,
                &coefficient_storage, {}));
            const auto full_report = test_geometry::take(full.report());
            const auto value_report = test_geometry::take(values.report());
            CHECK(value_report.used_steps < full_report.used_steps);
            CHECK(value_report.far_facets == full_report.far_facets);
            const auto ordinary = test_geometry::take(p::apply_shared_p1(
                full, *source, mu, std::nullopt, apply_steps, prepared.shared_workspace));
            const auto coefficient_velocity = test_geometry::take(p::apply_shared_p1(
                values, *source, mu, std::nullopt, apply_steps, prepared.shared_workspace));
            for (std::size_t axis = 0; axis < ordinary.velocity().size(); ++axis)
                same_complex(ordinary.velocity()[axis], coefficient_velocity.velocity()[axis]);
            std::array<p::p1::Interval, 9> ordinary_coefficients{}, value_coefficients{};
            REQUIRE(p::project_real_strength_transpose(full, *source, ordinary_coefficients,
                                                       apply_steps, prepared.shared_workspace, {}));
            REQUIRE(p::project_real_strength_transpose(values, *source, value_coefficients,
                                                       apply_steps, prepared.shared_workspace, {}));
            for (std::size_t i = 0; i < ordinary_coefficients.size(); ++i) {
                CHECK(p::p1::compare(ordinary_coefficients[i].lower, value_coefficients[i].lower) ==
                      0);
                CHECK(p::p1::compare(ordinary_coefficients[i].upper, value_coefficients[i].upper) ==
                      0);
            }
            rejected(p::apply_shared_p1(values, *source, mu, direction, apply_steps,
                                        prepared.shared_workspace),
                     core::FailureCode::InvalidInput);
            rejected(p::adjoint_shared_p1(values, *source, mu, dual, apply_steps,
                                          prepared.shared_workspace, coefficient_memory),
                     core::FailureCode::InvalidInput);
        }
        const auto restored = test_geometry::take(recycled(*source, coefficient_storage, limits()));
        REQUIRE(p::apply_shared_p1(restored, *source, mu, direction, apply_steps,
                                   prepared.shared_workspace));
    }

    core::ResourcePlan memory{2 * (near_plan.payload_bytes + far_plan.payload_bytes)};
    p::BoundWakeFacetStorage storage{memory};
    {
        const auto previous = [&] {
            const auto near_action = test_geometry::take(recycled(near_source, storage, limits()));
            return test_geometry::take(p::apply_shared_p1(near_action, near_source, mu, direction,
                                                          apply_steps, prepared.shared_workspace));
        }();
        const auto first_charge = memory.allocated_bytes();
        CHECK(first_charge < near_plan.payload_bytes);
        rejected(recycled(far_source, storage, limits()), core::FailureCode::InvalidInput);
        CHECK(memory.allocated_bytes() == first_charge);
        CHECK(previous.used_steps() > 0);
    }
    const auto first_charge = memory.allocated_bytes();
    CHECK(first_charge == serial.initial_payload_bytes);
    auto step_too_short = limits();
    step_too_short.steps = far_plan.reserved_steps;
    rejected(recycled(far_source, storage, step_too_short), core::FailureCode::ResourceExceeded);
    CHECK(memory.allocated_bytes() == first_charge);
    // NOLINTNEXTLINE(misc-const-correctness): MSVC 14.51 declares request_stop non-const.
    std::stop_source cancelled;
    REQUIRE(cancelled.request_stop());
    rejected(recycled(far_source, storage, limits(), cancelled.get_token()),
             core::FailureCode::Cancelled);
    CHECK(memory.allocated_bytes() == first_charge);
    std::uint64_t reuse_charge = 0;
    {
        const auto reused = test_geometry::take(recycled(far_source, storage, limits()));
        const auto report = test_geometry::take(reused.report());
        reuse_charge = memory.allocated_bytes() - first_charge;
        CAPTURE(first_charge, reuse_charge);
        CHECK(reuse_charge == report.plan.payload_bytes);
        CHECK(reuse_charge == serial.repeated_payload_bytes);
        CHECK(report.plan.reserved_steps == serial.repeated_reserved_steps);
        REQUIRE(far_plan.payload_bytes > reuse_charge);
        CHECK(far_plan.payload_bytes - reuse_charge >
              report.facets *
                  (p::LocalGeometryPacket::payload_bytes() + p::P1Basis::payload_bytes() +
                   p::ExactTriangleResult::payload_bytes()));
        CHECK(report.plan.reserved_steps > far_plan.reserved_steps);
        CHECK(report.far_facets == 1);

        core::ResourcePlan fresh_memory{far_plan.payload_bytes};
        const auto fresh = test_geometry::take(p::prepare_bound_wake_induction(
            far_source, maps, ids, {core::Precision::TwoWord, 1e-10},
            {core::Precision::TwoWord, 1e-10, 1e-10}, limits(), prepared.shared_workspace,
            *prepared.local_workspace, *prepared.kernel_workspace, fresh_memory));
        CHECK(report.used_steps < test_geometry::take(fresh.report()).used_steps);
        CAPTURE(fresh_memory.allocated_bytes());
        CHECK(fresh_memory.allocated_bytes() == far_plan.payload_bytes);
        const auto evaluate = [&](const p::SharedP1Action& action) {
            return test_geometry::take(p::apply_shared_p1(action, far_source, mu, direction,
                                                          apply_steps, prepared.shared_workspace));
        };
        const auto compare_adjoint = [&](const p::SharedP1Action& action) {
            return test_geometry::take(p::adjoint_shared_p1(action, far_source, mu, dual,
                                                            apply_steps, prepared.shared_workspace,
                                                            prepared.memory));
        };
        const auto original = evaluate(fresh);
        const auto current = evaluate(reused);
        const auto& original_jvp = original.jvp();
        const auto& current_jvp = current.jvp();
        if (!original_jvp || !current_jvp)
            FAIL("Missing simultaneous directional derivative");
        const auto original_dual = compare_adjoint(fresh);
        const auto current_dual = compare_adjoint(reused);
        for (std::size_t axis = 0; axis < 3; ++axis) {
            same_complex(original.velocity()[axis], current.velocity()[axis]);
            same_complex((*original_jvp)[axis], (*current_jvp)[axis]);
            const auto a = test_geometry::take(original_dual.free_query(axis));
            const auto b = test_geometry::take(current_dual.free_query(axis));
            if (!a || !b)
                FAIL("Missing free-query covector");
            same(*a, *b);
            for (std::size_t point = 0; point < fresh.point_count(); ++point)
                same(test_geometry::take(original_dual.point(point, axis)),
                     test_geometry::take(current_dual.point(point, axis)));
        }
        for (std::size_t i = 0; i < ids.size(); ++i)
            same_complex(test_geometry::take(original_dual.strength(i)),
                         test_geometry::take(current_dual.strength(i)));
        same(test_geometry::take(original_dual.radius()),
             test_geometry::take(current_dual.radius()));

        std::vector<p::p1::Interval> projected(3 * ids.size());
        rejected(p::project_real_strength_transpose(reused, far_source,
                                                    std::span{projected}.subspan(1), apply_steps,
                                                    prepared.shared_workspace, {}),
                 core::FailureCode::InvalidInput);
        rejected(p::project_real_strength_transpose(reused, far_source, projected, 1,
                                                    prepared.shared_workspace, {}),
                 core::FailureCode::ResourceExceeded);
        REQUIRE(p::project_real_strength_transpose(reused, far_source, projected, apply_steps,
                                                   prepared.shared_workspace, {}));
        for (std::size_t axis = 0; axis < 3; ++axis) {
            std::array<p::P1ComplexWords, 3> unit{};
            unit[axis].real = p::SingleWord{1};
            const auto full = test_geometry::take(
                p::adjoint_shared_p1(fresh, far_source, mu, unit, apply_steps,
                                     prepared.shared_workspace, prepared.memory));
            for (std::size_t i = 0; i < ids.size(); ++i) {
                const auto expected = decode(test_geometry::take(full.strength(i)).real);
                const auto& actual = projected[3 * i + axis];
                CHECK(ldexp(Big{actual.lower.mantissa}, actual.lower.power) ==
                      expected.center - expected.radius);
                CHECK(ldexp(Big{actual.upper.mantissa}, actual.upper.power) ==
                      expected.center + expected.radius);
            }
        }
    }
    const auto before_second_reuse = memory.allocated_bytes();
    const auto back = test_geometry::take(recycled(near_source, storage, limits()));
    CHECK(test_geometry::take(back.report()).far_facets == 0);
    CHECK(memory.allocated_bytes() - before_second_reuse == reuse_charge);

    core::ResourcePlan byte_short{first_charge + reuse_charge - 1};
    p::BoundWakeFacetStorage short_storage{byte_short};
    {
        const auto admitted = test_geometry::take(recycled(near_source, short_storage, limits()));
        CHECK(admitted.active());
    }
    rejected(recycled(far_source, short_storage, limits()), core::FailureCode::ResourceExceeded);
    CHECK(byte_short.allocated_bytes() == first_charge);

    core::ResourcePlan exact{first_charge + reuse_charge};
    p::BoundWakeFacetStorage exact_storage{exact};
    {
        const auto admitted = test_geometry::take(recycled(near_source, exact_storage, limits()));
        CHECK(admitted.active());
    }
    const auto completed = test_geometry::take(recycled(far_source, exact_storage, limits()));
    CHECK(completed.active());
    CHECK(exact.remaining_bytes() == 0);

    core::ResourcePlan retry_memory{3 * (near_plan.payload_bytes + far_plan.payload_bytes)};
    p::BoundWakeFacetStorage retry_storage{retry_memory};
    {
        const auto admitted = test_geometry::take(recycled(near_source, retry_storage, limits()));
        CHECK(admitted.active());
    }
    const auto before_refusal = retry_memory.allocated_bytes();
    rejected(p::prepare_bound_wake_induction_reusing_facets(
                 far_source, maps, ids, {core::Precision::TwoWord, 0},
                 {core::Precision::TwoWord, 1e-10, 1e-10}, limits(), prepared.shared_workspace,
                 *prepared.local_workspace, *prepared.kernel_workspace, retry_storage),
             core::FailureCode::PrecisionUnresolved);
    CHECK(retry_memory.allocated_bytes() > before_refusal);
    const auto before_retry = retry_memory.allocated_bytes();
    const auto retried = test_geometry::take(recycled(far_source, retry_storage, limits()));
    CHECK(test_geometry::take(retried.report()).far_facets == 1);
    CHECK(retry_memory.allocated_bytes() - before_retry == reuse_charge);

    // Keep the same IDs and diagnostic revisions but change an actual vertex.
    auto updated = geometry;
    updated.coordinates[8] = .5;
    Prepared updated_geometry{updated};
    const auto updated_query = updated_geometry.query(updated.query(12.));
    const p::P1ActionSource updated_source{updated_geometry.contact, updated_query, {3, 5}};
    const auto same_size_fresh = [&](p::BoundWakeFacetStorage& arena) {
        return p::prepare_bound_wake_induction_reusing_facets(
            updated_source, maps, ids, {core::Precision::TwoWord, 1e-10},
            {core::Precision::TwoWord, 1e-10, 1e-10}, limits(), updated_geometry.shared_workspace,
            *updated_geometry.local_workspace, *updated_geometry.kernel_workspace, arena);
    };
    core::ResourcePlan changed_memory{2 * (near_plan.payload_bytes + far_plan.payload_bytes)};
    p::BoundWakeFacetStorage changed_storage{changed_memory}, last_storage{changed_memory};
    {
        const auto first = test_geometry::take(recycled(near_source, changed_storage, limits()));
        CHECK(first.active());
    }
    const auto changed = test_geometry::take(same_size_fresh(changed_storage));
    const auto cold = test_geometry::take(same_size_fresh(last_storage));
    const auto warm_report = test_geometry::take(changed.report());
    const auto cold_report = test_geometry::take(cold.report());
    // Storage reuse has its own initialization visits. The new snapshot must
    // pay the same numerical preparation costs as the empty storage.
    CHECK(warm_report.used_steps >= cold_report.used_steps);
}

TEST_CASE("Serial bound midpoints exclude each new self edge and retain its wake alias",
          "[UW-D01][UW-D07][N4-06][negative][verification]") {
    GeometryFixture geometry;
    geometry.roles[0].role = p::AeroFaceRole::Bound;
    geometry.share_by_alias();
    Prepared prepared{geometry};
    const std::array ids{strength[0], strength[1], p::P1StrengthId::bound(geometry.faces[0].id)};
    const std::array maps{
        p::P1FaceStrengths{geometry.faces[0].id, p::P0BoundColumn{}},
        p::P1FaceStrengths{geometry.faces[1].id, std::array{ids[0], ids[1], ids[0]}}};
    const auto first = geometry.edges[0].id, second = geometry.edges[1].id;
    const auto midpoint = [&](p::EdgeId edge) {
        return prepared.query({p::QueryId{edge.value()}, p::SourceEdgePoint{edge, word(.5)}});
    };
    const auto first_query = midpoint(first), second_query = midpoint(second);
    const p::P1ActionSource first_source{prepared.contact, first_query, {3, 5}};
    const p::P1ActionSource second_source{prepared.contact, second_query, {3, 5}};
    const auto forecast = test_geometry::take(
        p::forecast_bound_wake_midpoint_serial_induction(first_source, maps, ids, limits()));
    REQUIRE(forecast.reusable);
    const auto cap = forecast.initial_payload_bytes + forecast.repeated_payload_bytes;
    CHECK(cap < 2 * forecast.ordinary.payload_bytes);
    const auto make = [&](const p::P1ActionSource& source, p::EdgeId edge, auto& target) {
        return p::prepare_bound_wake_midpoint_induction(
            source, edge, maps, ids, {core::Precision::TwoWord, 1e-10},
            {core::Precision::TwoWord, 1e-10, 1e-10}, limits(), prepared.shared_workspace,
            *prepared.local_workspace, *prepared.kernel_workspace, target);
    };
    const auto velocity = [&](const p::SharedP1Action& action, const p::P1ActionSource& source) {
        return test_geometry::take(p::apply_shared_p1(action, source, mu, std::nullopt, apply_steps,
                                                      prepared.shared_workspace));
    };

    core::ResourcePlan memory{cap};
    p::BoundWakeFacetStorage storage{memory};
    rejected(make(first_source, second, storage), core::FailureCode::InvalidInput);
    CHECK(memory.allocated_bytes() == 0);
    for (const auto& [edge, source] :
         {std::pair{first, &first_source}, std::pair{second, &second_source}}) {
        CAPTURE(edge.value());
        core::ResourcePlan fresh_memory{forecast.ordinary.payload_bytes};
        const auto fresh = test_geometry::take(make(*source, edge, fresh_memory));
        const auto reused = test_geometry::take(make(*source, edge, storage));
        const auto report = test_geometry::take(reused.report());
        REQUIRE(report.excluded_bound_occurrences == 1);
        CHECK(report.plan.reserved_steps == (edge == first ? forecast.ordinary.reserved_steps
                                                           : forecast.repeated_reserved_steps));
        const auto expected = velocity(fresh, *source);
        const auto observed = velocity(reused, *source);
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const auto original = decode(expected.velocity()[axis].real);
            const auto current = decode(observed.velocity()[axis].real);
            CHECK(current.center == original.center);
            CHECK(current.radius == original.radius);
        }
        if (edge == first) {
            auto wake_strengths = mu;
            wake_strengths.back() = {};
            const auto wake = test_geometry::take(
                p::apply_shared_p1(reused, *source, wake_strengths, std::nullopt, apply_steps,
                                   prepared.shared_workspace));
            CHECK(std::ranges::any_of(wake.velocity(), [](const auto& v) {
                const auto enclosure = decode(v.real);
                return abs(enclosure.center) > enclosure.radius;
            }));
        }
    }
    CHECK(memory.remaining_bytes() == 0);
    core::ResourcePlan short_memory{cap - 1};
    p::BoundWakeFacetStorage short_storage{short_memory};
    {
        const auto first_result = test_geometry::take(make(first_source, first, short_storage));
        CHECK(first_result.active());
    }
    rejected(make(second_source, second, short_storage), core::FailureCode::ResourceExceeded);
    CHECK(short_memory.allocated_bytes() == forecast.initial_payload_bytes);
}

TEST_CASE("Bound collocation may lie on a planar panel and its artificial diagonal",
          "[UW-D01][UW-D02][UW-D07][verification]") {
    for (const double radius : {0., .125}) {
        GeometryFixture geometry;
        geometry.radius = radius;
        geometry.coordinates[8] = 0;
        geometry.roles[0].role = p::AeroFaceRole::Bound;
        Prepared prepared{geometry};
        const std::array ids{p::P1StrengthId{5}, p::P1StrengthId::bound(geometry.faces[0].id)};
        const std::array bound_wake{
            p::P1FaceStrengths{geometry.faces[0].id, p::P0BoundColumn{}},
            p::P1FaceStrengths{geometry.faces[1].id, std::array{ids[0], ids[0], ids[0]}}};
        for (const double x : {.375, .5, 1.}) {
            CAPTURE(radius, x);
            const auto query =
                prepared.query({p::QueryId{29}, p::FreePoint{p::FreeSingle{{x, .5, 0.}}}});
            const p::P1ActionSource source{prepared.contact, query, {3, 5}};
            const auto forecast = test_geometry::take(
                p::forecast_bound_wake_induction(source, bound_wake, ids, limits()));
            core::ResourcePlan memory{forecast.payload_bytes};
            auto result = p::prepare_bound_wake_induction(
                source, bound_wake, ids, {core::Precision::TwoWord, 1e-10},
                {core::Precision::TwoWord, 1e-10, 1e-10}, limits(), prepared.shared_workspace,
                *prepared.local_workspace, *prepared.kernel_workspace, memory);
            if (x == 1.) {
                rejected(result, core::FailureCode::GeometryInvalid);
                continue;
            }
            const auto action = test_geometry::take(std::move(result));
            const auto value = test_geometry::take(
                p::apply_shared_p1(action, source, std::array{scalar(0), scalar(1)}, std::nullopt,
                                   apply_steps, prepared.shared_workspace));
            const auto expected = ring(geometry, 0., Big{0}, x, .5);
            for (std::size_t axis = 0; axis < expected.size(); ++axis) {
                const auto bound = decode(value.velocity()[axis].real);
                CHECK(abs(bound.center - expected[axis]) <= bound.radius + Big{"1e-45"});
            }
        }
    }
}

namespace {
struct TwoBoundPanels {
    GeometryFixture base;
    std::array<p::PointId, 6> points{p::PointId{1}, p::PointId{2}, p::PointId{3},
                                     p::PointId{4}, p::PointId{5}, p::PointId{6}};
    std::array<p::PointBinding, 6> bindings{test_geometry::binding(), test_geometry::binding(),
                                            test_geometry::binding(), test_geometry::binding(),
                                            test_geometry::binding(), test_geometry::binding()};
    std::array<double, 18> coordinates{0, 0, 0, 0,  1, 0,     1,  1, .25,
                                       1, 0, 0, -1, 0, -.125, -1, 1, .125};
    std::array<p::EdgeInput, 9> edges{base.edges[0],
                                      base.edges[1],
                                      base.edges[2],
                                      base.edges[3],
                                      base.edges[4],
                                      p::EdgeInput{p::EdgeId{15}, {points[0], points[4]}},
                                      p::EdgeInput{p::EdgeId{16}, {points[4], points[5]}},
                                      p::EdgeInput{p::EdgeId{17}, {points[5], points[1]}},
                                      p::EdgeInput{p::EdgeId{18}, {points[1], points[4]}}};
    std::array<p::QuadInput, 2> quads{
        base.quad[0], p::QuadInput{base.faces[1].id,
                                   {points[1], points[0], points[4], points[5]},
                                   {edges[0].id, edges[5].id, edges[6].id, edges[7].id},
                                   edges[8].id,
                                   {p::FacetId{202}, p::FacetId{203}}}};
    std::array<p::EdgeBindingInput, 9> owned{p::EdgeBindingInput{edges[0].id, p::CanonicalEdge{}},
                                             p::EdgeBindingInput{edges[1].id, p::CanonicalEdge{}},
                                             p::EdgeBindingInput{edges[2].id, p::CanonicalEdge{}},
                                             p::EdgeBindingInput{edges[3].id, p::CanonicalEdge{}},
                                             p::EdgeBindingInput{edges[4].id, p::CanonicalEdge{}},
                                             p::EdgeBindingInput{edges[5].id, p::CanonicalEdge{}},
                                             p::EdgeBindingInput{edges[6].id, p::CanonicalEdge{}},
                                             p::EdgeBindingInput{edges[7].id, p::CanonicalEdge{}},
                                             p::EdgeBindingInput{edges[8].id, p::CanonicalEdge{}}};
    std::array<p::FaceRoleInput, 2> roles{
        p::FaceRoleInput{base.faces[0].id, p::AeroFaceRole::Bound},
        p::FaceRoleInput{base.faces[1].id, p::AeroFaceRole::Bound}};
    [[nodiscard]] p::GeometryInput input() const {
        const p::SinglePositions xyz{coordinates};
        return {points,
                bindings,
                xyz,
                xyz,
                edges,
                base.faces,
                base.bands,
                {},
                quads,
                word(base.radius),
                test_geometry::context()};
    }
    [[nodiscard]] p::GlobalTopologyInput topology() const {
        return {roles, owned, {}};
    }
};
} // namespace

TEST_CASE("Two bound circulation columns share the signed edge kernel",
          "[UW-D01][UW-D02][UW-D07][UI-D05][verification]") {
    const TwoBoundPanels geometry;
    Prepared prepared{geometry};
    const auto query = prepared.query(geometry.base.query());
    const p::P1ActionSource source{prepared.contact, query, {3, 5}};
    const std::array ids{p::P1StrengthId::bound(geometry.base.faces[0].id),
                         p::P1StrengthId::bound(geometry.base.faces[1].id)};
    const std::array bound{p::P1FaceStrengths{geometry.base.faces[0].id, p::P0BoundColumn{}},
                           p::P1FaceStrengths{geometry.base.faces[1].id, p::P0BoundColumn{}}};
    auto resources = limits();
    const auto bound_steps = p::P1Workspace::maximum_bound_steps();
    REQUIRE(bound_steps < p::P1Workspace::maximum_steps());
    resources.kernel.max_steps = bound_steps;
    const auto forecast =
        test_geometry::take(p::forecast_bound_wake_induction(source, bound, ids, resources));
    resources.steps = forecast.reserved_steps;
    auto insufficient = resources;
    --insufficient.kernel.max_steps;
    rejected(p::forecast_bound_wake_induction(source, bound, ids, insufficient),
             core::FailureCode::ResourceExceeded);
    const auto action = test_geometry::take(p::prepare_bound_wake_induction(
        source, bound, ids, {core::Precision::TwoWord, 1e-10},
        {core::Precision::TwoWord, 1e-10, 1e-10}, resources, prepared.shared_workspace,
        *prepared.local_workspace, *prepared.kernel_workspace, prepared.memory));
    const auto report = test_geometry::take(action.report());
    CHECK(report.cancelled_analytic_pairs == 2);
    CHECK(report.used_steps <= forecast.reserved_steps);
    const auto first = ring(geometry.base, .75, Big{0});
    // Re-pack the second quad's vertices in its declared clockwise loop.
    auto second_geometry = geometry.base;
    for (std::size_t i = 0; i < 4; ++i) {
        constexpr std::array<std::size_t, 4> loop{1, 0, 4, 5};
        for (std::size_t axis = 0; axis < 3; ++axis)
            second_geometry.coordinates[3 * i + axis] = geometry.coordinates[3 * loop[i] + axis];
    }
    const auto second = ring(second_geometry, .75, Big{0});
    const auto values =
        std::array{p::P1ComplexWords{p::DoubleWord{1., 0x1p-60}, word(0)}, scalar(1)};
    const auto answer = test_geometry::take(p::apply_shared_p1(
        action, source, values, std::nullopt, apply_steps, prepared.shared_workspace));
    auto core_resources = core::ResourcePlan{1024ULL * 1024};
    const auto transpose = test_geometry::take(
        p::adjoint_shared_p1(action, source, values, {scalar(0), scalar(0), scalar(1)}, apply_steps,
                             prepared.shared_workspace, core_resources));
    for (std::size_t axis = 0; axis < first.size(); ++axis) {
        const auto velocity = decode(answer.velocity()[axis].real);
        const Big expected = first[axis] * (Big{1} + Big{0x1p-60}) + second[axis];
        CHECK(abs(velocity.center - expected) <= velocity.radius + Big{"1e-45"});
    }
    const auto first_gradient = decode(test_geometry::take(transpose.strength(0)).real);
    const auto second_gradient = decode(test_geometry::take(transpose.strength(1)).real);
    CHECK(abs(first_gradient.center - first[2]) <= first_gradient.radius + Big{"1e-45"});
    CHECK(abs(second_gradient.center - second[2]) <= second_gradient.radius + Big{"1e-45"});
}
