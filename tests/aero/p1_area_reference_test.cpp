#include "aero/p1_induction.hpp"
#include "geometry_snapshot_fixture.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <nlohmann/json.hpp>

#include <fstream>

namespace {
using namespace uvlas;
using namespace uvlas::test_geometry;
using namespace subject;
using Real = boost::multiprecision::cpp_bin_float_100;
using Complex = boost::multiprecision::cpp_complex_100;
using ComplexVector = std::array<Complex, 3>;
Real center(const ScaledEnclosure& value) {
    const auto words = value.center();
    Real result = 0;
    if (const auto* one = std::get_if<SingleWord>(&words))
        result = Real{one->value};
    else {
        const auto& pair = std::get<DoubleWord>(words);
        result = Real{pair.high} + Real{pair.low};
    }
    return boost::multiprecision::ldexp(result, value.exponent());
}
Real radius(const ScaledEnclosure& value) {
    return boost::multiprecision::ldexp(Real{value.radius()}, value.exponent());
}
QueryDefinition definition(const nlohmann::json& query, const TriangleFixture& fixture) {
    const auto kind = query.at("kind").get<std::string>();
    if (kind == "free")
        return FreePoint{FreeSingle{query.at("point").get<std::array<double, 3>>()}};
    if (kind == "face")
        return SourceFacePoint{FacetId{70}, SingleWord{query.at("u").get<double>()},
                               SingleWord{query.at("v").get<double>()}};
    if (kind == "edge")
        return SourceEdgePoint{fixture.edges.at(query.at("edge_index").get<std::size_t>()).id,
                               SingleWord{query.at("t").get<double>()}};
    return SourceVertex{fixture.points.at(query.at("index").get<std::size_t>())};
}
LocalGeometryPacket make_packet(const TriangleFixture& fixture, QueryDefinition query,
                                core::ResourcePlan& plan) {
    auto source_workspace = workspace(plan);
    auto local_workspace = LocalPacketWorkspace::create(plan);
    REQUIRE(local_workspace);
    auto geometry =
        take(prepare_owned_geometry(fixture.input(), limits(), *source_workspace, plan));
    auto resolved = take(resolve_query(geometry, {QueryId{77}, query}, geometry.context(), limits(),
                                       *source_workspace, plan));
    return take(prepare_local_packet(geometry, resolved, FacetId{70},
                                     {core::Precision::TwoWord, 1e-10}, limits(), **local_workspace,
                                     plan));
}
void close(const Real& actual, const Real& expected) {
    // UI adoption protocol: all individual basis components and all ambient columns.
    const Real error = abs(actual - expected);
    CHECK((error <= Real{"1e-10"} || error <= Real{"1e-8"} * abs(expected)));
}
Complex dot(const ComplexVector& a, const ComplexVector& b) {
    Complex result = 0;
    for (std::size_t i = 0; i < 3; ++i)
        result += a[i] * b[i];
    return result;
}
ComplexVector cross(const ComplexVector& a, const ComplexVector& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
ComplexVector subtract(const ComplexVector& a, const ComplexVector& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
struct ReferenceRule {
    std::array<Real, 12> nodes, weights;
};
ReferenceRule reference_rule() {
    ReferenceRule result;
    const Real pi = boost::multiprecision::acos(Real{-1});
    for (unsigned i = 0; i < 12; ++i) {
        Real z = cos(pi * (Real{i + 1} - Real{".25"}) / Real{"12.5"});
        for (unsigned iteration = 0; iteration < 32; ++iteration) {
            Real a = 1, b = z;
            for (unsigned k = 2; k <= 12; ++k) {
                const Real next = ((2 * k - 1) * z * b - (k - 1) * a) / k;
                a = b;
                b = next;
            }
            const Real derivative = 12 * (z * b - a) / (z * z - 1);
            z -= b / derivative;
        }
        Real a = 1, b = z;
        for (unsigned k = 2; k <= 12; ++k) {
            const Real next = ((2 * k - 1) * z * b - (k - 1) * a) / k;
            a = b;
            b = next;
        }
        const Real derivative = 12 * (z * b - a) / (z * z - 1);
        result.nodes[11 - i] = (z + 1) / 2;
        result.weights[11 - i] = 1 / ((1 - z * z) * derivative * derivative);
    }
    return result;
}
std::array<Complex, 9> fixed_area(const std::array<Complex, 13>& values,
                                  const ReferenceRule& rule) {
    const ComplexVector p{values[0], values[1], values[2]};
    const std::array<ComplexVector, 3> x{{{values[3], values[4], values[5]},
                                          {values[6], values[7], values[8]},
                                          {values[9], values[10], values[11]}}};
    const auto c = cross(subtract(x[1], x[0]), subtract(x[2], x[0]));
    const Complex jac = sqrt(dot(c, c));
    ComplexVector normal;
    for (std::size_t axis = 0; axis < 3; ++axis)
        normal[axis] = c[axis] / jac;
    const Real pi = boost::multiprecision::acos(Real{-1});
    std::array<Complex, 9> output{};
    for (std::size_t i = 0; i < 12; ++i)
        for (std::size_t j = 0; j < 12; ++j) {
            const Real u = rule.nodes[i], v = rule.nodes[j];
            const std::array<Real, 3> beta{(1 - u) * (1 - v), u, (1 - u) * v};
            ComplexVector point{};
            for (std::size_t axis = 0; axis < 3; ++axis)
                for (std::size_t vertex = 0; vertex < 3; ++vertex)
                    point[axis] += beta[vertex] * x[vertex][axis];
            const auto r = subtract(p, point);
            const Complex r2 = dot(r, r), rn = dot(r, normal), epsilon2 = values[12] * values[12];
            const Complex den = r2 + epsilon2;
            const Complex factor = jac * rule.weights[i] * rule.weights[j] * (1 - u) /
                                   (4 * pi * den * den * sqrt(den));
            for (std::size_t k = 0; k < 3; ++k)
                for (std::size_t axis = 0; axis < 3; ++axis)
                    output[3 * k + axis] +=
                        beta[k] * ((r2 - 2 * epsilon2) * normal[axis] - 3 * r[axis] * rn) * factor;
        }
    return output;
}
} // namespace
namespace {
void verify_area_fixture(const char* path, std::size_t expected_cases,
                         bool check_reference_center = false) {
    std::ifstream source{path};
    REQUIRE(source.good());
    nlohmann::json fixture;
    source >> fixture;
    REQUIRE(fixture.at("cases").size() == expected_cases);
    core::ResourcePlan plan{256ULL * 1024 * 1024};
    auto kernel = P1Workspace::create(plan);
    REQUIRE(kernel);
    for (const auto& input : fixture.at("cases")) {
        INFO(input.at("id").get<std::string>());
        TriangleFixture geometry;
        const auto vertices = input.at("vertices").get<std::array<std::array<double, 3>, 3>>();
        for (std::size_t vertex = 0; vertex < 3; ++vertex)
            for (std::size_t axis = 0; axis < 3; ++axis)
                geometry.reference[3 * vertex + axis] = vertices[vertex][axis];
        geometry.current = geometry.reference;
        geometry.radius = SingleWord{input.at("radius").get<double>()};
        auto local = make_packet(geometry, definition(input.at("query"), geometry), plan);
        auto actual = take(evaluate_p1_basis(
            local, {core::Precision::TwoWord, 1e-10, 1e-10},
            {numerics::linear::EstimateAdmission{28672}, UINT64_C(1000000000)}, **kernel, plan));
        for (std::size_t k = 0; k < 3; ++k)
            for (std::size_t axis = 0; axis < 3; ++axis) {
                INFO("basis " << k << " axis " << axis);
                close(center(take(actual.basis(k, axis)).get()),
                      Real{input.at("basis")[k][axis].get<std::string>()});
                // This observes containment of the independently refined centre,
                // not a rigorous bound on the remaining continuous quadrature error.
                if (check_reference_center) {
                    const auto& enclosure = take(actual.basis(k, axis)).get();
                    CHECK(abs(center(enclosure) -
                              Real{input.at("basis")[k][axis].get<std::string>()}) <=
                          radius(enclosure));
                }
                for (std::size_t j = 0; j < 13; ++j) {
                    INFO("variable " << j);
                    close(center(take(actual.derivative(k, axis, j)).get()),
                          Real{input.at("jacobian")[k][axis][j].get<std::string>()});
                    if (check_reference_center) {
                        const auto& enclosure = take(actual.derivative(k, axis, j)).get();
                        CHECK(abs(center(enclosure) -
                                  Real{input.at("jacobian")[k][axis][j].get<std::string>()}) <=
                              radius(enclosure));
                    }
                }
            }
    }
}

} // namespace
TEST_CASE("P1 all 126 components agree with independent area and ambient complex step",
          "[UI-D01][UI-D02][UI-D03][UI-D04][UI-D05][UI-D06][N4-06][verification]") {
    verify_area_fixture(UVLAS_P1_AREA_FIXTURE, 19);
}
TEST_CASE("P1 refined near and self area references cover all 126 components",
          "[UI-D02][UI-D03][UI-D04][UI-D06][N4-06][verification]") {
    verify_area_fixture(UVLAS_P1_AREA_PRECISE_FIXTURE, 2, true);
}
TEST_CASE("P1 far roundoff encloses independent high precision fixed Gauss12 jets",
          "[UI-D02][UI-D04][UI-D06][N4-06][verification]") {
    core::ResourcePlan plan{128ULL * 1024 * 1024};
    auto kernel = P1Workspace::create(plan);
    REQUIRE(kernel);
    TriangleFixture geometry;
    geometry.reference = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    geometry.current = geometry.reference;
    geometry.radius = SingleWord{.08};
    auto local = make_packet(geometry, FreePoint{FreeSingle{{5, 4, 6}}}, plan);
    auto actual = take(evaluate_p1_basis(
        local, {core::Precision::TwoWord, 1e-10, 1e-10},
        {numerics::linear::EstimateAdmission{28672}, UINT64_C(1000000000)}, **kernel, plan));
    const auto rule = reference_rule();
    const std::array<Complex, 13> values{
        Complex{5}, Complex{4}, Complex{6}, Complex{0}, Complex{0}, Complex{0},  Complex{1},
        Complex{0}, Complex{0}, Complex{0}, Complex{1}, Complex{0}, Complex{.08}};
    const auto basis = fixed_area(values, rule);
    for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const auto& enclosure = take(actual.basis(k, axis)).get();
            CHECK(abs(center(enclosure) - real(basis[3 * k + axis])) <= radius(enclosure));
        }
    const Real step{"1e-45"};
    for (std::size_t j = 0; j < 13; ++j) {
        auto perturbed = values;
        perturbed[j] += Complex{Real{0}, step};
        const auto differentiated = fixed_area(perturbed, rule);
        for (std::size_t k = 0; k < 3; ++k)
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const auto& enclosure = take(actual.derivative(k, axis, j)).get();
                CHECK(abs(center(enclosure) - imag(differentiated[3 * k + axis]) / step) <=
                      radius(enclosure));
            }
    }
}
