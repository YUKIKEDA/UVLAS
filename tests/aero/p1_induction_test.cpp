#include "aero/local_geometry_bounds_private.hpp"
#include "aero/p1_exact_product_private.hpp"
#include "aero/p1_gauss_private.hpp"
#include "aero/p1_induction.hpp"
#include "aero/p1_interval_private.hpp"
#include "aero/p1_quadrature_private.hpp"
#include "aero/p1_qualification_private.hpp"
#include "geometry_snapshot_fixture.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <nlohmann/json.hpp>

#include <fstream>
#include <thread>

namespace {
using namespace uvlas;
using namespace uvlas::test_geometry;
using namespace subject;
namespace arithmetic = uvlas::aero::detail::p1;
using Real = boost::multiprecision::cpp_bin_float_100;
Real scalar(const arithmetic::Bound& bound) {
    return boost::multiprecision::ldexp(bound.mantissa.convert_to<Real>(), bound.power);
}
boost::multiprecision::cpp_rational exact_bound(const arithmetic::Bound& b) {
    using Rational = boost::multiprecision::cpp_rational;
    const Rational value = b.mantissa.convert_to<Rational>();
    return b.power >= 0
               ? Rational{value * Rational{boost::multiprecision::cpp_int{1} << b.power}}
               : Rational{value / Rational{boost::multiprecision::cpp_int{1} << (-b.power)}};
}
void contains(const arithmetic::Interval& interval, const Real& value) {
    CHECK(scalar(interval.lower) <= value);
    CHECK(scalar(interval.upper) >= value);
}
double center(const ScaledEnclosure& value) {
    double result = 0;
    const auto words = value.center();
    if (const auto* one = std::get_if<SingleWord>(&words))
        result = one->value;
    else {
        const auto& pair = std::get<DoubleWord>(words);
        result = pair.high + pair.low;
    }
    return std::scalbn(result, value.exponent());
}
LocalGeometryPacket packet(const TriangleFixture& fixture, QueryDefinition query,
                           core::ResourcePlan& plan) {
    auto geometry_workspace = workspace(plan);
    auto packet_workspace = LocalPacketWorkspace::create(plan);
    REQUIRE(packet_workspace);
    auto geometry =
        take(prepare_owned_geometry(fixture.input(), limits(), *geometry_workspace, plan));
    auto resolved = take(resolve_query(geometry, {QueryId{44}, query}, geometry.context(), limits(),
                                       *geometry_workspace, plan));
    return take(prepare_local_packet(geometry, resolved, FacetId{70},
                                     {core::Precision::TwoWord, 1e-10}, limits(),
                                     **packet_workspace, plan));
}
P1Limits p1_limits() {
    return {numerics::linear::EstimateAdmission{28672}, UINT64_C(1000000000)};
}
} // namespace
TEST_CASE("P1 scalar primitives enclose independent higher precision values",
          "[UI-D06][N4-06][verification]") {
    arithmetic::Work work{0, UINT64_C(10000000), {}};
    const auto one = arithmetic::integer(1), two = arithmetic::integer(2);
    const auto pi = arithmetic::pi_constant(work);
    const auto ln2 =
        arithmetic::logarithm_series(arithmetic::divide(one, arithmetic::integer(3), work), work);
    contains(pi, boost::multiprecision::acos(Real{-1}));
    contains(ln2, boost::multiprecision::log(Real{2}));
    for (int i = 1; i <= 7; ++i) {
        const auto value = arithmetic::divide(arithmetic::integer(i), arithmetic::integer(3), work);
        contains(arithmetic::square_root(value, work), boost::multiprecision::sqrt(Real{i} / 3));
        contains(arithmetic::log_one_plus(value, ln2, work),
                 boost::multiprecision::log(1 + Real{i} / 3));
        contains(arithmetic::atan_two(value, two, pi, work),
                 boost::multiprecision::atan2(Real{i} / 3, Real{2}));
        contains(arithmetic::atan_two(value, arithmetic::integer(-2), pi, work),
                 boost::multiprecision::atan2(Real{i} / 3, Real{-2}));
    }
    CHECK_THROWS_AS(arithmetic::divide(one, arithmetic::integer(0), work),
                    arithmetic::ArithmeticFailure);
}
TEST_CASE("Product Gauss brackets match independently isolated exact Sturm roots",
          "[UI-D02][N4-06][verification]") {
    std::ifstream input{UVLAS_P1_GAUSS_FIXTURE};
    REQUIRE(input.good());
    nlohmann::json reference;
    input >> reference;
    const auto check = [&]<std::size_t N>(const std::array<arithmetic::GaussRoot, N>& table) {
        const auto& entries = reference.at("tables").at(N == 8 ? 0 : 1).at("entries");
        arithmetic::Work work{0, UINT64_C(1000000), {}};
        arithmetic::QuadratureRule product;
        arithmetic::build_rule(table, product, work);
        REQUIRE(entries.size() == N);
        for (std::size_t i = 0; i < N; ++i) {
            CHECK(table[i].lower == entries[i].at("node").at("lo").template get<std::string>());
            CHECK(table[i].upper == entries[i].at("node").at("hi").template get<std::string>());
            const auto& counts = entries[i].at("sturm_variations");
            CHECK(counts[0].template get<int>() - counts[1].template get<int>() == 1);
            const arithmetic::Integer low{
                entries[i].at("weight").at("lo").template get<std::string>()};
            const arithmetic::Integer high{
                entries[i].at("weight").at("hi").template get<std::string>()};
            CHECK(arithmetic::compare(product.weights[i].lower, {low, -256}) <= 0);
            CHECK(arithmetic::compare(product.weights[i].upper, {high, -256}) >= 0);
            CHECK(arithmetic::compare(product.weights[i].lower, {}) > 0);
        }
        std::array<arithmetic::Interval, N> powers;
        powers.fill(arithmetic::integer(1));
        for (unsigned degree = 0; degree < 2 * N; ++degree) {
            auto sum = arithmetic::integer(0);
            for (std::size_t i = 0; i < N; ++i) {
                sum = arithmetic::add(
                    sum, arithmetic::multiply(powers[i], product.weights[i], work), work);
                powers[i] = arithmetic::multiply(powers[i], product.nodes[i], work);
            }
            const boost::multiprecision::cpp_rational expected =
                boost::multiprecision::cpp_rational{1} / (degree + 1);
            CHECK(exact_bound(sum.lower) <= expected);
            CHECK(exact_bound(sum.upper) >= expected);
        }
    };
    check(arithmetic::gauss_8);
    check(arithmetic::gauss_12);
}
TEST_CASE("P1 primitive dyadic boundaries retain subgrid logarithms and quadrant limits",
          "[UI-D06][N4-06][verification]") {
    arithmetic::Work work{0, UINT64_C(10000000), {}};
    const auto zero = arithmetic::integer(0), one = arithmetic::integer(1);
    const auto pi = arithmetic::pi_constant(work);
    const auto ln2 =
        arithmetic::logarithm_series(arithmetic::divide(one, arithmetic::integer(3), work), work);
    contains(arithmetic::square_root(zero, work), Real{0});
    contains(arithmetic::square_root(arithmetic::integer(4), work), Real{2});
    for (const int power : {-9000, -1074, -256, 0, 256, 9000}) {
        const auto dyadic = arithmetic::scale(one, power);
        const Real reference = boost::multiprecision::ldexp(Real{1}, power);
        contains(arithmetic::square_root(dyadic, work), sqrt(reference));
        contains(arithmetic::log_one_plus(dyadic, ln2, work),
                 boost::multiprecision::log1p(reference));
        contains(arithmetic::atan_two(one, dyadic, pi, work),
                 boost::multiprecision::atan2(Real{1}, reference));
        contains(arithmetic::atan_two(one, arithmetic::negate(dyadic, work), pi, work),
                 boost::multiprecision::atan2(Real{1}, -reference));
    }
    contains(arithmetic::atan_two(one, zero, pi, work), boost::multiprecision::acos(Real{-1}) / 2);
    for (const double x : {std::nextafter(1., 0.), 1., std::nextafter(1., 2.),
                           std::nextafter(3., 0.), 3., std::nextafter(3., 4.)}) {
        ExactInteger decoded;
        local_decode(x, decoded);
        contains(arithmetic::log_one_plus(arithmetic::exact(decoded, -1074), ln2, work),
                 boost::multiprecision::log1p(Real{x}));
    }
}
TEST_CASE("P1 complete basis and derivative owner evaluates all branches",
          "[UI-D01][UI-D02][UI-D03][UI-D04][UI-D06][N4-02][verification]") {
    core::ResourcePlan plan{128ULL * 1024 * 1024};
    auto kernel = P1Workspace::create(plan);
    REQUIRE(kernel);
    TriangleFixture fixture;
    for (const auto query : {FreeSingle{{.75, 1, 2}}, FreeSingle{{.75, 1, 30}}}) {
        auto local = packet(fixture, FreePoint{query}, plan);
        auto result = take(evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10},
                                             p1_limits(), **kernel, plan));
        CHECK(take(result.steps()) <= take(result.work_plan()).maximum_steps);
        for (std::size_t k = 0; k < 3; ++k)
            for (std::size_t axis = 0; axis < 3; ++axis) {
                CHECK(std::isfinite(center(take(result.basis(k, axis)).get())));
                for (std::size_t j = 0; j < 13; ++j)
                    CHECK(std::isfinite(center(take(result.derivative(k, axis, j)).get())));
            }
        CHECK(take(result.snapshot_identity()) == take(take(result.geometry()).get().identity()));
    }
    fixture.radius = SingleWord{0};
    auto outside = packet(fixture, FreePoint{FreeSingle{{6, 6, 0}}}, plan);
    auto result = take(evaluate_p1_basis(outside, {core::Precision::TwoWord, 1e-10, 1e-10},
                                         p1_limits(), **kernel, plan));
    CHECK(take(result.branch()) == P1Branch::IdealCoplanarExterior);
}
namespace {
std::array<std::uint64_t, 504> owned_bits(const P1Basis& result) {
    std::array<std::uint64_t, 504> bits{};
    for (std::size_t component = 0; component < 9; ++component)
        for (std::size_t jet = 0; jet < 14; ++jet) {
            const auto& item =
                take(jet == 0 ? result.basis(component / 3, component % 3)
                              : result.derivative(component / 3, component % 3, jet - 1))
                    .get();
            const auto words = item.center();
            double high = 0, low = 0;
            if (const auto* single = std::get_if<SingleWord>(&words))
                high = single->value;
            else {
                const auto& pair = std::get<DoubleWord>(words);
                high = pair.high;
                low = pair.low;
            }
            const auto index = 4 * (14 * component + jet);
            bits[index] = std::bit_cast<std::uint64_t>(high);
            bits[index + 1] = std::bit_cast<std::uint64_t>(low);
            bits[index + 2] = std::bit_cast<std::uint64_t>(item.radius());
            bits[index + 3] =
                static_cast<std::uint64_t>(static_cast<std::int64_t>(item.exponent()));
        }
    return bits;
}
Real exact_center(const ScaledEnclosure& value) {
    Real result = 0;
    const auto words = value.center();
    if (const auto* single = std::get_if<SingleWord>(&words))
        result = Real{single->value};
    else {
        const auto& pair = std::get<DoubleWord>(words);
        result = Real{pair.high} + Real{pair.low};
    }
    return boost::multiprecision::ldexp(result, value.exponent());
}
} // namespace
TEST_CASE("P1 owns original geometry and survives failures moves and workspace reuse",
          "[UI-D04][UI-D06][N4-02][N4-03][N4-08][negative][verification]") {
    core::ResourcePlan setup{128ULL * 1024 * 1024};
    auto kernel = P1Workspace::create(setup);
    REQUIRE(kernel);
    const TriangleFixture fixture;
    auto local =
        packet(fixture, SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}, setup);
    const auto identity = take(local.snapshot_identity());
    auto result = take(evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10},
                                         p1_limits(), **kernel, setup));
    const auto saved = owned_bits(result);
    const auto admitted = take(result.work_plan());
    {
        core::ResourcePlan exact{P1Basis::payload_bytes()};
        auto value = evaluate_p1_basis(
            local, {core::Precision::TwoWord, 1e-10, 1e-10},
            {numerics::linear::EstimateAdmission{28672}, admitted.maximum_steps}, **kernel, exact);
        REQUIRE(value);
        CHECK(owned_bits(*value) == saved);
        CHECK(exact.allocated_bytes() == P1Basis::payload_bytes());
    }
    core::ResourcePlan short_result{P1Basis::payload_bytes() - 1};
    failure(evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10}, p1_limits(),
                              **kernel, short_result),
            core::FailureCode::ResourceExceeded);
    auto short_work = p1_limits();
    short_work.max_steps = admitted.maximum_steps - 1;
    failure(evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10}, short_work, **kernel,
                              setup),
            core::FailureCode::ResourceExceeded);
    failure(
        evaluate_p1_basis(local, {core::Precision::TwoWord, 0, 0}, p1_limits(), **kernel, setup),
        core::FailureCode::PrecisionUnresolved);
    failure(evaluate_p1_basis(local, {core::Precision::TwoWord, -1, 1e-10}, p1_limits(), **kernel,
                              setup),
            core::FailureCode::InvalidInput);
    failure(evaluate_p1_basis(
                local, {core::Precision::TwoWord, std::numeric_limits<double>::infinity(), 1e-10},
                p1_limits(), **kernel, setup),
            core::FailureCode::Nonfinite);
    failure(evaluate_p1_basis(local, {static_cast<core::Precision>(255), 1e-10, 1e-10}, p1_limits(),
                              **kernel, setup),
            core::FailureCode::PrecisionUnsupported);
    auto strict = p1_limits();
    strict.external_memory = numerics::linear::RequireProvenBound{28672};
    failure(
        evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10}, strict, **kernel, setup),
        core::FailureCode::Unsupported);
    // MSVC request_stop() is non-const; libstdc++ declares it const.
    // NOLINTNEXTLINE(misc-const-correctness)
    std::stop_source cancellation;
    cancellation.request_stop();
    failure(evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10}, p1_limits(),
                              **kernel, setup, cancellation.get_token()),
            core::FailureCode::Cancelled);
    CHECK(owned_bits(result) == saved);
    CHECK(take(local.snapshot_identity()) == identity);
    auto moved = std::move(result);
    // NOLINTNEXTLINE(bugprone-use-after-move): verify documented moved-from rejection.
    CHECK_FALSE(result.active());
    // NOLINTNEXTLINE(bugprone-use-after-move): verify documented moved-from rejection.
    failure(result.basis(0, 0), core::FailureCode::InvalidInput);
    failure(moved.basis(3, 0), core::FailureCode::InvalidInput);
    failure(moved.derivative(0, 0, 13), core::FailureCode::InvalidInput);
    CHECK(take(take(moved.geometry()).get().identity()) == identity);
    auto local_owner = std::move(local);
    // NOLINTNEXTLINE(bugprone-use-after-move): invalid moved-from input must fail by value.
    failure(evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10}, p1_limits(),
                              **kernel, setup),
            core::FailureCode::InvalidInput);
    CHECK(local_owner.active());
    CHECK(owned_bits(moved) == saved);
}
TEST_CASE("P1 huge common origins preserve low words and scaled derivatives",
          "[UI-D04][UI-D06][N4-06][verification]") {
    core::ResourcePlan plan{256ULL * 1024 * 1024};
    auto kernel = P1Workspace::create(plan);
    REQUIRE(kernel);
    const TriangleFixture flat;
    auto basic = packet(flat, SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}, plan);
    auto baseline = take(evaluate_p1_basis(basic, {core::Precision::TwoWord, 1e-10, 1e-10},
                                           p1_limits(), **kernel, plan));
    TriangleFixture translated;
    translated.two_word = true;
    translated.reference_low = translated.reference;
    translated.current_low = translated.current;
    translated.reference.fill(std::scalbn(1., 500));
    translated.current = translated.reference;
    auto shifted =
        packet(translated, SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}, plan);
    auto actual = take(evaluate_p1_basis(shifted, {core::Precision::TwoWord, 1e-10, 1e-10},
                                         p1_limits(), **kernel, plan));
    CHECK(owned_bits(actual) == owned_bits(baseline));
    for (const int scale : {-900, 900}) {
        TriangleFixture scaled;
        for (auto& value : scaled.reference)
            value = std::scalbn(value, scale);
        scaled.current = scaled.reference;
        scaled.radius = SingleWord{std::scalbn(.25, scale)};
        auto local =
            packet(scaled, SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}, plan);
        auto value = take(evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10},
                                            p1_limits(), **kernel, plan));
        for (std::size_t k = 0; k < 3; ++k)
            for (std::size_t a = 0; a < 3; ++a) {
                const Real normalized = boost::multiprecision::ldexp(
                    exact_center(take(value.basis(k, a)).get()), scale);
                const Real original = exact_center(take(baseline.basis(k, a)).get());
                CHECK(abs(normalized - original) < Real{"1e-28"});
                for (std::size_t j = 0; j < 13; ++j) {
                    const Real derivative = boost::multiprecision::ldexp(
                        exact_center(take(value.derivative(k, a, j)).get()), 2 * scale);
                    CHECK(abs(derivative - exact_center(take(baseline.derivative(k, a, j)).get())) <
                          Real{"1e-28"});
                }
            }
    }
    TriangleFixture tiny;
    tiny.radius = SingleWord{std::scalbn(1., -600)};
    auto low_radius = packet(tiny, FreePoint{FreeSingle{{.75, 1, 2}}}, plan);
    CHECK(evaluate_p1_basis(low_radius, {core::Precision::TwoWord, 1e-10, 1e-10}, p1_limits(),
                            **kernel, plan)
              .has_value());
}
TEST_CASE("P1 independent workspaces return byte-identical owning results",
          "[UI-D06][N4-08][N4-09][verification]") {
    core::ResourcePlan setup{128ULL * 1024 * 1024};
    const TriangleFixture fixture;
    auto local = packet(fixture, FreePoint{FreeSingle{{.75, 1, 2}}}, setup);
    std::array<std::optional<std::array<std::uint64_t, 504>>, 2> values;
    const auto execute = [&](std::size_t index) {
        core::ResourcePlan plan{128ULL * 1024 * 1024};
        auto workspace = P1Workspace::create(plan);
        if (!workspace)
            return;
        auto result = evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10},
                                        p1_limits(), **workspace, plan);
        if (!result)
            return;
        values[index] = owned_bits(*result);
    };
    std::thread first{execute, 0}, second{execute, 1};
    first.join();
    second.join();
    const auto& left = values[0];
    const auto& right = values[1];
    REQUIRE(left);
    REQUIRE(right);
    if (!left || !right)
        return;
    CHECK(*left == *right);
}
TEST_CASE("P1 allocation free exact products match an unbounded independent integer oracle",
          "[UI-D06][N4-02][N4-06][verification]") {
    using Big = boost::multiprecision::cpp_int;
    arithmetic::ExactProductScratch scratch;
    const auto check = [&](const Big& a, const Big& b) {
        const auto actual = arithmetic::exact_product(a.convert_to<ExactInteger>(),
                                                      b.convert_to<ExactInteger>(), scratch);
        CHECK(actual.convert_to<Big>() == a * b);
    };
    for (const int sign_left : {-1, 1})
        for (const int sign_right : {-1, 1}) {
            check(Big{sign_left} * ((Big{1} << 3200) - 1),
                  Big{sign_right} * ((Big{1} << 6400) - 1));
            check(Big{sign_left} * ((Big{1} << 19000) - 1),
                  Big{sign_right} * ((Big{1} << 6400) - 1));
            check(Big{sign_left} * ((Big{1} << 19000) - 1),
                  Big{sign_right} * ((Big{1} << 3200) - 1));
            check(Big{sign_left} * (Big{0xffffffffU} << 4096),
                  Big{sign_right} * (Big{0xffffffffU} << 2048));
        }
    check(Big{0}, (Big{1} << 32767) - 1);
    check((Big{1} << 32767) - 1, Big{1});
    CHECK_THROWS(arithmetic::exact_product(((Big{1} << 19000) - 1).convert_to<ExactInteger>(),
                                           ((Big{1} << 19000) - 1).convert_to<ExactInteger>(),
                                           scratch));
}
TEST_CASE("P1 interval operations include low dyadic terms and zero crossing squares",
          "[UI-D06][N4-06][verification]") {
    using Rational = boost::multiprecision::cpp_rational;
    const auto rational = [](const arithmetic::Bound& b) {
        const Rational value = b.mantissa.convert_to<Rational>();
        return b.power >= 0 ? value * Rational{boost::multiprecision::cpp_int{1} << b.power}
                            : value / Rational{boost::multiprecision::cpp_int{1} << (-b.power)};
    };
    arithmetic::Work work{0, UINT64_C(100000), {}};
    for (const int a : {-17, -1, 0, 1, 19})
        for (const int b : {-7, -1, 1, 11}) {
            const auto left =
                arithmetic::divide(arithmetic::integer(a), arithmetic::integer(3), work);
            const auto right =
                arithmetic::divide(arithmetic::integer(b), arithmetic::integer(5), work);
            const auto sum = arithmetic::add(left, right, work),
                       product = arithmetic::multiply(left, right, work),
                       quotient = arithmetic::divide(left, right, work);
            const Rational x = Rational{a} / 3, y = Rational{b} / 5;
            CHECK(rational(sum.lower) <= x + y);
            CHECK(rational(sum.upper) >= x + y);
            CHECK(rational(product.lower) <= x * y);
            CHECK(rational(product.upper) >= x * y);
            CHECK(rational(quotient.lower) <= x / y);
            CHECK(rational(quotient.upper) >= x / y);
        }
    const arithmetic::Interval crossing{{arithmetic::Integer{-1}, -9000},
                                        {arithmetic::Integer{1}, -9000}};
    const auto square = arithmetic::multiply(crossing, crossing, work);
    CHECK(square.lower.mantissa == 0);
    CHECK(square.upper.power < 0);
    const auto high = arithmetic::integer(1),
               low = arithmetic::scale(arithmetic::integer(1), -9000);
    const auto sum = arithmetic::add(high, low, work);
    const Rational exact =
        Rational{1} + Rational{1} / Rational{boost::multiprecision::cpp_int{1} << 9000};
    CHECK(rational(sum.lower) <= exact);
    CHECK(rational(sum.upper) >= exact);
}

TEST_CASE("P1 exact qualification compares full width including physical derivative scale",
          "[UI-D06][N4-06][negative][verification]") {
    arithmetic::ExactProductScratch product;
    LocalRootScratch scratch;
    const ExactInteger unit = ExactInteger{1} << 4296;
    CHECK_FALSE(arithmetic::full_width_at_most(.25, 0, .25, unit, false, product, scratch));
    CHECK(arithmetic::full_width_at_most(.25, 0, .5, unit, false, product, scratch));
    CHECK_FALSE(arithmetic::full_width_at_most(.25, 0, std::nextafter(.5, 0.), unit, false, product,
                                               scratch));
    const ExactInteger four = ExactInteger{4} << 4296;
    CHECK(arithmetic::full_width_at_most(.25, 0, 1., four, false, product, scratch));
    CHECK_FALSE(arithmetic::full_width_at_most(.25, 0, 1., four, true, product, scratch));
    CHECK(arithmetic::full_width_at_most(.25, 0, 2., four, true, product, scratch));
    for (const int scale : {-900, 900}) {
        const ExactInteger length = ExactInteger{1} << static_cast<unsigned>(4296 - 2 * scale);
        CHECK(arithmetic::full_width_at_most(.25, scale, .5, length, false, product, scratch));
        CHECK_FALSE(
            arithmetic::full_width_at_most(.25, scale, .25, length, false, product, scratch));
        CHECK(arithmetic::full_width_at_most(.25, 2 * scale, .5, length, true, product, scratch));
        CHECK_FALSE(
            arithmetic::full_width_at_most(.25, 2 * scale, .25, length, true, product, scratch));
    }
    CHECK(arithmetic::full_width_at_most(0, 0, 0, unit, false, product, scratch));
    CHECK_FALSE(arithmetic::full_width_at_most(std::numeric_limits<double>::denorm_min(), -9000, 0,
                                               unit, false, product, scratch));
}
TEST_CASE("P1 far qualification rejects unresolved values in every basis and ambient column",
          "[UI-D02][UI-D04][UI-D06][negative][verification]") {
    arithmetic::FarComponents low{}, high{};
    arithmetic::Work work{0, UINT64_C(10000000), {}};
    const auto unit = arithmetic::integer(1);
    for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t jet = 0; jet < 14; ++jet) {
            const auto index = 14 * (3 * k) + jet;
            low[index] = arithmetic::integer(1);
            CHECK_THROWS_AS(arithmetic::qualify_far(low, high, unit, 0, work),
                            arithmetic::ArithmeticFailure);
            low[index] = arithmetic::integer(0);
        }
    ExactInteger native;
    local_decode(1e-10, native);
    low[0] = arithmetic::exact(native, -1074);
    CHECK_NOTHROW(arithmetic::qualify_far(low, high, unit, 0, work));
    local_decode(std::nextafter(1e-10, std::numeric_limits<double>::infinity()), native);
    low[0] = arithmetic::exact(native, -1074);
    CHECK_THROWS_AS(arithmetic::qualify_far(low, high, unit, 0, work),
                    arithmetic::ArithmeticFailure);
    low[0] = arithmetic::integer(0);
    CHECK_NOTHROW(arithmetic::qualify_far(low, high, unit, 0, work));
}
namespace {
void near_value(double actual, double expected, double tolerance = 1e-10) {
    CHECK(std::abs(actual - expected) <= tolerance * std::max(1.0, std::abs(expected)));
}
P1Basis evaluate_fixture(const TriangleFixture& fixture, QueryDefinition query, P1Workspace& kernel,
                         core::ResourcePlan& plan) {
    auto local = packet(fixture, query, plan);
    return take(evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10}, p1_limits(),
                                  kernel, plan));
}
} // namespace
TEST_CASE("P1 orientation permutations preserve all ambient columns and every support",
          "[UI-D01][UI-D02][UI-D03][UI-D04][UI-D06][property][verification]") {
    core::ResourcePlan plan{256ULL * 1024 * 1024};
    auto kernel = P1Workspace::create(plan);
    REQUIRE(kernel);
    const TriangleFixture original;
    const FreePoint free{FreeSingle{{.75, 1, 2}}};
    auto baseline = evaluate_fixture(original, free, **kernel, plan);
    std::array<std::size_t, 3> order{0, 1, 2};
    do {
        TriangleFixture permuted;
        unsigned inversions = 0;
        for (std::size_t i = 0; i < 3; ++i) {
            permuted.triangles[0].points[i] = original.points[order[i]];
            for (std::size_t j = i + 1; j < 3; ++j)
                if (order[i] > order[j])
                    ++inversions;
            const auto a = order[i], b = order[(i + 1) % 3];
            const auto edge = (a + b == 1) ? 0U : ((a + b == 3) ? 1U : 2U);
            permuted.triangles[0].edges[i] = original.edges[edge].id;
        }
        const double sign = inversions % 2 == 0 ? 1 : -1;
        auto value = evaluate_fixture(permuted, free, **kernel, plan);
        for (std::size_t k = 0; k < 3; ++k)
            for (std::size_t axis = 0; axis < 3; ++axis) {
                near_value(center(take(value.basis(k, axis)).get()),
                           sign * center(take(baseline.basis(order[k], axis)).get()));
                for (std::size_t j = 0; j < 13; ++j) {
                    const auto mapped =
                        j < 3 || j == 12 ? j : 3 + 3 * order[(j - 3) / 3] + (j - 3) % 3;
                    near_value(center(take(value.derivative(k, axis, j)).get()),
                               sign *
                                   center(take(baseline.derivative(order[k], axis, mapped)).get()));
                }
            }
        const std::array<QueryDefinition, 7> supports{
            SourceVertex{original.points[0]},
            SourceVertex{original.points[1]},
            SourceVertex{original.points[2]},
            SourceEdgePoint{EdgeId{10}, SingleWord{.25}},
            SourceEdgePoint{EdgeId{20}, SingleWord{.25}},
            SourceEdgePoint{EdgeId{30}, SingleWord{.25}},
            SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}};
        for (const auto& support : supports) {
            auto local = packet(permuted, support, plan);
            CHECK(std::holds_alternative<VerifiedIncidence>(take(local.qualification())));
            auto incident = take(evaluate_p1_basis(local, {core::Precision::TwoWord, 1e-10, 1e-10},
                                                   p1_limits(), **kernel, plan));
            for (std::size_t k = 0; k < 3; ++k)
                for (std::size_t a = 0; a < 3; ++a) {
                    CHECK(std::isfinite(center(take(incident.basis(k, a)).get())));
                    for (std::size_t j = 0; j < 13; ++j)
                        CHECK(std::isfinite(center(take(incident.derivative(k, a, j)).get())));
                }
        }
    } while (std::next_permutation(order.begin(), order.end()));
}
TEST_CASE("P1 geometric covariance and exact far boundary retain full ambient action",
          "[UI-D02][UI-D03][UI-D04][UI-D06][property][verification]") {
    core::ResourcePlan plan{256ULL * 1024 * 1024};
    auto kernel = P1Workspace::create(plan);
    REQUIRE(kernel);
    const TriangleFixture original;
    const std::array<double, 3> point{.75, 1, 2};
    auto baseline = evaluate_fixture(original, FreePoint{FreeSingle{point}}, **kernel, plan);
    for (const double scale : {1e-3, 1., 1e3}) {
        TriangleFixture transformed;
        std::array<double, 3> p{};
        // Cyclic proper rotation plus translation; no branch-dependent derivative exists.
        for (std::size_t a = 0; a < 3; ++a) {
            p[a] = scale * point[(a + 1) % 3] + 8;
            for (std::size_t v = 0; v < 3; ++v)
                transformed.current[3 * v + a] = scale * original.current[3 * v + (a + 1) % 3] + 8;
        }
        transformed.reference = transformed.current;
        transformed.radius = SingleWord{scale * .25};
        auto actual = evaluate_fixture(transformed, FreePoint{FreeSingle{p}}, **kernel, plan);
        for (std::size_t k = 0; k < 3; ++k)
            for (std::size_t a = 0; a < 3; ++a) {
                near_value(scale * center(take(actual.basis(k, a)).get()),
                           center(take(baseline.basis(k, (a + 1) % 3)).get()), 1e-8);
                for (std::size_t j = 0; j < 13; ++j) {
                    const auto mapped = j == 12 ? j : 3 * (j / 3) + (j % 3 + 1) % 3;
                    near_value(scale * scale * center(take(actual.derivative(k, a, j)).get()),
                               center(take(baseline.derivative(k, (a + 1) % 3, mapped)).get()),
                               1e-8);
                }
            }
    }
    for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t a = 0; a < 3; ++a)
            for (std::size_t direction = 0; direction < 3; ++direction) {
                double translation = center(take(baseline.derivative(k, a, direction)).get());
                for (std::size_t v = 0; v < 3; ++v)
                    translation +=
                        center(take(baseline.derivative(k, a, 3 + 3 * v + direction)).get());
                near_value(translation, 0);
            }
    const std::array heights{std::nextafter(20., 0.), 20., std::nextafter(20., 21.)};
    for (std::size_t i = 0; i < 3; ++i) {
        auto value =
            evaluate_fixture(original, FreePoint{FreeSingle{{.75, 1, heights[i]}}}, **kernel, plan);
        CHECK(take(value.branch()) == (i == 2 ? P1Branch::FarGauss12 : P1Branch::Analytic));
        // The frozen branch contributes no derivative of the dispatch predicate.
        for (std::size_t k = 0; k < 3; ++k)
            for (std::size_t a = 0; a < 3; ++a) {
                double translation = center(take(value.derivative(k, a, 2)).get());
                for (std::size_t v = 0; v < 3; ++v)
                    translation += center(take(value.derivative(k, a, 5 + 3 * v)).get());
                near_value(translation, 0);
            }
    }
}
TEST_CASE("P1 constant strength recovers an independently integrated closed vortex ring",
          "[UI-D02][UI-D04][UI-D06][verification]") {
    core::ResourcePlan plan{128ULL * 1024 * 1024};
    auto kernel = P1Workspace::create(plan);
    REQUIRE(kernel);
    const TriangleFixture fixture;
    const std::array<double, 3> point{.75, 1, 2};
    auto value = evaluate_fixture(fixture, FreePoint{FreeSingle{point}}, **kernel, plan);
    std::array<double, 3> ring{};
    constexpr unsigned subdivisions = 4096;
    for (std::size_t edge = 0; edge < 3; ++edge) {
        std::array<double, 3> tangent{};
        for (std::size_t a = 0; a < 3; ++a)
            tangent[a] = fixture.current[3 * ((edge + 1) % 3) + a] - fixture.current[3 * edge + a];
        for (unsigned i = 0; i <= subdivisions; ++i) {
            const double t = static_cast<double>(i) / subdivisions;
            std::array<double, 3> r{};
            double r2 = .25 * .25;
            for (std::size_t a = 0; a < 3; ++a) {
                r[a] = point[a] - fixture.current[3 * edge + a] - t * tangent[a];
                r2 += r[a] * r[a];
            }
            const double weight = i == 0 || i == subdivisions ? 1 : (i % 2 == 0 ? 2 : 4);
            const double factor =
                -weight / (3 * subdivisions * 4 * std::acos(-1.) * r2 * std::sqrt(r2));
            for (std::size_t a = 0; a < 3; ++a)
                ring[a] += factor * (tangent[(a + 1) % 3] * r[(a + 2) % 3] -
                                     tangent[(a + 2) % 3] * r[(a + 1) % 3]);
        }
    }
    for (std::size_t a = 0; a < 3; ++a) {
        double sum = 0;
        for (std::size_t k = 0; k < 3; ++k)
            sum += center(take(value.basis(k, a)).get());
        near_value(sum, ring[a]);
    }
}
TEST_CASE("P1 material supports include every position occurrence in their chain rule",
          "[UI-D01][UI-D03][UI-D04][UI-D06][verification]") {
    core::ResourcePlan plan{512ULL * 1024 * 1024};
    auto kernel = P1Workspace::create(plan);
    REQUIRE(kernel);
    const TriangleFixture fixture;
    const std::array<QueryDefinition, 7> supports{
        SourceVertex{fixture.points[0]},
        SourceVertex{fixture.points[1]},
        SourceVertex{fixture.points[2]},
        SourceEdgePoint{EdgeId{10}, SingleWord{.25}},
        SourceEdgePoint{EdgeId{20}, SingleWord{.25}},
        SourceEdgePoint{EdgeId{30}, SingleWord{.25}},
        SourceFacePoint{FacetId{70}, SingleWord{.25}, SingleWord{.25}}};
    // Edge parameters use ascending PointId order, independently of facet orientation.
    const std::array<std::array<double, 3>, 7> beta{{{1, 0, 0},
                                                     {0, 1, 0},
                                                     {0, 0, 1},
                                                     {.25, .75, 0},
                                                     {0, .75, .25},
                                                     {.25, 0, .75},
                                                     {.5, .25, .25}}};
    constexpr double h = 0x1p-20;
    for (std::size_t support = 0; support < supports.size(); ++support) {
        auto baseline = evaluate_fixture(fixture, supports[support], **kernel, plan);
        for (std::size_t vertex = 0; vertex < 3; ++vertex)
            for (std::size_t direction = 0; direction < 3; ++direction) {
                TriangleFixture plus, minus;
                plus.current[3 * vertex + direction] += h;
                minus.current[3 * vertex + direction] -= h;
                auto upper = evaluate_fixture(plus, supports[support], **kernel, plan);
                auto lower = evaluate_fixture(minus, supports[support], **kernel, plan);
                for (std::size_t k = 0; k < 3; ++k)
                    for (std::size_t a = 0; a < 3; ++a) {
                        const double chain =
                            center(
                                take(baseline.derivative(k, a, 3 + 3 * vertex + direction)).get()) +
                            beta[support][vertex] *
                                center(take(baseline.derivative(k, a, direction)).get());
                        const double difference = (center(take(upper.basis(k, a)).get()) -
                                                   center(take(lower.basis(k, a)).get())) /
                                                  (2 * h);
                        near_value(chain, difference, 1e-8);
                    }
            }
    }
}
