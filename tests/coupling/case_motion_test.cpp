#include "coupling/case_motion_rotation_private.hpp"

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstddef>
#include <cstdint>

namespace {
using namespace uvlas;
namespace motion = coupling::detail::motion;
namespace m = motion::m;
using Dec = boost::multiprecision::cpp_dec_float_100;

motion::ScalarJet jet(std::int64_t value, std::int64_t first, std::int64_t second) {
    return {m::integer(value), m::integer(first), m::integer(second)};
}

motion::Matrix identity() {
    motion::Matrix result{};
    for (const auto index : {0U, 4U, 8U})
        result[index] = m::integer(1);
    return result;
}

void verify_exact(const m::Interval& actual, std::int64_t expected) {
    const auto exact = m::integer(expected);
    CHECK(m::compare(actual.lower, exact.lower) <= 0);
    CHECK(m::compare(actual.upper, exact.upper) >= 0);
    const auto tolerance = m::scale(m::integer(1), -180);
    m::Work work{0, 4, {}};
    CHECK(m::compare(actual.lower, m::subtract(exact, tolerance, work).lower) >= 0);
    CHECK(m::compare(actual.upper, m::add(exact, tolerance, work).upper) <= 0);
}

void verify_exact(const motion::ScalarJet& actual, const std::array<std::int64_t, 3>& expected) {
    verify_exact(actual.value, expected[0]);
    verify_exact(actual.first, expected[1]);
    verify_exact(actual.second, expected[2]);
}

using DecimalMatrix = std::array<Dec, 9>;
using DecimalJet = std::array<DecimalMatrix, 3>;

// Coefficients of an independent Taylor expansion in time, through t².
// Multiply each coefficient as a whole matrix before expanding exp(K).
DecimalJet product(const DecimalJet& a, const DecimalJet& b) {
    DecimalJet result{};
    for (std::size_t order = 0; order < result.size(); ++order)
        for (std::size_t left = 0; left <= order; ++left)
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t j = 0; j < 3; ++j)
                    for (std::size_t k = 0; k < 3; ++k)
                        result[order][3 * i + j] += a[left][3 * i + k] * b[order - left][3 * k + j];
    return result;
}

DecimalJet exponential(const DecimalJet& generator) {
    DecimalJet term{}, total{};
    for (const auto index : {0U, 4U, 8U})
        term[0][index] = total[0][index] = 1;
    for (std::size_t n = 1; n <= 144; ++n) {
        term = product(term, generator);
        for (std::size_t order = 0; order < term.size(); ++order)
            for (std::size_t entry = 0; entry < term[order].size(); ++entry) {
                term[order][entry] /= n;
                total[order][entry] += term[order][entry];
            }
    }
    return total;
}

DecimalMatrix skew(const std::array<Dec, 3>& entries) {
    const auto& [x, y, z] = entries;
    return {0, -z, y, z, 0, -x, -y, x, 0};
}

Dec decimal(const m::Bound& bound) {
    return ldexp(bound.mantissa.convert_to<Dec>(), bound.power);
}

void verify_oracle(const m::Interval& actual, const Dec& expected) {
    const Dec tolerance{"1e-70"}, oracle_rounding{"1e-90"};
    CHECK(decimal(actual.lower) <= expected + oracle_rounding);
    CHECK(decimal(actual.upper) >= expected - oracle_rounding);
    CHECK(decimal(actual.lower) >= expected - tolerance);
    CHECK(decimal(actual.upper) <= expected + tolerance);
}

} // namespace

TEST_CASE("Prescribed motion composes local rotation and translation with the world frame",
          "[IO4-01][N4-03][verification]") {
    const auto zero = m::integer(0);
    const auto one = m::integer(1);
    const motion::Matrix frame{one, zero, zero, zero, zero, m::integer(-1), zero, one, zero};
    const motion::Vector centre{m::integer(3), m::integer(4), m::integer(5)};
    const motion::VectorJet offset{{{}, jet(5, 7, 11), {}}};
    m::Work work{0, 20'000, {}};
    const auto pi = m::pi_constant(work);
    motion::VectorJet phi{{{}, {}, {m::scale(pi, -1), m::integer(2), m::integer(3)}}};
    const auto result = motion::compose_local(centre, frame, offset, phi, pi, work);
    verify_exact(result.centre[0], {3, 0, 0});
    verify_exact(result.centre[1], {4, 0, 0});
    verify_exact(result.centre[2], {10, 7, 11});
    verify_exact(result.rotation[0], {0, -2, -3});
    verify_exact(result.rotation[3], {0, 0, 0});
    verify_exact(result.rotation[6], {1, 0, -4});

    phi[2].value = pi;
    const auto half_turn = motion::compose_local(centre, frame, offset, phi, pi, work);
    verify_exact(half_turn.rotation[0], {-1, 0, 4});
    verify_exact(half_turn.rotation[3], {0, 0, 0});
    verify_exact(half_turn.rotation[6], {0, -2, -3});
}

TEST_CASE("Prescribed rotation retains derivatives at zero and under changing axes",
          "[IO4-01][N4-03][verification][negative]") {
    m::Work work{0, 20'000, {}};
    const auto pi = m::pi_constant(work);
    const motion::VectorJet zero_point{jet(0, 1, 2), jet(0, 2, -3), jet(0, 3, 4)};
    const auto at_zero = motion::compose_local({}, identity(), {}, zero_point, pi, work);
    verify_exact(at_zero.rotation[0], {1, 0, -13});
    verify_exact(at_zero.rotation[1], {0, -3, -2});
    verify_exact(at_zero.rotation[5], {0, -1, 4});

    const std::array<std::int64_t, 3> value{2, -3, 4}, first{1, 3, -2}, second{-2, 4, 1};
    for (const int power : {-3, 0}) {
        CAPTURE(power);
        motion::VectorJet phi{};
        std::array<Dec, 3> a{}, b{}, c{};
        for (std::size_t i = 0; i < phi.size(); ++i) {
            phi[i] = {m::scale(m::integer(value[i]), power), m::scale(m::integer(first[i]), -3),
                      m::scale(m::integer(second[i]), -3)};
            a[i] = ldexp(Dec{value[i]}, power);
            b[i] = Dec{first[i]} / 8;
            c[i] = Dec{second[i]} / 16;
        }
        const auto expected = exponential({skew(a), skew(b), skew(c)});
        const auto actual = motion::compose_local({}, identity(), {}, phi, pi, work);
        for (std::size_t entry = 0; entry < actual.rotation.size(); ++entry) {
            CAPTURE(entry);
            verify_oracle(actual.rotation[entry].value, expected[0][entry]);
            verify_oracle(actual.rotation[entry].first, expected[1][entry]);
            verify_oracle(actual.rotation[entry].second, 2 * expected[2][entry]);
        }
    }

    m::Work denied{0, 1, {}};
    REQUIRE_THROWS_AS(motion::compose_local({}, identity(), {}, zero_point, pi, denied),
                      m::ArithmeticFailure);
}
