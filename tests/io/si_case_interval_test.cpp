#include "si_program.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/rational_adaptor.hpp>
#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <limits>

namespace {
using namespace uvlas::io;
namespace si = uvlas::io::detail::si;
using Big = boost::multiprecision::cpp_int;
using Rational = boost::multiprecision::cpp_rational;
Rational exact(const si::Bound& x) {
    Rational value{x.mantissa.convert_to<Big>()};
    if (x.power >= 0)
        value *= Big{1} << static_cast<unsigned>(x.power);
    else
        value /= Big{1} << static_cast<unsigned>(-x.power);
    return value;
}
void contains(const si::Interval& x, const Rational& value) {
    CHECK(exact(x.lower) <= value);
    CHECK(exact(x.upper) >= value);
}
} // namespace
TEST_CASE("SI directed arithmetic encloses independent rational operations", "[IO4-01][property]") {
    si::Work work{0, 1000000, {}};
    si::ProgramScratch scratch;
    const auto third = si::divide(si::integer(1), si::integer(3), work);
    contains(third, Rational{1} / 3);
    const auto negative = si::divide(si::integer(-7), si::integer(11), work);
    contains(negative, Rational{-7} / 11);
    contains(si::add(third, negative, work), Rational{1} / 3 - Rational{7} / 11);
    contains(si::multiply(third, negative, work), Rational{-7} / 33);
    contains(si::divide(third, negative, work), Rational{-11} / 21);
    contains(si::programmed_power(si::integer(10), -400, scratch, work),
             Rational{1} / boost::multiprecision::pow(Big{10}, 400));
    contains(si::programmed_decimal({1, "1000000000000000000000000000001", -30}, scratch, work),
             Rational{boost::multiprecision::pow(Big{10}, 30) + 1} /
                 boost::multiprecision::pow(Big{10}, 30));
}
TEST_CASE("SI Machin pi encloses independent rational alternating sums", "[IO4-01][verification]") {
    si::Work work{0, 1000000, {}};
    si::ProgramScratch scratch;
    const auto atan = [](unsigned denominator) {
        Rational sum = 0;
        for (unsigned j = 0; j < 160; ++j) {
            const Rational term =
                Rational{1} /
                ((2 * j + 1) * boost::multiprecision::pow(Big{denominator}, 2 * j + 1));
            sum += (j & 1) ? -term : term;
        }
        return sum;
    };
    const Rational reference = 16 * atan(5) - 4 * atan(239);
    const auto pi = si::programmed_pi(scratch, work);
    contains(pi, reference);
    CHECK(exact(pi.upper) - exact(pi.lower) < Rational{1} / (Big{1} << 240));
    const auto output = si::center(pi, true, work);
    CHECK(exact(output.error) >= 0);
    CHECK(abs(reference - exact(output.center)) <= exact(output.error));
}

TEST_CASE("SI word encoding preserves ties subnormals and normalized carry", "[IO4-01][property]") {
    si::Work work{0, 1000000, {}};
    const si::Bound tie{(si::Integer{1} << 53U) + 1, -53};
    CHECK(std::bit_cast<std::uint64_t>(si::word(tie, false, true, work).word) ==
          0x3ff0000000000000ULL);
    CHECK(std::bit_cast<std::uint64_t>(si::word(tie, true, false, work).word) ==
          0x3ff0000000000001ULL);
    const si::Bound half_subnormal{1, -1075};
    CHECK(std::bit_cast<std::uint64_t>(si::word(half_subnormal, false, true, work).word) == 0);
    CHECK(std::bit_cast<std::uint64_t>(si::word(half_subnormal, true, false, work).word) == 1);
    CHECK(std::bit_cast<std::uint64_t>(
              si::word(si::negate(half_subnormal), false, false, work).word) ==
          0x8000000000000001ULL);
    const si::Bound near_two{(si::Integer{1} << 54U) - 1, -53};
    const auto binary = si::center({near_two, near_two}, false, work);
    const auto one = std::get<ScaledBinary64>(binary.value);
    CHECK(one.center == 1.0);
    CHECK(one.exponent2 == 1);
    CHECK(exact(binary.error) == Rational{1} / (Big{1} << 53U));
    const auto paired = si::center({near_two, near_two}, true, work);
    CHECK(exact(paired.center) == exact(near_two));
    CHECK(exact(paired.error) == 0);
}

TEST_CASE("SI programmed work counts match execution and reject before excess",
          "[IO4-01][negative]") {
    si::ProgramScratch scratch;
    for (const std::int64_t exponent : {0, 1, -1, 31, -31, 400, -400}) {
        const auto expected = si::power_work_bound(si::magnitude(exponent)) +
                              static_cast<std::uint64_t>(exponent < 0);
        si::Work exact_limit{0, expected, {}};
        (void)si::programmed_power(si::integer(10), exponent, scratch, exact_limit);
        CHECK(exact_limit.used == expected);
        if (expected != 0) {
            si::Work insufficient{0, expected - 1, {}};
            CHECK_THROWS_AS(si::programmed_power(si::integer(10), exponent, scratch, insufficient),
                            si::ArithmeticFailure);
            CHECK(insufficient.used <= insufficient.limit);
        }
    }
    si::Work pi_limit{0, si::pi_work, {}};
    (void)si::programmed_pi(scratch, pi_limit);
    CHECK(pi_limit.used == si::pi_work);
    si::Work pi_short{0, si::pi_work - 1, {}};
    CHECK_THROWS_AS(si::programmed_pi(scratch, pi_short), si::ArithmeticFailure);
    CHECK(pi_short.used <= pi_short.limit);
    constexpr auto minimum = std::numeric_limits<std::int64_t>::min();
    constexpr auto maximum = std::numeric_limits<std::int64_t>::max();
    CHECK(si::magnitude(minimum) == (std::uint64_t{1} << 63U));
    CHECK(si::subtract_exponent(minimum, minimum) == 0);
    CHECK_THROWS_AS(si::add_exponent(maximum, 1), si::ArithmeticFailure);
    CHECK_THROWS_AS(si::subtract_exponent(0, minimum), si::ArithmeticFailure);
}
