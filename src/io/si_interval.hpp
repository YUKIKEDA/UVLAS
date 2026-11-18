#pragma once

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <uvlas/io/si_case.hpp>

namespace uvlas::io::detail::si {
using Integer = boost::multiprecision::number<
    boost::multiprecision::cpp_int_backend<1024, 1024, boost::multiprecision::signed_magnitude,
                                           boost::multiprecision::checked, void>,
    boost::multiprecision::et_off>;
struct ArithmeticFailure {
    core::FailureCode code;
};
struct Work {
    std::uint64_t used = 0, limit = 0;
    std::stop_token stop;
    void step(std::uint64_t amount = 1) {
        if (stop.stop_requested())
            throw ArithmeticFailure{core::FailureCode::Cancelled};
        if (amount > limit - used)
            throw ArithmeticFailure{core::FailureCode::ResourceExceeded};
        used += amount;
    }
};
inline std::int64_t add_exponent(std::int64_t a, std::int64_t b) {
    if ((b > 0 && a > std::numeric_limits<std::int64_t>::max() - b) ||
        (b < 0 && a < std::numeric_limits<std::int64_t>::min() - b))
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    return a + b;
}
inline std::int64_t subtract_exponent(std::int64_t a, std::int64_t b) {
    if ((b > 0 && a < std::numeric_limits<std::int64_t>::min() + b) ||
        (b < 0 && a > std::numeric_limits<std::int64_t>::max() + b))
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    return a - b;
}
inline std::uint64_t magnitude(std::int64_t x) noexcept {
    return x < 0 ? static_cast<std::uint64_t>(-(x + 1)) + 1 : static_cast<std::uint64_t>(x);
}
inline Integer abs_integer(const Integer& x) {
    return x < 0 ? -x : x;
}
inline unsigned bits(const Integer& x) {
    return x == 0 ? 0U : static_cast<unsigned>(boost::multiprecision::msb(abs_integer(x))) + 1U;
}
struct Bound {
    Integer mantissa = 0;
    std::int64_t power = 0;
};
struct Interval {
    Bound lower, upper;
};
inline Bound negate(Bound x) {
    x.mantissa = -x.mantissa;
    return x;
}
inline Bound rounded(Integer value, std::int64_t power, bool upper) {
    while (bits(value) > 256U) {
        const unsigned shift = bits(value) - 256U;
        const bool negative = value < 0;
        const Integer positive = abs_integer(value);
        Integer quotient = positive >> shift;
        if ((quotient << shift) != positive && upper != negative)
            ++quotient;
        value = negative ? -quotient : quotient;
        power = add_exponent(power, shift);
    }
    return value == 0 ? Bound{} : Bound{std::move(value), power};
}
inline int compare(const Bound& a, const Bound& b) {
    if ((a.mantissa < 0) != (b.mantissa < 0))
        return a.mantissa < 0 ? -1 : 1;
    if (a.mantissa == 0 || b.mantissa == 0)
        return a.mantissa == b.mantissa ? 0 : (a.mantissa < b.mantissa ? -1 : 1);
    const auto at = add_exponent(a.power, bits(a.mantissa));
    const auto bt = add_exponent(b.power, bits(b.mantissa));
    if (at != bt)
        return ((at < bt) != (a.mantissa < 0)) ? -1 : 1;
    const auto common = std::min(a.power, b.power);
    const Integer left = abs_integer(a.mantissa)
                         << static_cast<unsigned>(subtract_exponent(a.power, common));
    const Integer right = abs_integer(b.mantissa)
                          << static_cast<unsigned>(subtract_exponent(b.power, common));
    const int result = left == right ? 0 : (left < right ? -1 : 1);
    return a.mantissa < 0 ? -result : result;
}
// Rescale onto a fixed grid without ever constructing a huge shifted integer.
inline Integer at_power(const Bound& value, std::int64_t power, bool upper) {
    if (value.mantissa == 0)
        return 0;
    const bool negative = value.mantissa < 0;
    const Integer positive = abs_integer(value.mantissa);
    Integer quotient;
    bool remainder = false;
    if (value.power >= power) {
        const auto shift = subtract_exponent(value.power, power);
        if (shift > 768 - bits(positive))
            throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
        quotient = positive << static_cast<unsigned>(shift);
    } else {
        // Difference may exceed int64; compare with a bounded offset first.
        const bool all_lost = power > add_exponent(value.power, bits(positive));
        if (all_lost) {
            quotient = 0;
            remainder = true;
        } else {
            const auto shift = static_cast<unsigned>(subtract_exponent(power, value.power));
            quotient = positive >> shift;
            remainder = (quotient << shift) != positive;
        }
    }
    if (remainder && upper != negative)
        ++quotient;
    return negative ? -quotient : quotient;
}
inline Bound sum_bound(const Bound& a, const Bound& b, bool upper) {
    if (a.mantissa == 0)
        return b;
    if (b.mantissa == 0)
        return a;
    const auto top =
        std::max(add_exponent(a.power, bits(a.mantissa)), add_exponent(b.power, bits(b.mantissa)));
    const auto common = add_exponent(top, -258);
    return rounded(at_power(a, common, upper) + at_power(b, common, upper), common, upper);
}
inline Bound product_bound(const Bound& a, const Bound& b, bool upper) {
    if (a.mantissa == 0 || b.mantissa == 0)
        return {};
    return rounded(a.mantissa * b.mantissa, add_exponent(a.power, b.power), upper);
}
inline Bound quotient_bound(const Bound& a, const Bound& b, bool upper) {
    if (b.mantissa == 0)
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    if (a.mantissa == 0)
        return {};
    const bool negative = (a.mantissa < 0) != (b.mantissa < 0);
    const Integer denominator = abs_integer(b.mantissa);
    const auto shift = 256U + bits(denominator) - bits(a.mantissa);
    const Integer enlarged = abs_integer(a.mantissa) << shift;
    Integer quotient = enlarged / denominator;
    if (enlarged % denominator != 0 && upper != negative)
        ++quotient;
    return rounded(negative ? -quotient : quotient,
                   subtract_exponent(subtract_exponent(a.power, b.power), shift), upper);
}
inline Interval integer(std::int64_t x) {
    const Bound b{Integer{x}, 0};
    return {b, b};
}
inline Interval add(const Interval& a, const Interval& b, Work& w) {
    w.step();
    return {sum_bound(a.lower, b.lower, false), sum_bound(a.upper, b.upper, true)};
}
inline Interval subtract(const Interval& a, const Interval& b, Work& w) {
    w.step();
    return {sum_bound(a.lower, negate(b.upper), false), sum_bound(a.upper, negate(b.lower), true)};
}
inline Interval multiply(const Interval& a, const Interval& b, Work& w) {
    w.step();
    std::array<Bound, 4> low{
        product_bound(a.lower, b.lower, false), product_bound(a.lower, b.upper, false),
        product_bound(a.upper, b.lower, false), product_bound(a.upper, b.upper, false)};
    std::array<Bound, 4> high{
        product_bound(a.lower, b.lower, true), product_bound(a.lower, b.upper, true),
        product_bound(a.upper, b.lower, true), product_bound(a.upper, b.upper, true)};
    const auto less = [](const Bound& x, const Bound& y) { return compare(x, y) < 0; };
    return {*std::min_element(low.begin(), low.end(), less),
            *std::max_element(high.begin(), high.end(), less)};
}
inline Interval divide(const Interval& a, const Interval& b, Work& w) {
    w.step();
    if (compare(b.lower, {}) <= 0 && compare(b.upper, {}) >= 0)
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    std::array<Bound, 4> low{
        quotient_bound(a.lower, b.lower, false), quotient_bound(a.lower, b.upper, false),
        quotient_bound(a.upper, b.lower, false), quotient_bound(a.upper, b.upper, false)};
    std::array<Bound, 4> high{
        quotient_bound(a.lower, b.lower, true), quotient_bound(a.lower, b.upper, true),
        quotient_bound(a.upper, b.lower, true), quotient_bound(a.upper, b.upper, true)};
    const auto less = [](const Bound& x, const Bound& y) { return compare(x, y) < 0; };
    return {*std::min_element(low.begin(), low.end(), less),
            *std::max_element(high.begin(), high.end(), less)};
}
inline Interval scale(Interval value, std::int64_t exponent) {
    if (value.lower.mantissa != 0)
        value.lower.power = add_exponent(value.lower.power, exponent);
    if (value.upper.mantissa != 0)
        value.upper.power = add_exponent(value.upper.power, exponent);
    return value;
}
inline Interval absolute(const Interval& x) {
    if (compare(x.lower, {}) >= 0)
        return x;
    if (compare(x.upper, {}) <= 0)
        return {negate(x.upper), negate(x.lower)};
    return {{}, compare(negate(x.lower), x.upper) > 0 ? negate(x.lower) : x.upper};
}
struct RoundedWord {
    double word = 0;
    Bound exact;
};
// Encode IEEE bits directly, independent of the current FP rounding mode and FTZ.
inline RoundedWord word(const Bound& value, bool upward, bool nearest, Work& w) {
    w.step();
    if (value.mantissa == 0)
        return {0.0, {}};
    const bool negative = value.mantissa < 0;
    const Integer positive = abs_integer(value.mantissa);
    const auto top = add_exponent(value.power, static_cast<std::int64_t>(bits(positive)) - 1);
    if (top > 1023)
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    const auto grid = std::max<std::int64_t>(add_exponent(top, -52), -1074);
    Integer q;
    if (value.power >= grid)
        q = positive << static_cast<unsigned>(subtract_exponent(value.power, grid));
    else {
        const auto shift = subtract_exponent(grid, value.power);
        q = shift >= bits(positive) ? Integer{0}
                                    : Integer{positive >> static_cast<unsigned>(shift)};
        const bool remainder =
            shift >= bits(positive) || ((q << static_cast<unsigned>(shift)) != positive);
        bool increment = remainder && upward != negative;
        if (nearest) {
            increment = false;
            if (shift <= static_cast<std::int64_t>(bits(positive))) {
                const Integer discarded = positive - (q << static_cast<unsigned>(shift));
                const Integer half = Integer{1} << static_cast<unsigned>(shift - 1);
                increment = discarded > half || (discarded == half && ((q & 1) != 0));
            }
        }
        if (increment)
            ++q;
    }
    if (q == 0)
        return {std::bit_cast<double>(negative ? std::uint64_t{1} << 63 : 0), {}};
    auto exponent = grid;
    if (bits(q) > 53) {
        q >>= 1;
        exponent = add_exponent(exponent, 1);
    }
    const auto word_top = add_exponent(exponent, static_cast<std::int64_t>(bits(q)) - 1);
    if (word_top > 1023)
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    std::uint64_t wire = negative ? std::uint64_t{1} << 63 : 0;
    if (word_top >= -1022) {
        const auto aligned = q << static_cast<unsigned>(53 - bits(q));
        wire |= static_cast<std::uint64_t>(word_top + 1023) << 52;
        wire |= aligned.convert_to<std::uint64_t>() & ((std::uint64_t{1} << 52) - 1);
    } else
        wire |= q.convert_to<std::uint64_t>();
    return {std::bit_cast<double>(wire), {negative ? -q : q, exponent}};
}
struct Center {
    SiEnclosure value;
    Bound center, error;
};
inline Bound exact_sum(const Bound& a, const Bound& b) {
    if (a.mantissa == 0)
        return b;
    if (b.mantissa == 0)
        return a;
    const auto common = std::min(a.power, b.power);
    return {at_power(a, common, false) + at_power(b, common, false), common};
}
inline Interval relative_to(Interval value, std::int64_t q) {
    if (value.lower.mantissa != 0)
        value.lower.power = subtract_exponent(value.lower.power, q);
    if (value.upper.mantissa != 0)
        value.upper.power = subtract_exponent(value.upper.power, q);
    return value;
}
inline Center center(const Interval& original, bool two, Work& w) {
    const Bound sum = exact_sum(original.lower, original.upper);
    const Bound midpoint{sum.mantissa, add_exponent(sum.power, -1)};
    std::int64_t q =
        midpoint.mantissa == 0
            ? 0
            : add_exponent(midpoint.power, static_cast<std::int64_t>(bits(midpoint.mantissa)) - 1);
    const auto scaled = relative_to(Interval{midpoint, midpoint}, q);
    auto hi = word(scaled.lower, false, true, w);
    auto low = word(exact_sum(scaled.lower, negate(hi.exact)), false, true, w);
    if (!two)
        low = {0.0, {}};
    Bound c = sum_bound(hi.exact, low.exact, false);
    if (c.mantissa != 0) {
        const auto shift = add_exponent(c.power, static_cast<std::int64_t>(bits(c.mantissa)) - 1);
        if (shift != 0) {
            q = add_exponent(q, shift);
            hi = word(scale(Interval{hi.exact, hi.exact}, subtract_exponent(0, shift)).lower, false,
                      true, w);
            low = word(scale(Interval{low.exact, low.exact}, subtract_exponent(0, shift)).lower,
                       false, true, w);
            c = sum_bound(hi.exact, low.exact, false);
        }
    }
    const auto interval = relative_to(original, q);
    const auto left = absolute(Interval{sum_bound(interval.lower, negate(c), false),
                                        sum_bound(interval.lower, negate(c), true)});
    const auto right = absolute(Interval{sum_bound(interval.upper, negate(c), false),
                                         sum_bound(interval.upper, negate(c), true)});
    const auto radius =
        word(compare(left.upper, right.upper) > 0 ? left.upper : right.upper, true, false, w);
    const auto error = scale(Interval{radius.exact, radius.exact}, q).upper;
    const auto actual_center = scale(Interval{c, c}, q).lower;
    return {two ? SiEnclosure{ScaledTwoWord{hi.word, low.word, radius.word, q}}
                : SiEnclosure{ScaledBinary64{hi.word, radius.word, q}},
            actual_center, error};
}
} // namespace uvlas::io::detail::si
