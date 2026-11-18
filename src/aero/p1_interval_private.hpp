#pragma once

#include "exact_triangle_predicates.hpp"

#include <algorithm>
#include <array>
#include <limits>

namespace uvlas::aero::detail::p1 {
using Integer = boost::multiprecision::number<
    boost::multiprecision::cpp_int_backend<1024, 1024, boost::multiprecision::signed_magnitude,
                                           boost::multiprecision::checked, void>,
    boost::multiprecision::et_off>;
struct ArithmeticFailure {
    core::FailureCode code;
};
struct Work {
    std::uint64_t used = 0;
    std::uint64_t limit = 0;
    std::stop_token stop;
    void step() {
        if (stop.stop_requested())
            throw ArithmeticFailure{core::FailureCode::Cancelled};
        if (used == limit)
            throw ArithmeticFailure{core::FailureCode::ResourceExceeded};
        ++used;
    }
};
inline int exponent(std::int64_t value) {
    if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    return static_cast<int>(value);
}
inline Integer magnitude(const Integer& x) {
    return x < 0 ? -x : x;
}
inline unsigned bits(const Integer& x) {
    return x == 0 ? 0U : static_cast<unsigned>(boost::multiprecision::msb(magnitude(x))) + 1U;
}
/// An exact dyadic endpoint. Arithmetic normalizes outward to 256 significant bits.
struct Bound {
    Integer mantissa = 0;
    int power = 0;
};
inline Bound rounded(Integer value, int power, bool upper) {
    const unsigned length = bits(value);
    if (length == 0)
        return {};
    if (length > 256U) {
        const unsigned shift = length - 256U;
        const bool negative = value < 0;
        const Integer positive = magnitude(value);
        Integer quotient = positive >> shift;
        const bool remainder = (quotient << shift) != positive;
        if (remainder && (upper != negative))
            ++quotient;
        value = negative ? -quotient : quotient;
        power = exponent(static_cast<std::int64_t>(power) + shift);
    }
    return {std::move(value), power};
}
inline Bound negated(Bound x) {
    x.mantissa = -x.mantissa;
    return x;
}
inline Integer at_power(const Bound& x, int power, bool upper) {
    const auto shift = static_cast<std::int64_t>(x.power) - power;
    if (shift >= 0) {
        if (shift > 768)
            throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
        const Integer shifted = magnitude(x.mantissa) << static_cast<unsigned>(shift);
        return x.mantissa < 0 ? -shifted : shifted;
    }
    const auto right = -shift;
    const bool negative = x.mantissa < 0;
    const Integer positive = magnitude(x.mantissa);
    Integer quotient =
        right >= bits(positive) ? Integer{0} : Integer{positive >> static_cast<unsigned>(right)};
    const bool remainder = right >= bits(positive)
                               ? positive != 0
                               : (quotient << static_cast<unsigned>(right)) != positive;
    if (remainder && (upper != negative))
        ++quotient;
    return negative ? -quotient : quotient;
}
inline int compare(const Bound& left, const Bound& right) {
    if ((left.mantissa < 0) != (right.mantissa < 0))
        return left.mantissa < 0 ? -1 : 1;
    if (left.mantissa == 0 || right.mantissa == 0)
        return left.mantissa == right.mantissa ? 0 : (left.mantissa < right.mantissa ? -1 : 1);
    const auto a = static_cast<std::int64_t>(bits(left.mantissa)) + left.power;
    const auto b = static_cast<std::int64_t>(bits(right.mantissa)) + right.power;
    if (a != b)
        return (a < b) != (left.mantissa < 0) ? -1 : 1;
    const int common = std::min(left.power, right.power);
    Integer l = magnitude(left.mantissa) << static_cast<unsigned>(left.power - common);
    Integer r = magnitude(right.mantissa) << static_cast<unsigned>(right.power - common);
    if (left.mantissa < 0)
        l = -l;
    if (right.mantissa < 0)
        r = -r;
    return l == r ? 0 : (l < r ? -1 : 1);
}
inline Bound sum_bound(const Bound& a, const Bound& b, bool upper) {
    if (a.mantissa == 0)
        return b;
    if (b.mantissa == 0)
        return a;
    const auto top = std::max(static_cast<std::int64_t>(bits(a.mantissa)) + a.power,
                              static_cast<std::int64_t>(bits(b.mantissa)) + b.power);
    const int common = exponent(top - 258);
    return rounded(at_power(a, common, upper) + at_power(b, common, upper), common, upper);
}
inline Bound product_bound(const Bound& a, const Bound& b, bool upper) {
    return rounded(a.mantissa * b.mantissa, exponent(static_cast<std::int64_t>(a.power) + b.power),
                   upper);
}
inline Bound quotient_bound(const Bound& a, const Bound& b, bool upper) {
    if (b.mantissa == 0)
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    if (a.mantissa == 0)
        return {};
    const bool negative = (a.mantissa < 0) != (b.mantissa < 0);
    const Integer numerator = magnitude(a.mantissa), denominator = magnitude(b.mantissa);
    const int shift = 256 - static_cast<int>(bits(numerator)) + static_cast<int>(bits(denominator));
    const Integer enlarged = numerator << static_cast<unsigned>(shift);
    Integer quotient = enlarged / denominator;
    if (enlarged % denominator != 0 && (upper != negative))
        ++quotient;
    return rounded(negative ? -quotient : quotient,
                   exponent(static_cast<std::int64_t>(a.power) - b.power - shift), upper);
}
struct Interval {
    Bound lower, upper;
};
inline Interval integer(std::int64_t value) {
    const Bound b{Integer{value}, 0};
    return {b, b};
}
inline Interval exact(const ExactInteger& value, int power) {
    if (value == 0)
        return integer(0);
    const bool negative = value < 0;
    const ExactInteger positive = negative ? -value : value;
    const unsigned length = static_cast<unsigned>(boost::multiprecision::msb(positive)) + 1U;
    const unsigned shift = length > 256U ? length - 256U : 0U;
    const ExactInteger head = positive >> shift;
    const bool remainder = (head << shift) != positive;
    Integer low = head.convert_to<Integer>(), high = low;
    if (remainder)
        ++high;
    const int e = exponent(static_cast<std::int64_t>(power) + shift);
    return negative ? Interval{{-high, e}, {-low, e}} : Interval{{low, e}, {high, e}};
}
inline Interval scale(Interval value, int power) {
    for (auto* b : {&value.lower, &value.upper})
        b->power = exponent(static_cast<std::int64_t>(b->power) + power);
    return value;
}
inline Interval add(const Interval& a, const Interval& b, Work& work) {
    work.step();
    return {sum_bound(a.lower, b.lower, false), sum_bound(a.upper, b.upper, true)};
}
inline Interval negate(const Interval& a, Work& work) {
    work.step();
    return {negated(a.upper), negated(a.lower)};
}
inline Interval subtract(const Interval& a, const Interval& b, Work& work) {
    work.step();
    return {sum_bound(a.lower, negated(b.upper), false),
            sum_bound(a.upper, negated(b.lower), true)};
}
inline Interval multiply(const Interval& a, const Interval& b, Work& work) {
    work.step();
    std::array<Bound, 4> low{
        product_bound(a.lower, b.lower, false), product_bound(a.lower, b.upper, false),
        product_bound(a.upper, b.lower, false), product_bound(a.upper, b.upper, false)};
    std::array<Bound, 4> high{
        product_bound(a.lower, b.lower, true), product_bound(a.lower, b.upper, true),
        product_bound(a.upper, b.lower, true), product_bound(a.upper, b.upper, true)};
    const auto less = [](const Bound& x, const Bound& y) { return compare(x, y) < 0; };
    const Bound lower = (&a == &b && compare(a.lower, {}) <= 0 && compare(a.upper, {}) >= 0)
                            ? Bound{}
                            : *std::min_element(low.begin(), low.end(), less);
    return {lower, *std::max_element(high.begin(), high.end(), less)};
}
inline Interval divide(const Interval& a, const Interval& b, Work& work) {
    work.step();
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
inline Bound root_bound(const Bound& input, bool upper, Work& work) {
    if (input.mantissa < 0)
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    if (input.mantissa == 0)
        return {};
    const auto k = static_cast<std::int64_t>(bits(input.mantissa)) - 1 + input.power;
    const auto t = k >= 0 ? k / 2 : -((-k + 1) / 2);
    const auto shift = static_cast<std::int64_t>(input.power) + 512 - 2 * t;
    const Integer a = input.mantissa << static_cast<unsigned>(shift);
    Integer root = 0;
    for (int bit = 256; bit >= 0; --bit) {
        work.step();
        const Integer candidate = root + (Integer{1} << static_cast<unsigned>(bit));
        if (candidate * candidate <= a)
            root = candidate;
    }
    if (upper && root * root != a)
        ++root;
    return rounded(root, exponent(t - 256), upper);
}
inline Interval square_root(const Interval& x, Work& work) {
    return {root_bound(x.lower, false, work), root_bound(x.upper, true, work)};
}
inline Interval absolute(const Interval& x) {
    if (compare(x.lower, {}) >= 0)
        return x;
    if (compare(x.upper, {}) <= 0)
        return {negated(x.upper), negated(x.lower)};
    const Bound high = compare(negated(x.lower), x.upper) > 0 ? negated(x.lower) : x.upper;
    return {{}, high};
}
/// All series work is explicit. No platform libm ulp assumption enters an enclosure.
inline Interval logarithm_series(const Interval& t, Work& work) {
    const Interval square = multiply(t, t, work);
    Interval power = t, total = integer(0);
    for (unsigned j = 0; j < 128U; ++j) {
        total = add(total, divide(power, integer(2 * j + 1), work), work);
        power = multiply(power, square, work);
    }
    total = multiply(integer(2), total, work);
    const Interval remainder =
        divide(multiply(integer(2), absolute(power), work),
               multiply(integer(257), subtract(integer(1), square, work), work), work);
    const Interval radius{{}, remainder.upper};
    return {subtract(total, radius, work).lower, add(total, radius, work).upper};
}
inline Interval log_endpoint(const Bound& x, const Interval& ln2, Work& work) {
    if (x.mantissa <= 0)
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    const auto q = static_cast<std::int64_t>(bits(x.mantissa)) - 1 + x.power;
    const Bound normalized{x.mantissa, exponent(static_cast<std::int64_t>(x.power) - q)};
    const Interval r{normalized, normalized};
    const Interval t = divide(subtract(r, integer(1), work), add(r, integer(1), work), work);
    return add(logarithm_series(t, work), multiply(integer(q), ln2, work), work);
}
inline Interval log_one_plus(const Interval& x, const Interval& ln2, Work& work) {
    const Interval y = add(integer(1), x, work);
    return {log_endpoint(y.lower, ln2, work).lower, log_endpoint(y.upper, ln2, work).upper};
}
inline Interval atan_series(const Interval& z, Work& work) {
    const Interval square = multiply(z, z, work);
    Interval power = z, total = integer(0);
    for (unsigned j = 0; j < 128U; ++j) {
        const Interval term = divide(power, integer(2 * j + 1), work);
        total = (j % 2U == 0U) ? add(total, term, work) : subtract(total, term, work);
        power = multiply(power, square, work);
    }
    const Interval remainder = divide(absolute(power), integer(257), work);
    const Interval radius{{}, remainder.upper};
    return {subtract(total, radius, work).lower, add(total, radius, work).upper};
}
inline Interval pi_constant(Work& work) {
    return subtract(
        multiply(integer(16), atan_series(divide(integer(1), integer(5), work), work), work),
        multiply(integer(4), atan_series(divide(integer(1), integer(239), work), work), work),
        work);
}
inline Interval atan_nonnegative(const Interval& r, Work& work) {
    const Interval z = divide(
        r, add(integer(1), square_root(add(integer(1), multiply(r, r, work), work), work), work),
        work);
    return multiply(integer(2), atan_series(z, work), work);
}
inline Interval atan2_point(const Bound& y, const Bound& x, const Interval& pi, Work& work) {
    if (y.mantissa <= 0)
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    const Interval yy{y, y},
        xx{negated(x).mantissa > 0 ? negated(x) : x, negated(x).mantissa > 0 ? negated(x) : x};
    Interval angle;
    if (compare(xx.upper, y) <= 0)
        angle = subtract(scale(pi, -1), atan_nonnegative(divide(xx, yy, work), work), work);
    else
        angle = atan_nonnegative(divide(yy, xx, work), work);
    return x.mantissa >= 0 ? angle : subtract(pi, angle, work);
}
inline Interval atan_two(const Interval& y, const Interval& x, const Interval& pi, Work& work) {
    if (y.lower.mantissa <= 0)
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    const std::array<Interval, 4> corners{
        atan2_point(y.lower, x.lower, pi, work), atan2_point(y.lower, x.upper, pi, work),
        atan2_point(y.upper, x.lower, pi, work), atan2_point(y.upper, x.upper, pi, work)};
    Bound low = corners[0].lower, high = corners[0].upper;
    for (const auto& corner : corners) {
        if (compare(corner.lower, low) < 0)
            low = corner.lower;
        if (compare(corner.upper, high) > 0)
            high = corner.upper;
    }
    return {low, high};
}
} // namespace uvlas::aero::detail::p1
