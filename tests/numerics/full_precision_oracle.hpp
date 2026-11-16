#pragma once

// Independent test-only dyadic arithmetic. No UVLAS residual/EFT helper or
// numerical backend participates in this oracle.
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace uvlas::test_support::full_precision {
using Rational = boost::multiprecision::cpp_rational;
using Complex = std::complex<double>;
struct Audit {
    Rational l2_squared_ratio;
    Rational rhs_norm_squared;
    Rational backward_upper;
};
struct ExactComplex {
    Rational real;
    Rational imag;
};
namespace detail {
using Integer = boost::multiprecision::number<
    boost::multiprecision::cpp_int_backend<1024, 1024, boost::multiprecision::signed_magnitude,
                                           boost::multiprecision::checked, void>,
    boost::multiprecision::et_off>;
using SquareInteger = boost::multiprecision::number<
    boost::multiprecision::cpp_int_backend<2048, 2048, boost::multiprecision::unsigned_magnitude,
                                           boost::multiprecision::checked, void>,
    boost::multiprecision::et_off>;
inline constexpr std::size_t maximum_order = 4096;
inline void require(bool condition, const char* message) {
    if (!condition)
        throw std::domain_error(message);
}
struct Scalar {
    std::uint64_t mantissa;
    int exponent;
    bool negative;
};
inline Scalar decode(double value) {
    const auto bits = std::bit_cast<std::uint64_t>(value);
    const auto magnitude = bits & UINT64_C(0x7fffffffffffffff);
    require((magnitude >> 52) != 0x7ff, "Nonfinite exact-oracle input");
    if (magnitude == 0)
        return {0, 0, false};
    require(magnitude >= std::bit_cast<std::uint64_t>(0x1p-200) &&
                magnitude <= std::bit_cast<std::uint64_t>(0x1p200),
            "Exact-oracle scalar lies outside [2^-200,2^200]");
    return {(magnitude & UINT64_C(0x000fffffffffffff)) | UINT64_C(0x0010000000000000),
            static_cast<int>(magnitude >> 52) - 1075, (bits >> 63) != 0};
}
inline Integer in_units(Scalar value, int unit_exponent) {
    if (value.mantissa == 0)
        return 0;
    const int shift = value.exponent - unit_exponent;
    require(shift >= 0, "Invalid dyadic oracle unit");
    Integer result{value.mantissa};
    result <<= static_cast<unsigned>(shift);
    return value.negative ? -result : result;
}
inline Integer absolute(const Integer& value) {
    return value < 0 ? -value : value;
}
struct Actual {
    Integer real;
    Integer imag;
    Integer maximum;
    bool negative_real;
    bool negative_imag;
    Actual(const Integer& re, const Integer& im)
        : real(absolute(re)),
          imag(absolute(im)),
          maximum(std::max(real, imag)),
          negative_real(re < 0),
          negative_imag(im < 0) {}
};
inline std::vector<Actual> actual_words(std::span<const Complex> high,
                                        std::span<const Complex> low) {
    require(high.size() <= maximum_order && (low.empty() || low.size() == high.size()),
            "Invalid exact-oracle word shape");
    std::vector<Actual> result;
    result.reserve(high.size());
    for (std::size_t i = 0; i < high.size(); ++i) {
        Integer re = in_units(decode(high[i].real()), -252);
        Integer im = in_units(decode(high[i].imag()), -252);
        if (!low.empty()) {
            re += in_units(decode(low[i].real()), -252);
            im += in_units(decode(low[i].imag()), -252);
        }
        // Form hi+lo before taking its absolute value. Cancellation is physical.
        result.emplace_back(re, im);
    }
    return result;
}
inline Integer magnitude_product(Scalar coefficient, const Integer& actual_magnitude) {
    if (coefficient.mantissa == 0 || actual_magnitude == 0)
        return 0;
    // Multiply by the at-most-53-bit significand before shifting. Multiplication
    // of a previously shifted coefficient by the full word would do unnecessary
    // wide-by-wide work for each of the 4096^2 entries.
    Integer result = actual_magnitude;
    result *= coefficient.mantissa;
    const int shift = coefficient.exponent + 252;
    require(shift >= 0, "Invalid exact product exponent");
    result <<= static_cast<unsigned>(shift);
    return result;
}
inline Integer signed_product(Scalar coefficient, const Integer& actual_magnitude,
                              bool actual_negative) {
    const Integer result = magnitude_product(coefficient, actual_magnitude);
    return coefficient.negative != actual_negative ? -result : result;
}
inline SquareInteger square(const Integer& value) {
    SquareInteger positive{absolute(value)};
    positive *= positive;
    return positive;
}
template <class IntegerType>
inline Rational dyadic(const IntegerType& value, unsigned unit_bits) {
    // Boost 1.92 mixed cpp_int division can retain 0/k, which does not compare
    // equal to integral zero. Keep the exact zero canonical before that operation.
    if (value == 0)
        return Rational{0};
    Rational result{boost::multiprecision::cpp_int{value}};
    result /= boost::multiprecision::cpp_int{1} << unit_bits;
    return result;
}
inline Rational ratio(const Integer& numerator, const Integer& denominator) {
    require(denominator > 0, "Nonpositive exact ratio denominator");
    if (numerator == 0)
        return Rational{0};
    Rational result{boost::multiprecision::cpp_int{numerator}};
    result /= boost::multiprecision::cpp_int{denominator};
    return result;
}
} // namespace detail

// Every nonzero input scalar (including the rounded effective coefficients) is
// checked in [2^-200,2^200]. Its integer significand exponent is >=-252. Actual
// words are exact multiples of 2^-252 and have magnitude <=2^201 after hi+lo.
// For n<=4096 a residual component is <2^415: its integer magnitude is <2^919
// in units 2^-504. The integer sum of all real/imag squares is <2^1851 in
// units 2^-1008 (at most 919 and 1851 bits respectively).
// The L1 residual numerator and denominator cross-products also fit 2048 bits.
// Checked fixed backends enforce these limits; no arbitrary-precision summation
// is used for the O(n^2) path. Dynamic rationals are constructed only afterward.
// Aeff uses the original packed column-major real base. Only its final column
// changes: real=round64(base+update.real), imag=update.imag. The update transpose
// is nonconjugated. Empty low means an exact zero low plane.
inline Audit audit_full_column(std::span<const double> base_colmajor,
                               std::span<const Complex> update, std::span<const Complex> rhs,
                               std::span<const Complex> high, std::span<const Complex> low = {}) {
    using namespace detail;
    const std::size_t n = update.size();
    require(n > 0 && n <= maximum_order && base_colmajor.size() == n * n && rhs.size() == n &&
                high.size() == n && (low.empty() || low.size() == n),
            "Invalid full-column oracle shape");
    require(std::fegetround() == FE_TONEAREST,
            "Oracle effective coefficients require nearest rounding");
    const auto x = actual_words(high, low);
    for (const auto value : update) {
        (void)decode(value.real());
        (void)decode(value.imag());
    }
    SquareInteger residual_squared = 0, rhs_squared = 0;
    Integer worst_numerator = 0, worst_denominator = 1;
    for (std::size_t row = 0; row < n; ++row) {
        const auto br = decode(rhs[row].real()), bi = decode(rhs[row].imag());
        Integer real = in_units(br, -504), imag = in_units(bi, -504);
        rhs_squared += square(real) + square(imag);
        Integer denominator = std::max(absolute(real), absolute(imag));
        for (std::size_t col = 0; col < n; ++col) {
            const double original = base_colmajor[row + col * n];
            const auto base = decode(original); // Validate even if the actual word is zero.
            Scalar ar = base, ai{0, 0, false};
            Scalar maximum = base;
            if (col + 1 == n) {
                const double effective_real = original + update[row].real();
                ar = decode(effective_real);
                ai = decode(update[row].imag());
                maximum = std::abs(effective_real) >= std::abs(update[row].imag()) ? ar : ai;
            }
            const auto& word = x[col];
            real -= signed_product(ar, word.real, word.negative_real);
            imag -= signed_product(ar, word.imag, word.negative_imag);
            if (ai.mantissa != 0) {
                real += signed_product(ai, word.imag, word.negative_imag);
                imag -= signed_product(ai, word.real, word.negative_real);
            }
            denominator += magnitude_product(maximum, word.maximum);
        }
        residual_squared += square(real) + square(imag);
        const Integer numerator = absolute(real) + absolute(imag);
        if (denominator == 0) {
            require(numerator == 0, "Zero original denominator with nonzero exact residual");
        } else {
            SquareInteger left{numerator}, right{worst_numerator};
            left *= SquareInteger{worst_denominator};
            right *= SquareInteger{denominator};
            if (left > right) {
                worst_numerator = numerator;
                worst_denominator = denominator;
            }
        }
    }
    Rational residual_ratio = 0;
    if (rhs_squared == 0) {
        require(residual_squared == 0, "Nonzero residual has infinite zero-RHS relative ratio");
    } else if (residual_squared != 0) {
        residual_ratio = Rational{boost::multiprecision::cpp_int{residual_squared}} /
                         boost::multiprecision::cpp_int{rhs_squared};
    }
    return {residual_ratio, dyadic(rhs_squared, 1008), ratio(worst_numerator, worst_denominator)};
}
inline Rational exact_rhs_norm_squared(std::span<const Complex> rhs) {
    using namespace detail;
    require(rhs.size() <= maximum_order, "Invalid RHS oracle shape");
    SquareInteger sum = 0;
    for (const auto value : rhs) {
        sum += square(in_units(decode(value.real()), -252));
        sum += square(in_units(decode(value.imag()), -252));
    }
    return dyadic(sum, 504);
}
// Return a positive binary64 lower bound on relative*||rhs||2. Big100 supplies
// a candidate only: the squared comparison below proves its direction exactly.
inline double absolute_allowance(std::span<const Complex> rhs, double relative) {
    using Big = boost::multiprecision::cpp_bin_float_100;
    detail::require(std::isfinite(relative) && relative > 0, "Invalid relative allowance");
    const Rational norm_squared = exact_rhs_norm_squared(rhs);
    detail::require(norm_squared > 0, "A positive RHS norm is required");
    const Rational target_squared = Rational{relative} * Rational{relative} * norm_squared;
    const Big target = target_squared.convert_to<Big>();
    double result = sqrt(target).convert_to<double>();
    detail::require(std::isfinite(result) && result > 0, "Unrepresentable positive allowance");
    while (Rational{result} * Rational{result} > target_squared) {
        result = std::nextafter(result, 0.0);
        detail::require(result > 0, "Unrepresentable downward allowance");
    }
    return result;
}
inline double conservative_l2_allowance(std::span<const Complex> rhs, double relative) {
    return absolute_allowance(rhs, relative);
}
inline double to_upper(const Rational& nonnegative) {
    if (boost::multiprecision::numerator(nonnegative) == 0)
        return 0;
    detail::require(nonnegative >= 0, "A nonnegative oracle upper bound is required");
    double result = nonnegative.convert_to<double>();
    detail::require(std::isfinite(result), "Unrepresentable oracle upper bound");
    if (result == 0)
        result = std::numeric_limits<double>::denorm_min();
    while (Rational{result} < nonnegative) {
        result = std::nextafter(result, std::numeric_limits<double>::infinity());
        detail::require(std::isfinite(result), "Overflow in upward oracle conversion");
    }
    return result;
}
// Plain nonconjugated complex dot. A conjugated observation must explicitly
// supply conjugated coefficients; no convention is inferred from the data.
inline ExactComplex exact_complex_dot(std::span<const Complex> coefficients,
                                      std::span<const Complex> high,
                                      std::span<const Complex> low = {}) {
    using namespace detail;
    require(coefficients.size() == high.size() && high.size() <= maximum_order,
            "Invalid exact dot shape");
    const auto x = actual_words(high, low);
    Integer real = 0, imag = 0;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        const auto ar = decode(coefficients[i].real()), ai = decode(coefficients[i].imag());
        real += signed_product(ar, x[i].real, x[i].negative_real);
        real -= signed_product(ai, x[i].imag, x[i].negative_imag);
        imag += signed_product(ar, x[i].imag, x[i].negative_imag);
        imag += signed_product(ai, x[i].real, x[i].negative_real);
    }
    return {dyadic(real, 504), dyadic(imag, 504)};
}
inline Rational squared_distance(const ExactComplex& value, Complex reference) {
    detail::require(std::isfinite(reference.real()) && std::isfinite(reference.imag()),
                    "Nonfinite frozen observation");
    const Rational real = value.real - Rational{reference.real()};
    const Rational imag = value.imag - Rational{reference.imag()};
    return real * real + imag * imag;
}
inline double observation_error_upper(std::span<const Complex> coefficients,
                                      std::span<const Complex> high, std::span<const Complex> low,
                                      Complex reference, double scale) {
    detail::require(std::isfinite(scale) && scale > 0 && std::isfinite(reference.real()) &&
                        std::isfinite(reference.imag()),
                    "Invalid frozen observation scale/value");
    const auto value = exact_complex_dot(coefficients, high, low);
    const Rational error =
        abs(value.real - Rational{reference.real()}) + abs(value.imag - Rational{reference.imag()});
    return to_upper(error / Rational{scale});
}

// Independent small Rational arithmetic checks the bit-decoded fixed-width
// oracle, including complex signs, actual-word cancellation and rounded Aeff.
// Call this from a dedicated Catch2 case; failures throw rather than disappear
// under NDEBUG. No product solver or production helper is used by this check.
inline void self_test() {
    const auto check = [](bool value) { detail::require(value, "Exact oracle self-test failed"); };
    const std::array<double, 4> base{2, -1, 1, 3};
    const std::array<Complex, 2> update{Complex{.25, 2}, Complex{-.5, -1}};
    const std::array<Complex, 2> rhs{Complex{3, .5}, Complex{-2, 3}};
    const std::array<Complex, 2> high{Complex{1, 1}, Complex{2, -1}};
    const std::array<Complex, 2> low{Complex{0x1p-80, -0x1p-81}, Complex{-0x1p-82, 0x1p-83}};
    const auto actual = audit_full_column(base, update, rhs, high, low);
    Rational rs = 0, bs = 0, worst = 0;
    for (std::size_t i = 0; i < 2; ++i) {
        Rational rr{rhs[i].real()}, ri{rhs[i].imag()};
        Rational den =
            std::max(Rational{std::abs(rhs[i].real())}, Rational{std::abs(rhs[i].imag())});
        bs += rr * rr + ri * ri;
        for (std::size_t j = 0; j < 2; ++j) {
            Rational ar{base[i + 2 * j]}, ai = 0;
            if (j == 1) {
                ar += Rational{update[i].real()}; // This literal's quarter sums are exact.
                ai = Rational{update[i].imag()};
            }
            const Rational xr = Rational{high[j].real()} + Rational{low[j].real()};
            const Rational xi = Rational{high[j].imag()} + Rational{low[j].imag()};
            rr -= ar * xr - ai * xi;
            ri -= ar * xi + ai * xr;
            const Rational am = std::max(Rational{abs(ar)}, Rational{abs(ai)});
            const Rational xm = std::max(Rational{abs(xr)}, Rational{abs(xi)});
            den += am * xm;
        }
        rs += rr * rr + ri * ri;
        const Rational ratio = (abs(rr) + abs(ri)) / den;
        if (ratio > worst)
            worst = ratio;
    }
    check(actual.rhs_norm_squared == bs && actual.l2_squared_ratio == rs / bs &&
          actual.backward_upper == worst);
    const std::array<double, 1> identity{1};
    const std::array<Complex, 1> zero{Complex{}}, one{Complex{1, 0}}, minus_one{Complex{-1, 0}};
    const auto cancelled = audit_full_column(identity, zero, one, one, minus_one);
    check(cancelled.l2_squared_ratio == 1 && cancelled.backward_upper == 1);
    const std::array<Complex, 1> edge_high{Complex{std::nextafter(0x1p-200, 1.0), 0}};
    const std::array<Complex, 1> edge_low{Complex{-0x1p-200, 0}};
    const Rational tiny = Rational{1} / (boost::multiprecision::cpp_int{1} << 252);
    const auto edge = audit_full_column(identity, zero, one, edge_high, edge_low);
    check(edge.l2_squared_ratio == (1 - tiny) * (1 - tiny));
    check(edge.backward_upper == (1 - tiny) / (1 + tiny));
    const std::array<Complex, 1> rounded_update{Complex{0x1p-53, 0x1p-80}};
    const std::array<Complex, 1> rounded_rhs{Complex{1, 0x1p-80}};
    check(audit_full_column(identity, rounded_update, rounded_rhs, one).l2_squared_ratio == 0);
    const std::array<Complex, 2> coefficients{Complex{1, 2}, Complex{3, -4}};
    const std::array<Complex, 2> dot_high{Complex{2, -1}, Complex{-1, 3}};
    const std::array<Complex, 2> dot_low{Complex{.5, .25}, Complex{-.25, .5}};
    const auto dot = exact_complex_dot(coefficients, dot_high, dot_low);
    check(dot.real == Rational{57} / 4 && dot.imag == Rational{79} / 4);
    check(observation_error_upper(coefficients, dot_high, dot_low, Complex{14.25, 19.75}, .2) == 0);
    check(exact_rhs_norm_squared(zero) == 0);
    const auto zero_dot = exact_complex_dot(one, zero);
    check(zero_dot.real == 0 && zero_dot.imag == 0);
    Rational mixed_zero{boost::multiprecision::cpp_int{0}};
    mixed_zero /= boost::multiprecision::cpp_int{1024};
    check(to_upper(mixed_zero) == 0);
    const std::array<Complex, 1> pythagorean{Complex{3, 4}};
    check(exact_rhs_norm_squared(pythagorean) == 25);
    check(absolute_allowance(pythagorean, .125) == .625);
    check(conservative_l2_allowance(pythagorean, 1) == 5);
    const Rational below_binary64 = Rational{1} / (boost::multiprecision::cpp_int{1} << 1100);
    const double upper = to_upper(below_binary64);
    check(upper > 0 && Rational{upper} >= below_binary64);
}
} // namespace uvlas::test_support::full_precision
