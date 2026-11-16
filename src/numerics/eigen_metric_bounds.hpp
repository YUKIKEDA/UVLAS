#pragma once

#include "compensated.hpp"
#include "linear_solve_detail.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <complex>
#include <limits>

namespace uvlas::numerics::eigen::detail {
using Bounds = linear::detail::MagnitudeBounds;
using core::Result;
using linear::detail::add_magnitudes;
using linear::detail::multiply_magnitudes;
constexpr auto metric_operation = core::Operation::ExtractEigen;

inline auto arithmetic_failure(core::Failure error) noexcept {
    // These are internal arithmetic operations, not producer callbacks. A finite
    // input whose intermediate overflows has an unresolved enclosure.
    error.operation = metric_operation;
    if (error.code == core::FailureCode::Nonfinite)
        error.code = core::FailureCode::PrecisionUnresolved;
    return std::unexpected(error);
}

template <StorageScalar Scalar>
Result<Bounds> magnitude(Scalar high, Scalar low = {}) noexcept {
    auto result = linear::detail::actual_magnitude_bounds(high, low);
    if (!result)
        return arithmetic_failure(result.error());
    return *result;
}

// Directed positive division follows the existing linear metric policy, without
// changing that solver's arithmetic or exporting its anonymous implementation.
inline double divide_bound(double numerator, double denominator, bool upper) noexcept {
    constexpr double infinity = std::numeric_limits<double>::infinity();
    if (numerator == 0)
        return 0;
    if (denominator == 0)
        return upper ? infinity : std::numeric_limits<double>::max();
    const double rounded = numerator / denominator;
    if (!std::isfinite(rounded))
        return upper ? infinity : std::numeric_limits<double>::max();
    if (rounded == 0)
        return upper ? std::numeric_limits<double>::denorm_min() : 0;
    const auto product = two_product(rounded, denominator);
    if (!product)
        return std::nextafter(rounded, upper ? infinity : 0.0);
    const bool same =
        std::bit_cast<std::uint64_t>(product->hi) == std::bit_cast<std::uint64_t>(numerator);
    if ((upper && (product->hi > numerator || (same && product->lo >= 0))) ||
        (!upper && (product->hi < numerator || (same && product->lo <= 0))))
        return rounded;
    return std::nextafter(rounded, upper ? infinity : 0.0);
}

inline Bounds divide(Bounds numerator, double denominator) noexcept {
    return {divide_bound(numerator.lower, denominator, false),
            divide_bound(numerator.upper, denominator, true)};
}
inline double power_scale(double value, int exponent, bool upper) noexcept {
    if (value == 0)
        return 0;
    const double rounded = std::scalbn(value, exponent);
    if (!std::isfinite(rounded))
        return upper ? std::numeric_limits<double>::infinity() : std::numeric_limits<double>::max();
    if (rounded == 0)
        return upper ? std::numeric_limits<double>::denorm_min() : 0;
    const double restored = std::scalbn(rounded, -exponent);
    if ((upper && restored < value) || (!upper && restored > value))
        return std::nextafter(rounded, upper ? std::numeric_limits<double>::infinity() : 0.0);
    return rounded;
}
// Up to three positive scale factors. Exponents remain small signed integers,
// while mantissa arithmetic stays normal even when the physical factors do not.
inline Bounds scale(Bounds value, double numerator, double first_denominator = 1,
                    double second_denominator = 1) noexcept {
    const auto endpoint = [&](double input, bool upper) {
        if (input == 0)
            return 0.0;
        int exponent = 0, next = 0;
        const double leading = std::frexp(input, &exponent);
        const double factor = std::frexp(numerator, &next);
        exponent += next;
        const auto product = multiply_magnitudes({leading, leading}, {factor, factor});
        double result = upper ? product.upper : product.lower;
        for (double denominator : {first_denominator, second_denominator}) {
            const double mantissa = std::frexp(denominator, &next);
            exponent -= next;
            result = divide_bound(result, mantissa, upper);
        }
        return power_scale(result, exponent, upper);
    };
    return {endpoint(value.lower, false), endpoint(value.upper, true)};
}
inline double lower_subtract(double a, double b) noexcept {
    if (a <= b)
        return 0;
    const auto value = two_sum(a, -b);
    if (!value)
        return 0;
    return value->lo < 0 ? std::nextafter(value->hi, 0.0) : value->hi;
}
inline double distance_from_one(Bounds norm) noexcept {
    double result = 0;
    for (double value : {norm.lower, norm.upper}) {
        const auto difference = two_sum(value, -1.0);
        if (!difference)
            return std::numeric_limits<double>::infinity();
        const auto bound = magnitude(difference->hi, difference->lo);
        if (!bound)
            return std::numeric_limits<double>::infinity();
        result = std::max(result, bound->upper);
    }
    return result;
}
inline Result<Bounds> norm_step(Bounds accumulated, Bounds component) noexcept {
    const auto lower = magnitude(std::complex<double>{accumulated.lower, component.lower});
    if (!lower)
        return std::unexpected(lower.error());
    const auto upper = magnitude(std::complex<double>{accumulated.upper, component.upper});
    if (!upper)
        return std::unexpected(upper.error());
    return Bounds{lower->lower, upper->upper};
}

template <StorageScalar Scalar>
struct Evaluation {
    Scalar value;
    double error;
};
template <StorageScalar Scalar>
class Sum {
public:
    Result<void> add(Scalar value, double sign = 1) noexcept {
        if constexpr (std::same_as<Scalar, double>) {
            return translate(real_.add(sign * value));
        } else {
            auto real = translate(real_.add(sign * value.real()));
            if (!real)
                return real;
            return translate(imaginary_.add(sign * value.imag()));
        }
    }
    Result<void> product(Scalar a, Scalar b, double sign = 1) noexcept {
        if constexpr (std::same_as<Scalar, double>) {
            return translate(real_.product(a, b, sign));
        } else {
            for (const auto& [left, right, weight] :
                 {std::array{a.real(), b.real(), sign}, std::array{a.imag(), b.imag(), -sign}}) {
                auto result = translate(real_.product(left, right, weight));
                if (!result)
                    return result;
            }
            auto first = translate(imaginary_.product(a.real(), b.imag(), sign));
            if (!first)
                return first;
            return translate(imaginary_.product(a.imag(), b.real(), sign));
        }
    }
    [[nodiscard]] Result<Evaluation<Scalar>> finish() const noexcept {
        const auto real = real_.result();
        if (!real)
            return arithmetic_failure(real.error());
        if constexpr (std::same_as<Scalar, double>) {
            return Evaluation<Scalar>{real->value, real->error};
        } else {
            const auto imaginary = imaginary_.result();
            if (!imaginary)
                return arithmetic_failure(imaginary.error());
            return Evaluation<Scalar>{
                {real->value, imaginary->value},
                add_magnitudes({real->error, real->error}, {imaginary->error, imaginary->error})
                    .upper};
        }
    }

private:
    static Result<void> translate(Result<void> result) noexcept {
        if (!result)
            return arithmetic_failure(result.error());
        return {};
    }
    uvlas::numerics::detail::CompensatedAccumulator real_, imaginary_;
};
} // namespace uvlas::numerics::eigen::detail
