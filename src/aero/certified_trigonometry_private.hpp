#pragma once

#include "p1_interval_private.hpp"

namespace uvlas::aero::detail::p1 {

enum class PhaseUnits : std::uint8_t { Radians, PiMultiple };
struct SinCos {
    Interval sine, cosine;
};

/// On [-1,1], the alternating Taylor terms are decreasing in magnitude.
/// Term 32 bounds the truncation of the first 32 terms of both series.
inline SinCos small_sin_cos(const Interval& reduced, Work& work) {
    const Interval square = multiply(reduced, reduced, work);
    const Interval minus_square = negate(square, work);
    Interval sine = reduced;
    Interval cosine = integer(1);
    Interval sine_term = sine;
    Interval cosine_term = cosine;
    for (std::int64_t i = 1; i <= 32; ++i) {
        const auto even = 2 * i;
        sine_term =
            divide(multiply(sine_term, minus_square, work), integer(even * (even + 1)), work);
        cosine_term =
            divide(multiply(cosine_term, minus_square, work), integer((even - 1) * even), work);
        if (i != 32) {
            sine = add(sine, sine_term, work);
            cosine = add(cosine, cosine_term, work);
        }
    }
    const Interval sine_radius{{}, absolute(sine_term).upper};
    const Interval cosine_radius{{}, absolute(cosine_term).upper};
    return {{subtract(sine, sine_radius, work).lower, add(sine, sine_radius, work).upper},
            {subtract(cosine, cosine_radius, work).lower, add(cosine, cosine_radius, work).upper}};
}

/// Any integer quadrant identifies the same functions. Pick one using an interval
/// endpoint and certify that the whole remaining argument is in the series range.
inline SinCos sine_cosine(const Interval& phase, PhaseUnits units, const Interval& pi, Work& work) {
    const auto half_pi = scale(pi, -1);
    const auto quadrants =
        units == PhaseUnits::PiMultiple ? scale(phase, 1) : divide(phase, half_pi, work);
    const auto shifted = add(quadrants, scale(integer(1), -1), work);
    const auto& lower = shifted.lower;
    if (static_cast<std::int64_t>(bits(lower.mantissa)) + lower.power > 256)
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    const auto quadrant = at_power(lower, 0, false);
    const Interval selected{{quadrant, 0}, {quadrant, 0}};
    const auto reduced = units == PhaseUnits::PiMultiple
                             ? multiply(subtract(phase, scale(selected, -1), work), pi, work)
                             : subtract(phase, multiply(selected, half_pi, work), work);
    if (compare(absolute(reduced).upper, integer(1).upper) > 0)
        throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
    auto values = small_sin_cos(reduced, work);
    const Integer index = (quadrant % 4 + 4) % 4;
    switch (index.convert_to<unsigned>()) {
    case 0:
        return values;
    case 1:
        return {values.cosine, negate(values.sine, work)};
    case 2:
        return {negate(values.sine, work), negate(values.cosine, work)};
    default:
        return {negate(values.cosine, work), values.sine};
    }
}

} // namespace uvlas::aero::detail::p1
