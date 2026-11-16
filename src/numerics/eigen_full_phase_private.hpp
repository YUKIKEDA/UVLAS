#pragma once

#include <boost/multiprecision/cpp_int.hpp>

#include <bit>
#include <complex>
#include <cstdint>
#include <uvlas/numerics/eigen_metrics.hpp>

namespace uvlas::numerics::eigen::detail {
// Finite binary64 x=M*2^e has at most 53 significand bits and -1074<=e<=971.
// Each |complex|^2 * D^2 term has four significands (212 bits) and exponent
// between -4296 and 3884. A common -4296 origin needs at most 8393 bits after
// adding real and imaginary squares. This fixed allocator-free backend has 8448.
using FullPhaseInteger = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
    8448, 8448, boost::multiprecision::unsigned_magnitude, boost::multiprecision::checked, void>>;
struct FullBinaryPart {
    std::uint64_t significand;
    int exponent;
};
inline FullBinaryPart full_binary_part(double value) noexcept {
    const auto bits = std::bit_cast<std::uint64_t>(value);
    const auto encoded = (bits >> 52) & 0x7ff;
    const auto fraction = bits & 0x000fffffffffffffULL;
    return {encoded == 0 ? fraction : fraction | 0x0010000000000000ULL,
            encoded == 0 ? -1074 : static_cast<int>(encoded) - 1075};
}
inline FullPhaseInteger full_phase_term(double value, double scale) {
    const auto a = full_binary_part(value), b = full_binary_part(scale);
    FullPhaseInteger product{a.significand};
    product *= a.significand;
    product *= b.significand;
    product *= b.significand;
    product <<= static_cast<unsigned>(2 * a.exponent + 2 * b.exponent + 4296);
    return product;
}
inline FullPhaseInteger full_phase_square(std::complex<double> value, double scale) {
    return full_phase_term(value.real(), scale) + full_phase_term(value.imag(), scale);
}
// Inputs are finite, scales positive. Compare exact physical/D or D*physical
// magnitudes without sqrt, rounding a division, or treating a near-tie as a tie.
inline bool full_phase_coordinate_larger(std::complex<double> a, double da, std::complex<double> b,
                                         double db, PairSide side) {
    return side == PairSide::Right ? full_phase_square(a, db) > full_phase_square(b, da)
                                   : full_phase_square(a, da) > full_phase_square(b, db);
}
} // namespace uvlas::numerics::eigen::detail
