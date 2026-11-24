#pragma once

#include "pg_rank_storage.hpp"

#include <bit>
#include <limits>

namespace uvlas::structure::detail {

// A product of binary64 words is an integer multiple of 2^-2148. It takes
// fewer than 4196 upper bits before accumulation. Four products per
// contracted Index, plus the error terms, fit in 4288 bits.
using VolumeInteger = RankInteger<4288>;
using VolumeMantissaProduct = RankInteger<128>;
inline constexpr std::uint64_t volume_sign = std::uint64_t{1} << 63U;
inline constexpr std::uint64_t volume_infinity = 0x7ff0000000000000ULL;
inline constexpr std::uint64_t volume_normal_bit = std::uint64_t{1} << 52U;
inline constexpr int volume_product_shift = 2148;
inline constexpr auto volume_operation = core::Operation::EvaluateStructure;
static_assert(std::numeric_limits<double>::is_iec559 && sizeof(double) == sizeof(std::uint64_t));

[[nodiscard]] inline bool volume_finite(double value) noexcept {
    return (std::bit_cast<std::uint64_t>(value) & ~volume_sign) < volume_infinity;
}
struct VolumeWord {
    std::uint64_t mantissa;
    int exponent;
    bool negative;
};
[[nodiscard]] inline VolumeWord volume_decode(double value) noexcept {
    const auto bits = std::bit_cast<std::uint64_t>(value);
    const auto exponent = static_cast<int>((bits >> 52U) & 0x7ffU);
    const auto fraction = bits & (volume_normal_bit - 1);
    return {exponent == 0 ? fraction : fraction | volume_normal_bit,
            exponent == 0 ? -1074 : exponent - 1075, (bits & volume_sign) != 0};
}

[[nodiscard]] inline VolumeInteger volume_lift(VolumeWord value) {
    if (value.mantissa == 0)
        return {};
    VolumeInteger exact = value.mantissa;
    exact <<= static_cast<unsigned>(value.exponent + volume_product_shift);
    if (value.negative)
        exact = -exact;
    return exact;
}

[[nodiscard]] inline VolumeInteger volume_product(VolumeWord a, VolumeWord b) {
    // Only the two 53-bit mantissas are multiplied. The wider fixed array
    // uses shifts and signed additions, with no large-multiply workspace.
    if (a.mantissa == 0 || b.mantissa == 0)
        return {};
    auto exact = VolumeInteger{VolumeMantissaProduct{a.mantissa} * b.mantissa};
    exact <<= static_cast<unsigned>(a.exponent + b.exponent + volume_product_shift);
    if (a.negative != b.negative)
        exact = -exact;
    return exact;
}

[[nodiscard]] inline core::Result<double> volume_round(const VolumeInteger& exact) {
    if (exact == 0)
        return 0;
    const bool negative = exact < 0;
    const VolumeInteger magnitude = negative ? -exact : exact;
    const auto bit = boost::multiprecision::msb(magnitude);
    auto shift = bit > 1126U ? bit - 52U : 1074U;
    auto mantissa = (magnitude >> shift).convert_to<std::uint64_t>();
    const auto discarded = magnitude - (VolumeInteger{mantissa} << shift);
    const auto halfway = VolumeInteger{1} << (shift - 1);
    if (discarded > halfway || (discarded == halfway && mantissa % 2 != 0))
        ++mantissa;
    if (mantissa == volume_normal_bit * 2) {
        mantissa >>= 1U;
        ++shift;
    }
    if (shift > 3119U)
        return std::unexpected(
            core::Failure{core::FailureCode::PrecisionUnresolved, volume_operation});
    const auto biased_exponent = mantissa < volume_normal_bit ? 0 : std::uint64_t{shift - 1073U};
    const auto fraction = mantissa & (volume_normal_bit - 1);
    const auto word = (negative ? volume_sign : 0) | (biased_exponent << 52U) | fraction;
    return std::bit_cast<double>(word);
}

struct VolumeProjection {
    double high, low;
};

[[nodiscard]] inline core::Result<VolumeProjection>
volume_finalize(VolumeInteger total, const VolumeInteger& reference_error,
                const VolumeInteger& allowed, bool two_words) {
    const auto high = volume_round(total);
    if (!high)
        return std::unexpected(high.error());
    total -= volume_lift(volume_decode(*high));
    const auto low = two_words ? volume_round(total) : core::Result<double>{0.0};
    if (!low)
        return std::unexpected(low.error());
    total -= volume_lift(volume_decode(*low));
    if (rank_abs(total) + reference_error > allowed)
        return std::unexpected(
            core::Failure{core::FailureCode::PrecisionUnresolved, volume_operation});
    return VolumeProjection{*high, *low};
}

} // namespace uvlas::structure::detail
