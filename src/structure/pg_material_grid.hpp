#pragma once

#include "pg_conditioning_factor.hpp"
#include "pg_reference_arithmetic.hpp"

namespace uvlas::structure::detail {

enum class GridEnd : std::uint8_t { Lower, Upper };

/// Outward rounding to the caller's unchanged dyadic assembly grid.
[[nodiscard]] inline ExactRational outward_material_grid(const ExactRational& value,
                                                         std::uint32_t bits, GridEnd endpoint,
                                                         ReferenceArithmetic& math) {
    math.cancel();
    if (value.sign() == 0)
        return {};
    const auto numerator_bits = partition_add(value.numerator_bits(), bits);
    const auto denominator_bits = value.denominator_bits();
    const auto width = std::max(numerator_bits, denominator_bits);
    using GridWide = RankInteger<32768>;
    if (width > std::numeric_limits<GridWide>::digits ||
        width > math.limits.algebra.maximum_work_bits)
        partition_fail(core::FailureCode::ResourceExceeded);
    const bool negative = value.sign() < 0;
    const bool away = negative ? endpoint == GridEnd::Lower : endpoint == GridEnd::Upper;
    constexpr auto short_width = std::numeric_limits<RankWide>::digits;
    constexpr auto limb_bits = std::numeric_limits<RankWord>::digits;
    math.charge(2 + numerator_bits / limb_bits + denominator_bits / limb_bits);
    if (width <= short_width) {
        const auto fraction = rank_fraction<RankWide>(value);
        const RankWide scaled = rank_abs(fraction.numerator) << bits;
        auto units = conditioning_divide(math, scaled, fraction.denominator);
        if (away && rank_product(math, units, fraction.denominator) != scaled)
            ++units;
        if (negative)
            units = -units;
        math.charge();
        return rank_materialize_dyadic(units, bits, math.leaf_limits, math.plan);
    }

    // The fixed wide path only needs shift and quotient/remainder. It uses
    // allocator-free integer slots, with the division charged for both operands.
    const auto fraction = rank_fraction<GridWide>(value);
    const GridWide scaled = rank_abs(fraction.numerator) << bits;
    math.charge(partition_mul(scaled.backend().size(), fraction.denominator.backend().size()));
    GridWide units, remainder;
    boost::multiprecision::divide_qr(scaled, fraction.denominator, units, remainder);
    if (away && remainder != 0)
        ++units;
    if (negative)
        units = -units;
    math.charge();
    return rank_materialize_dyadic(units, bits, math.leaf_limits, math.plan);
}

} // namespace uvlas::structure::detail
