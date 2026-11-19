#pragma once

#include "../structure/exact_polynomial_detail.hpp"
#include "si_case_detail.hpp"
#include "si_rational_plan.hpp"

#include <bit>

namespace uvlas::io::detail {
// These source-derived frames enumerate large local payloads, including return
// temporaries. They are conservative stack-payload envelopes, not process RSS.
// Fixed Boost division has remainder and product backends; 1024 bits is below
// the locked multiplication backend's 40-limb Karatsuba threshold.
struct SiBackendDivisionFrame {
    si::Integer remainder, product;
};
struct SiRoundingFrame {
    si::Integer parameter, positive, quotient, comparison, negated_result, bits_argument;
};
struct SiQuotientFrame {
    si::Integer numerator, denominator, enlarged, quotient, modulus;
    SiBackendDivisionFrame backend;
    SiRoundingFrame rounding;
    si::Bound returned;
};
struct SiIntervalPrimitiveFrame {
    std::array<si::Bound, 4> low, high;
    si::Interval returned;
    SiQuotientFrame quotient;
};
struct SiWordFrame {
    si::Integer positive, quotient, discarded, half, shifted, bits_argument;
    si::Bound returned;
    SiBackendDivisionFrame backend;
};
struct SiCenterFrame {
    si::Bound sum, midpoint, center, error, actual_center;
    si::Interval scaled, interval, left, right;
    si::RoundedWord high, low, radius;
    SiWordFrame word;
    SiRoundingFrame rounding;
};
struct SiSequenceFrame {
    si::Interval power_base, power_result;
    std::array<si::Interval, 2> power_inputs, power_outputs;
    si::Interval pi_x, pi_square, pi_term, pi_sum, pi_remainder, pi_first, pi_second;
    std::array<si::Interval, 4> pi_inputs;
    std::array<si::Interval, 2> pi_outputs;
    std::array<si::Interval, 4> unit_factors;
    si::Interval unit_result, parsed_value;
};
struct SiConversionFrame {
    si::Interval original, truth, scale, absolute, relative, budget, pi;
    si::Center candidate;
};
inline constexpr std::uint64_t si_fixed_stack_bytes =
    sizeof(SiConversionFrame) + sizeof(SiSequenceFrame) +
    std::max(sizeof(SiIntervalPrimitiveFrame), sizeof(SiCenterFrame)) +
    rational_plan::certificate_stack_bytes;
inline constexpr std::uint64_t si_token_work =
    si::power_work_bound(std::numeric_limits<std::uint64_t>::max()) + 2;
inline constexpr std::uint64_t si_center_work =
    2ULL * (3 + 2 + 1); // high/low/radius, two renormalization words, compare.

struct SiExactLeafPlan {
    std::uint64_t cumulative = 0, simultaneous = 0, work = 0, external = 0;
};
inline SiExactLeafPlan si_exact_leaf_plan(SiContext& c, const SiRecipeView& recipe) {
    const auto cap = c.limits.exact().maximum_work_bits;
    const auto bounded_add = [cap](std::uint64_t a, std::uint64_t b) {
        return a >= cap || b >= cap - a ? cap : a + b;
    };
    const auto bounded_multiply = [cap](std::uint64_t a, std::uint64_t b) {
        return a && b > cap / a ? cap : a * b;
    };
    std::uint64_t bits = 0, characters = 0, largest_character_buffer = 0, storages = 0;
    bool zero = false;
    const auto decimal = [&](SiDecimalView value) {
        zero = zero || value.sign == 0;
        bits = bounded_add(bits, bounded_add(bounded_multiply(value.digits.size(), 4),
                                             bounded_multiply(si::magnitude(value.exponent10), 4)));
        const auto count = si_decimal_token_size(value);
        if (!count)
            throw count.error();
        const auto text = *count;
        if (text <= c.limits.exact().maximum_input_characters) {
            characters = c.add(characters, text);
            largest_character_buffer = std::max(largest_character_buffer, text);
        }
        storages = c.add(storages, 1);
    };
    if (const auto* value = std::get_if<SiDecimalView>(&recipe.value))
        decimal(*value);
    else {
        const auto& r = std::get<SiRationalView>(recipe.value);
        decimal(r.numerator);
        decimal(r.denominator);
        storages = c.add(storages, 1);
    }
    if (zero)
        return {};
    if (recipe.pi_power)
        return {0, 0, c.add(characters, 1), 0};
    // Exact numerator+denominator bit lengths of each reduced adopted root factor.
    constexpr std::array<std::uint64_t, 6> length_bits{2, 11, 8, 11, 20, 20};
    constexpr std::array<std::uint64_t, 4> mass_bits{2, 11, 11, 53};
    constexpr std::array<std::uint64_t, 3> time_bits{2, 11, 7};
    const std::array factors{length_bits[static_cast<std::size_t>(recipe.units.length)],
                             mass_bits[static_cast<std::size_t>(recipe.units.mass)],
                             time_bits[static_cast<std::size_t>(recipe.units.time)]};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const auto exponent = static_cast<std::int64_t>(recipe.dimensions.exponents[axis]);
        const auto absolute = si::magnitude(exponent);
        if (!absolute)
            continue;
        bits = bounded_add(bits, bounded_multiply(absolute, factors[axis]));
        // unit parse, power seed, reciprocal, bit products/squares, final multiply.
        const auto operations = 2 + static_cast<std::uint64_t>(exponent < 0) +
                                static_cast<std::uint64_t>(std::popcount(absolute)) +
                                static_cast<std::uint64_t>(std::bit_width(absolute));
        storages = c.add(storages, operations);
    }
    bits = bounded_add(bits, 2);
    // Four backend reductions plus the optional fixed certificate and admission walk.
    // Each Euclidean chain has at most twice the operand bit length in divisions.
    const auto rational_cost = c.add(c.multiply(bits, 16), 128);
    const auto work = c.add(c.multiply(storages, rational_cost),
                            c.add(c.multiply(characters, 8), c.add(2ULL * (3 * 64 + 1), 24)));
    const auto cumulative =
        c.add(characters, c.multiply(storages, sizeof(structure::detail::RationalStorage)));
    // The nested exact chain retains at most eight RationalStorage leaves. A
    // canonical decimal's character buffer is live in only one parse at a time.
    const auto simultaneous =
        c.add(largest_character_buffer, c.multiply(std::min<std::uint64_t>(storages, 8),
                                                   sizeof(structure::detail::RationalStorage)));
    const auto limb_bytes = c.multiply(c.add(bits, 63) / 64, sizeof(std::uint64_t));
    // Two source rational components and five proof temporaries; zero is exact.
    return {cumulative, simultaneous, work, c.multiply(limb_bytes, 7)};
}
inline std::uint64_t si_decimal_work(SiContext& c, const DecimalRange& value) {
    return c.add(c.multiply(value.count, si::digit_program.work + 1),
                 c.add(si::power_work_bound(si::magnitude(value.exponent)), 2));
}
inline std::uint64_t si_quantity_work(SiContext& c, const SchemaNode& node, bool rational,
                                      SiUnits units) {
    auto work = rational ? c.add(c.add(si_decimal_work(c, node.rational_numerator),
                                       si_decimal_work(c, node.rational_denominator)),
                                 1)
                         : si_decimal_work(c, node.number);
    work = c.add(work, 4); // rational root multipliers, including the degree denominator.
    for (const auto exponent : node.dimensions.exponents) {
        if (exponent)
            work = c.add(work, c.add(si::power_work_bound(si::magnitude(exponent)), 2));
    }
    if (units.angle == SiAngleUnit::Degree && node.dimensions.exponents[3])
        work = c.add(work,
                     c.add(si::power_work_bound(si::magnitude(node.dimensions.exponents[3])), 2));
    return c.add(work, 1 + 3 + si_center_work +
                           2); // SI product, budget operations, centers and row visits.
}
inline std::uint64_t si_rebind_work(SiContext& c, const SchemaData& source) {
    const auto& input = source.snapshot.values()[0];
    std::uint64_t work = 0;
    for (std::uint64_t node = 0; node < source.nodes.size(); ++node) {
        c.work.step();
        const auto& typed = source.nodes.values()[static_cast<std::size_t>(node)];
        if (!typed.assigned || node_kind(input, node) != InputNodeKind::Object)
            continue;
        work = c.add(work, 1);
        const auto& object =
            std::get<ObjectValue>(input.nodes.values()[static_cast<std::size_t>(node)].value);
        for (const auto& field : schema_records()[typed.record_kind].fields) {
            auto member = object.members.head;
            while (member >= 0) {
                c.work.step();
                work = c.add(work, 1);
                const auto& value = input.members.values()[static_cast<std::size_t>(member)];
                const std::string_view key{input.text.values().data() + value.key.offset,
                                           value.key.size};
                if (key.size() == field.name.size()) {
                    c.work.step(key.size());
                    work = c.add(work, key.size());
                    if (key == field.name)
                        break;
                }
                member = value.next;
            }
        }
    }
    return work;
}
} // namespace uvlas::io::detail
