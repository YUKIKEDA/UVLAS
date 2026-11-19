#include "../structure/exact_polynomial_detail.hpp"
#include "si_case_detail.hpp"
#include "si_case_plan.hpp"

namespace uvlas::io::detail {
namespace {
using Big = structure::detail::Integer;
std::uint64_t big_bits(const Big& value) {
    return value == 0 ? 0 : static_cast<std::uint64_t>(boost::multiprecision::msb(value)) + 1;
}
si::RoundedWord nearest(SiContext& c, const Big& signed_value, std::int64_t power) {
    c.work.step();
    if (signed_value == 0)
        return {0.0, {}};
    const bool negative = signed_value < 0;
    const Big positive = negative ? -signed_value : signed_value;
    const auto width = big_bits(positive);
    if (width > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
        c.reject(SiReason::Overflow, core::FailureCode::ResourceExceeded);
    const auto top = si::add_exponent(power, static_cast<std::int64_t>(width) - 1);
    const auto grid = std::max<std::int64_t>(-1074, si::add_exponent(top, -52));
    si::Integer head;
    if (power >= grid) {
        c.work.step();
        head = positive.convert_to<si::Integer>()
               << static_cast<unsigned>(si::subtract_exponent(power, grid));
    } else {
        const auto shift = si::subtract_exponent(grid, power);
        if (static_cast<std::uint64_t>(shift) > std::numeric_limits<unsigned>::max())
            c.reject(SiReason::Overflow, core::FailureCode::ResourceExceeded);
        if (static_cast<std::uint64_t>(shift) >= width)
            head = 0;
        else {
            c.work.step();
            const Big prefix = positive >> static_cast<std::uint64_t>(shift);
            head = prefix.convert_to<si::Integer>();
        }
        c.work.step();
        const bool half =
            static_cast<std::uint64_t>(shift) <= width &&
            boost::multiprecision::bit_test(positive, static_cast<unsigned>(shift - 1));
        c.work.step();
        const bool sticky = static_cast<std::uint64_t>(shift) > width ||
                            static_cast<std::uint64_t>(boost::multiprecision::lsb(positive)) <
                                static_cast<std::uint64_t>(shift - 1);
        if (half && (sticky || ((head & 1) != 0)))
            ++head;
    }
    return si::word({negative ? -head : head, grid}, false, true, c.work);
}
Big at_original_grid(SiContext& c, const si::Bound& word, std::int64_t grid) {
    if (word.mantissa == 0)
        return 0;
    const auto shift = si::subtract_exponent(word.power, grid);
    // Nearest words cannot need source bits below an exact dyadic source grid.
    if (shift < 0) {
        const auto right = si::magnitude(shift);
        if (right > si::bits(word.mantissa))
            c.reject(SiReason::UnresolvedPrecision, core::FailureCode::PrecisionUnresolved);
        const si::Integer quotient = si::abs_integer(word.mantissa) >> static_cast<unsigned>(right);
        if ((quotient << static_cast<unsigned>(right)) != si::abs_integer(word.mantissa))
            c.reject(SiReason::UnresolvedPrecision, core::FailureCode::PrecisionUnresolved);
        c.work.step();
        return (word.mantissa < 0 ? -quotient : quotient).convert_to<Big>();
    }
    c.work.step(2);
    Big result = word.mantissa.convert_to<Big>();
    result <<= static_cast<std::uint64_t>(shift);
    return result;
}
} // namespace
std::optional<SiEnclosure> si_exact_zero_center(SiContext& c, const SiRecipeView& recipe,
                                                bool two) {
    c.work.step();
    const auto leaf = si_exact_leaf_plan(c, recipe);
    c.work.step(leaf.work);
    c.opaque_exact_work = c.add(c.opaque_exact_work, leaf.work);
    c.exact_scratch_peak = std::max(c.exact_scratch_peak, leaf.simultaneous);
    c.require(SiQuota::ScratchBytes,
              c.add(c.add(c.scratch, c.automatic_scratch), leaf.simultaneous));
    c.require(SiQuota::CumulativeBytes,
              c.add(c.resources.allocated_bytes() - c.start_charged, leaf.cumulative));
    if (leaf.cumulative > c.resources.remaining_bytes())
        c.reject(SiReason::QuotaExceeded, core::FailureCode::ResourceExceeded,
                 SiQuota::CumulativeBytes, leaf.cumulative);
    if (const auto* estimate =
            std::get_if<structure::ExactEstimateAdmission>(&c.limits.exact().external_memory)) {
        if (leaf.external > estimate->bytes)
            c.reject(SiReason::QuotaExceeded, core::FailureCode::ResourceExceeded);
    }
    auto exact = materialize_rational(recipe, c.limits.exact(), c.resources);
    c.require(SiQuota::CumulativeBytes, c.resources.allocated_bytes() - c.start_charged);
    if (!exact) {
        if (exact.error().code == core::FailureCode::PrecisionUnsupported)
            return std::nullopt;
        throw exact.error();
    }
    if (exact->sign() == 0)
        return SiEnclosure{ScaledBinary64{0, 0, 0}};
    const auto& value = structure::detail::ExactFactory::value(*exact);
    const auto& source_numerator = boost::multiprecision::numerator(value);
    const auto& denominator = boost::multiprecision::denominator(value);
    const auto denominator_bits = exact->denominator_bits();
    const auto numerator_bits = exact->numerator_bits();
    // Existing exact limits govern these bounded backend temporaries as well as parsing.
    const auto needed = c.add(c.multiply(std::max(numerator_bits, denominator_bits), 4), 128);
    if (needed > c.limits.exact().maximum_work_bits)
        c.reject(SiReason::QuotaExceeded, core::FailureCode::ResourceExceeded);
    c.work.step(2);
    if (static_cast<std::uint64_t>(boost::multiprecision::lsb(denominator)) + 1 != denominator_bits)
        return std::nullopt;
    if (numerator_bits > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()) ||
        denominator_bits > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
        c.reject(SiReason::Overflow, core::FailureCode::ResourceExceeded);
    const auto grid = 1 - static_cast<std::int64_t>(numerator_bits);
    auto q = si::subtract_exponent(static_cast<std::int64_t>(numerator_bits),
                                   static_cast<std::int64_t>(denominator_bits));
    const auto high = nearest(c, source_numerator, grid);
    c.work.step();
    Big residual = source_numerator - at_original_grid(c, high.exact, grid);
    if (residual == 0)
        return SiEnclosure{ScaledBinary64{high.word, 0, q}};
    auto low = si::RoundedWord{0.0, {}};
    if (two)
        low = nearest(c, residual, grid);
    c.work.step();
    residual -= at_original_grid(c, low.exact, grid);
    if (residual != 0)
        return std::nullopt;
    auto h = high, l = low;
    if (si::compare(high.exact, {si::Integer{2}, 0}) >= 0 ||
        si::compare(high.exact, {si::Integer{-2}, 0}) <= 0) {
        // A carry to two occurs only with a residual toward the interval interior.
        // Inspect the complete exact center rather than the high word alone.
        const bool nonzero_low = l.exact.mantissa != 0;
        const bool points_inward =
            nonzero_low && ((h.exact.mantissa < 0) != (l.exact.mantissa < 0));
        if (!points_inward) {
            q = si::add_exponent(q, 1);
            h = si::word(si::scale({h.exact, h.exact}, -1).lower, false, true, c.work);
            l = si::word(si::scale({l.exact, l.exact}, -1).lower, false, true, c.work);
        }
    }
    return two ? SiEnclosure{ScaledTwoWord{h.word, l.word, 0, q}}
               : SiEnclosure{ScaledBinary64{h.word, 0, q}};
}
} // namespace uvlas::io::detail
