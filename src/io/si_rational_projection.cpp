#include "si_rational_projection.hpp"

#include "si_interval.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <optional>
#include <type_traits>

namespace uvlas::io::detail {
namespace {
constexpr auto operation = core::Operation::EvaluateAero;
constexpr std::uint64_t endpoint_bits = 256;
using enum core::FailureCode;

struct Magnitude {
    std::span<const std::uint32_t> words;
    std::uint64_t width = 0, lowest = 0;

    [[nodiscard]] bool bit(std::uint64_t index, si::Work& work) const {
        work.step();
        return ((words[static_cast<std::size_t>(index / 32)] >> (index % 32)) & 1U) != 0;
    }
};

struct Rational {
    int sign = 0;
    Magnitude numerator, denominator;
};

Magnitude inspect(std::span<const std::uint32_t> words, std::uint64_t cap, si::Work& work) {
    work.step();
    if (words.empty())
        return {words};
    const auto bounded_cap =
        std::min(cap, static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()));
    const auto maximum_words = bounded_cap / 32 + (bounded_cap % 32 != 0 ? 1U : 0U);
    if (words.size() > maximum_words)
        throw si::ArithmeticFailure{ResourceExceeded};
    if (words.back() == 0)
        throw si::ArithmeticFailure{InvalidInput};
    const auto width = (static_cast<std::uint64_t>(words.size()) - 1) * 32 +
                       static_cast<std::uint64_t>(std::bit_width(words.back()));
    if (width > bounded_cap)
        throw si::ArithmeticFailure{ResourceExceeded};
    for (std::size_t i = 0; i < words.size(); ++i) {
        work.step();
        if (words[i] != 0)
            return {words, width,
                    static_cast<std::uint64_t>(i) * 32 +
                        static_cast<std::uint64_t>(std::countr_zero(words[i]))};
    }
    throw si::ArithmeticFailure{InvalidInput};
}

Rational inspect_fraction(SiRationalWordsView value, std::uint64_t cap, si::Work& work) {
    if (value.denominator.empty() ||
        (value.numerator.empty() ? value.sign != 0 : (value.sign != -1 && value.sign != 1)))
        throw si::ArithmeticFailure{InvalidInput};
    return {value.sign, inspect(value.numerator, cap, work), inspect(value.denominator, cap, work)};
}

si::Interval enclosure(const Magnitude& value, si::Work& work) {
    if (value.width == 0)
        return {};
    const auto shift = value.width > endpoint_bits ? value.width - endpoint_bits : 0;
    const auto first = static_cast<std::size_t>(shift / 32);
    si::Integer head{};
    for (std::size_t i = value.words.size(); i > first; --i) {
        work.step();
        head <<= 32;
        head += value.words[i - 1];
    }
    head >>= static_cast<unsigned>(shift % 32);
    const auto q = static_cast<std::int64_t>(shift);
    const auto discarded = value.lowest < shift;
    return {{head, q}, si::rounded(head + (discarded ? 1 : 0), q, true)};
}

si::Interval enclosure(const Rational& value, si::Work& work) {
    const auto numerator = enclosure(value.numerator, work);
    const auto denominator = enclosure(value.denominator, work);
    const auto result = si::divide(numerator, denominator, work);
    return value.sign < 0 ? si::Interval{si::negate(result.upper), si::negate(result.lower)}
                          : result;
}

SiEnclosure zero_enclosure(bool two_word) noexcept {
    if (two_word)
        return ScaledTwoWord{0.0, 0.0, 0.0, 0};
    return ScaledBinary64{0.0, 0.0, 0};
}

template <class Bit>
std::uint64_t significant(std::uint64_t first, std::uint64_t end, Bit bit) {
    std::uint64_t result = 0;
    for (auto i = end; i > first; --i)
        result = (result << 1U) | static_cast<std::uint64_t>(bit(i - 1));
    return result;
}

// Preserve an exactly representable low word even across a gap exceeding the interval precision.
std::optional<SiEnclosure> exact_dyadic(const Rational& value, bool two_word, si::Work& work) {
    work.step();
    if (value.sign == 0)
        return zero_enclosure(two_word);
    if (value.denominator.lowest != value.denominator.width - 1)
        return std::nullopt;

    const auto& magnitude = value.numerator;
    const auto bit = [&](auto index) { return magnitude.bit(index, work); };
    const auto top = magnitude.width - 1;
    const auto scale = static_cast<std::int64_t>(magnitude.width) -
                       static_cast<std::int64_t>(value.denominator.width);
    constexpr std::uint64_t significand_bits = 53;
    const auto cut = magnitude.width > significand_bits ? magnitude.width - significand_bits : 0;
    const auto discarded = magnitude.lowest < cut;
    if (discarded && (!two_word || top - magnitude.lowest > 1074))
        return std::nullopt;
    auto high = significant(cut, magnitude.width, bit);
    const bool carry =
        discarded && bit(cut - 1) && (magnitude.lowest < cut - 1 || ((high & 1U) != 0));
    if (carry)
        ++high;
    const auto encode = [&](std::uint64_t integer, std::uint64_t exponent, bool negative) {
        const si::Integer signed_word = negative ? -si::Integer{integer} : si::Integer{integer};
        return si::word({signed_word,
                         static_cast<std::int64_t>(exponent) - static_cast<std::int64_t>(top)},
                        false, true, work)
            .word;
    };
    const auto h = encode(high, cut, value.sign < 0);
    if (!two_word)
        return ScaledBinary64{h, 0.0, scale};
    if (!discarded)
        return ScaledTwoWord{h, 0.0, 0.0, scale};

    // With a high-word carry, the magnitude of the negative residual is the
    // two's complement of the discarded bits in their original grid.
    const auto residual = [&](auto index) {
        if (carry && index == magnitude.lowest)
            return true;
        return bit(index) != carry;
    };
    auto end = cut;
    while (end > magnitude.lowest && !residual(end - 1))
        --end;
    if (end - magnitude.lowest > significand_bits)
        return std::nullopt;
    const auto tail = significant(magnitude.lowest, end, residual);
    return ScaledTwoWord{h, encode(tail, magnitude.lowest, (value.sign < 0) != carry), 0.0, scale};
}

si::Interval exact_word(double word, si::Work& work) {
    work.step();
    const auto wire = std::bit_cast<std::uint64_t>(word);
    const auto exponent = (wire >> 52U) & 2047U;
    auto mantissa = wire & ((std::uint64_t{1} << 52U) - 1);
    if (exponent == 2047U)
        throw si::ArithmeticFailure{PrecisionUnresolved};
    if (exponent != 0)
        mantissa |= std::uint64_t{1} << 52U;
    const si::Integer signed_mantissa =
        (wire >> 63U) != 0 ? -si::Integer{mantissa} : si::Integer{mantissa};
    const si::Bound result{signed_mantissa,
                           exponent == 0 ? -1074 : static_cast<std::int64_t>(exponent) - 1075};
    return {result, result};
}

void qualify_candidate(SiEnclosure& center, const si::Interval& truth, const si::Bound& budget,
                       si::Work& work) {
    std::visit(
        [&](auto& candidate) {
            const auto actual = [&] {
                if constexpr (std::is_same_v<std::remove_cvref_t<decltype(candidate)>,
                                             ScaledBinary64>)
                    return exact_word(candidate.center, work);
                else
                    return si::add(exact_word(candidate.high, work),
                                   exact_word(candidate.low, work), work);
            }();
            // Compute the distance to the retained words, including alignment of both words.
            const auto scaled_truth = si::relative_to(truth, candidate.exponent2);
            const auto residual = si::absolute(si::subtract(scaled_truth, actual, work));
            const auto outward = si::word(residual.upper, true, false, work);
            const auto error = si::scale({outward.exact, outward.exact}, candidate.exponent2);
            work.step();
            if (si::compare(error.upper, budget) > 0)
                throw si::ArithmeticFailure{PrecisionUnresolved};
            candidate.radius = outward.word;
        },
        center);
}
} // namespace

core::Result<SiRationalProjection>
project_rational_words(SiRationalWordsView value, const SiRationalProjectionBudget& requested,
                       core::Precision representation, SiRationalProjectionLimits limits,
                       const std::stop_token& stop, core::FailureDetails details) noexcept {
    return core::guard_entry(
        operation, InvalidInput, details, [&]() -> core::Result<SiRationalProjection> {
            try {
                si::Work work{0, limits.steps, stop};
                work.step();
                if (representation != core::Precision::Binary64 &&
                    representation != core::Precision::TwoWord)
                    throw si::ArithmeticFailure{PrecisionUnsupported};
                const auto scale =
                    inspect_fraction(requested.scale, limits.input_integer_bits, work);
                const auto a =
                    inspect_fraction(requested.absolute, limits.input_integer_bits, work);
                const auto r =
                    inspect_fraction(requested.relative, limits.input_integer_bits, work);
                if (scale.sign <= 0 || a.sign < 0 || r.sign < 0)
                    throw si::ArithmeticFailure{InvalidInput};
                const auto source = inspect_fraction(value, limits.input_integer_bits, work);
                const bool two = representation == core::Precision::TwoWord;
                if (const auto exact = exact_dyadic(source, two, work)) {
                    work.step();
                    return SiRationalProjection{*exact, work.used};
                }
                const auto truth = enclosure(source, work);
                const auto abs_budget =
                    si::multiply(enclosure(scale, work), enclosure(a, work), work);
                const auto rel_budget = si::multiply(enclosure(r, work), si::absolute(truth), work);
                const auto budget = si::add(abs_budget, rel_budget, work);
                auto candidate = si::center(truth, two, work).value;
                qualify_candidate(candidate, truth, budget.lower, work);
                return SiRationalProjection{candidate, work.used};
            } catch (const si::ArithmeticFailure& error) {
                return std::unexpected(core::Failure{error.code, operation, details});
            } catch (const std::overflow_error&) {
                return std::unexpected(core::Failure{PrecisionUnresolved, operation, details});
            }
        });
}
} // namespace uvlas::io::detail
