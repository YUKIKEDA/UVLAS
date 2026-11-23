#include "si_case_detail.hpp"
#include "si_exact_words.hpp"
#include "si_rational_plan.hpp"

#include <charconv>

namespace uvlas::io::detail {
core::Result<std::uint64_t> si_decimal_token_size(SiDecimalView value) noexcept {
    if (value.sign == 0)
        return 0;
    std::array<char, 24> exponent{};
    const auto encoded =
        std::to_chars(exponent.data(), exponent.data() + exponent.size(), value.exponent10);
    const auto suffix =
        value.exponent10 ? static_cast<std::uint64_t>(encoded.ptr - exponent.data()) + 1 : 0;
    const auto extra = suffix + static_cast<std::uint64_t>(value.sign < 0);
    if (value.digits.size() > std::numeric_limits<std::uint64_t>::max() - extra)
        return std::unexpected(
            core::Failure{core::FailureCode::ResourceExceeded, core::Operation::PrepareModel});
    return value.digits.size() + extra;
}
namespace {
si::Interval token_interval(SiContext& c, std::string_view token, bool integer_only) {
    c.require(SiQuota::RequestTokenBytes, token.size());
    if (token.empty())
        c.reject(SiReason::InvalidExactInput);
    std::size_t position = 0;
    int sign = 1;
    if (token[position] == '+' || token[position] == '-') {
        sign = token[position] == '-' ? -1 : 1;
        if (++position == token.size())
            c.reject(SiReason::InvalidExactInput);
    }
    auto value = si::integer(0);
    bool dot = false;
    std::uint64_t digits = 0, fractional = 0, trailing = 0;
    for (; position < token.size() && token[position] != 'e' && token[position] != 'E';
         ++position) {
        c.work.step();
        const char ch = token[position];
        if (ch == '.' && !integer_only && !dot) {
            dot = true;
            continue;
        }
        if (ch < '0' || ch > '9')
            c.reject(SiReason::InvalidExactInput);
        digits = c.add(digits, 1);
        trailing = ch == '0' ? c.add(trailing, 1) : 0;
        if (dot)
            fractional = c.add(fractional, 1);
        const std::array operands{value, si::integer(ch - '0')};
        value = si::execute(si::digit_program, operands, c.program.values()[0], c.work)[0];
    }
    if (!digits ||
        fractional > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
        c.reject(SiReason::InvalidExactInput);
    c.require(SiQuota::RequestDigits, digits);
    const bool zero = value.lower.mantissa == 0 && value.upper.mantissa == 0;
    std::int64_t exponent = 0;
    if (position < token.size()) {
        if (integer_only)
            c.reject(SiReason::InvalidExactInput);
        ++position;
        bool negative = false;
        if (position < token.size() && (token[position] == '+' || token[position] == '-')) {
            negative = token[position] == '-';
            ++position;
        }
        if (position == token.size())
            c.reject(SiReason::InvalidExactInput);
        std::uint64_t magnitude = 0;
        const auto cap = negative
                             ? std::uint64_t{1} << 63
                             : static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
        for (; position < token.size(); ++position) {
            c.work.step();
            const char ch = token[position];
            if (ch < '0' || ch > '9')
                c.reject(SiReason::InvalidExactInput);
            const auto digit = static_cast<unsigned>(ch - '0');
            if (!zero && (magnitude > cap / 10 || (magnitude == cap / 10 && digit > cap % 10)))
                c.reject(SiReason::Overflow, core::FailureCode::ResourceExceeded,
                         SiQuota::ExponentMagnitude);
            if (!zero)
                magnitude = magnitude * 10 + digit;
        }
        exponent = negative ? (magnitude == (std::uint64_t{1} << 63)
                                   ? std::numeric_limits<std::int64_t>::min()
                                   : -static_cast<std::int64_t>(magnitude))
                            : static_cast<std::int64_t>(magnitude);
    }
    if (zero)
        return value;
    exponent = si::subtract_exponent(exponent, static_cast<std::int64_t>(fractional));
    if (trailing > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
        c.reject(SiReason::Overflow, core::FailureCode::ResourceExceeded);
    const auto normalized = si::add_exponent(exponent, static_cast<std::int64_t>(trailing));
    c.require(SiQuota::ExponentMagnitude, si::magnitude(normalized));
    if (sign < 0)
        value = {si::negate(value.upper), si::negate(value.lower)};
    return si::multiply(
        value, si::programmed_power(si::integer(10), exponent, c.program.values()[0], c.work),
        c.work);
}
} // namespace
si::Interval si_exact_interval(SiContext& c, const SiExactInput& input) {
    if (const auto* decimal = std::get_if<SiDecimalInput>(&input))
        return token_interval(c, decimal->token, false);
    const auto& rational = std::get<SiRationalInput>(input);
    const auto numerator = token_interval(c, rational.numerator, true);
    const auto denominator = token_interval(c, rational.denominator, true);
    if (si::compare(denominator.lower, {}) <= 0)
        c.reject(SiReason::InvalidExactInput);
    return si::divide(numerator, denominator, c.work);
}
si::Interval si_unit_interval(SiContext& c, const SiRecipeView& recipe, const si::Interval& pi) {
    constexpr std::array<std::int64_t, 6> length_n{1, 1, 1, 1000, 127, 381};
    constexpr std::array<std::int64_t, 6> length_d{1, 1000, 100, 1, 5000, 1250};
    constexpr std::array<std::int64_t, 4> mass_n{1, 1, 1000, 45359237};
    constexpr std::array<std::int64_t, 4> mass_d{1, 1000, 1, 100000000};
    constexpr std::array<std::int64_t, 3> time_n{1, 1, 60};
    constexpr std::array<std::int64_t, 3> time_d{1, 1000, 1};
    const auto l = static_cast<std::size_t>(recipe.units.length);
    const auto m = static_cast<std::size_t>(recipe.units.mass);
    const auto t = static_cast<std::size_t>(recipe.units.time);
    if (l >= length_n.size() || m >= mass_n.size() || t >= time_n.size() ||
        static_cast<unsigned>(recipe.units.angle) > static_cast<unsigned>(SiAngleUnit::Degree))
        c.reject(SiReason::InvalidExactInput);
    auto result = si::integer(1);
    const std::array factors{si::divide(si::integer(length_n[l]), si::integer(length_d[l]), c.work),
                             si::divide(si::integer(mass_n[m]), si::integer(mass_d[m]), c.work),
                             si::divide(si::integer(time_n[t]), si::integer(time_d[t]), c.work),
                             recipe.units.angle == SiAngleUnit::Degree
                                 ? si::divide(si::integer(1), si::integer(180), c.work)
                                 : si::integer(1)};
    for (std::size_t i = 0; i < 4; ++i) {
        const auto exponent = static_cast<std::int64_t>(recipe.dimensions.exponents[i]);
        c.require(SiQuota::UnitPower, si::magnitude(exponent));
        if (exponent)
            result = si::multiply(
                result, si::programmed_power(factors[i], exponent, c.program.values()[0], c.work),
                c.work);
    }
    if (recipe.pi_power)
        result = si::multiply(
            result, si::programmed_power(pi, recipe.pi_power, c.program.values()[0], c.work),
            c.work);
    return result;
}
} // namespace uvlas::io::detail

namespace uvlas::io {
namespace {
using structure::ExactRational;
core::Result<ExactRational> decimal_rational(SiDecimalView value,
                                             const structure::ExactLimits& limits,
                                             core::ResourcePlan& plan) {
    if (value.sign == 0)
        return ExactRational{};
    const auto required = detail::si_decimal_token_size(value);
    if (!required)
        return std::unexpected(required.error());
    const auto capacity = *required;
    if (capacity > static_cast<std::size_t>(std::numeric_limits<core::Index>::max()) ||
        capacity > limits.maximum_input_characters)
        return std::unexpected(
            core::Failure{core::FailureCode::ResourceExceeded, core::Operation::PrepareModel});
    auto text =
        plan.allocate<char>(static_cast<core::Index>(capacity), core::Operation::PrepareModel);
    if (!text)
        return std::unexpected(text.error());
    auto* begin = text->values().data();
    auto* next = begin;
    if (value.sign < 0)
        *next++ = '-';
    next = std::copy(value.digits.begin(), value.digits.end(), next);
    if (value.exponent10) {
        *next++ = 'e';
        const auto exponent = std::to_chars(next, begin + capacity, value.exponent10);
        if (exponent.ec != std::errc{})
            return std::unexpected(
                core::Failure{core::FailureCode::ResourceExceeded, core::Operation::PrepareModel});
        next = exponent.ptr;
    }
    return ExactRational::from_decimal({begin, static_cast<std::size_t>(next - begin)}, limits,
                                       plan);
}
core::Result<ExactRational> rational_power(ExactRational factor, std::int64_t exponent,
                                           const structure::ExactLimits& limits,
                                           core::ResourcePlan& plan) {
    auto result = ExactRational::from_decimal("1", limits, plan);
    if (!result)
        return std::unexpected(result.error());
    auto remaining = detail::si::magnitude(exponent);
    if (exponent < 0) {
        auto inverse = result->divide(factor, limits, plan);
        if (!inverse)
            return std::unexpected(inverse.error());
        factor = std::move(*inverse);
    }
    while (remaining) {
        if (remaining & 1U) {
            auto next = result->multiply(factor, limits, plan);
            if (!next)
                return std::unexpected(next.error());
            result = std::move(next);
        }
        remaining >>= 1;
        if (remaining) {
            auto next = factor.multiply(factor, limits, plan);
            if (!next)
                return std::unexpected(next.error());
            factor = std::move(*next);
        }
    }
    return result;
}
} // namespace
core::Result<ExactRational> materialize_rational(const SiRecipeView& recipe,
                                                 const structure::ExactLimits& limits,
                                                 core::ResourcePlan& plan) noexcept {
    return core::guard_entry(
        core::Operation::PrepareModel, core::FailureCode::InvalidInput, {},
        [&]() -> core::Result<ExactRational> {
            const auto check_decimal = [&](SiDecimalView value) {
                if (value.digits.size() > limits.maximum_input_characters)
                    throw core::Failure{core::FailureCode::ResourceExceeded,
                                        core::Operation::PrepareModel};
                if (value.digits.empty() || value.sign < -1 || value.sign > 1)
                    throw core::Failure{core::FailureCode::InvalidInput,
                                        core::Operation::PrepareModel};
                bool nonzero = false;
                for (const char digit : value.digits) {
                    if (digit < '0' || digit > '9')
                        throw core::Failure{core::FailureCode::InvalidInput,
                                            core::Operation::PrepareModel};
                    nonzero = nonzero || digit != '0';
                }
                if ((value.sign != 0) != nonzero)
                    throw core::Failure{core::FailureCode::InvalidInput,
                                        core::Operation::PrepareModel};
            };
            if (static_cast<unsigned>(recipe.units.length) >
                    static_cast<unsigned>(SiLengthUnit::Foot) ||
                static_cast<unsigned>(recipe.units.mass) >
                    static_cast<unsigned>(SiMassUnit::PoundMass) ||
                static_cast<unsigned>(recipe.units.time) >
                    static_cast<unsigned>(SiTimeUnit::Minute) ||
                static_cast<unsigned>(recipe.units.angle) >
                    static_cast<unsigned>(SiAngleUnit::Degree) ||
                recipe.pi_power != (recipe.units.angle == SiAngleUnit::Degree
                                        ? recipe.dimensions.exponents[3]
                                        : 0))
                return std::unexpected(
                    core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
            if (const auto* decimal = std::get_if<SiDecimalView>(&recipe.value))
                check_decimal(*decimal);
            else {
                const auto& r = std::get<SiRationalView>(recipe.value);
                check_decimal(r.numerator);
                check_decimal(r.denominator);
                if (r.denominator.sign <= 0 || r.numerator.exponent10 < 0 ||
                    r.denominator.exponent10 < 0)
                    return std::unexpected(core::Failure{core::FailureCode::InvalidInput,
                                                         core::Operation::PrepareModel});
            }
            const auto zero = std::visit(
                [](const auto& x) {
                    if constexpr (std::is_same_v<std::remove_cvref_t<decltype(x)>, SiDecimalView>)
                        return x.sign == 0;
                    else
                        return x.numerator.sign == 0;
                },
                recipe.value);
            if (zero)
                return ExactRational{};
            if (recipe.pi_power)
                return std::unexpected(core::Failure{core::FailureCode::PrecisionUnsupported,
                                                     core::Operation::PrepareModel});
            try {
                detail::rational_plan::admit(recipe, limits, plan);
            } catch (const core::Failure& failure) {
                return std::unexpected(failure);
            }
            core::Result<ExactRational> result = ExactRational{};
            if (const auto* decimal = std::get_if<SiDecimalView>(&recipe.value))
                result = decimal_rational(*decimal, limits, plan);
            else {
                const auto& r = std::get<SiRationalView>(recipe.value);
                auto n = decimal_rational(r.numerator, limits, plan);
                if (!n)
                    return std::unexpected(n.error());
                auto d = decimal_rational(r.denominator, limits, plan);
                if (!d)
                    return std::unexpected(d.error());
                result = n->divide(*d, limits, plan);
            }
            if (!result)
                return std::unexpected(result.error());
            const auto factors = detail::si_rational_units(recipe.units);
            for (std::size_t i = 0; i < 3; ++i) {
                if (!recipe.dimensions.exponents[i])
                    continue;
                auto base = ExactRational::from_decimal(factors[i].token, limits, plan);
                if (!base)
                    return std::unexpected(base.error());
                auto power =
                    rational_power(std::move(*base), recipe.dimensions.exponents[i], limits, plan);
                if (!power)
                    return std::unexpected(power.error());
                result = result->multiply(*power, limits, plan);
                if (!result)
                    return std::unexpected(result.error());
            }
            return result;
        });
}
} // namespace uvlas::io

namespace uvlas::io::detail {
namespace {
constexpr auto word_operation = core::Operation::PrepareModel;
constexpr std::array<std::uint32_t, 1> one_word{1};

constexpr std::uint64_t size_in_words(std::uint64_t bits) noexcept {
    return bits / 32 + (bits % 32 != 0 ? 1U : 0U);
}

bool copy_words(const structure::detail::IntegerBackend& magnitude, std::span<std::uint32_t> target,
                const std::stop_token& stop) noexcept {
    constexpr auto limb_bits = structure::detail::IntegerBackend::limb_bits;
    static_assert(limb_bits % 32 == 0);
    constexpr auto words_per_limb = limb_bits / 32;
    for (std::size_t i = 0; i < target.size(); ++i) {
        if (stop.stop_requested())
            return false;
        target[i] = static_cast<std::uint32_t>(magnitude.limbs()[i / words_per_limb] >>
                                               (32 * (i % words_per_limb)));
    }
    return true;
}
} // namespace

SiRationalWords::SiRationalWords(core::PlannedBuffer<std::uint32_t> words, std::size_t numerator,
                                 int sign) noexcept
    : words_(std::move(words)), numerator_(numerator), sign_(sign) {}

SiRationalWordsView SiRationalWords::view() const& noexcept {
    const auto data = words_.values();
    if (data.empty())
        return {0, {}, one_word};
    return {sign_, data.first(numerator_), data.subspan(numerator_)};
}

core::Result<SiRationalWords> materialize_rational_words(const SiRecipeView& recipe,
                                                         const structure::ExactLimits& limits,
                                                         core::ResourcePlan& memory,
                                                         const std::stop_token& stop) noexcept {
    return core::guard_entry(
        word_operation, core::FailureCode::InvalidInput, {},
        [&]() -> core::Result<SiRationalWords> {
            const auto cancelled = core::Failure{core::FailureCode::Cancelled, word_operation};
            if (stop.stop_requested())
                return std::unexpected(cancelled);
            auto exact = materialize_rational(recipe, limits, memory);
            if (!exact)
                return std::unexpected(exact.error());
            if (stop.stop_requested())
                return std::unexpected(cancelled);
            const auto sign = exact->sign();
            if (sign == 0)
                return SiRationalWords{};
            const auto numerator = size_in_words(exact->numerator_bits());
            const auto denominator = size_in_words(exact->denominator_bits());
            constexpr auto index_max =
                static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
            if (numerator > index_max || denominator > index_max - numerator)
                return std::unexpected(
                    core::Failure{core::FailureCode::ResourceExceeded, word_operation});
            auto words = memory.allocate<std::uint32_t>(
                static_cast<core::Index>(numerator + denominator), word_operation);
            if (!words)
                return std::unexpected(words.error());
            const auto& value = structure::detail::ExactFactory::value(*exact).backend();
            const auto storage = words->values();
            const auto split = static_cast<std::size_t>(numerator);
            if (!copy_words(value.num(), storage.first(split), stop) ||
                !copy_words(value.denom(), storage.subspan(split), stop) || stop.stop_requested())
                return std::unexpected(cancelled);
            return SiRationalWords{std::move(*words), split, sign};
        });
}
} // namespace uvlas::io::detail
