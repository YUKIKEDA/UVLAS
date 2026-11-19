#pragma once

#include "../structure/exact_polynomial_detail.hpp"
#include "si_case_detail.hpp"

namespace uvlas::io::detail {
struct SiRationalUnit {
    constexpr SiRationalUnit(std::string_view spelling, std::string_view coefficient,
                             std::int64_t scale) noexcept
        : token(spelling), digits(coefficient), exponent(scale) {}
    std::string_view token, digits;
    std::int64_t exponent;
};
inline std::array<SiRationalUnit, 3> si_rational_units(SiUnits units) {
    constexpr std::array<SiRationalUnit, 6> length{{{"1", "1", 0},
                                                    {"0.001", "0001", -3},
                                                    {"0.01", "001", -2},
                                                    {"1000", "1000", 0},
                                                    {"0.0254", "00254", -4},
                                                    {"0.3048", "03048", -4}}};
    constexpr std::array<SiRationalUnit, 4> mass{{{"1", "1", 0},
                                                  {"0.001", "0001", -3},
                                                  {"1000", "1000", 0},
                                                  {"0.45359237", "045359237", -8}}};
    constexpr std::array<SiRationalUnit, 3> time{
        {{"1", "1", 0}, {"0.001", "0001", -3}, {"60", "60", 0}}};
    return {length[static_cast<std::size_t>(units.length)],
            mass[static_cast<std::size_t>(units.mass)], time[static_cast<std::size_t>(units.time)]};
}
namespace rational_plan {
[[noreturn]] inline void refuse() {
    throw core::Failure{core::FailureCode::ResourceExceeded, core::Operation::PrepareModel};
}
inline std::uint64_t add(std::uint64_t a, std::uint64_t b) {
    if (b > std::numeric_limits<std::uint64_t>::max() - a)
        refuse();
    return a + b;
}
inline std::uint64_t multiply(std::uint64_t a, std::uint64_t b) {
    if (a && b > std::numeric_limits<std::uint64_t>::max() / a)
        refuse();
    return a * b;
}
struct Admission {
    const structure::ExactLimits& limits;
    bool exact;
    bool certified = true;
    std::uint64_t storages = 0;
    void stage(std::uint64_t work, std::uint64_t borrowed) {
        storages = add(storages, 1);
        const auto result = structure::detail::exact_rational_admission(limits, work, borrowed);
        if (!result) {
            if (result.error().code != core::FailureCode::ResourceExceeded)
                throw result.error();
            certified = false;
            if (exact)
                refuse();
        }
    }
    void value(std::uint64_t numerator, std::uint64_t denominator) {
        if (std::max(numerator, denominator) > limits.maximum_value_bits) {
            certified = false;
            if (exact)
                refuse();
        }
    }
    [[nodiscard]] std::uint64_t retained(std::uint64_t numerator, std::uint64_t denominator,
                                         std::uint64_t work) const {
        using Backend = structure::detail::IntegerBackend;
        constexpr auto width = Backend::limb_bits;
        constexpr auto inline_bits = width * Backend::internal_limb_count;
        // The same two guard limbs and fourfold growth envelope as backend admit.
        // Normalization can reduce size without releasing the pre-normalization capacity.
        const auto limbs = add(add(work / width, work % width != 0 ? 1 : 0), 2);
        const auto capacity =
            multiply(multiply(4, limbs), sizeof(boost::multiprecision::limb_type));
        return multiply(static_cast<std::uint64_t>(numerator > inline_bits) +
                            static_cast<std::uint64_t>(denominator > inline_bits),
                        capacity);
    }
    [[nodiscard]] std::uint64_t parse_work(SiDecimalView value, std::uint64_t characters) const {
        if (characters > limits.maximum_input_characters)
            refuse();
        return add(multiply(4, add(value.digits.size(), si::magnitude(value.exponent10))), 1);
    }
};
struct Extent {
    std::uint64_t numerator, denominator, borrowed;
};
struct Envelope : Admission {
    explicit Envelope(const structure::ExactLimits& configured) : Admission{configured, false} {}
    Extent parse(SiDecimalView input, std::uint64_t characters, std::uint64_t extra = 0) {
        const auto work = parse_work(input, characters);
        stage(work, extra);
        const auto coefficient = multiply(4, input.digits.size());
        const auto power = multiply(4, si::magnitude(input.exponent10));
        const auto n = input.exponent10 < 0 ? coefficient : add(coefficient, power);
        const auto d = input.exponent10 < 0 ? add(power, 1) : 1;
        value(n, d);
        return {n, d, retained(n, d, work)};
    }
    Extent binary(const Extent& a, const Extent& b, bool divide, std::uint64_t extra = 0) {
        const auto n = add(a.numerator, divide ? b.denominator : b.numerator);
        const auto d = add(a.denominator, divide ? b.numerator : b.denominator);
        const auto work = std::max({n, d, a.numerator, a.denominator, b.numerator, b.denominator});
        stage(work, add(extra, add(a.borrowed, b.borrowed)));
        value(n, d);
        return {n, d, retained(n, d, work)};
    }
};
using Fixed = si::Integer;
inline std::uint64_t bits(const Fixed& value) {
    return si::bits(value);
}
inline Fixed product(const Fixed& a, const Fixed& b) {
    if (a != 0 && si::abs_integer(b) > std::numeric_limits<Fixed>::max() / si::abs_integer(a))
        refuse();
    return a * b;
}
inline Fixed power(Fixed base, std::uint64_t exponent) {
    Fixed result = 1;
    while (exponent) {
        if (exponent & 1U)
            result = product(result, base);
        exponent >>= 1;
        if (exponent)
            base = product(base, base);
    }
    return result;
}
struct Fraction {
    Fraction(Fixed n, Fixed d, std::uint64_t retained_bytes)
        : numerator(std::move(n)), denominator(std::move(d)), borrowed(retained_bytes) {}
    Fixed numerator, denominator;
    std::uint64_t borrowed;
};
struct Certificate : Admission {
    explicit Certificate(const structure::ExactLimits& configured) : Admission{configured, true} {}
    Fraction normalized(Fixed n, Fixed d, std::uint64_t work) {
        const auto raw_n = bits(n), raw_d = bits(d);
        Fixed a = si::abs_integer(n), b = d;
        while (b != 0) {
            Fixed remainder = a % b;
            a = std::move(b);
            b = std::move(remainder);
        }
        n /= a;
        d /= a;
        value(bits(n), bits(d));
        return {std::move(n), std::move(d), retained(raw_n, raw_d, work)};
    }
    Fraction parse(SiDecimalView input, std::uint64_t characters, std::uint64_t extra = 0) {
        const auto work = parse_work(input, characters);
        stage(work, extra);
        Fixed coefficient = 0;
        for (const char c : input.digits) {
            const auto digit = static_cast<unsigned>(c - '0');
            if (coefficient > (std::numeric_limits<Fixed>::max() - digit) / 10)
                refuse();
            coefficient *= 10;
            coefficient += digit;
        }
        auto ten = power(Fixed{10}, si::magnitude(input.exponent10));
        Fixed denominator = 1;
        if (input.exponent10 < 0)
            denominator = std::move(ten);
        else
            coefficient = product(coefficient, ten);
        if (input.sign < 0)
            coefficient = -coefficient;
        return normalized(std::move(coefficient), std::move(denominator), work);
    }
    Fraction binary(const Fraction& a, const Fraction& b, bool divide, std::uint64_t extra = 0) {
        const auto an = bits(a.numerator), ad = bits(a.denominator);
        const auto bn = bits(b.numerator), bd = bits(b.denominator);
        const auto work =
            std::max({add(an, divide ? bd : bn), add(ad, divide ? bn : bd), an, ad, bn, bd});
        stage(work, add(extra, add(a.borrowed, b.borrowed)));
        auto n = product(a.numerator, divide ? b.denominator : b.numerator);
        auto d = product(a.denominator, divide ? b.numerator : b.denominator);
        if (d < 0) {
            n = -n;
            d = -d;
        }
        return normalized(std::move(n), std::move(d), work);
    }
};
// Large automatic payloads on the optional exact boundary-certificate path.
// Names follow walk -> power/parse -> binary -> normalize/product call nesting.
struct WalkFrame {
    Fraction source_n, source_d, source, result, unit, powered, returned;
};
struct PowerFrame {
    Fraction parameter, result, seed, product, square;
};
struct ParseFrame {
    Fixed coefficient, ten, denominator, power_base, power_result, power_product, returned;
};
struct NormalizeFrame {
    Fixed parameter_n, parameter_d, a, b, remainder, divided_n, divided_d;
};
struct ProductFrame {
    Fixed parameter_a, parameter_b, maximum, absolute_a, absolute_b, guard_quotient, returned;
};
struct BackendFrame {
    Fixed division_remainder, division_product;
};
inline constexpr std::uint64_t certificate_stack_bytes =
    sizeof(WalkFrame) + sizeof(PowerFrame) + sizeof(ParseFrame) + sizeof(NormalizeFrame) +
    sizeof(ProductFrame) + sizeof(BackendFrame);
template <class Model, class Value>
Value rational_power(Model& model, Value base, std::int64_t exponent, std::uint64_t outer) {
    auto result = model.parse(SiDecimalView{1, "1", 0}, 1, add(outer, base.borrowed));
    if (exponent < 0)
        base = model.binary(result, base, true, outer);
    auto remaining = si::magnitude(exponent);
    while (remaining) {
        if (remaining & 1U)
            result = model.binary(result, base, false, outer);
        remaining >>= 1;
        if (remaining)
            base = model.binary(base, base, false, add(outer, result.borrowed));
    }
    return result;
}
template <class Model>
void walk(Model& model, const SiRecipeView& recipe) {
    const auto parse = [&](SiDecimalView value, std::uint64_t extra = 0) {
        const auto size = si_decimal_token_size(value);
        if (!size)
            throw size.error();
        return model.parse(value, *size, extra);
    };
    const auto source = [&] {
        if (const auto* value = std::get_if<SiDecimalView>(&recipe.value))
            return parse(*value);
        const auto& value = std::get<SiRationalView>(recipe.value);
        const auto n = parse(value.numerator);
        const auto d = parse(value.denominator, n.borrowed);
        return model.binary(n, d, true);
    };
    auto result = source();
    const auto factors = si_rational_units(recipe.units);
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const auto exponent = recipe.dimensions.exponents[axis];
        if (!exponent)
            continue;
        const auto& unit = factors[axis];
        auto base = model.parse(SiDecimalView{1, unit.digits, unit.exponent}, unit.token.size(),
                                result.borrowed);
        auto powered = rational_power(model, std::move(base), exponent, result.borrowed);
        result = model.binary(result, powered, false);
    }
}
inline void admit(const SiRecipeView& recipe, const structure::ExactLimits& limits,
                  const core::ResourcePlan& resources) {
    Envelope envelope{limits};
    walk(envelope, recipe);
    std::uint64_t characters = 0;
    const auto count = [&](SiDecimalView value) {
        const auto size = si_decimal_token_size(value);
        if (!size)
            throw size.error();
        if (*size > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
            refuse();
        characters = add(characters, *size);
    };
    if (const auto* value = std::get_if<SiDecimalView>(&recipe.value))
        count(*value);
    else {
        const auto& rational = std::get<SiRationalView>(recipe.value);
        count(rational.numerator);
        count(rational.denominator);
    }
    const auto payload =
        add(characters, multiply(envelope.storages, sizeof(structure::detail::RationalStorage)));
    if (payload > resources.remaining_bytes())
        refuse();
    if (!envelope.certified) {
        // Tight caps are proved with allocator-free, checked 1024-bit integers.
        // Larger ambiguous witnesses remain a declined conservative admission.
        Certificate certificate{limits};
        walk(certificate, recipe);
    }
}
} // namespace rational_plan
} // namespace uvlas::io::detail
