#pragma once

#include "aero/local_geometry_bounds_private.hpp"
#include "aero/p1_shared_kernel.hpp"
#include "case_uniform_flow_projection.hpp"

#include <bit>
#include <ranges>

namespace uvlas::coupling::detail::case_components {
namespace a = aero::detail;
namespace p = a::pchip;
namespace m = a::p1;
using enum core::FailureCode;

[[noreturn]] inline void reject(core::FailureCode code) {
    throw m::ArithmeticFailure{code};
}

inline std::uint64_t bits(std::span<const std::uint32_t> limbs) {
    if (limbs.empty())
        return 0;
    if (limbs.size() > 1024)
        reject(ResourceExceeded);
    if (limbs.back() == 0)
        reject(InvalidInput);
    return 32U * (limbs.size() - 1U) + static_cast<std::uint64_t>(std::bit_width(limbs.back()));
}

inline void validate_words(p::RationalView value, std::uint64_t input_bits) {
    if (value.sign < -1 || value.sign > 1 || value.denominator.empty() ||
        value.numerator.empty() != (value.sign == 0) ||
        (value.sign == 0 && (value.denominator.size() != 1 || value.denominator.front() != 1)))
        reject(InvalidInput);
    if (bits(value.numerator) > input_bits || bits(value.denominator) > input_bits)
        reject(ResourceExceeded);
}

inline m::Interval rational(p::RationalView value, std::uint64_t input_bits, m::Work& work) {
    work.step();
    validate_words(value, input_bits);
    const auto assemble = [&](std::span<const std::uint32_t> limbs) {
        a::ExactInteger result = 0;
        for (const auto limb : limbs | std::views::reverse) {
            work.step();
            result <<= 32U;
            result += limb;
        }
        return result;
    };
    if (value.sign == 0)
        return m::integer(0);
    const auto n = assemble(value.numerator);
    return m::divide(m::exact(value.sign < 0 ? -n : n, 0), m::exact(assemble(value.denominator), 0),
                     work);
}

inline io::SiEnclosure project(const m::Interval& value, const CaseFlowComponentBudget& budget,
                               core::Precision representation, std::uint64_t input_bits,
                               a::LocalRootScratch& scratch, m::Work& work) {
    const auto scale = rational(budget.scale, input_bits, work);
    const auto absolute = rational(budget.absolute, input_bits, work);
    const auto relative = rational(budget.relative, input_bits, work);
    if (m::compare(scale.lower, {}) <= 0 || m::compare(absolute.lower, {}) < 0 ||
        m::compare(relative.lower, {}) < 0)
        reject(InvalidInput);
    const auto tolerance = m::add(m::multiply(scale, absolute, work),
                                  m::multiply(relative, m::absolute(value), work), work);
    work.step();
    const auto result = a::enclose_p1_value(value, representation, scratch);
    const auto error = m::scale(a::decode_p1_word(result.radius()), result.exponent());
    if (m::compare(error.upper, tolerance.lower) > 0)
        reject(PrecisionUnresolved);
    const auto center = result.center();
    if (const auto* single = std::get_if<a::SingleWord>(&center))
        return io::ScaledBinary64{single->value, result.radius(), result.exponent()};
    if (const auto* two = std::get_if<a::DoubleWord>(&center))
        return io::ScaledTwoWord{two->high, two->low, result.radius(), result.exponent()};
    reject(PrecisionUnsupported);
}

} // namespace uvlas::coupling::detail::case_components
