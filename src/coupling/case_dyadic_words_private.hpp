#pragma once

#include "aero/certified_trigonometry_private.hpp"
#include "aero/pchip_curve.hpp"

namespace uvlas::coupling::detail {
struct DyadicWords {
    std::array<std::uint32_t, 128> numerator{}, denominator{};
    std::size_t numerator_size = 0, denominator_size = 1;
    int sign = 0;

    DyadicWords() = default;

    DyadicWords(aero::detail::ExactInteger number, unsigned shift, aero::detail::p1::Work& work) {
        work.step();
        denominator[0] = 1;
        if (number == 0)
            return;
        sign = number < 0 ? -1 : 1;
        if (sign < 0)
            number = -number;
        const auto trim =
            std::min(shift, static_cast<unsigned>(boost::multiprecision::lsb(number)));
        number >>= trim;
        shift -= trim;
        const auto word = shift / 32U;
        if (word >= denominator.size())
            throw aero::detail::p1::ArithmeticFailure{core::FailureCode::ResourceExceeded};
        denominator[0] = 0;
        denominator[word] = std::uint32_t{1} << (shift % 32U);
        denominator_size = static_cast<std::size_t>(word) + 1;
        while (number != 0) {
            work.step();
            if (numerator_size == numerator.size())
                throw aero::detail::p1::ArithmeticFailure{core::FailureCode::ResourceExceeded};
            numerator[numerator_size++] = (number & 0xffffffffU).convert_to<std::uint32_t>();
            number >>= 32U;
        }
    }

    [[nodiscard]] aero::detail::pchip::RationalView view() const noexcept {
        return {sign, std::span{numerator}.first(numerator_size),
                std::span{denominator}.first(denominator_size)};
    }
};

} // namespace uvlas::coupling::detail
