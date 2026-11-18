#pragma once
#include "p1_interval_private.hpp"

namespace uvlas::aero::detail::p1 {
/// Fixed payload for exact source arithmetic. No Boost Karatsuba arena is allocated.
struct ExactProductScratch {
    std::array<std::uint32_t, 1024> left{}, right{}, output{};
};
/// Multiply fixed integers using 32-bit digits after removing exact powers of two.
/// Each inner product plus digit and carry fits exactly in uint64_t.
/// (2^32-1)^2 + (2^32-1) + (2^32-1) = 2^64-1.
inline ExactInteger exact_product(const ExactInteger& left, const ExactInteger& right,
                                  ExactProductScratch& scratch) {
    if (left == 0 || right == 0)
        return 0;
    const bool negative = (left < 0) != (right < 0);
    ExactInteger a = left < 0 ? -left : left, b = right < 0 ? -right : right;
    const auto ap = boost::multiprecision::lsb(a), bp = boost::multiprecision::lsb(b);
    a >>= ap;
    b >>= bp;
    const auto na = static_cast<std::size_t>(boost::multiprecision::msb(a) / 32 + 1);
    const auto nb = static_cast<std::size_t>(boost::multiprecision::msb(b) / 32 + 1);
    scratch.left.fill(0);
    scratch.right.fill(0);
    scratch.output.fill(0);
    boost::multiprecision::export_bits(a, scratch.left.begin(), 32, false);
    boost::multiprecision::export_bits(b, scratch.right.begin(), 32, false);
    for (std::size_t i = 0; i < na; ++i) {
        std::uint64_t carry = 0;
        for (std::size_t j = 0; j < nb; ++j) {
            if (i + j >= scratch.output.size())
                throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
            const std::uint64_t value =
                static_cast<std::uint64_t>(scratch.left[i]) * scratch.right[j] +
                scratch.output[i + j] + carry;
            scratch.output[i + j] = static_cast<std::uint32_t>(value);
            carry = value >> 32U;
        }
        std::size_t index = i + nb;
        while (carry != 0) {
            if (index >= scratch.output.size())
                throw ArithmeticFailure{core::FailureCode::PrecisionUnresolved};
            const std::uint64_t value = static_cast<std::uint64_t>(scratch.output[index]) + carry;
            scratch.output[index++] = static_cast<std::uint32_t>(value);
            carry = value >> 32U;
        }
    }
    const std::size_t count = std::min(na + nb, scratch.output.size());
    ExactInteger result = 0;
    for (std::size_t i = count; i > 0; --i) {
        result <<= 32U;
        result += scratch.output[i - 1];
    }
    result <<= ap + bp;
    return negative ? -result : result;
}
} // namespace uvlas::aero::detail::p1
