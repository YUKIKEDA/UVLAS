#pragma once
#include "local_geometry_packet.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>

namespace uvlas::aero::detail {
/// Optional, allocation-free clock; the target is checked at the use site.
class LocalClock {
public:
    explicit LocalClock(core::StageTiming* target) noexcept : target_(target) {
        if (target_ != nullptr)
            start_ = std::chrono::steady_clock::now();
    }
    ~LocalClock() {
        if (target_ != nullptr && start_) {
            const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                     std::chrono::steady_clock::now() - *start_)
                                     .count();
            ++target_->calls;
            target_->elapsed_nanoseconds += static_cast<std::uint64_t>(elapsed);
        }
    }
    LocalClock(const LocalClock&) = delete;
    LocalClock& operator=(const LocalClock&) = delete;

private:
    core::StageTiming* target_;
    std::optional<std::chrono::steady_clock::time_point> start_;
};
struct LocalRootScratch {
    std::array<ExactInteger, 23> integers;
};
struct LocalWork {
    std::uint64_t used, maximum;
    std::stop_token stop;
    [[nodiscard]] core::Result<void> step() noexcept {
        if (stop.stop_requested())
            return std::unexpected(
                core::Failure{core::FailureCode::Cancelled, core::Operation::EvaluateAero});
        if (used >= maximum)
            return std::unexpected(
                core::Failure{core::FailureCode::ResourceExceeded, core::Operation::EvaluateAero});
        ++used;
        return {};
    }
};
/// Positive integer comparison with separate exponents; no exponent-sized allocation.
inline int local_compare(const ExactInteger& a, int ae, const ExactInteger& b, int be) {
    if (a == 0 || b == 0)
        return a == b ? 0 : (a == 0 ? -1 : 1);
    const auto ab = static_cast<int>(boost::multiprecision::msb(a)) + 1;
    const auto bb = static_cast<int>(boost::multiprecision::msb(b)) + 1;
    const auto atop = static_cast<std::int64_t>(ab) + ae;
    const auto btop = static_cast<std::int64_t>(bb) + be;
    if (atop != btop)
        return atop < btop ? -1 : 1;
    if (ab < bb) {
        ExactInteger shifted = a;
        shifted <<= static_cast<unsigned>(bb - ab);
        return shifted == b ? 0 : (shifted < b ? -1 : 1);
    }
    ExactInteger shifted = b;
    shifted <<= static_cast<unsigned>(ab - bb);
    return a == shifted ? 0 : (a < shifted ? -1 : 1);
}
/// Decode a finite nonnegative word in units 2^-1074.
inline void local_decode(double word, ExactInteger& output) {
    const auto bits = std::bit_cast<std::uint64_t>(word);
    const auto exponent = (bits >> 52U) & 2047U;
    auto mantissa = bits & UINT64_C(0xfffffffffffff);
    if (exponent != 0)
        mantissa |= UINT64_C(0x10000000000000);
    output = mantissa;
    if (exponent != 0)
        output <<= static_cast<unsigned>(exponent - 1);
    if ((bits >> 63U) != 0)
        output = -output;
}
class LocalBounds {
    static core::Failure unresolved() noexcept {
        return {core::FailureCode::PrecisionUnresolved, core::Operation::EvaluateAero};
    }
    // floor of a positive dyadic value, and its exactly consumed integer units.
    static double floor_word(const ExactInteger& value, int exponent, ExactInteger& consumed) {
        consumed = 0;
        if (value == 0)
            return 0;
        const int bits = static_cast<int>(boost::multiprecision::msb(value)) + 1;
        const int shift = std::max({0, bits - 53, -1074 - exponent});
        if (shift >= bits)
            return 0;
        consumed = value >> static_cast<unsigned>(shift);
        const double word = std::scalbn(consumed.convert_to<double>(), exponent + shift);
        consumed <<= static_cast<unsigned>(shift);
        return word;
    }

public:
    /// Preserve dyadic low words even when the two significands are far apart.
    [[nodiscard]] static core::Result<ScaledEnclosure>
    dyadic(const ExactInteger& signed_value, int exponent, core::Precision precision,
           LocalRootScratch& scratch, bool uncertain = false) {
        auto& value = scratch.integers[8];
        auto& consumed = scratch.integers[9];
        auto& remainder = scratch.integers[10];
        value = signed_value < 0 ? -signed_value : signed_value;
        const bool negative = signed_value < 0;
        int scale = 0;
        if (value != 0) {
            const int top = static_cast<int>(boost::multiprecision::msb(value)) + exponent;
            const int bottom = static_cast<int>(boost::multiprecision::lsb(value)) + exponent;
            if (top > 1023 || bottom < -1074)
                scale = top;
        }
        const double high = floor_word(value, exponent - scale, consumed);
        remainder = value - consumed;
        double low = 0;
        if (precision == core::Precision::TwoWord) {
            low = floor_word(remainder, exponent - scale, consumed);
            remainder -= consumed;
        }
        if (uncertain)
            ++remainder; // Strict upper root endpoint, in the original dyadic units.
        double radius = floor_word(remainder, exponent - scale, consumed);
        if (consumed < remainder)
            radius = std::nextafter(radius, std::numeric_limits<double>::infinity());
        if (!std::isfinite(high) || !std::isfinite(low) || !std::isfinite(radius))
            return std::unexpected(unresolved());
        const double h = negative ? -high : high;
        const double l = negative ? -low : low;
        if (precision == core::Precision::TwoWord)
            return ScaledEnclosure{DoubleWord{h, l}, radius, scale};
        return ScaledEnclosure{SingleWord{h}, radius, scale};
    }
    /// Certified sqrt(N/D)*2^exponent, signed afterwards; fixed 257 exact trials.
    [[nodiscard]] static core::Result<ScaledEnclosure>
    root(const ExactInteger& numerator, const ExactInteger& denominator, int exponent,
         bool negative, core::Precision precision, LocalRootScratch& scratch, LocalWork& work,
         core::StageTiming* timing = nullptr) {
        const LocalClock clock{timing};
        auto& n = scratch.integers[0];
        auto& d = scratch.integers[1];
        auto& a = scratch.integers[2];
        auto& b = scratch.integers[3];
        auto& root = scratch.integers[4];
        auto& candidate = scratch.integers[5];
        auto& square = scratch.integers[6];
        auto& product = scratch.integers[7];
        n = numerator;
        d = denominator;
        if (n < 0 || d <= 0)
            return std::unexpected(unresolved());
        if (n == 0)
            return dyadic(n, 0, precision, scratch);
        int k = static_cast<int>(boost::multiprecision::msb(n)) -
                static_cast<int>(boost::multiprecision::msb(d));
        if (local_compare(n, 0, d, k) < 0)
            --k;
        const int t = k >= 0 ? k / 2 : -((-k + 1) / 2);
        const int shift = 512 - 2 * t;
        a = n;
        b = d;
        if (shift >= 0)
            a <<= static_cast<unsigned>(shift);
        else
            b <<= static_cast<unsigned>(-shift);
        root = 0;
        for (int bit = 256; bit >= 0; --bit) {
            if (auto check = work.step(); !check)
                return std::unexpected(check.error());
            candidate = 1;
            candidate <<= static_cast<unsigned>(bit);
            candidate += root;
            square = candidate * candidate;
            product = square * b;
            if (product <= a)
                root = candidate;
        }
        square = root * root;
        product = square * b;
        const bool uncertain = product != a;
        if (negative)
            root = -root;
        return dyadic(root, exponent + t - 256, precision, scratch, uncertain);
    }
    /// Compare full width 2r (not radius) against the requested original scale.
    [[nodiscard]] static bool width(const ScaledEnclosure& value, double tau,
                                    const ExactInteger& longest_squared, bool area, bool normal,
                                    LocalRootScratch& scratch) {
        auto& r = scratch.integers[0];
        auto& t = scratch.integers[1];
        auto& left = scratch.integers[2];
        auto& right = scratch.integers[3];
        auto& scale = scratch.integers[4];
        local_decode(value.radius(), r);
        local_decode(tau, t);
        if (normal)
            return local_compare(r, value.exponent() - 1073, t, -1074) <= 0;
        left = r * r;
        right = t * t;
        if (area) {
            scale = longest_squared * longest_squared;
            right *= scale;
        } else
            right *= longest_squared;
        return local_compare(left, 2 * value.exponent() - 2146, right, area ? -10740 : -6444) <= 0;
    }
};
} // namespace uvlas::aero::detail
