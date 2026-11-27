#pragma once
#include <cmath>
#include <limits>
#include <uvlas/numerics/precision.hpp>

namespace uvlas::numerics::detail {
using core::Result;

struct PerOperationEft {
    [[nodiscard]] Result<TwoWordValue> sum(double a, double b) const noexcept {
        return two_sum(a, b);
    }
    [[nodiscard]] Result<TwoWordValue> product(double a, double b) const noexcept {
        return two_product(a, b);
    }
};

/// A checked environment for a synchronous, callback-free computation on
/// the thread that obtained it. Validate again after changing the environment.
class CheckedEft {
public:
    CheckedEft(CheckedEft&&) noexcept = default;
    CheckedEft(const CheckedEft&) = delete;
    CheckedEft& operator=(CheckedEft&&) = delete;
    CheckedEft& operator=(const CheckedEft&) = delete;
    [[nodiscard]] static Result<CheckedEft> prepare() noexcept;
    [[nodiscard]] Result<TwoWordValue> sum(double a, double b) const noexcept;
    [[nodiscard]] Result<TwoWordValue> product(double a, double b) const noexcept;

private:
    CheckedEft() = default;
};

// Shared fixed-order compensated accumulation with an absolute enclosure of every
// discarded correction. No allocation or untracked precision conversion occurs here.
inline double outward_add(double a, double b) noexcept {
    if (a == 0) {
        return b;
    }
    if (b == 0) {
        return a;
    }
    return std::nextafter(a + b, std::numeric_limits<double>::infinity());
}
struct SumEvaluation {
    double value;
    double error;
};
class CompensatedAccumulator {
public:
    template <class Eft = PerOperationEft>
    Result<void> add(double value, const Eft& eft = {}) noexcept {
        const auto leading = eft.sum(sum_, value);
        if (!leading) {
            return std::unexpected(leading.error());
        }
        const auto trailing = eft.sum(correction_, leading->lo);
        if (!trailing) {
            return std::unexpected(trailing.error());
        }
        // Enclose rounding lost by compensation itself, even when the estimate becomes zero.
        error_ = outward_add(error_, std::abs(trailing->lo));
        sum_ = leading->hi;
        correction_ = trailing->hi;
        return {};
    }
    template <class Eft = PerOperationEft>
    Result<void> product(double a, double b, double sign, const Eft& eft = {}) noexcept {
        const auto value = eft.product(a, b);
        if (!value) {
            return std::unexpected(value.error());
        }
        const auto leading = add(sign * value->hi, eft);
        if (!leading) {
            return leading;
        }
        return add(sign * value->lo, eft);
    }
    template <class Eft = PerOperationEft>
    [[nodiscard]] Result<SumEvaluation> result(const Eft& eft = {}) const noexcept {
        const auto value = eft.sum(sum_, correction_);
        if (!value) {
            return std::unexpected(value.error());
        }
        return SumEvaluation{value->hi, outward_add(error_, std::abs(value->lo))};
    }

private:
    double sum_ = 0, correction_ = 0, error_ = 0;
};

inline Result<TwoWordValue> correct_two_word(double hi, double lo, double delta) noexcept {
    const auto first = two_sum(hi, delta);
    if (!first)
        return std::unexpected(first.error());
    const auto tail = two_sum(lo, first->lo);
    if (!tail)
        return std::unexpected(tail.error());
    const auto lead = two_sum(first->hi, tail->hi);
    if (!lead)
        return std::unexpected(lead.error());
    return two_sum(lead->hi, lead->lo + tail->lo);
}
template <StorageScalar Scalar>
Result<void> apply_two_word_correction(Scalar& hi, Scalar& lo, Scalar delta) noexcept {
    if constexpr (std::same_as<Scalar, double>) {
        auto value = correct_two_word(hi, lo, delta);
        if (!value)
            return std::unexpected(value.error());
        hi = value->hi;
        lo = value->lo;
    } else {
        auto r = correct_two_word(hi.real(), lo.real(), delta.real());
        if (!r)
            return std::unexpected(r.error());
        auto i = correct_two_word(hi.imag(), lo.imag(), delta.imag());
        if (!i)
            return std::unexpected(i.error());
        hi = {r->hi, i->hi};
        lo = {r->lo, i->lo};
    }
    return {};
}

} // namespace uvlas::numerics::detail
