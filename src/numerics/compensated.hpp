#pragma once
#include <cmath>
#include <limits>
#include <uvlas/numerics/precision.hpp>

namespace uvlas::numerics::detail {
using core::Result;
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
    Result<void> add(double value) noexcept {
        const auto leading = two_sum(sum_, value);
        if (!leading) {
            return std::unexpected(leading.error());
        }
        const auto trailing = two_sum(correction_, leading->lo);
        if (!trailing) {
            return std::unexpected(trailing.error());
        }
        // Enclose rounding lost by compensation itself, even when the estimate becomes zero.
        error_ = outward_add(error_, std::abs(trailing->lo));
        sum_ = leading->hi;
        correction_ = trailing->hi;
        return {};
    }
    Result<void> product(double a, double b, double sign) noexcept {
        const auto value = two_product(a, b);
        if (!value) {
            return std::unexpected(value.error());
        }
        const auto leading = add(sign * value->hi);
        if (!leading) {
            return leading;
        }
        return add(sign * value->lo);
    }
    [[nodiscard]] Result<SumEvaluation> result() const noexcept {
        const auto value = two_sum(sum_, correction_);
        if (!value) {
            return std::unexpected(value.error());
        }
        return SumEvaluation{value->hi, outward_add(error_, std::abs(value->lo))};
    }

private:
    double sum_ = 0, correction_ = 0, error_ = 0;
};

} // namespace uvlas::numerics::detail
