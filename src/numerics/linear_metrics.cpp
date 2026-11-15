#include "linear_solve_detail.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>

namespace uvlas::numerics::linear::detail {
namespace {
using core::FailureCode;
using core::Result;
constexpr double infinity = std::numeric_limits<double>::infinity();
auto fail(FailureCode code) noexcept {
    return std::unexpected(core::Failure{code, core::Operation::SolveLinear});
}
bool finite(double value) noexcept {
    return std::isfinite(value);
}
bool finite(std::complex<double> value) noexcept {
    return finite(value.real()) && finite(value.imag());
}
double adjacent(double value, bool up) noexcept {
    return std::nextafter(value, up ? infinity : 0.0);
}
// Nonnegative directed arithmetic. Exact zero and exact representable operations
// are retained, while overflow/underflow remain honest one-sided enclosures.
double add_bound(double a, double b, bool up) noexcept {
    if (a == 0) {
        return b;
    }
    if (b == 0) {
        return a;
    }
    if (!finite(a) || !finite(b)) {
        return up ? infinity : std::numeric_limits<double>::max();
    }
    const auto value = two_sum(a, b);
    if (!value) {
        return up ? infinity : std::numeric_limits<double>::max();
    }
    if ((up && value->lo > 0) || (!up && value->lo < 0)) {
        return adjacent(value->hi, up);
    }
    return value->hi;
}
double multiply_bound(double a, double b, bool up) noexcept {
    if (a == 0 || b == 0) {
        return 0;
    }
    const double rounded = a * b;
    if (!finite(rounded)) {
        return up ? infinity : std::numeric_limits<double>::max();
    }
    const auto value = two_product(a, b);
    if (!value) {
        // An exact product below the subnormal grid is not an EFT, but its
        // correctly-rounded binary64 result still has adjacent outward bounds.
        return adjacent(rounded, up);
    }
    if ((up && value->lo > 0) || (!up && value->lo < 0)) {
        return adjacent(rounded, up);
    }
    return rounded;
}
double divide_bound(double a, double b, bool up) noexcept {
    if (a == 0) {
        return 0;
    }
    if (b == 0) {
        return up ? infinity : std::numeric_limits<double>::max();
    }
    const double rounded = a / b;
    if (!finite(rounded)) {
        return up ? infinity : std::numeric_limits<double>::max();
    }
    if (rounded == 0) {
        return up ? std::numeric_limits<double>::denorm_min() : 0;
    }
    const auto product = two_product(rounded, b);
    if (!product) {
        return adjacent(rounded, up);
    }
    const bool same = std::bit_cast<std::uint64_t>(product->hi) == std::bit_cast<std::uint64_t>(a);
    // Determine which side of the exact quotient the rounded result occupies.
    if (up) {
        if (product->hi > a || (same && product->lo >= 0)) {
            return rounded;
        }
    } else {
        if (product->hi < a || (same && product->lo <= 0)) {
            return rounded;
        }
    }
    return adjacent(rounded, up);
}
double sqrt_bound(double value, bool up) noexcept {
    if (value == 0) {
        return 0;
    }
    const double root = std::sqrt(value);
    const double difference = std::fma(root, root, -value);
    if ((up && difference < 0) || (!up && difference > 0)) {
        return adjacent(root, up);
    }
    return root;
}
double hypot_bound(double a, double b, bool up) noexcept {
    const double scale = std::max(a, b), other = std::min(a, b);
    if (scale == 0) {
        return 0;
    }
    if (!finite(scale)) {
        return up ? infinity : std::numeric_limits<double>::max();
    }
    const double ratio = divide_bound(other, scale, up);
    const double square = multiply_bound(ratio, ratio, up);
    const double radicand = add_bound(1, square, up);
    // The normalized radicand is near [1,2], so sqrt/FMA certification avoids
    // both extreme physical squares and assumptions about std::hypot rounding.
    return multiply_bound(scale, sqrt_bound(radicand, up), up);
}
Result<MagnitudeBounds> real_magnitude(double high, double low) noexcept {
    if (!finite(high) || !finite(low)) {
        return fail(FailureCode::Nonfinite);
    }
    const auto sum = two_sum(high, low);
    if (!sum) {
        if (sum.error().code != FailureCode::Nonfinite) {
            return std::unexpected(sum.error());
        }
        return MagnitudeBounds{std::numeric_limits<double>::max(), infinity};
    }
    const double leading = std::abs(sum->hi);
    const double remainder = std::signbit(sum->hi) ? -sum->lo : sum->lo;
    return MagnitudeBounds{remainder < 0 ? adjacent(leading, false) : leading,
                           remainder > 0 ? adjacent(leading, true) : leading};
}
} // namespace

MagnitudeBounds add_magnitudes(MagnitudeBounds a, MagnitudeBounds b) noexcept {
    return {add_bound(a.lower, b.lower, false), add_bound(a.upper, b.upper, true)};
}
MagnitudeBounds multiply_magnitudes(MagnitudeBounds a, MagnitudeBounds b) noexcept {
    return {multiply_bound(a.lower, b.lower, false), multiply_bound(a.upper, b.upper, true)};
}

template <StorageScalar Scalar>
Result<MagnitudeBounds> actual_magnitude_bounds(Scalar high, Scalar low) noexcept {
    if constexpr (std::same_as<Scalar, double>) {
        return real_magnitude(high, low);
    } else {
        const auto real = real_magnitude(high.real(), low.real());
        if (!real) {
            return std::unexpected(real.error());
        }
        const auto imag = real_magnitude(high.imag(), low.imag());
        if (!imag) {
            return std::unexpected(imag.error());
        }
        return MagnitudeBounds{hypot_bound(real->lower, imag->lower, false),
                               hypot_bound(real->upper, imag->upper, true)};
    }
}

template <StorageScalar Scalar>
Result<ColumnMetrics>
evaluate_column_metrics(core::Index rows, std::span<const Scalar> residual_estimate,
                        std::span<const double> absolute_error_upper,
                        std::span<const double> backward_denominator_lower,
                        std::span<const Scalar> rhs_hi, std::span<const Scalar> rhs_lo,
                        std::span<const double> physical_residual_scales,
                        const Accuracy& accuracy) noexcept {
    const auto request = check_accuracy_inputs(accuracy);
    if (!request) {
        return std::unexpected(request.error());
    }
    if (rows < 0 || std::cmp_not_equal(residual_estimate.size(), rows) ||
        std::cmp_not_equal(absolute_error_upper.size(), rows) ||
        std::cmp_not_equal(backward_denominator_lower.size(), rows) ||
        std::cmp_not_equal(rhs_hi.size(), rows) ||
        (!rhs_lo.empty() && std::cmp_not_equal(rhs_lo.size(), rows)) ||
        std::cmp_not_equal(physical_residual_scales.size(), rows)) {
        return fail(FailureCode::InvalidInput);
    }
    ColumnMetrics result;
    for (std::size_t i = 0; i < residual_estimate.size(); ++i) {
        if (!finite(residual_estimate[i]) || !finite(rhs_hi[i]) ||
            (!rhs_lo.empty() && !finite(rhs_lo[i])) || !finite(physical_residual_scales[i]) ||
            !finite(backward_denominator_lower[i])) {
            return fail(FailureCode::Nonfinite);
        }
        if (!finite(absolute_error_upper[i])) {
            return fail(FailureCode::PrecisionUnresolved);
        }
        if (absolute_error_upper[i] < 0 || backward_denominator_lower[i] < 0 ||
            physical_residual_scales[i] <= 0) {
            return fail(FailureCode::InvalidInput);
        }
        const auto residual = actual_magnitude_bounds(residual_estimate[i], Scalar{});
        if (!residual) {
            return std::unexpected(residual.error());
        }
        const auto rhs = actual_magnitude_bounds(rhs_hi[i], rhs_lo.empty() ? Scalar{} : rhs_lo[i]);
        if (!rhs) {
            return std::unexpected(rhs.error());
        }
        const double numerator = add_bound(residual->upper, absolute_error_upper[i], true);
        // Lower allowances are mandatory: an upward-rounded tau could turn
        // an actually unmet criterion into a false certificate.
        const double allowance = add_bound(
            multiply_bound(accuracy.absolute_tolerance, physical_residual_scales[i], false),
            multiply_bound(accuracy.relative_tolerance, rhs->lower, false), false);
        const double normalized = divide_bound(numerator, allowance, true);
        const double backward = divide_bound(numerator, backward_denominator_lower[i], true);
        if (!finite(normalized) || !finite(backward)) {
            return fail(FailureCode::PrecisionUnresolved);
        }
        result.normalized_original_residual_upper =
            accuracy.norm == Norm::Infinity
                ? std::max(result.normalized_original_residual_upper, normalized)
                : hypot_bound(result.normalized_original_residual_upper, normalized, true);
        result.componentwise_backward_error_upper =
            std::max(result.componentwise_backward_error_upper, backward);
    }
    if (!finite(result.normalized_original_residual_upper)) {
        return fail(FailureCode::PrecisionUnresolved);
    }
    return result;
}
bool meets_accuracy(ColumnMetrics metrics, const Accuracy& accuracy) noexcept {
    return finite(metrics.normalized_original_residual_upper) &&
           finite(metrics.componentwise_backward_error_upper) &&
           metrics.normalized_original_residual_upper >= 0 &&
           metrics.componentwise_backward_error_upper >= 0 &&
           finite(accuracy.componentwise_backward_tolerance) &&
           accuracy.componentwise_backward_tolerance >= 0 &&
           metrics.normalized_original_residual_upper <= 1 &&
           metrics.componentwise_backward_error_upper <= accuracy.componentwise_backward_tolerance;
}
template Result<MagnitudeBounds> actual_magnitude_bounds(double, double) noexcept;
template Result<MagnitudeBounds> actual_magnitude_bounds(std::complex<double>,
                                                         std::complex<double>) noexcept;
template Result<ColumnMetrics>
evaluate_column_metrics(core::Index, std::span<const double>, std::span<const double>,
                        std::span<const double>, std::span<const double>, std::span<const double>,
                        std::span<const double>, const Accuracy&) noexcept;
template Result<ColumnMetrics>
evaluate_column_metrics(core::Index, std::span<const std::complex<double>>, std::span<const double>,
                        std::span<const double>, std::span<const std::complex<double>>,
                        std::span<const std::complex<double>>, std::span<const double>,
                        const Accuracy&) noexcept;
} // namespace uvlas::numerics::linear::detail
