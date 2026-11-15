#pragma once
#include <span>
#include <uvlas/numerics/linear_solve.hpp>

namespace uvlas::numerics::linear::detail {
struct MagnitudeBounds {
    double lower;
    double upper;
};
// Enclose abs(high+low), not abs(high)+abs(low). Nonfinite words return Nonfinite.
// A finite exact magnitude outside binary64 can have an infinite upper bound.
template <StorageScalar Scalar>
[[nodiscard]] core::Result<MagnitudeBounds> actual_magnitude_bounds(Scalar high,
                                                                    Scalar low) noexcept;
// All spans describe one packed logical column, with rows elements. rhs_lo alone
// may be empty (a binary64 RHS). Scales are positive finite; the same scales are
// broadcast across RHS columns by the caller. Error bounds and denominator lower
// bounds must be nonnegative. Physical words and scales reject nonfinite inputs.
// A valid but unmet criterion is returned as a diagnostic metric, not as success
// qualification. An unrepresentable conservative error upper is PrecisionUnresolved.
template <StorageScalar Scalar>
[[nodiscard]] core::Result<ColumnMetrics>
evaluate_column_metrics(core::Index rows, std::span<const Scalar> residual_estimate,
                        std::span<const double> absolute_error_upper,
                        std::span<const double> backward_denominator_lower,
                        std::span<const Scalar> rhs_hi, std::span<const Scalar> rhs_lo,
                        std::span<const double> physical_residual_scales,
                        const Accuracy& accuracy) noexcept;
// One common acceptance criterion: every column needs residual norm <=1 and its
// independently computed backward-error upper <= the requested backward budget.
[[nodiscard]] bool meets_accuracy(ColumnMetrics metrics, const Accuracy& accuracy) noexcept;
} // namespace uvlas::numerics::linear::detail
