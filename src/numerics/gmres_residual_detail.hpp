#pragma once

#include <span>
#include <uvlas/core/checked_dimensions.hpp>
#include <uvlas/core/failure.hpp>

namespace uvlas::numerics::linear::detail {
// Enclose b_hi+b_lo-(Ax_hi+its signed tail)-(Ax_lo+its signed tail), including both errors.
// All spans have rows entries. Output spans must not overlap any input or each other.
// Output is scratch: a failure may leave earlier rows written, never an owning candidate.
// Malformed sizes/negative errors are InvalidInput; nonfinite inputs are Nonfinite.
// Unrepresentable finite arithmetic/enclosures are PrecisionUnresolved.
[[nodiscard]] core::Result<void>
combine_gmres_residual(core::Index rows, std::span<const double> b_hi, std::span<const double> b_lo,
                       std::span<const double> action_hi, std::span<const double> action_lo,
                       std::span<const double> action_hi_tail,
                       std::span<const double> action_lo_tail, std::span<const double> error_hi,
                       std::span<const double> error_lo, std::span<double> residual_out,
                       std::span<double> error_out) noexcept;
} // namespace uvlas::numerics::linear::detail
