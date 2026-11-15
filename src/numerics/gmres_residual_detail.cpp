#include "gmres_residual_detail.hpp"

#include "compensated.hpp"
#include "linear_solve_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>

namespace uvlas::numerics::linear::detail {
namespace {
[[nodiscard]] auto fail(core::FailureCode code) noexcept {
    return std::unexpected(core::Failure{code, core::Operation::SolveLinear});
}
[[nodiscard]] auto arithmetic_failure(core::Failure failure) noexcept {
    // Every operand was already finite. Overflow of this evaluation does not
    // establish a nonfinite physical residual, which may still cancel exactly.
    if (failure.code == core::FailureCode::Nonfinite) {
        failure.code = core::FailureCode::PrecisionUnresolved;
    }
    return std::unexpected(failure);
}
} // namespace

core::Result<void>
combine_gmres_residual(core::Index rows, std::span<const double> b_hi, std::span<const double> b_lo,
                       std::span<const double> action_hi, std::span<const double> action_lo,
                       std::span<const double> action_hi_tail,
                       std::span<const double> action_lo_tail, std::span<const double> error_hi,
                       std::span<const double> error_lo, std::span<double> residual_out,
                       std::span<double> error_out) noexcept {
    if (!std::in_range<std::size_t>(rows)) {
        return fail(core::FailureCode::InvalidInput);
    }
    const auto size = static_cast<std::size_t>(rows);
    const std::array extents{b_hi.size(),      b_lo.size(),           action_hi.size(),
                             action_lo.size(), action_hi_tail.size(), action_lo_tail.size(),
                             error_hi.size(),  error_lo.size(),       residual_out.size(),
                             error_out.size()};
    if (std::ranges::any_of(extents, [size](std::size_t extent) { return extent != size; })) {
        return fail(core::FailureCode::InvalidInput);
    }
    // Validate the complete borrowed input before writing any scratch row.
    for (std::size_t i = 0; i < size; ++i) {
        if (!std::isfinite(b_hi[i]) || !std::isfinite(b_lo[i]) || !std::isfinite(action_hi[i]) ||
            !std::isfinite(action_lo[i]) || !std::isfinite(action_hi_tail[i]) ||
            !std::isfinite(action_lo_tail[i]) || !std::isfinite(error_hi[i]) ||
            !std::isfinite(error_lo[i])) {
            return fail(core::FailureCode::Nonfinite);
        }
        if (error_hi[i] < 0 || error_lo[i] < 0) {
            return fail(core::FailureCode::InvalidInput);
        }
    }
    for (std::size_t i = 0; i < size; ++i) {
        uvlas::numerics::detail::CompensatedAccumulator sum;
        for (const double term : std::array{b_hi[i], b_lo[i], -action_hi[i], -action_hi_tail[i],
                                            -action_lo[i], -action_lo_tail[i]}) {
            const auto added = sum.add(term);
            if (!added) {
                return arithmetic_failure(added.error());
            }
        }
        const auto value = sum.result();
        if (!value) {
            return arithmetic_failure(value.error());
        }
        if (!std::isfinite(value->error)) {
            return fail(core::FailureCode::PrecisionUnresolved);
        }
        const auto high_bound =
            add_magnitudes({value->error, value->error}, {error_hi[i], error_hi[i]});
        const auto total = add_magnitudes(high_bound, {error_lo[i], error_lo[i]});
        if (!std::isfinite(total.upper)) {
            return fail(core::FailureCode::PrecisionUnresolved);
        }
        residual_out[i] = value->value;
        error_out[i] = total.upper;
    }
    return {};
}
} // namespace uvlas::numerics::linear::detail
