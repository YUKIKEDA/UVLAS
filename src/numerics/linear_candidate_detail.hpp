#pragma once

#include "linear_solve_detail.hpp"

#include <cmath>
#include <complex>
#include <concepts>
#include <utility>
#include <variant>

namespace uvlas::numerics::linear::detail {
// A raw residual is usable for refinement, but its normalized conservative
// metric cannot be represented finitely. This state can never certify a result.
struct UnrepresentableMetrics {};
// Default construction is deliberately unqualified, including prepared scratch.
using CandidateMetrics = std::variant<UnrepresentableMetrics, ColumnMetrics>;
struct RefinementProgress {
    CandidateMetrics current;
    CandidateMetrics previous;
};

[[nodiscard]] inline bool candidate_meets_accuracy(const CandidateMetrics& candidate,
                                                   const Accuracy& accuracy) noexcept {
    const auto* metric = std::get_if<ColumnMetrics>(&candidate);
    return metric && meets_accuracy(*metric, accuracy);
}

// Call only for a binary64 candidate after at least one attempted refinement.
// Unrepresented ratios cannot demonstrate the adopted half-progress condition;
// conservatively promote without restarting the total refinement counter.
[[nodiscard]] inline bool needs_precision_promotion(const CandidateMetrics& current,
                                                    const CandidateMetrics& previous,
                                                    const Accuracy& accuracy) noexcept {
    if (candidate_meets_accuracy(current, accuracy))
        return false;
    const auto* now = std::get_if<ColumnMetrics>(&current);
    const auto* before = std::get_if<ColumnMetrics>(&previous);
    if (!now || !before)
        return true;
    return (now->normalized_original_residual_upper <= 1 ||
            now->normalized_original_residual_upper >
                0.5 * before->normalized_original_residual_upper) &&
           (now->componentwise_backward_error_upper <= accuracy.componentwise_backward_tolerance ||
            now->componentwise_backward_error_upper >
                0.5 * before->componentwise_backward_error_upper);
}

template <StorageScalar Scalar>
[[nodiscard]] bool candidate_scalar_finite(Scalar value) noexcept {
    if constexpr (std::same_as<Scalar, double>)
        return std::isfinite(value);
    else
        return std::isfinite(value.real()) && std::isfinite(value.imag());
}

// Only normalization/aggregation failure after a complete valid raw column is
// an unmet candidate. In particular, a nonfinite producer/accumulation error
// remains a genuine enclosure failure, even after an earlier row's ratio/0.
// No allocation, owning state or public certificate is created here.
template <StorageScalar Scalar>
[[nodiscard]] core::Result<CandidateMetrics>
assess_candidate_column_metrics(core::Index rows, std::span<const Scalar> residual_estimate,
                                std::span<const double> absolute_error_upper,
                                std::span<const double> backward_denominator_lower,
                                std::span<const Scalar> rhs_hi, std::span<const Scalar> rhs_lo,
                                std::span<const double> physical_residual_scales,
                                const Accuracy& accuracy) noexcept {
    const auto fail = [](core::FailureCode code) {
        return std::unexpected(core::Failure{code, core::Operation::SolveLinear});
    };
    const auto request = check_accuracy_inputs(accuracy);
    if (!request)
        return std::unexpected(request.error());
    if (rows < 0 || std::cmp_not_equal(residual_estimate.size(), rows) ||
        std::cmp_not_equal(absolute_error_upper.size(), rows) ||
        std::cmp_not_equal(backward_denominator_lower.size(), rows) ||
        std::cmp_not_equal(rhs_hi.size(), rows) ||
        (!rhs_lo.empty() && std::cmp_not_equal(rhs_lo.size(), rows)) ||
        std::cmp_not_equal(physical_residual_scales.size(), rows))
        return fail(core::FailureCode::InvalidInput);
    for (std::size_t i = 0; i < residual_estimate.size(); ++i) {
        if (!candidate_scalar_finite(residual_estimate[i]) || !candidate_scalar_finite(rhs_hi[i]) ||
            (!rhs_lo.empty() && !candidate_scalar_finite(rhs_lo[i])) ||
            !std::isfinite(physical_residual_scales[i]) ||
            !std::isfinite(backward_denominator_lower[i]))
            return fail(core::FailureCode::Nonfinite);
        if (!std::isfinite(absolute_error_upper[i]))
            return fail(core::FailureCode::PrecisionUnresolved);
        if (absolute_error_upper[i] < 0 || backward_denominator_lower[i] < 0 ||
            physical_residual_scales[i] <= 0)
            return fail(core::FailureCode::InvalidInput);
    }
    const auto metric = evaluate_column_metrics<Scalar>(
        rows, residual_estimate, absolute_error_upper, backward_denominator_lower, rhs_hi, rhs_lo,
        physical_residual_scales, accuracy);
    if (!metric) {
        if (metric.error().code == core::FailureCode::PrecisionUnresolved)
            return CandidateMetrics{UnrepresentableMetrics{}};
        return std::unexpected(metric.error());
    }
    return CandidateMetrics{*metric};
}
} // namespace uvlas::numerics::linear::detail
