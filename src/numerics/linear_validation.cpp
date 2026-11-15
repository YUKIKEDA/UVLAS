#include "linear_solve_detail.hpp"

#include <cmath>

namespace uvlas::numerics::linear::detail {
namespace {
auto fail(core::FailureCode code) noexcept {
    return std::unexpected(core::Failure{code, core::Operation::SolveLinear});
}
} // namespace
core::Result<void> check_context_inputs(const Context& context) noexcept {
    if (!std::isfinite(context.time)) {
        return fail(core::FailureCode::Nonfinite);
    }
    if (context.cache.representation != core::Precision::Binary64 &&
        context.cache.representation != core::Precision::TwoWord) {
        return fail(core::FailureCode::PrecisionUnsupported);
    }
    return {};
}
core::Result<void> check_scale_inputs(core::Index rows, const Scales& scales,
                                      const Context& expected) noexcept {
    const auto context = check_context_inputs(scales.context);
    if (!context) {
        return context;
    }
    if (rows < 0 || std::cmp_not_equal(scales.unknown.size(), rows) ||
        std::cmp_not_equal(scales.residual.size(), rows) ||
        std::cmp_not_equal(scales.unknown_units.size(), rows) ||
        std::cmp_not_equal(scales.residual_units.size(), rows)) {
        return fail(core::FailureCode::InvalidInput);
    }
    for (std::size_t i = 0; i < scales.unknown.size(); ++i) {
        if (!std::isfinite(scales.unknown[i]) || !std::isfinite(scales.residual[i])) {
            return fail(core::FailureCode::Nonfinite);
        }
        if (scales.unknown[i] <= 0 || scales.residual[i] <= 0) {
            return fail(core::FailureCode::InvalidInput);
        }
    }
    if (scales.context != expected) {
        return fail(core::FailureCode::StaleGeneration);
    }
    return {};
}
core::Result<void> check_capacity_inputs(core::Index rows, core::Index rhs_columns,
                                         core::Index update_rank, const Limits& limits) noexcept {
    if (rows < 0 || rhs_columns < 0 || update_rank < 0 || limits.maximum_order < 0 ||
        limits.maximum_rhs_columns < 0 || limits.maximum_update_rank < 0 ||
        limits.maximum_refinements > 5) {
        return fail(core::FailureCode::InvalidInput);
    }
    if (rows > limits.maximum_order || rhs_columns > limits.maximum_rhs_columns ||
        update_rank > limits.maximum_update_rank) {
        return fail(core::FailureCode::ResourceExceeded);
    }
    return {};
}
template <StorageScalar Scalar>
core::Result<void> check_solve_inputs(core::Index rows, core::Index update_rank,
                                      const Context& expected, const Rhs<Scalar>& rhs,
                                      const Scales& scales, const Precision<Scalar>& precision,
                                      const Accuracy& accuracy, const Limits& limits) noexcept {
    const auto context = check_context_inputs(expected);
    if (!context) {
        return context;
    }
    const auto rhs_context = check_context_inputs(rhs.context());
    if (!rhs_context) {
        return rhs_context;
    }
    const auto scale = check_scale_inputs(rows, scales, expected);
    if (!scale) {
        return scale;
    }
    const auto shape = rhs.values().shape();
    if (shape.rows != rows) {
        return fail(core::FailureCode::InvalidInput);
    }
    if (const auto* binary = std::get_if<Binary64<Scalar>>(&rhs.values().storage())) {
        const auto valid = binary->values.validate(core::Operation::SolveLinear);
        if (!valid) {
            return valid;
        }
    } else {
        // PrecisionBuffer's private, nonthrowing moves cannot make its variant valueless.
        const auto* pair = std::get_if<TwoWord<Scalar>>(&rhs.values().storage());
        const auto high = pair->hi.validate(core::Operation::SolveLinear);
        if (!high) {
            return high;
        }
        const auto low = pair->lo.validate(core::Operation::SolveLinear);
        if (!low) {
            return low;
        }
    }
    const auto request = check_accuracy_inputs(accuracy);
    if (!request) {
        return request;
    }
    const auto checks = check_consumer_inputs(precision.required_checks);
    if (!checks) {
        return checks;
    }
    if (!precision.required_checks.empty()) {
        const auto supported = precision.evaluator.supported_precision();
        const auto capability = require_precision(core::Precision::Binary64, supported);
        if (!capability) {
            return capability;
        }
        // Projection evaluation also reads the exact retained source, so a known
        // TwoWord RHS cannot bypass that capability requirement through its target tag.
        const auto initial = require_precision(rhs.values().precision(), supported);
        if (!initial) {
            return initial;
        }
    }
    if (precision.returned.valueless_by_exception()) {
        return fail(core::FailureCode::PrecisionUnsupported);
    }
    if (const auto* projection = std::get_if<Binary64Projection>(&precision.returned)) {
        if (projection->retained_generation == projection->target_generation) {
            return fail(core::FailureCode::InvalidInput);
        }
    }
    if (rhs.context() != expected || !core::can_reuse(expected.cache, rhs.cache_request())) {
        return fail(core::FailureCode::StaleGeneration);
    }
    return check_capacity_inputs(rows, shape.columns, update_rank, limits);
}
template core::Result<void> check_solve_inputs(core::Index, core::Index, const Context&,
                                               const Rhs<double>&, const Scales&,
                                               const Precision<double>&, const Accuracy&,
                                               const Limits&) noexcept;
template core::Result<void> check_solve_inputs(core::Index, core::Index, const Context&,
                                               const Rhs<std::complex<double>>&, const Scales&,
                                               const Precision<std::complex<double>>&,
                                               const Accuracy&, const Limits&) noexcept;
} // namespace uvlas::numerics::linear::detail
