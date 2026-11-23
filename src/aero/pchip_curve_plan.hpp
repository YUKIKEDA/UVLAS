#pragma once

#include "pchip_preparation_input.hpp"

namespace uvlas::aero::detail::pchip {

/// Conservative numerical preparation plan for one already-owned input.
/// Covers cumulative coefficient and boundary payload and the work of both forecasts
/// followed by both existing calls. Shared owner and allocator metadata are separate.
class CurveNumericsPlan {
public:
    [[nodiscard]] std::uint64_t payload_bytes() const noexcept {
        return bytes_;
    }
    [[nodiscard]] std::uint64_t reserved_steps() const noexcept {
        return work_;
    }
    [[nodiscard]] const CoefficientPlan& coefficients() const& noexcept {
        return coefficients_;
    }
    [[nodiscard]] const CoefficientPlan& coefficients() const&& = delete;
    [[nodiscard]] const BoundaryPlan& boundary() const& noexcept {
        return boundary_;
    }
    [[nodiscard]] const BoundaryPlan& boundary() const&& = delete;
    /// Derived coefficient inputs to the boundary may be wider than the source inputs.
    [[nodiscard]] RationalLimits boundary_limits() const noexcept {
        return boundary_limits_;
    }

private:
    friend struct boundary_detail::Builder;
    CoefficientPlan coefficients_;
    BoundaryPlan boundary_;
    RationalLimits boundary_limits_;
    std::uint64_t bytes_ = 0, work_ = 0;
};

/// Does not allocate or create coefficients. The plan is for this immutable input only.
[[nodiscard]] core::Result<CurveNumericsPlan>
forecast_curve_numerics(const PchipPreparationInput&, RationalLimits,
                        const std::stop_token& stop = {},
                        core::FailureDetails details = {}) noexcept;

} // namespace uvlas::aero::detail::pchip
