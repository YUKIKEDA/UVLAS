#pragma once

#include "case_motion_frame.hpp"
#include "case_motion_time_series.hpp"
#include "case_uniform_flow_projection.hpp"

#include <array>
#include <functional>
#include <span>

namespace uvlas::coupling::detail {
namespace case_motion_detail {
struct Builder;
}

struct CaseMotionLimits {
    CaseMotionFrameLimits frame;
    CaseMotionTimeLimits time_series;
    std::uint64_t maximum_interval_input_bits;
    std::uint64_t maximum_interval_work;
};

struct CaseMotionProjection {
    core::Precision representation;
    /// Per derivative order: three world centre components, then nine row-major
    /// rotation components. Each scale uses its component's native units.
    std::span<const CaseFlowComponentBudget> components;
};

class CaseMotionReady {
public:
    CaseMotionReady(CaseMotionReady&&) noexcept = default;
    CaseMotionReady& operator=(CaseMotionReady&&) noexcept = default;
    CaseMotionReady(const CaseMotionReady&) = delete;
    CaseMotionReady& operator=(const CaseMotionReady&) = delete;

    [[nodiscard]] std::span<const aero::detail::pchip::RationalView> events() const& noexcept;
    [[nodiscard]] std::span<const aero::detail::pchip::RationalView> events() const&& = delete;

private:
    friend struct case_motion_detail::Builder;
    CaseMotionReady(CaseMotionFrame, CaseMotionTimeSeries, motion::m::Interval) noexcept;

    CaseMotionFrame frame_;
    CaseMotionTimeSeries curves_;
    motion::m::Interval pi_;
};

/// Owns the qualified time and the projected world position/rotation jet.
class CaseMotionEvaluation {
public:
    CaseMotionEvaluation(CaseMotionEvaluation&&) noexcept;
    CaseMotionEvaluation& operator=(CaseMotionEvaluation&&) noexcept;
    CaseMotionEvaluation(const CaseMotionEvaluation&) = delete;
    CaseMotionEvaluation& operator=(const CaseMotionEvaluation&) = delete;

    [[nodiscard]] core::Result<io::SiEnclosure>
    centre(std::size_t axis, aero::detail::pchip::JetOrder) const noexcept;
    [[nodiscard]] core::Result<io::SiEnclosure>
    rotation(std::size_t row, std::size_t column, aero::detail::pchip::JetOrder) const noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const CaseMotionTimeJet>>
    time() const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const CaseMotionTimeJet>>
    time() const&& = delete;
    [[nodiscard]] std::uint64_t interval_work() const noexcept;

private:
    friend struct case_motion_detail::Builder;
    static constexpr std::size_t components_per_order = 12;
    using Values = std::array<io::SiEnclosure, 3 * components_per_order>;
    CaseMotionEvaluation(CaseMotionTimeJet, Values, std::size_t orders, std::uint64_t) noexcept;

    CaseMotionTimeJet time_;
    Values values_;
    std::size_t orders_;
    std::uint64_t interval_work_;
};

[[nodiscard]] core::Result<CaseMotionReady>
prepare_case_motion(const io::SiCase&, CaseMotionField, CaseMotionTimeWindow,
                    aero::detail::pchip::PchipInputStamp, aero::detail::pchip::CurveWorkspace&,
                    const CaseMotionLimits&, core::ResourcePlan&, const std::stop_token& = {},
                    core::FailureDetails = {}) noexcept;

[[nodiscard]] core::Result<CaseMotionEvaluation>
evaluate_case_motion(const CaseMotionReady&, const CaseMotionReady& current, CaseMotionTimeQuery,
                     const CaseMotionProjection&, aero::detail::pchip::CurveWorkspace&,
                     const CaseMotionLimits&, core::ResourcePlan&, const std::stop_token& = {},
                     core::FailureDetails = {}) noexcept;

} // namespace uvlas::coupling::detail
