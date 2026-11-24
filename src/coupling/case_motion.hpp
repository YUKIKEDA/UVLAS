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
    std::uint64_t maximum_interval_input_bits{};
    std::uint64_t maximum_interval_work{};
};

struct CaseMotionReferencePoint {
    std::array<aero::detail::pchip::RationalView, 3> world_metres;
    /// Three world components per requested derivative order.
    std::span<const CaseFlowComponentBudget> components;
};

struct CaseMotionProjection {
    CaseMotionProjection(core::Precision format, std::span<const CaseFlowComponentBudget> budgets,
                         std::span<const CaseMotionReferencePoint> points = {}) noexcept
        : representation(format), components(budgets), reference_points(points) {}
    core::Precision representation;
    /// Per derivative order: three world centre components, then nine row-major
    /// rotation components. Each scale uses its component's native units.
    std::span<const CaseFlowComponentBudget> components;
    std::span<const CaseMotionReferencePoint> reference_points{};
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

/// Stored exact inputs for this motion recipe. The version identifies the rigid
/// composition and formula-series rules; each PCHIP archive has its own version.
struct CaseMotionSnapshotView {
    std::uint64_t formula_version{};
    CaseMotionFrameSnapshotView frame;
    CaseMotionTimeSnapshotView time;
};

class CaseMotionSnapshot {
public:
    CaseMotionSnapshot(CaseMotionSnapshot&&) noexcept = default;
    CaseMotionSnapshot& operator=(CaseMotionSnapshot&&) noexcept = default;
    CaseMotionSnapshot(const CaseMotionSnapshot&) = delete;
    CaseMotionSnapshot& operator=(const CaseMotionSnapshot&) = delete;

    [[nodiscard]] core::Result<CaseMotionSnapshotView> view() const& noexcept;
    [[nodiscard]] core::Result<CaseMotionSnapshotView> view() const&& = delete;

private:
    friend struct case_motion_detail::Builder;
    CaseMotionSnapshot(CaseMotionFrameSnapshot&&, CaseMotionTimeSnapshot&&) noexcept;
    CaseMotionFrameSnapshot frame_;
    CaseMotionTimeSnapshot time_;
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
    [[nodiscard]] std::size_t point_count() const noexcept;
    [[nodiscard]] core::Result<io::SiEnclosure> point(std::size_t point, std::size_t axis,
                                                      aero::detail::pchip::JetOrder) const noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const CaseMotionTimeJet>>
    time() const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const CaseMotionTimeJet>>
    time() const&& = delete;
    [[nodiscard]] std::uint64_t interval_work() const noexcept;

private:
    friend struct case_motion_detail::Builder;
    static constexpr std::size_t components_per_order = 12;
    using Values = std::array<io::SiEnclosure, 3 * components_per_order>;
    CaseMotionEvaluation(CaseMotionTimeJet, const Values&, core::PlannedBuffer<io::SiEnclosure>,
                         std::size_t orders, std::uint64_t) noexcept;

    CaseMotionTimeJet time_;
    Values values_;
    core::PlannedBuffer<io::SiEnclosure> points_;
    std::size_t orders_;
    std::uint64_t interval_work_;
};

[[nodiscard]] core::Result<CaseMotionReady>
prepare_case_motion(const io::SiCase&, CaseMotionField, CaseMotionTimeWindow,
                    aero::detail::pchip::PchipInputStamp, aero::detail::pchip::CurveWorkspace&,
                    const CaseMotionLimits&, core::ResourcePlan&, const std::stop_token& = {},
                    core::FailureDetails = {}) noexcept;

[[nodiscard]] core::Result<CaseMotionSnapshot>
snapshot_case_motion(const CaseMotionReady&, aero::detail::pchip::CurveWorkspace&,
                     const CaseMotionLimits&, core::ResourcePlan&, const std::stop_token& = {},
                     core::FailureDetails = {}) noexcept;

/// Requalify saved inputs against a separately prepared current case. Saved
/// curves never become active until their actual coefficients have been checked.
[[nodiscard]] core::Result<CaseMotionReady>
restore_case_motion(const CaseMotionSnapshotView&, const CaseMotionReady& current,
                    aero::detail::pchip::CurveWorkspace&, const CaseMotionLimits&,
                    core::ResourcePlan&, const std::stop_token& = {},
                    core::FailureDetails = {}) noexcept;

[[nodiscard]] core::Result<CaseMotionEvaluation>
evaluate_case_motion(const CaseMotionReady&, const CaseMotionReady& current, CaseMotionTimeQuery,
                     const CaseMotionProjection&, aero::detail::pchip::CurveWorkspace&,
                     const CaseMotionLimits&, core::ResourcePlan&, const std::stop_token& = {},
                     core::FailureDetails = {}) noexcept;

} // namespace uvlas::coupling::detail
