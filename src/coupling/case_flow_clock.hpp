#pragma once

#include "case_uniform_flow_projection.hpp"
#include "io/si_exact_words.hpp"

namespace uvlas::coupling::detail {

using CaseClockRestart = decltype(std::declval<io::schema::TimeInterval>().restart.selected());

struct CaseFlowClockLimits {
    CaseUniformFlowLimits fields; ///< The arithmetic and field caps apply to each evaluation role.
    std::uint64_t maximum_intervals = 0;
    std::uint64_t maximum_steps_per_interval = 0;
};

struct CaseFlowClockQuery {
    CaseFlowClockQuery(aero::detail::pchip::PchipInputStamp input, std::uint64_t segment,
                       std::uint64_t node, aero::detail::pchip::RequestedSide requested) noexcept
        : stamp(input), interval(segment), step(node), side(requested) {}

    aero::detail::pchip::PchipInputStamp stamp;
    std::uint64_t interval, step; ///< Zero-based interval and node. The terminal node is included.
    aero::detail::pchip::RequestedSide side;
};

namespace case_clock_detail {
struct Builder;
}

class CaseClockedFlows {
public:
    CaseClockedFlows(CaseClockedFlows&&) noexcept = default;
    CaseClockedFlows& operator=(CaseClockedFlows&&) noexcept = default;
    CaseClockedFlows(const CaseClockedFlows&) = delete;
    CaseClockedFlows& operator=(const CaseClockedFlows&) = delete;

    [[nodiscard]] std::size_t interval_count() const noexcept {
        return intervals_.size();
    }
    [[nodiscard]] core::Result<std::uint64_t> interval_steps(std::uint64_t) const noexcept;

private:
    friend struct case_clock_detail::Builder;
    struct Interval {
        structure::ExactRational begin, end, step;
        std::uint64_t steps = 0;
        CaseClockRestart restart = CaseClockRestart::ItemInitial;
    };
    CaseClockedFlows(io::SchemaIdentity source, aero::detail::pchip::PchipInputStamp stamp,
                     core::PlannedBuffer<Interval>&& intervals, CaseUniformFlowReady&& ambient,
                     CasePrescribedWakeReady&& wake) noexcept;

    io::SchemaIdentity source_;
    aero::detail::pchip::PchipInputStamp stamp_;
    core::PlannedBuffer<Interval> intervals_;
    CaseUniformFlowReady ambient_;
    CasePrescribedWakeReady wake_;
};

/// Own the exact clock coordinate and both projected fields for the same query.
class CaseClockedFlowSample {
public:
    CaseClockedFlowSample(CaseClockedFlowSample&&) noexcept = default;
    CaseClockedFlowSample& operator=(CaseClockedFlowSample&&) noexcept = default;
    CaseClockedFlowSample(const CaseClockedFlowSample&) = delete;
    CaseClockedFlowSample& operator=(const CaseClockedFlowSample&) = delete;

    [[nodiscard]] io::SchemaIdentity source() const noexcept {
        return source_;
    }
    [[nodiscard]] CaseFlowClockQuery address() const noexcept {
        return address_;
    }
    [[nodiscard]] CaseClockRestart restart() const noexcept {
        return restart_;
    }
    [[nodiscard]] core::Result<aero::detail::pchip::RationalView> time() const& noexcept;
    [[nodiscard]] core::Result<aero::detail::pchip::RationalView> time() const&& = delete;
    [[nodiscard]] core::Result<aero::detail::pchip::RationalView> step_seconds() const& noexcept;
    [[nodiscard]] core::Result<aero::detail::pchip::RationalView> step_seconds() const&& = delete;
    [[nodiscard]] const ProjectedCaseUniformFlow& ambient() const& noexcept {
        return ambient_;
    }
    [[nodiscard]] const ProjectedCaseUniformFlow& ambient() const&& = delete;
    [[nodiscard]] const ProjectedCasePrescribedWake& prescribed() const& noexcept {
        return wake_;
    }
    [[nodiscard]] const ProjectedCasePrescribedWake& prescribed() const&& = delete;

private:
    friend struct case_clock_detail::Builder;
    CaseClockedFlowSample(io::SchemaIdentity source, CaseFlowClockQuery address,
                          CaseClockRestart restart, io::detail::SiRationalWords&& instant,
                          io::detail::SiRationalWords&& interval_step,
                          ProjectedCaseUniformFlow&& ambient,
                          ProjectedCasePrescribedWake&& wake) noexcept;

    io::SchemaIdentity source_;
    CaseFlowClockQuery address_;
    CaseClockRestart restart_;
    io::detail::SiRationalWords instant_, interval_step_;
    ProjectedCaseUniformFlow ambient_;
    ProjectedCasePrescribedWake wake_;
};

/// Read the checked transient grid and its no-gust ambient and prescribed wake fields.
[[nodiscard]] core::Result<CaseClockedFlows>
prepare_case_clocked_flows(const io::SiCase&, aero::detail::pchip::PchipInputStamp,
                           aero::detail::UniformFlowWorkspace&, const CaseFlowClockLimits&,
                           core::ResourcePlan&, const std::stop_token& = {},
                           core::FailureDetails = {}) noexcept;

/// Validate the independent current grid and both fields before publishing a pair.
/// The addressed interval owns its right-sided initial and left-sided terminal node.
[[nodiscard]] core::Result<CaseClockedFlowSample> evaluate_case_clocked_flows(
    const CaseClockedFlows&, const CaseClockedFlows& current, CaseFlowClockQuery,
    const CaseFlowProjectionRequest& ambient, const CaseFlowProjectionRequest& wake,
    aero::detail::UniformFlowWorkspace&, const CaseFlowClockLimits&, core::ResourcePlan&,
    const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

} // namespace uvlas::coupling::detail
