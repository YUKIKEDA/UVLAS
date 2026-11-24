#pragma once

#include "case_uniform_flow_projection.hpp"

namespace uvlas::coupling::detail {

struct CaseKnownFlowLimits {
    CaseUniformFlowLimits time_curves;    ///< Exact SI preparation and work caps per time curve.
    std::uint64_t maximum_input_bits = 0; ///< Canonical integer inputs to interval arithmetic.
    std::uint64_t interval_steps = 0; ///< Composition, phase reduction, projection and matching.
};

struct CaseKnownFlowQuery {
    CaseKnownFlowQuery(CaseUniformFlowQuery time_query,
                       const std::array<aero::detail::pchip::RationalView, 3>& x) noexcept
        : time(time_query), position(x) {}
    CaseUniformFlowQuery time;
    std::array<aero::detail::pchip::RationalView, 3> position; ///< Global Cartesian SI metres.
};

struct CaseKnownFlowProjection {
    CaseKnownFlowProjection(core::Precision precision,
                            const std::array<CaseFlowComponentBudget, 3>& speed,
                            const std::array<CaseFlowComponentBudget, 3>& acceleration,
                            const std::array<CaseFlowComponentBudget, 9>& gradient) noexcept
        : representation(precision),
          velocity(speed),
          time_partial(acceleration),
          spatial_partial(gradient) {}
    core::Precision representation;
    std::array<CaseFlowComponentBudget, 3> velocity, time_partial;
    std::array<CaseFlowComponentBudget, 9>
        spatial_partial; ///< Row-major velocity component / coordinate.
};

namespace known_flow_detail {
struct Prepared;
struct Saved;
struct Position;
struct Builder;
} // namespace known_flow_detail

class CaseKnownFlowReady {
public:
    CaseKnownFlowReady(CaseKnownFlowReady&&) noexcept;
    CaseKnownFlowReady& operator=(CaseKnownFlowReady&&) noexcept;
    CaseKnownFlowReady(const CaseKnownFlowReady&) = delete;
    CaseKnownFlowReady& operator=(const CaseKnownFlowReady&) = delete;
    ~CaseKnownFlowReady();

private:
    friend struct known_flow_detail::Builder;
    explicit CaseKnownFlowReady(core::PlannedBuffer<known_flow_detail::Prepared>&&) noexcept;
    core::PlannedBuffer<known_flow_detail::Prepared> data_;
};

class CaseKnownFlowSnapshot {
public:
    CaseKnownFlowSnapshot(CaseKnownFlowSnapshot&&) noexcept;
    CaseKnownFlowSnapshot& operator=(CaseKnownFlowSnapshot&&) noexcept;
    CaseKnownFlowSnapshot(const CaseKnownFlowSnapshot&) = delete;
    CaseKnownFlowSnapshot& operator=(const CaseKnownFlowSnapshot&) = delete;
    ~CaseKnownFlowSnapshot();

private:
    friend struct known_flow_detail::Builder;
    explicit CaseKnownFlowSnapshot(core::PlannedBuffer<known_flow_detail::Saved>&&) noexcept;
    core::PlannedBuffer<known_flow_detail::Saved> data_;
};

/// Projected value and physical partial derivatives at the same exact time and position.
class CaseKnownFlowJet {
public:
    CaseKnownFlowJet(CaseKnownFlowJet&&) noexcept;
    CaseKnownFlowJet& operator=(CaseKnownFlowJet&&) noexcept;
    CaseKnownFlowJet(const CaseKnownFlowJet&) = delete;
    CaseKnownFlowJet& operator=(const CaseKnownFlowJet&) = delete;
    ~CaseKnownFlowJet();
    [[nodiscard]] core::Result<aero::detail::pchip::PchipInputStamp> stamp() const noexcept;
    [[nodiscard]] core::Result<aero::detail::pchip::RationalView>
    position(std::size_t) const& noexcept;
    [[nodiscard]] core::Result<aero::detail::pchip::RationalView>
    position(std::size_t) const&& = delete;
    [[nodiscard]] core::Result<aero::detail::pchip::RationalView> time() const& noexcept;
    [[nodiscard]] core::Result<aero::detail::pchip::RationalView> time() const&& = delete;
    [[nodiscard]] aero::detail::pchip::RequestedSide requested_side() const noexcept;
    [[nodiscard]] core::Result<io::SiEnclosure> velocity(std::size_t) const noexcept;
    [[nodiscard]] core::Result<io::SiEnclosure> time_partial(std::size_t) const noexcept;
    [[nodiscard]] core::Result<io::SiEnclosure> spatial_partial(std::size_t,
                                                                std::size_t) const noexcept;
    [[nodiscard]] std::uint64_t interval_steps() const noexcept {
        return time_.jet.active() ? steps_ : 0;
    }

private:
    friend struct known_flow_detail::Builder;
    CaseKnownFlowJet(CaseUniformFlowEvaluation&&, const std::array<io::SiEnclosure, 15>&,
                     core::PlannedBuffer<known_flow_detail::Position>&&,
                     std::uint64_t steps) noexcept;
    CaseUniformFlowEvaluation time_;
    core::PlannedBuffer<known_flow_detail::Position> position_;
    std::array<io::SiEnclosure, 15> values_;
    std::uint64_t steps_;
};

/// Prepare the supported constant/polynomial/PCHIP ambient and optional analytic sine.
[[nodiscard]] core::Result<CaseKnownFlowReady>
prepare_case_known_flow(const io::SiCase&, aero::detail::pchip::PchipInputStamp,
                        std::optional<CaseUniformTimeWindow>, aero::detail::UniformFlowWorkspace&,
                        const CaseKnownFlowLimits&, core::ResourcePlan&,
                        const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

[[nodiscard]] core::Result<CaseKnownFlowJet> evaluate_case_known_flow(
    const CaseKnownFlowReady&, const CaseKnownFlowReady& current, const CaseKnownFlowQuery&,
    const CaseKnownFlowProjection&, aero::detail::UniformFlowWorkspace&, const CaseKnownFlowLimits&,
    core::ResourcePlan&, const std::stop_token& = {}, core::FailureDetails = {}) noexcept;
[[nodiscard]] core::Result<CaseKnownFlowSnapshot>
snapshot_case_known_flow(const CaseKnownFlowReady&, aero::detail::UniformFlowWorkspace&,
                         const CaseKnownFlowLimits&, core::ResourcePlan&,
                         const std::stop_token& = {}, core::FailureDetails = {}) noexcept;
[[nodiscard]] core::Result<CaseKnownFlowReady>
restore_case_known_flow(const CaseKnownFlowSnapshot&, const CaseKnownFlowReady& current,
                        aero::detail::UniformFlowWorkspace&, const CaseKnownFlowLimits&,
                        core::ResourcePlan&, const std::stop_token& = {},
                        core::FailureDetails = {}) noexcept;

} // namespace uvlas::coupling::detail
