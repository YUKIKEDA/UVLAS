#pragma once

#include "aero/uniform_flow.hpp"

#include <utility>
#include <uvlas/io/si_case.hpp>

namespace uvlas::coupling::detail {

struct CaseUniformFlowLimits {
    structure::ExactLimits conversion; ///< Existing SI materialization limits for each recipe.
    aero::detail::pchip::CurveLimits
        flow; ///< Aggregate exact uniform evaluation/preparation limits.
    /// SI recipes imported for one time curve, including typed coverage endpoints.
    /// Provider padding is accounted for by the storage budget.
    std::uint64_t maximum_exact_numbers = 0;
};

/// Required C1 window in exact SI seconds, supplied by the consumer for non-transient PCHIP.
struct CaseUniformTimeWindow {
    aero::detail::pchip::RationalView begin, end;
};

namespace case_flow_detail {
struct Builder;
struct PrescribedProjection;
} // namespace case_flow_detail

/// Owns the exact velocity provider prepared from the case's typed no-gust flow.
class CaseUniformFlowReady {
public:
    CaseUniformFlowReady(CaseUniformFlowReady&&) noexcept = default;
    CaseUniformFlowReady& operator=(CaseUniformFlowReady&&) noexcept = default;
    CaseUniformFlowReady(const CaseUniformFlowReady&) = delete;
    CaseUniformFlowReady& operator=(const CaseUniformFlowReady&) = delete;
    [[nodiscard]] const aero::detail::UniformFlowReady& uniform() const& noexcept;
    [[nodiscard]] const aero::detail::UniformFlowReady& uniform() const&& = delete;

private:
    friend struct case_flow_detail::Builder;
    explicit CaseUniformFlowReady(aero::detail::UniformFlowReady) noexcept;
    aero::detail::UniformFlowReady flow_;
};

/// Explicit query time in SI seconds and requested side. Stamp is issued by the input host.
struct CaseUniformFlowQuery {
    CaseUniformFlowQuery(aero::detail::pchip::PchipInputStamp identity,
                         aero::detail::pchip::RationalView instant,
                         aero::detail::pchip::RequestedSide side_requested =
                             aero::detail::pchip::RequestedSide::Automatic) noexcept
        : stamp(identity), time(instant), side(side_requested) {}

    aero::detail::pchip::PchipInputStamp stamp;
    aero::detail::pchip::RationalView time;
    aero::detail::pchip::RequestedSide side = aero::detail::pchip::RequestedSide::Automatic;
};

struct CaseUniformFlowEvaluation {
    aero::detail::ExactUniformFlowJet jet;
    std::uint64_t used_steps = 0; ///< Includes current-input comparison.
    aero::detail::pchip::RequestedSide requested_side =
        aero::detail::pchip::RequestedSide::Automatic;
};

/// Transient PCHIP uses the checked case intervals. Other PCHIP callers supply a C1 window.
[[nodiscard]] core::Result<CaseUniformFlowReady>
prepare_case_uniform_flow(const io::SiCase&, aero::detail::pchip::PchipInputStamp,
                          std::optional<CaseUniformTimeWindow>, aero::detail::UniformFlowWorkspace&,
                          const CaseUniformFlowLimits&, core::ResourcePlan&,
                          const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

/// The host supplies its current, independently prepared case input for every evaluation.
[[nodiscard]] core::Result<aero::detail::pchip::CurveAdmission>
forecast_case_uniform_flow(const CaseUniformFlowReady&, const CaseUniformFlowReady& current,
                           CaseUniformFlowQuery, aero::detail::pchip::CurveLimits,
                           const std::stop_token& = {}, core::FailureDetails = {}) noexcept;
[[nodiscard]] core::Result<CaseUniformFlowEvaluation>
evaluate_case_uniform_flow(const CaseUniformFlowReady&, const CaseUniformFlowReady& current,
                           CaseUniformFlowQuery, aero::detail::UniformFlowWorkspace&,
                           aero::detail::pchip::CurveLimits, core::ResourcePlan&,
                           const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

[[nodiscard]] core::Result<aero::detail::UniformFlowSnapshot>
snapshot_case_uniform_flow(const CaseUniformFlowReady&, aero::detail::UniformFlowWorkspace&,
                           aero::detail::pchip::CurveLimits, core::ResourcePlan&,
                           const std::stop_token& = {}, core::FailureDetails = {}) noexcept;
[[nodiscard]] core::Result<CaseUniformFlowReady>
restore_case_uniform_flow(const aero::detail::UniformFlowSnapshotView&,
                          const CaseUniformFlowReady& current, aero::detail::UniformFlowWorkspace&,
                          aero::detail::pchip::CurveLimits, core::ResourcePlan&,
                          const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

/// Own the independently specified wake-transport velocity curve.
class CasePrescribedWakeReady {
public:
    CasePrescribedWakeReady(CasePrescribedWakeReady&&) noexcept = default;
    CasePrescribedWakeReady& operator=(CasePrescribedWakeReady&&) noexcept = default;
    CasePrescribedWakeReady(const CasePrescribedWakeReady&) = delete;
    CasePrescribedWakeReady& operator=(const CasePrescribedWakeReady&) = delete;

private:
    friend struct case_flow_detail::Builder;
    explicit CasePrescribedWakeReady(CaseUniformFlowReady&& provider) noexcept
        : provider_(std::move(provider)) {}
    CaseUniformFlowReady provider_;
};

class CasePrescribedWakeEvaluation {
public:
    CasePrescribedWakeEvaluation(CasePrescribedWakeEvaluation&&) noexcept = default;
    CasePrescribedWakeEvaluation& operator=(CasePrescribedWakeEvaluation&&) noexcept = default;
    CasePrescribedWakeEvaluation(const CasePrescribedWakeEvaluation&) = delete;
    CasePrescribedWakeEvaluation& operator=(const CasePrescribedWakeEvaluation&) = delete;

    [[nodiscard]] const aero::detail::ExactUniformFlowJet& jet() const& noexcept {
        return value_.jet;
    }
    [[nodiscard]] const aero::detail::ExactUniformFlowJet& jet() const&& = delete;
    [[nodiscard]] aero::detail::pchip::RequestedSide requested_side() const noexcept {
        return value_.requested_side;
    }
    [[nodiscard]] std::uint64_t used_steps() const noexcept {
        return value_.jet.active() ? value_.used_steps : 0;
    }

private:
    friend struct case_flow_detail::Builder;
    friend struct case_flow_detail::PrescribedProjection;
    explicit CasePrescribedWakeEvaluation(CaseUniformFlowEvaluation&& value) noexcept
        : value_(std::move(value)) {}
    CaseUniformFlowEvaluation value_;
};

/// Owns the curve input and prepared coefficients, tagged as prescribed convection.
class CasePrescribedWakeSnapshot {
public:
    CasePrescribedWakeSnapshot(CasePrescribedWakeSnapshot&&) noexcept = default;
    CasePrescribedWakeSnapshot& operator=(CasePrescribedWakeSnapshot&&) noexcept = default;
    CasePrescribedWakeSnapshot(const CasePrescribedWakeSnapshot&) = delete;
    CasePrescribedWakeSnapshot& operator=(const CasePrescribedWakeSnapshot&) = delete;

private:
    friend struct case_flow_detail::Builder;
    explicit CasePrescribedWakeSnapshot(aero::detail::UniformFlowSnapshot&& snapshot) noexcept
        : snapshot_(std::move(snapshot)) {}
    aero::detail::UniformFlowSnapshot snapshot_;
};

/// Read WakeActive.convection=prescribed from the typed SI case. The velocity
/// curve is independent of the ambient FlowKnown/gust and has zero spatial derivative.
[[nodiscard]] core::Result<CasePrescribedWakeReady> prepare_case_prescribed_wake(
    const io::SiCase&, aero::detail::pchip::PchipInputStamp, std::optional<CaseUniformTimeWindow>,
    aero::detail::UniformFlowWorkspace&, const CaseUniformFlowLimits&, core::ResourcePlan&,
    const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

[[nodiscard]] core::Result<aero::detail::pchip::CurveAdmission>
forecast_case_prescribed_wake(const CasePrescribedWakeReady&,
                              const CasePrescribedWakeReady& current, CaseUniformFlowQuery,
                              aero::detail::pchip::CurveLimits, const std::stop_token& = {},
                              core::FailureDetails = {}) noexcept;
[[nodiscard]] core::Result<CasePrescribedWakeEvaluation> evaluate_case_prescribed_wake(
    const CasePrescribedWakeReady&, const CasePrescribedWakeReady& current, CaseUniformFlowQuery,
    aero::detail::UniformFlowWorkspace&, aero::detail::pchip::CurveLimits, core::ResourcePlan&,
    const std::stop_token& = {}, core::FailureDetails = {}) noexcept;
[[nodiscard]] core::Result<CasePrescribedWakeSnapshot>
snapshot_case_prescribed_wake(const CasePrescribedWakeReady&, aero::detail::UniformFlowWorkspace&,
                              aero::detail::pchip::CurveLimits, core::ResourcePlan&,
                              const std::stop_token& = {}, core::FailureDetails = {}) noexcept;
[[nodiscard]] core::Result<CasePrescribedWakeReady> restore_case_prescribed_wake(
    const CasePrescribedWakeSnapshot&, const CasePrescribedWakeReady& current,
    aero::detail::UniformFlowWorkspace&, aero::detail::pchip::CurveLimits, core::ResourcePlan&,
    const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

} // namespace uvlas::coupling::detail
