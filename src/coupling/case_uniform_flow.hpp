#pragma once

#include "aero/uniform_flow.hpp"

#include <uvlas/io/si_case.hpp>

namespace uvlas::coupling::detail {

struct CaseUniformFlowLimits {
    structure::ExactLimits conversion; ///< Existing SI materialization limits for each recipe.
    aero::detail::pchip::CurveLimits
        flow; ///< Aggregate exact uniform evaluation/preparation limits.
    std::uint64_t maximum_exact_numbers = 0; ///< Total SI recipes imported in one preparation.
};

/// Required C1 window in exact SI seconds, supplied by the consumer for non-transient PCHIP.
struct CaseUniformTimeWindow {
    aero::detail::pchip::RationalView begin, end;
};

namespace case_flow_detail {
struct Builder;
}

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

} // namespace uvlas::coupling::detail
