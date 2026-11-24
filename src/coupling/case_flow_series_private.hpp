#pragma once

#include "case_uniform_flow.hpp"

namespace uvlas::coupling::detail::case_flow_detail {

struct KnownSeries {
    CaseUniformFlowReady ambient;
    std::optional<CaseUniformFlowReady> amplitude;
};

[[nodiscard]] core::Result<KnownSeries>
prepare_known_series(const io::SiCase&, aero::detail::pchip::PchipInputStamp,
                     std::optional<CaseUniformTimeWindow>, aero::detail::UniformFlowWorkspace&,
                     const CaseUniformFlowLimits&, core::ResourcePlan&, const std::stop_token&,
                     core::FailureDetails) noexcept;

} // namespace uvlas::coupling::detail::case_flow_detail
