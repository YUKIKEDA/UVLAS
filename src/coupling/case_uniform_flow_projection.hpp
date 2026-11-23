#pragma once

#include "case_uniform_flow.hpp"

#include <uvlas/core/error_certificate.hpp>

namespace uvlas::coupling::detail {

/// The caller supplies the exact physical scale in each component's own unit.
struct CaseFlowComponentBudget {
    aero::detail::pchip::RationalView scale, absolute, relative;
};

struct CaseFlowProjectionLimits {
    std::uint64_t input_integer_bits = 0;
    std::uint64_t steps = 0;
};

struct CaseFlowProjectionRequest {
    CaseFlowProjectionRequest(core::Precision selected, std::array<CaseFlowComponentBudget, 3> u,
                              std::array<CaseFlowComponentBudget, 3> dt,
                              CaseFlowProjectionLimits bound) noexcept
        : representation(selected), velocity(u), time_partial(dt), limits(bound) {}

    core::Precision representation;
    std::array<CaseFlowComponentBudget, 3> velocity;     ///< XYZ, scales in m/s.
    std::array<CaseFlowComponentBudget, 3> time_partial; ///< XYZ, scales in m/s².
    CaseFlowProjectionLimits limits;
};

/// An owned finite center and outward radius for each jet component, with its exact parent.
class ProjectedCaseUniformFlow {
public:
    ProjectedCaseUniformFlow(ProjectedCaseUniformFlow&&) noexcept = default;
    ProjectedCaseUniformFlow& operator=(ProjectedCaseUniformFlow&&) noexcept = default;
    ProjectedCaseUniformFlow(const ProjectedCaseUniformFlow&) = delete;
    ProjectedCaseUniformFlow& operator=(const ProjectedCaseUniformFlow&) = delete;
    [[nodiscard]] const CaseUniformFlowEvaluation& source() const& noexcept;
    [[nodiscard]] const CaseUniformFlowEvaluation& source() const&& = delete;
    [[nodiscard]] core::Result<io::SiEnclosure> velocity(std::size_t component) const noexcept;
    [[nodiscard]] core::Result<io::SiEnclosure> time_partial(std::size_t component) const noexcept;
    [[nodiscard]] core::Result<io::SiEnclosure>
    spatial_partial(std::size_t component, std::size_t coordinate) const noexcept;
    [[nodiscard]] core::Precision representation() const noexcept;
    [[nodiscard]] std::uint64_t used_steps() const noexcept;

private:
    friend core::Result<ProjectedCaseUniformFlow>
    project_case_uniform_flow(CaseUniformFlowEvaluation&&, const CaseFlowProjectionRequest&,
                              const std::stop_token&, core::FailureDetails) noexcept;
    ProjectedCaseUniformFlow(CaseUniformFlowEvaluation&&, const std::array<io::SiEnclosure, 6>&,
                             core::Precision, std::uint64_t) noexcept;
    CaseUniformFlowEvaluation source_;
    std::array<io::SiEnclosure, 6> values_;
    core::Precision representation_;
    std::uint64_t work_;
};

/// Qualify the actual rounded words against E <= a*s + r*abs(exact value).
/// The exact parent is consumed only on success. The spatial partial is exactly zero.
[[nodiscard]] core::Result<ProjectedCaseUniformFlow>
project_case_uniform_flow(CaseUniformFlowEvaluation&&, const CaseFlowProjectionRequest&,
                          const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

} // namespace uvlas::coupling::detail
