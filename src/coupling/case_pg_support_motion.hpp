#pragma once

#include "case_flow_clock.hpp"

#include <uvlas/structure/pg_rigid_support_rows.hpp>

namespace uvlas::coupling::detail {

struct CasePgSupportProjection {
    core::Precision representation;
    std::span<const CaseFlowComponentBudget> motion;       ///< 12 components per order.
    std::span<const CaseFlowComponentBudget> probe;        ///< 3 components per order per probe.
    std::span<const CaseFlowComponentBudget> original_rhs; ///< 3 components per order per row.
};

struct CasePgSupportLimits {
    CaseClockedPrescribedLimits clock;
    std::uint64_t maximum_row_width_bits{}, maximum_row_work{}, maximum_original_rows{};
};

struct CasePgSupportSample {
    CaseClockedPrescribedSample inputs;
    structure::PgSparseChartIdentity chart;
    core::PlannedBuffer<io::SiEnclosure> original_rhs; ///< Derivative, original row, axis.
    std::size_t rows{}, orders{};
    std::uint64_t row_work{};
};

/// Evaluate the single shared Motion on the adopted clock, then all original support rows.
/// Keep actual enclosures and their radius with the clock/side and exact chart provenance.
[[nodiscard]] core::Result<CasePgSupportSample> evaluate_case_clocked_pg_support(
    const structure::PgRigidSupportRows&, const CaseClockedPrescribedInputs&,
    const CaseClockedPrescribedInputs& current, CaseFlowClockQuery,
    const std::array<aero::detail::pchip::RationalView, 3>& flow_position,
    const CaseKnownFlowProjection&, const std::optional<CaseFlowProjectionRequest>& wake,
    aero::detail::pchip::JetOrder, const CasePgSupportProjection&,
    aero::detail::UniformFlowWorkspace&, aero::detail::pchip::CurveWorkspace&,
    const CasePgSupportLimits&, core::ResourcePlan&, const std::stop_token& = {},
    core::FailureDetails = {}) noexcept;

} // namespace uvlas::coupling::detail
