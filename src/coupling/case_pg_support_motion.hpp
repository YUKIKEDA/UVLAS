#pragma once

#include "case_flow_clock.hpp"

#include <uvlas/structure/pg_chart_candidate.hpp>
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

struct CaseStructuralSupportLimits {
    CasePgSupportLimits evaluation;
    structure::PgRigidSupportLimits exact;
    structure::PgTraceLimits region_trace;
    std::uint64_t maximum_group_supports{}, maximum_group_work{};
};

namespace structural_support_detail {
struct Builder;
}

class CaseStructuralSupport {
public:
    CaseStructuralSupport(CaseStructuralSupport&&) noexcept = default;
    CaseStructuralSupport& operator=(CaseStructuralSupport&&) noexcept = default;
    CaseStructuralSupport(const CaseStructuralSupport&) = delete;
    CaseStructuralSupport& operator=(const CaseStructuralSupport&) = delete;

    [[nodiscard]] io::SchemaIdentity source() const noexcept {
        return clock_.source;
    }
    [[nodiscard]] const structure::PgRigidSupportRows& support() const& noexcept {
        return support_;
    }
    [[nodiscard]] const structure::PgRigidSupportRows& support() const&& = delete;

private:
    friend struct structural_support_detail::Builder;
    CaseStructuralSupport(case_clock_detail::Grid&& clock, CaseMotionReady&& motion,
                          structure::PgRigidSupportRows&& support) noexcept
        : clock_(std::move(clock)), motion_(std::move(motion)), support_(std::move(support)) {}

    case_clock_detail::Grid clock_;
    CaseMotionReady motion_;
    structure::PgRigidSupportRows support_;
};

class CaseStructuralSupportSample {
public:
    CaseStructuralSupportSample(CaseStructuralSupportSample&&) noexcept = default;
    CaseStructuralSupportSample& operator=(CaseStructuralSupportSample&&) noexcept = default;
    CaseStructuralSupportSample(const CaseStructuralSupportSample&) = delete;
    CaseStructuralSupportSample& operator=(const CaseStructuralSupportSample&) = delete;

    [[nodiscard]] const CaseReferenceClockSample& node() const& noexcept {
        return node_;
    }
    [[nodiscard]] const CaseReferenceClockSample& node() const&& = delete;
    [[nodiscard]] const CaseMotionEvaluation& motion() const& noexcept {
        return motion_;
    }
    [[nodiscard]] const CaseMotionEvaluation& motion() const&& = delete;
    [[nodiscard]] structure::PgSparseChartIdentity chart() const noexcept {
        return chart_;
    }
    [[nodiscard]] std::span<const io::SiEnclosure> original_rhs() const& noexcept {
        return original_rhs_.values();
    }
    [[nodiscard]] std::span<const io::SiEnclosure> original_rhs() const&& = delete;
    [[nodiscard]] io::SchemaIdentity source() const noexcept {
        return source_;
    }
    [[nodiscard]] std::span<const bool> moving_rows() const& noexcept {
        return moving_.values();
    }
    [[nodiscard]] std::span<const bool> moving_rows() const&& = delete;
    [[nodiscard]] std::size_t rows() const noexcept {
        return rows_;
    }
    [[nodiscard]] std::size_t orders() const noexcept {
        return orders_;
    }

private:
    friend struct structural_support_detail::Builder;
    CaseStructuralSupportSample(CaseReferenceClockSample&& node, CaseMotionEvaluation&& motion,
                                structure::PgSparseChartIdentity chart, io::SchemaIdentity source,
                                core::PlannedBuffer<io::SiEnclosure>&& rhs,
                                core::PlannedBuffer<bool>&& moving, std::size_t rows,
                                std::size_t orders) noexcept
        : node_(std::move(node)),
          motion_(std::move(motion)),
          chart_(chart),
          source_(source),
          original_rhs_(std::move(rhs)),
          moving_(std::move(moving)),
          rows_(rows),
          orders_(orders) {}

    CaseReferenceClockSample node_;
    CaseMotionEvaluation motion_;
    structure::PgSparseChartIdentity chart_;
    io::SchemaIdentity source_;
    core::PlannedBuffer<io::SiEnclosure> original_rhs_; ///< Derivative, original row, axis.
    core::PlannedBuffer<bool> moving_;
    std::size_t rows_, orders_;
};

/// Resolve actual SupportTraceConstraint IDs in this immutable case. All selected
/// supports must use the same exact Motion recipe and complete named-edge trace.
[[nodiscard]] core::Result<CaseStructuralSupport> prepare_case_structural_pg_support(
    const io::SiCase&, structure::PgSparseChart&&,
    std::span<const structure::PgConstraintPatchBinding>,
    std::span<const structure::PgConstraintId>, aero::detail::pchip::PchipInputStamp,
    core::Workspace&, aero::detail::pchip::CurveWorkspace&, const CaseStructuralSupportLimits&,
    core::ResourcePlan&, const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

[[nodiscard]] core::Result<CaseStructuralSupportSample> evaluate_case_structural_pg_support(
    const CaseStructuralSupport&, const CaseStructuralSupport& current, CaseFlowClockQuery,
    aero::detail::pchip::JetOrder, const CasePgSupportProjection&,
    aero::detail::pchip::CurveWorkspace&, const CasePgSupportLimits&, core::ResourcePlan&,
    const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

namespace structural_support_detail {
struct Retention;
}

class CaseStructuralRetainedRhs {
public:
    CaseStructuralRetainedRhs(CaseStructuralRetainedRhs&&) noexcept = default;
    CaseStructuralRetainedRhs& operator=(CaseStructuralRetainedRhs&&) noexcept = default;
    CaseStructuralRetainedRhs(const CaseStructuralRetainedRhs&) = delete;
    CaseStructuralRetainedRhs& operator=(const CaseStructuralRetainedRhs&) = delete;

    [[nodiscard]] structure::PgSparseChartIdentity chart() const noexcept {
        return chart_;
    }
    [[nodiscard]] const numerics::PrecisionBuffer<double>& centre() const& noexcept {
        return centre_;
    }
    [[nodiscard]] const numerics::PrecisionBuffer<double>& centre() const&& = delete;
    [[nodiscard]] std::span<const double> radii() const& noexcept {
        return radii_.values();
    }
    [[nodiscard]] std::span<const double> radii() const&& = delete;
    [[nodiscard]] double maximum_radius() const noexcept {
        return maximum_radius_;
    }

private:
    friend struct structural_support_detail::Retention;
    CaseStructuralRetainedRhs(structure::PgSparseChartIdentity chart,
                              numerics::PrecisionBuffer<double>&& centre,
                              core::PlannedBuffer<double>&& radii, double maximum_radius) noexcept
        : chart_(chart),
          centre_(std::move(centre)),
          radii_(std::move(radii)),
          maximum_radius_(maximum_radius) {}

    structure::PgSparseChartIdentity chart_;
    numerics::PrecisionBuffer<double> centre_;
    core::PlannedBuffer<double> radii_;
    double maximum_radius_;
};

[[nodiscard]] core::Result<CaseStructuralRetainedRhs>
retain_case_structural_pg_rhs(const CaseStructuralSupportSample&, aero::detail::pchip::JetOrder,
                              core::NumericalGeneration, std::uint64_t maximum_work,
                              core::ResourcePlan&, const std::stop_token& = {},
                              core::FailureDetails = {}) noexcept;

/// Add the true RHS radius to the numerical chart's original-CSR certificate.
[[nodiscard]] core::Result<structure::PgChartCandidateLift> lift_case_structural_pg_support(
    structure::PgChartCandidate&, structure::PgChartCandidateInput independent,
    const CaseStructuralRetainedRhs&, const structure::PgChartCandidateSolve&, core::Workspace&,
    core::ResourcePlan&, core::FailureDetails = {}) noexcept;

/// Certify the affine work with uncertainty in the true RHS at this derivative.
[[nodiscard]] core::Result<double> work_case_structural_pg_support(
    const structure::PgChartCandidate&, structure::PgChartCandidateInput force,
    structure::PgChartCandidateInput independent, const CaseStructuralRetainedRhs&,
    const structure::PgChartCandidateLift&, const structure::PgChartCandidateDual&,
    double maximum_absolute_work_error, std::uint64_t maximum_work, const std::stop_token& = {},
    core::FailureDetails = {}) noexcept;

struct CaseStructuralExactLift {
    core::PlannedBuffer<structure::PgCartesianExactPosition> original_rhs, coordinates;
    structure::ExactRational maximum_rhs_error; ///< Difference from the physical RHS, in row units.
};

/// Build a compatible rational RHS using one retained shared motion centre.
[[nodiscard]] core::Result<CaseStructuralExactLift> lift_case_structural_pg_support_exact(
    const CaseStructuralSupport&, const CaseStructuralSupportSample&, aero::detail::pchip::JetOrder,
    std::span<const structure::PgCartesianExactPosition> independent,
    const structure::PgRigidSupportLimits&, core::Workspace&, core::ResourcePlan&) noexcept;

} // namespace uvlas::coupling::detail
