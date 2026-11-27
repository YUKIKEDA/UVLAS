#pragma once

#include "case_motion_frame.hpp"
#include "case_pg_reference.hpp"
#include "case_pg_support_motion.hpp"

#include <uvlas/structure/pg_prepared_patch.hpp>
#include <uvlas/structure/pg_sparse_chart.hpp>
#include <variant>

namespace uvlas::coupling::detail {

/// One caller-supplied finite point list or positive product rule per named Shell in Case order.
struct CasePgRuleBinding {
    std::uint64_t shell;
    std::variant<std::span<const structure::PgDiscreteVolumePoint>, structure::PgVolumeProductRule>
        quadrature;
};

struct CasePgChartLimits {
    CasePgReferenceLimits reference;
    structure::PgRankPreparationLimits rank;
    structure::PgDiscreteGramGrid gram_grid{};
    structure::PgDiscreteGramLimits gram;
    structure::PgConditioningGeneration factor_grid;
    structure::PgConditioningLimits factor;
    structure::PgReferenceLimits patch;
    CaseMotionFrameLimits support_frame;
    structure::PgTraceLimits trace;
    structure::PgGlobalConstraintLimits constraints;
    structure::PgSparseChartLimits chart;
    std::uint64_t maximum_input_work{};
};

namespace case_pg_chart_detail {
struct Builder;
}

/// Synchronous borrow of a conditioned Shell from its owning Case preparation.
class CaseStaticPgPointPatch {
public:
    [[nodiscard]] numerics::PatchId shell() const noexcept {
        return shell_;
    }
    [[nodiscard]] io::SchemaIdentity source() const noexcept {
        return source_;
    }
    [[nodiscard]] const structure::PgPreparedPatch& patch() const noexcept {
        return patch_;
    }

private:
    friend class CaseComponentPgChart;
    friend struct case_pg_chart_detail::Builder;
    CaseStaticPgPointPatch(io::SchemaIdentity source, numerics::PatchId shell,
                           const structure::PgPreparedPatch& patch) noexcept
        : source_(source), shell_(shell), patch_(patch) {}
    io::SchemaIdentity source_;
    numerics::PatchId shell_;
    const structure::PgPreparedPatch& patch_;
};

class CaseComponentPgChart {
public:
    CaseComponentPgChart(CaseComponentPgChart&&) noexcept = default;
    CaseComponentPgChart& operator=(CaseComponentPgChart&&) noexcept = default;

    [[nodiscard]] std::size_t patch_count() const noexcept {
        return patches_.size();
    }
    [[nodiscard]] core::Result<const structure::PgPreparedPatch*>
        patch(core::Index) const& noexcept;
    [[nodiscard]] core::Result<const structure::PgPreparedPatch*>
        patch(core::Index) const&& = delete;
    [[nodiscard]] core::Result<CaseStaticPgPointPatch> point_patch(core::Index) const& noexcept;
    [[nodiscard]] core::Result<CaseStaticPgPointPatch> point_patch(core::Index) const&& = delete;
    [[nodiscard]] const structure::PgComponentSparseChart& chart() const& noexcept {
        return chart_;
    }
    [[nodiscard]] const structure::PgComponentSparseChart& chart() const&& = delete;

private:
    friend struct case_pg_chart_detail::Builder;
    CaseComponentPgChart(io::SchemaIdentity source,
                         core::PlannedBuffer<std::optional<structure::PgPreparedPatch>>&& patches,
                         structure::PgComponentSparseChart&& chart) noexcept
        : source_(source), patches_(std::move(patches)), chart_(std::move(chart)) {}
    io::SchemaIdentity source_;
    core::PlannedBuffer<std::optional<structure::PgPreparedPatch>> patches_;
    structure::PgComponentSparseChart chart_;
};

/// Prepare a static Case with full-depth traces and selected world components
/// of centre-surface point supports in the same PG bases.
[[nodiscard]] core::Result<CaseComponentPgChart>
prepare_case_component_pg_chart(const io::SiCase&, std::span<const CasePgRuleBinding>,
                                const CasePgChartLimits&, core::Workspace&, core::ResourcePlan&,
                                std::stop_token = {}) noexcept;

/// Synchronous borrow from a CasePgChart. Keep its owner live and do not move it while borrowed.
class CasePgShellReference {
public:
    [[nodiscard]] const structure::PgPreparedPatch& patch() const noexcept {
        return patch_;
    }
    [[nodiscard]] const structure::PgSparseChart& chart() const noexcept {
        return chart_;
    }
    [[nodiscard]] structure::PgConstraintPatchLayout columns() const noexcept {
        return columns_;
    }
    [[nodiscard]] structure::PgSparseChartIdentity chart_identity() const noexcept {
        return chart_id_;
    }
    [[nodiscard]] std::span<const structure::PgOwnedDiscreteVolumePoint>
    fixed_rule() const noexcept {
        return rule_;
    }

private:
    friend class CasePgChart;
    CasePgShellReference(const structure::PgPreparedPatch& patch,
                         const structure::PgSparseChart& chart,
                         structure::PgConstraintPatchLayout columns,
                         structure::PgSparseChartIdentity chart_id,
                         std::span<const structure::PgOwnedDiscreteVolumePoint> rule) noexcept
        : patch_(patch), chart_(chart), columns_(columns), chart_id_(chart_id), rule_(rule) {}

    const structure::PgPreparedPatch& patch_;
    const structure::PgSparseChart& chart_;
    structure::PgConstraintPatchLayout columns_;
    structure::PgSparseChartIdentity chart_id_;
    std::span<const structure::PgOwnedDiscreteVolumePoint> rule_;
};

/// Owns the exact same references and column order used in the chart.
class CasePgChart {
public:
    CasePgChart(CasePgChart&&) noexcept = default;
    CasePgChart& operator=(CasePgChart&&) noexcept = default;
    CasePgChart(const CasePgChart&) = delete;
    CasePgChart& operator=(const CasePgChart&) = delete;

    [[nodiscard]] std::size_t patch_count() const noexcept {
        return patches_.values().size();
    }
    [[nodiscard]] core::Result<const structure::PgPreparedPatch*>
    patch(core::Index index) const& noexcept;
    [[nodiscard]] core::Result<const structure::PgPreparedPatch*>
        patch(core::Index) const&& = delete;

    /// Resolve one physical Shell ID and reconcile its actual basis with the original columns.
    [[nodiscard]] core::Result<CasePgShellReference>
    shell_reference(numerics::PatchId, std::uint64_t maximum_patch_visits,
                    std::stop_token = {}) const& noexcept;
    [[nodiscard]] core::Result<CasePgShellReference>
    shell_reference(numerics::PatchId, std::uint64_t, std::stop_token = {}) const&& = delete;

    [[nodiscard]] const structure::PgSparseChart& chart() const& noexcept;
    [[nodiscard]] const structure::PgSparseChart& chart() const&& = delete;
    [[nodiscard]] const CaseStructuralSupport* motion() const& noexcept;
    [[nodiscard]] const CaseStructuralSupport* motion() const&& = delete;

    /// Attach or replace numerical preparation on the same exact chart, preserving it on failure.
    /// Do not retain borrows of this chart or its numerical candidate across successful
    /// preparation.
    [[nodiscard]] core::Result<void>
    prepare_numerical(const numerics::linear::Scales&, const numerics::linear::Scales&,
                      const structure::PgChartCandidatePreparation&,
                      core::ResourcePlan&) & noexcept;
    [[nodiscard]] const structure::PgChartCandidate* numerical() const& noexcept;
    [[nodiscard]] const structure::PgChartCandidate* numerical() const&& = delete;

private:
    friend struct case_pg_chart_detail::Builder;
    CasePgChart(core::PlannedBuffer<std::optional<structure::PgPreparedPatch>>&& patches,
                structure::PgSparseChart&& chart) noexcept
        : patches_(std::move(patches)),
          chart_(std::in_place_type<structure::PgChartOwner>, std::move(chart)) {}
    CasePgChart(core::PlannedBuffer<std::optional<structure::PgPreparedPatch>>&& patches,
                CaseStructuralSupport&& motion) noexcept
        : patches_(std::move(patches)), chart_(std::move(motion)) {}

    core::PlannedBuffer<std::optional<structure::PgPreparedPatch>> patches_;
    std::variant<structure::PgChartOwner, CaseStructuralSupport> chart_;
};

/// Prepare every PG shell, fixed full-position edge support and compatible
/// full-depth affine edge adhesive from one immutable Case. The finite rule
/// qualifies discrete conditioning for those retained points and weights.
[[nodiscard]] core::Result<CasePgChart> prepare_case_pg_chart(const io::SiCase&,
                                                              std::span<const CasePgRuleBinding>,
                                                              const CasePgChartLimits&,
                                                              core::Workspace&, core::ResourcePlan&,
                                                              std::stop_token = {}) noexcept;

/// Prepare the same fixed chart and bind every actual moving support to the Case clock.
[[nodiscard]] core::Result<CasePgChart> prepare_case_pg_motion_chart(
    const io::SiCase&, std::span<const CasePgRuleBinding>, const CasePgChartLimits&,
    aero::detail::pchip::PchipInputStamp, const CaseStructuralSupportLimits&, core::Workspace&,
    aero::detail::pchip::CurveWorkspace&, core::ResourcePlan&, std::stop_token = {}) noexcept;

} // namespace uvlas::coupling::detail
