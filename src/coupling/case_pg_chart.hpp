#pragma once

#include "case_motion_frame.hpp"
#include "case_pg_reference.hpp"
#include "case_pg_support_motion.hpp"

#include <uvlas/structure/pg_prepared_patch.hpp>
#include <uvlas/structure/pg_sparse_chart.hpp>

namespace uvlas::coupling::detail {

/// The finite points and weights are supplied for each named shell in Case order.
struct CasePgRuleBinding {
    std::uint64_t shell;
    std::span<const structure::PgDiscreteVolumePoint> points;
};

struct CasePgChartLimits {
    CasePgReferenceLimits reference;
    structure::PgRankPreparationLimits rank;
    structure::PgDiscreteGramGrid gram_grid;
    structure::PgDiscreteGramLimits gram;
    structure::PgConditioningGeneration factor_grid;
    structure::PgConditioningLimits factor;
    structure::PgReferenceLimits patch;
    CaseMotionFrameLimits support_frame;
    structure::PgTraceLimits trace;
    structure::PgGlobalConstraintLimits constraints;
    structure::PgSparseChartLimits chart;
    std::uint64_t maximum_input_work;
};

namespace case_pg_chart_detail {
struct Builder;
}

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

    [[nodiscard]] const structure::PgSparseChart& chart() const& noexcept;
    [[nodiscard]] const structure::PgSparseChart& chart() const&& = delete;
    [[nodiscard]] const CaseStructuralSupport* motion() const& noexcept;
    [[nodiscard]] const CaseStructuralSupport* motion() const&& = delete;

private:
    friend struct case_pg_chart_detail::Builder;
    CasePgChart(core::PlannedBuffer<std::optional<structure::PgPreparedPatch>>&& patches,
                structure::PgSparseChart&& chart) noexcept
        : patches_(std::move(patches)), chart_(std::move(chart)) {}

    core::PlannedBuffer<std::optional<structure::PgPreparedPatch>> patches_;
    // The fixed owner moves into the motion result only after successful qualification.
    structure::PgSparseChart chart_;
    std::optional<CaseStructuralSupport> motion_;
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
