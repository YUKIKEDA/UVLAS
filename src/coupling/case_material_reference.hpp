#pragma once

#include "case_material_chart.hpp"
#include "case_reference_graph.hpp"

#include <optional>

namespace uvlas::coupling::detail {

struct CaseMaterialReferenceCandidate {
    core::Index partition_record{};
    std::array<structure::ExactRational, 3> local;
};

/// The exact reference point and the two chained physical tangent vectors use
/// its single owned radicand. The original layer face remains explicit.
struct CaseMaterialReferenceSample {
    CaseMaterialChartPoint chart;
    structure::PgReferenceValue reference;
    std::optional<std::array<structure::PgReferenceVector, 2>> tangents; ///< Collocations only.
    core::Index triangle;
    std::array<structure::ExactRational, 3> triangle_coordinates; ///< u, v, layer depth.
    core::PlannedBuffer<CaseMaterialReferenceCandidate> alternate;
    std::size_t alternate_count{};
    structure::PgReferenceCallInfo call_info;
};

[[nodiscard]] core::Result<CaseMaterialReferenceSample>
evaluate_case_material_collocation_reference(const CaseMaterialChart&, core::Index chord_cell,
                                             core::Index span_cell,
                                             const structure::PgReferenceLimits&,
                                             core::ResourcePlan&) noexcept;

[[nodiscard]] core::Result<CaseMaterialReferenceSample>
evaluate_case_material_vertex_reference(const CaseMaterialChart&, core::Index chord_column,
                                        core::Index span_knot, core::Index owning_span_cell,
                                        const structure::PgReferenceLimits&,
                                        core::ResourcePlan&) noexcept;

/// Compare with the independent reference from the same owner and use the Case's
/// own physical tolerance. The squared-distance and normal-sign comparisons are exact.
[[nodiscard]] core::Result<CaseMaterialReferenceSample> qualify_case_material_collocation_reference(
    const CaseMaterialChart&, const CaseReferenceSurface&, core::Index chord_cell,
    core::Index span_cell, const structure::PgReferenceLimits&, core::ResourcePlan&) noexcept;
[[nodiscard]] core::Result<CaseMaterialReferenceSample> qualify_case_material_vertex_reference(
    const CaseMaterialChart&, const CaseReferenceSurface&, core::Index chord_column,
    core::Index span_knot, core::Index owning_span_cell, const structure::PgReferenceLimits&,
    core::ResourcePlan&) noexcept;

} // namespace uvlas::coupling::detail
