#pragma once

#include "external_planar_pg.hpp"

namespace uvlas::io::detail {

struct ExternalPgThicknessClaim {
    ExternalPgThicknessClaim(std::uint64_t tag, const structure::ExactRational& metres,
                             std::optional<core::Index> declaration = {},
                             std::optional<NastranBulkPosition> card = {}) noexcept
        : element(tag), value_metres(metres), sidecar_node(declaration), field(card) {}

    std::uint64_t element;
    const structure::ExactRational& value_metres;
    std::optional<core::Index> sidecar_node;
    std::optional<NastranBulkPosition> field;
};

enum class ExternalPgThicknessFailure : std::uint8_t { Source, Selection, Difference };

struct ExternalPgThicknessDiagnostic {
    ExternalPgThicknessFailure reason = ExternalPgThicknessFailure::Source;
    std::optional<std::uint64_t> element{};
    std::optional<core::Index> sidecar_node{};
    std::optional<NastranBulkPosition> field{};
    std::optional<core::Index> case_triangle{};
    std::optional<std::size_t> overlap_vertex{};
};

struct ExternalPgThicknessLimits {
    ExternalPgThicknessLimits(const structure::PgReferenceLimits& arithmetic_limits,
                              std::uint64_t claims) noexcept
        : arithmetic(arithmetic_limits), maximum_claims(claims) {}

    structure::PgReferenceLimits arithmetic;
    std::uint64_t maximum_claims;
};

struct ExternalPgThicknessMatch {
    std::uint64_t element{};
    ExternalPgWinding winding = ExternalPgWinding::Forward;
    structure::ExactRational value_metres;
    structure::ExactRational twice_chart_area;
    std::uint64_t pieces{};
    std::uint64_t vertices{};
    std::optional<core::Index> sidecar_node{};
    std::optional<NastranBulkPosition> field{};
};

struct ExternalPgThicknessReport {
    RootSnapshotLocation source;
    SourceDigest digest;
    structure::PgReferenceIdentity pg;
    core::PlannedBuffer<ExternalPgThicknessMatch> matches;
    structure::PgReferenceCallInfo work;
};

// Compare one external constant per declared original face on the exact common
// subdivision between that face's certified chart polygon and the prepared PG.
// Each declaration applies only on its own face. Zero-area contacts do not compare
// values from the adjacent material region. The caller's owners remain read-only.
[[nodiscard]] core::Result<ExternalPgThicknessReport> compare_external_pg_regional_thickness(
    const ExternalPlanarPg&, const structure::PgReferenceGeometry&,
    std::span<const ExternalPgThicknessClaim>, const ExternalPgThicknessLimits&,
    core::ResourcePlan&, std::optional<ExternalPgThicknessDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::io::detail
