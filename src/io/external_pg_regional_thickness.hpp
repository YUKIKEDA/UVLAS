#pragma once

#include "external_planar_pg.hpp"

namespace uvlas::io::detail {

struct ExternalPgCornerThicknessClaim {
    std::uint64_t node{};
    const structure::ExactRational* value_metres{};
    std::optional<NastranBulkPosition> field{};
    std::optional<NastranBulkPosition> scale_field{};
};

struct ExternalPgThicknessClaim {
    ExternalPgThicknessClaim() noexcept = default;
    ExternalPgThicknessClaim(std::uint64_t tag, const structure::ExactRational& metres,
                             std::optional<core::Index> declaration = {},
                             std::optional<NastranBulkPosition> card = {}) noexcept
        : element(tag), value_metres(&metres), sidecar_node(declaration), field(card) {}

    ExternalPgThicknessClaim(std::uint64_t tag,
                             std::span<const ExternalPgCornerThicknessClaim> samples,
                             std::optional<core::Index> declaration = {}) noexcept
        : element(tag), sidecar_node(declaration), corners(samples) {}

    std::uint64_t element{};
    const structure::ExactRational* value_metres{};
    std::optional<core::Index> sidecar_node{};
    std::optional<NastranBulkPosition> field{};
    std::span<const ExternalPgCornerThicknessClaim> corners{};
    std::optional<NastranPropertyId> property{};
};

enum class ExternalPgThicknessFailure : std::uint8_t { Source, Selection, Difference };

struct ExternalPgThicknessDiagnostic {
    ExternalPgThicknessFailure reason = ExternalPgThicknessFailure::Source;
    std::optional<std::uint64_t> element{};
    std::optional<core::Index> sidecar_node{};
    std::optional<NastranBulkPosition> field{};
    std::optional<core::Index> case_triangle{};
    std::optional<std::size_t> overlap_vertex{};
    std::optional<std::uint64_t> node{};
    std::optional<NastranBulkPosition> scale_field{};
    std::optional<NastranPropertyId> property{};
};

struct ExternalPgThicknessLimits {
    ExternalPgThicknessLimits(const structure::PgReferenceLimits& arithmetic_limits,
                              std::uint64_t claims, std::uint64_t vertices) noexcept
        : arithmetic(arithmetic_limits), maximum_claims(claims), case_vertices(vertices) {}

    structure::PgReferenceLimits arithmetic;
    std::uint64_t maximum_claims;
    std::uint64_t case_vertices;
};

struct ExternalPgCornerThickness {
    std::uint64_t node{};
    structure::ExactRational value_metres;
    std::optional<NastranBulkPosition> field{};
    std::optional<NastranBulkPosition> scale_field{};
};

struct ExternalPgNodalThickness {
    std::array<ExternalPgCornerThickness, 4> values;
    std::size_t count{};
};

struct ExternalPgThicknessMatch {
    std::uint64_t element{};
    ExternalPgWinding winding = ExternalPgWinding::Forward;
    std::variant<structure::ExactRational, ExternalPgNodalThickness> value_metres;
    structure::ExactRational twice_chart_area;
    std::uint64_t pieces{};
    std::uint64_t vertices{}; // Vertices of the clipped positive-area polygons.
    std::optional<core::Index> sidecar_node{};
    std::optional<NastranBulkPosition> field{};
    std::optional<NastranPropertyId> property{};
};

struct ExternalPgThicknessReport {
    RootSnapshotLocation source;
    SourceDigest digest;
    structure::PgReferenceIdentity pg;
    core::PlannedBuffer<ExternalPgThicknessMatch> matches;
    structure::PgReferenceCallInfo work;
};

// Reserve every staged uniform value and four slots per optional original-node
// face. The returned leaf limits retain that reservation while values are built.
[[nodiscard]] core::Result<structure::ExactLimits>
admit_external_pg_thickness_inputs(const structure::PgReferenceLimits&,
                                   std::uint64_t uniform_values, std::uint64_t nodal_faces,
                                   core::ResourcePlan&) noexcept;

// Compare a constant or original-node thickness on the exact common subdivision
// of the certified face and PG. Original triangles use affine interpolation and
// original quadrilaterals use their bilinear chart map without a diagonal split.
// Each declaration applies only on its positive-area region. All inputs are read-only.
[[nodiscard]] core::Result<ExternalPgThicknessReport> compare_external_pg_regional_thickness(
    const ExternalPlanarPg&, const structure::PgReferenceGeometry&,
    std::span<const ExternalPgThicknessClaim>, const ExternalPgThicknessLimits&,
    core::ResourcePlan&, std::optional<ExternalPgThicknessDiagnostic>* = nullptr) noexcept;

// The bounded Gmsh map certifies the same chart coverage. The thickness
// comparison on that chart remains exact and does not spend the geometry budget.
[[nodiscard]] core::Result<ExternalPgThicknessReport> compare_external_pg_regional_thickness(
    const GmshPgBounds&, const structure::PgReferenceGeometry&,
    std::span<const ExternalPgThicknessClaim>, const ExternalPgThicknessLimits&,
    core::ResourcePlan&, std::optional<ExternalPgThicknessDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::io::detail
