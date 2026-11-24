#pragma once

#include <uvlas/structure/pg_prepared_patch.hpp>

namespace uvlas::structure {

/// Part of an external edge. Along coordinates use the canonical vertex order.
struct PgTraceEdgeSpan {
    core::Index edge{};
    const ExactRational& first;
    const ExactRational& last;
    const ExactRational& along_first;
    const ExactRational& along_last;
    std::optional<core::Index> layer; ///< Empty selects the entire physical thickness.
};

struct PgTraceLimits {
    PgReferenceLimits arithmetic;
    std::uint64_t maximum_nonzeros = 0;
};

struct PgTraceDepth {
    ExactRational lower_first, upper_first, lower_last, upper_last;
};

/// One common planar interval, with physical depths affine in its local parameter.
struct PgTraceOverlap {
    core::Index left_edge{}, right_edge = -1;
    std::optional<core::Index> left_layer, right_layer;
    ExactRational first, last;
    ExactRational active_first, active_last; ///< Local [0,1] interval with positive overlap.
    PgTraceDepth left_depth, right_depth;
};

struct PgTraceOriginalRow {
    core::Index overlap{};
    PgCoefficientField field = PgCoefficientField::R;
    std::uint32_t power = 0;
    core::Index unique_row = -1;
    ExactRational scale; ///< Original row = scale * unique row. Zero has no representative.
    PgCartesianExactPosition prescribed;
};

struct PgTraceSparseEntry {
    core::Index column{}; ///< Left basis followed by right basis.
    ExactRational coefficient;
};

enum class PgTraceKind : std::uint8_t { Connection, Support };

struct PgTraceInfo {
    PgTraceKind kind = PgTraceKind::Connection;
    PgPreparedPatchIdentity left, right;
    core::Index left_columns{}, right_columns{};
    std::size_t overlaps = 0, original_rows = 0, unique_rows = 0, nonzeros = 0;
    PgReferenceCallInfo preparation;
};

struct PgTraceCoordinates {
    PgPreparedPatchIdentity left, right;
    std::span<const PgCartesianExactPosition> left_values, right_values;
};

/// Owns the original homogeneous full-position connection and its row equivalence.
class PgFullTrace {
public:
    PgFullTrace(PgFullTrace&&) noexcept = default;
    PgFullTrace& operator=(PgFullTrace&&) noexcept = default;
    PgFullTrace(const PgFullTrace&) = delete;
    PgFullTrace& operator=(const PgFullTrace&) = delete;

    /// Both lists tile the same material interval. Only positive physical depth is connected.
    [[nodiscard]] static core::Result<PgFullTrace>
    prepare(const PgPreparedPatch& left, std::span<const PgTraceEdgeSpan> left_spans,
            const PgPreparedPatch& right, std::span<const PgTraceEdgeSpan> right_spans,
            const PgTraceLimits&, core::Workspace&, core::ResourcePlan&) noexcept;

    /// Freeze a compatible nonzero prescribed position in this same prepared basis.
    [[nodiscard]] static core::Result<PgFullTrace>
    prepare_support(const PgPreparedPatch&, std::span<const PgTraceEdgeSpan>,
                    PgPreparedPatchIdentity, std::span<const PgCartesianExactPosition> prescribed,
                    const PgTraceLimits&, core::Workspace&, core::ResourcePlan&) noexcept;

    [[nodiscard]] core::Result<PgTraceInfo> info() const noexcept;
    [[nodiscard]] std::span<const PgTraceOverlap> overlaps() const& noexcept;
    [[nodiscard]] std::span<const PgTraceOverlap> overlaps() const&& = delete;
    [[nodiscard]] std::span<const PgTraceOriginalRow> original_rows() const& noexcept;
    [[nodiscard]] std::span<const PgTraceOriginalRow> original_rows() const&& = delete;
    [[nodiscard]] core::Result<std::span<const PgTraceSparseEntry>>
        unique_row(core::Index) const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgTraceSparseEntry>>
        unique_row(core::Index) const&& = delete;

    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    apply(PgTraceCoordinates, const PgReferenceLimits&, core::Workspace&,
          core::ResourcePlan&) const noexcept;
    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    residual(PgTraceCoordinates, const PgReferenceLimits&, core::Workspace&,
             core::ResourcePlan&) const noexcept;
    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    transpose(std::span<const PgCartesianExactPosition>, const PgReferenceLimits&, core::Workspace&,
              core::ResourcePlan&) const noexcept;
    /// Put each normalized dual on its first original representative.
    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    restore_dual(std::span<const PgCartesianExactPosition>, const PgReferenceLimits&,
                 core::Workspace&, core::ResourcePlan&) const noexcept;

private:
    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    apply_impl(PgTraceCoordinates, bool subtract_prescribed, const PgReferenceLimits&,
               core::Workspace&, core::ResourcePlan&) const noexcept;

    PgFullTrace(core::PlannedBuffer<PgTraceOverlap>, core::PlannedBuffer<PgTraceOriginalRow>,
                core::PlannedBuffer<std::size_t>, core::PlannedBuffer<PgTraceSparseEntry>,
                PgTraceInfo) noexcept;

    core::PlannedBuffer<PgTraceOverlap> overlaps_;
    core::PlannedBuffer<PgTraceOriginalRow> rows_;
    core::PlannedBuffer<std::size_t> offsets_;
    core::PlannedBuffer<PgTraceSparseEntry> coefficients_;
    PgTraceInfo info_;
};

} // namespace uvlas::structure
