#pragma once

#include <uvlas/numerics/storage.hpp>
#include <uvlas/structure/pg_full_trace.hpp>

namespace uvlas::structure {

namespace detail {
struct PgConstraintTag;
}
using PgConstraintId = core::Identifier<detail::PgConstraintTag>;

struct PgConstraintPatchBinding {
    PgConstraintPatchBinding(numerics::PatchId physical_id,
                             const PgPreparedPatch& reference) noexcept
        : id(physical_id), patch(reference) {}

    numerics::PatchId id;
    const PgPreparedPatch& patch;
};

/// Endpoints index the patch binding list supplied to the same preparation.
/// The ID is unique within the trace kind, preserving separate Joint and Support namespaces.
struct PgConstraintTraceBinding {
    PgConstraintTraceBinding(PgConstraintId physical_id, const PgFullTrace& constraint,
                             core::Index left, std::optional<core::Index> right) noexcept
        : id(physical_id), trace(constraint), left_patch(left), right_patch(right) {}

    PgConstraintId id;
    const PgFullTrace& trace;
    core::Index left_patch;
    std::optional<core::Index> right_patch; ///< Empty for a fixed support.
};

/// A static material point row in one prepared basis, in increasing local column order.
struct PgConstraintPointBinding {
    PgConstraintId id;
    PgPreparedPatchIdentity basis;
    core::Index patch;
    std::uint8_t world_components; ///< Bit 0/1/2 selects world X/Y/Z.
    std::span<const PgTraceSparseEntry> terms;
    const PgCartesianExactPosition& prescribed;
};

struct PgConstraintPatchLayout {
    numerics::PatchId id{0};
    PgPreparedPatchIdentity basis;
    core::Index first{}, count{}; ///< Contiguous coefficient triples; x/y/z interleave.
    bool operator==(const PgConstraintPatchLayout&) const = default;
};

struct PgConstraintSource {
    PgConstraintId id{0};
    PgTraceKind kind = PgTraceKind::Connection;
    std::uint8_t world_components = 7;
    numerics::PatchId left{0};
    std::optional<numerics::PatchId> right;
    core::Index first_row{}, rows{};
    core::Index first_overlap{}, overlaps{};
};

struct PgConstraintOriginalRow {
    core::Index source{}, source_row{}, source_overlap{}; ///< A point has no overlap (-1).
    PgCoefficientField field = PgCoefficientField::R;
    std::uint32_t power = 0;
    core::Index local_representative = -1;
    ExactRational local_scale;
    PgCartesianExactPosition prescribed; ///< Original right-hand side c.
};

struct PgGlobalConstraintLimits {
    PgReferenceLimits arithmetic;
    std::uint64_t maximum_patches{}, maximum_sources{}, maximum_rows{}, maximum_nonzeros{};
};

struct PgGlobalConstraintInfo {
    core::Index columns{}, rows{}, nonzeros{}; ///< One coefficient per participating triple.
    bool common_cartesian = true;              ///< All rows act on three components.
    PgReferenceCallInfo preparation;
};

struct PgGlobalConstraintCoordinates {
    std::span<const PgConstraintPatchLayout> layout;
    std::span<const PgCartesianExactPosition> values;
};

/// Original exact equation Cq=c. The CSR stores one coefficient for each identical
/// Cartesian action. A flat scalar coordinate of triple j and component a is 3*j+a.
class PgGlobalConstraints {
public:
    PgGlobalConstraints(PgGlobalConstraints&&) noexcept;
    PgGlobalConstraints& operator=(PgGlobalConstraints&&) noexcept;
    PgGlobalConstraints(const PgGlobalConstraints&) = delete;
    PgGlobalConstraints& operator=(const PgGlobalConstraints&) = delete;
    ~PgGlobalConstraints();

    /// Preserve input patch and constraint order and own all original rows, including zero rows.
    [[nodiscard]] static core::Result<PgGlobalConstraints>
    prepare(std::span<const PgConstraintPatchBinding>, std::span<const PgConstraintTraceBinding>,
            const PgGlobalConstraintLimits&, core::Workspace&, core::ResourcePlan&) noexcept;

    /// Append fixed point rows after the existing traces. Point and trace support
    /// IDs belong to the same namespace. A point owns only its selected components.
    [[nodiscard]] static core::Result<PgGlobalConstraints>
    prepare(std::span<const PgConstraintPatchBinding>, std::span<const PgConstraintTraceBinding>,
            std::span<const PgConstraintPointBinding>, const PgGlobalConstraintLimits&,
            core::Workspace&, core::ResourcePlan&) noexcept;

    [[nodiscard]] core::Result<PgGlobalConstraintInfo> info() const noexcept;
    [[nodiscard]] std::span<const PgConstraintPatchLayout> layout() const& noexcept;
    [[nodiscard]] std::span<const PgConstraintPatchLayout> layout() const&& = delete;
    [[nodiscard]] std::span<const PgConstraintSource> sources() const& noexcept;
    [[nodiscard]] std::span<const PgConstraintSource> sources() const&& = delete;
    [[nodiscard]] std::span<const PgTraceOverlap> overlaps() const& noexcept;
    [[nodiscard]] std::span<const PgTraceOverlap> overlaps() const&& = delete;
    [[nodiscard]] std::span<const PgConstraintOriginalRow> original_rows() const& noexcept;
    [[nodiscard]] std::span<const PgConstraintOriginalRow> original_rows() const&& = delete;
    [[nodiscard]] core::Result<std::span<const PgTraceSparseEntry>>
        row(core::Index) const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgTraceSparseEntry>>
        row(core::Index) const&& = delete;

    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    apply(PgGlobalConstraintCoordinates, const PgReferenceLimits&, core::Workspace&,
          core::ResourcePlan&) const noexcept;
    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    residual(PgGlobalConstraintCoordinates, const PgReferenceLimits&, core::Workspace&,
             core::ResourcePlan&) const noexcept;
    /// Lambda uses the original row order of this exact owner. Unselected components
    /// of point rows contribute zero, as in the forward action.
    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    transpose(std::span<const PgCartesianExactPosition>, const PgReferenceLimits&, core::Workspace&,
              core::ResourcePlan&) const noexcept;

private:
    struct Data;
    explicit PgGlobalConstraints(core::PlannedBuffer<Data>) noexcept;
    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    apply_impl(PgGlobalConstraintCoordinates, bool subtract_prescribed, const PgReferenceLimits&,
               core::Workspace&, core::ResourcePlan&) const noexcept;

    core::PlannedBuffer<Data> data_;
};

} // namespace uvlas::structure
