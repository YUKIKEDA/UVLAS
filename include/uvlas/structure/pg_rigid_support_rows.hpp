#pragma once

#include <uvlas/structure/pg_sparse_chart.hpp>

namespace uvlas::structure {

struct PgRigidSupportLimits {
    PgReferenceLimits arithmetic;
    /// Per generator check. All work shares the arithmetic cancellation token.
    PgReferenceLimits chart_verification;
    std::uint64_t maximum_original_rows{};
};

struct PgRigidSupportRow {
    bool moving = false;
    PgCartesianExactPosition reference_position; ///< C X in world metres.
    ExactRational translation; ///< C b, the weight of a constant Cartesian translation.
};

/// Original RHS maps for a single common rigid motion, qualified on its full affine span.
class PgRigidSupportRows {
public:
    PgRigidSupportRows(PgRigidSupportRows&&) noexcept = default;
    PgRigidSupportRows& operator=(PgRigidSupportRows&&) noexcept = default;
    PgRigidSupportRows(const PgRigidSupportRows&) = delete;
    PgRigidSupportRows& operator=(const PgRigidSupportRows&) = delete;

    /// Unlisted supports and connections keep their original constant RHS.
    /// Moves the chart only after compatibility of all generators is verified.
    [[nodiscard]] static core::Result<PgRigidSupportRows>
    prepare(PgSparseChart&&, std::span<const PgConstraintPatchBinding>,
            std::span<const PgConstraintId> moving, const PgRigidSupportLimits&, core::Workspace&,
            core::ResourcePlan&) noexcept;

    [[nodiscard]] const PgSparseChart& chart() const& noexcept;
    [[nodiscard]] const PgSparseChart& chart() const&& = delete;
    [[nodiscard]] std::span<const PgRigidSupportRow> rows() const& noexcept;
    [[nodiscard]] std::span<const PgRigidSupportRow> rows() const&& = delete;

private:
    PgRigidSupportRows(PgSparseChart&&, core::PlannedBuffer<PgRigidSupportRow>) noexcept;
    PgSparseChart chart_;
    core::PlannedBuffer<PgRigidSupportRow> rows_;
};

} // namespace uvlas::structure
