#pragma once

#include <cstddef>
#include <optional>
#include <uvlas/structure/pg_sparse_chart.hpp>

namespace uvlas::structure {

struct PgRigidSupportLimits {
    PgReferenceLimits arithmetic;
    /// Per generator check. All work shares the arithmetic cancellation token.
    PgReferenceLimits chart_verification;
    std::uint64_t maximum_original_rows{};
};

struct PgRigidSupportRow {
    std::optional<std::size_t> group{};
    PgCartesianExactPosition reference_position; ///< C X in world metres.
    ExactRational translation; ///< C b, the weight of a constant Cartesian translation.
};

/// Associate a physical support with a dense group ordinal.
struct PgRigidSupportBinding {
    PgConstraintId support{0};
    std::size_t group{};
};

/// Original RHS maps qualified separately for the affine span of each rigid motion.
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

    /// Each group is checked independently with zero entries on all other original rows.
    [[nodiscard]] static core::Result<PgRigidSupportRows>
    prepare_groups(PgSparseChart&&, std::span<const PgConstraintPatchBinding>,
                   std::span<const PgRigidSupportBinding>, std::size_t group_count,
                   const PgRigidSupportLimits&, core::Workspace&, core::ResourcePlan&) noexcept;

    [[nodiscard]] std::size_t group_count() const noexcept {
        return group_count_;
    }
    [[nodiscard]] const PgSparseChart& chart() const& noexcept;
    [[nodiscard]] const PgSparseChart& chart() const&& = delete;
    [[nodiscard]] std::span<const PgRigidSupportRow> rows() const& noexcept;
    [[nodiscard]] std::span<const PgRigidSupportRow> rows() const&& = delete;

private:
    PgRigidSupportRows(PgSparseChart&&, core::PlannedBuffer<PgRigidSupportRow>,
                       std::size_t group_count) noexcept;
    PgSparseChart chart_;
    core::PlannedBuffer<PgRigidSupportRow> rows_;
    std::size_t group_count_;
};

} // namespace uvlas::structure
