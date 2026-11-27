#pragma once

#include <uvlas/structure/pg_global_constraints.hpp>

namespace uvlas::structure {

struct PgSparseChartLimits {
    PgReferenceLimits arithmetic;
    std::uint64_t maximum_pivots{}, maximum_row_terms{}, maximum_factor_terms{};
    bool measure_preparation = false;
};

struct PgSparseChartPivot {
    core::Index original_row{}, dependent_column{}; ///< Original row and global coefficient triple.
};

struct PgSparseChartPreparation {
    std::uint64_t rank{}, original_nonzeros{}, lower_terms{}, upper_terms{};
    std::uint64_t pivot_capacity{}, work_row_capacity{}, factor_term_capacity{};
    PgReferenceCallInfo arithmetic; ///< Cumulative ResourcePlan payload, including scratch.
    std::optional<core::StageTiming> timing; ///< Only present if the caller enabled the clock.
};

struct PgSparseChartIdentity {
    std::array<std::byte, 32> bytes{};
    bool operator==(const PgSparseChartIdentity&) const = default;
};

struct PgSparseChartInput {
    PgSparseChartIdentity identity;
    std::span<const PgCartesianExactPosition> values;
};

struct PgSparseChartDual {
    core::PlannedBuffer<PgCartesianExactPosition> independent_force;
    core::PlannedBuffer<PgCartesianExactPosition> original_multipliers;
    PgReferenceCallInfo call_info;
};

/// Exact affine fixed-pivot chart with sparse factors shared by lift and transpose.
class PgSparseChart {
public:
    PgSparseChart(PgSparseChart&&) noexcept;
    PgSparseChart& operator=(PgSparseChart&&) noexcept;
    PgSparseChart(const PgSparseChart&) = delete;
    PgSparseChart& operator=(const PgSparseChart&) = delete;
    ~PgSparseChart();

    /// Take the original equations only after their fixed RHS is proved compatible.
    [[nodiscard]] static core::Result<PgSparseChart> prepare(PgGlobalConstraints&&,
                                                             const PgSparseChartLimits&,
                                                             core::Workspace&,
                                                             core::ResourcePlan&) noexcept;

    [[nodiscard]] core::Result<PgSparseChartIdentity> identity() const noexcept;
    [[nodiscard]] core::Result<PgSparseChartPreparation> preparation() const noexcept;
    [[nodiscard]] std::span<const PgSparseChartPivot> pivots() const& noexcept;
    [[nodiscard]] std::span<const PgSparseChartPivot> pivots() const&& = delete;
    [[nodiscard]] std::span<const core::Index> free_columns() const& noexcept;
    [[nodiscard]] std::span<const core::Index> free_columns() const&& = delete;
    [[nodiscard]] const PgGlobalConstraints& original() const& noexcept;
    [[nodiscard]] const PgGlobalConstraints& original() const&& = delete;

    /// q=psi(r), using the original stored fixed-support values.
    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    lift(PgSparseChartInput r, const PgReferenceLimits&, core::Workspace&,
         core::ResourcePlan&) const noexcept;
    /// Homogeneous derivative J*r with zero original RHS.
    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    tangent(PgSparseChartInput r, const PgReferenceLimits&, core::Workspace&,
            core::ResourcePlan&) const noexcept;
    /// Check and lift a separately owned exact RHS vector, such as a time coefficient
    /// or derivative, in original row order. All original equations must be compatible.
    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    lift_rhs(PgSparseChartInput r, PgSparseChartInput original_rhs, const PgReferenceLimits&,
             core::Workspace&, core::ResourcePlan&) const noexcept;

    /// Return J^T f and the original-row multipliers for the eliminated components.
    [[nodiscard]] core::Result<PgSparseChartDual> transpose(PgSparseChartInput force,
                                                            const PgReferenceLimits&,
                                                            core::Workspace&,
                                                            core::ResourcePlan&) const noexcept;

private:
    friend class PgComponentSparseChart;
    struct Data;
    PgSparseChart(PgGlobalConstraints&&, core::PlannedBuffer<Data>) noexcept;
    enum class RhsKind : std::uint8_t { Fixed, Zero, Explicit };
    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    lift_impl(PgSparseChartInput, PgSparseChartInput, RhsKind, const PgReferenceLimits&,
              core::Workspace&, core::ResourcePlan&) const noexcept;

    PgGlobalConstraints original_;
    core::PlannedBuffer<Data> data_;
};

struct PgComponentSparseInput {
    PgSparseChartIdentity identity;
    std::span<const ExactRational> values;
};

struct PgComponentSparseDual {
    core::PlannedBuffer<ExactRational> independent_force;
    core::PlannedBuffer<PgCartesianExactPosition> original_multipliers;
    PgReferenceCallInfo call_info;
};

struct PgComponentSparsePreparation {
    std::uint64_t factors{}, scalar_rank{}, lower_terms{}, upper_terms{};
    PgReferenceCallInfo arithmetic;
    std::optional<core::StageTiming> timing;
};

/// Static affine chart for Cartesian rows with component selection. Axes with
/// identical row sets share one exact sparse factor.
class PgComponentSparseChart {
public:
    PgComponentSparseChart(PgComponentSparseChart&&) noexcept;
    PgComponentSparseChart& operator=(PgComponentSparseChart&&) noexcept;
    ~PgComponentSparseChart();

    /// Pivot, row-work and factor-term limits apply to each distinct factor.
    /// The common arithmetic and ResourcePlan limits cover all factors.
    [[nodiscard]] static core::Result<PgComponentSparseChart> prepare(PgGlobalConstraints&&,
                                                                      const PgSparseChartLimits&,
                                                                      core::Workspace&,
                                                                      core::ResourcePlan&) noexcept;

    [[nodiscard]] core::Result<PgSparseChartIdentity> identity() const noexcept;
    [[nodiscard]] core::Result<PgComponentSparsePreparation> preparation() const noexcept;
    [[nodiscard]] const PgGlobalConstraints& original() const& noexcept;
    [[nodiscard]] const PgGlobalConstraints& original() const&& = delete;
    /// Increasing Cartesian-interleaved scalar indices, 3*column + axis.
    [[nodiscard]] std::span<const core::Index> free_columns() const& noexcept;
    [[nodiscard]] std::span<const core::Index> free_columns() const&& = delete;

    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    lift(PgComponentSparseInput, const PgReferenceLimits&, core::Workspace&,
         core::ResourcePlan&) const noexcept;
    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    tangent(PgComponentSparseInput, const PgReferenceLimits&, core::Workspace&,
            core::ResourcePlan&) const noexcept;
    [[nodiscard]] core::Result<PgComponentSparseDual> transpose(PgSparseChartInput force,
                                                                const PgReferenceLimits&,
                                                                core::Workspace&,
                                                                core::ResourcePlan&) const noexcept;

private:
    struct Data;
    PgComponentSparseChart(PgGlobalConstraints&&, core::PlannedBuffer<Data>) noexcept;
    [[nodiscard]] core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
    lift_impl(PgComponentSparseInput, bool homogeneous, const PgReferenceLimits&, core::Workspace&,
              core::ResourcePlan&) const noexcept;

    PgGlobalConstraints original_;
    core::PlannedBuffer<Data> data_;
};

} // namespace uvlas::structure
