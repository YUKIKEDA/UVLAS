#pragma once
#include <span>
#include <uvlas/numerics/linear_solve.hpp>

namespace uvlas::numerics::detail {
class CheckedEft;
}

namespace uvlas::numerics::linear::detail {
[[nodiscard]] core::Result<void> check_accuracy_inputs(const Accuracy& accuracy) noexcept;
[[nodiscard]] core::Result<void>
check_consumer_inputs(std::span<const ConsumerCheck> checks) noexcept;
[[nodiscard]] core::Result<void> check_context_inputs(const Context& context) noexcept;
[[nodiscard]] core::Result<void> check_scale_inputs(core::Index rows, const Scales& scales,
                                                    const Context& expected) noexcept;
[[nodiscard]] core::Result<void> check_capacity_inputs(core::Index rows, core::Index rhs_columns,
                                                       core::Index update_rank,
                                                       const Limits& limits) noexcept;
template <StorageScalar Scalar>
[[nodiscard]] core::Result<void>
check_solve_inputs(core::Index rows, core::Index update_rank, const Context& expected,
                   const Rhs<Scalar>& rhs, const Scales& scales, const Precision<Scalar>& precision,
                   const Accuracy& accuracy, const Limits& limits) noexcept;

struct MagnitudeBounds {
    double lower;
    double upper;
};
// Outward arithmetic for already-validated nonnegative magnitude intervals.
// A finite overflow is represented by max-double lower / infinity upper.
[[nodiscard]] MagnitudeBounds add_magnitudes(MagnitudeBounds a, MagnitudeBounds b) noexcept;
[[nodiscard]] MagnitudeBounds multiply_magnitudes(MagnitudeBounds a, MagnitudeBounds b) noexcept;
/// Accumulate one product of nonnegative lower bounds with the same downward
/// rounding as the lower component of multiply_magnitudes and add_magnitudes.
[[nodiscard]] double
add_magnitude_product_lower(double sum, double a, double b,
                            const uvlas::numerics::detail::CheckedEft&) noexcept;

// Enclose abs(high+low), not abs(high)+abs(low). Nonfinite words return Nonfinite.
// A finite exact magnitude outside binary64 can have an infinite upper bound.
template <StorageScalar Scalar>
[[nodiscard]] core::Result<MagnitudeBounds> actual_magnitude_bounds(Scalar high,
                                                                    Scalar low) noexcept;
// All spans describe one packed logical column, with rows elements. rhs_lo alone
// may be empty (a binary64 RHS). Scales are positive finite; the same scales are
// broadcast across RHS columns by the caller. Error bounds and denominator lower
// bounds must be nonnegative. Physical words and scales reject nonfinite inputs.
// A valid but unmet criterion is returned as a diagnostic metric, not as success
// qualification. An unrepresentable conservative error upper is PrecisionUnresolved.
template <StorageScalar Scalar>
[[nodiscard]] core::Result<ColumnMetrics>
evaluate_column_metrics(core::Index rows, std::span<const Scalar> residual_estimate,
                        std::span<const double> absolute_error_upper,
                        std::span<const double> backward_denominator_lower,
                        std::span<const Scalar> rhs_hi, std::span<const Scalar> rhs_lo,
                        std::span<const double> physical_residual_scales,
                        const Accuracy& accuracy) noexcept;
// One common acceptance criterion: every column needs residual norm <=1 and its
// independently computed backward-error upper <= the requested backward budget.
[[nodiscard]] bool meets_accuracy(ColumnMetrics metrics, const Accuracy& accuracy) noexcept;

// Prepared before entering solve/IR kernels. All result arrays have their own
// ResourcePlan allocations and cannot borrow the prepared factor/workspace.
struct MetadataBuffers {
    core::PlannedBuffer<ColumnMetrics> retained_columns;
    core::PlannedBuffer<ColumnMetrics> returned_columns;
    core::PlannedBuffer<double> unknown_scales;
    core::PlannedBuffer<double> residual_scales;
    core::PlannedBuffer<UnitDimensions> unknown_units;
    core::PlannedBuffer<UnitDimensions> residual_units;
    core::PlannedBuffer<std::optional<ConsumerProof>> consumers;
};
struct Finalization {
    Finalization(Context identity, core::RhsGeneration rhs_identity,
                 core::NumericalGeneration rhs_values, Accuracy request, BackendAdmission prepared,
                 BackendAdmission current, unsigned corrections, RetentionPath path,
                 ConditionDiagnostic diagnostic,
                 LinearSystemOrientation system = LinearSystemOrientation::Original) noexcept
        : context(identity),
          rhs(rhs_identity),
          rhs_buffer_generation(rhs_values),
          accuracy(request),
          preparation_admission(prepared),
          solve_admission(current),
          refinements(corrections),
          retention(path),
          condition(diagnostic),
          orientation(system) {}
    Context context;
    core::RhsGeneration rhs;
    core::NumericalGeneration rhs_buffer_generation;
    Accuracy accuracy;
    BackendAdmission preparation_admission;
    BackendAdmission solve_admission;
    unsigned refinements;
    RetentionPath retention;
    ConditionDiagnostic condition;
    LinearSystemOrientation orientation;
};
class SolutionFactory {
public:
    // Validate and copy all scale/unit metadata, validate necessary requirements,
    // and preallocate every diagnostic/proof slot. Called before the inner kernel.
    [[nodiscard]] static core::Result<MetadataBuffers>
    prepare_metadata(core::Index rows, core::Index rhs_columns, const Scales& scales,
                     std::span<const ConsumerCheck> checks,
                     core::ResourcePlan& result_plan) noexcept;
    // Backends fill per-column metrics from evaluate_column_metrics, then transfer
    // the exact immutable candidates here. This entry validates numerical gates,
    // executes each required evaluator on the exact retained/returned values,
    // checks certificate context/budget/representation, and publishes only achievement.
    // With no projection the returned metrics are copied from retained metrics.
    // Failure consumes only these local candidates/metadata, never an input parent.
    template <StorageScalar Scalar>
    [[nodiscard]] static core::Result<Solution<Scalar>>
    finalize(PrecisionBuffer<Scalar> retained, std::optional<PrecisionBuffer<Scalar>> projected,
             const Precision<Scalar>& precision, MetadataBuffers metadata,
             Finalization context) noexcept;
};

} // namespace uvlas::numerics::linear::detail
