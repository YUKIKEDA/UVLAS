#pragma once
#include <memory>
#include <optional>
#include <uvlas/numerics/linear_solve.hpp>

namespace uvlas::numerics::linear {
namespace detail {
struct DenseData;
struct UpdatedData;
} // namespace detail

/// Prepared real pivoted-LU factorization, immutable original coefficients and reusable workspace.
/// One solve may use this object at a time; concurrent/reentrant use is rejected.
class DensePrepared {
public:
    /// Transfer preparation ownership; the source becomes inactive.
    DensePrepared(DensePrepared&&) noexcept;
    /// Replace preparation ownership; previous borrowed diagnostics expire.
    DensePrepared& operator=(DensePrepared&&) noexcept;
    /// A factorization/workspace cannot be duplicated implicitly.
    DensePrepared(const DensePrepared&) = delete;
    /// Copy assignment is forbidden.
    DensePrepared& operator=(const DensePrepared&) = delete;
    /// Release preparation independently of previously returned solutions.
    ~DensePrepared();
    /// Validate original square coefficients, scale/context/capacity and backend admission,
    /// then copy the original operator and construct a binary64 partial-pivot LU.
    /// Owned arrays use workspace_plan; external Eigen payload follows its separate plan.
    /// All backend exceptions become Failure. No input is modified.
    [[nodiscard]] static core::Result<DensePrepared>
    prepare(const DenseOperator& op, const Scales& scales, const Limits& limits,
            core::ResourcePlan& workspace_plan) noexcept;
    /// Whether this object still owns preparation after moving.
    [[nodiscard]] bool active() const noexcept {
        return data_ != nullptr;
    }
    /// Copy the declared external-backend plan, or no value after moving.
    [[nodiscard]] std::optional<BackendMemoryPlan> memory_plan() const noexcept;

private:
    friend core::Result<Solution<double>>
    solve_linear(const DenseOperator&, const Rhs<double>&, const Scales&, const Precision<double>&,
                 const Accuracy&, const Limits&, DensePrepared&, core::ResourcePlan&) noexcept;
    explicit DensePrepared(std::unique_ptr<detail::DenseData> data) noexcept;
    std::unique_ptr<detail::DenseData> data_;
};

/// Prepared real-base LU plus unconjugated complex Woodbury/Schur update.
/// Owns original A0/U/V snapshots and reusable workspaces; returns no scratch borrow.
class UpdatedPrepared {
public:
    /// Transfer preparation ownership; the source becomes inactive.
    UpdatedPrepared(UpdatedPrepared&&) noexcept;
    /// Replace preparation ownership.
    UpdatedPrepared& operator=(UpdatedPrepared&&) noexcept;
    /// Copying is forbidden.
    UpdatedPrepared(const UpdatedPrepared&) = delete;
    /// Copy assignment is forbidden.
    UpdatedPrepared& operator=(const UpdatedPrepared&) = delete;
    /// Release factors independently of returned solutions.
    ~UpdatedPrepared();
    /// Prepare A0's real LU, Y=A0^-1 U and S=I+V^T Y without conjugation.
    /// The original effective operator is the fixed-order rounded coefficient matrix;
    /// qualification never substitutes the staged low-rank action for its residual.
    [[nodiscard]] static core::Result<UpdatedPrepared>
    prepare(const ColumnUpdatedOperator& op, const Scales& scales, const Limits& limits,
            core::ResourcePlan& workspace_plan) noexcept;
    /// Whether this object still owns its preparation.
    [[nodiscard]] bool active() const noexcept {
        return data_ != nullptr;
    }
    /// Copy the declared external-backend plan, or no value after moving.
    [[nodiscard]] std::optional<BackendMemoryPlan> memory_plan() const noexcept;

    /// Uncertified reciprocal kappa_1 estimate of the unscaled small S only.
    /// Uses Eigen 5 complex PartialPivLU::rcond; never a condition estimate for A_eff.
    /// A moved object, rank-zero update or unavailable estimate returns NotEvaluated.
    [[nodiscard]] ConditionDiagnostic schur_condition() const noexcept;

private:
    friend core::Result<Solution<std::complex<double>>>
    solve_linear(const ColumnUpdatedOperator&, const Rhs<std::complex<double>>&, const Scales&,
                 const Precision<std::complex<double>>&, const Accuracy&, const Limits&,
                 UpdatedPrepared&, core::ResourcePlan&) noexcept;
    explicit UpdatedPrepared(std::unique_ptr<detail::UpdatedData> data) noexcept;
    std::unique_ptr<detail::UpdatedData> data_;
};
} // namespace uvlas::numerics::linear
