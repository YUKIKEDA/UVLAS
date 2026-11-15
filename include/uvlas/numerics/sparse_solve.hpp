#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <uvlas/numerics/linear_solve.hpp>

namespace uvlas::numerics::linear {

/// Borrow one canonical real 64-bit CSC operator for synchronous preparation/solve.
/// Values, indices, their owner and context remain unchanged during the call.
class SparseOperator {
public:
    /// Bind a live source; preparation copies its actual coefficients and structure.
    SparseOperator(const SparseReal& values, Context context) noexcept
        : values_(values), context_(context) {}
    /// A temporary cannot supply the original-operator borrow.
    SparseOperator(SparseReal&&, Context) = delete;
    /// A const temporary cannot supply the original-operator borrow.
    SparseOperator(const SparseReal&&, Context) = delete;
    /// Borrow the original source while its owning matrix remains alive and unmoved.
    [[nodiscard]] const SparseReal& values() const noexcept {
        return values_;
    }
    /// Copy model, graph, matrix, scale, precision and physical evaluation identities.
    [[nodiscard]] Context context() const noexcept {
        return context_;
    }

private:
    const SparseReal& values_;
    Context context_;
};

/// Which property was actually proved for the complete original numerical matrix.
enum class SparseSpdStatus : std::uint8_t {
    PositiveDefinite,  ///< A rigorous shifted-factor qualification succeeded.
    NegativeCurvature, ///< An owned nonzero witness has strictly negative quadratic upper bound.
    ZeroMode,          ///< An owned nonzero witness has provably exact zero original matrix action.
    Uncertified,       ///< Available finite checks did not establish any of the preceding claims.
    Asymmetric         ///< Opposite original coefficients differ; they were not averaged.
};

/// Owns matrix-property diagnostics without claiming that a linear system was solved.
/// This report is not accepted as proof for another solve; preparation checks its actual input.
class SparseSpdDiagnostic {
public:
    /// Transfer any owned witness.
    SparseSpdDiagnostic(SparseSpdDiagnostic&&) noexcept = default;
    /// Replace diagnostics and invalidate former witness borrows.
    SparseSpdDiagnostic& operator=(SparseSpdDiagnostic&&) noexcept = default;
    /// Owned witnesses cannot be implicitly copied.
    SparseSpdDiagnostic(const SparseSpdDiagnostic&) = delete;
    /// Copy assignment is forbidden.
    SparseSpdDiagnostic& operator=(const SparseSpdDiagnostic&) = delete;
    /// Release witness storage independently of the input/workspace lifetime.
    ~SparseSpdDiagnostic() = default;
    /// Distinct proved or unresolved matrix property.
    [[nodiscard]] SparseSpdStatus status() const noexcept {
        return status_;
    }
    /// Context of the original matrix actually inspected.
    [[nodiscard]] Context context() const noexcept {
        return context_;
    }
    /// Conservative qualification shift, when that check was evaluated.
    [[nodiscard]] std::optional<double> qualification_shift() const noexcept {
        return shift_;
    }
    /// Upper bound on w^T*A*w for a certified negative witness; otherwise absent.
    [[nodiscard]] std::optional<double> quadratic_upper() const noexcept {
        return quadratic_;
    }
    /// Borrow a nonzero n-by-1 witness in original coordinates, or nullptr if none was proved.
    [[nodiscard]] const DenseReal* witness() const& noexcept {
        return witness_ ? &*witness_ : nullptr;
    }
    /// A temporary cannot provide a borrowed witness.
    [[nodiscard]] const DenseReal* witness() const&& = delete;

private:
    friend core::Result<SparseSpdDiagnostic>
    inspect_sparse_spd(const SparseOperator&, const Limits&, core::ResourcePlan&) noexcept;
    SparseSpdDiagnostic(SparseSpdStatus status, Context context, std::optional<double> shift,
                        std::optional<double> quadratic, std::optional<DenseReal> witness) noexcept
        : status_(status),
          context_(context),
          shift_(shift),
          quadratic_(quadratic),
          witness_(std::move(witness)) {}
    SparseSpdStatus status_;
    Context context_;
    std::optional<double> shift_;
    std::optional<double> quadratic_;
    std::optional<DenseReal> witness_;
};

/// Inspect exact symmetry, rigorous SPD qualification and independently checked witnesses.
/// A failed shifted Cholesky is Uncertified, not a proof of negative curvature.
/// Matrix/limits/representation/resource errors return Failure; the original never changes.
[[nodiscard]] core::Result<SparseSpdDiagnostic>
inspect_sparse_spd(const SparseOperator& original, const Limits& limits,
                   core::ResourcePlan& result_plan) noexcept;

/// Explicit factorization policy; no regularization or pseudoinverse is selected implicitly.
enum class SparseMethod : std::uint8_t {
    GeneralLu,          ///< COLAMD SparseLU for indefinite KKT and general dynamic tangents.
    PositiveDefiniteLlt ///< Inspect and prove SPD before factoring the original with sparse LLT.
};

namespace detail {
class SparsePreparedData;
}

/// Owns an immutable original snapshot, validated sparse factors and exclusive prepared scratch.
/// Factors and permutations are extracted at preparation; inner triangular solves allocate nothing.
class SparsePrepared {
public:
    /// Prepare original CSC, positive finite scales with SI units, factors and RHS capacity.
    /// Scale shape/context must match the operator; values and units are snapshotted and
    /// compared at every solve, so reusing a generation label cannot replace actual scales.
    /// Plan accounts UVLAS-owned arrays; external backend forecast/proof admission is separate.
    /// A general stock Eigen factor has no proven external allocation ceiling.
    [[nodiscard]] static core::Result<SparsePrepared> prepare(const SparseOperator&, const Scales&,
                                                              SparseMethod, const Limits&,
                                                              core::ResourcePlan&) noexcept;
    /// Transfer prepared ownership; the source becomes inactive.
    SparsePrepared(SparsePrepared&&) noexcept;
    /// Replace ownership, invalidating internal source borrows.
    SparsePrepared& operator=(SparsePrepared&&) noexcept;
    /// Prepared scratch cannot be shared by copying.
    SparsePrepared(const SparsePrepared&) = delete;
    /// Copy assignment is forbidden.
    SparsePrepared& operator=(const SparsePrepared&) = delete;
    /// Release factors independently of prior returned solutions.
    ~SparsePrepared();
    /// Whether this owner still contains a usable prepared backend.
    [[nodiscard]] bool active() const noexcept;
    /// Original context labels; inactive moved-from owners cannot solve even when these match.
    [[nodiscard]] Context context() const noexcept {
        return context_;
    }
    /// Explicit factorization policy used at preparation.
    [[nodiscard]] SparseMethod method() const noexcept {
        return method_;
    }

private:
    friend core::Result<Solution<double>>
    solve_linear(const SparseOperator&, const Rhs<double>&, const Scales&, const Precision<double>&,
                 const Accuracy&, const Limits&, SparsePrepared&, core::ResourcePlan&) noexcept;
    SparsePrepared(std::unique_ptr<detail::SparsePreparedData>, Context, SparseMethod) noexcept;
    std::unique_ptr<detail::SparsePreparedData> data_;
    Context context_;
    SparseMethod method_;
};

/// Solve the actual original system with prepared factors and fresh independent qualification.
/// Full generations AND actual CSC contents must match. Busy/moved-from scratch is rejected.
/// Owns returned values/diagnostics; inputs and previous solutions remain unchanged on failure.
/// Both RHS planes are retained. Stalled binary64 refinement promotes the owning solution
/// to two-word values; an explicit rounded projection is independently requalified.
[[nodiscard]] core::Result<Solution<double>> solve_linear(const SparseOperator&, const Rhs<double>&,
                                                          const Scales&, const Precision<double>&,
                                                          const Accuracy&, const Limits&,
                                                          SparsePrepared&,
                                                          core::ResourcePlan& result_plan) noexcept;

} // namespace uvlas::numerics::linear
