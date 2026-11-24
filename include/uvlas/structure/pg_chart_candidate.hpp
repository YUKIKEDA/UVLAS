#pragma once

#include <utility>
#include <uvlas/numerics/sparse_solve.hpp>
#include <uvlas/structure/pg_sparse_chart.hpp>
#include <variant>

namespace uvlas::structure {

/// Explicit provisional coefficient tolerance, before production physical qualification.
struct PgChartCandidatePreparation {
    PgReferenceLimits projection;
    double maximum_absolute_coefficient_error{};
    numerics::linear::Limits solver;
};

struct PgChartCandidateSolve {
    numerics::linear::Scales scales;
    numerics::linear::Precision<double> precision;
    numerics::linear::Accuracy accuracy;
    numerics::linear::Limits solver;
    core::RhsGeneration rhs_generation;
    core::NumericalGeneration rhs_buffer_generation, output_generation;
    double maximum_absolute_original_defect;
    std::uint64_t maximum_visits;
    std::stop_token cancellation;
};

/// Each matrix has rows in the owner layout and three columns ordered x, y, z.
struct PgChartCandidateInput {
    PgSparseChartIdentity identity;
    const numerics::PrecisionBuffer<double>& value;
};

struct PgChartCandidateLift {
    PgSparseChartIdentity identity;
    numerics::PrecisionBuffer<double> coordinates;
    double original_residual_upper;
};

struct PgChartCandidateDual {
    PgSparseChartIdentity identity;
    numerics::PrecisionBuffer<double> independent_force, original_multipliers;
    double original_balance_upper;
};

/// Bounds original rational arithmetic across one candidate and all its corrections.
struct PgChartRefinement {
    PgReferenceLimits original_arithmetic;
    std::uint64_t maximum_corrections{};
};

template <class Output>
struct PgChartRefined {
    Output output;
    PgReferenceCallInfo original_work;
    std::uint64_t corrections{}; ///< Additional solves using the original prepared factor.
};

/// Fixed sparse projection candidate plus certificates in the original coefficient layout.
class PgChartCandidate {
public:
    PgChartCandidate(PgChartCandidate&&) noexcept;
    PgChartCandidate& operator=(PgChartCandidate&&) noexcept;
    PgChartCandidate(const PgChartCandidate&) = delete;
    PgChartCandidate& operator=(const PgChartCandidate&) = delete;
    ~PgChartCandidate();

    /// Ownership transfers only after coefficient qualification and factorization succeed.
    [[nodiscard]] static core::Result<PgChartCandidate>
    prepare(PgSparseChart&&, const numerics::linear::Scales& forward,
            const numerics::linear::Scales& transpose, const PgChartCandidatePreparation&,
            core::ResourcePlan&) noexcept;

    [[nodiscard]] const PgSparseChart& exact() const& noexcept;
    [[nodiscard]] const PgSparseChart& exact() const&& = delete;

    [[nodiscard]] core::Result<PgChartCandidateLift>
    lift(PgChartCandidateInput independent, PgChartCandidateInput original_rhs,
         const PgChartCandidateSolve&, core::Workspace&, core::ResourcePlan&) const noexcept;

    [[nodiscard]] core::Result<PgChartCandidateDual> transpose(PgChartCandidateInput force,
                                                               const PgChartCandidateSolve&,
                                                               core::Workspace&,
                                                               core::ResourcePlan&) const noexcept;

    /// Qualify a retained TwoWord result directly with the original sparse rational rows.
    /// The inner solve requires retained precision and algebraic checks only.
    [[nodiscard]] core::Result<PgChartRefined<PgChartCandidateLift>>
    refine_lift(PgChartCandidateInput independent, PgChartCandidateInput original_rhs,
                const PgChartCandidateSolve&, const PgChartRefinement&, core::Workspace&,
                core::ResourcePlan&) const noexcept;

    [[nodiscard]] core::Result<PgChartRefined<PgChartCandidateDual>>
    refine_transpose(PgChartCandidateInput force, const PgChartCandidateSolve&,
                     const PgChartRefinement&, core::Workspace&,
                     core::ResourcePlan&) const noexcept;

    /// Check actual retained work: f^T*q = h^T*r + lambda^T*c.
    [[nodiscard]] core::Result<double>
    work_defect(PgChartCandidateInput force, PgChartCandidateInput independent,
                PgChartCandidateInput original_rhs, const PgChartCandidateLift&,
                const PgChartCandidateDual&, double maximum_absolute_error,
                std::uint64_t maximum_visits, std::stop_token) const noexcept;

private:
    friend class PgChartOwner;
    struct Data;
    PgChartCandidate(PgSparseChart&&, core::PlannedBuffer<Data>) noexcept;
    PgSparseChart exact_;
    core::PlannedBuffer<Data> data_;
};

/// Owns one exact chart throughout optional numerical preparation.
/// A successful preparation invalidates prior borrows of this chart and candidate.
class PgChartOwner {
public:
    explicit PgChartOwner(PgSparseChart&& exact) noexcept : storage_(std::move(exact)) {}
    PgChartOwner(PgChartOwner&&) noexcept = default;
    PgChartOwner& operator=(PgChartOwner&&) noexcept = default;
    PgChartOwner(const PgChartOwner&) = delete;
    PgChartOwner& operator=(const PgChartOwner&) = delete;

    [[nodiscard]] const PgSparseChart& exact() const& noexcept;
    [[nodiscard]] const PgSparseChart& exact() const&& = delete;
    [[nodiscard]] const PgChartCandidate* numerical() const& noexcept;
    [[nodiscard]] const PgChartCandidate* numerical() const&& = delete;

    /// Prepare from the current exact chart. Failure retains both exact and numerical owners.
    [[nodiscard]] core::Result<void> prepare_numerical(const numerics::linear::Scales&,
                                                       const numerics::linear::Scales&,
                                                       const PgChartCandidatePreparation&,
                                                       core::ResourcePlan&) & noexcept;

private:
    std::variant<PgSparseChart, PgChartCandidate> storage_;
};

} // namespace uvlas::structure
