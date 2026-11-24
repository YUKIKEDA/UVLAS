#pragma once

#include <uvlas/numerics/sparse_solve.hpp>
#include <uvlas/structure/pg_sparse_chart.hpp>

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
         const PgChartCandidateSolve&, core::Workspace&, core::ResourcePlan&) noexcept;

    [[nodiscard]] core::Result<PgChartCandidateDual> transpose(PgChartCandidateInput force,
                                                               const PgChartCandidateSolve&,
                                                               core::Workspace&,
                                                               core::ResourcePlan&) noexcept;

    /// Check actual retained work: f^T*q = h^T*r + lambda^T*c.
    [[nodiscard]] core::Result<double>
    work_defect(PgChartCandidateInput force, PgChartCandidateInput independent,
                PgChartCandidateInput original_rhs, const PgChartCandidateLift&,
                const PgChartCandidateDual&, double maximum_absolute_error,
                std::uint64_t maximum_visits, std::stop_token) const noexcept;

private:
    struct Data;
    PgChartCandidate(PgSparseChart&&, core::PlannedBuffer<Data>) noexcept;
    PgSparseChart exact_;
    core::PlannedBuffer<Data> data_;
};

} // namespace uvlas::structure
