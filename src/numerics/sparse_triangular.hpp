#pragma once

#include "eigen_views.hpp"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

#include <span>
#include <uvlas/core/resource_plan.hpp>

namespace uvlas::numerics::linear::detail {
using SparseLuFactor =
    Eigen::SparseLU<uvlas::numerics::detail::EigenSparse, Eigen::COLAMDOrdering<core::Index>>;

using SparseLltFactor = Eigen::SimplicialLLT<uvlas::numerics::detail::EigenSparse, Eigen::Lower,
                                             Eigen::AMDOrdering<core::Index>>;

// Immutable after preparation. L owns an explicit unit diagonal; U owns every
// nonzero pivot. Both original-index maps point INTO factor coordinate order.
struct SparseTriangularFactors {
    SparseReal lower;
    SparseReal upper;
    core::PlannedBuffer<core::Index> row_old_to_factor;
    core::PlannedBuffer<core::Index> column_old_to_factor;
    bool unit_lower = true;
};

// Extract before repeated solves. Eigen stores within-supernode U in mapL and
// off-supernode U in mapU; only strict-lower mapL entries belong to unit L.
[[nodiscard]] core::Result<SparseTriangularFactors>
extract_sparse_lu(const SparseLuFactor& factor, core::ResourcePlan& plan) noexcept;

// Extract a qualified original LLT with explicit non-unit L and its exact transpose U.
[[nodiscard]] core::Result<SparseTriangularFactors>
extract_sparse_llt(const SparseLltFactor& factor, core::ResourcePlan& plan) noexcept;

// For R=Eigen.rowsPermutation(), E=Eigen.colsPermutation(), R*A*E^-1=L*U.
// work[row_old_to_factor[i]]=b[i]; x[i]=work[column_old_to_factor[i]].
// All logical spans are packed column-major and mutually non-overlapping.
// Overlap returns InvalidInput before updates; no scratch/output allocation occurs.
[[nodiscard]] core::Result<void> solve_sparse_triangular(const SparseTriangularFactors& factors,
                                                         std::span<const double> rhs,
                                                         core::Index columns,
                                                         std::span<double> work,
                                                         std::span<double> output) noexcept;
} // namespace uvlas::numerics::linear::detail
