#include "sparse_spd.hpp"

#include "eigen_views.hpp"

#include <Eigen/OrderingMethods>
#include <Eigen/SparseCholesky>

#include <cmath>
#include <limits>

namespace uvlas::numerics::linear::detail {

core::Result<SparseSpdQualification>
qualify_sparse_positive_definite(const SparseReal& matrix, core::ResourcePlan& plan) noexcept {
    constexpr auto operation = core::Operation::SolveLinear;
    return core::guard_entry(
        operation, core::FailureCode::FactorizationFailed, {},
        [&]() -> core::Result<SparseSpdQualification> {
            const auto symmetric = sparse_exact_symmetry(matrix);
            if (!symmetric) {
                return std::unexpected(symmetric.error());
            }
            if (!*symmetric) {
                return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
            }
            const auto n = matrix.rows();
            if (n == 0) {
                return SparseSpdQualification{true, 0.0};
            }
            if (n > (core::Index{1} << 52) - 2) {
                return SparseSpdQualification{false, std::nullopt};
            }
            auto congruence = sparse_exact_congruence(matrix, plan);
            if (!congruence) {
                return std::unexpected(congruence.error());
            }
            const auto input =
                uvlas::numerics::detail::eigen_read_map(congruence->matrix, n, n, operation);
            if (!input) {
                return std::unexpected(input.error());
            }
            using Matrix = uvlas::numerics::detail::EigenSparse;
            Eigen::AMDOrdering<core::Index> ordering;
            Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, core::Index> permutation;
            ordering(input->template selfadjointView<Eigen::Lower>(), permutation);
            Matrix reordered = permutation * *input * permutation.transpose();
            reordered.makeCompressed();
            auto permuted = SparseReal::copy_csc(
                n, n, {reordered.outerIndexPtr(), static_cast<std::size_t>(n + 1)},
                {reordered.innerIndexPtr(), static_cast<std::size_t>(reordered.nonZeros())},
                {reordered.valuePtr(), static_cast<std::size_t>(reordered.nonZeros())}, plan,
                operation);
            if (!permuted) {
                return std::unexpected(permuted.error());
            }
            const auto bound = rump_sparse_shift(*permuted);
            if (!bound) {
                return std::unexpected(bound.error());
            }
            auto shifted = plan.allocate<double>(
                static_cast<core::Index>(permuted->values().size()), operation);
            if (!shifted) {
                return std::unexpected(shifted.error());
            }
            std::copy(permuted->values().begin(), permuted->values().end(),
                      shifted->values().begin());
            for (core::Index column = 0; column < n; ++column) {
                for (auto entry = permuted->column_offsets()[static_cast<std::size_t>(column)];
                     entry < permuted->column_offsets()[static_cast<std::size_t>(column + 1)];
                     ++entry) {
                    const auto index = static_cast<std::size_t>(entry);
                    if (permuted->row_indices()[index] == column) {
                        const double diagonal =
                            std::nextafter(shifted->values()[index] - *bound,
                                           -std::numeric_limits<double>::infinity());
                        if (!std::isfinite(diagonal) || diagonal <= 0) {
                            return SparseSpdQualification{false, *bound};
                        }
                        shifted->values()[index] = diagonal;
                    }
                }
            }
            auto candidate =
                SparseReal::copy_csc(n, n, permuted->column_offsets(), permuted->row_indices(),
                                     shifted->values(), plan, operation);
            if (!candidate) {
                return std::unexpected(candidate.error());
            }
            const auto mapped =
                uvlas::numerics::detail::eigen_read_map(*candidate, n, n, operation);
            if (!mapped) {
                return std::unexpected(mapped.error());
            }
            // Natural ordering preserves the exact envelope used for the bound.
            using Factor =
                Eigen::SimplicialLLT<Matrix, Eigen::Lower, Eigen::NaturalOrdering<core::Index>>;
            Factor factor;
            factor.compute(*mapped);
            if (factor.info() != Eigen::Success) {
                return SparseSpdQualification{false, *bound};
            }
            const Factor::MatrixL view = factor.matrixL();
            const Matrix& data = view.nestedExpression();
            for (core::Index column = 0; column < n; ++column) {
                bool diagonal = false;
                for (Matrix::InnerIterator entry(data, column); entry; ++entry) {
                    if (!std::isfinite(entry.value())) {
                        return SparseSpdQualification{false, *bound};
                    }
                    if (entry.row() == column) {
                        if (!(entry.value() > 0)) {
                            return SparseSpdQualification{false, *bound};
                        }
                        diagonal = true;
                    }
                }
                if (!diagonal) {
                    return SparseSpdQualification{false, *bound};
                }
            }
            return SparseSpdQualification{true, *bound};
        });
}

} // namespace uvlas::numerics::linear::detail
