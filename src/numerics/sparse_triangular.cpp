#include "sparse_triangular.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <type_traits>

namespace uvlas::numerics::linear::detail {
namespace {
constexpr auto operation = core::Operation::SolveLinear;
auto fail(core::FailureCode code) noexcept {
    return std::unexpected(core::Failure{code, operation});
}
core::Result<void> solve_transposed_columns(const SparseReal& factor, std::span<double> work,
                                            bool upper, bool unit_diagonal) noexcept {
    const auto row_indices = factor.row_indices();
    const auto values = factor.values();
    const auto offsets = factor.column_offsets();
    const auto n = work.size();
    for (std::size_t step = 0; step < n; ++step) {
        const auto column = upper ? step : n - 1 - step;
        const auto begin = static_cast<std::size_t>(offsets[column]);
        const auto end = static_cast<std::size_t>(offsets[column + 1]);
        if (begin == end) {
            return fail(core::FailureCode::FactorizationFailed);
        }
        const auto diagonal = upper ? end - 1 : begin;
        if (row_indices[diagonal] != static_cast<core::Index>(column)) {
            return fail(core::FailureCode::FactorizationFailed);
        }
        const double pivot = values[diagonal];
        if (!unit_diagonal && (pivot == 0 || !std::isfinite(pivot))) {
            return fail(core::FailureCode::RankUnresolved);
        }
        double value = work[column];
        for (auto slot = upper ? begin : begin + 1; slot < (upper ? diagonal : end); ++slot) {
            value -= values[slot] * work[static_cast<std::size_t>(row_indices[slot])];
            if (!std::isfinite(value)) {
                return fail(core::FailureCode::Nonfinite);
            }
        }
        if (!unit_diagonal) {
            value /= pivot;
        }
        if (!std::isfinite(value)) {
            return fail(core::FailureCode::Nonfinite);
        }
        work[column] = value;
    }
    return {};
}
} // namespace

core::Result<SparseTriangularFactors> extract_sparse_lu(const SparseLuFactor& factor,
                                                        core::ResourcePlan& plan) noexcept {
    return core::guard_entry(
        operation, core::FailureCode::FactorizationFailed, {},
        [&]() -> core::Result<SparseTriangularFactors> {
            if (factor.info() != Eigen::Success) {
                return fail(core::FailureCode::RankUnresolved);
            }
            const auto n = static_cast<core::Index>(factor.rows());
            using UpperView = decltype(factor.matrixU());
            const UpperView upper = factor.matrixU();
            using Super = std::remove_cvref_t<decltype(upper.m_mapL)>;
            using OffUpper = std::remove_cvref_t<decltype(upper.m_mapU)>;
            core::Index lower_count = n;
            core::Index upper_count = 0;
            const auto increment = [](core::Index& value) noexcept {
                if (value == std::numeric_limits<core::Index>::max()) {
                    return false;
                }
                ++value;
                return true;
            };
            for (core::Index column = 0; column < n; ++column) {
                bool diagonal = false;
                for (Super::InnerIterator entry(upper.m_mapL, column); entry; ++entry) {
                    if (!std::isfinite(entry.value())) {
                        return fail(core::FailureCode::Nonfinite);
                    }
                    if (entry.row() < 0 || entry.row() >= n) {
                        return fail(core::FailureCode::FactorizationFailed);
                    }
                    auto& count = entry.row() > column ? lower_count : upper_count;
                    if (!increment(count)) {
                        return fail(core::FailureCode::ResourceExceeded);
                    }
                    if (entry.row() == column) {
                        if (diagonal || entry.value() == 0) {
                            return fail(core::FailureCode::RankUnresolved);
                        }
                        diagonal = true;
                    }
                }
                if (!diagonal) {
                    return fail(core::FailureCode::RankUnresolved);
                }
                for (OffUpper::InnerIterator entry(upper.m_mapU, column); entry; ++entry) {
                    if (!std::isfinite(entry.value())) {
                        return fail(core::FailureCode::Nonfinite);
                    }
                    if (entry.row() < 0 || entry.row() >= column) {
                        return fail(core::FailureCode::FactorizationFailed);
                    }
                    if (!increment(upper_count)) {
                        return fail(core::FailureCode::ResourceExceeded);
                    }
                }
            }
            auto lower = plan.allocate<Triplet>(lower_count, operation);
            if (!lower) {
                return std::unexpected(lower.error());
            }
            auto u = plan.allocate<Triplet>(upper_count, operation);
            if (!u) {
                return std::unexpected(u.error());
            }
            auto row = plan.allocate<core::Index>(n, operation);
            if (!row) {
                return std::unexpected(row.error());
            }
            auto column_map = plan.allocate<core::Index>(n, operation);
            if (!column_map) {
                return std::unexpected(column_map.error());
            }
            auto seen = plan.allocate<unsigned char>(n, operation);
            if (!seen) {
                return std::unexpected(seen.error());
            }
            if (factor.rowsPermutation().size() != n || factor.colsPermutation().size() != n) {
                return fail(core::FailureCode::FactorizationFailed);
            }
            for (unsigned pass = 0; pass < 2; ++pass) {
                std::fill(seen->values().begin(), seen->values().end(), 0);
                const auto& permutation =
                    pass == 0 ? factor.rowsPermutation() : factor.colsPermutation();
                auto destination = pass == 0 ? row->values() : column_map->values();
                for (core::Index original = 0; original < n; ++original) {
                    const auto mapped = permutation.indices()[original];
                    if (mapped < 0 || mapped >= n ||
                        seen->values()[static_cast<std::size_t>(mapped)] != 0) {
                        return fail(core::FailureCode::FactorizationFailed);
                    }
                    seen->values()[static_cast<std::size_t>(mapped)] = 1;
                    destination[static_cast<std::size_t>(original)] = mapped;
                }
            }
            std::size_t l_position = 0;
            std::size_t u_position = 0;
            for (core::Index column = 0; column < n; ++column) {
                lower->values()[l_position++] = Triplet{column, column, 1};
                for (Super::InnerIterator entry(upper.m_mapL, column); entry; ++entry) {
                    if (entry.row() > column) {
                        lower->values()[l_position++] = Triplet{entry.row(), column, entry.value()};
                    } else {
                        u->values()[u_position++] = Triplet{entry.row(), column, entry.value()};
                    }
                }
                for (OffUpper::InnerIterator entry(upper.m_mapU, column); entry; ++entry) {
                    u->values()[u_position++] = Triplet{entry.row(), column, entry.value()};
                }
            }
            auto lower_matrix = SparseReal::from_triplets(n, n, lower->values(), plan, operation);
            if (!lower_matrix) {
                return std::unexpected(lower_matrix.error());
            }
            auto upper_matrix = SparseReal::from_triplets(n, n, u->values(), plan, operation);
            if (!upper_matrix) {
                return std::unexpected(upper_matrix.error());
            }
            return SparseTriangularFactors{std::move(*lower_matrix), std::move(*upper_matrix),
                                           std::move(*row), std::move(*column_map), true};
        });
}

core::Result<SparseTriangularFactors> extract_sparse_llt(const SparseLltFactor& factor,
                                                         core::ResourcePlan& plan) noexcept {
    return core::guard_entry(
        operation, core::FailureCode::FactorizationFailed, {},
        [&]() -> core::Result<SparseTriangularFactors> {
            if (factor.info() != Eigen::Success) {
                return fail(core::FailureCode::FactorizationFailed);
            }
            using Matrix = uvlas::numerics::detail::EigenSparse;
            const SparseLltFactor::MatrixL view = factor.matrixL();
            const Matrix& data = view.nestedExpression();
            const auto n = static_cast<core::Index>(data.rows());
            const auto count = static_cast<core::Index>(data.nonZeros());
            auto lower = plan.allocate<Triplet>(count, operation);
            if (!lower) {
                return std::unexpected(lower.error());
            }
            auto upper = plan.allocate<Triplet>(count, operation);
            if (!upper) {
                return std::unexpected(upper.error());
            }
            auto row = plan.allocate<core::Index>(n, operation);
            if (!row) {
                return std::unexpected(row.error());
            }
            auto column_map = plan.allocate<core::Index>(n, operation);
            if (!column_map) {
                return std::unexpected(column_map.error());
            }
            auto seen = plan.allocate<unsigned char>(n, operation);
            if (!seen) {
                return std::unexpected(seen.error());
            }
            const auto& permutation = factor.permutationP();
            if (permutation.size() != 0 && permutation.size() != n) {
                return fail(core::FailureCode::FactorizationFailed);
            }
            for (core::Index original = 0; original < n; ++original) {
                const auto mapped =
                    permutation.size() == 0 ? original : permutation.indices()[original];
                if (mapped < 0 || mapped >= n ||
                    seen->values()[static_cast<std::size_t>(mapped)] != 0) {
                    return fail(core::FailureCode::FactorizationFailed);
                }
                seen->values()[static_cast<std::size_t>(mapped)] = 1;
                row->values()[static_cast<std::size_t>(original)] = mapped;
                column_map->values()[static_cast<std::size_t>(original)] = mapped;
            }
            std::size_t position = 0;
            for (core::Index column = 0; column < n; ++column) {
                bool diagonal = false;
                for (Matrix::InnerIterator entry(data, column); entry; ++entry) {
                    if (!std::isfinite(entry.value())) {
                        return fail(core::FailureCode::Nonfinite);
                    }
                    if (entry.row() < column || entry.row() >= n) {
                        return fail(core::FailureCode::FactorizationFailed);
                    }
                    if (entry.row() == column) {
                        if (diagonal || entry.value() <= 0) {
                            return fail(core::FailureCode::RankUnresolved);
                        }
                        diagonal = true;
                    }
                    lower->values()[position] = Triplet{entry.row(), column, entry.value()};
                    upper->values()[position] = Triplet{column, entry.row(), entry.value()};
                    ++position;
                }
                if (!diagonal) {
                    return fail(core::FailureCode::FactorizationFailed);
                }
            }
            auto l = SparseReal::from_triplets(n, n, lower->values(), plan, operation);
            if (!l) {
                return std::unexpected(l.error());
            }
            auto u = SparseReal::from_triplets(n, n, upper->values(), plan, operation);
            if (!u) {
                return std::unexpected(u.error());
            }
            return SparseTriangularFactors{std::move(*l), std::move(*u), std::move(*row),
                                           std::move(*column_map), false};
        });
}

core::Result<void> solve_sparse_triangular(const SparseTriangularFactors& factors,
                                           std::span<const double> rhs, core::Index columns,
                                           std::span<double> work, std::span<double> output,
                                           LinearSystemOrientation orientation) noexcept {
    const auto n = factors.lower.rows();
    const auto count = core::checked_extent_product(std::array{n, columns}, operation);
    if (!count) {
        return std::unexpected(count.error());
    }
    if (factors.lower.columns() != n || factors.upper.rows() != n || factors.upper.columns() != n ||
        factors.row_old_to_factor.size() != static_cast<std::size_t>(n) ||
        factors.column_old_to_factor.size() != static_cast<std::size_t>(n) ||
        rhs.size() != static_cast<std::size_t>(*count) || output.size() != rhs.size() ||
        work.size() < rhs.size()) {
        return fail(core::FailureCode::InvalidInput);
    }
    const auto overlaps = [](std::span<const double> a, std::span<const double> b) noexcept {
        if (a.empty() || b.empty()) {
            return false;
        }
        const std::less<const double*> before;
        return before(a.data(), b.data() + b.size()) && before(b.data(), a.data() + a.size());
    };
    const auto logical_work = work.first(rhs.size());
    if (overlaps(rhs, logical_work) || overlaps(output, logical_work) || overlaps(rhs, output)) {
        return fail(core::FailureCode::InvalidInput);
    }
    if (n == 0 || columns == 0) {
        return {};
    }
    const bool transposed = orientation == LinearSystemOrientation::Transpose;
    const auto rows = factors.row_old_to_factor.values();
    const auto columns_map = factors.column_old_to_factor.values();
    const auto rhs_map = transposed ? columns_map : rows;
    const auto output_map = transposed ? rows : columns_map;
    const auto l_offsets = factors.lower.column_offsets();
    const auto l_rows = factors.lower.row_indices();
    const auto l_values = factors.lower.values();
    const auto u_offsets = factors.upper.column_offsets();
    const auto u_rows = factors.upper.row_indices();
    const auto u_values = factors.upper.values();
    for (core::Index rhs_column = 0; rhs_column < columns; ++rhs_column) {
        const auto start = static_cast<std::size_t>(rhs_column * n);
        for (core::Index original = 0; original < n; ++original) {
            const auto index = start + static_cast<std::size_t>(original);
            if (!std::isfinite(rhs[index])) {
                return fail(core::FailureCode::Nonfinite);
            }
            work[start + static_cast<std::size_t>(rhs_map[static_cast<std::size_t>(original)])] =
                rhs[index];
        }
        if (transposed) {
            const auto column_work = work.subspan(start, static_cast<std::size_t>(n));
            const auto upper = solve_transposed_columns(factors.upper, column_work, true, false);
            if (!upper) {
                return upper;
            }
            const auto lower =
                solve_transposed_columns(factors.lower, column_work, false, factors.unit_lower);
            if (!lower) {
                return lower;
            }
        } else {
            for (core::Index column = 0; column < n; ++column) {
                const auto begin = l_offsets[static_cast<std::size_t>(column)];
                const auto end = l_offsets[static_cast<std::size_t>(column + 1)];
                if (begin == end || l_rows[static_cast<std::size_t>(begin)] != column) {
                    return fail(core::FailureCode::FactorizationFailed);
                }
                auto& value = work[start + static_cast<std::size_t>(column)];
                if (!factors.unit_lower) {
                    value /= l_values[static_cast<std::size_t>(begin)];
                }
                if (!std::isfinite(value)) {
                    return fail(core::FailureCode::Nonfinite);
                }
                for (auto entry = begin + 1; entry < end; ++entry) {
                    const auto index = static_cast<std::size_t>(entry);
                    auto& target = work[start + static_cast<std::size_t>(l_rows[index])];
                    target -= l_values[index] * value;
                    if (!std::isfinite(target)) {
                        return fail(core::FailureCode::Nonfinite);
                    }
                }
            }
            for (auto column = n; column > 0;) {
                --column;
                const auto begin = u_offsets[static_cast<std::size_t>(column)];
                const auto end = u_offsets[static_cast<std::size_t>(column + 1)];
                if (begin == end || u_rows[static_cast<std::size_t>(end - 1)] != column) {
                    return fail(core::FailureCode::FactorizationFailed);
                }
                const double pivot = u_values[static_cast<std::size_t>(end - 1)];
                if (pivot == 0 || !std::isfinite(pivot)) {
                    return fail(core::FailureCode::RankUnresolved);
                }
                auto& value = work[start + static_cast<std::size_t>(column)];
                value /= pivot;
                if (!std::isfinite(value)) {
                    return fail(core::FailureCode::Nonfinite);
                }
                for (auto entry = begin; entry < end - 1; ++entry) {
                    const auto index = static_cast<std::size_t>(entry);
                    auto& target = work[start + static_cast<std::size_t>(u_rows[index])];
                    target -= u_values[index] * value;
                    if (!std::isfinite(target)) {
                        return fail(core::FailureCode::Nonfinite);
                    }
                }
            }
        }
        for (core::Index original = 0; original < n; ++original) {
            const auto index = static_cast<std::size_t>(original);
            output[start + index] = work[start + static_cast<std::size_t>(output_map[index])];
        }
    }
    return {};
}
} // namespace uvlas::numerics::linear::detail
