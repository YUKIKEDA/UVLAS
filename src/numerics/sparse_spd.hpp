#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <uvlas/numerics/precision.hpp>

namespace uvlas::numerics::linear::detail {

// Nonnegative outward arithmetic. Exact zero terms stay zero; a positive
// underflow is enclosed by the next representable value rather than discarded.
[[nodiscard]] inline double sparse_up_add(double a, double b) noexcept {
    if (a == 0) {
        return b;
    }
    if (b == 0) {
        return a;
    }
    return std::nextafter(a + b, std::numeric_limits<double>::infinity());
}
[[nodiscard]] inline double sparse_up_multiply(double a, double b) noexcept {
    if (a == 0 || b == 0) {
        return 0;
    }
    return std::nextafter(a * b, std::numeric_limits<double>::infinity());
}

// Rump (2006), Corollary 2.4 and section 3 bound II. The caller verifies exact
// symmetry/congruence first. This computes only an error upper bound, never an
// SPD certificate: a separate downward-shifted Cholesky check is still required.
[[nodiscard]] inline core::Result<double> rump_sparse_shift(const SparseReal& matrix) noexcept {
    constexpr auto operation = core::Operation::SolveLinear;
    const auto reject = [](core::FailureCode code) -> core::Result<double> {
        return std::unexpected(core::Failure{code, operation});
    };
    if (matrix.rows() != matrix.columns()) {
        return reject(core::FailureCode::InvalidInput);
    }
    const auto order = matrix.rows();
    if (order == 0) {
        return 0.0;
    }
    // All subsequent dimension conversions are exact, with (t+2)*u < 1/2.
    constexpr std::int64_t maximum_order = (std::int64_t{1} << 52) - 2;
    if (order > maximum_order) {
        return reject(core::FailureCode::PrecisionUnresolved);
    }
    const auto environment = two_sum(0, 0);
    if (!environment) {
        return std::unexpected(environment.error());
    }
    constexpr double unit = 0x1p-53;
    const auto offsets = matrix.column_offsets();
    const auto rows = matrix.row_indices();
    const auto values = matrix.values();
    double shift = 0;
    double maximum_diagonal = 0;
    for (core::Index column = 0; column < order; ++column) {
        const auto start = offsets[static_cast<std::size_t>(column)];
        const auto finish = offsets[static_cast<std::size_t>(column + 1)];
        if (start == finish) {
            return reject(core::FailureCode::RankUnresolved);
        }
        double diagonal = 0;
        for (auto entry = start; entry < finish; ++entry) {
            const auto index = static_cast<std::size_t>(entry);
            if (!std::isfinite(values[index])) {
                return reject(core::FailureCode::Nonfinite);
            }
            if (rows[index] == column) {
                diagonal = values[index];
            }
        }
        if (diagonal <= 0) {
            return reject(core::FailureCode::RankUnresolved);
        }
        maximum_diagonal = std::max(maximum_diagonal, diagonal);
        // Envelope width in the PERMUTED matrix, not the number of stored entries.
        // Including an explicitly stored zero can only enlarge this bound.
        const auto first = rows[static_cast<std::size_t>(start)];
        const auto width = column - first;
        const auto count = width + 2;
        const double ku = static_cast<double>(count) * unit;
        if (!(ku < 0.5)) {
            return reject(core::FailureCode::PrecisionUnresolved);
        }
        const double denominator = std::nextafter(1 - 2 * ku, 0.0);
        if (!(denominator > 0)) {
            return reject(core::FailureCode::PrecisionUnresolved);
        }
        const double ratio =
            std::nextafter(ku / denominator, std::numeric_limits<double>::infinity());
        shift = sparse_up_add(shift, sparse_up_multiply(ratio, diagonal));
    }
    const double n = static_cast<double>(order);
    const double n_eta = sparse_up_multiply(n, std::numeric_limits<double>::denorm_min());
    // Expand 3*n*(2*n+maxdiag)*eta before evaluating: the original parenthesized
    // sum can overflow although this underflow correction has a finite upper bound.
    const double underflow =
        sparse_up_add(sparse_up_multiply(sparse_up_multiply(6, n), n_eta),
                      sparse_up_multiply(sparse_up_multiply(3, n_eta), maximum_diagonal));
    shift = sparse_up_add(shift, underflow);
    if (!std::isfinite(shift)) {
        return reject(core::FailureCode::PrecisionUnresolved);
    }
    return shift;
}

// Numeric symmetry is checked on the complete original CSC, including absent
// zero entries. The caller never averages triangles or discards a mismatch.
[[nodiscard]] inline core::Result<bool> sparse_exact_symmetry(const SparseReal& matrix) noexcept {
    if (matrix.rows() != matrix.columns()) {
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::SolveLinear});
    }
    const auto offsets = matrix.column_offsets();
    const auto rows = matrix.row_indices();
    const auto values = matrix.values();
    for (core::Index column = 0; column < matrix.columns(); ++column) {
        for (auto entry = offsets[static_cast<std::size_t>(column)];
             entry < offsets[static_cast<std::size_t>(column + 1)]; ++entry) {
            const auto index = static_cast<std::size_t>(entry);
            const auto row = rows[index];
            const auto begin = rows.begin() + offsets[static_cast<std::size_t>(row)];
            const auto end = rows.begin() + offsets[static_cast<std::size_t>(row + 1)];
            const auto opposite = std::lower_bound(begin, end, column);
            const double other = opposite != end && *opposite == column
                                     ? values[static_cast<std::size_t>(opposite - rows.begin())]
                                     : 0;
            if (values[index] != other) {
                return false;
            }
        }
    }
    return true;
}

struct SparseCongruence {
    SparseReal matrix;
    bool scaled;
};

// Prepare a separate validation matrix. Every power-of-two transformation must
// reverse exactly; otherwise preserve the complete original coefficients instead.
[[nodiscard]] inline core::Result<SparseCongruence>
sparse_exact_congruence(const SparseReal& matrix, core::ResourcePlan& plan) noexcept {
    constexpr auto operation = core::Operation::SolveLinear;
    if (matrix.rows() != matrix.columns()) {
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    }
    auto exponents = plan.allocate<int>(matrix.rows(), operation);
    if (!exponents) {
        return std::unexpected(exponents.error());
    }
    const auto offsets = matrix.column_offsets();
    const auto rows = matrix.row_indices();
    const auto values = matrix.values();
    for (core::Index column = 0; column < matrix.columns(); ++column) {
        double diagonal = 0;
        for (auto entry = offsets[static_cast<std::size_t>(column)];
             entry < offsets[static_cast<std::size_t>(column + 1)]; ++entry) {
            const auto index = static_cast<std::size_t>(entry);
            if (rows[index] == column) {
                diagonal = values[index];
            }
        }
        if (!(diagonal > 0)) {
            return std::unexpected(core::Failure{core::FailureCode::RankUnresolved, operation});
        }
        int exponent = 0;
        std::frexp(diagonal, &exponent);
        int half = exponent / 2;
        if (exponent < 0 && exponent % 2 != 0) {
            --half;
        }
        exponents->values()[static_cast<std::size_t>(column)] = -half;
    }
    bool exact = true;
    for (core::Index column = 0; column < matrix.columns() && exact; ++column) {
        for (auto entry = offsets[static_cast<std::size_t>(column)];
             entry < offsets[static_cast<std::size_t>(column + 1)]; ++entry) {
            const auto index = static_cast<std::size_t>(entry);
            const auto power = exponents->values()[static_cast<std::size_t>(column)] +
                               exponents->values()[static_cast<std::size_t>(rows[index])];
            const double scaled = std::scalbn(values[index], power);
            if (!std::isfinite(scaled) || std::scalbn(scaled, -power) != values[index]) {
                exact = false;
                break;
            }
        }
    }
    if (!exact) {
        auto copy = SparseReal::copy_csc(matrix.rows(), matrix.columns(), offsets, rows, values,
                                         plan, operation);
        if (!copy) {
            return std::unexpected(copy.error());
        }
        return SparseCongruence{std::move(*copy), false};
    }
    auto scaled = plan.allocate<double>(static_cast<core::Index>(values.size()), operation);
    if (!scaled) {
        return std::unexpected(scaled.error());
    }
    for (core::Index column = 0; column < matrix.columns(); ++column) {
        for (auto entry = offsets[static_cast<std::size_t>(column)];
             entry < offsets[static_cast<std::size_t>(column + 1)]; ++entry) {
            const auto index = static_cast<std::size_t>(entry);
            const auto power = exponents->values()[static_cast<std::size_t>(column)] +
                               exponents->values()[static_cast<std::size_t>(rows[index])];
            scaled->values()[index] = std::scalbn(values[index], power);
        }
    }
    auto copy = SparseReal::copy_csc(matrix.rows(), matrix.columns(), offsets, rows,
                                     scaled->values(), plan, operation);
    if (!copy) {
        return std::unexpected(copy.error());
    }
    return SparseCongruence{std::move(*copy), true};
}

struct SparseSpdQualification {
    bool certified;
    std::optional<double> shift;
};

// The public preparation entry must admit the backend memory plan BEFORE this
// helper constructs Eigen ordering/factor storage. This returns only qualification
// diagnostics; it never solves a shifted system or changes the original matrix.
[[nodiscard]] core::Result<SparseSpdQualification>
qualify_sparse_positive_definite(const SparseReal& matrix, core::ResourcePlan& plan) noexcept;

} // namespace uvlas::numerics::linear::detail
