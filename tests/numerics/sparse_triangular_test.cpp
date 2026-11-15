#include "sparse_triangular.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <type_traits>
#include <utility>
#include <uvlas/numerics/sparse_solve.hpp>
#include <uvlas/numerics/storage.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;
constexpr auto operation = core::Operation::SolveLinear;

template <class T>
concept BorrowsSparseWitness = requires(T&& value) { std::forward<T>(value).witness(); };
static_assert(!std::is_constructible_v<SparseOperator, SparseReal&&, Context>);
static_assert(!std::is_copy_constructible_v<SparsePrepared>);
static_assert(BorrowsSparseWitness<SparseSpdDiagnostic&>);
static_assert(!BorrowsSparseWitness<SparseSpdDiagnostic>);
static_assert(!BorrowsSparseWitness<const SparseSpdDiagnostic>);

void check_solve(core::Index n, std::span<const Triplet> triplets, std::span<const double> rhs,
                 std::span<const double> expected, bool require_both_permutations) {
    core::ResourcePlan plan{1000000};
    auto matrix = SparseReal::from_triplets(n, n, triplets, plan, operation);
    if (!matrix) {
        FAIL("Reference CSC");
        return;
    }
    const auto mapped = uvlas::numerics::detail::eigen_read_map(*matrix, n, n, operation);
    if (!mapped) {
        FAIL("Reference map");
        return;
    }
    linear::detail::SparseLuFactor factor;
    factor.compute(*mapped);
    REQUIRE(factor.info() == Eigen::Success);
    auto extracted = linear::detail::extract_sparse_lu(factor, plan);
    if (!extracted) {
        FAIL("Factor extraction");
        return;
    }
    bool row_changed = false, column_changed = false;
    for (core::Index i = 0; i < n; ++i) {
        row_changed |= extracted->row_old_to_factor.values()[static_cast<std::size_t>(i)] != i;
        column_changed |=
            extracted->column_old_to_factor.values()[static_cast<std::size_t>(i)] != i;
    }
    if (require_both_permutations) {
        CHECK(row_changed);
        CHECK(column_changed);
    }
    auto work = plan.allocate<double>(static_cast<core::Index>(rhs.size()), operation);
    auto result = plan.allocate<double>(static_cast<core::Index>(rhs.size()), operation);
    if (!work || !result) {
        FAIL("Reference workspace");
        return;
    }
    const auto columns = static_cast<core::Index>(rhs.size()) / n;
    std::copy(rhs.begin(), rhs.end(), work->values().begin());
    auto work_input_alias = linear::detail::solve_sparse_triangular(
        *extracted, work->values(), columns, work->values(), result->values());
    REQUIRE_FALSE(work_input_alias);
    CHECK(work_input_alias.error().code == core::FailureCode::InvalidInput);
    CHECK(std::equal(work->values().begin(), work->values().end(), rhs.begin(), rhs.end()));
    auto work_output_alias = linear::detail::solve_sparse_triangular(
        *extracted, rhs, columns, work->values(), work->values());
    REQUIRE_FALSE(work_output_alias);
    CHECK(work_output_alias.error().code == core::FailureCode::InvalidInput);
    CHECK(std::equal(work->values().begin(), work->values().end(), rhs.begin(), rhs.end()));
    std::copy(rhs.begin(), rhs.end(), result->values().begin());
    auto input_output_alias = linear::detail::solve_sparse_triangular(
        *extracted, result->values(), columns, work->values(), result->values());
    REQUIRE_FALSE(input_output_alias);
    CHECK(input_output_alias.error().code == core::FailureCode::InvalidInput);
    CHECK(std::equal(result->values().begin(), result->values().end(), rhs.begin(), rhs.end()));
    auto solved = linear::detail::solve_sparse_triangular(*extracted, rhs, columns, work->values(),
                                                          result->values());
    REQUIRE(solved.has_value());
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        CHECK(std::abs(result->values()[i] - expected[i]) <= 1e-12);
    }
    // The stock solve is an extra dependency-specific contrast; literal RHS/solution
    // above independently fix the row/column permutation directions.
    using Dense = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    const Eigen::Map<const Dense> input{rhs.data(), n, columns};
    const Dense stock = factor.solve(input);
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        CHECK(std::abs(result->values()[i] - stock.data()[i]) <= 1e-12);
    }
}
} // namespace

TEST_CASE("Extracted sparse LU solves a literal indefinite KKT with multiple RHS", "[N4-05]") {
    const std::array<Triplet, 6> entries{
        {{0, 0, 2}, {2, 0, 1}, {1, 1, 3}, {2, 1, 1}, {0, 2, 1}, {1, 2, 1}}};
    const std::array rhs{4.0, -1.0, 0.0, 8.0, -2.0, 0.0};
    const std::array expected{1.0, -1.0, 2.0, 2.0, -2.0, 4.0};
    check_solve(3, entries, rhs, expected, false);
}

TEST_CASE("Sparse row and column permutations are both nontrivial in a literal system", "[N4-05]") {
    const std::array<Triplet, 11> entries{{{3, 0, 1},
                                           {1, 0, 3},
                                           {2, 1, 4},
                                           {4, 1, 2},
                                           {1, 2, 1},
                                           {4, 2, 6},
                                           {0, 3, 2},
                                           {3, 3, 5},
                                           {0, 4, 1},
                                           {2, 4, 1},
                                           {4, 4, 7}}};
    const std::array rhs{-3.0, 6.0, -3.0, -19.0, 49.0};
    const std::array expected{1.0, -2.0, 3.0, -4.0, 5.0};
    check_solve(5, entries, rhs, expected, true);
}

TEST_CASE("Sparse triangular solve keeps a tiny nonzero pivot", "[N4-05]") {
    const double tiny = std::ldexp(1.0, -500);
    const std::array<Triplet, 2> entries{{{0, 0, 1}, {1, 1, tiny}}};
    const std::array rhs{1.0, tiny};
    const std::array expected{1.0, 1.0};
    check_solve(2, entries, rhs, expected, false);
}
