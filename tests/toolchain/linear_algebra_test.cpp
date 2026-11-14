// These exercise the pinned dependencies, not the future N4-05 product solvers.
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <catch2/catch_test_macros.hpp>
#include <unsupported/Eigen/IterativeSolvers>

#include <complex>
#include <cstdint>

TEST_CASE("pivoted LU solves real and complex nonsymmetric systems", "[toolchain][N4-01]") {
    Eigen::Matrix3d matrix;
    matrix << 0, 2, 0, 1, 3, 0, 0, 0, -4;
    const Eigen::Vector3d expected(1, 2, 3);
    const Eigen::Vector3d rhs(4, 7, -12);
    const Eigen::PartialPivLU<Eigen::Matrix3d> factor(matrix);
    const Eigen::Vector3d solution = factor.solve(rhs);
    // N4-05 small-system known solution and original residual, fixed 1e-12.
    CHECK((solution - expected).norm() <= 1e-12);
    CHECK((matrix * solution - rhs).norm() <= 1e-12);

    Eigen::Matrix3cd complex_matrix = matrix.cast<std::complex<double>>();
    complex_matrix(1, 1) += std::complex<double>(0, 1);
    const Eigen::Vector3cd complex_rhs(4, std::complex<double>(7, 2), -12);
    const Eigen::PartialPivLU<Eigen::Matrix3cd> complex_factor(complex_matrix);
    const Eigen::Vector3cd complex_solution = complex_factor.solve(complex_rhs);
    CHECK((complex_solution - expected.cast<std::complex<double>>()).norm() <= 1e-12);
    CHECK((complex_matrix * complex_solution - complex_rhs).norm() <= 1e-12);
}

TEST_CASE("64 bit CSC SparseLU with COLAMD solves an indefinite KKT", "[toolchain][N4-01]") {
    using Sparse = Eigen::SparseMatrix<double, Eigen::ColMajor, std::int64_t>;
    STATIC_REQUIRE(sizeof(Sparse::StorageIndex) == 8);
    Eigen::Matrix3d dense;
    dense << 2, 0, 1, 0, 3, 1, 1, 1, 0;
    Sparse matrix = dense.sparseView();
    matrix.makeCompressed();
    const Eigen::Vector3d expected(1, 2, -1);
    const Eigen::Vector3d rhs(1, 5, 3);
    Eigen::SparseLU<Sparse, Eigen::COLAMDOrdering<std::int64_t>> factor;
    factor.compute(matrix);
    REQUIRE(factor.info() == Eigen::Success);
    const Eigen::Vector3d solution = factor.solve(rhs);
    REQUIRE(factor.info() == Eigen::Success);
    CHECK((solution - expected).norm() <= 1e-12);
    CHECK((matrix * solution - rhs).norm() <= 1e-12);
}

TEST_CASE("dense and sparse LLT solve a known SPD system", "[toolchain][N4-01]") {
    // Symmetric with positive leading minors 4 and 11 (Sylvester criterion).
    Eigen::Matrix2d dense;
    dense << 4, 1, 1, 3;
    const Eigen::Vector2d rhs(6, 7);
    const Eigen::Vector2d expected(1, 2);
    const Eigen::LLT<Eigen::Matrix2d> factor(dense);
    REQUIRE(factor.info() == Eigen::Success);
    const Eigen::Vector2d solution = factor.solve(rhs);
    CHECK((solution - expected).norm() <= 1e-12);
    CHECK((dense * solution - rhs).norm() <= 1e-12);
    using Sparse = Eigen::SparseMatrix<double, Eigen::ColMajor, std::int64_t>;
    const Sparse matrix = dense.sparseView();
    const Eigen::SimplicialLLT<Sparse> sparse_factor(matrix);
    REQUIRE(sparse_factor.info() == Eigen::Success);
    const Eigen::Vector2d sparse_solution = sparse_factor.solve(rhs);
    CHECK((sparse_solution - expected).norm() <= 1e-12);
    CHECK((matrix * sparse_solution - rhs).norm() <= 1e-12);
}

TEST_CASE("GMRES with a fixed identity preconditioner solves a nonsymmetric system",
          "[toolchain][N4-01]") {
    // IterativeSolverBase's default dense Ref starts at 0x0, requiring a
    // dynamic matrix type even when this fixture has a fixed small size.
    Eigen::MatrixXd matrix(3, 3);
    matrix << 4, 1, 0, 0, 3, 1, 1, 0, 2;
    const Eigen::Vector3d rhs(6, 9, 7);
    const Eigen::Vector3d expected(1, 2, 3);
    Eigen::GMRES<Eigen::MatrixXd, Eigen::IdentityPreconditioner> solver;
    solver.set_restart(3);
    solver.setMaxIterations(12);
    solver.setTolerance(1e-14);
    solver.compute(matrix);
    const Eigen::Vector3d solution = solver.solve(rhs);
    REQUIRE(solver.info() == Eigen::Success);
    CHECK((solution - expected).norm() <= 1e-12);
    CHECK((matrix * solution - rhs).norm() <= 1e-12);
}
