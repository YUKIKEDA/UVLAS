#include <Eigen/Core>
#include <Spectra/GenEigsSolver.h>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <complex>

namespace {
// No assembled matrix is handed to Spectra. The dominant block has roots +/-5i.
struct RealOperator {
    using Scalar = double;
    [[nodiscard]] static Eigen::Index rows() {
        return 8;
    }
    [[nodiscard]] static Eigen::Index cols() {
        return 8;
    }
    static void perform_op(const double* input, double* output) {
        output[0] = -5 * input[1];
        output[1] = 5 * input[0];
        for (Eigen::Index i = 2; i < 8; ++i) {
            output[i] = (static_cast<double>(i) / 10) * input[i];
        }
    }
};

struct ComplexOperator {
    using Scalar = std::complex<double>;
    [[nodiscard]] static Eigen::Index rows() {
        return 6;
    }
    [[nodiscard]] static Eigen::Index cols() {
        return 6;
    }
    static void perform_op(const Scalar* input, Scalar* output) {
        // Upper triangular, with distinct eigenvalues known without a solver.
        constexpr std::array<Scalar, 6> diagonal{
            {{8, 2}, {6, -1}, {0, 4}, {2, 0}, {1, 0}, {.5, 0}}};
        for (Eigen::Index i = 0; i < 6; ++i) {
            output[i] = diagonal[static_cast<std::size_t>(i)] * input[i];
            if (i + 1 < 6) {
                output[i] += Scalar(0, .25) * input[i + 1];
            }
        }
    }
};
} // namespace

TEST_CASE("real operator Arnoldi returns the complex conjugate dominant pair",
          "[toolchain][N4-01]") {
    RealOperator op;
    Spectra::GenEigsSolver<RealOperator> solver(op, 2, 8);
    const Eigen::Vector<double, 8> initial = Eigen::Vector<double, 8>::LinSpaced(1, 8);
    solver.init(initial.data());
    REQUIRE(solver.compute(Spectra::SortRule::LargestMagn, 100, 1e-13) == 2);
    REQUIRE(solver.info() == Spectra::CompInfo::Successful);
    const auto values = solver.eigenvalues();
    const auto vectors = solver.eigenvectors();
    REQUIRE(values.size() == 2);
    for (Eigen::Index j = 0; j < 2; ++j) {
        CHECK(std::abs(values[j].real()) <= 1e-12);
        CHECK(std::abs(std::abs(values[j].imag()) - 5) <= 1e-12);
        Eigen::Vector<std::complex<double>, 8> applied;
        applied[0] = -5.0 * vectors(1, j);
        applied[1] = 5.0 * vectors(0, j);
        for (Eigen::Index i = 2; i < 8; ++i) {
            applied[i] = (static_cast<double>(i) / 10) * vectors(i, j);
        }
        CHECK((applied - values[j] * vectors.col(j)).norm() <= 1e-12);
    }
    CHECK(std::abs(values[0] - std::conj(values[1])) <= 1e-12);
}

TEST_CASE("complex operator Arnoldi solves a non Hermitian triangular operator",
          "[toolchain][N4-01]") {
    ComplexOperator op;
    Spectra::GenEigsSolver<ComplexOperator> solver(op, 2, 6);
    const Eigen::VectorXcd initial = Eigen::VectorXcd::Ones(6);
    solver.init(initial.data());
    REQUIRE(solver.compute(Spectra::SortRule::LargestMagn, 100, 1e-13) == 2);
    REQUIRE(solver.info() == Spectra::CompInfo::Successful);
    const auto values = solver.eigenvalues();
    const auto vectors = solver.eigenvectors();
    REQUIRE(values.size() == 2);
    CHECK(std::abs(values[0] - std::complex<double>(8, 2)) <= 1e-12);
    CHECK(std::abs(values[1] - std::complex<double>(6, -1)) <= 1e-12);
    // Check original components, separately from perform_op's diagonal array.
    for (Eigen::Index j = 0; j < 2; ++j) {
        CHECK(std::abs(std::complex<double>(8, 2) * vectors(0, j) +
                       std::complex<double>(0, .25) * vectors(1, j) - values[j] * vectors(0, j)) <=
              1e-12);
        CHECK(std::abs(std::complex<double>(6, -1) * vectors(1, j) +
                       std::complex<double>(0, .25) * vectors(2, j) - values[j] * vectors(1, j)) <=
              1e-12);
        CHECK(vectors.col(j).tail(4).norm() <= 1e-12);
    }
}
