#include "numerics/eigen_views.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <limits>

namespace {
using namespace uvlas::core;
using namespace uvlas::numerics;
namespace adapter = uvlas::numerics::detail;
constexpr Operation op = Operation::SolveLinear;
template <class T>
void invalid(const Result<T>& result) {
    REQUIRE_FALSE(result);
    CHECK(result.error().code == FailureCode::InvalidInput);
}
} // namespace
TEST_CASE("Private Eigen dense maps preserve padding and share storage", "[integration][N4-02]") {
    ResourcePlan plan{1024};
    const std::array<double, 6> source{1, 2, 99, 3, 4, 88};
    auto matrix = DenseReal::copy({2, 2, 3}, source, plan, op);
    REQUIRE(matrix);
    const auto mapped = adapter::eigen_read_map(*matrix, 2, 2, op);
    REQUIRE(mapped);
    CHECK(mapped->data() == matrix->values().data());
    CHECK((*mapped)(1, 1) == 4);
    CHECK(mapped->outerStride() == 3);
    invalid(adapter::eigen_read_map(*matrix, 1, 2, op));
    matrix->values()[1] = std::numeric_limits<double>::infinity();
    const auto bad = adapter::eigen_read_map(*matrix, 2, 2, op);
    REQUIRE_FALSE(bad);
    CHECK(bad.error().code == FailureCode::Nonfinite);
    const std::array<std::complex<double>, 2> values{{{1, 2}, {3, -4}}};
    auto complex = DenseComplex::copy({1, 2, 1}, values, plan, op);
    REQUIRE(complex);
    const auto complex_map = adapter::eigen_read_map(*complex, 1, 2, op);
    REQUIRE(complex_map);
    CHECK((*complex_map)(0, 1) == std::complex<double>{3, -4});
    const auto empty = DenseReal::copy({0, 3, 0}, {}, plan, op);
    REQUIRE(empty);
    REQUIRE(adapter::eigen_read_map(*empty, 0, 3, op));
}
TEST_CASE("Eigen rejects scalar layout stride and shape mismatches", "[negative][N4-02]") {
    Eigen::Matrix<double, 2, 2, Eigen::ColMajor> matrix;
    matrix << 1, 3, 2, 4;
    const auto borrowed = adapter::borrow_eigen_storage<double>(matrix, {2, 2, 2}, op);
    REQUIRE(borrowed);
    CHECK(borrowed->data() == matrix.data());
    CHECK((*borrowed)[1] == 2);
    invalid(adapter::borrow_eigen_storage<double>(matrix, {1, 2, 2}, op));
    const Eigen::Matrix<double, 2, 2, Eigen::RowMajor> row_major = matrix;
    invalid(adapter::borrow_eigen_storage<double>(row_major, {2, 2, 2}, op));
    const Eigen::Matrix<float, 2, 2, Eigen::ColMajor> floats = matrix.cast<float>();
    invalid(adapter::borrow_eigen_storage<double>(floats, {2, 2, 2}, op));
    std::array<double, 8> storage{1, 0, 2, 0, 3, 0, 4, 0};
    const Eigen::Map<const adapter::EigenMatrix<double>, Eigen::Unaligned,
                     Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>
        strided{storage.data(), 2, 2, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{4, 2}};
    invalid(adapter::borrow_eigen_storage<double>(strided, {2, 2, 4}, op));
    matrix(0, 0) = std::numeric_limits<double>::quiet_NaN();
    const auto bad = adapter::borrow_eigen_storage<double>(matrix, {2, 2, 2}, op);
    REQUIRE_FALSE(bad);
    CHECK(bad.error().code == FailureCode::Nonfinite);
}
TEST_CASE("Eigen complex views reject either nonfinite part", "[negative][N4-02]") {
    Eigen::Matrix<std::complex<double>, 2, 1> matrix;
    matrix << std::complex<double>{2, -3}, std::complex<double>{4, 5};
    const auto valid = adapter::borrow_eigen_storage<std::complex<double>>(matrix, {2, 1, 2}, op);
    REQUIRE(valid);
    CHECK(valid->data() == matrix.data());
    CHECK((*valid)[0] == std::complex<double>{2, -3});
    CHECK((*valid)[1] == std::complex<double>{4, 5});
    const double nonfinite = std::numeric_limits<double>::infinity();
    for (const auto value :
         {std::complex<double>{nonfinite, 1}, std::complex<double>{1, nonfinite}}) {
        matrix(0, 0) = value;
        const auto bad = adapter::borrow_eigen_storage<std::complex<double>>(matrix, {2, 1, 2}, op);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().code == FailureCode::Nonfinite);
    }
}
TEST_CASE("Private Eigen sparse maps borrow canonical int64 CSC arrays", "[integration][N4-02]") {
    ResourcePlan plan{1024};
    const std::array<Triplet, 2> entries{{{1, 0, 3}, {0, 1, 5}}};
    const auto matrix = SparseReal::from_triplets(2, 2, entries, plan, op);
    REQUIRE(matrix);
    const auto mapped = adapter::eigen_read_map(*matrix, 2, 2, op);
    REQUIRE(mapped);
    CHECK(mapped->valuePtr() == matrix->values().data());
    CHECK(mapped->coeff(1, 0) == 3);
    CHECK(mapped->coeff(0, 1) == 5);
    invalid(adapter::eigen_read_map(*matrix, 1, 2, op));
    const auto empty = SparseReal::from_triplets(0, 0, {}, plan, op);
    REQUIRE(empty);
    REQUIRE(adapter::eigen_read_map(*empty, 0, 0, op));
}
