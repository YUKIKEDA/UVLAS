// Numeric canary (QA-03): fixed floating-point work whose results are stored
// bit for bit per platform. A change of compiler, flags, instruction set,
// standard library or dependency that changes any of them is caught here,
// before it reaches a solver result. Output: one "name value" line each, the
// value in hexadecimal floating point, preceded by the platform key.

#include "uvlas/core/build_info.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include <cmath>
#include <cstdint>
#include <exception>
#include <format>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

void emit(std::string_view name, double value) {
    std::cout << std::format("{} {:a}\n", name, value);
}

// Error-free product (N4-06): p + e equals a * b exactly.
void two_product() {
    const double a = 1.0 + std::ldexp(1.0, -30);
    const double b = 1.0 - std::ldexp(1.0, -31);
    const double p = a * b;
    emit("two_product.p", p);
    emit("two_product.e", std::fma(a, b, -p));
}

// A sum that cancels: plain and compensated (Neumaier).
void sums() {
    std::vector<double> terms;
    for (int k = 1; k <= 2000; ++k) {
        const double x = 1.0 / static_cast<double>(k);
        terms.push_back((k % 2 == 0) ? -x * 1.0e8 : x * 1.0e8);
        terms.push_back(std::sin(static_cast<double>(k)) * 1.0e-8);
    }
    double plain = 0.0;
    double sum = 0.0;
    double correction = 0.0;
    for (const double x : terms) {
        plain += x;
        const double t = sum + x;
        correction += (std::abs(sum) >= std::abs(x)) ? (sum - t) + x : (x - t) + sum;
        sum = t;
    }
    emit("sum.plain", plain);
    emit("sum.compensated", sum + correction);
}

// Dense LU with partial pivoting on a mildly ill-conditioned matrix.
void dense_solve() {
    constexpr Eigen::Index n = 12;
    Eigen::MatrixXd a(n, n);
    Eigen::VectorXd b(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        b(i) = 1.0 + 0.125 * static_cast<double>(i);
        for (Eigen::Index j = 0; j < n; ++j) {
            a(i, j) = 1.0 / static_cast<double>(i + j + 1) + ((i == j) ? 1.0e-3 : 0.0);
        }
    }
    const Eigen::VectorXd x = a.partialPivLu().solve(b);
    emit("dense.x0", x(0));
    emit("dense.x_last", x(n - 1));
    emit("dense.residual_norm", (a * x - b).norm());
}

// Sparse LU on a small nonsymmetric band matrix.
void sparse_solve() {
    constexpr Eigen::Index n = 40;
    std::vector<Eigen::Triplet<double>> entries;
    for (Eigen::Index i = 0; i < n; ++i) {
        entries.emplace_back(i, i, 4.0 + 0.01 * static_cast<double>(i));
        if (i + 1 < n) {
            entries.emplace_back(i, i + 1, -1.0);
            entries.emplace_back(i + 1, i, -1.5);
        }
        if (i + 7 < n) {
            entries.emplace_back(i, i + 7, 0.25);
        }
    }
    Eigen::SparseMatrix<double> a(n, n);
    a.setFromTriplets(entries.begin(), entries.end());
    a.makeCompressed();
    const Eigen::VectorXd b = Eigen::VectorXd::LinSpaced(n, 1.0, 2.0);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
    lu.compute(a);
    const Eigen::VectorXd x = lu.solve(b);
    emit("sparse.x0", x(0));
    emit("sparse.x_last", x(n - 1));
    emit("sparse.residual_norm", (a * x - b).norm());
}

// Standard library functions are not correctly rounded and differ by vendor.
void math_library() {
    const double x = 0.7853981633974483;
    emit("libm.sin", std::sin(x));
    emit("libm.cos", std::cos(1.0e6));
    emit("libm.exp", std::exp(-12.5));
    emit("libm.log", std::log(3.0));
    emit("libm.pow", std::pow(1.0000001, 12345.678));
    emit("libm.atan2", std::atan2(1.0e-9, -1.0));
    emit("libm.hypot", std::hypot(3.0e-200, 4.0e-200));
    emit("libm.cbrt", std::cbrt(2.0));
    emit("libm.sqrt", std::sqrt(2.0));
}

// A 100-bit computation rounded to double: the reference path of N4-06.
void multiprecision() {
    using Big = boost::multiprecision::cpp_bin_float_100;
    const Big third = Big(1) / Big(3);
    emit("big.third", static_cast<double>(third));
    emit("big.sqrt2", static_cast<double>(sqrt(Big(2))));
    emit("big.cancel", static_cast<double>((Big(1) + Big("1e-40")) - Big(1)));
}

int run() {
    try {
        std::cout << "platform " << uvlas::build_info().platform_key << '\n';
        two_product();
        sums();
        dense_solve();
        sparse_solve();
        math_library();
        multiprecision();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

} // namespace

// No exception leaves the process, not even one thrown while reporting another.
int main() {
    try {
        return run();
    } catch (...) {
        return 1;
    }
}
