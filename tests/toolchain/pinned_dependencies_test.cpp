// Checks that the headers reaching the compiler are the versions named in
// vcpkg.json and that the build keeps the floating-point
// model the adopted numerics rely on (N4-01, N4-06). This is a toolchain check,
// not a verification of any solver path.

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Spectra/Util/Version.h>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/version.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdint>
#include <limits>

TEST_CASE("pinned dependency versions reach the compiler", "[toolchain][N4-01]") {
    // Eigen 5 keeps EIGEN_WORLD_VERSION at 3; the release is MAJOR.MINOR.PATCH.
    CHECK(EIGEN_MAJOR_VERSION == 5);
    CHECK(EIGEN_MINOR_VERSION == 0);
    CHECK(EIGEN_PATCH_VERSION == 1);
    CHECK(SPECTRA_MAJOR_VERSION == 1);
    CHECK(SPECTRA_MINOR_VERSION == 2);
    CHECK(SPECTRA_PATCH_VERSION == 0);
    CHECK(BOOST_VERSION == 109200);
    CHECK(NLOHMANN_JSON_VERSION_MAJOR == 3);
    CHECK(NLOHMANN_JSON_VERSION_MINOR == 12);
    CHECK(NLOHMANN_JSON_VERSION_PATCH == 0);
}

TEST_CASE("numerical build policy is in effect", "[toolchain][N4-01][N4-02]") {
#ifndef EIGEN_MPL2_ONLY
    FAIL("EIGEN_MPL2_ONLY is not defined");
#endif
    // N4-02: indices are signed 64-bit on every supported platform.
    STATIC_REQUIRE(sizeof(Eigen::Index) == 8);
    STATIC_REQUIRE(std::numeric_limits<Eigen::Index>::is_signed);
    STATIC_REQUIRE(std::numeric_limits<double>::is_iec559);
}

TEST_CASE("fma gives the exact product error", "[toolchain][N4-06]") {
    // Error-free two-product: p + e == a * b exactly. A build that contracts or
    // reorders floating-point expressions breaks the compensated residuals.
    const volatile double a = 1.0 + std::ldexp(1.0, -30);
    const volatile double b = 1.0 - std::ldexp(1.0, -31);
    const double p = a * b;
    const double e = std::fma(a, b, -p);

    using Big = boost::multiprecision::cpp_bin_float_100;
    const Big exact = Big(static_cast<double>(a)) * Big(static_cast<double>(b));
    CHECK(Big(p) + Big(e) == exact);
    CHECK(e != 0.0);
}
