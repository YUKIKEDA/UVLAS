#include "linear_solve_detail.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <span>
#include <utility>

namespace {
using namespace uvlas;
using namespace uvlas::numerics::linear;
using Big = boost::multiprecision::cpp_bin_float_100;
using Complex = std::complex<double>;
namespace metrics = uvlas::numerics::linear::detail;

Accuracy accuracy(Norm norm = Norm::Infinity, double absolute = 1, double relative = 0,
                  double backward = 1) {
    return {norm, absolute, relative, backward, core::ObservationId{1}, core::ObservationId{2}};
}

template <typename Scalar, std::size_t Size>
auto column(const std::array<Scalar, Size>& residual, const std::array<double, Size>& error,
            const std::array<double, Size>& denominator, const std::array<Scalar, Size>& rhs,
            const std::array<double, Size>& scales, const Accuracy& request,
            std::span<const Scalar> low = {}) {
    return metrics::evaluate_column_metrics<Scalar>(static_cast<core::Index>(Size), residual, error,
                                                    denominator, rhs, low, scales, request);
}

void encloses(double lower, double upper, const Big& exact) {
    REQUIRE(std::isfinite(lower));
    REQUIRE(std::isfinite(upper));
    CHECK(lower >= 0);
    CHECK(Big(lower) <= exact);
    CHECK(Big(upper) >= exact);
    // Fixed envelope-width check: bounds cannot pass by returning an arbitrary large number.
    const Big allowance = abs(exact) * Big(64 * std::numeric_limits<double>::epsilon()) +
                          Big(8 * std::numeric_limits<double>::denorm_min());
    CHECK(Big(upper) - Big(lower) <= allowance);
}
} // namespace

TEST_CASE("Actual two-word magnitude encloses cancellation and extreme exponents",
          "[verification][negative][N4-06]") {
    const std::array pairs{std::pair{1.0, -1.0},
                           std::pair{1.0, std::ldexp(1.0, -54)},
                           std::pair{1.0, -std::nextafter(1.0, 0.0)},
                           std::pair{std::ldexp(1.0, 1000), -std::ldexp(1.0, 999)},
                           std::pair{std::numeric_limits<double>::denorm_min(), 0.0},
                           std::pair{-0.0, 0.0}};
    for (const auto& [hi, lo] : pairs) {
        const auto bound = metrics::actual_magnitude_bounds(hi, lo);
        REQUIRE(bound);
        encloses(bound->lower, bound->upper, abs(Big(hi) + Big(lo)));
    }
    const auto cancellation = metrics::actual_magnitude_bounds(1.0, -1.0);
    REQUIRE(cancellation);
    CHECK(cancellation->lower == 0);
    CHECK(cancellation->upper == 0);

    const std::array complex_pairs{std::pair{Complex{3, 4}, Complex{0, 0}},
                                   std::pair{Complex{1, 1}, Complex{-1, -1}},
                                   std::pair{Complex{1, 0}, Complex{-1, std::ldexp(1.0, -100)}},
                                   std::pair{Complex{1.e308, 1.e308}, Complex{0, 0}},
                                   std::pair{Complex{std::numeric_limits<double>::denorm_min(),
                                                     std::numeric_limits<double>::denorm_min()},
                                             Complex{0, 0}}};
    for (const auto& [hi, lo] : complex_pairs) {
        const Big real = Big(hi.real()) + Big(lo.real());
        const Big imag = Big(hi.imag()) + Big(lo.imag());
        const auto bound = metrics::actual_magnitude_bounds(hi, lo);
        REQUIRE(bound);
        encloses(bound->lower, bound->upper, sqrt(real * real + imag * imag));
    }
}

TEST_CASE("Residual criteria use model scales and actual RHS rather than high word only",
          "[verification][negative][N4-05][N4-06]") {
    const std::array residual{0.25, 0.125};
    const std::array error{0.0, 0.0};
    const std::array denominator{2.0, 2.0};
    const std::array rhs{1.0, 2.0};
    const std::array low{-1.0, -2.0};
    const std::array scales{0.5, 0.25};
    const auto request = accuracy(Norm::Infinity, 0.5, 1, 0.2);
    const auto value =
        column(residual, error, denominator, rhs, scales, request, std::span<const double>{low});
    REQUIRE(value);
    // Actual b is exactly zero: tau=(1/4,1/8), so the normalized maximum is one.
    CHECK(value->normalized_original_residual_upper >= 1);
    CHECK(value->normalized_original_residual_upper <=
          1 + 64 * std::numeric_limits<double>::epsilon());
    CHECK(value->componentwise_backward_error_upper >= 0.125);
    CHECK(value->componentwise_backward_error_upper <=
          0.125 + 64 * std::numeric_limits<double>::epsilon());
}

TEST_CASE("Residual and backward criteria are independent and include evaluation error",
          "[verification][negative][N4-05]") {
    const std::array residual{0.125};
    const std::array error{0.125};
    const std::array denominator{0.5};
    const std::array rhs{0.0};
    const std::array scales{1.0};
    const auto request = accuracy(Norm::Infinity, 1, 0, 0.4);
    const auto value = column(residual, error, denominator, rhs, scales, request);
    REQUIRE(value);
    CHECK(value->normalized_original_residual_upper >= 0.25);
    CHECK(value->normalized_original_residual_upper < 0.3);
    CHECK(value->componentwise_backward_error_upper >= 0.5);
    CHECK(value->componentwise_backward_error_upper < 0.6);
    CHECK_FALSE(metrics::meets_accuracy(*value, request));

    auto synthetic = *value;
    synthetic.normalized_original_residual_upper = 1;
    synthetic.componentwise_backward_error_upper = 0.4;
    CHECK(metrics::meets_accuracy(synthetic, request));
    synthetic.normalized_original_residual_upper = std::nextafter(1.0, 2.0);
    CHECK_FALSE(metrics::meets_accuracy(synthetic, request));
    synthetic.normalized_original_residual_upper = 0;
    synthetic.componentwise_backward_error_upper = std::nextafter(0.4, 1.0);
    CHECK_FALSE(metrics::meets_accuracy(synthetic, request));
}

TEST_CASE("L2 criteria combine components without extreme squares or max-norm substitution",
          "[verification][negative][N4-05]") {
    const std::array error{0.0, 0.0};
    const std::array denominator{1.e308, 1.e308};
    const std::array rhs{0.0, 0.0};
    const std::array scales{1.0, 1.0};
    const auto request = accuracy(Norm::L2, 5, 0, 1);
    const auto value = column(std::array{3.0, 4.0}, error, denominator, rhs, scales, request);
    REQUIRE(value);
    CHECK(value->normalized_original_residual_upper >= 1);
    CHECK(value->normalized_original_residual_upper <
          1 + 64 * std::numeric_limits<double>::epsilon());

    const auto large = column(std::array{6.e307, 8.e307}, error, denominator, rhs, scales,
                              accuracy(Norm::L2, 1.e308));
    REQUIRE(large);
    CHECK(large->normalized_original_residual_upper >= 0.99999999999999);
    CHECK(large->normalized_original_residual_upper <= 1.00000000000001);

    const auto tiny =
        column(std::array{1.e-200, 1.e-200}, error, denominator, rhs, scales, accuracy(Norm::L2));
    REQUIRE(tiny);
    CHECK(tiny->normalized_original_residual_upper > 1.e-200);
    CHECK(tiny->normalized_original_residual_upper < 1.5e-200);
}

TEST_CASE("Complex residual uses modulus rather than a single component",
          "[verification][negative][N4-05]") {
    const std::array residual{Complex{3, 4}};
    const std::array error{0.0};
    const std::array denominator{10.0};
    const std::array rhs{Complex{0, 0}};
    const std::array scales{10.0};
    const auto value = column(residual, error, denominator, rhs, scales, accuracy());
    REQUIRE(value);
    CHECK(value->normalized_original_residual_upper >= 0.5);
    CHECK(value->normalized_original_residual_upper < 0.50000000000001);
    CHECK(value->componentwise_backward_error_upper >= 0.5);
    CHECK(value->componentwise_backward_error_upper < 0.50000000000001);
}

TEST_CASE("Zero allowances require proven zero and positive quotient underflow stays positive",
          "[verification][negative][N4-05][N4-06]") {
    const std::array zero{0.0};
    const std::array scales{1.0};
    const auto exact = column(zero, zero, zero, zero, scales, accuracy(Norm::Infinity, 0, 0, 0));
    REQUIRE(exact);
    CHECK(exact->normalized_original_residual_upper == 0);
    CHECK(exact->componentwise_backward_error_upper == 0);
    CHECK(metrics::meets_accuracy(*exact, accuracy(Norm::Infinity, 0, 0, 0)));

    const double least = std::numeric_limits<double>::denorm_min();
    const auto uncertain =
        column(zero, std::array{least}, scales, zero, scales, accuracy(Norm::Infinity, 0, 0, 0));
    CHECK_FALSE(
        (uncertain && metrics::meets_accuracy(*uncertain, accuracy(Norm::Infinity, 0, 0, 0))));
    const auto underflow = column(std::array{least}, zero, std::array{2.0}, zero, scales,
                                  accuracy(Norm::Infinity, 2, 0, 0));
    REQUIRE(underflow);
    CHECK(underflow->normalized_original_residual_upper > 0);
    CHECK(underflow->componentwise_backward_error_upper > 0);
    CHECK_FALSE(metrics::meets_accuracy(*underflow, accuracy(Norm::Infinity, 2, 0, 0)));

    // The exact allowance is half the least binary64 value; rounding it upward would falsely pass.
    const auto tiny_allowance = column(std::array{least}, zero, scales, zero, std::array{0.5},
                                       accuracy(Norm::Infinity, least));
    CHECK_FALSE((tiny_allowance &&
                 metrics::meets_accuracy(*tiny_allowance, accuracy(Norm::Infinity, least))));
}

TEST_CASE("Metric helpers reject malformed spans and nonfinite or negative inputs",
          "[negative][N4-05][N4-06]") {
    const std::array zero{0.0};
    const std::array one{1.0};
    const std::array minus{-1.0};
    const std::array infinity{std::numeric_limits<double>::infinity()};
    CHECK_FALSE(column(zero, minus, one, zero, one, accuracy()));
    CHECK_FALSE(column(zero, zero, minus, zero, one, accuracy()));
    CHECK_FALSE(column(zero, zero, one, zero, zero, accuracy()));
    CHECK_FALSE(column(infinity, zero, one, zero, one, accuracy()));
    CHECK_FALSE(column(zero, infinity, one, zero, one, accuracy()));
    CHECK_FALSE(column(zero, zero, one, infinity, one, accuracy()));
    CHECK_FALSE(column(zero, zero, one, zero, one, accuracy(Norm::Infinity, -1)));
    CHECK_FALSE(column(zero, zero, one, zero, one, accuracy(static_cast<Norm>(255))));
    CHECK_FALSE(
        metrics::evaluate_column_metrics<double>(1, zero, zero, one, zero, {}, {}, accuracy()));
    CHECK_FALSE(metrics::actual_magnitude_bounds(std::numeric_limits<double>::infinity(), 0.0));
    CHECK_FALSE(metrics::actual_magnitude_bounds(
        Complex{0, 0}, Complex{0, std::numeric_limits<double>::quiet_NaN()}));
}

TEST_CASE("Relative complex tolerance uses actual modulus instead of component sum",
          "[verification][negative][N4-05][N4-06]") {
    const auto request = accuracy(Norm::Infinity, 0, 1, 1);
    const auto value = column(std::array{Complex{6, 0}}, std::array{0.0}, std::array{10.0},
                              std::array{Complex{3, 4}}, std::array{1.0}, request);
    REQUIRE(value);
    // |b|=5, whereas |Re(b)|+|Im(b)|=7 would falsely qualify this residual.
    CHECK(value->normalized_original_residual_upper >= 1.19999999999999);
    CHECK(value->normalized_original_residual_upper <= 1.20000000000001);
    CHECK_FALSE(metrics::meets_accuracy(*value, request));
}

TEST_CASE("Relative allowance includes a nonzero low RHS word",
          "[verification][negative][N4-05][N4-06]") {
    const std::array low{1.0};
    const auto value =
        column(std::array{0.5}, std::array{0.0}, std::array{1.0}, std::array{1.0}, std::array{1.0},
               accuracy(Norm::Infinity, 0, 0.5), std::span<const double>{low});
    REQUIRE(value);
    // Actual b=2 gives tau=1; ignoring low gives tau=1/2, and ignoring rtol gives zero.
    CHECK(value->normalized_original_residual_upper >= 0.5);
    CHECK(value->normalized_original_residual_upper < 0.50000000000001);
    CHECK(metrics::meets_accuracy(*value, accuracy(Norm::Infinity, 0, 0.5)));
}
