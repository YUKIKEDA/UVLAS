#include "compensated.hpp"
#include "linear_solve_detail.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <complex>
#include <limits>
#include <span>
#include <stdexcept>
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

TEST_CASE("A lower-only backward denominator retains the directed product and sum",
          "[verification][N4-05][N4-06]") {
    const auto eft = uvlas::numerics::detail::CheckedEft::prepare();
    REQUIRE(eft);
    constexpr double maximum = std::numeric_limits<double>::max();
    constexpr double subnormal = std::numeric_limits<double>::denorm_min();
    const std::array inputs{std::array{0., 2., 0.},
                            std::array{0., subnormal, .5},
                            std::array{0., maximum, 2.},
                            std::array{maximum, maximum, 1.},
                            std::array{.1, std::nextafter(1., 2.), .2},
                            std::array{1., 2., 3.}};
    for (const auto& [base, a, b] : inputs) {
        const double lower = metrics::add_magnitude_product_lower(base, a, b, *eft);
        const auto both =
            metrics::add_magnitudes({base, base}, metrics::multiply_magnitudes({a, a}, {b, b}));
        CHECK(std::bit_cast<std::uint64_t>(lower) == std::bit_cast<std::uint64_t>(both.lower));
        const Big exact = Big(base) + Big(a) * Big(b);
        CHECK(Big(lower) <= exact);
    }
}

TEST_CASE("Thread-local EFT qualification retains arithmetic words and failures",
          "[verification][negative][N4-05][N4-06]") {
    const auto eft = uvlas::numerics::detail::CheckedEft::prepare();
    REQUIRE(eft);
    const auto same = [](const auto& checked, const auto& full) {
        REQUIRE(checked.has_value() == full.has_value());
        if (!full) {
            CHECK(checked.error().code == full.error().code);
            return;
        }
        CHECK(std::bit_cast<std::uint64_t>(checked->hi) == std::bit_cast<std::uint64_t>(full->hi));
        CHECK(std::bit_cast<std::uint64_t>(checked->lo) == std::bit_cast<std::uint64_t>(full->lo));
    };
    constexpr auto tiny = std::numeric_limits<double>::denorm_min();
    constexpr auto maximum = std::numeric_limits<double>::max();
    const std::array pairs{std::pair{1., 0x1p-53},
                           std::pair{-12.25, 8.2},
                           std::pair{-0., +0.},
                           std::pair{tiny, .5},
                           std::pair{maximum, 2.},
                           std::pair{maximum, maximum},
                           std::pair{std::numeric_limits<double>::infinity(), 1.}};
    for (const auto& [a, b] : pairs) {
        same(eft->sum(a, b), uvlas::numerics::two_sum(a, b));
        same(eft->product(a, b), uvlas::numerics::two_product(a, b));
    }

    const auto previous = std::fegetround();
    REQUIRE(previous == FE_TONEAREST);
    REQUIRE(std::fesetround(FE_DOWNWARD) == 0);
    const auto rejected = uvlas::numerics::detail::CheckedEft::prepare();
    const auto restored = std::fesetround(previous);
    REQUIRE(restored == 0);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::PrecisionUnsupported);
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

namespace {
using ExactInteger = boost::multiprecision::cpp_int;
using ExactRational = boost::multiprecision::cpp_rational;
struct ExactSqrtInterval {
    ExactSqrtInterval(ExactRational low, ExactRational high)
        : lower(std::move(low)), upper(std::move(high)) {}
    ExactRational lower, upper;
};
ExactRational exact_power_two(int exponent) {
    // Oracle inputs are finite binary64 words. The wider signed type also
    // prevents negation overflow when obtaining a negative exponent's magnitude.
    const auto magnitude =
        exponent < 0 ? -static_cast<long long>(exponent) : static_cast<long long>(exponent);
    if (!std::in_range<unsigned>(magnitude)) {
        throw std::overflow_error("oracle exponent exceeds integer shift range");
    }
    const ExactInteger power = ExactInteger(1) << static_cast<unsigned>(magnitude);
    if (exponent < 0) {
        return ExactRational(1) / ExactRational(power);
    }
    return ExactRational(power);
}
ExactSqrtInterval sqrt_interval_256(const ExactRational& squared) {
    if (squared < 0) {
        throw std::invalid_argument("negative squared oracle magnitude");
    }
    if (squared == 0) {
        return {0, 0};
    }
    const ExactInteger numerator = boost::multiprecision::numerator(squared);
    const ExactInteger denominator = boost::multiprecision::denominator(squared);
    const auto numerator_bit = boost::multiprecision::msb(numerator);
    const auto denominator_bit = boost::multiprecision::msb(denominator);
    if (!std::in_range<int>(numerator_bit) || !std::in_range<int>(denominator_bit)) {
        throw std::overflow_error("oracle radicand exponent exceeds signed range");
    }
    int exponent = static_cast<int>(numerator_bit) - static_cast<int>(denominator_bit);
    if (squared < exact_power_two(exponent)) {
        --exponent;
    }
    // Exact power-of-two normalization makes the reference interval narrower
    // than a binary64 ulp even for physical minsub and near-overflow magnitudes.
    int half_exponent = exponent / 2;
    if (exponent % 2 < 0) {
        --half_exponent;
    }
    const ExactRational normalized = squared / exact_power_two(2 * half_exponent);
    const ExactInteger scale = ExactInteger(1) << 256;
    const ExactInteger scaled_numerator =
        boost::multiprecision::numerator(normalized) * scale * scale;
    const ExactInteger scaled_denominator = boost::multiprecision::denominator(normalized);
    const ExactInteger radicand = scaled_numerator / scaled_denominator;
    const ExactInteger root = boost::multiprecision::sqrt(radicand);
    const ExactRational low = ExactRational(root) / ExactRational(scale);
    ExactRational high = ExactRational(root + 1) / ExactRational(scale);
    if (root * root * scaled_denominator == scaled_numerator) {
        high = low;
    }
    const auto physical_scale = exact_power_two(half_exponent);
    return {low * physical_scale, high * physical_scale};
}
void encloses_exact_sqrt(double lower, double upper, const ExactRational& squared) {
    const auto reference = sqrt_interval_256(squared);
    REQUIRE(std::isfinite(lower));
    REQUIRE(std::isfinite(upper));
    REQUIRE(lower >= 0);
    REQUIRE(upper >= lower);
    // Check the original rational inequalities as well as the independently
    // computed 256-bit bracket; no product sqrt or metric code builds the oracle.
    CHECK(ExactRational(lower) * ExactRational(lower) <= squared);
    CHECK(ExactRational(upper) * ExactRational(upper) >= squared);
    CHECK(ExactRational(lower) <= reference.lower);
    CHECK(ExactRational(upper) >= reference.upper);
    const ExactRational allowance =
        64 * ExactRational(std::numeric_limits<double>::epsilon()) * reference.lower +
        8 * ExactRational(std::numeric_limits<double>::denorm_min());
    // This is the existing binary64 envelope limit, not a 2^-256 product-width
    // requirement. The latter precision belongs to the independent oracle only.
    CHECK(ExactRational(upper) - ExactRational(lower) <= allowance);
}
} // namespace

TEST_CASE("256-bit rational sqrt oracle certifies product two-word complex magnitude bounds",
          "[verification][negative][N4-06][N4-10]") {
    const ExactRational original = ExactRational(3) / 2;
    const auto interval = sqrt_interval_256(original);
    CHECK(interval.lower * interval.lower <= original);
    CHECK(interval.upper * interval.upper >= original);
    CHECK(interval.upper - interval.lower == exact_power_two(-256));

    const double tiny = std::numeric_limits<double>::denorm_min();
    const std::array fixtures{std::pair{Complex{1, .5}, Complex{}},
                              std::pair{Complex{3, 4}, Complex{}},
                              std::pair{Complex{1, 1}, Complex{-1, -1}},
                              std::pair{Complex{1, 0}, Complex{0, 0x1p-1000}},
                              std::pair{Complex{1, 0}, Complex{-1, 0x1p-1000}},
                              std::pair{Complex{0x1p1000, 0x1p1000}, Complex{0x1p946, -0x1p945}},
                              std::pair{Complex{tiny, tiny}, Complex{}},
                              std::pair{Complex{std::numeric_limits<double>::min(), tiny},
                                        Complex{-std::numeric_limits<double>::min(), 0}},
                              std::pair{Complex{-0., 0.}, Complex{}}};
    for (std::size_t i = 0; i < fixtures.size(); ++i) {
        CAPTURE(i);
        const auto& [high, low] = fixtures[i];
        const ExactRational real = ExactRational(high.real()) + ExactRational(low.real());
        const ExactRational imaginary = ExactRational(high.imag()) + ExactRational(low.imag());
        const ExactRational squared = real * real + imaginary * imaginary;
        const auto actual = metrics::actual_magnitude_bounds(high, low);
        REQUIRE(actual.has_value());
        encloses_exact_sqrt(actual->lower, actual->upper, squared);
    }
}

TEST_CASE("256-bit sqrt three halves reference certifies the actual product L2 residual path",
          "[verification][negative][N4-05][N4-06][N4-10]") {
    const ExactRational original = ExactRational(3) / 2;
    const auto reference = sqrt_interval_256(original);
    REQUIRE(reference.upper - reference.lower == exact_power_two(-256));
    // Six exact .5 entries have squared norm 3/2. Scaling both model scales
    // and inputs by an exact power of two must preserve that same original
    // normalized criterion, including subnormal and very large physical inputs.
    for (const int exponent : {-1073, -1000, 0, 1000, 1023}) {
        CAPTURE(exponent);
        const double scale = std::scalbn(1., exponent);
        std::array<double, 6> residual{}, scales{}, zero{};
        residual.fill(scale / 2);
        scales.fill(scale);
        for (const double absolute : {1., 2.}) {
            CAPTURE(absolute);
            const auto request = accuracy(Norm::L2, absolute, 0, 1);
            const auto actual = column(residual, zero, scales, zero, scales, request);
            REQUIRE(actual.has_value());
            const ExactRational squared =
                original / (ExactRational(absolute) * ExactRational(absolute));
            const auto expected = sqrt_interval_256(squared);
            const auto upper = actual->normalized_original_residual_upper;
            REQUIRE(std::isfinite(upper));
            CHECK(ExactRational(upper) * ExactRational(upper) >= squared);
            CHECK(ExactRational(upper) >= expected.upper);
            const ExactRational allowance =
                64 * ExactRational(std::numeric_limits<double>::epsilon()) * expected.lower +
                8 * ExactRational(std::numeric_limits<double>::denorm_min());
            CHECK(ExactRational(upper) - expected.lower <= allowance);
            CHECK(metrics::meets_accuracy(*actual, request) == (absolute == 2));
        }
    }
}
