#include <boost/multiprecision/cpp_bin_float.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <complex>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <uvlas/numerics/precision.hpp>
#include <vector>
#if defined(_M_X64) || defined(__x86_64__)
#include <xmmintrin.h>
#endif

namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using Big = boost::multiprecision::cpp_bin_float_100;
using C = std::complex<double>;
constexpr auto op = core::Operation::SolveLinear;
constexpr auto generation = core::NumericalGeneration{17};
template <StorageScalar S>
DenseMatrix<S> matrix(DenseShape shape, std::initializer_list<S> values, core::ResourcePlan& plan) {
    auto result =
        DenseMatrix<S>::copy(shape, std::span<const S>{values.begin(), values.size()}, plan, op);
    REQUIRE(result.has_value());
    return std::move(*result);
}
template <StorageScalar S>
class Verifier final : public ProjectionVerifier<S> {
public:
    mutable int calls = 0;
    bool unavailable = false;
    bool unsupported = false;
    bool throwing = false;
    bool require_exact_first = false;
    core::Result<ProjectionAssessment> verify(const PrecisionBuffer<S>& source,
                                              const DenseMatrix<S>& candidate,
                                              const ProjectionContext& context) const override {
        ++calls;
        CHECK(context.source == source.generation());
        CHECK(context.target == core::NumericalGeneration{18});
        CHECK(context.consumer == core::ScopeId{31});
        if (throwing) {
            throw std::runtime_error("test evaluator failure");
        }
        if (unsupported) {
            return std::unexpected(core::Failure{core::FailureCode::PrecisionUnsupported, op});
        }
        if (require_exact_first) {
            if constexpr (std::same_as<S, double>) {
                const auto& pair = std::get<TwoWord<S>>(source.storage());
                if (Big(pair.hi.values()[0]) + Big(pair.lo.values()[0]) !=
                    Big(candidate.values()[0])) {
                    return std::unexpected(
                        core::Failure{core::FailureCode::PrecisionUnresolved, op});
                }
            }
        }
        return unavailable ? ProjectionAssessment::Unevaluated : ProjectionAssessment::Verified;
    }
};
template <StorageScalar S>
auto project(const PrecisionBuffer<S>& source, const DenseMatrix<S>& a, const DenseMatrix<S>& b,
             std::span<const double> scales, double tolerance, const Verifier<S>& verifier,
             core::ResourcePlan& plan) {
    return project_binary64(source, a, b, scales, tolerance, core::NumericalGeneration{18},
                            core::ScopeId{31}, verifier, plan);
}
} // namespace

TEST_CASE("EFT preserves exact sums and products including subnormal boundaries",
          "[verification][N4-06]") {
    const std::array pairs{std::pair{1.0, std::ldexp(1.0, -54)},
                           std::pair{1.e16, 1.0},
                           std::pair{1.0 + std::ldexp(1.0, -27), 1.0 - std::ldexp(1.0, -27)},
                           std::pair{std::ldexp(1.0, -1022), std::ldexp(1.0, -52)},
                           std::pair{std::numeric_limits<double>::denorm_min(), 1.0},
                           std::pair{-0.0, 2.0}};
    for (const auto [a, b] : pairs) {
        const auto sum = two_sum(a, b);
        REQUIRE(sum);
        CHECK(Big(sum->hi) + Big(sum->lo) == Big(a) + Big(b));
        const auto product = two_product(a, b);
        REQUIRE(product);
        CHECK(Big(product->hi) + Big(product->lo) == Big(a) * Big(b));
    }
    for (int i = 1; i <= 40; ++i) {
        const double a = std::ldexp(1.0 + static_cast<double>(i) / 64, -300 + i);
        const double b = std::ldexp(1.0 - static_cast<double>(i) / 128, 200 - i);
        const auto value = two_product(a, b);
        REQUIRE(value);
        CHECK(Big(value->hi) + Big(value->lo) == Big(a) * Big(b));
    }
    // A 334-bit reference cannot resolve this exponent gap; check the exact components.
    for (const double tiny :
         {std::numeric_limits<double>::denorm_min(), -std::numeric_limits<double>::denorm_min()}) {
        const auto extreme = two_sum(1.0, tiny);
        REQUIRE(extreme);
        CHECK(extreme->hi == 1.0);
        CHECK(extreme->lo == tiny);
    }
    CHECK(std::signbit(two_product(-0.0, 2.0)->hi));
    const auto tiny = two_product(std::numeric_limits<double>::denorm_min(), 0.5);
    REQUIRE_FALSE(tiny);
    CHECK(tiny.error().code == core::FailureCode::PrecisionUnresolved);
    const auto overflow = two_product(std::numeric_limits<double>::max(), 2);
    REQUIRE_FALSE(overflow);
    CHECK(overflow.error().code == core::FailureCode::Nonfinite);
    for (const double value :
         {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        CHECK_FALSE(two_product(value, 1));
        CHECK_FALSE(two_product(1, value));
        CHECK_FALSE(two_sum(value, 1));
        CHECK_FALSE(two_sum(1, value));
    }
    const int previous = std::fegetround();
    REQUIRE(std::fesetround(FE_DOWNWARD) == 0);
    const auto refused = two_product(1, 1);
    const auto refused_sum = two_sum(1, 1);
    REQUIRE(std::fesetround(previous) == 0);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().code == core::FailureCode::PrecisionUnsupported);
    REQUIRE_FALSE(refused_sum);
    CHECK(refused_sum.error().code == core::FailureCode::PrecisionUnsupported);
}
TEST_CASE("Precision buffers retain padded shapes finite planes and promoted corrections",
          "[unit][negative][N4-02][N4-06]") {
    core::ResourcePlan plan{1024ULL * 1024};
    auto source =
        PrecisionBuffer<double>::binary64(matrix<double>({1, 1, 2}, {1, 9}, plan), generation);
    REQUIRE(source);
    const auto delta = matrix<double>({1, 1, 2}, {std::ldexp(1.0, -54), 0}, plan);
    auto promoted = add_correction(*source, delta, core::NumericalGeneration{18}, plan);
    REQUIRE(promoted);
    CHECK(promoted->precision() == core::Precision::TwoWord);
    const auto& pair = std::get<TwoWord<double>>(promoted->storage());
    CHECK(Big(pair.hi.values()[0]) + Big(pair.lo.values()[0]) == Big(1) + pow(Big(2), -54));
    auto again = add_correction(*promoted, delta, core::NumericalGeneration{19}, plan);
    REQUIRE(again);
    const auto& second = std::get<TwoWord<double>>(again->storage());
    CHECK(Big(second.hi.values()[0]) + Big(second.lo.values()[0]) == Big(1) + pow(Big(2), -53));
    CHECK(std::get<Binary64<double>>(source->storage()).values.values()[0] == 1);
    CHECK(pair.hi.shape().leading_dimension == 2);
    CHECK(pair.lo.values()[1] == 0);
    CHECK_FALSE(add_correction(*source, delta, generation, plan));
    auto wrong = PrecisionBuffer<double>::two_word(
        matrix<double>({1, 1, 1}, {1}, plan), matrix<double>({1, 1, 2}, {0, 0}, plan), generation);
    REQUIRE_FALSE(wrong);
    CHECK(wrong.error().code == core::FailureCode::InvalidInput);
    for (const C value : {C{std::numeric_limits<double>::infinity(), 0},
                          C{0, std::numeric_limits<double>::infinity()}}) {
        for (const bool lower : {false, true}) {
            auto hi = matrix<C>({1, 1, 1}, {C{1, 2}}, plan);
            auto lo = matrix<C>({1, 1, 1}, {C{0, 0}}, plan);
            (lower ? lo.values() : hi.values())[0] = value;
            auto invalid = PrecisionBuffer<C>::two_word(std::move(hi), std::move(lo), generation);
            REQUIRE_FALSE(invalid);
            CHECK(invalid.error().code == core::FailureCode::Nonfinite);
        }
    }
    for (const auto precision : {core::Precision::Binary64, core::Precision::TwoWord}) {
        CHECK(require_precision(precision, SupportedPrecision::Binary64AndTwoWord));
    }
    CHECK(require_precision(core::Precision::Binary64, SupportedPrecision::Binary64Only));
    const auto unsupported =
        require_precision(core::Precision::TwoWord, SupportedPrecision::Binary64Only);
    REQUIRE_FALSE(unsupported);
    CHECK(unsupported.error().code == core::FailureCode::PrecisionUnsupported);
    CHECK_FALSE(
        require_precision(static_cast<core::Precision>(250), SupportedPrecision::Binary64Only));
    CHECK_FALSE(require_precision(core::Precision::Binary64, static_cast<SupportedPrecision>(250)));
}
TEST_CASE("Compensated original residual consumes complex cross terms and both words",
          "[verification][N4-06]") {
    core::ResourcePlan plan{1024ULL * 1024};
    const auto a = matrix<C>(
        {2, 2, 3}, {C{1, 2}, C{-3, 4}, C{999, 999}, C{5, -6}, C{7, 8}, C{999, 999}}, plan);
    auto x = PrecisionBuffer<C>::two_word(
        matrix<C>({2, 1, 3}, {C{1, 2}, C{-2, 3}, C{99, 99}}, plan),
        matrix<C>({2, 1, 3}, {C{0x1p-54, -0x1p-55}, C{0x1p-53, 0x1p-52}, C{99, 99}}, plan),
        generation);
    REQUIRE(x);
    const auto b = matrix<C>({2, 1, 3}, {C{4, 5}, C{-6, 7}, C{888, 888}}, plan);
    const auto residual = original_residual(a, *x, b, plan);
    REQUIRE(residual);
    const auto& pair = std::get<TwoWord<C>>(x->storage());
    for (std::size_t row = 0; row < 2; ++row) {
        Big real = b.values()[row].real(), imag = b.values()[row].imag();
        for (std::size_t j = 0; j < 2; ++j) {
            const C coefficient = a.values()[3 * j + row];
            const Big xr = Big(pair.hi.values()[j].real()) + Big(pair.lo.values()[j].real());
            const Big xi = Big(pair.hi.values()[j].imag()) + Big(pair.lo.values()[j].imag());
            real -= Big(coefficient.real()) * xr - Big(coefficient.imag()) * xi;
            imag -= Big(coefficient.real()) * xi + Big(coefficient.imag()) * xr;
        }
        // N4-05 small-system 1e-12 absolute tolerance, checked independently with binary100.
        CHECK(abs(Big(residual->values()[row].real()) - real) < Big("1e-12"));
        CHECK(abs(Big(residual->values()[row].imag()) - imag) < Big("1e-12"));
    }
    CHECK(residual->values()[2] == C(888, 888));
    auto empty = PrecisionBuffer<double>::binary64(DenseReal{}, generation);
    REQUIRE(empty);
    const auto zero = original_residual(DenseReal{}, *empty, DenseReal{}, plan);
    REQUIRE(zero);
    CHECK(zero->values().empty());
    const auto bad =
        original_residual(matrix<double>({1, 1, 1}, {1}, plan), *empty, DenseReal{}, plan);
    REQUIRE_FALSE(bad);
    CHECK(bad.error().code == core::FailureCode::InvalidInput);
}
TEST_CASE("Projection qualifies only its actual residual observation and consumer",
          "[verification][negative][N4-06]") {
    core::ResourcePlan plan{1024ULL * 1024};
    auto x = PrecisionBuffer<double>::two_word(matrix<double>({1, 1, 1}, {1}, plan),
                                               matrix<double>({1, 1, 1}, {0}, plan), generation);
    REQUIRE(x);
    const auto a = matrix<double>({1, 1, 1}, {2}, plan), b = matrix<double>({1, 1, 1}, {2}, plan);
    const std::array scales{1.0};
    Verifier<double> verifier;
    auto accepted = project(*x, a, b, scales, 0, verifier, plan);
    REQUIRE(accepted);
    CHECK(verifier.calls == 1);
    CHECK(accepted->context().source == generation);
    CHECK(accepted->values().values()[0] == 1);
    verifier.unavailable = true;
    auto unavailable = project(*x, a, b, scales, 0, verifier, plan);
    REQUIRE_FALSE(unavailable);
    CHECK(unavailable.error().code == core::FailureCode::PrecisionUnresolved);
    verifier.unavailable = false;
    verifier.unsupported = true;
    auto unsupported = project(*x, a, b, scales, 0, verifier, plan);
    REQUIRE_FALSE(unsupported);
    CHECK(unsupported.error().code == core::FailureCode::PrecisionUnsupported);
    verifier.unsupported = false;
    verifier.throwing = true;
    auto thrown = project(*x, a, b, scales, 0, verifier, plan);
    REQUIRE_FALSE(thrown);
    CHECK(thrown.error().code == core::FailureCode::PrecisionUnresolved);
    verifier.throwing = false;
    core::ResourcePlan none{0};
    auto limited = project(*x, a, b, scales, 0, verifier, none);
    REQUIRE_FALSE(limited);
    CHECK(limited.error().code == core::FailureCode::ResourceExceeded);
    CHECK(std::get<TwoWord<double>>(x->storage()).hi.values()[0] == 1);
    auto sensitive =
        PrecisionBuffer<double>::two_word(matrix<double>({1, 1, 1}, {1}, plan),
                                          matrix<double>({1, 1, 1}, {0x1p-54}, plan), generation);
    REQUIRE(sensitive);
    const auto zero = matrix<double>({1, 1, 1}, {0}, plan);
    verifier.require_exact_first = true;
    auto lost = project(*sensitive, zero, zero, scales, 0, verifier, plan);
    REQUIRE_FALSE(lost);
    CHECK(lost.error().code == core::FailureCode::PrecisionUnresolved);
}
TEST_CASE("Projection conservatively rejects compensation loss and underflowed normalization",
          "[negative][verification][N4-06]") {
    core::ResourcePlan plan{1024ULL * 1024};
    const Verifier<double> verifier;
    const std::array scales{1.0};
    const auto a = matrix<double>({1, 5, 1}, {-1, -1e-16, -1e-300, 1, 1e-16}, plan);
    auto x = PrecisionBuffer<double>::binary64(matrix<double>({5, 1, 5}, {1, 1, 1, 1, 1}, plan),
                                               generation);
    REQUIRE(x);
    const auto b = matrix<double>({1, 1, 1}, {0}, plan);
    const auto diagnostic = original_residual(a, *x, b, plan);
    REQUIRE(diagnostic);
    CHECK(diagnostic->values()[0] == 0);
    const auto rejected = project(*x, a, b, scales, 0, verifier, plan);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::PrecisionUnresolved);
    CHECK(verifier.calls == 0);
    // Independent oracle resolves cancellation in two groups; binary100 alone cannot retain
    // 1e-300 beside 1. The exact cancelling pairs leave the unchanged third term.
    const Big exact = (Big(1) - Big(1)) + (Big(1e-16) - Big(1e-16)) + Big(1e-300);
    CHECK(exact > 0);
    auto singleton =
        PrecisionBuffer<double>::binary64(matrix<double>({1, 1, 1}, {0}, plan), generation);
    REQUIRE(singleton);
    const auto identity = matrix<double>({1, 1, 1}, {1}, plan);
    const auto tiny = matrix<double>({1, 1, 1}, {std::numeric_limits<double>::denorm_min()}, plan);
    const std::array huge_scale{std::numeric_limits<double>::max()};
    const auto tiny_rejected = project(*singleton, identity, tiny, huge_scale, 0, verifier, plan);
    REQUIRE_FALSE(tiny_rejected);
    CHECK(tiny_rejected.error().code == core::FailureCode::PrecisionUnresolved);
    const auto accepted = project(*singleton, identity, tiny, huge_scale,
                                  std::numeric_limits<double>::denorm_min(), verifier, plan);
    REQUIRE(accepted);
}
TEST_CASE("Projection threshold comparisons preserve exact boundaries across extreme scales",
          "[verification][negative][N4-06]") {
    core::ResourcePlan plan{1024ULL * 1024};
    const Verifier<double> verifier;
    auto x = PrecisionBuffer<double>::binary64(matrix<double>({1, 1, 1}, {0}, plan), generation);
    REQUIRE(x);
    const auto a = matrix<double>({1, 1, 1}, {0}, plan);
    const auto b = matrix<double>({1, 1, 1}, {1}, plan);
    const std::array scales{1.0};
    const auto equal = project(*x, a, b, scales, 1, verifier, plan);
    REQUIRE(equal);
    const auto below = project(*x, a, b, scales, std::nextafter(1.0, 0.0), verifier, plan);
    REQUIRE_FALSE(below);
    CHECK(below.error().code == core::FailureCode::PrecisionUnresolved);
    const auto above = project(*x, a, b, scales, std::nextafter(1.0, 2.0), verifier, plan);
    REQUIRE(above);
    for (const int exponent : {-900, 900}) {
        const double unit = std::ldexp(1.0, exponent);
        const auto rhs = matrix<double>({1, 1, 1}, {unit}, plan);
        const std::array units{unit};
        CHECK(project(*x, a, rhs, units, 1, verifier, plan));
        CHECK_FALSE(project(*x, a, rhs, units, std::nextafter(1.0, 0.0), verifier, plan));
    }
    // The exact product is below its rounded high word; comparing high words alone is unsafe.
    const double t = 1.0 + 0x1p-27, scale = 1.0 - 0x1p-27;
    CHECK(Big(t) * Big(scale) < Big(1));
    const std::array tricky{scale};
    CHECK_FALSE(project(*x, a, b, tricky, t, verifier, plan));
}
TEST_CASE("Precision transformations ignore finite physical padding and reject true overflow",
          "[negative][N4-02][N4-06]") {
    core::ResourcePlan plan{1024ULL * 1024};
    const double largest = std::numeric_limits<double>::max();
    auto x = PrecisionBuffer<double>::two_word(matrix<double>({1, 1, 2}, {1, largest}, plan),
                                               matrix<double>({1, 1, 2}, {0, largest}, plan),
                                               generation);
    REQUIRE(x);
    const auto correction = matrix<double>({1, 1, 2}, {0, largest}, plan);
    auto updated = add_correction(*x, correction, core::NumericalGeneration{18}, plan);
    REQUIRE(updated);
    const auto& pair = std::get<TwoWord<double>>(updated->storage());
    CHECK(pair.hi.values()[1] == largest);
    CHECK(pair.lo.values()[1] == largest);
    const auto a = matrix<double>({1, 1, 1}, {1}, plan), b = matrix<double>({1, 1, 1}, {1}, plan);
    const std::array scales{1.0};
    const Verifier<double> verifier;
    auto projected = project(*x, a, b, scales, 0, verifier, plan);
    REQUIRE(projected);
    CHECK(projected->values().values()[1] == largest);
    auto overflow =
        PrecisionBuffer<double>::two_word(matrix<double>({1, 1, 1}, {largest}, plan),
                                          matrix<double>({1, 1, 1}, {largest}, plan), generation);
    REQUIRE(overflow);
    const auto refused = project(*overflow, a, b, scales, 1, verifier, plan);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().code == core::FailureCode::Nonfinite);
    auto complex_source =
        PrecisionBuffer<C>::binary64(matrix<C>({1, 1, 1}, {C{1, 2}}, plan), generation);
    REQUIRE(complex_source);
    const auto delta = matrix<C>({1, 1, 1}, {C{0x1p-54, 0x1p-53}}, plan);
    const auto promoted =
        add_correction(*complex_source, delta, core::NumericalGeneration{18}, plan);
    REQUIRE(promoted);
    const auto& cp = std::get<TwoWord<C>>(promoted->storage());
    CHECK(Big(cp.hi.values()[0].real()) + Big(cp.lo.values()[0].real()) ==
          Big(1) + pow(Big(2), -54));
    CHECK(Big(cp.hi.values()[0].imag()) + Big(cp.lo.values()[0].imag()) ==
          Big(2) + pow(Big(2), -53));
}
TEST_CASE("Projection checks invalid scales generations shapes and empty consumers",
          "[negative][N4-06]") {
    core::ResourcePlan plan{1024ULL * 1024};
    Verifier<double> verifier;
    auto x = PrecisionBuffer<double>::binary64(matrix<double>({1, 1, 1}, {1}, plan), generation);
    REQUIRE(x);
    const auto a = matrix<double>({1, 1, 1}, {1}, plan), b = matrix<double>({1, 1, 1}, {1}, plan);
    for (const double scale : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN()}) {
        const std::array scales{scale};
        const auto result = project(*x, a, b, scales, 1e-12, verifier, plan);
        REQUIRE_FALSE(result);
        CHECK(result.error().code == (std::isfinite(scale) ? core::FailureCode::InvalidInput
                                                           : core::FailureCode::Nonfinite));
    }
    const std::array scale{1.0};
    CHECK_FALSE(project(*x, a, b, {}, 1e-12, verifier, plan));
    CHECK_FALSE(project(*x, a, b, scale, -1, verifier, plan));
    CHECK_FALSE(
        project_binary64(*x, a, b, scale, 1e-12, generation, core::ScopeId{31}, verifier, plan));
    auto empty = PrecisionBuffer<double>::binary64(DenseReal{}, generation);
    REQUIRE(empty);
    auto accepted = project(*empty, DenseReal{}, DenseReal{}, {}, 0, verifier, plan);
    REQUIRE(accepted);
    CHECK(verifier.calls == 1);
    verifier.unavailable = true;
    CHECK_FALSE(project(*empty, DenseReal{}, DenseReal{}, {}, 0, verifier, plan));
}
TEST_CASE("Actual 4096-grid projection witness is rejected by the product residual",
          "[verification][negative][N4-06]") {
    std::ifstream input(UVLAS_PRECISION_WITNESS_PATH);
    REQUIRE(input.is_open());
    const auto fixture = nlohmann::json::parse(input);
    REQUIRE(fixture.at("m").get<int>() == 4096);
    REQUIRE(fixture.at("row").get<int>() == 4095);
    REQUIRE(fixture.at("column").get<int>() == 1);
    REQUIRE(fixture.at("input_hashes_equal").get<bool>());
    core::ResourcePlan plan{8 * 1024ULL * 1024};
    const auto read = [&](std::string_view key, DenseShape shape) {
        std::vector<C> values;
        for (const auto& pair : fixture.at(key)) {
            values.emplace_back(pair.at(0).get<double>(), pair.at(1).get<double>());
        }
        auto result = DenseComplex::copy(shape, values, plan, op);
        REQUIRE(result);
        return std::move(*result);
    };
    const auto a = read("a", {1, 4096, 1}), b = read("b", {1, 1, 1});
    auto hi = read("hi", {4096, 1, 4096}), lo = read("lo", {4096, 1, 4096});
    Big retained_real = b.values()[0].real(), retained_imag = b.values()[0].imag();
    Big projected_real = retained_real, projected_imag = retained_imag;
    for (std::size_t i = 0; i < 4096; ++i) {
        const Big ar = a.values()[i].real(), ai = a.values()[i].imag();
        const Big xr = Big(hi.values()[i].real()) + Big(lo.values()[i].real());
        const Big xi = Big(hi.values()[i].imag()) + Big(lo.values()[i].imag());
        retained_real -= ar * xr - ai * xi;
        retained_imag -= ar * xi + ai * xr;
        const C rounded = hi.values()[i] + lo.values()[i];
        projected_real -= ar * Big(rounded.real()) - ai * Big(rounded.imag());
        projected_imag -= ar * Big(rounded.imag()) + ai * Big(rounded.real());
    }
    const double scale = fixture.at("scale").get<double>();
    // Frozen R&D criterion: ||r||2/||b||2 <= 1e-11. One row above that same
    // denominator proves rejection, without claiming a product full-system acceptance.
    CHECK(sqrt(projected_real * projected_real + projected_imag * projected_imag) / Big(scale) >
          Big("1e-11"));
    CHECK(sqrt(retained_real * retained_real + retained_imag * retained_imag) / Big(scale) <
          Big("1e-18"));
    auto source = PrecisionBuffer<C>::two_word(std::move(hi), std::move(lo), generation);
    REQUIRE(source);
    const auto residual = original_residual(a, *source, b, plan);
    REQUIRE(residual);
    // Decimal-selected-row agreement fixed at 1e-18 in the original experiment protocol.
    CHECK(abs(Big(residual->values()[0].real()) - retained_real) / Big(scale) < Big("1e-18"));
    CHECK(abs(Big(residual->values()[0].imag()) - retained_imag) / Big(scale) < Big("1e-18"));
    const Verifier<C> verifier;
    const std::array scales{scale};
    const auto projection = project(*source, a, b, scales, 1e-11, verifier, plan);
    REQUIRE_FALSE(projection);
    CHECK(projection.error().code == core::FailureCode::PrecisionUnresolved);
    CHECK(verifier.calls == 0);
}
namespace {
template <class T>
concept TemporaryProjectionBorrow = requires(T&& value) { std::move(value).values(); };
template <class T>
concept TemporaryPrecisionBorrow = requires(T&& value) { std::move(value).storage(); };
static_assert(!TemporaryProjectionBorrow<QualifiedProjection<double>>);
static_assert(!TemporaryPrecisionBorrow<PrecisionBuffer<double>>);
static_assert(!std::is_default_constructible_v<QualifiedProjection<double>>);
static_assert(!std::is_constructible_v<QualifiedProjection<double>, DenseReal, ProjectionContext>);
static_assert(!std::is_constructible_v<QualifiedProjection<C>, DenseComplex, ProjectionContext>);
static_assert(!std::is_constructible_v<PrecisionBuffer<double>,
                                       std::variant<Binary64<double>, TwoWord<double>>,
                                       core::NumericalGeneration>);
static_assert(!std::is_convertible_v<PrecisionBuffer<double>, double>);
static_assert(std::is_same_v<decltype(std::declval<const QualifiedProjection<double>&>().values()),
                             const DenseReal&>);
} // namespace
#if defined(_M_X64) || defined(__x86_64__)
TEST_CASE("EFT refuses flush-to-zero and denormal-as-zero modes", "[negative][N4-06]") {
    const unsigned previous = _mm_getcsr();
    for (const unsigned unsupported : {0x8000U, 0x40U}) {
        _mm_setcsr(previous | unsupported);
        const auto sum = two_sum(1, 1), product = two_product(1, 1);
        _mm_setcsr(previous);
        REQUIRE_FALSE(sum);
        CHECK(sum.error().code == core::FailureCode::PrecisionUnsupported);
        REQUIRE_FALSE(product);
        CHECK(product.error().code == core::FailureCode::PrecisionUnsupported);
    }
}
#endif

TEST_CASE("Projection rejects malformed residual inputs before allocating", "[negative][N4-06]") {
    core::ResourcePlan inputs{1024ULL * 1024};
    auto source =
        PrecisionBuffer<double>::binary64(matrix<double>({1, 1, 1}, {1}, inputs), generation);
    REQUIRE(source);
    const auto wrong_shape = matrix<double>({2, 1, 2}, {1, 1}, inputs);
    const auto identity = matrix<double>({1, 1, 1}, {1}, inputs);
    auto a = matrix<double>({1, 1, 1}, {1}, inputs);
    auto b = matrix<double>({1, 1, 1}, {1}, inputs);
    const Verifier<double> verifier;
    const std::array scales{1.0};
    core::ResourcePlan zero{0};
    const auto shape = project(*source, wrong_shape, b, scales, 0, verifier, zero);
    REQUIRE_FALSE(shape);
    CHECK(shape.error().code == core::FailureCode::InvalidInput);
    CHECK(zero.allocated_bytes() == 0);
    a.values()[0] = std::numeric_limits<double>::quiet_NaN();
    const auto bad_a = project(*source, a, b, scales, 0, verifier, zero);
    REQUIRE_FALSE(bad_a);
    CHECK(bad_a.error().code == core::FailureCode::Nonfinite);
    CHECK(zero.allocated_bytes() == 0);
    b.values()[0] = std::numeric_limits<double>::infinity();
    const auto bad_b = project(*source, identity, b, scales, 0, verifier, zero);
    REQUIRE_FALSE(bad_b);
    CHECK(bad_b.error().code == core::FailureCode::Nonfinite);
    CHECK(zero.allocated_bytes() == 0);
    CHECK(verifier.calls == 0);
}
