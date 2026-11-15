#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <uvlas/numerics/dense_solve.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;
using Big = boost::multiprecision::cpp_bin_float_100;
using Complex = std::complex<double>;
constexpr auto operation = core::Operation::SolveLinear;
Context context() {
    return {{core::ModelFingerprint{}, core::ModelGeneration{1}, core::GeometryGeneration{2},
             core::ChartGeneration{3}, core::SparsityGeneration{4}, core::NumericalGeneration{5},
             core::ScaleGeneration{6}, core::PrecisionGeneration{7}, core::Precision::Binary64,
             core::WakePolicyGeneration{8}, core::TimeIntervalGeneration{9},
             core::BackendGeneration{10}},
            core::TrialId{11},
            core::ScopeId{12},
            core::ScopeId{13},
            0.25};
}
template <StorageScalar T>
class NoChecks final : public ConsumerEvaluator<T> {
public:
    [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
        return SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<T>&, const PrecisionBuffer<T>&, const Context&,
           const ConsumerCheck& check, core::CertificateContext c) const override {
        FAIL_CHECK("algebraic solve must not call an absent physical requirement");
        return core::ErrorCertificate::unevaluated(c, check.budget);
    }
};

template <StorageScalar T, std::size_t N>
DenseMatrix<T> matrix(core::Index rows, core::Index columns, const std::array<T, N>& values,
                      core::ResourcePlan& plan) {
    auto value = DenseMatrix<T>::copy({rows, columns, rows}, values, plan, operation);
    REQUIRE(value);
    return std::move(*value);
}
template <StorageScalar T>
const DenseMatrix<T>& high(const PrecisionBuffer<T>& values) {
    if (auto* b = std::get_if<Binary64<T>>(&values.storage()))
        return b->values;
    return std::get<TwoWord<T>>(values.storage()).hi;
}
template <StorageScalar T>
T low(const PrecisionBuffer<T>& values, std::size_t i) {
    if (auto* b = std::get_if<TwoWord<T>>(&values.storage()))
        return b->lo.values()[i];
    return {};
}
struct Setup {
    core::ResourcePlan plan{1 << 24};
    std::array<double, 3> s{1, 2, 4};
    std::array<UnitDimensions, 3> units{};
    [[nodiscard]] Scales scales(std::size_t n) const {
        return {std::span{s}.first(n), std::span{s}.first(n), std::span{units}.first(n),
                std::span{units}.first(n), context()};
    }
    [[nodiscard]] Limits limits(core::Index n, core::Index k = 1, core::Index r = 0) const {
        return {n, k, r, 5, EstimateAdmission{1 << 24}};
    }
    [[nodiscard]] Accuracy accuracy(double tolerance = 1e-12) const {
        return {Norm::Infinity,         tolerance, 0, tolerance, core::ObservationId{20},
                core::ObservationId{21}};
    }
};
// Direct fixed-size elimination in a separate 100-decimal-digit arithmetic domain.
// Row-major literal A and b are independent of Eigen, production LU and residual evaluation.
template <class Wide>
std::array<Wide, 3> reference(std::array<std::array<Wide, 3>, 3> a, std::array<Wide, 3> b) {
    for (std::size_t k = 0; k < 3; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < 3; ++i)
            if (abs(a[i][k]) > abs(a[pivot][k]))
                pivot = i;
        std::swap(a[k], a[pivot]);
        std::swap(b[k], b[pivot]);
        for (std::size_t i = k + 1; i < 3; ++i) {
            Wide ratio = a[i][k] / a[k][k];
            for (std::size_t j = k; j < 3; ++j)
                a[i][j] -= ratio * a[k][j];
            b[i] -= ratio * b[k];
        }
    }
    std::array<Wide, 3> x{};
    for (std::size_t end = 3; end > 0;) {
        auto i = --end;
        Wide value = b[i];
        for (std::size_t j = i + 1; j < 3; ++j)
            value -= a[i][j] * x[j];
        x[i] = value / a[i][i];
    }
    return x;
}
} // namespace
TEST_CASE("Dense pivoted LU agrees with an independent high precision original solve",
          "[verification][N4-05][N4-06]") {
    Setup f;
    auto a = matrix(3, 3, std::array{0., 1., 3., 2., -2., -1., 1., -3., 2.}, f.plan);
    auto b = matrix(3, 1, std::array{-1., -4., 11.}, f.plan);
    auto rhs = PrecisionBuffer<double>::binary64(std::move(b), core::NumericalGeneration{30});
    REQUIRE(rhs);
    const NoChecks<double> evaluator;
    auto scales = f.scales(3);
    const DenseOperator op{a, context()};
    auto prepared = DensePrepared::prepare(op, scales, f.limits(3), f.plan);
    REQUIRE(prepared);
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{31}}, {}, evaluator};
    auto solved = solve_linear(op, Rhs{*rhs, context(), core::RhsGeneration{40}}, scales, precision,
                               f.accuracy(), f.limits(3), *prepared, f.plan);
    REQUIRE(solved);
    auto expected = reference<Big>({{{0, 2, 1}, {1, -2, -3}, {3, -1, 2}}}, {-1, -4, 11});
    for (std::size_t i = 0; i < 3; ++i)
        CHECK(abs(Big(high(solved->values()).values()[i]) + Big(low(solved->values(), i)) -
                  expected[i]) < Big("1e-13"));
    CHECK(std::holds_alternative<ReciprocalOneNormEstimate>(solved->certificate().condition()));
}
TEST_CASE("Dense stagnation promotes and preserves the correction low word",
          "[verification][negative][N4-05][N4-06]") {
    Setup f;
    auto a = matrix(1, 1, std::array{3.}, f.plan);
    auto b = matrix(1, 1, std::array{1.}, f.plan);
    auto rhs = PrecisionBuffer<double>::binary64(std::move(b), core::NumericalGeneration{30});
    REQUIRE(rhs);
    auto scales = f.scales(1);
    const NoChecks<double> evaluator;
    const DenseOperator op{a, context()};
    auto prepared = DensePrepared::prepare(op, scales, f.limits(1), f.plan);
    REQUIRE(prepared);
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{31}}, {}, evaluator};
    auto solved = solve_linear(op, Rhs{*rhs, context(), core::RhsGeneration{40}}, scales, precision,
                               f.accuracy(1e-30), f.limits(1), *prepared, f.plan);
    REQUIRE(solved);
    CHECK(solved->values().precision() == core::Precision::TwoWord);
    CHECK(solved->certificate().retention_path() == RetentionPath::PromotedTwoWord);
    CHECK(solved->certificate().refinements() <= 5);
    CHECK(abs(Big(high(solved->values()).values()[0]) + Big(low(solved->values(), 0)) -
              Big(1) / 3) < Big("1e-31"));
    const Precision<double> projection{
        Binary64Projection{core::NumericalGeneration{31}, core::NumericalGeneration{32}},
        {},
        evaluator};
    auto rejected = solve_linear(op, Rhs{*rhs, context(), core::RhsGeneration{40}}, scales,
                                 projection, f.accuracy(1e-30), f.limits(1), *prepared, f.plan);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::PrecisionUnresolved);
}
TEST_CASE("A tiny nonzero pivot is not a rank or condition threshold",
          "[verification][negative][N4-05][N4-06]") {
    Setup f;
    auto a = matrix(2, 2, std::array{1., 0., 0., 1e-300}, f.plan);
    auto b = matrix(2, 1, std::array{1., 1e-300}, f.plan);
    auto rhs = PrecisionBuffer<double>::binary64(std::move(b), core::NumericalGeneration{30});
    REQUIRE(rhs);
    auto scales = f.scales(2);
    const NoChecks<double> evaluator;
    const DenseOperator op{a, context()};
    auto prepared = DensePrepared::prepare(op, scales, f.limits(2), f.plan);
    REQUIRE(prepared);
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{31}}, {}, evaluator};
    auto solved = solve_linear(op, Rhs{*rhs, context(), core::RhsGeneration{40}}, scales, precision,
                               f.accuracy(), f.limits(2), *prepared, f.plan);
    REQUIRE(solved);
    CHECK(high(solved->values()).values()[0] == 1);
    CHECK(high(solved->values()).values()[1] == 1);
}
TEST_CASE("Complex rank-two solve uses transpose without conjugation on both RHS columns",
          "[verification][N4-05][N4-06]") {
    Setup f;
    auto a = matrix(3, 3, std::array{0., 1., 3., 2., -2., -1., 1., -3., 2.}, f.plan);
    auto u = matrix(3, 2,
                    std::array{Complex{0, 1. / 16}, Complex{1. / 16, -2. / 16},
                               Complex{3. / 16, 1. / 16}, Complex{1. / 16, 1. / 16},
                               Complex{2. / 16, 0}, Complex{0, -1. / 16}},
                    f.plan);
    auto v = matrix(3, 2,
                    std::array{Complex{1, 1}, Complex{0, 0}, Complex{.5, 0}, Complex{0, 0},
                               Complex{1, -1}, Complex{-1, 0}},
                    f.plan);
    auto b = matrix(3, 2,
                    std::array{Complex{-81. / 32, 43. / 8}, Complex{-7. / 8, -41. / 32},
                               Complex{303. / 32, -51. / 32}, Complex{57. / 8, -127. / 32},
                               Complex{-321. / 32, 115. / 8}, Complex{-139. / 32, -79. / 32}},
                    f.plan);
    auto rhs = PrecisionBuffer<Complex>::binary64(std::move(b), core::NumericalGeneration{30});
    REQUIRE(rhs);
    const NoChecks<Complex> evaluator;
    auto scales = f.scales(3);
    const ColumnUpdatedOperator op{a, u, v, context()};
    auto limits = f.limits(3, 2, 2);
    auto prepared = UpdatedPrepared::prepare(op, scales, limits, f.plan);
    REQUIRE(prepared);
    const Precision<Complex> precision{
        RetainedPrecision{core::NumericalGeneration{31}}, {}, evaluator};
    auto solved = solve_linear(op, Rhs{*rhs, context(), core::RhsGeneration{40}}, scales, precision,
                               f.accuracy(), limits, *prepared, f.plan);
    REQUIRE(solved);
    const std::array expected{Complex{1, 1},  Complex{-2, 3}, Complex{2, -1},
                              Complex{-1, 2}, Complex{3, 0},  Complex{1, -4}};
    for (std::size_t i = 0; i < expected.size(); ++i)
        CHECK(std::abs(high(solved->values()).values()[i] + low(solved->values(), i) -
                       expected[i]) < 1e-12);
    using Wide = boost::multiprecision::cpp_complex_100;
    // Literal effective A and RHS are dyadic rationals, independent of the low-rank action.
    const std::array<std::array<Wide, 3>, 3> effective{
        {{Wide{Big(-1) / 16, Big(1) / 16}, Wide{Big(17) / 8, Big(0)},
          Wide{Big(30) / 32, Big(-1) / 32}},
         {Wide{Big(19) / 16, Big(-1) / 16}, Wide{Big(-15) / 8, Big(-1) / 8},
          Wide{Big(-99) / 32, Big(-2) / 32}},
         {Wide{Big(25) / 8, Big(2) / 8}, Wide{Big(-17) / 16, Big(-1) / 16},
          Wide{Big(67) / 32, Big(3) / 32}}}};
    const std::array<std::array<Wide, 3>, 2> loads{
        {{Wide{Big(-81) / 32, Big(43) / 8}, Wide{Big(-7) / 8, Big(-41) / 32},
          Wide{Big(303) / 32, Big(-51) / 32}},
         {Wide{Big(57) / 8, Big(-127) / 32}, Wide{Big(-321) / 32, Big(115) / 8},
          Wide{Big(-139) / 32, Big(-79) / 32}}}};
    for (std::size_t c = 0; c < 2; ++c) {
        const auto direct = reference<Wide>(effective, loads[c]);
        for (std::size_t i = 0; i < 3; ++i) {
            const auto actual =
                high(solved->values()).values()[c * 3 + i] + low(solved->values(), c * 3 + i);
            CHECK(abs(Wide{Big(actual.real()), Big(actual.imag())} - direct[i]) < Big("1e-12"));
            CHECK(abs(direct[i] - Wide{Big(expected[c * 3 + i].real()),
                                       Big(expected[c * 3 + i].imag())}) < Big("1e-90"));
        }
    }
    CHECK(std::holds_alternative<NotEvaluated>(solved->certificate().condition()));
}
TEST_CASE("Dense singularity, zero loads and empty order have distinct outcomes",
          "[negative][verification][N4-05][N4-06]") {
    Setup f;
    const NoChecks<double> evaluator;
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{31}}, {}, evaluator};
    auto singular = matrix(2, 2, std::array{1., 2., 2., 4.}, f.plan);
    const DenseOperator bad{singular, context()};
    auto factor = DensePrepared::prepare(bad, f.scales(2), f.limits(2), f.plan);
    REQUIRE_FALSE(factor);
    CHECK(factor.error().code == core::FailureCode::RankUnresolved);
    auto a = matrix(2, 2, std::array{2., 1., 0., 3.}, f.plan);
    auto zero = matrix(2, 1, std::array{0., 0.}, f.plan);
    auto rhs = PrecisionBuffer<double>::binary64(std::move(zero), core::NumericalGeneration{30});
    REQUIRE(rhs);
    const DenseOperator op{a, context()};
    auto prepared = DensePrepared::prepare(op, f.scales(2), f.limits(2), f.plan);
    REQUIRE(prepared);
    auto solved = solve_linear(op, Rhs{*rhs, context(), core::RhsGeneration{40}}, f.scales(2),
                               precision, f.accuracy(0), f.limits(2), *prepared, f.plan);
    REQUIRE(solved);
    for (auto value : high(solved->values()).values())
        CHECK(value == 0);
    auto empty = matrix(0, 0, std::array<double, 0>{}, f.plan);
    auto empty_rhs = matrix(0, 2, std::array<double, 0>{}, f.plan);
    auto buffer =
        PrecisionBuffer<double>::binary64(std::move(empty_rhs), core::NumericalGeneration{30});
    REQUIRE(buffer);
    const DenseOperator empty_op{empty, context()};
    auto empty_limits = Limits{0, 2, 0, 5, RequireProvenBound{0}};
    auto empty_prepared = DensePrepared::prepare(empty_op, f.scales(0), empty_limits, f.plan);
    REQUIRE(empty_prepared);
    auto empty_solution =
        solve_linear(empty_op, Rhs{*buffer, context(), core::RhsGeneration{40}}, f.scales(0),
                     precision, f.accuracy(0), empty_limits, *empty_prepared, f.plan);
    REQUIRE(empty_solution);
    CHECK(empty_solution->values().shape().rows == 0);
    CHECK(empty_solution->values().shape().columns == 2);
}
TEST_CASE("Singular complex Schur is rejected while rank-zero uses the real base",
          "[negative][verification][N4-05][N4-06]") {
    Setup f;
    auto a = matrix(1, 1, std::array{1.}, f.plan);
    auto u = matrix(1, 1, std::array{Complex{-1, 0}}, f.plan);
    auto v = matrix(1, 1, std::array{Complex{1, 0}}, f.plan);
    auto bad = UpdatedPrepared::prepare(ColumnUpdatedOperator{a, u, v, context()}, f.scales(1),
                                        f.limits(1, 1, 1), f.plan);
    REQUIRE_FALSE(bad);
    CHECK(bad.error().code == core::FailureCode::RankUnresolved);
    auto zero_u = matrix(1, 0, std::array<Complex, 0>{}, f.plan);
    auto zero_v = matrix(1, 0, std::array<Complex, 0>{}, f.plan);
    auto b = matrix(1, 1, std::array{Complex{2, 3}}, f.plan);
    auto rhs = PrecisionBuffer<Complex>::binary64(std::move(b), core::NumericalGeneration{30});
    REQUIRE(rhs);
    const NoChecks<Complex> evaluator;
    const Precision<Complex> precision{
        RetainedPrecision{core::NumericalGeneration{31}}, {}, evaluator};
    const ColumnUpdatedOperator op{a, zero_u, zero_v, context()};
    auto prepared = UpdatedPrepared::prepare(op, f.scales(1), f.limits(1), f.plan);
    REQUIRE(prepared);
    auto result = solve_linear(op, Rhs{*rhs, context(), core::RhsGeneration{40}}, f.scales(1),
                               precision, f.accuracy(), f.limits(1), *prepared, f.plan);
    REQUIRE(result);
    CHECK(high(result->values()).values()[0] == Complex{2, 3});
}
TEST_CASE("An initially retained low RHS word is preserved in real and complex solves",
          "[verification][N4-05][N4-06]") {
    Setup f;
    auto a = matrix(1, 1, std::array{1.}, f.plan);
    auto hi = matrix(1, 1, std::array{1.}, f.plan);
    auto lo = matrix(1, 1, std::array{std::ldexp(1., -60)}, f.plan);
    auto rhs = PrecisionBuffer<double>::two_word(std::move(hi), std::move(lo),
                                                 core::NumericalGeneration{30});
    REQUIRE(rhs);
    const NoChecks<double> evaluator;
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{31}}, {}, evaluator};
    const DenseOperator op{a, context()};
    auto prepared = DensePrepared::prepare(op, f.scales(1), f.limits(1), f.plan);
    REQUIRE(prepared);
    auto result = solve_linear(op, Rhs{*rhs, context(), core::RhsGeneration{40}}, f.scales(1),
                               precision, f.accuracy(1e-30), f.limits(1), *prepared, f.plan);
    REQUIRE(result);
    CHECK(low(result->values(), 0) == std::ldexp(1., -60));
    CHECK(result->certificate().retention_path() == RetentionPath::InitiallyTwoWord);
    auto zero_u = matrix(1, 0, std::array<Complex, 0>{}, f.plan),
         zero_v = matrix(1, 0, std::array<Complex, 0>{}, f.plan);
    auto ch = matrix(1, 1, std::array{Complex{1, 2}}, f.plan),
         cl = matrix(1, 1, std::array{Complex{std::ldexp(1., -60), -std::ldexp(1., -59)}}, f.plan);
    auto cb = PrecisionBuffer<Complex>::two_word(std::move(ch), std::move(cl),
                                                 core::NumericalGeneration{30});
    REQUIRE(cb);
    const NoChecks<Complex> ce;
    const Precision<Complex> cp{RetainedPrecision{core::NumericalGeneration{31}}, {}, ce};
    const ColumnUpdatedOperator cop{a, zero_u, zero_v, context()};
    auto prepared_c = UpdatedPrepared::prepare(cop, f.scales(1), f.limits(1), f.plan);
    REQUIRE(prepared_c);
    auto cr = solve_linear(cop, Rhs{*cb, context(), core::RhsGeneration{40}}, f.scales(1), cp,
                           f.accuracy(1e-30), f.limits(1), *prepared_c, f.plan);
    REQUIRE(cr);
    CHECK(low(cr->values(), 0) == Complex{std::ldexp(1., -60), -std::ldexp(1., -59)});
}
TEST_CASE("Dense input validation precedes a zero allocation budget and retains parents",
          "[negative][N4-05][N4-06]") {
    Setup f;
    auto a = matrix(1, 1, std::array{1.}, f.plan);
    core::ResourcePlan none{0};
    a.values()[0] = std::numeric_limits<double>::quiet_NaN();
    auto invalid =
        DensePrepared::prepare(DenseOperator{a, context()}, f.scales(1), f.limits(1), none);
    REQUIRE_FALSE(invalid);
    CHECK(invalid.error().code == core::FailureCode::Nonfinite);
    CHECK(none.allocated_bytes() == 0);
    CHECK(std::isnan(a.values()[0]));
    a.values()[0] = 1;
    auto u = matrix(1, 1, std::array{Complex{1, 2}}, f.plan),
         v = matrix(1, 1, std::array{Complex{3, 4}}, f.plan);
    for (auto* source : {&u, &v})
        for (const bool imaginary : {false, true}) {
            const auto saved = source->values()[0];
            source->values()[0] =
                imaginary ? Complex{saved.real(), std::numeric_limits<double>::infinity()}
                          : Complex{std::numeric_limits<double>::infinity(), saved.imag()};
            auto rejected = UpdatedPrepared::prepare(ColumnUpdatedOperator{a, u, v, context()},
                                                     f.scales(1), f.limits(1, 1, 1), none);
            REQUIRE_FALSE(rejected);
            CHECK(rejected.error().code == core::FailureCode::Nonfinite);
            CHECK(none.allocated_bytes() == 0);
            source->values()[0] = saved;
        }
    auto rectangular = matrix(1, 2, std::array{1., 2.}, f.plan);
    auto bad_shape = DensePrepared::prepare(DenseOperator{rectangular, context()}, f.scales(1),
                                            f.limits(1), none);
    REQUIRE_FALSE(bad_shape);
    CHECK(bad_shape.error().code == core::FailureCode::InvalidInput);
    CHECK(none.allocated_bytes() == 0);
}
TEST_CASE("Dense mathematical rows exclude finite storage padding",
          "[verification][N4-02][N4-05][N4-06]") {
    Setup f;
    const std::array ac{3., 1e300}, bc{1., -1e300};
    auto a = DenseReal::copy({1, 1, 2}, ac, f.plan, operation);
    auto b = DenseReal::copy({1, 1, 2}, bc, f.plan, operation);
    REQUIRE(a);
    REQUIRE(b);
    auto rhs = PrecisionBuffer<double>::binary64(std::move(*b), core::NumericalGeneration{30});
    REQUIRE(rhs);
    const NoChecks<double> evaluator;
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{31}}, {}, evaluator};
    const DenseOperator op{*a, context()};
    auto prepared = DensePrepared::prepare(op, f.scales(1), f.limits(1), f.plan);
    REQUIRE(prepared);
    auto result = solve_linear(op, Rhs{*rhs, context(), core::RhsGeneration{40}}, f.scales(1),
                               precision, f.accuracy(1e-30), f.limits(1), *prepared, f.plan);
    REQUIRE(result);
    CHECK(high(result->values()).shape().leading_dimension == 2);
    CHECK(high(result->values()).values()[1] == -1e300);
    CHECK(low(result->values(), 1) == 0);
    CHECK(a->values()[1] == 1e300);
}
TEST_CASE("Original residual stagnation promotes even while a coordinate keeps moving",
          "[verification][negative][N4-05][N4-06]") {
    Setup f;
    auto a = matrix(3, 3, std::array{7., -2., 6., 5., -5., 5., 8., 8., 6.}, f.plan);
    auto b = matrix(3, 1, std::array{-8., -7., -5.}, f.plan);
    auto rhs = PrecisionBuffer<double>::binary64(std::move(b), core::NumericalGeneration{30});
    REQUIRE(rhs);
    const NoChecks<double> evaluator;
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{31}}, {}, evaluator};
    const DenseOperator op{a, context()};
    auto prepared = DensePrepared::prepare(op, f.scales(3), f.limits(3), f.plan);
    REQUIRE(prepared);
    auto result = solve_linear(op, Rhs{*rhs, context(), core::RhsGeneration{40}}, f.scales(3),
                               precision, f.accuracy(1e-29), f.limits(3), *prepared, f.plan);
    REQUIRE(result);
    CHECK(result->certificate().retention_path() == RetentionPath::PromotedTwoWord);
    CHECK(result->certificate().refinements() <= 5);
    const std::array<Big, 3> truth{Big(-3), Big(13) / 5, Big(0)};
    for (std::size_t i = 0; i < 3; ++i)
        CHECK(abs(Big(high(result->values()).values()[i]) + Big(low(result->values(), i)) -
                  truth[i]) < Big("1e-29"));
    CHECK(result->certificate().columns()[0].scaled_correction_upper > 0);
}
TEST_CASE("Small Schur condition is a separate unscaled diagnostic", "[verification][N4-05]") {
    Setup f;
    auto a = matrix(2, 2, std::array{1., 0., 0., 1.}, f.plan);
    auto u = matrix(2, 2, std::array{Complex{1}, Complex{0}, Complex{0}, Complex{2}}, f.plan);
    auto v = matrix(2, 2, std::array{Complex{1}, Complex{0}, Complex{0}, Complex{1}}, f.plan);
    auto prepared = UpdatedPrepared::prepare(ColumnUpdatedOperator{a, u, v, context()}, f.scales(2),
                                             f.limits(2, 1, 2), f.plan);
    REQUIRE(prepared);
    const auto estimate = prepared->schur_condition();
    REQUIRE(std::holds_alternative<ReciprocalOneNormEstimate>(estimate));
    CHECK(std::abs(std::get<ReciprocalOneNormEstimate>(estimate).value - 2. / 3) < 1e-14);
    auto moved = std::move(*prepared);
    CHECK(std::holds_alternative<NotEvaluated>(prepared->schur_condition()));
    CHECK(std::holds_alternative<ReciprocalOneNormEstimate>(moved.schur_condition()));
}
