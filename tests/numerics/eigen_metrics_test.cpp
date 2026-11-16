#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <uvlas/numerics/eigen_metrics.hpp>
#include <variant>

namespace {
using namespace uvlas;
namespace eig = numerics::eigen;
using Complex = std::complex<double>;
using Big = boost::multiprecision::cpp_bin_float_100;
using Exact = boost::multiprecision::cpp_rational;

[[nodiscard]] eig::MetricContext context() {
    return {{{core::ModelFingerprint{}, core::ModelGeneration{1}, core::GeometryGeneration{2},
              core::ChartGeneration{3}, core::SparsityGeneration{4}, core::NumericalGeneration{5},
              core::ScaleGeneration{6}, core::PrecisionGeneration{7}, core::Precision::Binary64,
              core::WakePolicyGeneration{8}, core::TimeIntervalGeneration{9},
              core::BackendGeneration{10}},
             core::TrialId{11},
             core::ScopeId{12},
             core::ScopeId{13},
             .25},
            {core::AcceptedStateId{20}, core::AcceptedStateId{21}},
            core::NumericalGeneration{40}};
}
[[nodiscard]] eig::MetricRequirement requirement(double budget, std::uint64_t id) {
    return {budget, core::ObservationId{id}, core::ScopeId{id + 100}};
}
[[nodiscard]] eig::PairMetricRequest pair_request() {
    return {requirement(1.e-12, 50), requirement(1.e-12, 51), eig::NotRequested{}};
}
[[nodiscard]] eig::BasisMetricRequest basis_request() {
    return {requirement(1.e-12, 50), requirement(1.e-12, 52), eig::NotRequested{}};
}
[[nodiscard]] eig::GeneralizedMetricRequest generalized_request() {
    return {requirement(1.e-12, 50), requirement(1.e-12, 51), requirement(1.e-12, 52),
            eig::NotRequested{}};
}
template <numerics::StorageScalar Scalar>
[[nodiscard]] Scalar scalar(double real, double imaginary = 0) {
    if constexpr (std::is_same_v<Scalar, double>)
        return real;
    else
        return {real, imaginary};
}
template <numerics::StorageScalar Scalar, std::size_t Size>
[[nodiscard]] eig::WordView<Scalar> word(const std::array<Scalar, Size>& high) {
    auto value = eig::WordView<Scalar>::make(high, {}, core::Precision::Binary64);
    REQUIRE(value);
    return *value;
}
template <numerics::StorageScalar Scalar, std::size_t Size>
[[nodiscard]] eig::WordView<Scalar> word(const std::array<Scalar, Size>& high,
                                         const std::array<Scalar, Size>& low) {
    auto value = eig::WordView<Scalar>::make(high, low, core::Precision::TwoWord);
    REQUIRE(value);
    return *value;
}
[[nodiscard]] double upper(const core::ErrorCertificate& proof) {
    const auto* achieved = std::get_if<core::Achieved>(&proof.assessment());
    if (achieved == nullptr)
        throw std::runtime_error("A successful metric must own achievement");
    return achieved->upper_bound;
}
void verify_proof(const core::ErrorCertificate& proof, eig::MetricRequirement request,
                  core::Precision representation) {
    CHECK(proof.context().observation == request.observation);
    CHECK(proof.context().scope == request.scope);
    CHECK(proof.context().numerical_generation == core::NumericalGeneration{40});
    CHECK(proof.context().representation == representation);
    CHECK(std::bit_cast<std::uint64_t>(proof.budget()) ==
          std::bit_cast<std::uint64_t>(request.budget));
    CHECK(std::isfinite(upper(proof)));
    CHECK(upper(proof) >= 0);
    CHECK(upper(proof) <= request.budget);
}
void encloses_tightly(const core::ErrorCertificate& proof, const Big& exact) {
    const Big bound{upper(proof)};
    CHECK(bound >= exact);
    // Fixed independent accuracy check: an arbitrary large conservative bound is insufficient.
    CHECK(bound <=
          exact * (1 + Big("1e-12")) + Big(64) * Big(std::numeric_limits<double>::denorm_min()));
}
void failure_context(const eig::EigenFailure& failure) {
    CHECK(failure.cause.operation == core::Operation::ExtractEigen);
    CHECK(failure.cause.details.model == core::ModelGeneration{1});
    CHECK(failure.cause.details.trial == core::TrialId{11});
    CHECK(failure.cause.details.numerical_generation == core::NumericalGeneration{40});
    CHECK(failure.cause.details.parent == core::AcceptedStateId{20});
    CHECK(failure.cause.details.last_accepted == core::AcceptedStateId{21});
}

template <numerics::StorageScalar Scalar>
void verify_optional_backward() {
    const double epsilon = std::ldexp(1., -80);
    const std::array vector{scalar<Scalar>(1), scalar<Scalar>(epsilon, epsilon)};
    const std::array value{scalar<Scalar>(5)};
    const std::array applied{scalar<Scalar>(5), scalar<Scalar>(epsilon, epsilon)};
    const std::array<double, 2> errors{}, components{1, 1};
    const auto original = vector;
    const Big exact_residual = [epsilon]() -> Big {
        if constexpr (std::is_same_v<Scalar, Complex>)
            return Big(4) * Big(epsilon) * sqrt(Big(2));
        else
            return Big(4) * Big(epsilon);
    }();
    // A=diag(5,1), lambda=5. Exact r=(0,-4*epsilon), or its complex counterpart.
    // beta_second=4/6=2/3 for every nonzero epsilon, however small the norm residual.
    const Exact exact_beta = Exact{2} / 3;
    for (const auto side : {eig::PairSide::Right, eig::PairSide::Left}) {
        auto request = pair_request();
        const auto norm_only = eig::assess_pair_metrics(
            word(vector), word(value), eig::ActionValues<Scalar>{word(applied), errors},
            eig::MetricScales{components, 1}, side, request, context());
        REQUIRE(norm_only);
        CHECK(std::holds_alternative<eig::NotRequested>(norm_only->backward));
        verify_proof(norm_only->residual, request.residual, core::Precision::Binary64);
        verify_proof(norm_only->normalization, request.normalization, core::Precision::Binary64);
        encloses_tightly(norm_only->residual, exact_residual);
        request.backward = requirement(1.e-12, 53);
        const auto refused = eig::assess_pair_metrics(
            word(vector), word(value), eig::ActionValues<Scalar>{word(applied), errors},
            eig::MetricScales{components, 1}, side, request, context());
        REQUIRE_FALSE(refused);
        CHECK(refused.error().cause.code == core::FailureCode::PrecisionUnresolved);
        CHECK(refused.error().cause.details.problem == core::ProblemId{core::ObservationId{53}});
        failure_context(refused.error());
        request.backward = requirement(.75, 53);
        const auto accepted = eig::assess_pair_metrics(
            word(vector), word(value), eig::ActionValues<Scalar>{word(applied), errors},
            eig::MetricScales{components, 1}, side, request, context());
        REQUIRE(accepted);
        const auto* beta = std::get_if<core::ErrorCertificate>(&accepted->backward);
        if (beta == nullptr) {
            FAIL_CHECK("The requested backward criterion must own its proof");
            return;
        }
        verify_proof(*beta, requirement(.75, 53), core::Precision::Binary64);
        CHECK(Exact{upper(*beta)} >= exact_beta);
        encloses_tightly(*beta, Big(2) / 3);
        for (std::size_t i = 0; i < vector.size(); ++i)
            CHECK(vector[i] == original[i]);
        // A later refusal/qualification must not overwrite an earlier owned proof.
        encloses_tightly(norm_only->residual, exact_residual);
    }
}

template <numerics::StorageScalar Scalar>
void verify_scales_and_side() {
    const double epsilon = std::ldexp(1., -50);
    const std::array components{2., 4.};
    const std::array<double, 2> errors{};
    const auto request = pair_request();
    const std::array value{scalar<Scalar>(6)};
    const std::array complex_value{scalar<Scalar>(2, 3)};
    for (const auto side : {eig::PairSide::Right, eig::PairSide::Left}) {
        const bool right = side == eig::PairSide::Right;
        const std::array vector{scalar<Scalar>(0), scalar<Scalar>(right ? 4. : .25)};
        const auto lambda = std::is_same_v<Scalar, double> ? value : complex_value;
        Scalar second;
        if constexpr (std::is_same_v<Scalar, double>)
            second = right ? 24. : 1.5;
        else
            second = right ? Complex{8, 12} : Complex{.5, -.75};
        const std::array applied{scalar<Scalar>(right ? 2 * epsilon : epsilon / 2), second};
        // One literal physical matrix can produce both actions:
        // A00=7, A01=epsilon/2, A10=2*epsilon, A11=the supplied lambda.
        // D^-1 right and D left both give rho=(epsilon/3,0), with vector infinity norm 1.
        const auto result = eig::assess_pair_metrics(
            word(vector), word(lambda), eig::ActionValues<Scalar>{word(applied), errors},
            eig::MetricScales{components, 3}, side, request, context());
        REQUIRE(result);
        encloses_tightly(result->residual, Big(epsilon) / 3);
        CHECK(result->vector_infinity_lower <= 1);
        CHECK(result->vector_infinity_lower > .999999999999);
    }
}
} // namespace

TEST_CASE("Eigen norm qualification does not impose an unrequested componentwise criterion",
          "[verification][negative][N4-07][N4-08]") {
    verify_optional_backward<double>();
    verify_optional_backward<Complex>();
}
TEST_CASE("Eigen pair metrics use physical side scales and the conjugated left eigenvalue",
          "[verification][N4-07]") {
    verify_scales_and_side<double>();
    verify_scales_and_side<Complex>();
}
TEST_CASE("Eigen pair metrics consume cancellation in vector eigenvalue and action low words",
          "[verification][negative][N4-06][N4-07]") {
    // Actual phi=(0,1), lambda=1, Aphi=(0,1). Every high-only interpretation is wrong.
    const std::array high{1., 1.}, low{-1., 0.};
    const std::array lambda_hi{5.}, lambda_lo{-4.};
    const std::array action_hi{1., 2.}, action_lo{-1., -1.};
    const std::array<double, 2> errors{}, components{1, 1};
    const auto result = eig::assess_pair_metrics(
        word(high, low), word(lambda_hi, lambda_lo),
        eig::ActionValues<double>{word(action_hi, action_lo), errors},
        eig::MetricScales{components, 1}, eig::PairSide::Right, pair_request(), context());
    REQUIRE(result);
    verify_proof(result->residual, pair_request().residual, core::Precision::TwoWord);
    verify_proof(result->normalization, pair_request().normalization, core::Precision::TwoWord);
    CHECK(high[0] == 1);
    CHECK(low[0] == -1);
    CHECK(lambda_hi[0] == 5);
    CHECK(lambda_lo[0] == -4);
    CHECK(action_lo[1] == -1);
}
TEST_CASE("Eigen original producer uncertainty remains required under a zero residual budget",
          "[verification][negative][N4-07]") {
    const double epsilon = std::ldexp(1., -80);
    const std::array vector{1., 0.}, applied{5., 0.}, components{1., 1.};
    const std::array lambda{5.};
    const std::array errors{epsilon, 0.};
    auto request = pair_request();
    const auto bounded = eig::assess_pair_metrics(
        word(vector), word(lambda), eig::ActionValues<double>{word(applied), errors},
        eig::MetricScales{components, 1}, eig::PairSide::Right, request, context());
    REQUIRE(bounded);
    encloses_tightly(bounded->residual, Big(epsilon));
    request.residual.budget = 0;
    const auto refused = eig::assess_pair_metrics(
        word(vector), word(lambda), eig::ActionValues<double>{word(applied), errors},
        eig::MetricScales{components, 1}, eig::PairSide::Right, request, context());
    REQUIRE_FALSE(refused);
    CHECK(refused.error().cause.code == core::FailureCode::PrecisionUnresolved);
    CHECK(refused.error().cause.details.problem == core::ProblemId{core::ObservationId{50}});
}
TEST_CASE("Eigen Schur metrics qualify a noncontiguous invariant space without losing conjugation",
          "[verification][negative][N4-07]") {
    for (const double delta : {0., std::ldexp(1., -40)}) {
        // A=[[1,1,1],[0,3,0],[0,0,1+delta]], selected space span(e0,e2).
        const std::array q{1., 0., 0., 0., 0., 1.};
        const std::array t{1., 0., 1., 1. + delta};
        const std::array aq{1., 0., 0., 1., 0., 1. + delta};
        const std::array<double, 6> action_error{};
        const std::array<double, 4> t_error{};
        const auto result = eig::assess_schur_metrics(
            3, 2, word(q), eig::ActionValues<double>{word(t), t_error},
            eig::ActionValues<double>{word(aq), action_error}, basis_request(), context());
        REQUIRE(result);
        verify_proof(result->residual, basis_request().residual, core::Precision::Binary64);
        verify_proof(result->orthogonality, basis_request().orthogonality,
                     core::Precision::Binary64);
        CHECK(std::holds_alternative<eig::NotRequested>(result->backward));
        // Rotate only the first basis vector by i: Q*T still equals the same original A*Q.
        const std::array cq{Complex{0, 1}, Complex{}, Complex{}, Complex{}, Complex{}, Complex{1}};
        const std::array ct{Complex{1}, Complex{}, Complex{0, -1}, Complex{1. + delta}};
        const std::array caq{Complex{0, 1}, Complex{}, Complex{},
                             Complex{1},    Complex{}, Complex{1. + delta}};
        const auto complex_result = eig::assess_schur_metrics(
            3, 2, word(cq), eig::ActionValues<Complex>{word(ct), t_error},
            eig::ActionValues<Complex>{word(caq), action_error}, basis_request(), context());
        REQUIRE(complex_result);
        verify_proof(complex_result->orthogonality, basis_request().orthogonality,
                     core::Precision::Binary64);
    }
}
TEST_CASE("Eigen generalized metrics keep the original mass inner product and both action words",
          "[verification][N4-06][N4-07]") {
    // M=[[2,2],[2,5]], K=[[-2,-2],[-2,10]]. The exact eigendirections
    // e0 and (-1,1) have M squared norms 2 and 3, not unit Euclidean norms.
    const double u = 1 / std::sqrt(2.), v = 1 / std::sqrt(3.);
    const std::array phi{u, 0., -v, v};
    const std::array lambda{-1., 4.};
    const std::array kh{-2 * u, -2 * u, 0., 12 * v}, mh{2 * u, 2 * u, 0., 3 * v};
    const std::array kl{0., 0., 0., std::fma(12., v, -kh[3])};
    const std::array ml{0., 0., 0., std::fma(3., v, -mh[3])};
    const std::array<double, 4> errors{};
    const std::array components{2., 4.};
    const auto request = generalized_request();
    const auto result = eig::assess_generalized_metrics(
        2, word(phi), word(lambda), eig::ActionValues<double>{word(kh, kl), errors},
        eig::ActionValues<double>{word(mh, ml), errors}, eig::MetricScales{components, 3}, 1.,
        request, context());
    REQUIRE(result);
    verify_proof(result->residual, request.residual, core::Precision::Binary64);
    verify_proof(result->normalization, request.normalization, core::Precision::Binary64);
    verify_proof(result->orthogonality, request.orthogonality, core::Precision::Binary64);
    const Big first = abs(Big(2) * Big(u) * Big(u) - 1);
    const Big second = abs(Big(3) * Big(v) * Big(v) - 1);
    CHECK(Big(upper(result->normalization)) >= first);
    CHECK(Big(upper(result->normalization)) >= second);
    CHECK(Big(upper(result->orthogonality)) >= first);
    CHECK(Big(upper(result->orthogonality)) >= second);
    CHECK(std::holds_alternative<eig::NotRequested>(result->backward));
}
TEST_CASE("Eigen basis emptiness is qualified without inventing a nonzero vector norm",
          "[unit][negative][N4-07]") {
    const std::array<double, 0> empty{};
    const auto empty_action = eig::ActionValues<double>{word(empty), empty};
    for (const auto rows : {core::Index{0}, core::Index{4}}) {
        const auto basis = eig::assess_schur_metrics(rows, 0, word(empty), empty_action,
                                                     empty_action, basis_request(), context());
        REQUIRE(basis);
        CHECK(upper(basis->residual) == 0);
        CHECK(upper(basis->orthogonality) == 0);
    }
    const auto generalized = eig::assess_generalized_metrics(
        0, word(empty), word(empty), empty_action, empty_action, eig::MetricScales{empty, 1}, 1.,
        generalized_request(), context());
    REQUIRE(generalized);
    CHECK(upper(generalized->residual) == 0);
    CHECK(upper(generalized->normalization) == 0);
    CHECK(upper(generalized->orthogonality) == 0);
    const std::array lambda{1.};
    const auto pair = eig::assess_pair_metrics(word(empty), word(lambda), empty_action,
                                               eig::MetricScales{empty, 1}, eig::PairSide::Right,
                                               pair_request(), context());
    REQUIRE_FALSE(pair);
    CHECK(pair.error().cause.code == core::FailureCode::InvalidInput);
}

namespace {
template <class Metrics>
void rejected(const eig::Result<Metrics>& result, core::FailureCode expected) {
    REQUIRE_FALSE(result);
    CHECK(result.error().cause.code == expected);
    failure_context(result.error());
}
} // namespace
TEST_CASE("Eigen pair metrics reject all invalid raw bounds scales requests and contexts first",
          "[negative][N4-07][N4-08]") {
    enum class Site : std::uint8_t {
        NegativeError,
        NanError,
        InfiniteError,
        ZeroComponent,
        NegativeComponent,
        NanComponent,
        InfiniteComponent,
        ZeroEigenScale,
        NegativeEigenScale,
        NanEigenScale,
        InfiniteEigenScale,
        NegativeResidual,
        NanResidual,
        InfiniteResidual,
        NegativeNormalization,
        NanNormalization,
        InfiniteNormalization,
        NegativeBackward,
        NanBackward,
        InfiniteBackward,
        UnknownSide,
        NanTime,
        UnknownRepresentation
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (const auto site : {Site::NegativeError,         Site::NanError,
                            Site::InfiniteError,         Site::ZeroComponent,
                            Site::NegativeComponent,     Site::NanComponent,
                            Site::InfiniteComponent,     Site::ZeroEigenScale,
                            Site::NegativeEigenScale,    Site::NanEigenScale,
                            Site::InfiniteEigenScale,    Site::NegativeResidual,
                            Site::NanResidual,           Site::InfiniteResidual,
                            Site::NegativeNormalization, Site::NanNormalization,
                            Site::InfiniteNormalization, Site::NegativeBackward,
                            Site::NanBackward,           Site::InfiniteBackward,
                            Site::UnknownSide,           Site::NanTime,
                            Site::UnknownRepresentation}) {
        const std::array vector{1., 0.};
        const std::array lambda{5.};
        // The first row already misses its norm budget. Later raw invalidity must
        // still be reported, rather than hidden by early qualification failure.
        const std::array applied{6., 0.};
        std::array<double, 2> error{}, components{1, 1};
        double eigen_scale = 1;
        auto request = pair_request();
        auto identity = context();
        auto side = eig::PairSide::Right;
        auto expected = core::FailureCode::InvalidInput;
        switch (site) {
        case Site::NegativeError:
            error[1] = -1;
            break;
        case Site::NanError:
            error[1] = nan;
            expected = core::FailureCode::Nonfinite;
            break;
        case Site::InfiniteError:
            error[1] = inf;
            expected = core::FailureCode::Nonfinite;
            break;
        case Site::ZeroComponent:
            components[1] = 0;
            break;
        case Site::NegativeComponent:
            components[1] = -1;
            break;
        case Site::NanComponent:
            components[1] = nan;
            expected = core::FailureCode::Nonfinite;
            break;
        case Site::InfiniteComponent:
            components[1] = inf;
            expected = core::FailureCode::Nonfinite;
            break;
        case Site::ZeroEigenScale:
            eigen_scale = 0;
            break;
        case Site::NegativeEigenScale:
            eigen_scale = -1;
            break;
        case Site::NanEigenScale:
            eigen_scale = nan;
            expected = core::FailureCode::Nonfinite;
            break;
        case Site::InfiniteEigenScale:
            eigen_scale = inf;
            expected = core::FailureCode::Nonfinite;
            break;
        case Site::NegativeResidual:
            request.residual.budget = -1;
            break;
        case Site::NanResidual:
            request.residual.budget = nan;
            expected = core::FailureCode::Nonfinite;
            break;
        case Site::InfiniteResidual:
            request.residual.budget = inf;
            expected = core::FailureCode::Nonfinite;
            break;
        case Site::NegativeNormalization:
            request.normalization.budget = -1;
            break;
        case Site::NanNormalization:
            request.normalization.budget = nan;
            expected = core::FailureCode::Nonfinite;
            break;
        case Site::InfiniteNormalization:
            request.normalization.budget = inf;
            expected = core::FailureCode::Nonfinite;
            break;
        case Site::NegativeBackward:
            request.backward = requirement(-1, 53);
            break;
        case Site::NanBackward:
            request.backward = requirement(nan, 53);
            expected = core::FailureCode::Nonfinite;
            break;
        case Site::InfiniteBackward:
            request.backward = requirement(inf, 53);
            expected = core::FailureCode::Nonfinite;
            break;
        case Site::UnknownSide:
            side = static_cast<eig::PairSide>(255);
            break;
        case Site::NanTime:
            identity.system.time = nan;
            expected = core::FailureCode::Nonfinite;
            break;
        case Site::UnknownRepresentation:
            identity.system.cache.representation = static_cast<core::Precision>(255);
            expected = core::FailureCode::PrecisionUnsupported;
            break;
        }
        const auto result = eig::assess_pair_metrics(
            word(vector), word(lambda), eig::ActionValues<double>{word(applied), error},
            eig::MetricScales{components, eigen_scale}, side, request, identity);
        rejected(result, expected);
        CHECK(vector == std::array{1., 0.});
        CHECK(applied == std::array{6., 0.});
    }
}
TEST_CASE("Eigen pair metric dimensions and zero candidates fail without hiding their category",
          "[negative][N4-07][N4-08]") {
    const std::array vector{1., 0.}, applied{5., 0.}, components{1., 1.};
    const std::array lambda{5.};
    const std::array<double, 2> errors{}, zeros{};
    const std::array<double, 0> empty{};
    const std::array short_values{1.};
    const std::array long_lambda{5., 1.};
    const auto run = [&](eig::WordView<double> v, eig::WordView<double> l, eig::WordView<double> a,
                         std::span<const double> e, std::span<const double> d) {
        return eig::assess_pair_metrics(v, l, eig::ActionValues<double>{a, e},
                                        eig::MetricScales{d, 1}, eig::PairSide::Right,
                                        pair_request(), context());
    };
    rejected(run(word(vector), word(empty), word(applied), errors, components),
             core::FailureCode::InvalidInput);
    rejected(run(word(vector), word(long_lambda), word(applied), errors, components),
             core::FailureCode::InvalidInput);
    rejected(run(word(vector), word(lambda), word(short_values), errors, components),
             core::FailureCode::InvalidInput);
    rejected(run(word(vector), word(lambda), word(applied), std::span{errors}.first(1), components),
             core::FailureCode::InvalidInput);
    rejected(run(word(vector), word(lambda), word(applied), errors, std::span{components}.first(1)),
             core::FailureCode::InvalidInput);
    rejected(run(word(zeros), word(lambda), word(zeros), errors, components),
             core::FailureCode::PrecisionUnresolved);
}
TEST_CASE("Eigen Schur shape triangularity and rank guards precede false qualification",
          "[negative][N4-07][N4-08]") {
    const std::array q{1., 0., 0., 0., 0., 1.}, aq{1., 0., 0., 1., 0., 2.};
    const std::array t{1., 0., 1., 2.};
    const std::array<double, 6> error{};
    const std::array<double, 4> t_error{};
    const auto run = [&](core::Index n, core::Index k, eig::WordView<double> basis,
                         eig::WordView<double> triangular, std::span<const double> te,
                         eig::WordView<double> action, std::span<const double> ae) {
        return eig::assess_schur_metrics(n, k, basis, eig::ActionValues<double>{triangular, te},
                                         eig::ActionValues<double>{action, ae}, basis_request(),
                                         context());
    };
    const std::array short_q{1., 0., 0., 0., 1.};
    const std::array short_t{1., 0., 1.};
    rejected(run(3, 2, word(short_q), word(t), t_error, word(aq), error),
             core::FailureCode::InvalidInput);
    rejected(run(3, 2, word(q), word(short_t), t_error, word(aq), error),
             core::FailureCode::InvalidInput);
    rejected(run(3, 2, word(q), word(t), std::span{t_error}.first(3), word(aq), error),
             core::FailureCode::InvalidInput);
    rejected(run(3, 2, word(q), word(t), t_error, word(short_q), error),
             core::FailureCode::InvalidInput);
    rejected(run(3, 2, word(q), word(t), t_error, word(aq), std::span{error}.first(5)),
             core::FailureCode::InvalidInput);
    const std::array too_wide{1., 1.};
    const std::array<double, 2> two_error{};
    rejected(run(1, 2, word(too_wide), word(t), t_error, word(too_wide), two_error),
             core::FailureCode::InvalidInput);
    const std::array<double, 4> t_low{0., std::ldexp(1., -80), 0., 0.};
    rejected(run(3, 2, word(q), word(t, t_low), t_error, word(aq), error),
             core::FailureCode::InvalidInput);
    std::array bad_error = t_error;
    bad_error[3] = std::numeric_limits<double>::quiet_NaN();
    rejected(run(3, 2, word(q), word(t), bad_error, word(aq), error), core::FailureCode::Nonfinite);
    std::array uncertain_lower = t_error;
    uncertain_lower[1] = 1.e-30;
    rejected(run(3, 2, word(q), word(t), uncertain_lower, word(aq), error),
             core::FailureCode::PrecisionUnresolved);

    // A=I, duplicated unit columns: residual zero and any loose budget cannot prove rank two.
    const std::array dependent{1., 0., 1., 0.}, identity{1., 0., 0., 1.};
    auto loose = basis_request();
    loose.orthogonality.budget = 2;
    const auto rank = eig::assess_schur_metrics(
        2, 2, word(dependent), eig::ActionValues<double>{word(identity), t_error},
        eig::ActionValues<double>{word(dependent), t_error}, loose, context());
    rejected(rank, core::FailureCode::RankUnresolved);
    // Each Gram component q^2=9*2^1020 is finite; its two-term row sum exceeds max-double.
    const double large = std::ldexp(3., 510);
    const std::array large_q{large, 0., large, 0.};
    const std::array<double, 4> zero{};
    const auto overflow = eig::assess_schur_metrics(
        2, 2, word(large_q), eig::ActionValues<double>{word(zero), t_error},
        eig::ActionValues<double>{word(zero), t_error}, loose, context());
    rejected(overflow, core::FailureCode::PrecisionUnresolved);
}
TEST_CASE("Eigen generalized metric shape mass-scale and full-rank requirements are explicit",
          "[negative][N4-07][N4-08]") {
    const std::array identity{1., 0., 0., 1.};
    const std::array lambda{1., 1.}, components{1., 1.};
    const std::array<double, 4> error{};
    const auto run = [&](eig::WordView<double> basis, eig::WordView<double> values,
                         eig::WordView<double> stiffness, std::span<const double> ke,
                         eig::WordView<double> mass, std::span<const double> me,
                         double mass_scale) {
        return eig::assess_generalized_metrics(
            2, basis, values, eig::ActionValues<double>{stiffness, ke},
            eig::ActionValues<double>{mass, me}, eig::MetricScales{components, 1}, mass_scale,
            generalized_request(), context());
    };
    const std::array short_basis{1., 0., 0.};
    const std::array short_lambda{1.};
    rejected(run(word(short_basis), word(lambda), word(identity), error, word(identity), error, 1),
             core::FailureCode::InvalidInput);
    rejected(
        run(word(identity), word(short_lambda), word(identity), error, word(identity), error, 1),
        core::FailureCode::InvalidInput);
    rejected(run(word(identity), word(lambda), word(short_basis), error, word(identity), error, 1),
             core::FailureCode::InvalidInput);
    rejected(run(word(identity), word(lambda), word(identity), std::span{error}.first(3),
                 word(identity), error, 1),
             core::FailureCode::InvalidInput);
    rejected(run(word(identity), word(lambda), word(identity), error, word(short_basis), error, 1),
             core::FailureCode::InvalidInput);
    rejected(run(word(identity), word(lambda), word(identity), error, word(identity),
                 std::span{error}.first(3), 1),
             core::FailureCode::InvalidInput);
    for (const double mass_scale : {0., -1., std::numeric_limits<double>::quiet_NaN(),
                                    std::numeric_limits<double>::infinity()}) {
        const auto expected = std::isfinite(mass_scale) ? core::FailureCode::InvalidInput
                                                        : core::FailureCode::Nonfinite;
        rejected(run(word(identity), word(lambda), word(identity), error, word(identity), error,
                     mass_scale),
                 expected);
    }
    std::array bad_error = error;
    bad_error[3] = std::numeric_limits<double>::infinity();
    rejected(run(word(identity), word(lambda), word(identity), error, word(identity), bad_error, 1),
             core::FailureCode::Nonfinite);
    const std::array dependent{1., 0., 1., 0.};
    auto loose = generalized_request();
    loose.orthogonality.budget = 2;
    const auto rank = eig::assess_generalized_metrics(
        2, word(dependent), word(lambda), eig::ActionValues<double>{word(dependent), error},
        eig::ActionValues<double>{word(dependent), error}, eig::MetricScales{components, 1}, 1,
        loose, context());
    rejected(rank, core::FailureCode::RankUnresolved);
    // M=I, K=0, lambda=0: each diagonal Gram entry is finite, but the
    // complete row-sum bound is not. Loose finite budgets isolate that failure.
    const double large = std::ldexp(3., 510);
    const std::array large_basis{large, 0., large, 0.};
    const std::array<double, 4> zero{};
    const std::array<double, 2> zero_lambda{};
    loose.normalization.budget = std::numeric_limits<double>::max();
    loose.orthogonality.budget = std::numeric_limits<double>::max();
    const auto overflow = eig::assess_generalized_metrics(
        2, word(large_basis), word(zero_lambda), eig::ActionValues<double>{word(zero), error},
        eig::ActionValues<double>{word(large_basis), error}, eig::MetricScales{components, 1}, 1,
        loose, context());
    rejected(overflow, core::FailureCode::PrecisionUnresolved);
}
