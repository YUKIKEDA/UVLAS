#include "eigen_full_phase_private.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <uvlas/numerics/eigen_full.hpp>
#include <vector>
namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::eigen;
using Complex = std::complex<double>;
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
FullRequest request(core::NumericalGeneration generation = core::NumericalGeneration{42}) {
    return {BasisMetricRequest{{1e-12, core::ObservationId{1}, core::ScopeId{1}},
                               {1e-12, core::ObservationId{2}, core::ScopeId{1}},
                               NotRequested{}},
            linear::RetainedPrecision{generation},
            FullLimits{96,
                       96,
                       {10000, 10000, 0, 10},
                       linear::EstimateAdmission{128 * std::uint64_t{1024} * 1024}},
            {},
            nullptr,
            false};
}
struct Scales {
    explicit Scales(std::size_t n) : d(n, 1), units(n), aunits(n * n) {}
    double lambda_scale = 1;
    std::vector<double> d;
    std::vector<linear::UnitDimensions> units, aunits;
    [[nodiscard]] FullScales view() const {
        return {d, lambda_scale, units, {}, aunits};
    }
};
template <StorageScalar Scalar>
PrecisionBuffer<Scalar> input(core::Index n, std::span<const Scalar> values,
                              core::ResourcePlan& plan) {
    auto matrix = DenseMatrix<Scalar>::copy({n, n, n}, values, plan, core::Operation::ExtractEigen);
    if (!matrix)
        throw std::runtime_error("test matrix input failed");
    auto result =
        PrecisionBuffer<Scalar>::binary64(std::move(*matrix), core::NumericalGeneration{5});
    if (!result)
        throw std::runtime_error("test precision input failed");
    return std::move(*result);
}
} // namespace
TEST_CASE("Full Schur extracts all real roots without clipping the tiny negative root",
          "[unit][N4-07][N4-06][EIG-01][EIG-04][EIG-05]") {
    core::ResourcePlan plan{64 * std::uint64_t{1024} * 1024};
    std::array<double, 36> a{};
    a[1] = 2;
    a[6] = -2;
    a[2 + 2 * 6] = 3;
    a[2 + 3 * 6] = 1;
    a[3 + 3 * 6] = -1;
    a[5 + 5 * 6] = -0x1p-60;
    auto original = input<double>(6, a, plan);
    Scales const scales{6};
    const auto result = full_spectrum(original, context(), {}, scales.view(), request(), plan);
    if (!result) {
        CAPTURE(result.error().cause.code);
        FAIL("full Schur failed");
    }
    const auto v = result->view();
    REQUIRE(v);
    if (!v.has_value())
        throw std::runtime_error("missing expected v");
    REQUIRE(v->eigenvalues.size() == 6);
    const std::array expected{Complex{-1, 0}, Complex{-0x1p-60, 0}, Complex{0, -2},
                              Complex{0, 0},  Complex{0, 2},        Complex{3, 0}};
    std::array<bool, 6> used{};
    for (const auto value : v->eigenvalues) {
        std::size_t match = expected.size();
        for (std::size_t i = 0; i < expected.size(); ++i)
            if (!used[i] && std::abs(value - expected[i]) <= 1e-12) {
                match = i;
                break;
            }
        REQUIRE(match < expected.size());
        used[match] = true;
    }
    CHECK(std::count(v->eigenvalues.begin(), v->eigenvalues.end(), Complex{}) == 1);
    CHECK(v->eigenvalues[1].real() < 0);
    CHECK(std::abs(v->eigenvalues[1].real() + 0x1p-60) <= 0x1p-60 * 1e-12);
    REQUIRE(result->basis_proof());
    CHECK(result->diagnostics().converged_count == 6);
    CHECK(result->physical_coverage() == FullPhysicalCoverage::Unresolved);
}

TEST_CASE("Noncontiguous Schur clusters preserve Jordan multiplicity and independent left spaces",
          "[unit][N4-07][N4-06][EIG-04][EIG-05][EIG-07][EIG-13]") {
    for (const double delta : {0.0, 0x1p-40}) {
        core::ResourcePlan plan{64 * std::uint64_t{1024} * 1024};
        const std::array<double, 9> a{1, 0, 0, 1, 3, 0, 1, 0, 1 + delta};
        auto original = input<double>(3, a, plan);
        Scales const scales{3};
        auto full = full_spectrum(original, context(), {}, scales.view(), request(), plan);
        REQUIRE(full);
        const std::array<core::Index, 2> indices{0, 1};
        auto selected = FullSelection::make(*full, indices, plan);
        REQUIRE(selected);
        FullClusterRequest const wanted{
            request(core::NumericalGeneration{43}), {1e-12, 1e-12}, true};
        auto cluster = full_cluster(*full, *selected, wanted, plan);
        if (!cluster) {
            CAPTURE(delta, cluster.error().cause.code,
                    cluster.error().cause.details.problem.has_value(),
                    cluster.error().diagnostics.extraction_attempts,
                    cluster.error().diagnostics.forward_applications,
                    cluster.error().diagnostics.adjoint_applications);
            FAIL_CHECK("cluster failed");
            continue;
        }
        const auto view = cluster->view();
        REQUIRE(view);
        if (!view.has_value())
            throw std::runtime_error("missing expected view");
        CHECK(view->columns == 2);
        REQUIRE(view->right.size() == 6);
        REQUIRE(view->left.size() == 6);
        const std::array<double, 9> right_projector{1, 0, 0, 0, 0, 0, 0, 0, 1};
        const std::array<double, 9> left_projector{.8, -.4, 0, -.4, .2, 0, 0, 0, 1};
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t i = 0; i < 3; ++i) {
                Complex right{}, left{};
                for (std::size_t k = 0; k < 2; ++k) {
                    right += view->right[i + k * 3] * std::conj(view->right[j + k * 3]);
                    left += view->left[i + k * 3] * std::conj(view->left[j + k * 3]);
                }
                CHECK(std::abs(right - right_projector[i + j * 3]) <= 1e-12);
                CHECK(std::abs(left - left_projector[i + j * 3]) <= 1e-12);
            }
        REQUIRE(cluster->basis_proof());
        REQUIRE(cluster->left_basis_proof());
        CHECK(view->left_triangular.size() == 4);
    }
}

namespace {
using Big = boost::multiprecision::cpp_bin_float_100;
struct BC {
    Big r, i;
};
BC big(Complex z) {
    return {Big{z.real()}, Big{z.imag()}};
}
BC add(const BC& a, const BC& b) {
    return {a.r + b.r, a.i + b.i};
}
BC minus(const BC& a, const BC& b) {
    return {a.r - b.r, a.i - b.i};
}
BC product(const BC& a, const BC& b) {
    return {a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r};
}
BC quotient(const BC& a, const Big& b) {
    return {a.r / b, a.i / b};
}
BC conjugate(const BC& a) {
    return {a.r, -a.i};
}
Big absolute(const BC& a) {
    return boost::multiprecision::sqrt(a.r * a.r + a.i * a.i);
}
void original_schur_reference(const FullPacketView& v, std::span<const Complex> original,
                              const Scales& scales, std::span<const core::Index> mapping,
                              bool left = false) {
    const auto n = static_cast<std::size_t>(v.rows), k = static_cast<std::size_t>(v.columns);
    const auto q = left ? v.left : v.right;
    const auto t = left ? v.left_triangular : v.triangular;
    Big residual = 0, qnorm = 0, gram = 0;
    for (std::size_t i = 0; i < n; ++i) {
        Big row = 0, qr = 0;
        for (std::size_t j = 0; j < k; ++j) {
            BC difference{};
            for (std::size_t m = 0; m < n; ++m) {
                const BC coefficient =
                    left ? conjugate(big(original[m + i * n])) : big(original[i + m * n]);
                const Big numerator = left ? Big{scales.d[i]} : Big{scales.d[m]};
                const Big denominator =
                    (left ? Big{scales.d[m]} : Big{scales.d[i]}) * Big{scales.lambda_scale};
                difference = add(
                    difference, product(quotient(product(coefficient, {numerator, 0}), denominator),
                                        big(q[m + j * n])));
            }
            for (std::size_t m = 0; m < k; ++m) {
                BC coefficient = big(t[m + j * k]);
                if (m == j) {
                    const auto found =
                        std::find(mapping.begin(), mapping.end(), static_cast<core::Index>(j));
                    REQUIRE(found != mapping.end());
                    const auto index = static_cast<std::size_t>(found - mapping.begin());
                    coefficient =
                        quotient(big(left ? std::conj(v.eigenvalues[index]) : v.eigenvalues[index]),
                                 Big{scales.lambda_scale});
                }
                difference = minus(difference, product(big(q[i + m * n]), coefficient));
            }
            row += absolute(difference);
            qr += absolute(big(q[i + j * n]));
        }
        residual = std::max(residual, row);
        qnorm = std::max(qnorm, qr);
    }
    for (std::size_t j = 0; j < k; ++j) {
        Big row = 0;
        for (std::size_t m = 0; m < k; ++m) {
            BC value{};
            for (std::size_t i = 0; i < n; ++i)
                value = add(value, product(conjugate(big(q[i + j * n])), big(q[i + m * n])));
            if (j == m)
                value.r -= 1;
            row += absolute(value);
        }
        gram = std::max(gram, row);
    }
    if (k > 0) {
        CHECK(residual / qnorm <= Big{"1e-12"});
        CHECK(gram <= Big{"1e-12"});
    }
}
void root_matching(std::span<const Complex> actual, std::span<const Complex> expected) {
    REQUIRE(actual.size() == expected.size());
    std::vector<bool> used(expected.size(), false);
    for (const auto z : actual) {
        std::size_t match = expected.size();
        for (std::size_t i = 0; i < expected.size(); ++i)
            if (!used[i] &&
                std::abs(z - expected[i]) <= 1e-12 * std::max(1.0, std::abs(expected[i]))) {
                match = i;
                break;
            }
        REQUIRE(match < expected.size());
        used[match] = true;
    }
}
} // namespace

TEST_CASE(
    "Complex full roots and nonbinary scales are qualified against the original literal matrix",
    "[verification][N4-07][N4-06][EIG-02][EIG-03][EIG-05][EIG-13]") {
    core::ResourcePlan plan{64 * std::uint64_t{1024} * 1024};
    std::array<Complex, 36> a{};
    const std::array diag{Complex{7, 1}, Complex{3, 0},  Complex{2, 0},
                          Complex{1, 0}, Complex{.5, 0}, Complex{.25, 0}};
    for (std::size_t i = 0; i < 6; ++i)
        a[i + i * 6] = diag[i];
    a[6] = Complex{0, 1};
    auto original = input<Complex>(6, a, plan);
    Scales scales{6};
    scales.d = {3, .2, 5, .5, 7, 11};
    scales.lambda_scale = 7;
    const auto result = full_spectrum(original, context(), {}, scales.view(), request(), plan);
    if (!result) {
        CAPTURE(result.error().cause.code);
        FAIL_CHECK("complex spectrum failed");
        return;
    }
    const auto v = result->view();
    REQUIRE(v);
    if (!v.has_value())
        throw std::runtime_error("missing expected v");
    root_matching(v->eigenvalues, diag);
    original_schur_reference(*v, a, scales, result->diagonal_mapping());
    CHECK(v->eigenvalues.back().imag() > 0);
    CHECK(std::abs(v->eigenvalues.back() - Complex{7, 1}) <= 1e-12);
}

TEST_CASE(
    "Full selected pairs retain physical scales left correspondence and rigorous independence",
    "[verification][N4-07][N4-06][EIG-02][EIG-03][EIG-04][EIG-05][EIG-13]") {
    core::ResourcePlan plan{64 * std::uint64_t{1024} * 1024};
    const std::array<double, 4> a{2, 0, 1, 5};
    auto original = input<double>(2, a, plan);
    Scales scales{2};
    scales.d = {3, .2};
    scales.lambda_scale = 7;
    auto full = full_spectrum(original, context(), {}, scales.view(), request(), plan);
    REQUIRE(full);
    const std::array<core::Index, 2> indices{0, 1};
    auto selected = FullSelection::make(*full, indices, plan);
    REQUIRE(selected);
    FullPairRequest const wanted{
        request(core::NumericalGeneration{43}).output(),
        PairMetricRequest{{1e-12, core::ObservationId{31}, core::ScopeId{1}},
                          {1e-12, core::ObservationId{32}, core::ScopeId{1}},
                          NotRequested{}},
        {0, 1e-14},
        true,
        true};
    auto pairs = full_pairs(*full, *selected, wanted, plan);
    if (!pairs) {
        CAPTURE(pairs.error().cause.code);
        FAIL_CHECK("left/right pairs failed");
        return;
    }
    const auto v = pairs->view();
    REQUIRE(v);
    if (!v.has_value())
        throw std::runtime_error("missing expected v");
    REQUIRE(v->right.size() == 4);
    REQUIRE(v->left.size() == 4);
    CHECK(v->triangular.empty());
    CHECK(std::abs(v->right[1]) <= 1e-12);
    CHECK(std::abs(v->right[3] - 3.0 * v->right[2]) <= 1e-12);
    CHECK(std::abs(v->left[0] + 3.0 * v->left[1]) <= 1e-12);
    CHECK(std::abs(v->left[2]) <= 1e-12);
    const double expected = std::sqrt(2026.0) / 45.0;
    for (core::Index j = 0; j < 2; ++j) {
        const auto proof = pairs->pair_proof(j);
        REQUIRE(proof);
        REQUIRE(proof->left);
        REQUIRE(proof->sensitivity);
        if (!proof->sensitivity.has_value())
            throw std::runtime_error("missing sensitivity");
        const auto sensitivity = *proof->sensitivity;
        REQUIRE(sensitivity.estimate);
        if (!sensitivity.estimate.has_value())
            throw std::runtime_error("missing estimate");
        CHECK(std::abs(*sensitivity.estimate - expected) <= 1e-12);
        CHECK_FALSE(sensitivity.condition_unresolved);
    }
    const auto independence = pairs->independence_upper();
    REQUIRE(independence);
    if (!independence.has_value())
        throw std::runtime_error("missing expected independence");
    CHECK(*independence < 1);
    CHECK_FALSE(pairs->pair_proof(-1));
    CHECK_FALSE(pairs->pair_proof(2));
}

TEST_CASE("A defective Schur packet does not pretend to have a complete eigenvector basis",
          "[negative][N4-07][EIG-05][EIG-10]") {
    core::ResourcePlan plan{64 * std::uint64_t{1024} * 1024};
    const std::array<double, 4> a{1, 0, 1, 1};
    auto original = input<double>(2, a, plan);
    Scales const scales{2};
    auto full = full_spectrum(original, context(), {}, scales.view(), request(), plan);
    REQUIRE(full);
    const std::array<core::Index, 2> indices{0, 1};
    auto selected = FullSelection::make(*full, indices, plan);
    REQUIRE(selected);
    FullPairRequest const wanted{
        request(core::NumericalGeneration{43}).output(),
        PairMetricRequest{{1e-12, core::ObservationId{31}, core::ScopeId{1}},
                          {1e-12, core::ObservationId{32}, core::ScopeId{1}},
                          NotRequested{}},
        {0, 0},
        true,
        true};
    const auto pairs = full_pairs(*full, *selected, wanted, plan);
    REQUIRE_FALSE(pairs);
    CHECK(pairs.error().cause.code == core::FailureCode::RankUnresolved);
    REQUIRE(full->basis_proof());
}

TEST_CASE("Full PG24 roots use the frozen original P and independent SciPy spectrum",
          "[verification][N4-07][N4-06][EIG-03][EIG-05][EIG-13]") {
    std::ifstream file{UVLAS_FULL_FIXTURE_PATH};
    REQUIRE(file.good());
    nlohmann::json fixture;
    file >> fixture;
    CHECK(fixture.at("source_record").at("id") == "numerical-foundation-adoption-audit");
    CHECK(fixture.at("source_record").at("sha256") ==
          "e1c9ca1f6103e0651c416bc77eea6280d434a1934970f4779e93c26124ff3ad8");
    constexpr std::size_t n = 24;
    REQUIRE(fixture.at("step_values").size() == n * n);
    std::array<double, n * n> a{};
    std::array<Complex, n * n> complex_a{};
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            a[i + j * n] = fixture.at("step_values").at(i * n + j).get<double>();
            complex_a[i + j * n] = a[i + j * n];
        }
    std::vector<Complex> oracle;
    for (const auto& pair : fixture.at("eigenvalues"))
        oracle.emplace_back(pair[0].get<double>(), pair[1].get<double>());
    REQUIRE(oracle.size() == n);
    core::ResourcePlan plan{64 * std::uint64_t{1024} * 1024};
    auto original = input<double>(24, a, plan);
    Scales const scales{n};
    const auto result = full_spectrum(original, context(), {}, scales.view(), request(), plan);
    if (!result) {
        CAPTURE(result.error().cause.code);
        FAIL_CHECK("PG24 full spectrum failed");
        return;
    }
    const auto v = result->view();
    REQUIRE(v);
    if (!v.has_value())
        throw std::runtime_error("missing expected v");
    root_matching(v->eigenvalues, oracle);
    original_schur_reference(*v, complex_a, scales, result->diagonal_mapping());
    CHECK(result->diagnostics().converged_count == 24);
}

TEST_CASE("Full phase timing is owned and disabled mode leaves numerical bytes unchanged",
          "[unit][N4-07][N4-09][EIG-11][EIG-12]") {
    core::ResourcePlan plan{64 * std::uint64_t{1024} * 1024};
    const std::array<double, 4> a{2, 0, 1, 5};
    auto original = input<double>(2, a, plan);
    Scales const scales{2};
    auto disabled = request();
    auto enabled = request();
    enabled.timing = FullTimingMode::Enabled;
    auto off = full_spectrum(original, context(), {}, scales.view(), disabled, plan);
    auto on = full_spectrum(original, context(), {}, scales.view(), enabled, plan);
    REQUIRE(off);
    REQUIRE(on);
    CHECK(off->clock_reads() == 0);
    CHECK_FALSE(off->timings());
    CHECK(on->clock_reads() == 6);
    const auto timings = on->timings();
    REQUIRE(timings);
    if (!timings.has_value())
        throw std::runtime_error("missing expected timings");
    CHECK(timings->preparation.calls == 1);
    CHECK(timings->backend.calls == 1);
    CHECK(timings->qualification.calls == 1);
    const auto ov = off->view(), iv = on->view();
    REQUIRE(ov);
    if (!ov.has_value())
        throw std::runtime_error("missing expected ov");
    REQUIRE(iv);
    if (!iv.has_value())
        throw std::runtime_error("missing expected iv");
    for (const auto& pair :
         {std::pair{ov->eigenvalues, iv->eigenvalues}, std::pair{ov->right, iv->right},
          std::pair{ov->triangular, iv->triangular}}) {
        REQUIRE(pair.first.size() == pair.second.size());
        for (std::size_t i = 0; i < pair.first.size(); ++i) {
            CHECK(std::bit_cast<std::uint64_t>(pair.first[i].real()) ==
                  std::bit_cast<std::uint64_t>(pair.second[i].real()));
            CHECK(std::bit_cast<std::uint64_t>(pair.first[i].imag()) ==
                  std::bit_cast<std::uint64_t>(pair.second[i].imag()));
        }
    }
    auto invalid = request();
    invalid.timing = static_cast<FullTimingMode>(255);
    const auto bad = full_spectrum(original, context(), {}, scales.view(), invalid, plan);
    REQUIRE_FALSE(bad);
    CHECK(bad.error().cause.code == core::FailureCode::InvalidInput);
}

TEST_CASE("Full result capacity is independent from original order capacity",
          "[negative][N4-07][N4-09][EIG-01][EIG-11]") {
    core::ResourcePlan plan{64 * std::uint64_t{1024} * 1024};
    const std::array<double, 4> a{2, 0, 0, 5};
    auto original = input<double>(2, a, plan);
    Scales const scales{2};
    for (const bool output : {false, true}) {
        auto limited = request();
        if (output)
            limited.limits.maximum_results = 1;
        else
            limited.limits.maximum_order = 1;
        core::ResourcePlan empty{0};
        const auto bad = full_spectrum(original, context(), {}, scales.view(), limited, empty);
        REQUIRE_FALSE(bad);
        CHECK(bad.error().cause.code == core::FailureCode::ResourceExceeded);
        CHECK(empty.allocated_bytes() == 0);
    }
    auto full = full_spectrum(original, context(), {}, scales.view(), request(), plan);
    REQUIRE(full);
    const std::array<core::Index, 1> index{0};
    auto selection = FullSelection::make(*full, index, plan);
    REQUIRE(selection);
    FullClusterRequest derived{request(core::NumericalGeneration{43}), {0, 0}, false};
    derived.common.limits.maximum_results = 0;
    core::ResourcePlan empty{0};
    const auto bad = full_cluster(*full, *selection, derived, empty);
    REQUIRE_FALSE(bad);
    CHECK(bad.error().cause.code == core::FailureCode::ResourceExceeded);
    CHECK(empty.allocated_bytes() == 0);
}

TEST_CASE("Physical phase ordering is exact across the complete binary64 exponent range",
          "[unit][N4-07][N4-06][EIG-04][EIG-13]") {
    using Rational = boost::multiprecision::cpp_rational;
    const double maximum = std::numeric_limits<double>::max(),
                 tiny = std::numeric_limits<double>::denorm_min();
    struct Entry {
        Entry(Complex x, double sx, Complex y, double sy) : a(x), da(sx), b(y), db(sy) {}
        Complex a;
        double da;
        Complex b;
        double db;
    };
    const std::array entries{Entry{{maximum, maximum}, maximum, {maximum, 0}, maximum},
                             Entry{{maximum, tiny}, 1, {maximum, 0}, 1},
                             Entry{{0, -0.0}, 3, {-0.0, 0}, .2},
                             Entry{{3, 0}, 3, {.2, 0}, .2},
                             Entry{{1, 0}, 2, {1, 0}, 1},
                             Entry{{tiny, tiny}, maximum, {tiny, 0}, maximum},
                             Entry{{1, 1}, 3, {1, 1}, 3}};
    auto square = [](Complex z) {
        return Rational(z.real()) * Rational(z.real()) + Rational(z.imag()) * Rational(z.imag());
    };
    for (const auto& e : entries)
        for (const auto side : {PairSide::Right, PairSide::Left}) {
            const Rational a = side == PairSide::Right
                                   ? square(e.a) / (Rational(e.da) * Rational(e.da))
                                   : square(e.a) * Rational(e.da) * Rational(e.da);
            const Rational b = side == PairSide::Right
                                   ? square(e.b) / (Rational(e.db) * Rational(e.db))
                                   : square(e.b) * Rational(e.db) * Rational(e.db);
            CHECK(uvlas::numerics::eigen::detail::full_phase_coordinate_larger(e.a, e.da, e.b, e.db,
                                                                               side) == (a > b));
            CHECK(uvlas::numerics::eigen::detail::full_phase_coordinate_larger(e.b, e.db, e.a, e.da,
                                                                               side) == (b > a));
        }
}
TEST_CASE("Empty clusters retain whether a left basis was requested",
          "[unit][N4-07][N4-08][EIG-01][EIG-05][EIG-10]") {
    core::ResourcePlan plan{64 * std::uint64_t{1024} * 1024};
    const std::array<double, 4> a{2, 0, 0, 5};
    auto original = input<double>(2, a, plan);
    Scales const scales{2};
    auto full = full_spectrum(original, context(), {}, scales.view(), request(), plan);
    REQUIRE(full);
    const std::span<const core::Index> indices{};
    auto selection = FullSelection::make(*full, indices, plan);
    REQUIRE(selection);
    for (const bool left : {false, true}) {
        FullClusterRequest const req{request(core::NumericalGeneration{43}), {0, 0}, left};
        auto result = full_cluster(*full, *selection, req, plan);
        REQUIRE(result);
        const auto view = result->view();
        REQUIRE(view);
        if (!view.has_value())
            throw std::runtime_error("missing expected view");
        CHECK(view->left_requested == left);
        CHECK((result->left_basis_proof() != nullptr) == left);
        CHECK(result->diagnostics().forward_applications == 0);
        CHECK(result->diagnostics().adjoint_applications == 0);
    }
    const auto metadata = full->scales();
    REQUIRE(metadata);
    if (!metadata.has_value())
        throw std::runtime_error("missing expected metadata");
    CHECK(metadata->eigenvalue_scale == 1);
    CHECK(metadata->components[0] == 1);
    CHECK(full->original_coefficients().size() == 4);
    const auto limits = full->limits();
    REQUIRE(limits);
    if (!limits)
        throw std::runtime_error("missing limits");
    CHECK(limits->maximum_results == 96);
    REQUIRE(full->backend_admission());
    CHECK(full->schur_iteration_limit() == 60);
    CHECK(full->source_backend() == FullBackend::RealSchurComplexified);
    CHECK_FALSE(full->association_policy());
}

TEST_CASE("A noncontiguous Jordan cluster needs a genuine right unitary exchange",
          "[verification][N4-07][EIG-05][EIG-07][EIG-13]") {
    core::ResourcePlan plan{64 * std::uint64_t{1024} * 1024};
    const std::array<double, 9> a{1, 0, 0, 1, 3, 0, 2, 1, 1};
    auto original = input<double>(3, a, plan);
    Scales const scales{3};
    auto full = full_spectrum(original, context(), {}, scales.view(), request(), plan);
    REQUIRE(full);
    const std::array<core::Index, 2> indices{0, 1};
    auto selection = FullSelection::make(*full, indices, plan);
    REQUIRE(selection);
    FullClusterRequest const wanted{request(core::NumericalGeneration{43}), {1e-12, 1e-12}, true};
    auto cluster = full_cluster(*full, *selection, wanted, plan);
    REQUIRE(cluster);
    const auto view = cluster->view();
    REQUIRE(view);
    if (!view.has_value())
        throw std::runtime_error("missing expected view");
    const std::array<double, 9> rp{1, 0, 0, 0, .2, -.4, 0, -.4, .8},
        lp{.8, -.4, 0, -.4, .2, 0, 0, 0, 1};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            Complex r{}, l{};
            for (std::size_t k = 0; k < 2; ++k) {
                r += view->right[i + 3 * k] * std::conj(view->right[j + 3 * k]);
                l += view->left[i + 3 * k] * std::conj(view->left[j + 3 * k]);
            }
            CHECK(std::abs(r - rp[i + 3 * j]) <= 1e-12);
            CHECK(std::abs(l - lp[i + 3 * j]) <= 1e-12);
        }
}
TEST_CASE("Pair phases are chosen from actual physical complex vectors after nonpower scaling",
          "[verification][N4-07][EIG-02][EIG-03][EIG-04][EIG-13]") {
    core::ResourcePlan plan{64 * std::uint64_t{1024} * 1024};
    const std::array<Complex, 4> a{Complex{1, 0}, Complex{0, 0}, Complex{0, 30}, Complex{3, 0}};
    auto original = input<Complex>(2, a, plan);
    Scales scales{2};
    scales.d = {3, .2};
    auto full = full_spectrum(original, context(), {}, scales.view(), request(), plan);
    REQUIRE(full);
    const std::array<core::Index, 2> indices{0, 1};
    auto selection = FullSelection::make(*full, indices, plan);
    REQUIRE(selection);
    FullPairRequest const wanted{request(core::NumericalGeneration{43}).output(),
                                 {{1e-12, core::ObservationId{1}, core::ScopeId{1}},
                                  {1e-12, core::ObservationId{2}, core::ScopeId{1}},
                                  NotRequested{}},
                                 {0, 1e-14},
                                 true,
                                 true};
    auto pairs = full_pairs(*full, *selection, wanted, plan);
    REQUIRE(pairs);
    const auto view = pairs->view();
    REQUIRE(view);
    if (!view.has_value())
        throw std::runtime_error("missing expected view");
    using R = boost::multiprecision::cpp_rational;
    for (const bool left : {false, true})
        for (std::size_t col = 0; col < 2; ++col) {
            const auto values = left ? view->left : view->right;
            std::size_t pivot = 0;
            R best = 0;
            for (std::size_t i = 0; i < 2; ++i) {
                const auto z = values[i + col * 2];
                R norm = R(z.real()) * R(z.real()) + R(z.imag()) * R(z.imag());
                const R d = R(scales.d[i]) * R(scales.d[i]);
                if (left)
                    norm *= d;
                else
                    norm /= d;
                if (norm > best) {
                    best = norm;
                    pivot = i;
                }
            }
            CHECK(values[pivot + col * 2].real() >= 0);
            CHECK(values[pivot + col * 2].imag() == 0);
        }
}

TEST_CASE("A nonzero scaled coefficient is never silently replaced by an underflowed zero",
          "[negative][N4-07][N4-06][EIG-02][EIG-03]") {
    core::ResourcePlan plan{std::uint64_t{1024} * 1024};
    const std::array<double, 1> a{std::numeric_limits<double>::denorm_min()};
    auto original = input<double>(1, a, plan);
    Scales scales{1};
    scales.lambda_scale = 2;
    const auto result = full_spectrum(original, context(), {}, scales.view(), request(), plan);
    REQUIRE_FALSE(result);
    CHECK(result.error().cause.code == core::FailureCode::PrecisionUnresolved);
}
