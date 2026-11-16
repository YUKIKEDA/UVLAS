#include "eigen_arnoldi_cluster.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <bit>
#include <cmath>
#include <complex>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <uvlas/numerics/eigen_arnoldi.hpp>
#include <uvlas/numerics/eigen_full.hpp>
#include <vector>

namespace {
using namespace uvlas;
namespace eig = numerics::eigen;
using Complex = std::complex<double>;
// Catch2's assertion macros do not establish optional engagement for static
// analysis. Borrow lvalues, or copy temporaries, only after a real C++ check.
template <class T>
[[nodiscard]] const T& checked(const std::optional<T>& value) {
    if (!value) {
        throw std::runtime_error("required optional test value is absent");
    }
    return *value;
}
template <class T>
[[nodiscard]] T checked(std::optional<T>&& value) {
    if (!value) {
        throw std::runtime_error("required optional test value is absent");
    }
    return std::move(*value);
}

[[nodiscard]] eig::Context context() {
    return {{core::ModelFingerprint{}, core::ModelGeneration{1}, core::GeometryGeneration{2},
             core::ChartGeneration{3}, core::SparsityGeneration{4}, core::NumericalGeneration{5},
             core::ScaleGeneration{6}, core::PrecisionGeneration{7}, core::Precision::Binary64,
             core::WakePolicyGeneration{8}, core::TimeIntervalGeneration{9},
             core::BackendGeneration{10}},
            core::TrialId{11},
            core::ScopeId{12},
            core::ScopeId{13},
            .25};
}
struct NoConsumer final : eig::ArnoldiConsumerEvaluator {
    [[nodiscard]] numerics::SupportedPrecision supported_precision() const noexcept override {
        return numerics::SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const eig::ArnoldiCandidateView&, const eig::ArnoldiCandidateView&,
           const numerics::linear::ConsumerCheck&, core::CertificateContext) const override {
        throw std::runtime_error("an empty consumer list must never invoke its evaluator");
    }
};
struct Rotation final : eig::ActionProvider<double> {
    core::Index n = 8;
    mutable std::uint64_t calls = 0;
    unsigned fault = 0;
    std::uint64_t fault_at = 1;
    mutable std::uint64_t generation = 21;
    std::function<void()> on_first_action;
    double uncertainty = 0;
    [[nodiscard]] core::Index order() const noexcept override {
        return n;
    }
    [[nodiscard]] eig::ActionStamp stamp() const noexcept override {
        return {context(), core::NumericalGeneration{generation}, {}};
    }
    [[nodiscard]] numerics::SupportedPrecision
    supported_precision(eig::ActionKind) const noexcept override {
        return numerics::SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        return {};
    }
    [[nodiscard]] core::Result<void> apply(eig::ActionKind kind, eig::WordView<double> x,
                                           std::span<double> hi, std::span<double> lo,
                                           std::span<double> error) const override {
        ++calls;
        if (calls == 1 && on_first_action) {
            on_first_action();
        }
        for (core::Index i = 0; i < n; ++i) {
            const auto row = static_cast<std::size_t>(i);
            const auto column = static_cast<std::size_t>(i == 0 ? 1 : (i == 1 ? 0 : i));
            double coefficient = i == 0 ? -5 : (i == 1 ? 5 : static_cast<double>(i) / 10);
            if (kind == eig::ActionKind::Adjoint && i < 2) {
                coefficient = -coefficient;
            }
            const auto high = x.high()[column];
            const auto low = x.precision() == core::Precision::TwoWord ? x.low()[column] : 0;
            hi[row] = coefficient * high;
            lo[row] = std::fma(coefficient, high, -hi[row]) + coefficient * low;
            // Two products are exact here except possible subnormal rounding.
            error[row] = uncertainty + 4 * std::numeric_limits<double>::denorm_min();
        }
        if (calls == fault_at) {
            if (fault == 1) {
                throw std::bad_alloc{};
            }
            if (fault == 2) {
                throw std::runtime_error("known producer failure");
            }
            if (fault == 3) {
                core::FailureDetails details;
                details.resources = core::ResourceUsage{core::ResourceKind::Bytes, 900, 800};
                return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved,
                                                     core::Operation::CertifyError, details});
            }
            if (fault == 4) {
                error.back() = std::numeric_limits<double>::quiet_NaN();
            }
            if (fault == 5) {
                error.back() = -1;
            }
            if (fault == 6) {
                ++generation;
            }
        }
        return {};
    }
};
[[nodiscard]] eig::ArnoldiPairRequest accuracy() {
    return {{{1e-12, core::ObservationId{31}, core::ScopeId{32}},
             {1e-12, core::ObservationId{33}, core::ScopeId{32}},
             eig::NotRequested{}}};
}
template <class Scalar>
[[nodiscard]] eig::ArnoldiRequest<Scalar> request(std::span<const Scalar> seed, core::Index nev,
                                                  core::Index ncv, const NoConsumer& evaluator) {
    const auto words = eig::WordView<Scalar>::make(seed, {}, core::Precision::Binary64);
    if (!words) {
        throw std::runtime_error("invalid independent fixture seed");
    }
    return {nev,
            ncv,
            *words,
            1e-14,
            100,
            eig::ArnoldiLargestMagnitude{},
            accuracy(),
            eig::ArnoldiSides::RightOnly,
            {1e-12, 1e-12},
            numerics::linear::RetainedPrecision{core::NumericalGeneration{40}},
            {},
            {},
            evaluator,
            eig::ArnoldiTimingMode::Disabled};
}
template <class Scalar>
[[nodiscard]] eig::Result<eig::ArnoldiSpectrum>
extract(const eig::ActionProvider<Scalar>& provider, const eig::ArnoldiScales& fixed,
        const eig::ArnoldiRequest<Scalar>& wanted, const eig::ArnoldiLimits& capacities,
        core::ResourcePlan& workspace_plan, core::ResourcePlan& result_plan) {
    auto workspace = eig::ArnoldiWorkspace::create(provider.order(), capacities.maximum_returned,
                                                   workspace_plan);
    if (!workspace) {
        return std::unexpected(
            eig::complete_failure(workspace.error(), fixed.context, wanted.state));
    }
    return eig::arnoldi_spectrum(provider, fixed, wanted, capacities, *workspace, result_plan);
}
[[nodiscard]] eig::ArnoldiLimits limits(core::Index n, core::Index k, core::Index ncv) {
    return {
        n, k, ncv, 100, {10000, 10000, 10000, 4}, numerics::linear::EstimateAdmission{100000000}};
}
} // namespace

TEST_CASE("Arnoldi qualifies the original real conjugate pair",
          "[N4-07][EIG-01][EIG-08][EIG-13][verification]") {
    const Rotation op;
    const std::array<double, 8> seed{1, 2, 3, 4, 5, 6, 7, 8};
    const std::array<double, 8> scales{1, 1, 1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 8> units{};
    const NoConsumer consumer;
    auto wanted = request<double>(seed, 2, 8, consumer);
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    core::ResourcePlan workspace(100000000), result(100000000);
    auto extracted = extract(op, fixed, wanted, limits(8, 3, 8), workspace, result);
    INFO("failure code " << (extracted ? 255U
                                       : static_cast<unsigned>(extracted.error().cause.code)));
    REQUIRE(extracted.has_value());
    const auto values = extracted->values();
    REQUIRE(values.has_value());
    REQUIRE(checked(values).eigenvalues.shape.rows == 2);
    CHECK(std::abs(checked(values).eigenvalues.values[0] - Complex{0, -5}) <= 1e-12);
    CHECK(std::abs(checked(values).eigenvalues.values[1] - Complex{0, 5}) <= 1e-12);
    const auto* pairs = std::get_if<eig::ArnoldiPairView>(&checked(values).right);
    REQUIRE(pairs != nullptr);
    for (std::size_t j = 0; j < 2; ++j) {
        const auto lambda = checked(values).eigenvalues.values[j];
        const auto phi = pairs->vectors.values.subspan(8 * j, 8);
        CHECK(std::abs(-5.0 * phi[1] - lambda * phi[0]) <= 1e-12);
        CHECK(std::abs(5.0 * phi[0] - lambda * phi[1]) <= 1e-12);
        for (std::size_t i = 2; i < 8; ++i) {
            CHECK(std::abs((static_cast<double>(i) / 10 - lambda) * phi[i]) <= 1e-12);
        }
        REQUIRE(extracted->pair_metric(eig::ArnoldiSide::Right, static_cast<core::Index>(j))
                    .has_value());
    }
    const auto counts = extracted->diagnostics();
    REQUIRE(counts.has_value());
    CHECK(checked(counts).attempts.forward_applications == op.calls);
    CHECK(std::abs(checked(counts).initial_norm - std::sqrt(204.0)) <= 1e-12);
    // Independent Python hashlib oracle: domain, little-endian u64 tags/length,
    // then the exact IEEE binary64 bits of 1..8. No product hash helper is reused.
    const std::array<std::byte, 32> expected_seed_hash{
        std::byte{0xa3}, std::byte{0x9b}, std::byte{0x3f}, std::byte{0xff}, std::byte{0xc9},
        std::byte{0x41}, std::byte{0x71}, std::byte{0xd5}, std::byte{0xf2}, std::byte{0xe3},
        std::byte{0xc5}, std::byte{0x06}, std::byte{0x91}, std::byte{0xfa}, std::byte{0xb8},
        std::byte{0x94}, std::byte{0xe7}, std::byte{0x6b}, std::byte{0xe9}, std::byte{0xa7},
        std::byte{0x21}, std::byte{0x2b}, std::byte{0x69}, std::byte{0x09}, std::byte{0x78},
        std::byte{0xe4}, std::byte{0x09}, std::byte{0xed}, std::byte{0x7e}, std::byte{0x26},
        std::byte{0xc6}, std::byte{0xfd}};
    CHECK(checked(counts).initial_sha256 == expected_seed_hash);
}

TEST_CASE("Arnoldi completes a conjugate boundary at the maximum legal nev",
          "[N4-07][EIG-08][negative][verification]") {
    Rotation op;
    op.n = 3;
    const std::array<double, 3> seed{1, 2, 3}, scales{1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 3> units{};
    const NoConsumer consumer;
    auto wanted = request<double>(seed, 1, 3, consumer);
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    core::ResourcePlan workspace(1000000), result(1000000);
    auto extracted = extract(op, fixed, wanted, limits(3, 2, 3), workspace, result);
    REQUIRE(extracted.has_value());
    const auto values = extracted->values();
    REQUIRE(values.has_value());
    CHECK(checked(values).eigenvalues.values.size() == 2);
    CHECK(std::abs(checked(values).eigenvalues.values[0] - Complex{0, -5}) <= 1e-12);
    CHECK(std::abs(checked(values).eigenvalues.values[1] - Complex{0, 5}) <= 1e-12);
    CHECK_FALSE(extracted->pair_metric(eig::ArnoldiSide::Right, 2).has_value());
    core::ResourcePlan workspace2(1000000), result2(1000000);
    const auto calls = op.calls;
    auto capacity = extract(op, fixed, wanted, limits(3, 1, 3), workspace2, result2);
    REQUIRE_FALSE(capacity.has_value());
    CHECK(capacity.error().cause.code == core::FailureCode::ResourceExceeded);
    CHECK(op.calls == calls);
}

namespace {
struct Jordan final : eig::ActionProvider<double> {
    double second = 5;
    [[nodiscard]] core::Index order() const noexcept override {
        return 8;
    }
    [[nodiscard]] eig::ActionStamp stamp() const noexcept override {
        return {context(), core::NumericalGeneration{21}, {}};
    }
    [[nodiscard]] numerics::SupportedPrecision
    supported_precision(eig::ActionKind) const noexcept override {
        return numerics::SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        return {};
    }
    [[nodiscard]] core::Result<void> apply(eig::ActionKind kind, eig::WordView<double> x,
                                           std::span<double> hi, std::span<double> lo,
                                           std::span<double> error) const override {
        for (std::size_t i = 0; i < 8; ++i) {
            const auto diagonal = i == 0 ? 5.0 : (i == 1 ? second : static_cast<double>(i) / 10);
            const auto trailing = x.low().empty() ? 0 : x.low()[i];
            hi[i] = diagonal * x.high()[i];
            lo[i] = diagonal * trailing;
            error[i] = 8 * std::numeric_limits<double>::epsilon() * diagonal *
                           (std::abs(x.high()[i]) + std::abs(trailing)) +
                       16 * std::numeric_limits<double>::denorm_min();
        }
        const auto row = kind == eig::ActionKind::Adjoint ? 1U : 0U;
        const auto column = 1U - row;
        hi[row] += x.high()[column];
        if (!x.low().empty()) {
            lo[row] += x.low()[column];
        }
        error[row] += 8 * std::numeric_limits<double>::epsilon() *
                      (std::abs(x.high()[column]) + std::abs(hi[row]));
        return {};
    }
};
} // namespace
TEST_CASE("Arnoldi lifts a defective selected space from the active Krylov basis",
          "[N4-07][EIG-01][EIG-08][EIG-13][verification]") {
    const Jordan op;
    const std::array<double, 8> seed{1, 2, 3, 4, 5, 6, 7, 8}, scales{1, 1, 1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 8> units{};
    const NoConsumer consumer;
    auto wanted = request<double>(seed, 2, 8, consumer);
    wanted.qualification =
        eig::ArnoldiSubspaceRequest{{{1e-12, core::ObservationId{31}, core::ScopeId{32}},
                                     {1e-12, core::ObservationId{33}, core::ScopeId{32}},
                                     eig::NotRequested{}}};
    wanted.association = {1e-6, 1e-6};
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    core::ResourcePlan workspace(10000000), result(10000000);
    auto extracted = extract(op, fixed, wanted, limits(8, 3, 8), workspace, result);
    INFO("failure code " << (extracted ? 255U
                                       : static_cast<unsigned>(extracted.error().cause.code)));
    REQUIRE(extracted.has_value());
    const auto values = extracted->values();
    REQUIRE(values.has_value());
    const auto* space = std::get_if<eig::ArnoldiSubspaceView>(&checked(values).right);
    REQUIRE(space != nullptr);
    REQUIRE(space->basis.shape.columns == 2);
    // The exact invariant projector is diag(1,1,0,...,0), independent of
    // the nonunique Schur columns of the defective repeated eigenvalue.
    for (std::size_t i = 0; i < 8; ++i) {
        for (std::size_t j = 0; j < 8; ++j) {
            Complex projected{};
            for (std::size_t k = 0; k < 2; ++k) {
                projected +=
                    space->basis.values[k * 8 + i] * std::conj(space->basis.values[k * 8 + j]);
            }
            CHECK(std::abs(projected - Complex{i == j && i < 2 ? 1.0 : 0.0, 0}) <= 1e-12);
        }
    }
    REQUIRE(extracted->subspace_metric(eig::ArnoldiSide::Right).has_value());
    CHECK(extracted->sensitivities().empty());
}

namespace {
template <class Scalar>
struct LiteralAction final : eig::ActionProvider<Scalar> {
    core::Index n;
    std::vector<Scalar> coefficients;
    mutable std::uint64_t calls = 0;
    LiteralAction(core::Index order, std::vector<Scalar> matrix)
        : n(order), coefficients(std::move(matrix)) {}
    [[nodiscard]] core::Index order() const noexcept override {
        return n;
    }
    [[nodiscard]] eig::ActionStamp stamp() const noexcept override {
        return {context(), core::NumericalGeneration{21}, {}};
    }
    [[nodiscard]] numerics::SupportedPrecision
    supported_precision(eig::ActionKind) const noexcept override {
        return numerics::SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        return {};
    }
    [[nodiscard]] core::Result<void> apply(eig::ActionKind kind, eig::WordView<Scalar> input,
                                           std::span<Scalar> hi, std::span<Scalar> lo,
                                           std::span<double> error) const override {
        ++calls;
        for (core::Index i = 0; i < n; ++i) {
            Scalar sum{};
            double magnitude = 0;
            std::size_t terms = 0;
            for (core::Index j = 0; j < n; ++j) {
                auto coefficient = coefficients[static_cast<std::size_t>(
                    kind == eig::ActionKind::Adjoint ? j * n + i : i * n + j)];
                if constexpr (std::same_as<Scalar, Complex>) {
                    if (kind == eig::ActionKind::Adjoint) {
                        coefficient = std::conj(coefficient);
                    }
                }
                if (coefficient != Scalar{}) {
                    ++terms;
                }
                const auto value = input.high()[static_cast<std::size_t>(j)];
                const auto tail =
                    input.low().empty() ? Scalar{} : input.low()[static_cast<std::size_t>(j)];
                sum += coefficient * value;
                sum += coefficient * tail;
                magnitude += std::abs(coefficient) * (std::abs(value) + std::abs(tail));
            }
            hi[static_cast<std::size_t>(i)] = sum;
            lo[static_cast<std::size_t>(i)] = {};
            // A conservative gamma for both planes and complex multiply/adds.
            // Adding an exact zero product adds no rounding. Thus the count is
            // the actual nonzero terms, not matrix width. For the order-96
            // diagonal each row has one multiplication; the earlier width-96
            // bound refused qualification despite accurate candidates. The
            // product enclosure is independently checked below, with budget unchanged.
            error[static_cast<std::size_t>(i)] = std::nextafter(
                16 * static_cast<double>(terms) * std::numeric_limits<double>::epsilon() *
                        magnitude +
                    16 * static_cast<double>(n) * std::numeric_limits<double>::denorm_min(),
                std::numeric_limits<double>::infinity());
        }
        return {};
    }
};
} // namespace
TEST_CASE("Arnoldi keeps a complex nonconjugating low-rank original action",
          "[N4-07][EIG-01][EIG-08][EIG-13][verification]") {
    std::vector<Complex> matrix(36);
    const std::array<Complex, 6> diagonal{{{7, 1}, {3, 0}, {2, 0}, {1, 0}, {.5, 0}, {.25, 0}}};
    for (std::size_t i = 0; i < 6; ++i) {
        matrix[6 * i + i] = diagonal[i];
    }
    // diag(8,3,2,1,.5,.25) + (i e0) ((1+i)e0+e1)^T.
    // A conjugating update would give 9+i instead of the literal 7+i.
    matrix[1] = {0, 1};
    const LiteralAction<Complex> op(6, std::move(matrix));
    const std::array<Complex, 6> seed{{{1, 1}, {2, -1}, {3, 1}, {4, -1}, {5, 1}, {6, -1}}};
    const std::array<double, 6> scales{1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 6> units{};
    const NoConsumer consumer;
    auto wanted = request<Complex>(seed, 2, 6, consumer);
    wanted.sides = eig::ArnoldiSides::Both;
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    core::ResourcePlan workspace(10000000), result(10000000);
    auto extracted = extract(op, fixed, wanted, limits(6, 2, 6), workspace, result);
    INFO("failure code " << (extracted ? 255U
                                       : static_cast<unsigned>(extracted.error().cause.code)));
    REQUIRE(extracted.has_value());
    const auto values = extracted->values();
    REQUIRE(values.has_value());
    CHECK(std::abs(checked(values).eigenvalues.values[0] - Complex{3, 0}) <= 1e-12);
    CHECK(std::abs(checked(values).eigenvalues.values[1] - Complex{7, 1}) <= 1e-12);
    REQUIRE(extracted->pair_metric(eig::ArnoldiSide::Left, 0).has_value());
    REQUIRE(extracted->pair_metric(eig::ArnoldiSide::Left, 1).has_value());
}
TEST_CASE("Arnoldi honors explicit 24 and 64 on order 96",
          "[N4-07][EIG-01][EIG-08][EIG-13][verification]") {
    std::vector<double> matrix(std::size_t{96} * 96), seed(96), scales(96, 1);
    std::vector<numerics::linear::UnitDimensions> units(96);
    for (std::size_t i = 0; i < 96; ++i) {
        matrix[96 * i + i] = 96 - static_cast<double>(i);
        seed[i] = static_cast<double>(i + 1);
    }
    const LiteralAction<double> op(96, std::move(matrix));
    const NoConsumer consumer;
    auto wanted = request<double>(seed, 24, 64, consumer);
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    core::ResourcePlan workspace(100000000), result(100000000);
    auto extracted = extract(op, fixed, wanted, limits(96, 25, 64), workspace, result);
    INFO("failure code " << (extracted ? 255U
                                       : static_cast<unsigned>(extracted.error().cause.code)));
    REQUIRE(extracted.has_value());
    const auto values = extracted->values();
    REQUIRE(values.has_value());
    REQUIRE(checked(values).eigenvalues.values.size() == 24);
    for (std::size_t i = 0; i < 24; ++i) {
        CHECK(std::abs(checked(values).eigenvalues.values[i] -
                       Complex{73 + static_cast<double>(i), 0}) <= 1e-12);
    }
}
TEST_CASE("Arnoldi extracts the stored original PG24 action with explicit six and twenty",
          "[N4-07][EIG-01][EIG-08][EIG-13][verification]") {
    std::ifstream file(std::string(UVLAS_ARNOLDI_REFERENCE_ROOT) + "/fixture.json");
    REQUIRE(file.good());
    const auto fixture = nlohmann::json::parse(file);
    const auto matrix = fixture.at("step_values").get<std::vector<double>>();
    REQUIRE(matrix.size() == std::size_t{24} * 24);
    const LiteralAction<double> op(24, matrix);
    std::vector<double> seed(24), scales(24, 1);
    std::vector<numerics::linear::UnitDimensions> units(24);
    for (std::size_t i = 0; i < 24; ++i) {
        seed[i] = static_cast<double>(i + 1);
    }
    const NoConsumer consumer;
    auto wanted = request<double>(seed, 6, 20, consumer);
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    core::ResourcePlan workspace(100000000), result(100000000);
    auto extracted = extract(op, fixed, wanted, limits(24, 7, 20), workspace, result);
    INFO("failure code " << (extracted ? 255U
                                       : static_cast<unsigned>(extracted.error().cause.code)));
    REQUIRE(extracted.has_value());
    const auto values = extracted->values();
    REQUIRE(values.has_value());
    REQUIRE(checked(values).eigenvalues.values.size() == 6);
    std::vector<Complex> reference;
    for (const auto& root : fixture.at("eigenvalues")) {
        reference.emplace_back(root.at(0).get<double>(), root.at(1).get<double>());
    }
    std::sort(reference.begin(), reference.end(),
              [](auto a, auto b) { return std::abs(a) > std::abs(b); });
    reference.resize(6);
    std::array<bool, 6> used{};
    for (const auto lambda : checked(values).eigenvalues.values) {
        std::size_t best = 6;
        double distance = 1e-12;
        for (std::size_t j = 0; j < 6; ++j) {
            if (!used[j] && std::abs(lambda - reference[j]) <= distance) {
                best = j;
                distance = std::abs(lambda - reference[j]);
            }
        }
        REQUIRE(best < 6);
        used[best] = true;
    }
    // The independent frozen SciPy roots above remain the oracle. In addition,
    // compare the product complete method and this explicit Arnoldi subset.
    std::vector<double> packed(std::size_t{24} * 24);
    for (std::size_t j = 0; j < 24; ++j) {
        for (std::size_t i = 0; i < 24; ++i) {
            packed[j * 24 + i] = matrix[i * 24 + j];
        }
    }
    std::vector<numerics::linear::UnitDimensions> coefficient_units(std::size_t{24} * 24);
    core::ResourcePlan full_plan(100000000);
    auto dense =
        numerics::DenseReal::copy({24, 24, 24}, packed, full_plan, core::Operation::ExtractEigen);
    REQUIRE(dense.has_value());
    const auto original =
        numerics::PrecisionBuffer<double>::binary64(std::move(*dense), context().cache.matrix);
    REQUIRE(original.has_value());
    const eig::FullRequest full_request{
        {{1e-12, core::ObservationId{300}, core::ScopeId{301}},
         {1e-12, core::ObservationId{302}, core::ScopeId{301}},
         eig::NotRequested{}},
        numerics::linear::RetainedPrecision{core::NumericalGeneration{42}},
        {24, 24, {10000, 0, 0, 1}, numerics::linear::EstimateAdmission{100000000}},
        {},
        nullptr};
    const auto complete = eig::full_spectrum(
        *original, context(), {}, eig::FullScales{scales, 1, units, {}, coefficient_units},
        full_request, full_plan);
    REQUIRE(complete.has_value());
    const auto complete_view = complete->view();
    REQUIRE(complete_view.has_value());
    std::array<bool, 24> full_used{};
    for (const auto lambda : checked(values).eigenvalues.values) {
        std::size_t match = 24;
        for (std::size_t j = 0; j < 24; ++j) {
            if (!full_used[j] &&
                std::abs(lambda - checked(complete_view).eigenvalues[j]) <= 1e-12) {
                match = j;
                break;
            }
        }
        REQUIRE(match < 24);
        full_used[match] = true;
    }
    core::ResourcePlan other_workspace(100000000), other_result(100000000);
    auto unclipped = request<double>(seed, 24, 64, consumer);
    const auto rejected =
        extract(op, fixed, unclipped, limits(24, 25, 64), other_workspace, other_result);
    REQUIRE_FALSE(rejected.has_value());
    CHECK(rejected.error().cause.code == core::FailureCode::InvalidInput);
}

TEST_CASE("Literal diagonal producer bounds contain exact rational products",
          "[N4-07][EIG-03][verification]") {
    std::vector<double> matrix(std::size_t{96} * 96), input(96), hi(96), lo(96), error(96);
    for (std::size_t i = 0; i < 96; ++i) {
        matrix[96 * i + i] = 96 - static_cast<double>(i);
        input[i] = std::ldexp(1 + static_cast<double>(i) / 97, -10);
    }
    const LiteralAction<double> provider(96, std::move(matrix));
    const auto words = eig::WordView<double>::make(input, {}, core::Precision::Binary64);
    REQUIRE(words.has_value());
    REQUIRE(provider.apply(eig::ActionKind::Forward, *words, hi, lo, error).has_value());
    using Rational = boost::multiprecision::cpp_rational;
    for (std::size_t i = 0; i < 96; ++i) {
        const Rational exact = Rational(96 - i) * Rational(input[i]);
        const Rational reported = Rational(hi[i]) + Rational(lo[i]);
        CHECK(abs(exact - reported) <= Rational(error[i]));
    }
}

namespace {
struct Consumer final : eig::ArnoldiConsumerEvaluator {
    mutable std::size_t calls = 0;
    mutable bool exact_shape = true;
    bool refuse_second = false, refuse_projection = false;
    double fixed_second_upper = 0;
    const Rotation* mutate_stamp = nullptr;
    [[nodiscard]] numerics::SupportedPrecision supported_precision() const noexcept override {
        return numerics::SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const eig::ArnoldiCandidateView& retained, const eig::ArnoldiCandidateView& candidate,
           const numerics::linear::ConsumerCheck& check,
           core::CertificateContext expected) const override {
        ++calls;
        exact_shape = exact_shape && candidate.eigenvalues.values.size() == 2;
        if (const auto* pair = std::get_if<eig::ArnoldiPairView>(&candidate.right)) {
            exact_shape = exact_shape && pair->vectors.values.size() == 16;
            for (auto value : pair->vectors.values) {
                exact_shape =
                    exact_shape && std::isfinite(value.real()) && std::isfinite(value.imag());
            }
        }
        const bool refused = (refuse_second && calls % 2 == 0) ||
                             (refuse_projection && candidate.generation != retained.generation);
        if (mutate_stamp) {
            ++mutate_stamp->generation;
        }
        const auto upper = fixed_second_upper > 0 && calls % 2 == 0
                               ? fixed_second_upper
                               : (refused ? 2 * check.budget : 0);
        return core::ErrorCertificate::assess(expected, check.budget, upper);
    }
};
} // namespace
TEST_CASE("Arnoldi owns retained and projected consumer qualifications",
          "[N4-07][EIG-10][EIG-12][negative][verification]") {
    const Rotation op;
    const std::array<double, 8> seed{1, 2, 3, 4, 5, 6, 7, 8}, scales{1, 1, 1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 8> units{};
    const std::array<numerics::linear::ConsumerCheck, 2> checks{
        {{numerics::linear::CheckKind::Observation, core::ObservationId{100}, core::ScopeId{101},
          1e-12},
         {numerics::linear::CheckKind::Reaction, core::ObservationId{102}, core::ScopeId{103},
          1e-12}}};
    const NoConsumer unused;
    auto base = request<double>(seed, 2, 8, unused);
    Consumer evaluator;
    const auto words = eig::WordView<double>::make(seed, {}, core::Precision::Binary64);
    REQUIRE(words.has_value());
    auto wanted = eig::ArnoldiRequest<double>{
        base.nev,
        base.ncv,
        *words,
        base.backend_tolerance,
        base.backend_iterations,
        base.selection,
        base.qualification,
        base.sides,
        base.association,
        numerics::linear::Binary64Projection{core::NumericalGeneration{40},
                                             core::NumericalGeneration{41}},
        base.state,
        checks,
        evaluator,
        base.timing};
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    core::ResourcePlan workspace(10000000), result(10000000);
    auto extracted = extract(op, fixed, wanted, limits(8, 3, 8), workspace, result);
    REQUIRE(extracted.has_value());
    CHECK(evaluator.calls == 4);
    CHECK(evaluator.exact_shape);
    REQUIRE(extracted->values().has_value());
    REQUIRE(extracted->retained_values().has_value());
    CHECK(checked(extracted->values()).generation == core::NumericalGeneration{41});
    CHECK(checked(extracted->retained_values()).generation == core::NumericalGeneration{40});
    REQUIRE(extracted->consumer_certificate(1).has_value());
    CHECK(checked(extracted->consumer_certificate(1)).context().numerical_generation ==
          core::NumericalGeneration{41});
    evaluator.calls = 0;
    evaluator.refuse_projection = true;
    core::ResourcePlan workspace2(10000000), result2(10000000);
    const auto refused = extract(op, fixed, wanted, limits(8, 3, 8), workspace2, result2);
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error().cause.code == core::FailureCode::PrecisionUnresolved);
    CHECK(evaluator.calls == 3);
    CHECK(checked(extracted->values()).generation == core::NumericalGeneration{41});
    auto moved = std::move(*extracted);
    CHECK_FALSE(extracted->active());
    CHECK_FALSE(extracted->values().has_value());
    CHECK_FALSE(extracted->metadata().has_value());
    CHECK(extracted->sensitivities().empty());
    CHECK(moved.active());
}
TEST_CASE("Arnoldi preserves first provider failures and global action quotas",
          "[N4-07][EIG-07][EIG-11][EIG-12][negative][verification]") {
    const std::array<double, 8> seed{1, 2, 3, 4, 5, 6, 7, 8}, scales{1, 1, 1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 8> units{};
    const NoConsumer consumer;
    auto wanted = request<double>(seed, 2, 8, consumer);
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    const std::array expected{
        core::FailureCode::ResourceExceeded,    core::FailureCode::FactorizationFailed,
        core::FailureCode::PrecisionUnresolved, core::FailureCode::Nonfinite,
        core::FailureCode::InvalidInput,        core::FailureCode::StaleGeneration};
    for (unsigned mutation = 1; mutation <= 6; ++mutation) {
        Rotation op;
        op.fault = mutation;
        op.fault_at = 3;
        core::ResourcePlan workspace(10000000), result(10000000);
        const auto failed = extract(op, fixed, wanted, limits(8, 3, 8), workspace, result);
        REQUIRE_FALSE(failed.has_value());
        CHECK(failed.error().cause.code == expected[mutation - 1]);
        CHECK(op.calls == 3);
        CHECK(failed.error().diagnostics.forward_applications == 3);
        if (mutation == 3) {
            CHECK(failed.error().cause.operation == core::Operation::CertifyError);
            REQUIRE(failed.error().cause.details.resources.has_value());
            CHECK(*failed.error().cause.details.resources ==
                  core::ResourceUsage{core::ResourceKind::Bytes, 900, 800});
        }
    }
    const Rotation limited;
    auto capacities = limits(8, 3, 8);
    capacities.actions.forward = 3;
    core::ResourcePlan workspace(10000000), result(10000000);
    const auto quota = extract(limited, fixed, wanted, capacities, workspace, result);
    REQUIRE_FALSE(quota.has_value());
    CHECK(quota.error().cause.code == core::FailureCode::ResourceExceeded);
    CHECK(limited.calls == 3);
    CHECK(quota.error().diagnostics.forward_applications == 3);
}

namespace {
struct RotationInverse final : eig::ActionProvider<Complex> {
    Complex sigma;
    bool left;
    mutable std::uint64_t calls = 0;
    bool fail = false;
    RotationInverse(Complex shift, bool adjoint) : sigma(shift), left(adjoint) {}
    [[nodiscard]] core::Index order() const noexcept override {
        return 8;
    }
    [[nodiscard]] eig::ActionStamp stamp() const noexcept override {
        return {context(), core::NumericalGeneration{left ? 61U : 60U}, sigma};
    }
    [[nodiscard]] numerics::SupportedPrecision
    supported_precision(eig::ActionKind) const noexcept override {
        return numerics::SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        return {};
    }
    [[nodiscard]] core::Result<void> apply(eig::ActionKind, eig::WordView<Complex> x,
                                           std::span<Complex> hi, std::span<Complex> lo,
                                           std::span<double> error) const override {
        ++calls;
        if (fail) {
            core::FailureDetails details;
            details.resources = core::ResourceUsage{core::ResourceKind::Bytes, 70, 60};
            return std::unexpected(core::Failure{core::FailureCode::LinearUnconverged,
                                                 core::Operation::SolveLinear, details});
        }
        const auto x0 = x.high()[0] + (x.low().empty() ? Complex{} : x.low()[0]);
        const auto x1 = x.high()[1] + (x.low().empty() ? Complex{} : x.low()[1]);
        const auto determinant = sigma * sigma + 25.0;
        const auto off = left ? -5.0 : 5.0;
        hi[0] = (-sigma * x0 + off * x1) / determinant;
        hi[1] = (-off * x0 - sigma * x1) / determinant;
        const auto bound = 128 * std::numeric_limits<double>::epsilon() * (std::abs(sigma) + 5) *
                               (std::abs(x0) + std::abs(x1)) / std::abs(determinant) +
                           64 * std::numeric_limits<double>::denorm_min();
        error[0] = std::nextafter(bound, std::numeric_limits<double>::infinity());
        error[1] = error[0];
        for (std::size_t i = 2; i < 8; ++i) {
            const auto value = x.high()[i] + (x.low().empty() ? Complex{} : x.low()[i]);
            hi[i] = value / (static_cast<double>(i) / 10 - sigma);
            error[i] =
                std::nextafter(64 * std::numeric_limits<double>::epsilon() * std::abs(hi[i]) +
                                   64 * std::numeric_limits<double>::denorm_min(),
                               std::numeric_limits<double>::infinity());
        }
        std::fill(lo.begin(), lo.end(), Complex{});
        return {};
    }
};
struct RotationInverseAccess final : eig::ArnoldiInverseAccess {
    RotationInverse right, left;
    mutable std::uint64_t certificates = 0;
    mutable std::optional<core::ErrorCertificate> first;
    bool replay = false, throw_qualification = false;
    explicit RotationInverseAccess(Complex sigma)
        : right(sigma, false), left(std::conj(sigma), true) {}
    [[nodiscard]] const eig::ActionProvider<Complex>&
    provider(eig::ArnoldiSide side) const override {
        return side == eig::ArnoldiSide::Right ? right : left;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    qualify(eig::ArnoldiSide side, eig::WordView<Complex> rhs, eig::WordView<Complex> solution,
            std::span<const double>, const eig::ActionStamp& stamp,
            eig::MetricRequirement requirement, core::CertificateContext expected) const override {
        ++certificates;
        if (throw_qualification) {
            throw std::runtime_error("known inverse qualifier exception");
        }
        if (replay && first) {
            return *first;
        }
        const auto sigma = side == eig::ArnoldiSide::Right ? right.sigma : left.sigma;
        if (stamp.shift != sigma) {
            return std::unexpected(
                core::Failure{core::FailureCode::StaleGeneration, core::Operation::SolveLinear});
        }
        double worst = 0;
        for (std::size_t i = 0; i < 8; ++i) {
            const auto current = solution.high()[i] + solution.low()[i];
            Complex action;
            if (i < 2) {
                const auto other = solution.high()[1 - i] + solution.low()[1 - i];
                const auto coefficient =
                    (i == 0 ? -5.0 : 5.0) * (side == eig::ArnoldiSide::Right ? 1.0 : -1.0);
                action = coefficient * other - sigma * current;
            } else {
                action = (static_cast<double>(i) / 10 - sigma) * current;
            }
            const auto upper = std::abs(rhs.high()[i] - action) +
                               128 * std::numeric_limits<double>::epsilon() *
                                   (std::abs(rhs.high()[i]) + std::abs(action) + std::abs(current));
            worst = std::max(worst, upper);
        }
        auto proof = core::ErrorCertificate::assess(expected, requirement.budget, worst);
        if (proof && !first) {
            first = *proof;
        }
        return proof;
    }
};
} // namespace
TEST_CASE("Arnoldi maps a complex shift to the real original pair and distinct inverse candidates",
          "[N4-07][EIG-09][EIG-11][negative][verification]") {
    const Rotation original;
    RotationInverseAccess inverse(Complex{1, 4});
    const std::array<double, 8> seed{1, 2, 3, 4, 5, 6, 7, 8}, scales{1, 2, .5, 4, 1, 2, .5, 4};
    const std::array<numerics::linear::UnitDimensions, 8> units{};
    const NoConsumer consumer;
    auto wanted = request<double>(seed, 1, 8, consumer);
    wanted.selection.emplace<eig::ArnoldiShiftInvert>(
        eig::ArnoldiShiftInvert{{1, 4},
                                inverse,
                                {1e-11, core::ObservationId{200}, core::ScopeId{201}},
                                core::NumericalGeneration{1000}});
    wanted.sides = eig::ArnoldiSides::Both;
    const eig::ArnoldiScales fixed{{scales, 2}, units, {}, context()};
    core::ResourcePlan workspace(10000000), result(10000000);
    auto extracted = extract(original, fixed, wanted, limits(8, 2, 8), workspace, result);
    INFO("failure code " << (extracted ? 255U
                                       : static_cast<unsigned>(extracted.error().cause.code)));
    REQUIRE(extracted.has_value());
    const auto values = extracted->values();
    REQUIRE(values.has_value());
    REQUIRE(checked(values).eigenvalues.values.size() == 2);
    CHECK(std::abs(checked(values).eigenvalues.values[0] - Complex{0, -5}) <= 1e-12);
    CHECK(std::abs(checked(values).eigenvalues.values[1] - Complex{0, 5}) <= 1e-12);
    CHECK(extracted->inverse_certificate_count() == inverse.right.calls + inverse.left.calls);
    for (std::size_t i = 0; i < extracted->inverse_certificate_count(); ++i) {
        const auto proof = extracted->inverse_certificate(static_cast<core::Index>(i));
        REQUIRE(proof.has_value());
        CHECK(checked(proof).context().numerical_generation == core::NumericalGeneration{1000 + i});
    }
    inverse.left.fail = true;
    core::ResourcePlan workspace2(10000000), result2(10000000);
    const auto right_before = inverse.right.calls, left_before = inverse.left.calls;
    const auto failed = extract(original, fixed, wanted, limits(8, 2, 8), workspace2, result2);
    REQUIRE_FALSE(failed.has_value());
    CHECK(failed.error().cause.code == core::FailureCode::LinearUnconverged);
    CHECK(failed.error().cause.operation == core::Operation::SolveLinear);
    CHECK(failed.error().diagnostics.inverse_applications ==
          (inverse.right.calls - right_before) + (inverse.left.calls - left_before));
    REQUIRE(failed.error().cause.details.resources.has_value());
    CHECK(*failed.error().cause.details.resources ==
          core::ResourceUsage{core::ResourceKind::Bytes, 70, 60});
    inverse.left.fail = false;
    inverse.replay = true;
    core::ResourcePlan workspace3(10000000), result3(10000000);
    const auto replayed = extract(original, fixed, wanted, limits(8, 2, 8), workspace3, result3);
    REQUIRE_FALSE(replayed.has_value());
    CHECK(replayed.error().cause.code == core::FailureCode::PrecisionUnresolved);
}

namespace {
struct ReentrantAction final : eig::ActionProvider<double> {
    Rotation base;
    std::function<void()> on_order;
    mutable bool visited = false;
    mutable std::optional<core::Failure> nested_failure;
    [[nodiscard]] core::Index order() const noexcept override {
        if (!visited) {
            visited = true;
            if (on_order) {
                on_order();
            }
        }
        return 8;
    }
    [[nodiscard]] eig::ActionStamp stamp() const noexcept override {
        return base.stamp();
    }
    [[nodiscard]] numerics::SupportedPrecision
    supported_precision(eig::ActionKind kind) const noexcept override {
        return base.supported_precision(kind);
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        if (nested_failure) {
            return std::unexpected(*nested_failure);
        }
        return {};
    }
    [[nodiscard]] core::Result<void> apply(eig::ActionKind kind, eig::WordView<double> x,
                                           std::span<double> hi, std::span<double> lo,
                                           std::span<double> error) const override {
        return base.apply(kind, x, hi, lo, error);
    }
};
} // namespace
TEST_CASE(
    "Arnoldi leases actual workspace before any provider metadata and reuses it independently",
    "[N4-07][EIG-12][negative][verification]") {
    const std::array<double, 8> seed{1, 2, 3, 4, 5, 6, 7, 8}, scales{1, 1, 1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 8> units{};
    const NoConsumer consumer;
    const auto wanted = request<double>(seed, 2, 8, consumer);
    const auto capacities = limits(8, 3, 8);
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    core::ResourcePlan scratch_plan(10000000), other_scratch(10000000), result(10000000),
        nested_result(10000000), other_result(10000000);
    auto workspace = eig::ArnoldiWorkspace::create(8, 3, scratch_plan);
    auto independent = eig::ArnoldiWorkspace::create(8, 3, other_scratch);
    REQUIRE(workspace.has_value());
    REQUIRE(independent.has_value());
    const auto prepared_payload = scratch_plan.allocated_bytes();
    ReentrantAction op;
    Rotation separate;
    std::optional<eig::ArnoldiSpectrum> independent_solution;
    bool nested_rejected = false, other_succeeded = false;
    op.on_order = [&] {
        auto nested =
            eig::arnoldi_spectrum(op, fixed, wanted, capacities, *workspace, nested_result);
        nested_rejected = !nested && nested.error().cause.code == core::FailureCode::InvalidInput;
        if (!nested) {
            op.nested_failure = nested.error().cause;
        }
        auto other =
            eig::arnoldi_spectrum(separate, fixed, wanted, capacities, *independent, other_result);
        other_succeeded = other.has_value();
        if (other) {
            independent_solution.emplace(std::move(*other));
        }
    };
    const auto failed = eig::arnoldi_spectrum(op, fixed, wanted, capacities, *workspace, result);
    REQUIRE_FALSE(failed.has_value());
    CHECK(failed.error().cause.code == core::FailureCode::InvalidInput);
    CHECK(nested_rejected);
    CHECK(other_succeeded);
    CHECK(op.base.calls == 0);
    CHECK(nested_result.allocated_bytes() == 0);
    REQUIRE(independent_solution.has_value());
    if (!independent_solution) {
        return;
    }
    const auto previous = independent_solution->values();
    REQUIRE(previous.has_value());
    const std::array<Complex, 2> before{checked(previous).eigenvalues.values[0],
                                        checked(previous).eigenvalues.values[1]};
    const Rotation clean;
    core::ResourcePlan result2(10000000);
    auto reused = eig::arnoldi_spectrum(clean, fixed, wanted, capacities, *workspace, result2);
    REQUIRE(reused.has_value());
    CHECK(scratch_plan.allocated_bytes() == prepared_payload);
    CHECK(checked(previous).eigenvalues.values[0] == before[0]);
    CHECK(checked(previous).eigenvalues.values[1] == before[1]);
    CHECK(checked(reused->values()).eigenvalues.values.data() !=
          checked(previous).eigenvalues.values.data());
    auto moved = std::move(*workspace);
    CHECK(workspace->order() == 0);
    CHECK(moved.order() == 8);
    core::ResourcePlan result3(10000000);
    const auto inactive =
        eig::arnoldi_spectrum(clean, fixed, wanted, capacities, *workspace, result3);
    REQUIRE_FALSE(inactive.has_value());
    CHECK(inactive.error().cause.code == core::FailureCode::InvalidInput);
    CHECK(result3.allocated_bytes() == 0);
}
TEST_CASE("Arnoldi validates dimensions precision scales generations and admission before "
          "numerical actions",
          "[N4-07][EIG-07][EIG-11][EIG-12][negative][verification]") {
    const Rotation op;
    const std::array<double, 8> seed{1, 2, 3, 4, 5, 6, 7, 8};
    std::array<double, 8> scales{1, 1, 1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 8> units{};
    const NoConsumer consumer;
    const auto wanted = request<double>(seed, 2, 8, consumer);
    core::ResourcePlan scratch_plan(10000000);
    auto workspace = eig::ArnoldiWorkspace::create(8, 3, scratch_plan);
    REQUIRE(workspace.has_value());
    const std::array expected{
        core::FailureCode::InvalidInput,        core::FailureCode::Nonfinite,
        core::FailureCode::InvalidInput,        core::FailureCode::StaleGeneration,
        core::FailureCode::Unsupported,         core::FailureCode::ResourceExceeded,
        core::FailureCode::PrecisionUnsupported};
    for (unsigned mutation = 0; mutation < expected.size(); ++mutation) {
        auto copy = wanted;
        auto capacities = limits(8, 3, 8);
        scales.fill(1);
        if (mutation == 0) {
            copy.ncv = copy.nev + 1;
        }
        if (mutation == 1) {
            scales[3] = std::numeric_limits<double>::infinity();
        }
        if (mutation == 2) {
            scales[3] = 0;
        }
        if (mutation == 3) {
            copy.returned = numerics::linear::Binary64Projection{core::NumericalGeneration{40},
                                                                 core::NumericalGeneration{40}};
        }
        if (mutation == 4) {
            capacities.backend_memory = numerics::linear::RequireProvenBound{100000000};
        }
        if (mutation == 5) {
            capacities.maximum_order = 7;
        }
        if (mutation == 6) {
            const std::array<double, 8> low{};
            const auto words = eig::WordView<double>::make(seed, low, core::Precision::TwoWord);
            REQUIRE(words.has_value());
            copy.initial = *words;
            const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
            core::ResourcePlan result(10000000);
            const auto failed =
                eig::arnoldi_spectrum(op, fixed, copy, capacities, *workspace, result);
            REQUIRE_FALSE(failed.has_value());
            CHECK(failed.error().cause.code == expected[mutation]);
            CHECK(op.calls == 0);
            continue;
        }
        const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
        core::ResourcePlan result(10000000);
        const auto failed = eig::arnoldi_spectrum(op, fixed, copy, capacities, *workspace, result);
        REQUIRE_FALSE(failed.has_value());
        CHECK(failed.error().cause.code == expected[mutation]);
        CHECK(op.calls == 0);
    }
}

TEST_CASE("Arnoldi exact coordinate conversion rejects subnormal division loss",
          "[N4-07][EIG-02][negative][verification]") {
    Complex hi{}, lo{};
    const auto rejected = numerics::eigen::detail::arnoldi_physical_input(
        Complex{std::numeric_limits<double>::denorm_min(), 0}, 1.5, eig::ArnoldiSide::Left, hi, lo);
    REQUIRE_FALSE(rejected.has_value());
    CHECK(rejected.error().code == core::FailureCode::PrecisionUnresolved);
    const auto exact = numerics::eigen::detail::arnoldi_physical_input(
        Complex{3, 6}, 1.5, eig::ArnoldiSide::Left, hi, lo);
    REQUIRE(exact.has_value());
    CHECK(hi + lo == Complex{2, 4});
    const auto product = numerics::eigen::detail::arnoldi_physical_input(
        Complex{.1, .2}, 1.5, eig::ArnoldiSide::Right, hi, lo);
    REQUIRE(product.has_value());
    using Rational = boost::multiprecision::cpp_rational;
    CHECK(Rational(hi.real()) + Rational(lo.real()) == Rational(.1) * Rational(1.5));
    CHECK(Rational(hi.imag()) + Rational(lo.imag()) == Rational(.2) * Rational(1.5));
}
TEST_CASE("Arnoldi rejects original residual budgets above the adopted ceiling",
          "[N4-07][EIG-03][negative][verification]") {
    const Rotation op;
    const std::array<double, 8> seed{1, 2, 3, 4, 5, 6, 7, 8}, scales{1, 1, 1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 8> units{};
    const NoConsumer consumer;
    auto wanted = request<double>(seed, 2, 8, consumer);
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    for (const bool subspace : {false, true}) {
        const eig::MetricRequirement loose{
            std::nextafter(1e-10, std::numeric_limits<double>::infinity()), core::ObservationId{31},
            core::ScopeId{32}};
        if (subspace) {
            wanted.qualification = eig::ArnoldiSubspaceRequest{
                {loose, {1e-12, core::ObservationId{33}, core::ScopeId{32}}, eig::NotRequested{}}};
        } else {
            wanted.qualification = eig::ArnoldiPairRequest{
                {loose, {1e-12, core::ObservationId{33}, core::ScopeId{32}}, eig::NotRequested{}}};
        }
        core::ResourcePlan scratch(10000000), result(10000000);
        const auto failed = extract(op, fixed, wanted, limits(8, 3, 8), scratch, result);
        REQUIRE_FALSE(failed.has_value());
        CHECK(failed.error().cause.code == core::FailureCode::InvalidInput);
        CHECK(op.calls == 0);
    }
}

TEST_CASE("Arnoldi normalizes extreme working seeds and timing leaves every value bit unchanged",
          "[N4-07][EIG-08][EIG-11][verification]") {
    const Rotation op;
    std::array<double, 8> seed{1, 2, 3, 4, 5, 6, 7, 8};
    const std::array<double, 8> scales{1, 1, 1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 8> units{};
    const NoConsumer consumer;
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    core::ResourcePlan scratch(10000000);
    auto workspace = eig::ArnoldiWorkspace::create(8, 3, scratch);
    REQUIRE(workspace.has_value());
    std::optional<eig::ArnoldiSpectrum> baseline;
    for (const auto exponent : {0, 900, -900}) {
        for (std::size_t i = 0; i < seed.size(); ++i) {
            seed[i] = std::scalbn(static_cast<double>(i + 1), exponent);
        }
        auto wanted = request<double>(seed, 2, 8, consumer);
        wanted.timing =
            exponent == 900 ? eig::ArnoldiTimingMode::Enabled : eig::ArnoldiTimingMode::Disabled;
        core::ResourcePlan result(10000000);
        auto extracted =
            eig::arnoldi_spectrum(op, fixed, wanted, limits(8, 3, 8), *workspace, result);
        REQUIRE(extracted.has_value());
        const auto metadata = extracted->metadata();
        REQUIRE(metadata.has_value());
        for (std::size_t i = 0; i < seed.size(); ++i) {
            CHECK(checked(metadata).initial_vector.values[i] == Complex{seed[i], 0});
        }
        const auto measured = extracted->measurements();
        REQUIRE(measured.has_value());
        CHECK(std::holds_alternative<eig::ArnoldiTimings>(*measured) == (exponent == 900));
        if (baseline) {
            const auto current = extracted->values(), previous = baseline->values();
            REQUIRE(current.has_value());
            REQUIRE(previous.has_value());
            for (std::size_t i = 0; i < 2; ++i) {
                CHECK(std::bit_cast<std::uint64_t>(checked(current).eigenvalues.values[i].real()) ==
                      std::bit_cast<std::uint64_t>(checked(previous).eigenvalues.values[i].real()));
                CHECK(std::bit_cast<std::uint64_t>(checked(current).eigenvalues.values[i].imag()) ==
                      std::bit_cast<std::uint64_t>(checked(previous).eigenvalues.values[i].imag()));
            }
        } else {
            baseline.emplace(std::move(*extracted));
        }
    }
}

TEST_CASE("Arnoldi sorts subspace roots without permuting a non diagonal Schur basis",
          "[N4-07][EIG-05][EIG-08][EIG-10][verification]") {
    Jordan op;
    op.second = 6;
    const std::array<double, 8> seed{1, 2, 3, 4, 5, 6, 7, 8}, scales{1, 1, 1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 8> units{};
    const NoConsumer consumer;
    auto wanted = request<double>(seed, 2, 8, consumer);
    wanted.qualification =
        eig::ArnoldiSubspaceRequest{{{1e-12, core::ObservationId{31}, core::ScopeId{32}},
                                     {1e-12, core::ObservationId{33}, core::ScopeId{32}},
                                     eig::NotRequested{}}};
    wanted.sides = eig::ArnoldiSides::Both;
    wanted.returned = numerics::linear::Binary64Projection{core::NumericalGeneration{40},
                                                           core::NumericalGeneration{41}};
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    core::ResourcePlan scratch(10000000), result(10000000);
    const auto extracted = extract(op, fixed, wanted, limits(8, 3, 8), scratch, result);
    INFO("failure code " << (extracted ? 255U
                                       : static_cast<unsigned>(extracted.error().cause.code)));
    REQUIRE(extracted.has_value());
    const auto values = extracted->values();
    REQUIRE(values.has_value());
    CHECK(std::abs(checked(values).eigenvalues.values[0] - Complex{5, 0}) <= 1e-12);
    CHECK(std::abs(checked(values).eigenvalues.values[1] - Complex{6, 0}) <= 1e-12);
    const auto* right = std::get_if<eig::ArnoldiSubspaceView>(&checked(values).right);
    const auto* left = std::get_if<eig::ArnoldiSubspaceView>(&checked(values).left);
    REQUIRE(right != nullptr);
    REQUIRE(left != nullptr);
    for (const auto* side : {right, left}) {
        REQUIRE(side->diagonal_columns.size() == 2);
        CHECK(side->diagonal_columns[0] != side->diagonal_columns[1]);
        CHECK(std::abs(side->triangular.values[2]) > .9);
        for (std::size_t j = 0; j < 2; ++j) {
            const auto column = static_cast<std::size_t>(side->diagonal_columns[j]);
            REQUIRE(column < 2);
            const auto actual = side == right
                                    ? side->triangular.values[column * 2 + column]
                                    : std::conj(side->triangular.values[column * 2 + column]);
            CHECK(std::abs(actual - checked(values).eigenvalues.values[j]) <= 1e-12);
        }
        for (std::size_t j = 0; j < 2; ++j) {
            for (std::size_t i = 0; i < 8; ++i) {
                const auto q = side->basis.values[j * 8 + i];
                Complex applied = (i == 0 ? 5.0 : (i == 1 ? 6.0 : static_cast<double>(i) / 10)) * q;
                if (side == right && i == 0) {
                    applied += side->basis.values[j * 8 + 1];
                }
                if (side == left && i == 1) {
                    applied += side->basis.values[j * 8];
                }
                Complex restricted{};
                for (std::size_t k = 0; k < 2; ++k) {
                    restricted +=
                        side->basis.values[k * 8 + i] * side->triangular.values[j * 2 + k];
                }
                CHECK(std::abs(applied - restricted) <= 1e-12);
            }
        }
    }
    const auto metadata = extracted->metadata();
    REQUIRE(metadata.has_value());
    CHECK(checked(metadata).backend.eigen_version == std::array<unsigned, 3>{5, 0, 1});
    CHECK(checked(metadata).backend.spectra_version == std::array<unsigned, 3>{1, 2, 0});
    CHECK(checked(metadata).original_stamp == op.stamp());
    CHECK(checked(metadata).settings.nev == 2);
    CHECK(checked(metadata).settings.ncv == 8);
}

TEST_CASE("Arnoldi validates reserved inverse generation ranges and qualifier failures",
          "[N4-07][EIG-09][EIG-12][negative][verification]") {
    const Rotation original;
    RotationInverseAccess inverse(Complex{1, 4});
    const std::array<double, 8> seed{1, 2, 3, 4, 5, 6, 7, 8}, scales{1, 1, 1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 8> units{};
    const NoConsumer consumer;
    auto wanted = request<double>(seed, 1, 8, consumer);
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    const std::array<std::uint64_t, 6> starts{5,  21, 40,
                                              60, 61, std::numeric_limits<std::uint64_t>::max()};
    for (std::size_t i = 0; i < starts.size(); ++i) {
        wanted.selection.emplace<eig::ArnoldiShiftInvert>(
            eig::ArnoldiShiftInvert{{1, 4},
                                    inverse,
                                    {1e-11, core::ObservationId{200}, core::ScopeId{201}},
                                    core::NumericalGeneration{starts[i]}});
        wanted.sides = eig::ArnoldiSides::Both;
        core::ResourcePlan scratch(10000000), result(10000000);
        const auto failed = extract(original, fixed, wanted, limits(8, 2, 8), scratch, result);
        REQUIRE_FALSE(failed.has_value());
        CHECK(failed.error().cause.code == (i == starts.size() - 1
                                                ? core::FailureCode::ResourceExceeded
                                                : core::FailureCode::StaleGeneration));
        CHECK(original.calls == 0);
        CHECK(inverse.right.calls == 0);
        CHECK(inverse.left.calls == 0);
    }
    wanted.selection.emplace<eig::ArnoldiShiftInvert>(
        eig::ArnoldiShiftInvert{{1, 4},
                                inverse,
                                {1e-11, core::ObservationId{200}, core::ScopeId{201}},
                                core::NumericalGeneration{1000}});
    inverse.throw_qualification = true;
    core::ResourcePlan scratch(10000000), result(10000000);
    const auto failure = extract(original, fixed, wanted, limits(8, 2, 8), scratch, result);
    REQUIRE_FALSE(failure.has_value());
    CHECK(failure.error().cause.code == core::FailureCode::PrecisionUnresolved);
    CHECK(failure.error().diagnostics.inverse_applications == 1);
    CHECK(inverse.right.calls == 1);
    inverse.throw_qualification = false;
    wanted.selection.emplace<eig::ArnoldiShiftInvert>(eig::ArnoldiShiftInvert{
        {1, 4},
        inverse,
        {1e-11, core::ObservationId{200}, core::ScopeId{201}},
        core::NumericalGeneration{std::numeric_limits<std::uint64_t>::max()}});
    auto zero_limits = limits(8, 2, 8);
    zero_limits.actions.inverse = 0;
    core::ResourcePlan scratch2(10000000), result2(10000000);
    const auto no_capacity = extract(original, fixed, wanted, zero_limits, scratch2, result2);
    REQUIRE_FALSE(no_capacity.has_value());
    CHECK(no_capacity.error().cause.code == core::FailureCode::ResourceExceeded);
    CHECK(no_capacity.error().diagnostics.inverse_applications == 0);
    CHECK(inverse.right.calls == 1);
}

namespace {
struct DiagonalInverse final : eig::ActionProvider<Complex> {
    std::array<Complex, 6> diagonal;
    Complex sigma;
    bool left;
    double multiplier = 1;
    DiagonalInverse(std::array<Complex, 6> values, Complex shift, bool adjoint)
        : diagonal(values), sigma(shift), left(adjoint) {}
    [[nodiscard]] core::Index order() const noexcept override {
        return 6;
    }
    [[nodiscard]] eig::ActionStamp stamp() const noexcept override {
        return {context(), core::NumericalGeneration{left ? 61U : 60U}, sigma};
    }
    [[nodiscard]] numerics::SupportedPrecision
    supported_precision(eig::ActionKind) const noexcept override {
        return numerics::SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        return {};
    }
    [[nodiscard]] core::Result<void> apply(eig::ActionKind, eig::WordView<Complex> x,
                                           std::span<Complex> hi, std::span<Complex> lo,
                                           std::span<double> error) const override {
        for (std::size_t i = 0; i < 6; ++i) {
            const auto coefficient = (left ? std::conj(diagonal[i]) : diagonal[i]) - sigma;
            const auto input = x.high()[i] + (x.low().empty() ? Complex{} : x.low()[i]);
            const auto inverse = input / coefficient;
            hi[i] = multiplier * inverse;
            lo[i] = {};
            // This deliberately includes the full bias in the negative control,
            // rather than pretending that the biased callback is the exact inverse.
            error[i] = std::nextafter(
                (std::abs(multiplier - 1) + 128 * std::numeric_limits<double>::epsilon()) *
                        std::abs(inverse) +
                    128 * std::numeric_limits<double>::denorm_min(),
                std::numeric_limits<double>::infinity());
        }
        return {};
    }
};
struct DiagonalInverseAccess final : eig::ArnoldiInverseAccess {
    DiagonalInverse right, left;
    DiagonalInverseAccess(std::array<Complex, 6> values, Complex shift)
        : right(values, shift, false), left(values, std::conj(shift), true) {}
    [[nodiscard]] const eig::ActionProvider<Complex>&
    provider(eig::ArnoldiSide side) const override {
        return side == eig::ArnoldiSide::Right ? right : left;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    qualify(eig::ArnoldiSide side, eig::WordView<Complex> rhs, eig::WordView<Complex> solution,
            std::span<const double>, const eig::ActionStamp&, eig::MetricRequirement requirement,
            core::CertificateContext expected) const override {
        const auto& provider = side == eig::ArnoldiSide::Right ? right : left;
        double worst = 0;
        for (std::size_t i = 0; i < 6; ++i) {
            const auto coefficient =
                (provider.left ? std::conj(provider.diagonal[i]) : provider.diagonal[i]) -
                provider.sigma;
            const auto value = solution.high()[i] + solution.low()[i];
            const auto applied = coefficient * value;
            const auto residual = std::abs(rhs.high()[i] - applied) +
                                  128 * std::numeric_limits<double>::epsilon() *
                                      (std::abs(rhs.high()[i]) + std::abs(applied));
            worst = std::max(worst, residual);
        }
        return core::ErrorCertificate::assess(expected, requirement.budget, worst);
    }
};
} // namespace
TEST_CASE("Arnoldi uses original complex eigenvalues after shifted extraction",
          "[N4-07][EIG-09][verification]") {
    const std::array<Complex, 6> diagonal{{{8, 2}, {6, -1}, {0, 4}, {2, 0}, {1, 0}, {.5, 0}}};
    std::vector<Complex> matrix(36);
    for (std::size_t i = 0; i < 6; ++i) {
        matrix[6 * i + i] = diagonal[i];
    }
    const LiteralAction<Complex> op(6, std::move(matrix));
    const DiagonalInverseAccess inverse(diagonal, Complex{5, -1});
    const std::array<Complex, 6> seed{{{1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1}, {6, 1}}};
    const std::array<double, 6> scales{1, 2, .5, 4, 1, 2};
    const std::array<numerics::linear::UnitDimensions, 6> units{};
    const NoConsumer consumer;
    auto wanted = request<Complex>(seed, 1, 6, consumer);
    wanted.selection.emplace<eig::ArnoldiShiftInvert>(
        eig::ArnoldiShiftInvert{{5, -1},
                                inverse,
                                {1e-11, core::ObservationId{200}, core::ScopeId{201}},
                                core::NumericalGeneration{1000}});
    wanted.sides = eig::ArnoldiSides::Both;
    const eig::ArnoldiScales fixed{{scales, 2}, units, {}, context()};
    core::ResourcePlan scratch(10000000), result(10000000);
    const auto extracted = extract(op, fixed, wanted, limits(6, 1, 6), scratch, result);
    INFO("failure code " << (extracted ? 255U
                                       : static_cast<unsigned>(extracted.error().cause.code)));
    REQUIRE(extracted.has_value());
    const auto values = extracted->values();
    REQUIRE(values.has_value());
    CHECK(std::abs(checked(values).eigenvalues.values[0] - Complex{6, -1}) <= 1e-12);
    const auto bound = extracted->inverse_assessment(0);
    REQUIRE(bound.has_value());
    CHECK(checked(bound).native_output_error_upper > 0);
    CHECK(std::isfinite(checked(bound).native_output_error_upper));
    CHECK(checked(bound).candidate == core::NumericalGeneration{1000});
}
TEST_CASE("Arnoldi never replaces selected shifted roots with a different original restriction",
          "[N4-07][EIG-08][EIG-09][negative][verification]") {
    const std::array<Complex, 6> diagonal{{{8, 0}, {4, 0}, {2, 0}, {1, 0}, {.5, 0}, {.25, 0}}};
    std::vector<Complex> matrix(36);
    for (std::size_t i = 0; i < 6; ++i) {
        matrix[6 * i + i] = diagonal[i];
    }
    const LiteralAction<Complex> op(6, std::move(matrix));
    DiagonalInverseAccess inverse(diagonal, Complex{});
    inverse.right.multiplier = 2;
    const std::array<Complex, 6> seed{{{1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}}};
    const std::array<double, 6> scales{1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 6> units{};
    const NoConsumer consumer;
    auto wanted = request<Complex>(seed, 1, 6, consumer);
    wanted.selection.emplace<eig::ArnoldiShiftInvert>(
        eig::ArnoldiShiftInvert{{},
                                inverse,
                                {10, core::ObservationId{200}, core::ScopeId{201}},
                                core::NumericalGeneration{1000}});
    wanted.qualification =
        eig::ArnoldiSubspaceRequest{{{1e-12, core::ObservationId{31}, core::ScopeId{32}},
                                     {1e-12, core::ObservationId{33}, core::ScopeId{32}},
                                     eig::NotRequested{}}};
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    core::ResourcePlan scratch(10000000), result(10000000);
    const auto failed = extract(op, fixed, wanted, limits(6, 1, 6), scratch, result);
    REQUIRE_FALSE(failed.has_value());
    CHECK(failed.error().cause.code == core::FailureCode::PrecisionUnresolved);
    // Native theta=8 would map to .125; the original restriction is .25.
    // Its invariant-space residual alone is tiny, so the set-matching gate is essential.
}

TEST_CASE("Arnoldi cannot double count a real mode through an imaginary rounding tail",
          "[N4-07][EIG-08][negative][verification]") {
    const double epsilon = 1e-15;
    const std::array<Complex, 2> nearly_real{{{1, 0}, {0, epsilon}}},
        nearly_conjugate{{{1, 0}, {0, -epsilon}}};
    const std::array<double, 2> scales{1, 1};
    const auto rejected = eig::detail::arnoldi_conjugate_independence(
        nearly_real, nearly_conjugate, scales, eig::ArnoldiSide::Right);
    REQUIRE_FALSE(rejected.has_value());
    CHECK(rejected.error().code == core::FailureCode::PrecisionUnresolved);
    const auto normalization = std::sqrt(.5);
    const std::array<Complex, 2> pair{{{normalization, 0}, {0, normalization}}},
        conjugate{{{normalization, 0}, {0, -normalization}}};
    CHECK(eig::detail::arnoldi_conjugate_independence(pair, conjugate, scales,
                                                      eig::ArnoldiSide::Right)
              .has_value());
    const std::array<Complex, 6> diagonal{{{8, 0}, {4, 0}, {2, 0}, {1, 0}, {.5, 0}, {.25, 0}}};
    std::vector<double> matrix(36);
    for (std::size_t i = 0; i < 6; ++i) {
        matrix[6 * i + i] = diagonal[i].real();
    }
    const LiteralAction<double> op(6, std::move(matrix));
    const DiagonalInverseAccess inverse(diagonal, Complex{.125, .125});
    const std::array<double, 6> seed{1, 2, 3, 4, 5, 6}, components{1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 6> units{};
    const NoConsumer consumer;
    auto wanted = request<double>(seed, 1, 6, consumer);
    wanted.selection.emplace<eig::ArnoldiShiftInvert>(
        eig::ArnoldiShiftInvert{{.125, .125},
                                inverse,
                                {1e-11, core::ObservationId{200}, core::ScopeId{201}},
                                core::NumericalGeneration{1000}});
    const eig::ArnoldiScales fixed{{components, 1}, units, {}, context()};
    core::ResourcePlan scratch(10000000), result(10000000);
    const auto extracted = extract(op, fixed, wanted, limits(6, 2, 6), scratch, result);
    if (extracted) {
        const auto values = extracted->values();
        REQUIRE(values.has_value());
        REQUIRE(checked(values).eigenvalues.values.size() == 1);
        CHECK(checked(values).eigenvalues.values[0].imag() == 0);
        CHECK(std::abs(checked(values).eigenvalues.values[0].real() - .25) <= 1e-12);
    } else {
        CHECK(extracted.error().cause.code == core::FailureCode::PrecisionUnresolved);
    }
}

TEST_CASE("Arnoldi uses one owning snapshot for values scales units and required budgets",
          "[N4-07][EIG-10][EIG-12][negative][verification]") {
    Rotation provider;
    std::array<double, 8> seed{1, 2, 3, 4, 5, 6, 7, 8}, components{1, 1, 1, 1, 1, 1, 1, 1};
    std::array<numerics::linear::UnitDimensions, 8> units{};
    std::array<numerics::linear::ConsumerCheck, 2> checks{
        {{numerics::linear::CheckKind::Observation, core::ObservationId{100}, core::ScopeId{101},
          1e-12},
         {numerics::linear::CheckKind::Reaction, core::ObservationId{102}, core::ScopeId{103},
          1e-12}}};
    Consumer evaluator;
    const auto words = eig::WordView<double>::make(seed, {}, core::Precision::Binary64);
    REQUIRE(words.has_value());
    eig::ArnoldiRequest<double> wanted{
        2,
        8,
        *words,
        1e-14,
        100,
        eig::ArnoldiLargestMagnitude{},
        accuracy(),
        eig::ArnoldiSides::RightOnly,
        {1e-12, 1e-12},
        numerics::linear::RetainedPrecision{core::NumericalGeneration{40}},
        {},
        checks,
        evaluator,
        eig::ArnoldiTimingMode::Disabled};
    eig::ArnoldiScales fixed{{components, 1}, units, {}, context()};
    auto capacities = limits(8, 3, 8);
    provider.on_first_action = [&] {
        components.fill(std::numeric_limits<double>::min());
        seed.fill(0);
        units[0][0] = 1;
        wanted.returned = numerics::linear::RetainedPrecision{core::NumericalGeneration{777}};
        wanted.backend_iterations = 0;
        checks[1].budget = 1;
        checks[1].observation = core::ObservationId{999};
        capacities.actions.forward = 0;
        fixed.context.trial = core::TrialId{888};
    };
    core::ResourcePlan scratch(10000000), result(10000000);
    const auto extracted = extract(provider, fixed, wanted, capacities, scratch, result);
    REQUIRE(extracted.has_value());
    const auto metadata = extracted->metadata();
    CHECK(checked(metadata).components[0] == 1);
    CHECK(checked(metadata).component_units[0][0] == 0);
    CHECK(checked(metadata).initial_vector.values[0] == Complex{1, 0});
    CHECK(checked(metadata).context.trial == core::TrialId{11});
    CHECK(checked(metadata).settings.backend_iteration_limit == 100);
    CHECK(checked(metadata).limits.actions.forward == 10000);
    CHECK(checked(extracted->values()).generation == core::NumericalGeneration{40});
    const auto second = extracted->required_check(1);
    REQUIRE(second.has_value());
    CHECK(checked(second).budget == 1e-12);
    CHECK(checked(second).observation == core::ObservationId{102});
    provider.calls = 0;
    components.fill(1);
    units[0][0] = 0;
    fixed.context = context();
    capacities = limits(8, 3, 8);
    for (std::size_t i = 0; i < seed.size(); ++i) {
        seed[i] = static_cast<double>(i + 1);
    }
    wanted.returned = numerics::linear::RetainedPrecision{core::NumericalGeneration{40}};
    wanted.backend_iterations = 100;
    checks[1].budget = 1e-12;
    checks[1].observation = core::ObservationId{102};
    evaluator.calls = 0;
    evaluator.fixed_second_upper = 2e-12;
    core::ResourcePlan scratch2(10000000), result2(10000000);
    const auto refused = extract(provider, fixed, wanted, capacities, scratch2, result2);
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error().cause.code == core::FailureCode::PrecisionUnresolved);
    REQUIRE(refused.error().cause.details.problem.has_value());
    CHECK(*refused.error().cause.details.problem == core::ProblemId{core::ObservationId{102}});
    CHECK(refused.error().cause.details.numerical_generation == core::NumericalGeneration{40});
    CHECK(refused.error().cause.details.trial == core::TrialId{11});
    CHECK(checked(extracted->values()).generation == core::NumericalGeneration{40});
}
TEST_CASE("Arnoldi keeps producer uncertainty and the final frozen stamp gate",
          "[N4-07][EIG-03][EIG-12][negative][verification]") {
    const std::array<double, 8> seed{1, 2, 3, 4, 5, 6, 7, 8}, components{1, 1, 1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 8> units{};
    const NoConsumer unused;
    const auto wanted = request<double>(seed, 2, 8, unused);
    const eig::ArnoldiScales fixed{{components, 1}, units, {}, context()};
    Rotation uncertain;
    uncertain.uncertainty = 1e-3;
    core::ResourcePlan scratch(10000000), result(10000000);
    const auto refused = extract(uncertain, fixed, wanted, limits(8, 3, 8), scratch, result);
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error().cause.code == core::FailureCode::PrecisionUnresolved);
    const Rotation mutating;
    Consumer evaluator;
    evaluator.mutate_stamp = &mutating;
    const std::array<numerics::linear::ConsumerCheck, 1> checks{
        {{numerics::linear::CheckKind::Observation, core::ObservationId{100}, core::ScopeId{101},
          1e-12}}};
    const eig::ArnoldiRequest<double> checked_request{wanted.nev,
                                                      wanted.ncv,
                                                      wanted.initial,
                                                      wanted.backend_tolerance,
                                                      wanted.backend_iterations,
                                                      wanted.selection,
                                                      wanted.qualification,
                                                      wanted.sides,
                                                      wanted.association,
                                                      wanted.returned,
                                                      wanted.state,
                                                      checks,
                                                      evaluator,
                                                      wanted.timing};
    core::ResourcePlan scratch2(10000000), result2(10000000);
    const auto stale =
        extract(mutating, fixed, checked_request, limits(8, 3, 8), scratch2, result2);
    REQUIRE_FALSE(stale.has_value());
    CHECK(stale.error().cause.code == core::FailureCode::StaleGeneration);
}

namespace {
[[nodiscard]] std::size_t exact_phase_maximum(std::span<const Complex> values,
                                              std::span<const double> scales,
                                              eig::ArnoldiSide side) {
    using R = boost::multiprecision::cpp_rational;
    const auto energy = [&](std::size_t i) {
        const R real{values[i].real()}, imaginary{values[i].imag()}, d{scales[i]};
        const R square = real * real + imaginary * imaginary;
        if (side == eig::ArnoldiSide::Right) {
            return R{square / (d * d)};
        }
        return R{square * d * d};
    };
    std::size_t maximum = 0;
    for (std::size_t i = 1; i < values.size(); ++i) {
        if (energy(i) > energy(maximum)) {
            maximum = i;
        }
    }
    return maximum;
}
} // namespace
TEST_CASE("Arnoldi canonical phase inspects the actual returned near tie",
          "[N4-07][EIG-04][negative][verification]") {
    const std::array<std::array<double, 2>, 2> policies{{{1, 1}, {3, .2}}};
    for (const auto& scales : policies) {
        for (const auto side : {eig::ArnoldiSide::Right, eig::ArnoldiSide::Left}) {
            for (const auto tail : {0.0, std::ldexp(1.0, -30)}) {
                std::array<Complex, 2> values{{{-1, 0}, {1, tail}}}, backup{}, candidate{};
                const auto normalized =
                    eig::detail::arnoldi_normalize(values, scales, side, backup, candidate);
                const bool rounding_cycle_boundary = scales[0] == 1 && tail != 0;
                if (rounding_cycle_boundary && !normalized) {
                    // For a=round(1/sqrt(2)), either anchor transfers the tiny
                    // imaginary component to its rival; no actual maximum is its
                    // own anchor. Another libm may return a valid canonical value.
                    CHECK(normalized.error().code == core::FailureCode::PrecisionUnresolved);
                    continue;
                }
                REQUIRE(normalized.has_value());
                const auto maximum = exact_phase_maximum(values, scales, side);
                CHECK(values[maximum].imag() == 0);
                CHECK(values[maximum].real() >= 0);
            }
        }
    }
    const double tail = std::ldexp(1.0, -30);
    std::vector<Complex> matrix(36);
    const std::array<double, 6> diagonal{2, 1, .8, .6, .4, .2};
    for (std::size_t i = 0; i < 6; ++i) {
        matrix[6 * i + i] = diagonal[i];
    }
    matrix[6] = {-1, -tail};
    const LiteralAction<Complex> provider(6, std::move(matrix));
    const std::array<Complex, 6> seed{{{1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1}, {6, 1}}};
    const std::array<double, 6> scales{1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 6> units{};
    const NoConsumer consumer;
    auto wanted = request<Complex>(seed, 1, 6, consumer);
    wanted.sides = eig::ArnoldiSides::Both;
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    core::ResourcePlan scratch(10000000), result(10000000);
    const auto extracted = extract(provider, fixed, wanted, limits(6, 1, 6), scratch, result);
    REQUIRE(extracted.has_value());
    const auto view = extracted->values();
    const auto* right = std::get_if<eig::ArnoldiPairView>(&checked(view).right);
    const auto* left = std::get_if<eig::ArnoldiPairView>(&checked(view).left);
    REQUIRE(right != nullptr);
    REQUIRE(left != nullptr);
    for (const auto* side : {right, left}) {
        const auto maximum =
            exact_phase_maximum(side->vectors.values, scales,
                                side == right ? eig::ArnoldiSide::Right : eig::ArnoldiSide::Left);
        CHECK(side->vectors.values[maximum].imag() == 0);
        CHECK(side->vectors.values[maximum].real() >= 0);
    }
}
