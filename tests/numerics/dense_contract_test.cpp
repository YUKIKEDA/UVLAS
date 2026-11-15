#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>
#include <uvlas/numerics/dense_solve.hpp>

namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;
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
class Consumer final : public ConsumerEvaluator<double> {
public:
    std::function<void()> hook;
    bool throwing = false;
    [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
        return SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<double>&, const PrecisionBuffer<double>&, const Context&,
           const ConsumerCheck& request, core::CertificateContext expected) const override {
        if (throwing) {
            throw std::runtime_error("consumer failure");
        }
        if (hook) {
            hook();
        }
        return core::ErrorCertificate::assess(expected, request.budget, 0);
    }
};
struct Inputs {
    core::ResourcePlan source_plan{16384};
    core::ResourcePlan workspace_plan{65536};
    core::ResourcePlan result_plan{65536};
    std::array<double, 2> scale_values{1, 1};
    std::array<UnitDimensions, 2> units{UnitDimensions{}, UnitDimensions{}};
    Consumer consumer;
    DenseReal matrix() {
        // Column-major A=[[2,1],[0,3]], x=[1,2], b=[4,6].
        const std::array coefficients{2.0, 0.0, 1.0, 3.0};
        auto value =
            DenseReal::copy({2, 2, 2}, coefficients, source_plan, core::Operation::SolveLinear);
        REQUIRE(value);
        return std::move(*value);
    }
    PrecisionBuffer<double> rhs() {
        const std::array values{4.0, 6.0};
        auto matrix = DenseReal::copy({2, 1, 2}, values, source_plan, core::Operation::SolveLinear);
        REQUIRE(matrix);
        auto value =
            PrecisionBuffer<double>::binary64(std::move(*matrix), core::NumericalGeneration{30});
        REQUIRE(value);
        return std::move(*value);
    }
    [[nodiscard]] Scales scales(Context identity = context()) const {
        return {scale_values, scale_values, units, units, identity};
    }
    [[nodiscard]] Limits limits() const {
        return {2, 1, 1, 5, EstimateAdmission{1024ULL * 1024ULL}};
    }
    [[nodiscard]] Accuracy accuracy() const {
        return {Norm::Infinity,         1.e-12, 0, 1.e-12, core::ObservationId{20},
                core::ObservationId{21}};
    }
    DensePrepared prepare(const DenseReal& matrix) {
        auto value = DensePrepared::prepare(DenseOperator{matrix, context()}, scales(), limits(),
                                            workspace_plan);
        REQUIRE(value);
        return std::move(*value);
    }
    auto solve(const DenseReal& matrix, const PrecisionBuffer<double>& rhs, DensePrepared& prepared,
               std::span<const ConsumerCheck> requirements = {}) {
        return solve_linear(DenseOperator{matrix, context()},
                            Rhs<double>{rhs, context(), core::RhsGeneration{1}}, scales(),
                            Precision<double>{RetainedPrecision{core::NumericalGeneration{40}},
                                              requirements, consumer},
                            accuracy(), limits(), prepared, result_plan);
    }
};
std::span<const double> values(const Solution<double>& solution) {
    REQUIRE(solution.values().precision() == core::Precision::Binary64);
    return std::get<Binary64<double>>(solution.values().storage()).values.values();
}
constexpr std::array requirement{
    ConsumerCheck{CheckKind::Reaction, core::ObservationId{50}, core::ScopeId{51}, 1.e-6}};
} // namespace

TEST_CASE("Dense result survives source factor and plan destruction and later solve reuse",
          "[integration][N4-05][N4-06]") {
    auto surviving = [] {
        Inputs input;
        const auto matrix = input.matrix();
        const auto rhs = input.rhs();
        auto prepared = input.prepare(matrix);
        return input.solve(matrix, rhs, prepared);
    }();
    REQUIRE(surviving);
    CHECK(std::abs(values(*surviving)[0] - 1) <= 1.e-12);
    CHECK(std::abs(values(*surviving)[1] - 2) <= 1.e-12);
    CHECK(surviving->certificate().unknown_scales().size() == 2);

    Inputs input;
    const auto matrix = input.matrix();
    const auto rhs = input.rhs();
    auto prepared = input.prepare(matrix);
    auto first = input.solve(matrix, rhs, prepared);
    REQUIRE(first);
    const auto first_values = values(*first);
    const std::array snapshot{first_values[0], first_values[1]};
    auto second = input.solve(matrix, rhs, prepared);
    REQUIRE(second);
    CHECK(values(*second).data() != first_values.data());
    CHECK(std::ranges::equal(first_values, snapshot));
}

TEST_CASE(
    "Dense cache rejects changed actual coefficients even when context and address are reused",
    "[negative][N4-03][N4-05]") {
    Inputs input;
    auto matrix = input.matrix();
    const auto rhs = input.rhs();
    auto prepared = input.prepare(matrix);
    auto previous = input.solve(matrix, rhs, prepared);
    REQUIRE(previous);
    matrix.values()[0] = 4;
    const auto rejected = input.solve(matrix, rhs, prepared);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::StaleGeneration);
    CHECK(matrix.values()[0] == 4);
    CHECK(std::abs(values(*previous)[0] - 1) <= 1.e-12);
    matrix.values()[0] = 2;
    CHECK(input.solve(matrix, rhs, prepared));
}

TEST_CASE("Dense cached identity cannot be relabelled coherently across all fresh inputs",
          "[negative][N4-03]") {
    Inputs input;
    const auto matrix = input.matrix();
    const auto rhs = input.rhs();
    auto prepared = input.prepare(matrix);
    auto changed = context();
    changed.cache.fingerprint.bytes[0] = std::byte{1};
    const auto rejected = solve_linear(
        DenseOperator{matrix, changed}, Rhs<double>{rhs, changed, core::RhsGeneration{99}},
        input.scales(changed),
        Precision<double>{RetainedPrecision{core::NumericalGeneration{40}}, {}, input.consumer},
        input.accuracy(), input.limits(), prepared, input.result_plan);
    REQUIRE_FALSE(rejected);
    CHECK(rejected.error().code == core::FailureCode::StaleGeneration);
}

TEST_CASE("Dense preparation and each reuse enforce explicit backend admission and output budget",
          "[negative][N4-09]") {
    Inputs input;
    const auto matrix = input.matrix();
    const auto rhs = input.rhs();
    auto prepared = input.prepare(matrix);
    auto strict = input.limits();
    strict.backend_memory = RequireProvenBound{1024ULL * 1024ULL};
    const auto unsupported = solve_linear(
        DenseOperator{matrix, context()}, Rhs<double>{rhs, context(), core::RhsGeneration{1}},
        input.scales(),
        Precision<double>{RetainedPrecision{core::NumericalGeneration{40}}, {}, input.consumer},
        input.accuracy(), strict, prepared, input.result_plan);
    REQUIRE_FALSE(unsupported);
    CHECK(unsupported.error().code == core::FailureCode::Unsupported);
    auto tiny = input.limits();
    tiny.backend_memory = EstimateAdmission{0};
    const auto exceeded = DensePrepared::prepare(DenseOperator{matrix, context()}, input.scales(),
                                                 tiny, input.workspace_plan);
    REQUIRE_FALSE(exceeded);
    CHECK(exceeded.error().code == core::FailureCode::ResourceExceeded);
    CHECK_FALSE(exceeded.error().details.resources);
    core::ResourcePlan none{0};
    const auto no_output = solve_linear(
        DenseOperator{matrix, context()}, Rhs<double>{rhs, context(), core::RhsGeneration{1}},
        input.scales(),
        Precision<double>{RetainedPrecision{core::NumericalGeneration{40}}, {}, input.consumer},
        input.accuracy(), input.limits(), prepared, none);
    REQUIRE_FALSE(no_output);
    CHECK(no_output.error().code == core::FailureCode::ResourceExceeded);
    CHECK(input.solve(matrix, rhs, prepared));
}

TEST_CASE("Moved preparation rejects use and the new owner remains usable", "[negative][N4-05]") {
    Inputs input;
    const auto matrix = input.matrix();
    const auto rhs = input.rhs();
    auto prepared = input.prepare(matrix);
    auto moved = std::move(prepared);
    // NOLINTNEXTLINE(bugprone-use-after-move): inactive queries are the documented move contract.
    CHECK_FALSE(prepared.active());
    CHECK_FALSE(prepared.memory_plan());
    CHECK(moved.active());
    CHECK_FALSE(input.solve(matrix, rhs, prepared));
    CHECK(input.solve(matrix, rhs, moved));
}

TEST_CASE("Dense callback failure and same-workspace reentry preserve reusable preparation",
          "[negative][N4-03][N4-05][N4-06]") {
    Inputs input;
    const auto matrix = input.matrix();
    const auto rhs = input.rhs();
    auto prepared = input.prepare(matrix);
    input.consumer.throwing = true;
    CHECK_FALSE(input.solve(matrix, rhs, prepared, requirement));
    input.consumer.throwing = false;
    bool rejected_inner = false;
    input.consumer.hook = [&] {
        const Consumer quiet;
        core::ResourcePlan output{65536};
        const auto inner = solve_linear(
            DenseOperator{matrix, context()}, Rhs<double>{rhs, context(), core::RhsGeneration{1}},
            input.scales(),
            Precision<double>{RetainedPrecision{core::NumericalGeneration{41}}, {}, quiet},
            input.accuracy(), input.limits(), prepared, output);
        rejected_inner = !inner && inner.error().code == core::FailureCode::InvalidInput;
    };
    auto outer = input.solve(matrix, rhs, prepared, requirement);
    REQUIRE(outer);
    CHECK(rejected_inner);
    input.consumer.hook = {};
    CHECK(input.solve(matrix, rhs, prepared));
    const std::array original_matrix{2.0, 0.0, 1.0, 3.0};
    const std::array original_rhs{4.0, 6.0};
    CHECK(std::ranges::equal(matrix.values(), original_matrix));
    CHECK(std::ranges::equal(std::get<Binary64<double>>(rhs.storage()).values.values(),
                             original_rhs));
}

TEST_CASE("A qualified total does not qualify a distinct local load observation",
          "[verification][negative][N4-05][N4-06]") {
    class LinearOutputs final : public ConsumerEvaluator<double> {
    public:
        mutable unsigned calls = 0;
        [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
            return SupportedPrecision::Binary64Only;
        }
        [[nodiscard]] core::Result<core::ErrorCertificate>
        assess(const PrecisionBuffer<double>&, const PrecisionBuffer<double>& candidate,
               const Context&, const ConsumerCheck& check,
               core::CertificateContext expected) const override {
            ++calls;
            const auto& x = std::get<Binary64<double>>(candidate.storage()).values;
            // This dyadic fixture makes both the sum and subtraction exact in binary64.
            const bool total = check.observation == core::ObservationId{80};
            const double quantity = total ? x.values()[0] + x.values()[1] : x.values()[0];
            const double reference = total ? 2.0 : 1.0;
            return core::ErrorCertificate::assess(expected, check.budget,
                                                  std::abs(quantity - reference));
        }
    };
    Inputs input;
    auto matrix = input.matrix();
    std::ranges::copy(std::array{1.0, 0.0, 0.0, 1.0}, matrix.values().begin());
    const double perturbation = std::ldexp(1.0, -20);
    const std::array loads{1 + perturbation, 1 - perturbation};
    auto b = DenseReal::copy({2, 1, 2}, loads, input.source_plan, core::Operation::SolveLinear);
    REQUIRE(b);
    auto rhs = PrecisionBuffer<double>::binary64(std::move(*b), core::NumericalGeneration{30});
    REQUIRE(rhs);
    auto prepared = input.prepare(matrix);
    const LinearOutputs evaluator;
    const std::array observations{
        ConsumerCheck{CheckKind::Observation, core::ObservationId{80}, core::ScopeId{81}, 1.e-8},
        ConsumerCheck{CheckKind::Observation, core::ObservationId{82}, core::ScopeId{83}, 1.e-8}};
    const auto solve_with = [&](std::span<const ConsumerCheck> required) {
        return solve_linear(DenseOperator{matrix, context()},
                            Rhs<double>{*rhs, context(), core::RhsGeneration{1}}, input.scales(),
                            Precision<double>{RetainedPrecision{core::NumericalGeneration{40}},
                                              required, evaluator},
                            input.accuracy(), input.limits(), prepared, input.result_plan);
    };
    const auto total_only = solve_with(std::span{observations}.first(1));
    REQUIRE(total_only);
    CHECK(evaluator.calls == 1);
    CHECK(std::abs(values(*total_only)[0] + values(*total_only)[1] - 2) <= 1.e-12);
    evaluator.calls = 0;
    const auto local_required = solve_with(observations);
    REQUIRE_FALSE(local_required);
    CHECK(evaluator.calls == 2);
    CHECK(local_required.error().code == core::FailureCode::PrecisionUnresolved);
    CHECK(total_only->certificate().consumer_count() == 1);
}
