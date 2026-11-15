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
    CHECK(local_required.error().details.problem == core::ProblemId{core::ObservationId{82}});
    CHECK(local_required.error().details.numerical_generation == core::NumericalGeneration{40});
    CHECK(total_only->certificate().consumer_count() == 1);
}
namespace {
template <class Value>
void check_known_failure(const core::Result<Value>& result) {
    REQUIRE_FALSE(result);
    CHECK(result.error().details.model == core::ModelGeneration{1});
    CHECK(result.error().details.trial == core::TrialId{11});
    CHECK(result.error().details.numerical_generation == core::NumericalGeneration{5});
}
template <StorageScalar Scalar>
class DiagnosticConsumer final : public ConsumerEvaluator<Scalar> {
public:
    core::FailureDetails supplied{};
    [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
        return SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<Scalar>&, const PrecisionBuffer<Scalar>&, const Context&,
           const ConsumerCheck&, core::CertificateContext) const override {
        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved,
                                             core::Operation::SolveLinear, supplied});
    }
};
} // namespace
TEST_CASE("Every dense public boundary completes known failure identities and preserves resource "
          "diagnostics",
          "[negative][N4-05][N4-08]") {
    Inputs input;
    const auto a = input.matrix();
    const auto b = input.rhs();
    core::ResourcePlan none{0};
    const DenseOperator real_op{a, context()};
    const auto prepare_real = DensePrepared::prepare(real_op, input.scales(), input.limits(), none);
    check_known_failure(prepare_real);
    REQUIRE(prepare_real.error().details.resources);
    CHECK(*prepare_real.error().details.resources ==
          core::ResourceUsage{core::ResourceKind::Bytes, 4 * sizeof(double), 0});
    auto real_prepared = input.prepare(a);
    const auto solve_real = solve_linear(
        real_op, Rhs{b, context(), core::RhsGeneration{1}}, input.scales(),
        Precision<double>{RetainedPrecision{core::NumericalGeneration{40}}, {}, input.consumer},
        input.accuracy(), input.limits(), real_prepared, none);
    check_known_failure(solve_real);
    REQUIRE(solve_real.error().details.resources);
    CHECK(*solve_real.error().details.resources ==
          core::ResourceUsage{core::ResourceKind::Bytes, sizeof(ColumnMetrics), 0});
    using Complex = std::complex<double>;
    const std::array<Complex, 2> zero{Complex{}, Complex{}}, load{Complex{4}, Complex{6}};
    auto u = DenseComplex::copy({2, 1, 2}, zero, input.source_plan, core::Operation::SolveLinear);
    auto v = DenseComplex::copy({2, 1, 2}, zero, input.source_plan, core::Operation::SolveLinear);
    auto raw_rhs =
        DenseComplex::copy({2, 1, 2}, load, input.source_plan, core::Operation::SolveLinear);
    REQUIRE(u);
    REQUIRE(v);
    REQUIRE(raw_rhs);
    auto complex_rhs =
        PrecisionBuffer<Complex>::binary64(std::move(*raw_rhs), core::NumericalGeneration{30});
    REQUIRE(complex_rhs);
    const ColumnUpdatedOperator updated{a, *u, *v, context()};
    const auto prepare_complex =
        UpdatedPrepared::prepare(updated, input.scales(), input.limits(), none);
    check_known_failure(prepare_complex);
    REQUIRE(prepare_complex.error().details.resources);
    CHECK(*prepare_complex.error().details.resources ==
          core::ResourceUsage{core::ResourceKind::Bytes, 4 * sizeof(double), 0});
    auto complex_prepared =
        UpdatedPrepared::prepare(updated, input.scales(), input.limits(), input.workspace_plan);
    REQUIRE(complex_prepared);
    const DiagnosticConsumer<Complex> complex_consumer;
    const auto solve_complex = solve_linear(
        updated, Rhs{*complex_rhs, context(), core::RhsGeneration{1}}, input.scales(),
        Precision<Complex>{RetainedPrecision{core::NumericalGeneration{40}}, {}, complex_consumer},
        input.accuracy(), input.limits(), *complex_prepared, none);
    check_known_failure(solve_complex);
    REQUIRE(solve_complex.error().details.resources);
    CHECK(*solve_complex.error().details.resources ==
          core::ResourceUsage{core::ResourceKind::Bytes, sizeof(ColumnMetrics), 0});
    CHECK(none.allocated_bytes() == 0);
    // Consumer failures can identify a different numerical observation and explicit resource unit.
    for (const bool all_known : {false, true}) {
        core::FailureDetails supplied{};
        supplied.numerical_generation = core::NumericalGeneration{333};
        supplied.resources =
            core::ResourceUsage{core::ResourceKind::OperatorApplications, 501, 101};
        supplied.parent = core::AcceptedStateId{44};
        supplied.problem = core::ScopeId{61};
        if (all_known) {
            supplied.model = core::ModelGeneration{77};
            supplied.trial = core::TrialId{88};
        }
        auto expected = supplied;
        if (!all_known) {
            expected.model = context().cache.model;
            expected.trial = context().trial;
        }
        DiagnosticConsumer<double> evaluator;
        evaluator.supplied = supplied;
        const auto real_failure =
            solve_linear(real_op, Rhs{b, context(), core::RhsGeneration{1}}, input.scales(),
                         Precision<double>{RetainedPrecision{core::NumericalGeneration{40}},
                                           requirement, evaluator},
                         input.accuracy(), input.limits(), real_prepared, input.result_plan);
        REQUIRE_FALSE(real_failure);
        CHECK(real_failure.error().details == expected);
        DiagnosticConsumer<Complex> complex_evaluator;
        complex_evaluator.supplied = supplied;
        const auto complex_failure = solve_linear(
            updated, Rhs{*complex_rhs, context(), core::RhsGeneration{1}}, input.scales(),
            Precision<Complex>{RetainedPrecision{core::NumericalGeneration{40}}, requirement,
                               complex_evaluator},
            input.accuracy(), input.limits(), *complex_prepared, input.result_plan);
        REQUIRE_FALSE(complex_failure);
        CHECK(complex_failure.error().details == expected);
    }
}
