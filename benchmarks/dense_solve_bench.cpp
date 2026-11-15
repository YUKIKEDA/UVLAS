#include <benchmark/benchmark.h>
#include <uvlas/numerics/dense_solve.hpp>
#include <vector>
namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;
constexpr auto operation = core::Operation::SolveLinear;
class AlgebraicOnly final : public ConsumerEvaluator<double> {
public:
    [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
        return SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<double>&, const PrecisionBuffer<double>&, const Context&,
           const ConsumerCheck& check, core::CertificateContext expected) const override {
        return core::ErrorCertificate::unevaluated(expected, check.budget);
    }
};
Context identity() {
    return {{core::ModelFingerprint{}, core::ModelGeneration{1}, core::GeometryGeneration{1},
             core::ChartGeneration{1}, core::SparsityGeneration{1}, core::NumericalGeneration{1},
             core::ScaleGeneration{1}, core::PrecisionGeneration{1}, core::Precision::Binary64,
             core::WakePolicyGeneration{1}, core::TimeIntervalGeneration{1},
             core::BackendGeneration{1}},
            core::TrialId{1},
            core::ScopeId{1},
            core::ScopeId{2},
            0};
}
void dense_solve(benchmark::State& state, bool promote, bool prepare_only = false) {
    const auto n = state.range(0);
    core::ResourcePlan inputs{1 << 24};
    std::vector<double> coefficients(static_cast<std::size_t>(n * n), 0),
        loads(static_cast<std::size_t>(n), 1), scale(static_cast<std::size_t>(n), 1);
    std::vector<UnitDimensions> units(static_cast<std::size_t>(n));
    for (core::Index i = 0; i < n; ++i)
        coefficients[static_cast<std::size_t>(i * n + i)] = promote ? 3 : 4;
    auto matrix = DenseReal::copy({n, n, n}, coefficients, inputs, operation);
    auto b = DenseReal::copy({n, 1, n}, loads, inputs, operation);
    if (!matrix || !b) {
        state.SkipWithError("fixture allocation failed");
        return;
    }
    auto rhs = PrecisionBuffer<double>::binary64(std::move(*b), core::NumericalGeneration{2});
    if (!rhs) {
        state.SkipWithError("rhs failed");
        return;
    }
    const Scales scales{scale, scale, units, units, identity()};
    const Limits limits{n, 1, 0, 5, EstimateAdmission{1 << 24}};
    const DenseOperator op{*matrix, identity()};
    auto prepared = DensePrepared::prepare(op, scales, limits, inputs);
    if (!prepared) {
        state.SkipWithError("preparation failed");
        return;
    }
    const AlgebraicOnly evaluator;
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{3}}, {}, evaluator};
    const Accuracy accuracy{
        Norm::Infinity,          promote ? 1e-30 : 1e-12, 0,
        promote ? 1e-30 : 1e-12, core::ObservationId{1},  core::ObservationId{2}};
    if (prepare_only) {
        for (auto iteration : state) {
            (void)iteration;
            core::ResourcePlan workspace{1 << 24};
            auto candidate = DensePrepared::prepare(op, scales, limits, workspace);
            if (!candidate) {
                state.SkipWithError("preparation failed");
                break;
            }
            benchmark::DoNotOptimize(candidate);
        }
        state.counters["order"] = static_cast<double>(n);
        state.SetLabel(
            "synthetic capacity baseline; snapshots, workspace, LU and condition estimate");
        return;
    }
    for (auto iteration : state) {
        (void)iteration;
        core::ResourcePlan output{1 << 24};
        auto result = solve_linear(op, Rhs{*rhs, identity(), core::RhsGeneration{1}}, scales,
                                   precision, accuracy, limits, *prepared, output);
        if (!result) {
            state.SkipWithError("solve failed");
            break;
        }
        benchmark::DoNotOptimize(result);
    }
    state.counters["order"] = static_cast<double>(n);
    state.counters["rhs_columns"] = 1;
    state.counters["output_allocation_included"] = 1;
    state.counters["two_word"] = promote ? 1 : 0;
    state.SetLabel("synthetic capacity baseline; prepared factor; original residual and owning "
                   "result included");
}
void DenseSolveBinary64(benchmark::State& state) {
    dense_solve(state, false);
}
void DenseSolvePromotedTwoWord(benchmark::State& state) {
    dense_solve(state, true);
}
void DensePreparation(benchmark::State& state) {
    dense_solve(state, false, true);
}
BENCHMARK(DensePreparation)->Arg(16)->Arg(64);
BENCHMARK(DenseSolveBinary64)->Arg(16)->Arg(64);
BENCHMARK(DenseSolvePromotedTwoWord)->Arg(16)->Arg(64);
} // namespace
