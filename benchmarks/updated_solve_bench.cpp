#include <array>
#include <benchmark/benchmark.h>
#include <complex>
#include <utility>
#include <uvlas/numerics/dense_solve.hpp>
#include <vector>

namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;
using Complex = std::complex<double>;
constexpr auto operation = core::Operation::SolveLinear;

class AlgebraicOnly final : public ConsumerEvaluator<Complex> {
public:
    [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
        return SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<Complex>&, const PrecisionBuffer<Complex>&, const Context&,
           const ConsumerCheck& check, core::CertificateContext expected) const override {
        return core::ErrorCertificate::unevaluated(expected, check.budget);
    }
};
[[nodiscard]] Context identity() {
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
void updated(benchmark::State& state, bool prepare_only) {
    const core::Index n = state.range(0);
    const auto count = core::checked_extent_product(std::array<core::Index, 2>{n, n}, operation);
    if (!count || n <= 0) {
        state.SkipWithError("invalid fixture dimensions");
        return;
    }
    core::ResourcePlan inputs{1 << 24};
    std::vector<double> coefficients(static_cast<std::size_t>(*count), 0);
    std::vector<double> scales_storage(static_cast<std::size_t>(n), 1);
    std::vector<UnitDimensions> units(static_cast<std::size_t>(n));
    // A0=4I, U=i/32 and V=(1+i)/16 give a dense, non-conjugated rank-one update.
    std::vector<Complex> u_storage(static_cast<std::size_t>(n), Complex{0, 1.0 / 32});
    std::vector<Complex> v_storage(static_cast<std::size_t>(n), Complex{1.0 / 16, 1.0 / 16});
    std::vector<Complex> rhs_storage(static_cast<std::size_t>(n), Complex{1, 1.0 / 8});
    for (core::Index i = 0; i < n; ++i) {
        coefficients[static_cast<std::size_t>(i * n + i)] = 4;
    }
    auto base = DenseReal::copy({n, n, n}, coefficients, inputs, operation);
    auto u = DenseComplex::copy({n, 1, n}, u_storage, inputs, operation);
    auto v = DenseComplex::copy({n, 1, n}, v_storage, inputs, operation);
    auto b = DenseComplex::copy({n, 1, n}, rhs_storage, inputs, operation);
    if (!base || !u || !v || !b) {
        state.SkipWithError("fixture allocation failed");
        return;
    }
    auto rhs = PrecisionBuffer<Complex>::binary64(std::move(*b), core::NumericalGeneration{2});
    if (!rhs) {
        state.SkipWithError("rhs failed");
        return;
    }
    const Scales scales{scales_storage, scales_storage, units, units, identity()};
    const Limits limits{n, 1, 1, 5, EstimateAdmission{1 << 24}};
    const ColumnUpdatedOperator op{*base, *u, *v, identity()};
    auto prepared = UpdatedPrepared::prepare(op, scales, limits, inputs);
    if (!prepared) {
        state.SkipWithError("updated preparation failed");
        return;
    }
    const AlgebraicOnly evaluator;
    const Precision<Complex> precision{
        RetainedPrecision{core::NumericalGeneration{3}}, {}, evaluator};
    const Accuracy accuracy{Norm::Infinity,        1.e-12, 0, 1.e-12, core::ObservationId{1},
                            core::ObservationId{2}};
    unsigned refinements = 0;
    for (auto iteration : state) {
        (void)iteration;
        core::ResourcePlan output{1 << 24};
        if (prepare_only) {
            auto result = UpdatedPrepared::prepare(op, scales, limits, output);
            if (!result) {
                state.SkipWithError("updated preparation failed");
                break;
            }
            benchmark::DoNotOptimize(result);
        } else {
            auto result = solve_linear(op, Rhs{*rhs, identity(), core::RhsGeneration{1}}, scales,
                                       precision, accuracy, limits, *prepared, output);
            if (!result) {
                state.SkipWithError("updated solve failed");
                break;
            }
            refinements = result->certificate().refinements();
            benchmark::DoNotOptimize(result);
        }
    }
    state.counters["order"] = static_cast<double>(n);
    state.counters["rhs_columns"] = 1;
    state.counters["update_rank"] = 1;
    state.counters["refinements"] = static_cast<double>(refinements);
    state.counters["output_allocation_included"] = prepare_only ? 0 : 1;
    state.SetLabel(
        prepare_only
            ? "synthetic capacity baseline; original snapshots, workspace, real LU and complex "
              "Schur"
            : "synthetic capacity baseline; prepared factors, original residual and owning output");
}
void UpdatedPreparation(benchmark::State& state) {
    updated(state, true);
}
void UpdatedSolve(benchmark::State& state) {
    updated(state, false);
}
BENCHMARK(UpdatedPreparation)->Arg(16)->Arg(64);
BENCHMARK(UpdatedSolve)->Arg(16)->Arg(64);
} // namespace
