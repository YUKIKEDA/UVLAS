#include <algorithm>
#include <benchmark/benchmark.h>
#include <cmath>
#include <limits>
#include <uvlas/numerics/gmres_solve.hpp>
#include <vector>
namespace {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;
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
// Synthetic tridiagonal capacity baseline; it represents no physical coupled case.
class Tridiagonal final : public GmresOperator {
public:
    explicit Tridiagonal(core::Index n) : n_(n) {}
    [[nodiscard]] core::Index order() const noexcept override {
        return n_;
    }
    [[nodiscard]] Context context() const noexcept override {
        return identity();
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        return {};
    }
    [[nodiscard]] core::Result<void> apply_word(std::span<const double> x, std::span<double> out,
                                                std::span<double> tail,
                                                std::span<double> error) const override {
        for (std::size_t i = 0; i < x.size(); ++i) {
            const double left = i ? x[i - 1] * 0.25 : 0,
                         right = i + 1 < x.size() ? x[i + 1] * 0.25 : 0;
            out[i] = (x[i] - left) - right;
            tail[i] = 0;
            error[i] = 8 * std::numeric_limits<double>::epsilon() *
                       (std::abs(x[i]) + std::abs(left) + std::abs(right));
        }
        return {};
    }
    [[nodiscard]] core::Result<void> absolute_action_lower(std::span<const double> x,
                                                           std::span<double> out) const override {
        for (std::size_t i = 0; i < x.size(); ++i) {
            const double left = i ? std::nextafter(x[i - 1] * 0.25, 0.0) : 0,
                         right = i + 1 < x.size() ? std::nextafter(x[i + 1] * 0.25, 0.0) : 0;
            out[i] = std::nextafter(std::nextafter(x[i] + left, 0.0) + right, 0.0);
        }
        return {};
    }

private:
    core::Index n_;
};
class IdentityP final : public FixedPreconditioner {
public:
    explicit IdentityP(core::Index n) : n_(n) {}
    [[nodiscard]] core::Index order() const noexcept override {
        return n_;
    }
    [[nodiscard]] FixedPreconditionerIdentity identity() const noexcept override {
        return {::identity(), ::identity().cache};
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        return {};
    }
    [[nodiscard]] core::Result<void> apply(std::span<const double> r,
                                           std::span<double> x) const override {
        std::copy(r.begin(), r.end(), x.begin());
        return {};
    }

private:
    core::Index n_;
};
class NoChecks final : public ConsumerEvaluator<double> {
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
void gmres(benchmark::State& state, bool prepare_only) {
    const auto n = state.range(0);
    const Tridiagonal a{n};
    const IdentityP p{n};
    const std::vector<double> scale(static_cast<std::size_t>(n), 1);
    const std::vector<UnitDimensions> units(static_cast<std::size_t>(n));
    const Scales scales{scale, scale, units, units, identity()};
    const GmresLimits limits{n, 1, EstimateAdmission{1 << 26}};
    core::ResourcePlan input{1 << 26};
    auto matrix = DenseReal::copy({n, 1, n}, scale, input, core::Operation::SolveLinear);
    if (!matrix) {
        state.SkipWithError("input allocation failed");
        return;
    }
    auto rhs = PrecisionBuffer<double>::binary64(std::move(*matrix), core::NumericalGeneration{1});
    if (!rhs) {
        state.SkipWithError("RHS failed");
        return;
    }
    const NoChecks consumer;
    const Precision<double> precision{
        RetainedPrecision{core::NumericalGeneration{2}}, {}, consumer};
    const Accuracy accuracy{Norm::Infinity,        1e-11, 0, 1e-11, core::ObservationId{1},
                            core::ObservationId{2}};
    auto prepared = GmresPrepared::prepare(a, p, scales, limits, input);
    if (!prepared) {
        state.SkipWithError("prepare failed");
        return;
    }
    for (auto iteration : state) {
        (void)iteration;
        core::ResourcePlan payload{1 << 26};
        if (prepare_only) {
            auto result = GmresPrepared::prepare(a, p, scales, limits, payload);
            if (!result) {
                state.SkipWithError("prepare failed");
                break;
            }
            benchmark::DoNotOptimize(result->active());
        } else {
            auto result = solve_linear(a, p, Rhs<double>{*rhs, identity(), core::RhsGeneration{1}},
                                       scales, precision, accuracy, limits, *prepared, payload);
            if (!result) {
                state.SkipWithError("solve failed");
                break;
            }
            benchmark::DoNotOptimize(result->diagnostics().operator_actions);
            state.counters["operator_actions"] =
                static_cast<double>(result->diagnostics().operator_actions);
        }
    }
    state.counters["order"] = static_cast<double>(n);
    state.counters["restart"] = 64;
    state.counters["operator_quota"] = 512;
    state.counters["rhs_columns"] = 1;
}
BENCHMARK_CAPTURE(gmres, Preparation, true)->Arg(96)->Arg(1024)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(gmres, Solve, false)->Arg(96)->Arg(1024)->Unit(benchmark::kMillisecond);
} // namespace
