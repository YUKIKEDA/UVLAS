#pragma once
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <span>
#include <utility>
#include <uvlas/numerics/sparse_solve.hpp>
#include <vector>

namespace sparse_test {
using namespace uvlas;
using namespace uvlas::numerics;
using namespace uvlas::numerics::linear;
constexpr auto operation = core::Operation::SolveLinear;

template <class T>
T take(core::Result<T> value) {
    if (!value) {
        FAIL("Unexpected failure code " << static_cast<unsigned>(value.error().code));
    }
    return std::move(value).value();
}
inline Context context() {
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
inline Limits limits(core::Index n = 16, core::Index rhs = 4) {
    return {n, rhs, 0, 5, EstimateAdmission{std::uint64_t{512} * 1024 * 1024}};
}
inline Accuracy accuracy(Norm norm = Norm::Infinity) {
    return {norm, 1e-12, 0, 1e-12, core::ObservationId{20}, core::ObservationId{21}};
}
class Evaluator final : public ConsumerEvaluator<double> {
public:
    mutable unsigned calls = 0;
    bool reject = false;
    mutable bool wrong_generation = false;
    [[nodiscard]] SupportedPrecision supported_precision() const noexcept override {
        return SupportedPrecision::Binary64Only;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const PrecisionBuffer<double>&, const PrecisionBuffer<double>& candidate, const Context&,
           const ConsumerCheck& request, core::CertificateContext expected) const override {
        ++calls;
        const auto* binary = std::get_if<Binary64<double>>(&candidate.storage());
        if (binary == nullptr) {
            return std::unexpected(
                core::Failure{core::FailureCode::PrecisionUnsupported, operation});
        }
        double upper = 0;
        const auto shape = binary->values.shape();
        for (core::Index column = 0; column < shape.columns; ++column) {
            const auto start = static_cast<std::size_t>(column * shape.leading_dimension);
            const auto values = binary->values.values();
            if (request.kind == CheckKind::Constraint) {
                upper = std::max(upper, std::abs(values[start] + values[start + 1]));
            } else if (request.kind == CheckKind::Reaction) {
                upper = std::max(
                    upper, std::abs(-values[start + 2] + 2 * static_cast<double>(column + 1)));
            }
        }
        if (reject) {
            upper = request.budget + 1;
        }
        if (wrong_generation) {
            expected.numerical_generation = core::NumericalGeneration{999};
        }
        return core::ErrorCertificate::assess(expected, request.budget, upper);
    }
};
struct Fixture {
    core::ResourcePlan plan{std::uint64_t{512} * 1024 * 1024};
    std::vector<double> unknown, residual;
    std::vector<UnitDimensions> unknown_units, residual_units;
    Evaluator evaluator;
    explicit Fixture(core::Index n)
        : unknown(static_cast<std::size_t>(n), 1),
          residual(static_cast<std::size_t>(n), 1),
          unknown_units(static_cast<std::size_t>(n)),
          residual_units(static_cast<std::size_t>(n)) {}
    Scales scales(Context identity = context()) const {
        return {unknown, residual, unknown_units, residual_units, identity};
    }
    Precision<double> precision(std::span<const ConsumerCheck> checks = {}) const {
        return {RetainedPrecision{core::NumericalGeneration{100}}, checks, evaluator};
    }
    SparseReal matrix(core::Index n, std::span<const Triplet> entries) {
        return take(SparseReal::from_triplets(n, n, entries, plan, operation));
    }
    PrecisionBuffer<double> rhs(core::Index n, core::Index columns,
                                std::span<const double> values) {
        auto matrix = take(DenseReal::copy({n, columns, n}, values, plan, operation));
        return take(
            PrecisionBuffer<double>::binary64(std::move(matrix), core::NumericalGeneration{40}));
    }
};
inline std::span<const double> values(const Solution<double>& result) {
    return std::get<Binary64<double>>(result.values().storage()).values.values();
}
inline const std::array<Triplet, 6> kkt{
    {{0, 0, 2}, {2, 0, 1}, {1, 1, 3}, {2, 1, 1}, {0, 2, 1}, {1, 2, 1}}};
} // namespace sparse_test
