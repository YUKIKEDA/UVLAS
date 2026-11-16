// Fixed actual action-Arnoldi results, printed without decimal rounding (QA-03).
#include <array>
#include <complex>
#include <exception>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <uvlas/core/build_info.hpp>
#include <uvlas/numerics/eigen_arnoldi.hpp>

namespace {
using namespace uvlas;
namespace eig = numerics::eigen;
[[nodiscard]] eig::Context context() {
    return {{core::ModelFingerprint{}, core::ModelGeneration{1}, core::GeometryGeneration{2},
             core::ChartGeneration{3}, core::SparsityGeneration{4}, core::NumericalGeneration{5},
             core::ScaleGeneration{6}, core::PrecisionGeneration{7}, core::Precision::Binary64,
             core::WakePolicyGeneration{8}, core::TimeIntervalGeneration{9},
             core::BackendGeneration{10}},
            core::TrialId{11},
            core::ScopeId{12},
            core::ScopeId{13},
            0};
}
class Rotation final : public eig::ActionProvider<double> {
public:
    [[nodiscard]] core::Index order() const noexcept override {
        return 8;
    }
    [[nodiscard]] eig::ActionStamp stamp() const noexcept override {
        return {context(), core::NumericalGeneration{21}, {}};
    }
    [[nodiscard]] numerics::SupportedPrecision
    supported_precision(eig::ActionKind) const noexcept override {
        return numerics::SupportedPrecision::Binary64Only;
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        return {};
    }
    [[nodiscard]] core::Result<void> apply(eig::ActionKind kind, eig::WordView<double> x,
                                           std::span<double> high, std::span<double> low,
                                           std::span<double> error) const override {
        constexpr std::array<double, 8> diagonal{0, 0, .5, .25, .125, .0625, .03125, .015625};
        for (std::size_t i = 0; i < 8; ++i) {
            const auto source = i == 0 ? 1 : (i == 1 ? 0 : i);
            double coefficient = i == 0 ? -4 : (i == 1 ? 4 : diagonal[i]);
            if (i < 2 && kind == eig::ActionKind::Adjoint) {
                coefficient = -coefficient;
            }
            const auto product = numerics::two_product(coefficient, x.high()[source]);
            if (!product) {
                return std::unexpected(product.error());
            }
            high[i] = product->hi;
            low[i] = product->lo;
            error[i] = 0;
        }
        return {};
    }
};
class NoConsumer final : public eig::ArnoldiConsumerEvaluator {
public:
    [[nodiscard]] numerics::SupportedPrecision supported_precision() const noexcept override {
        return numerics::SupportedPrecision::Binary64Only;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const eig::ArnoldiCandidateView&, const eig::ArnoldiCandidateView&,
           const numerics::linear::ConsumerCheck&, core::CertificateContext) const override {
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::ExtractEigen});
    }
};
void emit(std::string_view name, const core::ErrorCertificate& proof) {
    const auto* achieved = std::get_if<core::Achieved>(&proof.assessment());
    if (achieved == nullptr) {
        throw std::runtime_error("unqualified fixed Arnoldi result");
    }
    std::cout << std::format("{} {:a}\n", name, achieved->upper_bound);
}
void evaluate() {
    const Rotation provider;
    const NoConsumer evaluator;
    const std::array<double, 8> seed{1, 2, 3, 4, 5, 6, 7, 8}, components{1, 1, 1, 1, 1, 1, 1, 1};
    const std::array<numerics::linear::UnitDimensions, 8> units{};
    const auto initial = eig::WordView<double>::make(seed, {}, core::Precision::Binary64);
    if (!initial) {
        throw std::runtime_error("invalid fixed Arnoldi seed");
    }
    const eig::ArnoldiPairRequest qualification{
        {{1e-12, core::ObservationId{31}, core::ScopeId{32}},
         {1e-12, core::ObservationId{33}, core::ScopeId{32}},
         eig::NotRequested{}}};
    const eig::ArnoldiRequest<double> request{
        2,
        8,
        *initial,
        1e-14,
        100,
        eig::ArnoldiLargestMagnitude{},
        qualification,
        eig::ArnoldiSides::RightOnly,
        {1e-12, 1e-12},
        numerics::linear::RetainedPrecision{core::NumericalGeneration{40}},
        {},
        {},
        evaluator,
        eig::ArnoldiTimingMode::Disabled};
    const eig::ArnoldiScales scales{{components, 1}, units, {}, context()};
    const eig::ArnoldiLimits limits{
        8, 3, 8, 100, {1000, 0, 0, 1}, numerics::linear::EstimateAdmission{10000000}};
    core::ResourcePlan scratch_plan(1000000), result_plan(1000000);
    auto workspace = eig::ArnoldiWorkspace::create(8, 3, scratch_plan);
    if (!workspace) {
        throw std::runtime_error("fixed Arnoldi workspace failed");
    }
    const auto result =
        eig::arnoldi_spectrum(provider, scales, request, limits, *workspace, result_plan);
    if (!result) {
        throw std::runtime_error("fixed Arnoldi extraction failed");
    }
    const auto values = result->values();
    if (!values || values->eigenvalues.values.size() != 2) {
        throw std::runtime_error("fixed Arnoldi result shape changed");
    }
    for (core::Index i = 0; i < 2; ++i) {
        const auto root = values->eigenvalues.values[static_cast<std::size_t>(i)];
        std::cout << std::format("root{}.real {:a}\nroot{}.imag {:a}\n", i, root.real(), i,
                                 root.imag());
        const auto proof = result->pair_metric(eig::ArnoldiSide::Right, i);
        if (!proof) {
            throw std::runtime_error("fixed Arnoldi proof missing");
        }
        emit(std::format("root{}.residual", i), proof->residual);
        emit(std::format("root{}.normalization", i), proof->normalization);
    }
}
} // namespace
int main() {
    try {
        std::cout << "platform " << uvlas::build_info().platform_key << '\n';
        evaluate();
        return 0;
    } catch (const std::exception& error) {
        try {
            std::cerr << error.what() << '\n';
        } catch (...) {
            return 1;
        }
        return 1;
    } catch (...) {
        return 1;
    }
}
