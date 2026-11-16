// Fixed product extraction output in exact hexadecimal form (QA-03).
#include <array>
#include <exception>
#include <format>
#include <iostream>
#include <stdexcept>
#include <uvlas/core/build_info.hpp>
#include <uvlas/numerics/eigen_generalized.hpp>

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
            0.};
}
void emit(std::string_view name, const core::ErrorCertificate& proof) {
    const auto* value = std::get_if<core::Achieved>(&proof.assessment());
    if (value == nullptr)
        throw std::runtime_error("fixed generalized candidate was not qualified");
    std::cout << std::format("{} {:a}\n", name, value->upper_bound);
}
void evaluate() {
    constexpr double delta = 0x1p-60;
    constexpr auto op = core::Operation::ExtractEigen;
    core::ResourcePlan input(10000), plan(1000000);
    const std::array stiffness{-delta, 0., 0., 0., 0., 0., 0., 0., 16 * delta};
    const std::array mass{1., 0., 0., 0., 4., 0., 0., 0., 16.};
    auto k = numerics::DenseReal::copy({3, 3, 3}, stiffness, input, op);
    auto m = numerics::DenseReal::copy({3, 3, 3}, mass, input, op);
    if (!k || !m)
        throw std::runtime_error("fixed input allocation failed");
    const std::array d{1., 1., 1.};
    const std::array<numerics::linear::UnitDimensions, 3> units{};
    const eig::GeneralizedRequest request{
        {{d, 1.}, 1., units, {}, {}, context()},
        {{1.e-12, core::ObservationId{20}, core::ScopeId{21}},
         {1.e-12, core::ObservationId{22}, core::ScopeId{23}},
         {1.e-12, core::ObservationId{24}, core::ScopeId{25}},
         eig::NotRequested{}},
        {core::Precision::Binary64,
         numerics::linear::RetainedPrecision{core::NumericalGeneration{40}},
         {},
         nullptr},
        {3, 3, {3, 0, 0, 1}, numerics::linear::EstimateAdmission{1000000}},
        {},
        false};
    auto result =
        eig::symmetric_generalized_spectrum({*k, context()}, {*m, context()}, request, plan);
    if (!result)
        throw std::runtime_error("fixed generalized extraction failed");
    for (std::size_t i = 0; i < 3; ++i) {
        std::cout << std::format("root.{} {:a}\n", i, result->eigenvalues()[i]);
        std::cout << std::format("basis.diagonal.{} {:a}\n", i,
                                 result->basis().values()[i * 3 + i]);
    }
    emit("residual", result->metrics().residual);
    emit("normalization", result->metrics().normalization);
    emit("orthogonality", result->metrics().orthogonality);
}
int run() {
    try {
        std::cout << "platform " << build_info().platform_key << '\n';
        evaluate();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
} // namespace
int main() {
    try {
        return run();
    } catch (...) {
        return 1;
    }
}
