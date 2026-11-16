// Fixed product qualification results, recorded without decimal rounding (QA-03).
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <uvlas/core/build_info.hpp>
#include <uvlas/numerics/eigen_metrics.hpp>
#include <variant>

namespace {
using namespace uvlas;
namespace eig = numerics::eigen;

[[nodiscard]] eig::MetricContext context() {
    return {{{core::ModelFingerprint{}, core::ModelGeneration{1}, core::GeometryGeneration{2},
              core::ChartGeneration{3}, core::SparsityGeneration{4}, core::NumericalGeneration{5},
              core::ScaleGeneration{6}, core::PrecisionGeneration{7}, core::Precision::Binary64,
              core::WakePolicyGeneration{8}, core::TimeIntervalGeneration{9},
              core::BackendGeneration{10}},
             core::TrialId{11},
             core::ScopeId{12},
             core::ScopeId{13},
             0.0},
            {},
            core::NumericalGeneration{14}};
}
[[nodiscard]] eig::MetricRequirement requirement(std::uint64_t id) {
    return {1e-12, core::ObservationId{id}, core::ScopeId{20}};
}
template <std::size_t N>
[[nodiscard]] eig::WordView<double> view(const std::array<double, N>& values) {
    auto result = eig::WordView<double>::make(values, {}, core::Precision::Binary64);
    if (!result)
        throw std::runtime_error("invalid fixed words");
    return *result;
}
void emit(std::string_view name, const core::ErrorCertificate& certificate) {
    const auto* value = std::get_if<core::Achieved>(&certificate.assessment());
    if (value == nullptr)
        throw std::runtime_error("unqualified fixed candidate");
    std::cout << std::format("{} {:a}\n", name, value->upper_bound);
}
void evaluate() {
    const auto identity = context();
    const std::array scales{1., 1.}, zero{0., 0.};
    const std::array vector{1., std::ldexp(1., -80)};
    const std::array lambda{5.};
    const std::array applied{5., std::ldexp(1., -80)};
    const eig::PairMetricRequest pair_request{requirement(1), requirement(2), eig::NotRequested{}};
    auto pair =
        eig::assess_pair_metrics(view(vector), view(lambda), {view(applied), zero}, {scales, 1.},
                                 eig::PairSide::Right, pair_request, identity);
    if (!pair)
        throw std::runtime_error("fixed pair qualification failed");
    emit("pair.residual", pair->residual);
    emit("pair.normalization", pair->normalization);
    std::cout << std::format("pair.norm_lower {:a}\n", pair->vector_infinity_lower);
    const eig::PairMetricRequest strict_request{requirement(1), requirement(2), requirement(3)};
    auto refused =
        eig::assess_pair_metrics(view(vector), view(lambda), {view(applied), zero}, {scales, 1.},
                                 eig::PairSide::Right, strict_request, identity);
    if (refused || refused.error().cause.code != core::FailureCode::PrecisionUnresolved)
        throw std::runtime_error("componentwise negative control was accepted");

    const std::array basis{1., 0., 0., 1.}, triangular{2., 0., 0., 3.};
    const std::array errors{0., 0., 0., 0.};
    const eig::BasisMetricRequest basis_request{requirement(4), requirement(5),
                                                eig::NotRequested{}};
    auto schur = eig::assess_schur_metrics(2, 2, view(basis), {view(triangular), errors},
                                           {view(triangular), errors}, basis_request, identity);
    if (!schur)
        throw std::runtime_error("fixed Schur qualification failed");
    emit("schur.residual", schur->residual);
    emit("schur.orthogonality", schur->orthogonality);

    // M=diag(1,4), K=diag(-1,16): the physical modes are diag(1,1/2).
    const std::array modes{1., 0., 0., .5};
    const std::array roots{-1., 4.};
    const std::array stiffness{-1., 0., 0., 8.}, mass{1., 0., 0., 2.};
    const eig::GeneralizedMetricRequest generalized_request{requirement(6), requirement(7),
                                                            requirement(8), eig::NotRequested{}};
    auto generalized = eig::assess_generalized_metrics(
        2, view(modes), view(roots), {view(stiffness), errors}, {view(mass), errors}, {scales, 1.},
        1., generalized_request, identity);
    if (!generalized)
        throw std::runtime_error("fixed generalized qualification failed");
    emit("generalized.residual", generalized->residual);
    emit("generalized.normalization", generalized->normalization);
    emit("generalized.orthogonality", generalized->orthogonality);
}
int run() {
    try {
        std::cout << "platform " << uvlas::build_info().platform_key << '\n';
        evaluate();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
} // namespace

// Reporting a failure can itself throw; no exception may escape the process boundary.
int main() {
    try {
        return run();
    } catch (...) {
        return 1;
    }
}
