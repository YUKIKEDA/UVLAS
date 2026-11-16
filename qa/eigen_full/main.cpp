// Product complete-spectrum and derived-result baseline, with no rounded decimal values.
#include <array>
#include <exception>
#include <format>
#include <iostream>
#include <stdexcept>
#include <uvlas/core/build_info.hpp>
#include <uvlas/numerics/eigen_full.hpp>
#include <variant>
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
[[nodiscard]] eig::FullRequest request(core::NumericalGeneration gen) {
    return {{{1e-12, core::ObservationId{1}, core::ScopeId{1}},
             {1e-12, core::ObservationId{2}, core::ScopeId{1}},
             eig::NotRequested{}},
            numerics::linear::RetainedPrecision{gen},
            {3, 3, {32, 32, 0, 2}, numerics::linear::EstimateAdmission{std::uint64_t{1024} * 1024}},
            {},
            nullptr,
            false};
}
[[nodiscard]] double bound(const core::ErrorCertificate& certificate) {
    const auto* a = std::get_if<core::Achieved>(&certificate.assessment());
    if (!a)
        throw std::runtime_error("fixed result not qualified");
    return a->upper_bound;
}
void emit(const char* name, double value) {
    std::cout << std::format("{} {:a}\n", name, value);
}
void evaluate() {
    core::ResourcePlan plan{std::uint64_t{1024} * 1024};
    const std::array<double, 9> a{-0x1p-60, 0, 0, 0, 2, 0, 0, 0, 3};
    auto matrix = numerics::DenseReal::copy({3, 3, 3}, a, plan, core::Operation::ExtractEigen);
    if (!matrix)
        throw std::runtime_error("fixed input");
    auto input = numerics::PrecisionBuffer<double>::binary64(std::move(*matrix),
                                                             core::NumericalGeneration{5});
    if (!input)
        throw std::runtime_error("fixed precision");
    const std::array d{1., 1., 1.};
    const std::array<numerics::linear::UnitDimensions, 3> units{};
    const std::array<numerics::linear::UnitDimensions, 9> au{};
    auto full = eig::full_spectrum(*input, context(), {}, {d, 1, units, {}, au},
                                   request(core::NumericalGeneration{40}), plan);
    if (!full)
        throw std::runtime_error("fixed spectrum");
    const auto view = full->view();
    if (!view || !full->basis_proof())
        throw std::runtime_error("fixed spectrum proof");
    emit("root.negative", view->eigenvalues[0].real());
    emit("root.two", view->eigenvalues[1].real());
    emit("root.three", view->eigenvalues[2].real());
    emit("schur.residual", bound(full->basis_proof()->residual));
    emit("schur.gram", bound(full->basis_proof()->orthogonality));
    const std::array<core::Index, 1> indices{1};
    auto selection = eig::FullSelection::make(*full, indices, plan);
    if (!selection)
        throw std::runtime_error("fixed selection");
    eig::FullPairRequest const pr{request(core::NumericalGeneration{41}).output(),
                                  {{1e-12, core::ObservationId{3}, core::ScopeId{1}},
                                   {1e-12, core::ObservationId{4}, core::ScopeId{1}},
                                   eig::NotRequested{}},
                                  {0, 0},
                                  true,
                                  false};
    auto pairs = eig::full_pairs(*full, *selection, pr, plan);
    if (!pairs)
        throw std::runtime_error("fixed pairs");
    const auto proof = pairs->pair_proof(0);
    if (!proof)
        throw std::runtime_error("fixed pair proof");
    emit("pair.normalization", bound(proof->right.normalization));
    const auto pv = pairs->view();
    if (!pv)
        throw std::runtime_error("fixed pair view");
    emit("pair.component", pv->right[1].real());
    eig::FullClusterRequest const cr{request(core::NumericalGeneration{42}), {0, 0}, true};
    auto cluster = eig::full_cluster(*full, *selection, cr, plan);
    if (!cluster || !cluster->left_basis_proof())
        throw std::runtime_error("fixed cluster");
    emit("cluster.left_residual", bound(cluster->left_basis_proof()->residual));
}
int run() {
    try {
        std::cout << "platform " << uvlas::build_info().platform_key << '\n';
        evaluate();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
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
