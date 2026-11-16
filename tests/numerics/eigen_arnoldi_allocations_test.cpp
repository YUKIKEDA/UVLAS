#include "../support/allocation_audit.hpp"
#include "eigen_arnoldi_engine.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <complex>
#include <stdexcept>

// The independent high-precision oracle lives in another executable. This target
// links only recompiled audited product TUs and known header-only dependencies.
namespace {
using namespace uvlas;
namespace eig = numerics::eigen;
using Complex = std::complex<double>;
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
template <class Scalar>
class Diagonal final : public eig::ActionProvider<Scalar> {
public:
    bool inverse = false;
    [[nodiscard]] core::Index order() const noexcept override {
        return 4;
    }
    [[nodiscard]] eig::ActionStamp stamp() const noexcept override {
        return {context(), core::NumericalGeneration{inverse ? 60U : 21U}, {}};
    }
    [[nodiscard]] numerics::SupportedPrecision
    supported_precision(eig::ActionKind) const noexcept override {
        return numerics::SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        return {};
    }
    [[nodiscard]] core::Result<void> apply(eig::ActionKind, eig::WordView<Scalar> input,
                                           std::span<Scalar> high, std::span<Scalar> low,
                                           std::span<double> error) const override {
        const auto coefficient = inverse ? .5 : 2.0;
        for (std::size_t i = 0; i < input.size(); ++i) {
            high[i] = coefficient * input.high()[i];
            low[i] = input.low().empty() ? Scalar{} : coefficient * input.low()[i];
            error[i] = 0;
        }
        return {};
    }
};
class Inverse final : public eig::ArnoldiInverseAccess {
public:
    Diagonal<Complex> action;
    Inverse() {
        action.inverse = true;
    }
    [[nodiscard]] const eig::ActionProvider<Complex>& provider(eig::ArnoldiSide) const override {
        return action;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    qualify(eig::ArnoldiSide, eig::WordView<Complex>, eig::WordView<Complex>,
            std::span<const double>, const eig::ActionStamp&, eig::MetricRequirement request,
            core::CertificateContext expected) const override {
        return core::ErrorCertificate::assess(expected, request.budget, 0);
    }
};
class NoConsumer final : public eig::ArnoldiConsumerEvaluator {
public:
    [[nodiscard]] numerics::SupportedPrecision supported_precision() const noexcept override {
        return numerics::SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const eig::ArnoldiCandidateView&, const eig::ArnoldiCandidateView&,
           const numerics::linear::ConsumerCheck&, core::CertificateContext) const override {
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::ExtractEigen});
    }
};

template <class Original, class Native, bool Shift>
void verify_adapter() {
    const Diagonal<Original> provider;
    const Inverse inverse;
    const NoConsumer consumer;
    const std::array<Original, 4> seed{Original{1}, Original{2}, Original{3}, Original{4}};
    const std::array<double, 4> scales{1, 2, .5, 4};
    const std::array<numerics::linear::UnitDimensions, 4> units{};
    const auto words = eig::WordView<Original>::make(seed, {}, core::Precision::Binary64);
    REQUIRE(words.has_value());
    const eig::ArnoldiPairRequest criteria{{{1e-12, core::ObservationId{31}, core::ScopeId{32}},
                                            {1e-12, core::ObservationId{33}, core::ScopeId{32}},
                                            eig::NotRequested{}}};
    const auto selection = [&]() -> eig::ArnoldiSelection {
        if constexpr (Shift) {
            return eig::ArnoldiShiftInvert{{},
                                           inverse,
                                           {1e-12, core::ObservationId{50}, core::ScopeId{51}},
                                           core::NumericalGeneration{1000}};
        } else {
            return eig::ArnoldiLargestMagnitude{};
        }
    }();
    const auto metadata = [&]() -> eig::ArnoldiSelectionMetadata {
        if constexpr (Shift) {
            return eig::ArnoldiShiftMetadata{{},
                                             inverse.action.stamp(),
                                             eig::NotRequested{},
                                             {1e-12, core::ObservationId{50}, core::ScopeId{51}},
                                             core::NumericalGeneration{1000},
                                             8};
        } else {
            return eig::ArnoldiLargestMagnitude{};
        }
    }();
    const eig::ArnoldiRequest<Original> wanted{
        1,
        3,
        *words,
        1e-14,
        100,
        selection,
        criteria,
        eig::ArnoldiSides::RightOnly,
        {1e-12, 1e-12},
        numerics::linear::RetainedPrecision{core::NumericalGeneration{40}},
        {},
        {},
        consumer,
        eig::ArnoldiTimingMode::Disabled};
    const eig::ArnoldiScales fixed{{scales, 1}, units, {}, context()};
    const eig::ArnoldiLimits limits{
        4, 2, 3, 100, {8, 8, 8, 1}, numerics::linear::EstimateAdmission{10000000}};
    eig::detail::ArnoldiData data{
        context(),
        provider.stamp(),
        {},
        wanted.returned,
        wanted.association,
        metadata,
        {{1, 2, 0},
         {5, 0, 1},
         std::same_as<Original, double> ? eig::ArnoldiDomain::Real : eig::ArnoldiDomain::Complex,
         std::same_as<Native, double> ? eig::ArnoldiDomain::Real : eig::ArnoldiDomain::Complex},
        {},
        1,
        false,
        false,
        {1, 3, 1e-14, 100, wanted.sides, wanted.qualification, wanted.timing},
        limits};
    data.order = 4;
    core::ResourcePlan scratch_plan(1000000), result_plan(1000000);
    auto workspace = eig::ArnoldiWorkspace::create(4, 2, scratch_plan);
    REQUIRE(workspace.has_value());
    auto* scratch = eig::detail::ArnoldiWorkspaceAccess::data(*workspace);
    REQUIRE(scratch != nullptr);
    auto session = eig::ActionSession<Original>::make(provider, limits.actions);
    REQUIRE(session.has_value());
    eig::detail::ArnoldiRunner<Original> runner{
        provider, fixed, wanted, limits, data, std::move(*session), *scratch};
    if constexpr (Shift) {
        auto actions = eig::ActionSession<Complex>::make(inverse.action, limits.actions);
        REQUIRE(actions.has_value());
        runner.inverse_right.emplace(std::move(*actions));
        auto proofs = result_plan.allocate<std::optional<eig::ArnoldiInverseAssessment>>(
            8, core::Operation::ExtractEigen);
        REQUIRE(proofs.has_value());
        data.inverse_assessments = std::move(*proofs);
    }
    const eig::detail::ArnoldiOperation<Original, Native> operation{runner,
                                                                    eig::ArnoldiSide::Right};
    const std::array<Native, 4> input{Native{1}, Native{2}, Native{3}, Native{4}};
    std::array<Native, 4> output{};
    test_support::AllocationAudit audit(test_support::strict_kernel_backends,
                                        test_support::AuditThreadDomain::CallingThread);
    operation.perform_op(input.data(), output.data());
    const auto result = audit.finish();
    REQUIRE_FALSE(runner.first.has_value());
    CHECK(result.status == test_support::AllocationAuditStatus::Clean);
    CHECK(result.cpp_allocation_attempts == 0);
    for (std::size_t i = 0; i < input.size(); ++i) {
        CHECK(output[i] == (Shift ? .5 : 2.0) * input[i]);
    }
}
} // namespace
TEST_CASE("Arnoldi production real complex and shifted adapters allocate no memory",
          "[N4-07][EIG-11][negative][verification]") {
    verify_adapter<double, double, false>();
    verify_adapter<Complex, Complex, false>();
    verify_adapter<double, Complex, true>();
}
