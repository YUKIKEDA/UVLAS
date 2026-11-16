#pragma once

#include "eigen_arnoldi_data.hpp"
#include "eigen_arnoldi_native.hpp"
#include "eigen_metric_bounds.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <type_traits>

namespace uvlas::numerics::eigen::detail {
using ArnoldiCode = core::FailureCode;
inline constexpr auto arnoldi_operation = core::Operation::ExtractEigen;
[[nodiscard]] inline bool arnoldi_finite(ArnoldiComplex z) noexcept {
    return std::isfinite(z.real()) && std::isfinite(z.imag());
}
[[nodiscard]] inline bool arnoldi_certificate(const core::ErrorCertificate& proof,
                                              core::CertificateContext expected,
                                              double budget) noexcept {
    return proof.context() == expected &&
           std::bit_cast<std::uint64_t>(proof.budget()) == std::bit_cast<std::uint64_t>(budget) &&
           std::holds_alternative<core::Achieved>(proof.assessment());
}
template <class T>
[[nodiscard]] core::Result<void> arnoldi_allocate(core::PlannedBuffer<T>& destination,
                                                  core::Index count,
                                                  core::ResourcePlan& plan) noexcept {
    auto allocation = plan.allocate<T>(count, arnoldi_operation);
    if (!allocation) {
        return std::unexpected(allocation.error());
    }
    destination = std::move(*allocation);
    return {};
}

// All spans in this adapter are allocated before Spectra starts. It never creates
// a vector/matrix, grows a container, or keeps a borrowed provider beyond extraction.
template <StorageScalar Original>
struct ArnoldiRunner {
    const ActionProvider<Original>& original;
    ArnoldiScales scales;
    const ArnoldiRequest<Original>& request;
    const ArnoldiLimits& limits;
    ArnoldiData& data;
    ActionSession<Original> session;
    std::optional<ActionSession<ArnoldiComplex>> inverse_right, inverse_left;
    std::optional<EigenFailure> first;
    core::PlannedBuffer<double>& real_input;
    core::PlannedBuffer<double>& real_input_low;
    core::PlannedBuffer<double>& real_high;
    core::PlannedBuffer<double>& real_low;
    core::PlannedBuffer<double>& real_error;
    core::PlannedBuffer<ArnoldiComplex>& input;
    core::PlannedBuffer<ArnoldiComplex>& input_low;
    core::PlannedBuffer<ArnoldiComplex>& high;
    core::PlannedBuffer<ArnoldiComplex>& low;
    core::PlannedBuffer<double>& error;
    core::PlannedBuffer<ArnoldiComplex>& cluster_action;
    core::PlannedBuffer<double>& cluster_error;
    core::PlannedBuffer<double>& triangular_error;
    core::PlannedBuffer<ArnoldiComplex>& candidate_values;
    core::PlannedBuffer<ArnoldiComplex>& candidate_vectors;

    ArnoldiRunner(const ActionProvider<Original>& provider, const ArnoldiScales& fixed,
                  const ArnoldiRequest<Original>& wanted, const ArnoldiLimits& capacities,
                  ArnoldiData& owned, ActionSession<Original>&& actions,
                  ArnoldiWorkspaceData& scratch) noexcept
        : original(provider),
          scales(fixed),
          request(wanted),
          limits(capacities),
          data(owned),
          session(std::move(actions)),
          real_input(scratch.real_input),
          real_input_low(scratch.real_input_low),
          real_high(scratch.real_high),
          real_low(scratch.real_low),
          real_error(scratch.real_error),
          input(scratch.input),
          input_low(scratch.input_low),
          high(scratch.high),
          low(scratch.low),
          error(scratch.error),
          cluster_action(scratch.cluster_action),
          cluster_error(scratch.cluster_error),
          triangular_error(scratch.triangular_error),
          candidate_values(scratch.candidate_values),
          candidate_vectors(scratch.candidate_vectors) {}
    [[nodiscard]] AttemptDiagnostics counters() const noexcept {
        auto result = session.diagnostics();
        result.requested_count = static_cast<std::uint64_t>(request.nev);
        result.converged_count = data.diagnostics.attempts.converged_count;
        result.backend_iterations = data.diagnostics.attempts.backend_iterations;
        result.inverse_applications = 0;
        if (inverse_right) {
            result.inverse_applications = inverse_right->diagnostics().inverse_applications;
        }
        if (inverse_left) {
            result.inverse_applications =
                result.inverse_applications.value_or(0) +
                inverse_left->diagnostics().inverse_applications.value_or(0);
        }
        return result;
    }
    [[nodiscard]] EigenFailure failure(core::Failure cause) const noexcept {
        return complete_failure(cause, scales.context, request.state, counters());
    }
    void record_failure(ArnoldiCode code) noexcept {
        if (!first) {
            first = failure(core::Failure{code, arnoldi_operation});
        }
    }
    [[nodiscard]] ArnoldiResult<void> fail(ArnoldiCode code) noexcept {
        record_failure(code);
        return std::unexpected(*first);
    }
    [[nodiscard]] ArnoldiResult<void> remember(EigenFailure failed) noexcept {
        if (!first) {
            // Session counters describe one provider; translate that internal domain
            // to the whole extraction without altering the producer cause.
            first = failure(std::move(failed.cause));
        }
        return std::unexpected(*first);
    }
    // Real providers receive the two components separately, including a zero
    // component. Those are two actual original/adjoint applications, not one.
    [[nodiscard]] ArnoldiResult<void>
    original_complex(ArnoldiSide side, std::span<const ArnoldiComplex> x,
                     std::span<ArnoldiComplex> hi, std::span<ArnoldiComplex> tail,
                     std::span<double> upper,
                     std::span<const ArnoldiComplex> input_tail = {}) noexcept {
        if (first) {
            return std::unexpected(*first);
        }
        if (!input_tail.empty() &&
            std::ranges::all_of(input_tail, [](auto value) { return value == ArnoldiComplex{}; })) {
            input_tail = {};
        }
        const auto kind = side == ArnoldiSide::Right ? ActionKind::Forward : ActionKind::Adjoint;
        if constexpr (std::same_as<Original, ArnoldiComplex>) {
            auto words = WordView<Original>::make(x, input_tail,
                                                  input_tail.empty() ? core::Precision::Binary64
                                                                     : core::Precision::TwoWord);
            if (!words) {
                return remember(failure(words.error()));
            }
            const auto action = session.apply(original, kind, *words, hi, tail, upper);
            if (!action) {
                return remember(action.error());
            }
        } else {
            for (unsigned component = 0; component < 2; ++component) {
                for (std::size_t i = 0; i < x.size(); ++i) {
                    real_input.values()[i] = component == 0 ? x[i].real() : x[i].imag();
                    if (!input_tail.empty()) {
                        real_input_low.values()[i] =
                            component == 0 ? input_tail[i].real() : input_tail[i].imag();
                    }
                }
                auto words = WordView<double>::make(
                    real_input.values(),
                    input_tail.empty() ? std::span<const double>{} : real_input_low.values(),
                    input_tail.empty() ? core::Precision::Binary64 : core::Precision::TwoWord);
                if (!words) {
                    return remember(failure(words.error()));
                }
                const auto action = session.apply(original, kind, *words, real_high.values(),
                                                  real_low.values(), real_error.values());
                if (!action) {
                    return remember(action.error());
                }
                for (std::size_t i = 0; i < x.size(); ++i) {
                    if (component == 0) {
                        hi[i] = {real_high.values()[i], 0};
                        tail[i] = {real_low.values()[i], 0};
                        upper[i] = real_error.values()[i];
                    } else {
                        hi[i].imag(real_high.values()[i]);
                        tail[i].imag(real_low.values()[i]);
                        const auto sum = upper[i] + real_error.values()[i];
                        upper[i] =
                            sum == 0 ? 0
                                     : std::nextafter(sum, std::numeric_limits<double>::infinity());
                        if (!std::isfinite(upper[i])) {
                            return fail(ArnoldiCode::PrecisionUnresolved);
                        }
                    }
                }
            }
        }
        return {};
    }
    [[nodiscard]] ArnoldiResult<void> original_real(ArnoldiSide side, std::span<const double> x,
                                                    std::span<double> y) noexcept {
        static_assert(std::same_as<Original, double>);
        auto words = WordView<double>::make(x, {}, core::Precision::Binary64);
        if (!words) {
            return remember(failure(words.error()));
        }
        const auto action = session.apply(
            original, side == ArnoldiSide::Right ? ActionKind::Forward : ActionKind::Adjoint,
            *words, real_high.values(), real_low.values(), real_error.values());
        if (!action) {
            return remember(action.error());
        }
        for (std::size_t i = 0; i < x.size(); ++i) {
            y[i] = real_high.values()[i] + real_low.values()[i];
            if (!std::isfinite(y[i])) {
                return fail(ArnoldiCode::PrecisionUnresolved);
            }
        }
        return {};
    }
    [[nodiscard]] ArnoldiResult<void> inverse(ArnoldiSide side) {
        const auto* policy = std::get_if<ArnoldiShiftInvert>(&request.selection);
        if (policy == nullptr) {
            return fail(ArnoldiCode::InvalidInput);
        }
        auto& actions = side == ArnoldiSide::Right ? inverse_right : inverse_left;
        if (!actions) {
            return fail(ArnoldiCode::InvalidInput);
        }
        if (counters().inverse_applications.value_or(0) >= limits.actions.inverse) {
            core::FailureDetails details;
            details.resources = core::ResourceUsage{
                core::ResourceKind::OperatorApplications,
                limits.actions.inverse == std::numeric_limits<std::uint64_t>::max()
                    ? limits.actions.inverse
                    : limits.actions.inverse + 1,
                limits.actions.inverse};
            return remember(
                failure(core::Failure{ArnoldiCode::ResourceExceeded, arnoldi_operation, details}));
        }
        const auto& provider = policy->inverse.provider(side);
        const auto* metadata = std::get_if<ArnoldiShiftMetadata>(&data.selection);
        if (metadata == nullptr) {
            return fail(ArnoldiCode::InvalidInput);
        }
        const auto* left_stamp = std::get_if<ActionStamp>(&metadata->left_inverse);
        const auto expected_stamp =
            side == ArnoldiSide::Right
                ? metadata->right_inverse
                : (left_stamp != nullptr ? *left_stamp : metadata->right_inverse);
        if (provider.stamp() != expected_stamp) {
            return fail(ArnoldiCode::StaleGeneration);
        }
        const auto rhs =
            WordView<ArnoldiComplex>::make(input.values(), {}, core::Precision::Binary64);
        if (!rhs) {
            return remember(failure(rhs.error()));
        }
        const auto action = actions->apply(provider, ActionKind::Inverse, *rhs, high.values(),
                                           low.values(), error.values());
        if (!action) {
            return remember(action.error());
        }
        const auto solution =
            WordView<ArnoldiComplex>::make(high.values(), low.values(), core::Precision::TwoWord);
        if (!solution) {
            return remember(failure(solution.error()));
        }
        // A separate exact-inverse-output enclosure is required from the provider.
        // The linear proof here is necessary but is not substituted for that bound.
        const core::CertificateContext expected{
            policy->linear_accuracy.observation, policy->linear_accuracy.scope,
            core::NumericalGeneration{policy->first_inverse_generation.value() +
                                      counters().inverse_applications.value_or(0) - 1},
            core::Precision::TwoWord};
        const auto before_qualification = provider.stamp();
        const auto proof =
            core::guard_entry(arnoldi_operation, ArnoldiCode::PrecisionUnresolved, {}, [&] {
                return policy->inverse.qualify(side, *rhs, *solution, error.values(),
                                               before_qualification, policy->linear_accuracy,
                                               expected);
            });
        if (!proof) {
            return remember(failure(proof.error()));
        }
        if (provider.stamp() != before_qualification) {
            return fail(ArnoldiCode::StaleGeneration);
        }
        if (!arnoldi_certificate(*proof, expected, policy->linear_accuracy.budget)) {
            return fail(ArnoldiCode::PrecisionUnresolved);
        }
        double native_error = 0;
        for (std::size_t i = 0; i < error.size(); ++i) {
            const auto enclosed =
                side == ArnoldiSide::Right
                    ? scale({error.values()[i], error.values()[i]}, scales.metric.eigenvalue_scale,
                            scales.metric.components[i])
                    : scale(scale({error.values()[i], error.values()[i]},
                                  scales.metric.eigenvalue_scale),
                            scales.metric.components[i]);
            if (!std::isfinite(enclosed.upper)) {
                return fail(ArnoldiCode::PrecisionUnresolved);
            }
            native_error = std::max(native_error, enclosed.upper);
        }
        if (data.inverse_count >= data.inverse_assessments.size()) {
            return fail(ArnoldiCode::ResourceExceeded);
        }
        data.inverse_assessments.values()[data.inverse_count++] =
            ArnoldiInverseAssessment{side, expected.numerical_generation, native_error, *proof};
        return {};
    }
    template <StorageScalar Native>
    void perform(ArnoldiSide side, const Native* x, Native* y) noexcept {
        std::fill_n(y, static_cast<std::size_t>(data.order), Native{});
        if (first) {
            return;
        }
        // These rounded coordinates belong only to candidate generation. The
        // exact stored physical candidate is evaluated again without this rounding
        // during final original-system qualification.
        try {
            const auto n = static_cast<std::size_t>(data.order);
            for (std::size_t i = 0; i < n; ++i) {
                const auto scale = scales.metric.components[i];
                input.values()[i] = side == ArnoldiSide::Right ? ArnoldiComplex{x[i]} * scale
                                                               : ArnoldiComplex{x[i]} / scale;
                if (!arnoldi_finite(input.values()[i])) {
                    record_failure(ArnoldiCode::PrecisionUnresolved);
                    return;
                }
            }
            if (std::holds_alternative<ArnoldiShiftInvert>(request.selection)) {
                if (!inverse(side)) {
                    return;
                }
                for (std::size_t i = 0; i < n; ++i) {
                    const auto value =
                        (high.values()[i] + low.values()[i]) * scales.metric.eigenvalue_scale;
                    const auto scaled = side == ArnoldiSide::Right
                                            ? value / scales.metric.components[i]
                                            : value * scales.metric.components[i];
                    if (!arnoldi_finite(scaled)) {
                        record_failure(ArnoldiCode::PrecisionUnresolved);
                        return;
                    }
                    if constexpr (std::same_as<Native, double>) {
                        y[i] = scaled.real();
                    } else {
                        y[i] = scaled;
                    }
                }
            } else {
                if constexpr (std::same_as<Native, double> && std::same_as<Original, double>) {
                    for (std::size_t i = 0; i < n; ++i) {
                        real_input.values()[i] = input.values()[i].real();
                    }
                    if (!original_real(side, real_input.values(), std::span<double>{y, n})) {
                        return;
                    }
                    for (std::size_t i = 0; i < n; ++i) {
                        y[i] = (side == ArnoldiSide::Right ? y[i] / scales.metric.components[i]
                                                           : y[i] * scales.metric.components[i]) /
                               scales.metric.eigenvalue_scale;
                        if (!std::isfinite(y[i])) {
                            record_failure(ArnoldiCode::PrecisionUnresolved);
                            return;
                        }
                    }
                } else {
                    if (!original_complex(side, input.values(), high.values(), low.values(),
                                          error.values())) {
                        return;
                    }
                    for (std::size_t i = 0; i < n; ++i) {
                        const auto value = high.values()[i] + low.values()[i];
                        const auto scaled =
                            (side == ArnoldiSide::Right ? value / scales.metric.components[i]
                                                        : value * scales.metric.components[i]) /
                            scales.metric.eigenvalue_scale;
                        if (!arnoldi_finite(scaled)) {
                            record_failure(ArnoldiCode::PrecisionUnresolved);
                            return;
                        }
                        if constexpr (std::same_as<Native, double>) {
                            y[i] = scaled.real();
                        } else {
                            y[i] = scaled;
                        }
                    }
                }
            }
        } catch (const std::bad_alloc&) {
            record_failure(ArnoldiCode::ResourceExceeded);
        } catch (...) {
            record_failure(ArnoldiCode::FactorizationFailed);
        }
    }
};
template <StorageScalar Original, StorageScalar Native>
struct ArnoldiOperation {
    using Scalar = Native;
    ArnoldiRunner<Original>& runner;
    ArnoldiSide side;
    [[nodiscard]] Eigen::Index rows() const noexcept {
        return runner.data.order;
    }
    [[nodiscard]] Eigen::Index cols() const noexcept {
        return runner.data.order;
    }
    void perform_op(const Native* x, Native* y) const noexcept {
        runner.template perform<Native>(side, x, y);
    }
};
} // namespace uvlas::numerics::eigen::detail
