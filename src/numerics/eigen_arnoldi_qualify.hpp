#pragma once
#include "eigen_arnoldi_build.hpp"

namespace uvlas::numerics::eigen::detail {
template <StorageScalar Scalar>
[[nodiscard]] ArnoldiResult<void> arnoldi_qualify_pairs(ArnoldiRunner<Scalar>& runner,
                                                        ArnoldiPacket& packet,
                                                        core::NumericalGeneration generation) {
    const auto* requirements = std::get_if<ArnoldiPairRequest>(&runner.request.qualification);
    if (requirements == nullptr) {
        return runner.fail(ArnoldiCode::InvalidInput);
    }
    const auto n = static_cast<std::size_t>(runner.data.order);
    const MetricContext identity{runner.scales.context, runner.request.state, generation};
    for (const auto side : {ArnoldiSide::Right, ArnoldiSide::Left}) {
        if (side == ArnoldiSide::Left && !runner.data.left_requested) {
            continue;
        }
        const auto& storage = side == ArnoldiSide::Right ? packet.right : packet.left;
        auto& proofs = side == ArnoldiSide::Right ? packet.right_pairs : packet.left_pairs;
        for (core::Index j = 0; j < packet.count; ++j) {
            const auto column = storage.values().subspan(static_cast<std::size_t>(j) * n, n);
            if (auto r = runner.original_complex(side, column, runner.high.values(),
                                                 runner.low.values(), runner.error.values());
                !r) {
                return r;
            }
            const auto vector =
                WordView<ArnoldiComplex>::make(column, {}, core::Precision::Binary64);
            const auto value = WordView<ArnoldiComplex>::make(
                packet.eigenvalues.values().subspan(static_cast<std::size_t>(j), 1), {},
                core::Precision::Binary64);
            const auto action = WordView<ArnoldiComplex>::make(
                runner.high.values(), runner.low.values(), core::Precision::TwoWord);
            if (!vector) {
                return runner.remember(runner.failure(vector.error()));
            }
            if (!value) {
                return runner.remember(runner.failure(value.error()));
            }
            if (!action) {
                return runner.remember(runner.failure(action.error()));
            }
            const auto proof = assess_pair_metrics(
                *vector, *value, ActionValues<ArnoldiComplex>{*action, runner.error.values()},
                runner.scales.metric, side == ArnoldiSide::Right ? PairSide::Right : PairSide::Left,
                requirements->accuracy, identity);
            if (!proof) {
                return runner.remember(proof.error());
            }
            proofs.values()[static_cast<std::size_t>(j)] = *proof;
        }
    }
    return {};
}

template <StorageScalar Scalar>
[[nodiscard]] ArnoldiResult<void>
arnoldi_consumers(ArnoldiRunner<Scalar>& runner, ArnoldiPacket& packet,
                  core::NumericalGeneration generation, std::size_t offset) {
    const auto retained = arnoldi_view(runner.data, runner.data.retained,
                                       retained_generation(runner.request.returned));
    const auto candidate = arnoldi_view(runner.data, packet, generation);
    if (runner.data.required_checks.size() != 0) {
        auto precision = require_precision(core::Precision::Binary64,
                                           runner.request.evaluator.supported_precision());
        if (!precision) {
            auto cause = precision.error();
            cause.operation = arnoldi_operation;
            return runner.remember(runner.failure(cause));
        }
    }
    const auto fail_consumer = [&](core::Failure cause,
                                   const linear::ConsumerCheck& check) -> ArnoldiResult<void> {
        if (!cause.details.problem) {
            cause.details.problem = check.observation;
        }
        if (!cause.details.numerical_generation) {
            cause.details.numerical_generation = generation;
        }
        return runner.remember(runner.failure(cause));
    };
    std::size_t i = 0;
    for (const auto& owned : runner.data.required_checks.values()) {
        if (!owned) {
            return runner.fail(ArnoldiCode::InvalidInput);
        }
        const auto& check = *owned;
        const core::CertificateContext expected{check.observation, check.scope, generation,
                                                core::Precision::Binary64};
        const auto proof =
            core::guard_entry(arnoldi_operation, ArnoldiCode::PrecisionUnresolved, {}, [&] {
                return runner.request.evaluator.assess(retained, candidate, check, expected);
            });
        if (!proof) {
            return fail_consumer(proof.error(), check);
        }
        if (!arnoldi_certificate(*proof, expected, check.budget)) {
            return fail_consumer(core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation},
                                 check);
        }
        runner.data.consumer_proofs.values()[offset + i] = *proof;
        ++i;
    }
    runner.data.consumer_count = i;
    return {};
}

template <StorageScalar Scalar>
[[nodiscard]] ArnoldiResult<void> arnoldi_qualify_basis(ArnoldiRunner<Scalar>& runner,
                                                        ArnoldiPacket& packet,
                                                        core::NumericalGeneration generation) {
    const auto* requested = std::get_if<ArnoldiSubspaceRequest>(&runner.request.qualification);
    if (requested == nullptr) {
        return runner.fail(ArnoldiCode::InvalidInput);
    }
    const auto n = static_cast<std::size_t>(runner.data.order),
               k = static_cast<std::size_t>(packet.count);
    for (const auto side : {ArnoldiSide::Right, ArnoldiSide::Left}) {
        if (side == ArnoldiSide::Left && !runner.data.left_requested) {
            continue;
        }
        const auto& basis = side == ArnoldiSide::Right ? packet.right : packet.left;
        const auto& triangular = side == ArnoldiSide::Right ? packet.right_t : packet.left_t;
        for (std::size_t j = 0; j < k; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const auto physical = arnoldi_physical_input(
                    basis.values()[j * n + i], runner.scales.metric.components[i], side,
                    runner.input.values()[i], runner.input_low.values()[i]);
                if (!physical) {
                    return runner.remember(runner.failure(physical.error()));
                }
            }
            if (auto r = runner.original_complex(side, runner.input.values(), runner.high.values(),
                                                 runner.low.values(), runner.error.values(),
                                                 runner.input_low.values());
                !r) {
                return r;
            }
            for (std::size_t i = 0; i < n; ++i) {
                const auto converted = arnoldi_scaled_action(
                    runner.high.values()[i], runner.low.values()[i], runner.error.values()[i],
                    runner.scales.metric.components[i], runner.scales.metric.eigenvalue_scale, side,
                    runner.cluster_action.values()[j * n + i],
                    runner.cluster_error.values()[j * n + i]);
                if (!converted) {
                    return runner.remember(runner.failure(converted.error()));
                }
            }
        }
        std::fill(runner.triangular_error.values().begin(), runner.triangular_error.values().end(),
                  0.0);
        // Physical lambda is itself rounded. Relate that actual public value to
        // the corresponding original restriction diagonal, without discarding it.
        const auto& mapping =
            side == ArnoldiSide::Right ? packet.right_mapping : packet.left_mapping;
        for (std::size_t j = 0; j < k; ++j) {
            const auto diagonal = static_cast<std::size_t>(mapping.values()[j]);
            if (diagonal >= k) {
                return runner.fail(ArnoldiCode::InvalidInput);
            }
            auto physical = packet.eigenvalues.values()[j];
            if (side == ArnoldiSide::Left) {
                physical = std::conj(physical);
            }
            Sum<ArnoldiComplex> difference;
            const auto added = difference.add(physical);
            const auto product =
                difference.product(triangular.values()[diagonal * k + diagonal],
                                   ArnoldiComplex{runner.scales.metric.eigenvalue_scale}, -1);
            if (!added) {
                return runner.remember(runner.failure(added.error()));
            }
            if (!product) {
                return runner.remember(runner.failure(product.error()));
            }
            const auto evaluated = difference.finish();
            if (!evaluated) {
                return runner.remember(runner.failure(evaluated.error()));
            }
            const auto center = magnitude(evaluated->value);
            if (!center) {
                return runner.remember(runner.failure(center.error()));
            }
            const auto upper = add_magnitudes(*center, {evaluated->error, evaluated->error}).upper;
            runner.triangular_error.values()[diagonal * k + diagonal] =
                divide_bound(upper, runner.scales.metric.eigenvalue_scale, true);
        }
        const auto q = WordView<ArnoldiComplex>::make(basis.values().first(n * k), {},
                                                      core::Precision::Binary64);
        const auto t = WordView<ArnoldiComplex>::make(triangular.values().first(k * k), {},
                                                      core::Precision::Binary64);
        const auto aq = WordView<ArnoldiComplex>::make(runner.cluster_action.values().first(n * k),
                                                       {}, core::Precision::Binary64);
        if (!q) {
            return runner.remember(runner.failure(q.error()));
        }
        if (!t) {
            return runner.remember(runner.failure(t.error()));
        }
        if (!aq) {
            return runner.remember(runner.failure(aq.error()));
        }
        const auto proof = assess_schur_metrics(
            runner.data.order, packet.count, *q,
            ActionValues<ArnoldiComplex>{*t, runner.triangular_error.values().first(k * k)},
            ActionValues<ArnoldiComplex>{*aq, runner.cluster_error.values().first(n * k)},
            requested->accuracy,
            MetricContext{runner.scales.context, runner.request.state, generation});
        if (!proof) {
            return runner.remember(proof.error());
        }
        (side == ArnoldiSide::Right ? packet.right_basis : packet.left_basis) = *proof;
    }
    return {};
}

template <StorageScalar Scalar>
[[nodiscard]] ArnoldiResult<void> arnoldi_qualify(ArnoldiRunner<Scalar>& runner,
                                                  ArnoldiPacket& packet,
                                                  core::NumericalGeneration generation) {
    return runner.data.subspace_requested ? arnoldi_qualify_basis(runner, packet, generation)
                                          : arnoldi_qualify_pairs(runner, packet, generation);
}

template <StorageScalar Scalar>
inline void arnoldi_sensitivity(ArnoldiRunner<Scalar>& runner,
                                const ArnoldiPacket& packet) noexcept {
    if (runner.data.subspace_requested) {
        return;
    }
    for (core::Index j = 0; j < packet.count; ++j) {
        auto& output = runner.data.sensitivities.values()[static_cast<std::size_t>(j)];
        if (!runner.data.left_requested) {
            output = NotRequested{};
            continue;
        }
        output = ArnoldiConditionUnresolved{};
        Sum<ArnoldiComplex> dot;
        bool valid = true;
        double right_norm = 0, left_norm = 0;
        for (core::Index i = 0; i < runner.data.order; ++i) {
            const auto offset = static_cast<std::size_t>(j * runner.data.order + i);
            const auto right = packet.right.values()[offset] /
                               runner.scales.metric.components[static_cast<std::size_t>(i)];
            const auto left = packet.left.values()[offset] *
                              runner.scales.metric.components[static_cast<std::size_t>(i)];
            if (!arnoldi_finite(right) || !arnoldi_finite(left)) {
                valid = false;
                break;
            }
            right_norm = std::hypot(right_norm, std::abs(right));
            left_norm = std::hypot(left_norm, std::abs(left));
            if (!dot.product(std::conj(packet.left.values()[offset]),
                             packet.right.values()[offset])) {
                valid = false;
                break;
            }
        }
        if (!valid) {
            continue;
        }
        const auto evaluated = dot.finish();
        if (!evaluated) {
            continue;
        }
        const auto center = magnitude(evaluated->value);
        if (!center || !(lower_subtract(center->lower, evaluated->error) > 0)) {
            continue;
        }
        const auto estimate = (right_norm * left_norm) / std::abs(evaluated->value);
        if (std::isfinite(estimate) && estimate >= 0) {
            output = ArnoldiConditionEstimate{estimate};
        }
    }
}

// Stable association is deliberately a partial matching: ambiguity rejects instead
// of silently changing the requested algebraic set or reusing one left column twice.
template <StorageScalar Scalar>
[[nodiscard]] ArnoldiResult<void>
arnoldi_match_left(ArnoldiRunner<Scalar>& runner, const ArnoldiPacket& left, ArnoldiPacket& target,
                   std::span<core::Index> assigned) noexcept {
    if (left.count != target.count) {
        return runner.fail(ArnoldiCode::EigenCoverageUnresolved);
    }
    std::fill(assigned.begin(), assigned.end(), core::Index{-1});
    for (core::Index j = 0; j < target.count; ++j) {
        core::Index match = -1;
        for (core::Index k = 0; k < left.count; ++k) {
            bool used = false;
            for (core::Index p = 0; p < j; ++p) {
                used = used || assigned[static_cast<std::size_t>(p)] == k;
            }
            if (!used &&
                arnoldi_near(target.eigenvalues.values()[static_cast<std::size_t>(j)],
                             left.eigenvalues.values()[static_cast<std::size_t>(k)],
                             runner.request.association, runner.scales.metric.eigenvalue_scale)) {
                if (match >= 0) {
                    if (runner.data.subspace_requested) {
                        continue;
                    }
                    return runner.fail(ArnoldiCode::PrecisionUnresolved);
                }
                match = k;
            }
        }
        if (match < 0) {
            return runner.fail(ArnoldiCode::PrecisionUnresolved);
        }
        assigned[static_cast<std::size_t>(j)] = match;
        if (runner.data.subspace_requested) {
            target.left_mapping.values()[static_cast<std::size_t>(j)] =
                left.left_mapping.values()[static_cast<std::size_t>(match)];
        }
        for (core::Index i = 0; i < runner.data.order; ++i) {
            target.left.values()[static_cast<std::size_t>(j * runner.data.order + i)] =
                left.left.values()[static_cast<std::size_t>(
                    (runner.data.subspace_requested ? j : match) * runner.data.order + i)];
        }
    }
    if (runner.data.subspace_requested) {
        const auto count = static_cast<std::size_t>(target.count);
        std::ranges::copy(left.left_t.values().first(count * count),
                          target.left_t.values().begin());
    }
    return {};
}
} // namespace uvlas::numerics::eigen::detail
