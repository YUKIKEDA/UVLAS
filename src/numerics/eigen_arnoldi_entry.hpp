#pragma once
#include "eigen_arnoldi_qualify.hpp"

namespace uvlas::numerics::eigen::detail {
template <StorageScalar Scalar>
[[nodiscard]] ArnoldiResult<std::unique_ptr<ArnoldiData>>
arnoldi_compute(const ActionProvider<Scalar>& provider, const ArnoldiScales& input_scales,
                const ArnoldiRequest<Scalar>& input_request, const ArnoldiLimits& input_limits,
                ArnoldiWorkspaceData& scratch, core::ResourcePlan& result_plan) {
    const ArnoldiRequest<Scalar> request = input_request;
    const ArnoldiLimits limits = input_limits;
    const ArnoldiScales scales = input_scales;
    const auto fail = [&](core::Failure cause) -> ArnoldiResult<std::unique_ptr<ArnoldiData>> {
        return std::unexpected(complete_failure(cause, scales.context, request.state));
    };
    const auto initial_stamp = provider.stamp();
    const auto valid = arnoldi_validate(provider, scales, request, limits);
    if (!valid) {
        return fail(valid.error());
    }
    if (provider.stamp() != initial_stamp) {
        return fail(core::Failure{ArnoldiCode::StaleGeneration, arnoldi_operation});
    }
    const auto n = provider.order();
    const core::Index capacity = request.nev + (std::same_as<Scalar, double> ? 1 : 0);
    if (capacity > limits.maximum_returned) {
        return fail(core::Failure{ArnoldiCode::ResourceExceeded, arnoldi_operation});
    }
    if (scratch.order != n || scratch.capacity < capacity) {
        return fail(core::Failure{ArnoldiCode::InvalidInput, arnoldi_operation});
    }
    const bool left = request.sides == ArnoldiSides::Both;
    const bool subspace = std::holds_alternative<ArnoldiSubspaceRequest>(request.qualification);
    const bool shifted = std::holds_alternative<ArnoldiShiftInvert>(request.selection);
    const bool projected = std::holds_alternative<linear::Binary64Projection>(request.returned);
    ArnoldiSelectionMetadata selection{ArnoldiLargestMagnitude{}};
    if (const auto* policy = std::get_if<ArnoldiShiftInvert>(&request.selection)) {
        const auto right_stamp = policy->inverse.provider(ArnoldiSide::Right).stamp();
        std::variant<NotRequested, ActionStamp> left_stamp{NotRequested{}};
        if (left) {
            left_stamp = policy->inverse.provider(ArnoldiSide::Left).stamp();
        }
        if (right_stamp.original != scales.context || right_stamp.shift != policy->shift) {
            return fail(core::Failure{ArnoldiCode::StaleGeneration, arnoldi_operation});
        }
        if (const auto* actual = std::get_if<ActionStamp>(&left_stamp)) {
            if (actual->original != scales.context || actual->shift != std::conj(policy->shift)) {
                return fail(core::Failure{ArnoldiCode::StaleGeneration, arnoldi_operation});
            }
        }
        if (limits.actions.inverse > 0) {
            const auto first = policy->first_inverse_generation.value();
            const auto last = first + limits.actions.inverse - 1;
            const auto collision = [&](core::NumericalGeneration id) {
                return id.value() >= first && id.value() <= last;
            };
            if (collision(right_stamp.transformation)) {
                return fail(core::Failure{ArnoldiCode::StaleGeneration, arnoldi_operation});
            }
            if (const auto* actual = std::get_if<ActionStamp>(&left_stamp)) {
                if (collision(actual->transformation)) {
                    return fail(core::Failure{ArnoldiCode::StaleGeneration, arnoldi_operation});
                }
            }
        }
        selection = ArnoldiShiftMetadata{policy->shift,
                                         right_stamp,
                                         left_stamp,
                                         policy->linear_accuracy,
                                         policy->first_inverse_generation,
                                         limits.actions.inverse};
    }
    const ArnoldiBackend versions{
        {SPECTRA_MAJOR_VERSION, SPECTRA_MINOR_VERSION, SPECTRA_PATCH_VERSION},
        {EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION, EIGEN_PATCH_VERSION},
        std::same_as<Scalar, double> ? ArnoldiDomain::Real : ArnoldiDomain::Complex,
        std::same_as<Scalar, double> && !shifted ? ArnoldiDomain::Real : ArnoldiDomain::Complex};
    auto data = std::make_unique<ArnoldiData>(
        scales.context, initial_stamp, request.state, request.returned, request.association,
        selection, versions, scales.eigenvalue_units, scales.metric.eigenvalue_scale, left,
        subspace,
        ArnoldiSettings{request.nev, request.ncv, request.backend_tolerance,
                        request.backend_iterations, request.sides, request.qualification,
                        request.timing},
        limits);
    data->order = n;
    data->projected = projected;
    data->diagnostics.requested_nev = request.nev;
    data->diagnostics.requested_ncv = request.ncv;
    ArnoldiTimings measured{};
    const bool timed = request.timing == ArnoldiTimingMode::Enabled;
    // Checked matrix/byte products precede both external backend construction and
    // every numerical callback. The forecast is deliberately not a malloc cap.
    const std::array<PayloadExtent, 5> external{
        {{n, request.ncv, n, sizeof(ArnoldiComplex) * 8},
         {request.ncv, request.ncv, request.ncv, sizeof(ArnoldiComplex) * 20},
         {n, capacity, n, sizeof(ArnoldiComplex) * 8},
         {capacity, capacity, capacity, sizeof(ArnoldiComplex) * 16},
         {n, 1, n, sizeof(ArnoldiComplex) * 16}}};
    const auto estimated = checked_payload_bytes(external);
    if (!estimated) {
        return fail(estimated.error());
    }
    if (left && *estimated > std::numeric_limits<std::uint64_t>::max() / 2) {
        return fail(core::Failure{ArnoldiCode::ResourceExceeded, arnoldi_operation});
    }
    const auto admission = linear::admit_backend_memory(
        linear::detail::BackendPlanFactory::make(*estimated * (left ? 2U : 1U)),
        limits.backend_memory);
    if (!admission) {
        auto cause = admission.error();
        // The shared admission model is internal here, not a producer callback.
        cause.operation = arnoldi_operation;
        return fail(cause);
    }
    data->admission = *admission;
    {
        const ArnoldiPhase snapshot_phase(timed, measured.preparation);
        if (auto r = arnoldi_allocate(data->components, n, result_plan); !r) {
            return fail(r.error());
        }
        if (auto r = arnoldi_allocate(data->units, n, result_plan); !r) {
            return fail(r.error());
        }
        if (auto r = arnoldi_allocate(data->initial, n, result_plan); !r) {
            return fail(r.error());
        }
        if (!std::in_range<core::Index>(request.required_checks.size()) ||
            request.required_checks.size() >
                static_cast<std::size_t>(std::numeric_limits<core::Index>::max() / 2)) {
            return fail(core::Failure{ArnoldiCode::ResourceExceeded, arnoldi_operation});
        }
        if (auto r = arnoldi_allocate(data->required_checks,
                                      static_cast<core::Index>(request.required_checks.size()),
                                      result_plan);
            !r) {
            return fail(r.error());
        }
        for (std::size_t i = 0; i < request.required_checks.size(); ++i) {
            data->required_checks.values()[i] = request.required_checks[i];
        }
        std::ranges::copy(scales.metric.components, data->components.values().begin());
        std::ranges::copy(scales.component_units, data->units.values().begin());
        arnoldi_seed(*data, request.initial);
    }
    auto actions = ActionSession<Scalar>::make(provider, limits.actions, request.state);
    if (!actions) {
        return std::unexpected(actions.error());
    }
    if (provider.stamp() != initial_stamp) {
        return fail(core::Failure{ArnoldiCode::StaleGeneration, arnoldi_operation});
    }
    ArnoldiRunner<Scalar> runner(provider, scales, request, limits, *data, std::move(*actions),
                                 scratch);
    runner.scales = ArnoldiScales{{data->components.values(), data->eigenvalue_scale},
                                  data->units.values(),
                                  data->eigenvalue_units,
                                  data->context};
    const auto setup_failure =
        [&](core::Failure cause) -> ArnoldiResult<std::unique_ptr<ArnoldiData>> {
        return std::unexpected(runner.failure(cause));
    };
    try {
        auto& left_candidates = scratch.left_candidates;
        auto& assigned = scratch.assigned;
        left_candidates.count = 0;
        {
            const ArnoldiPhase preparation(timed, measured.preparation);
            if (auto r = arnoldi_packet_allocate(data->retained, n, capacity, left, subspace,
                                                 result_plan);
                !r) {
                return setup_failure(r.error());
            }
            if (projected) {
                if (auto r = arnoldi_packet_allocate(data->projection, n, capacity, left, subspace,
                                                     result_plan);
                    !r) {
                    return setup_failure(r.error());
                }
            }
            if (auto r =
                    arnoldi_allocate(data->sensitivities, subspace ? 0 : capacity, result_plan);
                !r) {
                return setup_failure(r.error());
            }
            if (auto r = arnoldi_allocate(data->consumer_proofs,
                                          static_cast<core::Index>(request.required_checks.size()) *
                                              (projected ? 2 : 1),
                                          result_plan);
                !r) {
                return setup_failure(r.error());
            }
            if (shifted) {
                if (!std::in_range<core::Index>(limits.actions.inverse)) {
                    return setup_failure(
                        core::Failure{ArnoldiCode::ResourceExceeded, arnoldi_operation});
                }
                if (auto r = arnoldi_allocate(data->inverse_assessments,
                                              static_cast<core::Index>(limits.actions.inverse),
                                              result_plan);
                    !r) {
                    return setup_failure(r.error());
                }
            }
            if (const auto* policy = std::get_if<ArnoldiShiftInvert>(&request.selection)) {
                auto right_inverse = ActionSession<ArnoldiComplex>::make(
                    policy->inverse.provider(ArnoldiSide::Right), limits.actions, request.state);
                if (!right_inverse) {
                    return std::unexpected(right_inverse.error());
                }
                runner.inverse_right.emplace(std::move(*right_inverse));
                if (left) {
                    auto left_inverse = ActionSession<ArnoldiComplex>::make(
                        policy->inverse.provider(ArnoldiSide::Left), limits.actions, request.state);
                    if (!left_inverse) {
                        return std::unexpected(left_inverse.error());
                    }
                    runner.inverse_left.emplace(std::move(*left_inverse));
                }
            }
        }
        const auto extract = [&](ArnoldiSide side, ArnoldiPacket& packet) -> ArnoldiResult<void> {
            if constexpr (std::same_as<Scalar, double>) {
                if (!shifted) {
                    return arnoldi_extract_side<Scalar, double>(
                        runner, side, packet, measured.krylov, timed, measured.subspace);
                }
            }
            return arnoldi_extract_side<Scalar, ArnoldiComplex>(
                runner, side, packet, measured.krylov, timed, measured.subspace);
        };
        if (auto r = extract(ArnoldiSide::Right, data->retained); !r) {
            return std::unexpected(r.error());
        }
        arnoldi_sort(data->retained, n, subspace);
        if (left) {
            if (auto r = extract(ArnoldiSide::Left, left_candidates); !r) {
                return std::unexpected(r.error());
            }
            if (auto r =
                    arnoldi_match_left(runner, left_candidates, data->retained, assigned.values());
                !r) {
                return std::unexpected(r.error());
            }
        }
        {
            const ArnoldiPhase qualification(timed, measured.qualification);
            if (auto r =
                    arnoldi_qualify(runner, data->retained, retained_generation(request.returned));
                !r) {
                return std::unexpected(r.error());
            }
            if (auto r = arnoldi_consumers(runner, data->retained,
                                           retained_generation(request.returned), 0);
                !r) {
                return std::unexpected(r.error());
            }
            if (const auto* projection =
                    std::get_if<linear::Binary64Projection>(&request.returned)) {
                data->projection.count = data->retained.count;
                const auto copy_actual = [](const auto& source, auto& target, std::size_t count) {
                    std::ranges::copy(source.values().first(count), target.values().begin());
                };
                const auto count = static_cast<std::size_t>(data->retained.count);
                copy_actual(data->retained.eigenvalues, data->projection.eigenvalues, count);
                copy_actual(data->retained.right, data->projection.right,
                            count * static_cast<std::size_t>(n));
                if (left) {
                    copy_actual(data->retained.left, data->projection.left,
                                count * static_cast<std::size_t>(n));
                }
                if (subspace) {
                    copy_actual(data->retained.right_t, data->projection.right_t, count * count);
                    copy_actual(data->retained.right_mapping, data->projection.right_mapping,
                                count);
                    if (left) {
                        copy_actual(data->retained.left_t, data->projection.left_t, count * count);
                        copy_actual(data->retained.left_mapping, data->projection.left_mapping,
                                    count);
                    }
                }
                if (auto r =
                        arnoldi_qualify(runner, data->projection, projection->target_generation);
                    !r) {
                    return std::unexpected(r.error());
                }
                if (auto r =
                        arnoldi_consumers(runner, data->projection, projection->target_generation,
                                          request.required_checks.size());
                    !r) {
                    return std::unexpected(r.error());
                }
            }
        }
        // The opaque sources must still denote the frozen system after every consumer.
        if (provider.stamp() != initial_stamp) {
            return setup_failure(core::Failure{ArnoldiCode::StaleGeneration, arnoldi_operation});
        }
        const auto frozen = core::guard_entry(arnoldi_operation, ArnoldiCode::FactorizationFailed,
                                              {}, [&] { return provider.validate_frozen_state(); });
        if (!frozen) {
            return setup_failure(frozen.error());
        }
        if (provider.stamp() != initial_stamp) {
            return setup_failure(core::Failure{ArnoldiCode::StaleGeneration, arnoldi_operation});
        }
        if (const auto* policy = std::get_if<ArnoldiShiftInvert>(&request.selection)) {
            const auto* metadata = std::get_if<ArnoldiShiftMetadata>(&data->selection);
            if (metadata == nullptr) {
                return setup_failure(core::Failure{ArnoldiCode::InvalidInput, arnoldi_operation});
            }
            for (const auto side : {ArnoldiSide::Right, ArnoldiSide::Left}) {
                if (side == ArnoldiSide::Left && !left) {
                    continue;
                }
                const auto& inverse = policy->inverse.provider(side);
                const auto* left_stamp = std::get_if<ActionStamp>(&metadata->left_inverse);
                const auto expected =
                    side == ArnoldiSide::Right
                        ? metadata->right_inverse
                        : (left_stamp != nullptr ? *left_stamp : metadata->right_inverse);
                if (inverse.stamp() != expected) {
                    return setup_failure(
                        core::Failure{ArnoldiCode::StaleGeneration, arnoldi_operation});
                }
                const auto checked =
                    core::guard_entry(arnoldi_operation, ArnoldiCode::PrecisionUnresolved, {},
                                      [&] { return inverse.validate_frozen_state(); });
                if (!checked) {
                    return setup_failure(checked.error());
                }
                if (inverse.stamp() != expected) {
                    return setup_failure(
                        core::Failure{ArnoldiCode::StaleGeneration, arnoldi_operation});
                }
            }
        }
        arnoldi_sensitivity(runner, data->projected ? data->projection : data->retained);
        data->diagnostics.returned_count = data->retained.count;
        data->diagnostics.attempts = runner.counters();
        if (timed) {
            data->measurements = measured;
        }
        return data;
    } catch (const std::bad_alloc&) {
        if (runner.first) {
            return std::unexpected(*runner.first);
        }
        return setup_failure(core::Failure{ArnoldiCode::ResourceExceeded, arnoldi_operation});
    } catch (...) {
        if (runner.first) {
            return std::unexpected(*runner.first);
        }
        return setup_failure(core::Failure{ArnoldiCode::FactorizationFailed, arnoldi_operation});
    }
}
} // namespace uvlas::numerics::eigen::detail
