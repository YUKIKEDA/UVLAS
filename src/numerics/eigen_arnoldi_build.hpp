#pragma once
#include "backend_memory_detail.hpp"
#include "eigen_arnoldi_extract.hpp"

#include <boost/hash2/sha2.hpp>

namespace uvlas::numerics::eigen::detail {
[[nodiscard]] inline core::Result<void> arnoldi_packet_allocate(ArnoldiPacket& packet,
                                                                core::Index n, core::Index capacity,
                                                                bool left, bool subspace,
                                                                core::ResourcePlan& plan) noexcept {
    const std::array dimensions{n, capacity};
    const auto entries = core::checked_extent_product(dimensions, arnoldi_operation);
    if (!entries) {
        return std::unexpected(entries.error());
    }
    const std::array square{capacity, capacity};
    const auto small = core::checked_extent_product(square, arnoldi_operation);
    if (!small) {
        return std::unexpected(small.error());
    }
    if (auto r = arnoldi_allocate(packet.eigenvalues, capacity, plan); !r) {
        return r;
    }
    if (auto r = arnoldi_allocate(packet.right, *entries, plan); !r) {
        return r;
    }
    if (left) {
        if (auto r = arnoldi_allocate(packet.left, *entries, plan); !r) {
            return r;
        }
    }
    if (subspace) {
        if (auto r = arnoldi_allocate(packet.right_mapping, capacity, plan); !r) {
            return r;
        }
        if (left) {
            if (auto r = arnoldi_allocate(packet.left_mapping, capacity, plan); !r) {
                return r;
            }
        }
        if (auto r = arnoldi_allocate(packet.right_t, *small, plan); !r) {
            return r;
        }
        if (left) {
            if (auto r = arnoldi_allocate(packet.left_t, *small, plan); !r) {
                return r;
            }
        }
    } else {
        if (auto r = arnoldi_allocate(packet.right_pairs, capacity, plan); !r) {
            return r;
        }
        if (left) {
            if (auto r = arnoldi_allocate(packet.left_pairs, capacity, plan); !r) {
                return r;
            }
        }
    }
    return {};
}

template <StorageScalar Scalar>
[[nodiscard]] core::Result<void>
arnoldi_validate(const ActionProvider<Scalar>& provider, const ArnoldiScales& scales,
                 const ArnoldiRequest<Scalar>& request, const ArnoldiLimits& limits) noexcept {
    const auto invalid = [](ArnoldiCode code) {
        return std::unexpected(core::Failure{code, arnoldi_operation});
    };
    const auto n = provider.order();
    if (n < 3 || request.nev < 1 || request.nev > n - 2 || request.ncv < request.nev + 2 ||
        request.ncv > n || std::cmp_not_equal(request.initial.size(), n) ||
        std::cmp_not_equal(scales.metric.components.size(), n) ||
        std::cmp_not_equal(scales.component_units.size(), n) || !std::in_range<Eigen::Index>(n) ||
        !std::in_range<Eigen::Index>(request.backend_iterations) ||
        request.backend_iterations == 0) {
        return invalid(ArnoldiCode::InvalidInput);
    }
    if (request.initial.precision() != core::Precision::Binary64) {
        return invalid(ArnoldiCode::PrecisionUnsupported);
    }
    const auto original_stamp = provider.stamp();
    if (!std::isfinite(scales.context.time) || !std::isfinite(original_stamp.original.time) ||
        !arnoldi_finite(original_stamp.shift)) {
        return invalid(ArnoldiCode::Nonfinite);
    }
    if (scales.context != original_stamp.original) {
        return invalid(ArnoldiCode::StaleGeneration);
    }
    if ((request.sides != ArnoldiSides::RightOnly && request.sides != ArnoldiSides::Both) ||
        (request.timing != ArnoldiTimingMode::Disabled &&
         request.timing != ArnoldiTimingMode::Enabled)) {
        return invalid(ArnoldiCode::InvalidInput);
    }
    if (!std::isfinite(request.backend_tolerance) ||
        !std::isfinite(scales.metric.eigenvalue_scale) ||
        !std::isfinite(request.association.absolute) ||
        !std::isfinite(request.association.relative)) {
        return invalid(ArnoldiCode::Nonfinite);
    }
    if (!(request.backend_tolerance > 0) || !(scales.metric.eigenvalue_scale > 0) ||
        request.association.absolute < 0 || request.association.relative < 0) {
        return invalid(ArnoldiCode::InvalidInput);
    }
    for (const auto scale : scales.metric.components) {
        if (!std::isfinite(scale)) {
            return invalid(ArnoldiCode::Nonfinite);
        }
        if (!(scale > 0)) {
            return invalid(ArnoldiCode::InvalidInput);
        }
    }
    const auto seed = WordView<Scalar>::make(request.initial.high(), {}, core::Precision::Binary64);
    if (!seed) {
        return std::unexpected(seed.error());
    }
    double seed_norm = 0;
    for (const auto value : seed->high()) {
        seed_norm = std::hypot(seed_norm, std::abs(value));
    }
    if (!std::isfinite(seed_norm)) {
        return invalid(ArnoldiCode::PrecisionUnresolved);
    }
    if (!(seed_norm > 0)) {
        return invalid(ArnoldiCode::InvalidInput);
    }
    const auto requirement = [&](MetricRequirement wanted) -> core::Result<void> {
        if (!std::isfinite(wanted.budget)) {
            return invalid(ArnoldiCode::Nonfinite);
        }
        if (wanted.budget < 0) {
            return invalid(ArnoldiCode::InvalidInput);
        }
        return {};
    };
    const auto check_accuracy = [&](const auto& choice) -> core::Result<void> {
        const auto& accuracy = choice.accuracy;
        if (auto r = requirement(accuracy.residual); !r) {
            return r;
        }
        // EIG-03/N4-07 mandate original norm residual at most 1e-10.
        if (accuracy.residual.budget > 1e-10) {
            return invalid(ArnoldiCode::InvalidInput);
        }
        if constexpr (std::same_as<std::decay_t<decltype(choice)>, ArnoldiPairRequest>) {
            if (auto r = requirement(accuracy.normalization); !r) {
                return r;
            }
        } else {
            if (auto r = requirement(accuracy.orthogonality); !r) {
                return r;
            }
        }
        if (const auto* beta = std::get_if<MetricRequirement>(&accuracy.backward)) {
            return requirement(*beta);
        }
        if (!std::holds_alternative<NotRequested>(accuracy.backward)) {
            return invalid(ArnoldiCode::InvalidInput);
        }
        return {};
    };
    core::Result<void> requests;
    if (const auto* pair = std::get_if<ArnoldiPairRequest>(&request.qualification)) {
        requests = check_accuracy(*pair);
    } else if (const auto* basis = std::get_if<ArnoldiSubspaceRequest>(&request.qualification)) {
        requests = check_accuracy(*basis);
    } else {
        return invalid(ArnoldiCode::InvalidInput);
    }
    if (!requests) {
        return requests;
    }
    for (const auto& check : request.required_checks) {
        if (!std::isfinite(check.budget)) {
            return invalid(ArnoldiCode::Nonfinite);
        }
        if (check.budget < 0 || (check.kind != linear::CheckKind::Observation &&
                                 check.kind != linear::CheckKind::Constraint &&
                                 check.kind != linear::CheckKind::Reaction)) {
            return invalid(ArnoldiCode::InvalidInput);
        }
    }
    if (const auto* projection = std::get_if<linear::Binary64Projection>(&request.returned)) {
        if (projection->retained_generation == projection->target_generation) {
            return invalid(ArnoldiCode::StaleGeneration);
        }
    } else if (!std::holds_alternative<linear::RetainedPrecision>(request.returned)) {
        return invalid(ArnoldiCode::PrecisionUnsupported);
    }
    if (const auto* shift = std::get_if<ArnoldiShiftInvert>(&request.selection)) {
        if (!arnoldi_finite(shift->shift)) {
            return invalid(ArnoldiCode::Nonfinite);
        }
        if (auto r = requirement(shift->linear_accuracy); !r) {
            return r;
        }
        if (limits.actions.inverse > 0) {
            const auto first = shift->first_inverse_generation.value();
            if (limits.actions.inverse - 1 > std::numeric_limits<std::uint64_t>::max() - first) {
                return invalid(ArnoldiCode::ResourceExceeded);
            }
            const auto last = first + limits.actions.inverse - 1;
            const auto collision = [&](core::NumericalGeneration id) {
                return id.value() >= first && id.value() <= last;
            };
            if (collision(scales.context.cache.matrix) ||
                collision(retained_generation(request.returned)) ||
                collision(provider.stamp().transformation)) {
                return invalid(ArnoldiCode::StaleGeneration);
            }
            if (const auto* projection =
                    std::get_if<linear::Binary64Projection>(&request.returned)) {
                if (collision(projection->target_generation)) {
                    return invalid(ArnoldiCode::StaleGeneration);
                }
            }
        }
    } else if (!std::holds_alternative<ArnoldiLargestMagnitude>(request.selection)) {
        return invalid(ArnoldiCode::InvalidInput);
    }
    if (limits.maximum_order < n || limits.maximum_ncv < request.ncv ||
        request.backend_iterations > limits.maximum_backend_iterations ||
        limits.maximum_returned < request.nev || limits.maximum_order < 0 ||
        limits.maximum_ncv < 0 || limits.maximum_returned < 0) {
        return invalid(ArnoldiCode::ResourceExceeded);
    }
    return {};
}

inline void arnoldi_hash_u64(boost::hash2::sha2_256& hash, std::uint64_t value) noexcept {
    std::array<unsigned char, 8> bytes{};
    for (auto& byte : bytes) {
        byte = static_cast<unsigned char>(value & 255U);
        value >>= 8U;
    }
    hash.update(bytes.data(), bytes.size());
}
template <StorageScalar Scalar>
inline void arnoldi_seed(ArnoldiData& data, WordView<Scalar> seed) noexcept {
    boost::hash2::sha2_256 hash;
    constexpr std::array<unsigned char, 8> domain{'U', 'V', 'A', 'R', 'N', '0', '0', '1'};
    hash.update(domain.data(), domain.size());
    arnoldi_hash_u64(hash, std::same_as<Scalar, double> ? 0 : 1);
    arnoldi_hash_u64(hash, static_cast<std::uint64_t>(seed.precision()));
    arnoldi_hash_u64(hash, seed.size());
    double norm = 0;
    for (std::size_t i = 0; i < seed.size(); ++i) {
        const ArnoldiComplex value{seed.high()[i]};
        data.initial.values()[i] = value;
        norm = std::hypot(norm, std::abs(value));
        arnoldi_hash_u64(hash, std::bit_cast<std::uint64_t>(value.real()));
        if constexpr (!std::same_as<Scalar, double>) {
            arnoldi_hash_u64(hash, std::bit_cast<std::uint64_t>(value.imag()));
        }
    }
    data.diagnostics.initial_norm = norm;
    const auto digest = hash.result();
    for (std::size_t i = 0; i < data.diagnostics.initial_sha256.size(); ++i) {
        data.diagnostics.initial_sha256[i] = static_cast<std::byte>(digest[i]);
    }
}
} // namespace uvlas::numerics::eigen::detail
