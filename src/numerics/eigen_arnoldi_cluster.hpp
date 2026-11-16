#pragma once
#include "eigen_arnoldi_extract.hpp"
#include "eigen_metric_bounds.hpp"

namespace uvlas::numerics::eigen::detail {
// Construct an exact two-word physical input for a dimensionless basis column.
// Division can require infinitely many bits. In that case no opaque action
// forward bound is invented for the discarded input tail.
[[nodiscard]] inline core::Result<void> arnoldi_physical_input(ArnoldiComplex q, double d,
                                                               ArnoldiSide side, ArnoldiComplex& hi,
                                                               ArnoldiComplex& lo) noexcept {
    for (unsigned lane = 0; lane < 2; ++lane) {
        const auto value = lane == 0 ? q.real() : q.imag();
        double leading = 0, trailing = 0;
        if (side == ArnoldiSide::Right) {
            const auto product = two_product(value, d);
            if (!product) {
                return std::unexpected(
                    core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
            }
            leading = product->hi;
            trailing = product->lo;
        } else {
            leading = value / d;
            const auto remainder = std::fma(-leading, d, value);
            trailing = remainder / d;
            Sum<double> identity;
            const auto added = identity.add(value);
            const auto head = identity.product(leading, d, -1);
            const auto tail = identity.product(trailing, d, -1);
            const auto equality = identity.finish();
            if (!added || !head || !tail || !equality || equality->value != 0 ||
                equality->error != 0 || !std::isfinite(leading) || !std::isfinite(trailing)) {
                return std::unexpected(
                    core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
            }
        }
        if (lane == 0) {
            hi.real(leading);
            lo.real(trailing);
        } else {
            hi.imag(leading);
            lo.imag(trailing);
        }
    }
    return {};
}

// Transform an enclosed physical action to the original dimensionless action.
// Roundoff of each explicit scaling is included, never absorbed into the budget.
[[nodiscard]] inline core::Result<void> arnoldi_scaled_action(ArnoldiComplex hi, ArnoldiComplex lo,
                                                              double error, double d,
                                                              double lambda_scale, ArnoldiSide side,
                                                              ArnoldiComplex& value,
                                                              double& upper) noexcept {
    const auto unscaled = hi + lo;
    const auto coordinate = side == ArnoldiSide::Right ? unscaled / d : unscaled * d;
    value = coordinate / lambda_scale;
    if (!arnoldi_finite(unscaled) || !arnoldi_finite(coordinate) || !arnoldi_finite(value)) {
        return std::unexpected(core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
    }
    const auto producer = side == ArnoldiSide::Right ? scale({error, error}, 1, d, lambda_scale)
                                                     : scale({error, error}, d, lambda_scale);
    const auto high_magnitude = magnitude(hi), low_magnitude = magnitude(lo);
    if (!high_magnitude) {
        return std::unexpected(high_magnitude.error());
    }
    if (!low_magnitude) {
        return std::unexpected(low_magnitude.error());
    }
    const auto magnitude_input = add_magnitudes(*high_magnitude, *low_magnitude);
    const auto normalized = side == ArnoldiSide::Right ? scale(magnitude_input, 1, d, lambda_scale)
                                                       : scale(magnitude_input, d, lambda_scale);
    constexpr auto epsilon = std::numeric_limits<double>::epsilon();
    // Each complex scaling is real and diagonal. Six rounded real operations,
    // L1 bounding the real/imaginary errors, fit in this conservative gamma.
    const auto rounding = multiply_magnitudes(normalized, {16 * epsilon, 16 * epsilon});
    upper = add_magnitudes(producer, rounding).upper;
    upper = add_magnitudes({upper, upper}, {16 * std::numeric_limits<double>::denorm_min(),
                                            16 * std::numeric_limits<double>::denorm_min()})
                .upper;
    if (!std::isfinite(upper)) {
        return std::unexpected(core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
    }
    return {};
}

template <StorageScalar Scalar>
[[nodiscard]] ArnoldiResult<void>
arnoldi_apply_basis(ArnoldiRunner<Scalar>& runner, ArnoldiSide side, const Eigen::MatrixXcd& basis,
                    Eigen::MatrixXcd& action, Eigen::MatrixXd& error) {
    for (Eigen::Index j = 0; j < basis.cols(); ++j) {
        for (Eigen::Index i = 0; i < basis.rows(); ++i) {
            auto converted = arnoldi_physical_input(
                basis(i, j), runner.scales.metric.components[static_cast<std::size_t>(i)], side,
                runner.input.values()[static_cast<std::size_t>(i)],
                runner.input_low.values()[static_cast<std::size_t>(i)]);
            if (!converted) {
                return runner.remember(runner.failure(converted.error()));
            }
        }
        if (auto r = runner.original_complex(side, runner.input.values(), runner.high.values(),
                                             runner.low.values(), runner.error.values(),
                                             runner.input_low.values());
            !r) {
            return r;
        }
        for (Eigen::Index i = 0; i < basis.rows(); ++i) {
            auto transformed = arnoldi_scaled_action(
                runner.high.values()[static_cast<std::size_t>(i)],
                runner.low.values()[static_cast<std::size_t>(i)],
                runner.error.values()[static_cast<std::size_t>(i)],
                runner.scales.metric.components[static_cast<std::size_t>(i)],
                runner.scales.metric.eigenvalue_scale, side, action(i, j), error(i, j));
            if (!transformed) {
                return runner.remember(runner.failure(transformed.error()));
            }
        }
    }
    return {};
}

// Prefix ordering is obtained by adjacent unitary Schur exchanges from fixed
// Eigen's matrix-functions implementation. Every selected root remains present.
[[nodiscard]] inline bool arnoldi_reorder(Eigen::MatrixXcd& q, Eigen::MatrixXcd& t,
                                          const Eigen::VectorXcd& wanted, ArnoldiAssociation policy,
                                          double physical_scale) {
    Eigen::VectorXi permutation(t.rows());
    permutation.setConstant(-1);
    for (Eigen::Index j = 0; j < wanted.size(); ++j) {
        Eigen::Index match = -1;
        double distance = std::numeric_limits<double>::infinity();
        for (Eigen::Index k = 0; k < t.rows(); ++k) {
            if (permutation[k] >= 0) {
                continue;
            }
            const auto delta = std::abs(t(k, k) - wanted[j]);
            if (arnoldi_near(t(k, k) * physical_scale, wanted[j] * physical_scale, policy,
                             physical_scale) &&
                delta < distance) {
                match = k;
                distance = delta;
            }
        }
        if (match < 0) {
            return false;
        }
        permutation[match] = static_cast<int>(j);
    }
    Eigen::Index next = wanted.size();
    for (Eigen::Index k = 0; k < permutation.size(); ++k) {
        if (permutation[k] < 0) {
            permutation[k] = static_cast<int>(next++);
        }
    }
    Eigen::internal::matrix_function_permute_schur(permutation, q, t);
    return true;
}

template <StorageScalar Scalar, class NativeSolver>
[[nodiscard]] ArnoldiResult<void> arnoldi_cluster(ArnoldiRunner<Scalar>& runner, ArnoldiSide side,
                                                  NativeSolver& native, ArnoldiPacket& packet,
                                                  core::StageTiming& timing, bool timed) {
    const ArnoldiPhase measured(timed, timing);
    const auto m = native.actual_dimension();
    const auto k = packet.count;
    if (m < k || m > runner.request.ncv) {
        return runner.fail(ArnoldiCode::RankUnresolved);
    }
    Eigen::MatrixXcd reduced =
        native.hessenberg_storage().topLeftCorner(m, m).template cast<ArnoldiComplex>();
    Eigen::ComplexSchur<Eigen::MatrixXcd> small(reduced);
    if (small.info() != Eigen::Success) {
        return runner.fail(ArnoldiCode::EigenCoverageUnresolved);
    }
    Eigen::MatrixXcd rotations = small.matrixU(), triangular = small.matrixT();
    Eigen::VectorXcd wanted(k);
    Eigen::VectorXcd selected_physical(k);
    for (core::Index j = 0; j < k; ++j) {
        auto lambda = packet.eigenvalues.values()[static_cast<std::size_t>(j)];
        selected_physical[j] = lambda;
        if (side == ArnoldiSide::Left) {
            lambda = std::conj(lambda);
        }
        if (const auto* shift = std::get_if<ArnoldiShiftInvert>(&runner.request.selection)) {
            const auto sigma = side == ArnoldiSide::Right ? shift->shift : std::conj(shift->shift);
            if (lambda == sigma) {
                return runner.fail(ArnoldiCode::PrecisionUnresolved);
            }
            wanted[j] = runner.scales.metric.eigenvalue_scale / (lambda - sigma);
        } else {
            wanted[j] = lambda / runner.scales.metric.eigenvalue_scale;
        }
    }
    if (!arnoldi_reorder(rotations, triangular, wanted, runner.request.association, 1)) {
        return runner.fail(ArnoldiCode::PrecisionUnresolved);
    }
    Eigen::MatrixXcd basis =
        native.active_basis_storage().leftCols(m).template cast<ArnoldiComplex>() *
        rotations.leftCols(k);
    Eigen::MatrixXcd original_action(runner.data.order, k);
    Eigen::MatrixXd original_error(runner.data.order, k);
    if (auto r = arnoldi_apply_basis(runner, side, basis, original_action, original_error); !r) {
        return r;
    }
    // The restriction is reconstructed from the original action, never from the
    // transformed Ritz values. A second Schur makes that restriction triangular.
    reduced = basis.adjoint() * original_action;
    small.compute(reduced);
    if (small.info() != Eigen::Success) {
        return runner.fail(ArnoldiCode::EigenCoverageUnresolved);
    }
    rotations = small.matrixU();
    triangular = small.matrixT();
    Eigen::VectorXi matched(k);
    matched.setConstant(-1);
    for (core::Index j = 0; j < k; ++j) {
        const auto physical =
            (side == ArnoldiSide::Right ? triangular(j, j) : std::conj(triangular(j, j))) *
            runner.scales.metric.eigenvalue_scale;
        core::Index match = -1;
        for (core::Index l = 0; l < k; ++l) {
            if (matched[l] < 0 &&
                arnoldi_near(physical, selected_physical[l], runner.request.association,
                             runner.scales.metric.eigenvalue_scale)) {
                match = l;
                break;
            }
        }
        if (match < 0) {
            return runner.fail(ArnoldiCode::PrecisionUnresolved);
        }
        matched[match] = static_cast<int>(j);
    }
    basis = (basis * rotations).eval();
    auto& output = side == ArnoldiSide::Right ? packet.right : packet.left;
    auto& t = side == ArnoldiSide::Right ? packet.right_t : packet.left_t;
    auto& mapping = side == ArnoldiSide::Right ? packet.right_mapping : packet.left_mapping;
    for (core::Index j = 0; j < k; ++j) {
        for (core::Index i = 0; i < runner.data.order; ++i) {
            output.values()[static_cast<std::size_t>(j * runner.data.order + i)] = basis(i, j);
        }
        auto physical = triangular(j, j) * runner.scales.metric.eigenvalue_scale;
        if (side == ArnoldiSide::Left) {
            physical = std::conj(physical);
        }
        if (!arnoldi_finite(physical)) {
            return runner.fail(ArnoldiCode::PrecisionUnresolved);
        }
        packet.eigenvalues.values()[static_cast<std::size_t>(j)] = physical;
        mapping.values()[static_cast<std::size_t>(j)] = j;
        for (core::Index i = 0; i < k; ++i) {
            t.values()[static_cast<std::size_t>(j * k + i)] =
                i <= j ? triangular(i, j) : ArnoldiComplex{};
        }
    }
    return {};
}
} // namespace uvlas::numerics::eigen::detail
