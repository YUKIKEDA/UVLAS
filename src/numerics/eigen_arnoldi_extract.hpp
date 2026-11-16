#pragma once
#include "eigen_arnoldi_engine.hpp"

#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>

#include <chrono>

namespace uvlas::numerics::eigen::detail {

// Phase timers have no clock access at all in disabled mode.
class ArnoldiPhase {
public:
    ArnoldiPhase(bool enabled, core::StageTiming& destination) noexcept
        : enabled_(enabled), destination_(destination) {
        if (enabled_) {
            start_ = std::chrono::steady_clock::now();
        }
    }
    ~ArnoldiPhase() noexcept {
        if (enabled_) {
            const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                     std::chrono::steady_clock::now() - start_)
                                     .count();
            ++destination_.calls;
            destination_.elapsed_nanoseconds += static_cast<std::uint64_t>(elapsed);
        }
    }
    ArnoldiPhase(const ArnoldiPhase&) = delete;
    ArnoldiPhase& operator=(const ArnoldiPhase&) = delete;

private:
    bool enabled_;
    core::StageTiming& destination_;
    std::chrono::steady_clock::time_point start_{};
};

[[nodiscard]] inline bool arnoldi_near(ArnoldiComplex a, ArnoldiComplex b,
                                       ArnoldiAssociation policy, double physical_scale) noexcept {
    if (!arnoldi_finite(a) || !arnoldi_finite(b) || !(physical_scale > 0) ||
        !std::isfinite(physical_scale)) {
        return false;
    }
    if (policy.absolute == 0 && policy.relative == 0) {
        return a == b;
    }
    const auto divisor = std::max({std::abs(a.real()), std::abs(a.imag()), std::abs(b.real()),
                                   std::abs(b.imag()), physical_scale});
    const auto x = a / divisor, y = b / divisor;
    const auto distance = std::abs(x - y);
    const auto allowance = policy.absolute * (physical_scale / divisor) +
                           policy.relative * std::max(std::abs(x), std::abs(y));
    // Normalize physical eigenvalues and the eigenvalue scale together, avoiding
    // overflow in mu=lambda/s_lambda, a-b, and a+b*max(abs(mu)).
    if (distance == 0 && a != b && allowance == 0) {
        return false;
    }
    return distance <= allowance;
}

// Sorting has bounded, explicit accesses; all payload and column order move
// together. It does not allocate and never examines capacity beyond count.
inline void arnoldi_sort(ArnoldiPacket& packet, core::Index n, bool subspace) noexcept {
    for (core::Index i = 1; i < packet.count; ++i) {
        for (core::Index j = i; j > 0; --j) {
            auto& a = packet.eigenvalues.values()[static_cast<std::size_t>(j - 1)];
            auto& b = packet.eigenvalues.values()[static_cast<std::size_t>(j)];
            if (a.real() < b.real() || (a.real() == b.real() && a.imag() <= b.imag())) {
                break;
            }
            std::swap(a, b);
            if (subspace) {
                std::swap(packet.right_mapping.values()[static_cast<std::size_t>(j - 1)],
                          packet.right_mapping.values()[static_cast<std::size_t>(j)]);
                continue;
            }
            for (core::Index row = 0; row < n; ++row) {
                std::swap(packet.right.values()[static_cast<std::size_t>((j - 1) * n + row)],
                          packet.right.values()[static_cast<std::size_t>(j * n + row)]);
            }
        }
    }
}

[[nodiscard]] core::Result<void> arnoldi_normalize(std::span<ArnoldiComplex> vector,
                                                   std::span<const double> scales, ArnoldiSide side,
                                                   std::span<ArnoldiComplex> backup,
                                                   std::span<ArnoldiComplex> candidate);

[[nodiscard]] inline core::Result<void> arnoldi_dimensionless_coordinate(ArnoldiComplex physical,
                                                                         double d, ArnoldiSide side,
                                                                         ArnoldiComplex& center,
                                                                         double& error) noexcept {
    error = 0;
    for (unsigned lane = 0; lane < 2; ++lane) {
        const auto input = lane == 0 ? physical.real() : physical.imag();
        double value = 0, upper = 0;
        if (side == ArnoldiSide::Left) {
            const auto product = two_product(input, d);
            if (!product) {
                return std::unexpected(
                    core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
            }
            value = product->hi;
            upper = std::abs(product->lo);
        } else {
            value = input / d;
            if (!std::isfinite(value)) {
                return std::unexpected(
                    core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
            }
            Sum<double> difference;
            const auto first = difference.add(input), second = difference.product(value, d, -1);
            if (!first || !second) {
                return std::unexpected(
                    core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
            }
            const auto evaluated = difference.finish();
            if (!evaluated) {
                return std::unexpected(
                    core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
            }
            const auto magnitude_value = magnitude(evaluated->value);
            if (!magnitude_value) {
                return std::unexpected(magnitude_value.error());
            }
            upper = divide_bound(
                add_magnitudes(*magnitude_value, {evaluated->error, evaluated->error}).upper, d,
                true);
        }
        if (lane == 0) {
            center.real(value);
        } else {
            center.imag(value);
        }
        error = add_magnitudes({error, error}, {upper, upper}).upper;
    }
    if (!std::isfinite(error)) {
        return std::unexpected(core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
    }
    return {};
}

// Only a certified Gram error below one establishes an independent conjugate
// group. A tiny imaginary rounding tail must not duplicate one real eigenmode.
[[nodiscard]] inline core::Result<void>
arnoldi_conjugate_independence(std::span<const ArnoldiComplex> a, std::span<const ArnoldiComplex> b,
                               std::span<const double> scales, ArnoldiSide side) noexcept {
    Sum<double> norm_a, norm_b;
    Sum<ArnoldiComplex> dot;
    double error_a = 0, error_b = 0, error_dot = 0;
    const auto add_error = [](double& total, double contribution) {
        total = add_magnitudes({total, total}, {contribution, contribution}).upper;
    };
    for (std::size_t i = 0; i < a.size(); ++i) {
        ArnoldiComplex x{}, y{};
        double xerror = 0, yerror = 0;
        const auto converted_a = arnoldi_dimensionless_coordinate(a[i], scales[i], side, x, xerror);
        const auto converted_b = arnoldi_dimensionless_coordinate(b[i], scales[i], side, y, yerror);
        if (!converted_a) {
            return converted_a;
        }
        if (!converted_b) {
            return converted_b;
        }
        for (const auto value : {x.real(), x.imag()}) {
            if (!norm_a.product(value, value)) {
                return std::unexpected(
                    core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
            }
        }
        for (const auto value : {y.real(), y.imag()}) {
            if (!norm_b.product(value, value)) {
                return std::unexpected(
                    core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
            }
        }
        if (!dot.product(std::conj(x), y)) {
            return std::unexpected(
                core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
        }
        const auto xm = magnitude(x), ym = magnitude(y);
        if (!xm || !ym) {
            return std::unexpected(
                core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
        }
        const Bounds xe{xerror, xerror}, ye{yerror, yerror};
        add_error(error_a, multiply_magnitudes(*xm, xe).upper);
        add_error(error_a, multiply_magnitudes(*xm, xe).upper);
        add_error(error_a, multiply_magnitudes(xe, xe).upper);
        add_error(error_b, multiply_magnitudes(*ym, ye).upper);
        add_error(error_b, multiply_magnitudes(*ym, ye).upper);
        add_error(error_b, multiply_magnitudes(ye, ye).upper);
        add_error(error_dot, multiply_magnitudes(*xm, ye).upper);
        add_error(error_dot, multiply_magnitudes(*ym, xe).upper);
        add_error(error_dot, multiply_magnitudes(xe, ye).upper);
    }
    if (!norm_a.add(-1) || !norm_b.add(-1)) {
        return std::unexpected(core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
    }
    const auto left = norm_a.finish(), right = norm_b.finish();
    const auto overlap = dot.finish();
    if (!left || !right || !overlap) {
        return std::unexpected(core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
    }
    const auto lm = magnitude(left->value), rm = magnitude(right->value),
               dm = magnitude(overlap->value);
    if (!lm || !rm || !dm) {
        return std::unexpected(core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
    }
    add_error(error_a, left->error);
    add_error(error_b, right->error);
    add_error(error_dot, overlap->error);
    const auto diagonal_a = add_magnitudes(*lm, {error_a, error_a}),
               diagonal_b = add_magnitudes(*rm, {error_b, error_b});
    const auto off_diagonal = add_magnitudes(*dm, {error_dot, error_dot});
    const auto gram = std::max(add_magnitudes(diagonal_a, off_diagonal).upper,
                               add_magnitudes(diagonal_b, off_diagonal).upper);
    if (!(gram < 1)) {
        return std::unexpected(core::Failure{ArnoldiCode::PrecisionUnresolved, arnoldi_operation});
    }
    return {};
}

template <StorageScalar Scalar, class NativeSolver>
[[nodiscard]] ArnoldiResult<void> arnoldi_cluster(ArnoldiRunner<Scalar>&, ArnoldiSide,
                                                  NativeSolver&, ArnoldiPacket&, core::StageTiming&,
                                                  bool);

template <StorageScalar Original, StorageScalar Native>
[[nodiscard]] ArnoldiResult<void>
arnoldi_extract_side(ArnoldiRunner<Original>& runner, ArnoldiSide side, ArnoldiPacket& destination,
                     core::StageTiming& phase, bool timed, core::StageTiming& subspace_phase) {
    auto started = runner.session.start_attempt();
    if (!started) {
        return runner.remember(started.error());
    }
    ArnoldiOperation<Original, Native> operation{runner, side};
    ArnoldiNative<decltype(operation)> solver(operation, runner.request.nev, runner.request.ncv);
    // The exact explicit seed is converted only when a real original problem has
    // a complex shift. No backend random initializer is ever used.
    Eigen::Matrix<Native, Eigen::Dynamic, 1> seed(runner.data.order);
    for (Eigen::Index i = 0; i < seed.size(); ++i) {
        if constexpr (std::same_as<Native, double>) {
            seed[i] = runner.data.initial.values()[static_cast<std::size_t>(i)].real();
        } else {
            seed[i] = runner.data.initial.values()[static_cast<std::size_t>(i)];
        }
    }
    double largest_seed = 0;
    for (Eigen::Index i = 0; i < seed.size(); ++i) {
        largest_seed = std::max(largest_seed, std::abs(seed[i]));
    }
    double working_norm = 0;
    for (Eigen::Index i = 0; i < seed.size(); ++i) {
        seed[i] /= largest_seed;
        working_norm = std::hypot(working_norm, std::abs(seed[i]));
    }
    for (Eigen::Index i = 0; i < seed.size(); ++i) {
        seed[i] /= working_norm;
    }
    // Spectra's internal norm may square its inputs. Only this working copy is
    // safely normalized; the exact caller bits and original norm remain owned.
    Eigen::Index converged = 0;
    {
        const ArnoldiPhase measured(timed, phase);
        try {
            solver.init(seed.data());
            converged =
                solver.compute(Spectra::SortRule::LargestMagn,
                               static_cast<Eigen::Index>(runner.request.backend_iterations),
                               runner.request.backend_tolerance, Spectra::SortRule::LargestMagn);
        } catch (...) {
            // A failed provider can make Spectra reject its dummy action later.
            // That secondary exception must never replace the first actual cause.
            if (runner.first) {
                return std::unexpected(*runner.first);
            }
            throw;
        }
    }
    runner.data.diagnostics.backend_reported_operations +=
        static_cast<std::uint64_t>(solver.num_operations());
    runner.data.diagnostics.attempts.backend_iterations =
        runner.data.diagnostics.attempts.backend_iterations.value_or(0) +
        static_cast<std::uint64_t>(solver.num_iterations());
    runner.data.diagnostics.attempts.converged_count =
        runner.data.diagnostics.attempts.converged_count.value_or(0) +
        static_cast<std::uint64_t>(converged);
    if (runner.first) {
        *runner.first = complete_failure(*runner.first, runner.scales.context, runner.request.state,
                                         runner.counters());
        return std::unexpected(*runner.first);
    }
    if (side == ArnoldiSide::Right) {
        runner.data.diagnostics.right_krylov_dimension = solver.actual_dimension();
    } else {
        runner.data.diagnostics.left_krylov_dimension = solver.actual_dimension();
    }
    if (solver.info() != Spectra::CompInfo::Successful || converged < runner.request.nev) {
        return runner.fail(ArnoldiCode::EigenCoverageUnresolved);
    }
    const auto values = solver.eigenvalues();
    const auto count = static_cast<core::Index>(values.size());
    if (count != runner.request.nev) {
        return runner.fail(ArnoldiCode::EigenCoverageUnresolved);
    }
    Eigen::MatrixXcd vectors;
    if (!runner.data.subspace_requested) {
        vectors = solver.eigenvectors();
        if (vectors.rows() != runner.data.order || vectors.cols() != count) {
            return runner.fail(ArnoldiCode::EigenCoverageUnresolved);
        }
    }
    auto& output = side == ArnoldiSide::Right ? destination.right : destination.left;
    // Keep this complete stage chain: native values -> original physical values
    // -> normalized actual phase -> conjugate closure -> original qualification.
    // A subspace uses active V/H, never the rank of individual Ritz vectors.
    for (core::Index j = 0; j < count; ++j) {
        auto physical = values[j];
        if (const auto* shift = std::get_if<ArnoldiShiftInvert>(&runner.request.selection)) {
            if (physical == ArnoldiComplex{}) {
                return runner.fail(ArnoldiCode::PrecisionUnresolved);
            }
            physical = (side == ArnoldiSide::Right ? shift->shift : std::conj(shift->shift)) +
                       runner.scales.metric.eigenvalue_scale / physical;
        } else {
            physical *= runner.scales.metric.eigenvalue_scale;
        }
        if (side == ArnoldiSide::Left) {
            physical = std::conj(physical);
        }
        if (!arnoldi_finite(physical)) {
            return runner.fail(ArnoldiCode::PrecisionUnresolved);
        }
        destination.eigenvalues.values()[static_cast<std::size_t>(j)] = physical;
        if (!runner.data.subspace_requested) {
            for (core::Index i = 0; i < runner.data.order; ++i) {
                const auto scale = runner.scales.metric.components[static_cast<std::size_t>(i)];
                output.values()[static_cast<std::size_t>(j * runner.data.order + i)] =
                    side == ArnoldiSide::Right ? vectors(i, j) * scale : vectors(i, j) / scale;
            }
            const auto n = static_cast<std::size_t>(runner.data.order);
            auto column = output.values().subspan(static_cast<std::size_t>(j) * n, n);
            const auto normalized = arnoldi_normalize(column, runner.scales.metric.components, side,
                                                      runner.input.values(), runner.high.values());
            if (!normalized) {
                return runner.remember(runner.failure(normalized.error()));
            }
        }
    }
    destination.count = count;
    if constexpr (std::same_as<Original, double>) {
        const auto n = static_cast<std::size_t>(runner.data.order);
        std::ranges::copy(destination.eigenvalues.values().first(static_cast<std::size_t>(count)),
                          runner.candidate_values.values().begin());
        std::ranges::copy(output.values().first(n * static_cast<std::size_t>(count)),
                          runner.candidate_vectors.values().begin());
        // Each native candidate belongs to at most one selected original group.
        // Preserve every native selected root, completing its entire group.
        // If closure exceeds nev+1, reject rather than discard a selected root.
        std::fill(runner.real_error.values().begin(), runner.real_error.values().end(), 0.0);
        destination.count = 0;
        for (core::Index j = 0; j < count; ++j) {
            if (runner.real_error.values()[static_cast<std::size_t>(j)] != 0) {
                continue;
            }
            const auto value = runner.candidate_values.values()[static_cast<std::size_t>(j)];
            if (std::cmp_greater_equal(destination.count, destination.eigenvalues.size())) {
                return runner.fail(ArnoldiCode::EigenCoverageUnresolved);
            }
            const auto first = destination.count++;
            destination.eigenvalues.values()[static_cast<std::size_t>(first)] = value;
            std::ranges::copy(
                runner.candidate_vectors.values().subspan(static_cast<std::size_t>(j) * n, n),
                output.values().subspan(static_cast<std::size_t>(first) * n, n).begin());
            runner.real_error.values()[static_cast<std::size_t>(j)] = 1;
            if (value.imag() == 0) {
                continue;
            }
            core::Index partner = -1;
            for (core::Index l = 0; l < count; ++l) {
                const auto candidate =
                    runner.candidate_values.values()[static_cast<std::size_t>(l)];
                if (runner.real_error.values()[static_cast<std::size_t>(l)] == 0 &&
                    candidate.imag() != 0 &&
                    std::signbit(candidate.imag()) != std::signbit(value.imag()) &&
                    arnoldi_near(std::conj(value), candidate, runner.request.association,
                                 runner.scales.metric.eigenvalue_scale)) {
                    partner = l;
                    break;
                }
            }
            if (std::cmp_greater_equal(destination.count, destination.eigenvalues.size())) {
                return runner.fail(ArnoldiCode::EigenCoverageUnresolved);
            }
            const auto second = destination.count++;
            destination.eigenvalues.values()[static_cast<std::size_t>(second)] =
                partner >= 0 ? runner.candidate_values.values()[static_cast<std::size_t>(partner)]
                             : std::conj(value);
            for (std::size_t i = 0; i < n; ++i) {
                output.values()[static_cast<std::size_t>(second) * n + i] =
                    partner >= 0
                        ? runner.candidate_vectors
                              .values()[static_cast<std::size_t>(partner) * n + i]
                        : std::conj(output.values()[static_cast<std::size_t>(first) * n + i]);
            }
            if (partner >= 0) {
                runner.real_error.values()[static_cast<std::size_t>(partner)] = 1;
            }
            if (!runner.data.subspace_requested) {
                const auto independent = arnoldi_conjugate_independence(
                    output.values().subspan(static_cast<std::size_t>(first) * n, n),
                    output.values().subspan(static_cast<std::size_t>(second) * n, n),
                    runner.scales.metric.components, side);
                if (!independent) {
                    return runner.remember(runner.failure(independent.error()));
                }
            }
        }
    }
    if (runner.data.subspace_requested) {
        return arnoldi_cluster(runner, side, solver, destination, subspace_phase, timed);
    }
    return {};
}
} // namespace uvlas::numerics::eigen::detail

#include "eigen_arnoldi_cluster.hpp"
