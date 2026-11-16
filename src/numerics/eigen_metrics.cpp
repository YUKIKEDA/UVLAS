#include "eigen_metric_bounds.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <uvlas/numerics/eigen_metrics.hpp>

namespace uvlas::numerics::eigen {
namespace {
using core::FailureCode;
using detail::Bounds;
constexpr auto operation = core::Operation::ExtractEigen;

auto fail(FailureCode code) noexcept {
    return std::unexpected(core::Failure{code, operation});
}
auto fail(const MetricContext& context, FailureCode code) noexcept {
    core::FailureDetails details;
    details.numerical_generation = context.values_generation;
    return std::unexpected(
        complete_failure(core::Failure{code, operation, details}, context.system, context.state));
}
template <class T>
auto forward_failure(const MetricContext& context, const core::Result<T>& result) noexcept {
    auto cause = result.error();
    if (!cause.details.numerical_generation)
        cause.details.numerical_generation = context.values_generation;
    return std::unexpected(complete_failure(cause, context.system, context.state));
}
template <StorageScalar Scalar>
Scalar low(WordView<Scalar> view, std::size_t index) noexcept {
    return view.low().empty() ? Scalar{} : view.low()[index];
}
template <StorageScalar Scalar>
Scalar conjugate(Scalar value) noexcept {
    if constexpr (std::same_as<Scalar, double>)
        return value;
    else
        return std::conj(value);
}
template <StorageScalar Scalar>
core::Result<void> validate(WordView<Scalar> view) noexcept {
    // A borrowed view may have been made earlier; validate the actual call inputs.
    const auto checked = WordView<Scalar>::make(view.high(), view.low(), view.precision());
    if (!checked)
        return std::unexpected(checked.error());
    return {};
}
template <StorageScalar Scalar>
core::Result<void> validate(ActionValues<Scalar> action) noexcept {
    const auto words = validate(action.values);
    if (!words)
        return words;
    if (action.absolute_error_upper.size() != action.values.size())
        return fail(FailureCode::InvalidInput);
    for (const double error : action.absolute_error_upper) {
        // Nonfinite supplied bounds are invalid external data. Overflow of a
        // bound computed below is instead an unresolved internal enclosure.
        if (!std::isfinite(error))
            return fail(FailureCode::Nonfinite);
        if (error < 0)
            return fail(FailureCode::InvalidInput);
    }
    return {};
}
core::Result<void> validate(MetricRequirement requirement) noexcept {
    if (!std::isfinite(requirement.budget))
        return fail(FailureCode::Nonfinite);
    if (requirement.budget < 0)
        return fail(FailureCode::InvalidInput);
    return {};
}
core::Result<void> validate(const BackwardRequirement& request) noexcept {
    if (request.valueless_by_exception())
        return fail(FailureCode::InvalidInput);
    const auto* required = std::get_if<MetricRequirement>(&request);
    return required ? validate(*required) : core::Result<void>{};
}
core::Result<void> validate(MetricScales scales, std::size_t rows) noexcept {
    if (scales.components.size() != rows)
        return fail(FailureCode::InvalidInput);
    if (!std::isfinite(scales.eigenvalue_scale))
        return fail(FailureCode::Nonfinite);
    if (scales.eigenvalue_scale <= 0)
        return fail(FailureCode::InvalidInput);
    for (const double scale : scales.components) {
        if (!std::isfinite(scale))
            return fail(FailureCode::Nonfinite);
        if (scale <= 0)
            return fail(FailureCode::InvalidInput);
    }
    return {};
}
core::Result<void> validate(MetricRequirement residual, MetricRequirement other,
                            const BackwardRequirement& backward) noexcept {
    for (const auto& requirement : {residual, other}) {
        const auto checked = validate(requirement);
        if (!checked)
            return checked;
    }
    if (residual.budget > 1e-10)
        return fail(FailureCode::InvalidInput);
    return validate(backward);
}
core::Result<double> finite_upper(double value) noexcept {
    if (!std::isfinite(value) || value < 0)
        return fail(FailureCode::PrecisionUnresolved);
    return value;
}
core::Result<core::ErrorCertificate> certificate(MetricRequirement requirement, double upper,
                                                 core::Precision precision,
                                                 const MetricContext& context) noexcept {
    const auto bounded = finite_upper(upper);
    if (!bounded)
        return std::unexpected(bounded.error());
    if (upper > requirement.budget) {
        core::FailureDetails details;
        details.problem = requirement.observation;
        details.numerical_generation = context.values_generation;
        return std::unexpected(core::Failure{FailureCode::PrecisionUnresolved, operation, details});
    }
    return core::ErrorCertificate::assess(
        {requirement.observation, requirement.scope, context.values_generation, precision},
        requirement.budget, upper);
}
core::Result<BackwardCertificate> backward_certificate(const BackwardRequirement& request,
                                                       double upper, core::Precision precision,
                                                       const MetricContext& context) noexcept {
    const auto* required = std::get_if<MetricRequirement>(&request);
    if (!required)
        return BackwardCertificate{NotRequested{}};
    const auto proof = certificate(*required, upper, precision, context);
    if (!proof)
        return std::unexpected(proof.error());
    return BackwardCertificate{*proof};
}
template <StorageScalar Scalar>
core::Result<void> add_words(detail::Sum<Scalar>& sum, WordView<Scalar> words,
                             std::size_t index) noexcept {
    auto high = sum.add(words.high()[index]);
    if (!high)
        return high;
    return sum.add(low(words, index));
}
template <StorageScalar Scalar>
core::Result<void> product_words(detail::Sum<Scalar>& sum, WordView<Scalar> a, std::size_t ai,
                                 WordView<Scalar> b, std::size_t bi, double sign,
                                 bool conjugate_first = false) noexcept {
    const auto transform = [&](Scalar value) { return conjugate_first ? conjugate(value) : value; };
    for (const Scalar left : {a.high()[ai], low(a, ai)}) {
        for (const Scalar right : {b.high()[bi], low(b, bi)}) {
            const auto product = sum.product(transform(left), right, sign);
            if (!product)
                return product;
        }
    }
    return {};
}
template <StorageScalar Scalar>
core::Result<Bounds> word_magnitude(WordView<Scalar> view, std::size_t index) noexcept {
    return detail::magnitude(view.high()[index], low(view, index));
}
template <StorageScalar Scalar>
core::Result<double> residual_upper(detail::Sum<Scalar>& sum, double action_error) noexcept {
    const auto evaluated = sum.finish();
    if (!evaluated)
        return std::unexpected(evaluated.error());
    const auto absolute = detail::magnitude(evaluated->value);
    if (!absolute)
        return std::unexpected(absolute.error());
    return finite_upper(detail::add_magnitudes(
                            *absolute, detail::add_magnitudes({evaluated->error, evaluated->error},
                                                              {action_error, action_error}))
                            .upper);
}
template <StorageScalar Scalar>
core::Result<double> action_lower(ActionValues<Scalar> action, std::size_t index) noexcept {
    const auto absolute = word_magnitude(action.values, index);
    if (!absolute)
        return std::unexpected(absolute.error());
    return detail::lower_subtract(absolute->lower, action.absolute_error_upper[index]);
}
template <StorageScalar Scalar>
core::Result<Bounds> scaled_vector_magnitude(WordView<Scalar> vector, std::size_t index,
                                             MetricScales scales, PairSide side) noexcept {
    const auto absolute = word_magnitude(vector, index);
    if (!absolute)
        return std::unexpected(absolute.error());
    return side == PairSide::Right ? detail::scale(*absolute, 1, scales.components[index])
                                   : detail::scale(*absolute, scales.components[index]);
}
template <StorageScalar Scalar>
core::Precision representation(WordView<Scalar> a, WordView<Scalar> b) noexcept {
    return a.precision() == core::Precision::TwoWord || b.precision() == core::Precision::TwoWord
               ? core::Precision::TwoWord
               : core::Precision::Binary64;
}
} // namespace

template <StorageScalar Scalar>
Result<PairMetrics> assess_pair_metrics(WordView<Scalar> vector, WordView<Scalar> eigenvalue,
                                        ActionValues<Scalar> action, MetricScales scales,
                                        PairSide side, const PairMetricRequest& request,
                                        const MetricContext& context) noexcept {
    return guard_entry<PairMetrics>(
        context.system, context.state, {}, [&]() -> Result<PairMetrics> {
            if (!std::isfinite(context.system.time))
                return fail(context, FailureCode::Nonfinite);
            if (context.system.cache.representation != core::Precision::Binary64 &&
                context.system.cache.representation != core::Precision::TwoWord)
                return fail(context, FailureCode::PrecisionUnsupported);
            if (vector.size() == 0 || eigenvalue.size() != 1 ||
                action.values.size() != vector.size() ||
                (side != PairSide::Right && side != PairSide::Left))
                return fail(context, FailureCode::InvalidInput);
            for (const auto& words : {vector, eigenvalue}) {
                const auto checked = validate(words);
                if (!checked)
                    return forward_failure(context, checked);
            }
            const auto checked_action = validate(action);
            if (!checked_action)
                return forward_failure(context, checked_action);
            const auto checked_scales = validate(scales, vector.size());
            if (!checked_scales)
                return forward_failure(context, checked_scales);
            const auto checked_request =
                validate(request.residual, request.normalization, request.backward);
            if (!checked_request)
                return forward_failure(context, checked_request);

            double infinity_lower = 0, residual = 0, backward = 0;
            Bounds norm{0, 0};
            const auto lambda_absolute = word_magnitude(eigenvalue, 0);
            if (!lambda_absolute)
                return forward_failure(context, lambda_absolute);
            const bool need_backward = std::holds_alternative<MetricRequirement>(request.backward);
            for (std::size_t i = 0; i < vector.size(); ++i) {
                const auto coordinate = scaled_vector_magnitude(vector, i, scales, side);
                if (!coordinate)
                    return forward_failure(context, coordinate);
                infinity_lower = std::max(infinity_lower, coordinate->lower);
                const auto next_norm = detail::norm_step(norm, *coordinate);
                if (!next_norm)
                    return forward_failure(context, next_norm);
                norm = *next_norm;
                detail::Sum<Scalar> sum;
                auto added = add_words(sum, action.values, i);
                if (!added)
                    return forward_failure(context, added);
                added = product_words(sum, eigenvalue, 0, vector, i, -1, side == PairSide::Left);
                if (!added)
                    return forward_failure(context, added);
                const auto upper = residual_upper(sum, action.absolute_error_upper[i]);
                if (!upper)
                    return forward_failure(context, upper);
                const Bounds physical{*upper, *upper};
                const auto scaled =
                    side == PairSide::Right
                        ? detail::scale(physical, 1, scales.components[i], scales.eigenvalue_scale)
                        : detail::scale(physical, scales.components[i], scales.eigenvalue_scale);
                residual = std::max(residual, scaled.upper);
                if (need_backward) {
                    const auto original = action_lower(action, i);
                    if (!original)
                        return forward_failure(context, original);
                    const auto vector_absolute = word_magnitude(vector, i);
                    if (!vector_absolute)
                        return forward_failure(context, vector_absolute);
                    const double denominator =
                        detail::add_magnitudes(
                            {*original, *original},
                            detail::multiply_magnitudes(*lambda_absolute, *vector_absolute))
                            .lower;
                    backward = std::max(backward, detail::divide_bound(*upper, denominator, true));
                }
            }
            if (!(infinity_lower > 0))
                return fail(context, FailureCode::PrecisionUnresolved);
            const auto precision = representation(vector, eigenvalue);
            const auto residual_proof =
                certificate(request.residual, detail::divide_bound(residual, infinity_lower, true),
                            precision, context);
            if (!residual_proof)
                return forward_failure(context, residual_proof);
            const auto norm_proof = certificate(
                request.normalization, detail::distance_from_one(norm), precision, context);
            if (!norm_proof)
                return forward_failure(context, norm_proof);
            const auto backward_proof =
                backward_certificate(request.backward, backward, precision, context);
            if (!backward_proof)
                return forward_failure(context, backward_proof);
            return PairMetrics{*residual_proof, *norm_proof, *backward_proof, infinity_lower};
        });
}

template <StorageScalar Scalar>
Result<BasisMetrics>
assess_schur_metrics(core::Index rows, core::Index columns, WordView<Scalar> basis,
                     ActionValues<Scalar> triangular, ActionValues<Scalar> action,
                     const BasisMetricRequest& request, const MetricContext& context) noexcept {
    return guard_entry<BasisMetrics>(
        context.system, context.state, {}, [&]() -> Result<BasisMetrics> {
            if (!std::isfinite(context.system.time))
                return fail(context, FailureCode::Nonfinite);
            if (context.system.cache.representation != core::Precision::Binary64 &&
                context.system.cache.representation != core::Precision::TwoWord)
                return fail(context, FailureCode::PrecisionUnsupported);
            if (rows < 0 || columns < 0 || columns > rows)
                return fail(context, FailureCode::InvalidInput);
            const auto count = core::checked_extent_product(std::array{rows, columns}, operation);
            if (!count)
                return forward_failure(context, count);
            const auto square =
                core::checked_extent_product(std::array{columns, columns}, operation);
            if (!square)
                return forward_failure(context, square);
            if (std::cmp_not_equal(basis.size(), *count) || action.values.size() != basis.size() ||
                std::cmp_not_equal(triangular.values.size(), *square))
                return fail(context, FailureCode::InvalidInput);
            const auto checked_basis = validate(basis);
            if (!checked_basis)
                return forward_failure(context, checked_basis);
            for (const auto& candidate : {triangular, action}) {
                const auto checked = validate(candidate);
                if (!checked)
                    return forward_failure(context, checked);
            }
            const auto checked_request =
                validate(request.residual, request.orthogonality, request.backward);
            if (!checked_request)
                return forward_failure(context, checked_request);
            const auto n = static_cast<std::size_t>(rows), k = static_cast<std::size_t>(columns);
            for (std::size_t j = 0; j < k; ++j) {
                for (std::size_t i = j + 1; i < k; ++i) {
                    const auto absolute = word_magnitude(triangular.values, i + j * k);
                    if (!absolute)
                        return forward_failure(context, absolute);
                    if (absolute->upper != 0 || triangular.absolute_error_upper[i + j * k] != 0)
                        return fail(context, FailureCode::PrecisionUnresolved);
                }
            }
            double basis_norm_lower = 0, residual = 0, backward = 0, gram = 0;
            const bool need_backward = std::holds_alternative<MetricRequirement>(request.backward);
            for (std::size_t i = 0; i < n; ++i) {
                double row_lower = 0, row_error = 0;
                for (std::size_t j = 0; j < k; ++j) {
                    const auto coordinate = word_magnitude(basis, i + j * n);
                    if (!coordinate)
                        return forward_failure(context, coordinate);
                    row_lower = detail::add_magnitudes({row_lower, row_lower}, *coordinate).lower;
                    detail::Sum<Scalar> sum;
                    auto added = add_words(sum, action.values, i + j * n);
                    if (!added)
                        return forward_failure(context, added);
                    double propagated = action.absolute_error_upper[i + j * n], denominator = 0;
                    if (need_backward) {
                        const auto lower = action_lower(action, i + j * n);
                        if (!lower)
                            return forward_failure(context, lower);
                        denominator = *lower;
                    }
                    for (std::size_t inner = 0; inner < k; ++inner) {
                        added = product_words(sum, basis, i + inner * n, triangular.values,
                                              inner + j * k, -1);
                        if (!added)
                            return forward_failure(context, added);
                        const auto absolute = word_magnitude(basis, i + inner * n);
                        if (!absolute)
                            return forward_failure(context, absolute);
                        const double error = triangular.absolute_error_upper[inner + j * k];
                        propagated = detail::add_magnitudes(
                                         {propagated, propagated},
                                         detail::multiply_magnitudes(*absolute, {error, error}))
                                         .upper;
                        if (need_backward) {
                            const auto t = action_lower(triangular, inner + j * k);
                            if (!t)
                                return forward_failure(context, t);
                            denominator = detail::add_magnitudes(
                                              {denominator, denominator},
                                              detail::multiply_magnitudes(*absolute, {*t, *t}))
                                              .lower;
                        }
                    }
                    const auto upper = residual_upper(sum, propagated);
                    if (!upper)
                        return forward_failure(context, upper);
                    row_error =
                        detail::add_magnitudes({row_error, row_error}, {*upper, *upper}).upper;
                    if (need_backward)
                        backward =
                            std::max(backward, detail::divide_bound(*upper, denominator, true));
                }
                basis_norm_lower = std::max(basis_norm_lower, row_lower);
                residual = std::max(residual, row_error);
            }
            for (std::size_t i = 0; i < k; ++i) {
                double row_error = 0;
                for (std::size_t j = 0; j < k; ++j) {
                    detail::Sum<Scalar> sum;
                    for (std::size_t inner = 0; inner < n; ++inner) {
                        const auto added =
                            product_words(sum, basis, inner + i * n, basis, inner + j * n, 1, true);
                        if (!added)
                            return forward_failure(context, added);
                    }
                    if (i == j) {
                        const auto diagonal = sum.add(Scalar{-1});
                        if (!diagonal)
                            return forward_failure(context, diagonal);
                    }
                    const auto upper = residual_upper(sum, 0);
                    if (!upper)
                        return forward_failure(context, upper);
                    row_error =
                        detail::add_magnitudes({row_error, row_error}, {*upper, *upper}).upper;
                }
                gram = std::max(gram, row_error);
            }
            // An unrepresentable enclosure is not a finite rank diagnostic.
            if (!std::isfinite(gram))
                return fail(context, FailureCode::PrecisionUnresolved);
            if (k != 0 && !(gram < 1))
                return fail(context, FailureCode::RankUnresolved);
            if (k != 0 && !(basis_norm_lower > 0))
                return fail(context, FailureCode::PrecisionUnresolved);
            const auto precision = representation(basis, triangular.values);
            const auto residual_proof =
                certificate(request.residual,
                            k == 0 ? 0 : detail::divide_bound(residual, basis_norm_lower, true),
                            precision, context);
            if (!residual_proof)
                return forward_failure(context, residual_proof);
            const auto gram_proof = certificate(request.orthogonality, gram, precision, context);
            if (!gram_proof)
                return forward_failure(context, gram_proof);
            const auto backward_proof =
                backward_certificate(request.backward, backward, precision, context);
            if (!backward_proof)
                return forward_failure(context, backward_proof);
            return BasisMetrics{*residual_proof, *gram_proof, *backward_proof};
        });
}

Result<GeneralizedMetrics> assess_generalized_metrics(core::Index rows, WordView<double> basis,
                                                      WordView<double> eigenvalues,
                                                      ActionValues<double> stiffness_action,
                                                      ActionValues<double> mass_action,
                                                      MetricScales scales, double mass_scale,
                                                      const GeneralizedMetricRequest& request,
                                                      const MetricContext& context) noexcept {
    return guard_entry<GeneralizedMetrics>(
        context.system, context.state, {}, [&]() -> Result<GeneralizedMetrics> {
            if (!std::isfinite(context.system.time))
                return fail(context, FailureCode::Nonfinite);
            if (context.system.cache.representation != core::Precision::Binary64 &&
                context.system.cache.representation != core::Precision::TwoWord)
                return fail(context, FailureCode::PrecisionUnsupported);
            if (rows < 0)
                return fail(context, FailureCode::InvalidInput);
            const auto square = core::checked_extent_product(std::array{rows, rows}, operation);
            if (!square)
                return forward_failure(context, square);
            if (std::cmp_not_equal(basis.size(), *square) ||
                std::cmp_not_equal(eigenvalues.size(), rows) ||
                stiffness_action.values.size() != basis.size() ||
                mass_action.values.size() != basis.size())
                return fail(context, FailureCode::InvalidInput);
            for (const auto& words : {basis, eigenvalues}) {
                const auto checked = validate(words);
                if (!checked)
                    return forward_failure(context, checked);
            }
            for (const auto& action : {stiffness_action, mass_action}) {
                const auto checked = validate(action);
                if (!checked)
                    return forward_failure(context, checked);
            }
            const auto n = static_cast<std::size_t>(rows);
            const auto checked_scales = validate(scales, n);
            if (!checked_scales)
                return forward_failure(context, checked_scales);
            if (!std::isfinite(mass_scale))
                return fail(context, FailureCode::Nonfinite);
            if (mass_scale <= 0)
                return fail(context, FailureCode::InvalidInput);
            const auto checked_request =
                validate(request.residual, request.orthogonality, request.backward);
            if (!checked_request)
                return forward_failure(context, checked_request);
            const auto checked_normalization = validate(request.normalization);
            if (!checked_normalization)
                return forward_failure(context, checked_normalization);

            double residual = 0, backward = 0, gram = 0, normalization = 0;
            const bool need_backward = std::holds_alternative<MetricRequirement>(request.backward);
            for (std::size_t j = 0; j < n; ++j) {
                double infinity_lower = 0, column_error = 0;
                const auto lambda = word_magnitude(eigenvalues, j);
                if (!lambda)
                    return forward_failure(context, lambda);
                for (std::size_t i = 0; i < n; ++i) {
                    const auto coordinate = word_magnitude(basis, i + j * n);
                    if (!coordinate)
                        return forward_failure(context, coordinate);
                    infinity_lower = std::max(
                        infinity_lower, detail::scale(*coordinate, 1, scales.components[i]).lower);
                    detail::Sum<double> sum;
                    auto added = add_words(sum, stiffness_action.values, i + j * n);
                    if (!added)
                        return forward_failure(context, added);
                    added = product_words(sum, eigenvalues, j, mass_action.values, i + j * n, -1);
                    if (!added)
                        return forward_failure(context, added);
                    const double propagated =
                        detail::add_magnitudes(
                            {stiffness_action.absolute_error_upper[i + j * n],
                             stiffness_action.absolute_error_upper[i + j * n]},
                            detail::multiply_magnitudes(
                                *lambda, {mass_action.absolute_error_upper[i + j * n],
                                          mass_action.absolute_error_upper[i + j * n]}))
                            .upper;
                    const auto upper = residual_upper(sum, propagated);
                    if (!upper)
                        return forward_failure(context, upper);
                    column_error =
                        std::max(column_error, detail::scale({*upper, *upper}, scales.components[i],
                                                             mass_scale, scales.eigenvalue_scale)
                                                   .upper);
                    if (need_backward) {
                        const auto k_lower = action_lower(stiffness_action, i + j * n);
                        if (!k_lower)
                            return forward_failure(context, k_lower);
                        const auto m_lower = action_lower(mass_action, i + j * n);
                        if (!m_lower)
                            return forward_failure(context, m_lower);
                        const double denominator =
                            detail::add_magnitudes(
                                {*k_lower, *k_lower},
                                detail::multiply_magnitudes(*lambda, {*m_lower, *m_lower}))
                                .lower;
                        backward =
                            std::max(backward, detail::divide_bound(*upper, denominator, true));
                    }
                }
                if (!(infinity_lower > 0))
                    return fail(context, FailureCode::PrecisionUnresolved);
                residual =
                    std::max(residual, detail::divide_bound(column_error, infinity_lower, true));
            }
            for (std::size_t i = 0; i < n; ++i) {
                double row_error = 0;
                for (std::size_t j = 0; j < n; ++j) {
                    detail::Sum<double> sum;
                    double propagated = 0;
                    for (std::size_t inner = 0; inner < n; ++inner) {
                        const auto added = product_words(sum, basis, inner + i * n,
                                                         mass_action.values, inner + j * n, 1);
                        if (!added)
                            return forward_failure(context, added);
                        const auto coordinate = word_magnitude(basis, inner + i * n);
                        if (!coordinate)
                            return forward_failure(context, coordinate);
                        const double error = mass_action.absolute_error_upper[inner + j * n];
                        propagated = detail::add_magnitudes(
                                         {propagated, propagated},
                                         detail::multiply_magnitudes(*coordinate, {error, error}))
                                         .upper;
                    }
                    if (i == j) {
                        const auto diagonal = sum.add(-mass_scale);
                        if (!diagonal)
                            return forward_failure(context, diagonal);
                    }
                    const auto upper = residual_upper(sum, propagated);
                    if (!upper)
                        return forward_failure(context, upper);
                    const double scaled = detail::scale({*upper, *upper}, 1, mass_scale).upper;
                    row_error =
                        detail::add_magnitudes({row_error, row_error}, {scaled, scaled}).upper;
                    if (i == j)
                        normalization = std::max(normalization, scaled);
                }
                gram = std::max(gram, row_error);
            }
            if (!std::isfinite(gram))
                return fail(context, FailureCode::PrecisionUnresolved);
            if (n != 0 && !(gram < 1))
                return fail(context, FailureCode::RankUnresolved);
            const auto precision = representation(basis, eigenvalues);
            const auto residual_proof = certificate(request.residual, residual, precision, context);
            if (!residual_proof)
                return forward_failure(context, residual_proof);
            const auto norm_proof =
                certificate(request.normalization, normalization, precision, context);
            if (!norm_proof)
                return forward_failure(context, norm_proof);
            const auto gram_proof = certificate(request.orthogonality, gram, precision, context);
            if (!gram_proof)
                return forward_failure(context, gram_proof);
            const auto backward_proof =
                backward_certificate(request.backward, backward, precision, context);
            if (!backward_proof)
                return forward_failure(context, backward_proof);
            return GeneralizedMetrics{*residual_proof, *norm_proof, *gram_proof, *backward_proof};
        });
}

template Result<PairMetrics> assess_pair_metrics<double>(WordView<double>, WordView<double>,
                                                         ActionValues<double>, MetricScales,
                                                         PairSide, const PairMetricRequest&,
                                                         const MetricContext&) noexcept;
template Result<PairMetrics> assess_pair_metrics<std::complex<double>>(
    WordView<std::complex<double>>, WordView<std::complex<double>>,
    ActionValues<std::complex<double>>, MetricScales, PairSide, const PairMetricRequest&,
    const MetricContext&) noexcept;
template Result<BasisMetrics> assess_schur_metrics<double>(core::Index, core::Index,
                                                           WordView<double>, ActionValues<double>,
                                                           ActionValues<double>,
                                                           const BasisMetricRequest&,
                                                           const MetricContext&) noexcept;
template Result<BasisMetrics> assess_schur_metrics<std::complex<double>>(
    core::Index, core::Index, WordView<std::complex<double>>, ActionValues<std::complex<double>>,
    ActionValues<std::complex<double>>, const BasisMetricRequest&, const MetricContext&) noexcept;
} // namespace uvlas::numerics::eigen
