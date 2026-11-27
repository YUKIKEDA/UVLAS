#include "compensated.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <uvlas/numerics/precision.hpp>

namespace uvlas::numerics {
namespace {
using core::Failure;
using core::FailureCode;
using core::Index;
using core::Operation;
using core::Result;
constexpr auto operation = Operation::SolveLinear;
auto fail(FailureCode code) noexcept {
    return std::unexpected(Failure{code, operation});
}
bool finite(double value) noexcept {
    return std::isfinite(value);
}
bool finite(std::complex<double> value) noexcept {
    return finite(value.real()) && finite(value.imag());
}
bool same_shape(DenseShape a, DenseShape b) noexcept {
    return a.rows == b.rows && a.columns == b.columns && a.leading_dimension == b.leading_dimension;
}
bool supported_environment() noexcept {
    if (std::fegetround() != FE_TONEAREST) {
        return false;
    }
    // Volatile probes cover flush-to-zero results and denormal-as-zero inputs.
    const volatile double normal = std::numeric_limits<double>::min();
    const volatile double half = 0.5;
    const volatile double subnormal = std::numeric_limits<double>::denorm_min();
    const volatile double one = 1.0;
    const double produced = normal * half;
    const double consumed = subnormal * one;
    return std::bit_cast<std::uint64_t>(produced) == 0x0008000000000000ULL &&
           std::bit_cast<std::uint64_t>(consumed) == 1ULL;
}
int lowest_bit(double value) noexcept {
    // For a nonzero finite value this is the exponent of its least significant 1.
    // Multiplication below that bit cannot be represented exactly even by two doubles.
    const auto bits = std::bit_cast<std::uint64_t>(value) & 0x7fffffffffffffffULL;
    const auto exponent = static_cast<int>(bits >> 52);
    const auto significand = (bits & 0xfffffffffffffULL) | (exponent == 0 ? 0ULL : (1ULL << 52));
    return (exponent == 0 ? -1074 : exponent - 1075) + std::countr_zero(significand);
}
Result<TwoWordValue> sum_finite(double a, double b) noexcept {
    const double sum = a + b;
    if (!finite(sum))
        return fail(FailureCode::Nonfinite);
    const double part_b = sum - a;
    const double error = (a - (sum - part_b)) + (b - part_b);
    if (!finite(error))
        return fail(FailureCode::Nonfinite);
    return TwoWordValue{sum, error};
}
Result<TwoWordValue> product_finite(double a, double b) noexcept {
    const double product = a * b;
    if (!finite(product))
        return fail(FailureCode::Nonfinite);
    if (a != 0 && b != 0 && lowest_bit(a) + lowest_bit(b) < -1074)
        return fail(FailureCode::PrecisionUnresolved);
    return TwoWordValue{product, std::fma(a, b, -product)};
}
using detail::outward_add;
using CompensatedSum = detail::CompensatedAccumulator;
// Compare with a downward enclosure of tolerance*scale without forming its extreme exponent.
bool within_budget(double upper, double scale, double tolerance) noexcept {
    if (!finite(upper)) {
        return false;
    }
    if (upper == 0) {
        return true;
    }
    if (tolerance == 0) {
        return false;
    }
    int se = 0, te = 0, ue = 0, be = 0;
    const double sf = std::frexp(scale, &se), tf = std::frexp(tolerance, &te);
    const auto product = two_product(sf, tf);
    if (!product) {
        return false;
    }
    const double lower = product->lo < 0 ? std::nextafter(product->hi, 0.0) : product->hi;
    const double bf = std::frexp(lower, &be), uf = std::frexp(upper, &ue);
    be += se + te;
    return ue < be || (ue == be && uf <= bf);
}
template <StorageScalar Scalar>
const DenseMatrix<Scalar>& high(const PrecisionBuffer<Scalar>& value) noexcept {
    if (const auto* binary = std::get_if<Binary64<Scalar>>(&value.storage())) {
        return binary->values;
    }
    // Neither alternative can throw on move, so this private variant never becomes valueless.
    return std::get_if<TwoWord<Scalar>>(&value.storage())->hi;
}
template <StorageScalar Scalar>
const DenseMatrix<Scalar>* low(const PrecisionBuffer<Scalar>& value) noexcept {
    const auto* pair = std::get_if<TwoWord<Scalar>>(&value.storage());
    return pair == nullptr ? nullptr : &pair->lo;
}
template <StorageScalar Scalar>
Result<void> residual_product(CompensatedSum& real, CompensatedSum& imag, Scalar a,
                              Scalar x) noexcept {
    if constexpr (std::same_as<Scalar, double>) {
        return real.product(a, x, -1);
    } else {
        const auto rr = real.product(a.real(), x.real(), -1);
        if (!rr) {
            return rr;
        }
        const auto ii = real.product(a.imag(), x.imag(), 1);
        if (!ii) {
            return ii;
        }
        const auto ri = imag.product(a.real(), x.imag(), -1);
        if (!ri) {
            return ri;
        }
        return imag.product(a.imag(), x.real(), -1);
    }
}
} // namespace

Result<TwoWordValue> two_sum(double a, double b) noexcept {
    if (!finite(a) || !finite(b)) {
        return fail(FailureCode::Nonfinite);
    }
    if (!supported_environment()) {
        return fail(FailureCode::PrecisionUnsupported);
    }
    return sum_finite(a, b);
}
Result<TwoWordValue> two_product(double a, double b) noexcept {
    if (!finite(a) || !finite(b)) {
        return fail(FailureCode::Nonfinite);
    }
    if (!supported_environment()) {
        return fail(FailureCode::PrecisionUnsupported);
    }
    return product_finite(a, b);
}

Result<detail::CheckedEft> detail::CheckedEft::prepare() noexcept {
    if (!supported_environment())
        return fail(FailureCode::PrecisionUnsupported);
    return CheckedEft{};
}
Result<TwoWordValue> detail::CheckedEft::sum(double a, double b) const noexcept {
    if (!finite(a) || !finite(b))
        return fail(FailureCode::Nonfinite);
    return sum_finite(a, b);
}
Result<TwoWordValue> detail::CheckedEft::product(double a, double b) const noexcept {
    if (!finite(a) || !finite(b))
        return fail(FailureCode::Nonfinite);
    return product_finite(a, b);
}

template <StorageScalar Scalar>
PrecisionBuffer<Scalar>::PrecisionBuffer(std::variant<Binary64<Scalar>, TwoWord<Scalar>> storage,
                                         core::NumericalGeneration generation) noexcept
    : storage_(std::move(storage)), generation_(generation) {}
template <StorageScalar Scalar>
Result<PrecisionBuffer<Scalar>>
PrecisionBuffer<Scalar>::binary64(DenseMatrix<Scalar> values,
                                  core::NumericalGeneration generation) noexcept {
    const auto valid = values.validate(operation);
    if (!valid) {
        return std::unexpected(valid.error());
    }
    return PrecisionBuffer{Binary64<Scalar>{std::move(values)}, generation};
}
template <StorageScalar Scalar>
Result<PrecisionBuffer<Scalar>>
PrecisionBuffer<Scalar>::two_word(DenseMatrix<Scalar> hi, DenseMatrix<Scalar> lo,
                                  core::NumericalGeneration generation) noexcept {
    const auto valid_hi = hi.validate(operation);
    if (!valid_hi) {
        return std::unexpected(valid_hi.error());
    }
    const auto valid_lo = lo.validate(operation);
    if (!valid_lo) {
        return std::unexpected(valid_lo.error());
    }
    if (!same_shape(hi.shape(), lo.shape())) {
        return fail(FailureCode::InvalidInput);
    }
    return PrecisionBuffer{TwoWord<Scalar>{std::move(hi), std::move(lo)}, generation};
}
template <StorageScalar Scalar>
core::Precision PrecisionBuffer<Scalar>::precision() const noexcept {
    return std::holds_alternative<Binary64<Scalar>>(storage_) ? core::Precision::Binary64
                                                              : core::Precision::TwoWord;
}
template <StorageScalar Scalar>
core::NumericalGeneration PrecisionBuffer<Scalar>::generation() const noexcept {
    return generation_;
}
template <StorageScalar Scalar>
DenseShape PrecisionBuffer<Scalar>::shape() const noexcept {
    return high(*this).shape();
}
template <StorageScalar Scalar>
const std::variant<Binary64<Scalar>, TwoWord<Scalar>>&
PrecisionBuffer<Scalar>::storage() const& noexcept {
    return storage_;
}
Result<void> require_precision(core::Precision actual, SupportedPrecision supported) noexcept {
    if ((actual != core::Precision::Binary64 && actual != core::Precision::TwoWord) ||
        (supported != SupportedPrecision::Binary64Only &&
         supported != SupportedPrecision::Binary64AndTwoWord)) {
        return fail(FailureCode::PrecisionUnsupported);
    }
    if (actual == core::Precision::TwoWord && supported == SupportedPrecision::Binary64Only) {
        return fail(FailureCode::PrecisionUnsupported);
    }
    return {};
}
template <StorageScalar Scalar>
Result<PrecisionBuffer<Scalar>>
add_correction(const PrecisionBuffer<Scalar>& source, const DenseMatrix<Scalar>& correction,
               core::NumericalGeneration generation, core::ResourcePlan& plan) noexcept {
    const auto valid = correction.validate(operation);
    if (!valid) {
        return std::unexpected(valid.error());
    }
    if (!same_shape(source.shape(), correction.shape()) || generation == source.generation()) {
        return fail(FailureCode::InvalidInput);
    }
    auto hi = DenseMatrix<Scalar>::copy(source.shape(), high(source).values(), plan, operation);
    if (!hi) {
        return std::unexpected(hi.error());
    }
    const auto* source_low = low(source);
    auto lo = DenseMatrix<Scalar>::copy(
        source.shape(), source_low == nullptr ? correction.values() : source_low->values(), plan,
        operation);
    if (!lo) {
        return std::unexpected(lo.error());
    }
    if (source_low == nullptr) {
        std::fill(lo->values().begin(), lo->values().end(), Scalar{});
    }
    for (Index column = 0; column < source.shape().columns; ++column) {
        for (Index row = 0; row < source.shape().rows; ++row) {
            const auto i =
                static_cast<std::size_t>(column * source.shape().leading_dimension + row);
            if constexpr (std::same_as<Scalar, double>) {
                const auto value = detail::correct_two_word(hi->values()[i], lo->values()[i],
                                                            correction.values()[i]);
                if (!value) {
                    return std::unexpected(value.error());
                }
                hi->values()[i] = value->hi;
                lo->values()[i] = value->lo;
            } else {
                const auto real = detail::correct_two_word(
                    hi->values()[i].real(), lo->values()[i].real(), correction.values()[i].real());
                if (!real) {
                    return std::unexpected(real.error());
                }
                const auto imag = detail::correct_two_word(
                    hi->values()[i].imag(), lo->values()[i].imag(), correction.values()[i].imag());
                if (!imag) {
                    return std::unexpected(imag.error());
                }
                hi->values()[i] = {real->hi, imag->hi};
                lo->values()[i] = {real->lo, imag->lo};
            }
        }
    }
    return PrecisionBuffer<Scalar>::two_word(std::move(*hi), std::move(*lo), generation);
}
template <StorageScalar Scalar>
Result<DenseMatrix<Scalar>>
evaluate_residual(const DenseMatrix<Scalar>& a, const PrecisionBuffer<Scalar>& x,
                  const DenseMatrix<Scalar>& b, core::ResourcePlan& plan,
                  std::span<double> bounds) noexcept {
    const auto valid_a = a.validate(operation);
    if (!valid_a) {
        return std::unexpected(valid_a.error());
    }
    const auto valid_b = b.validate(operation);
    if (!valid_b) {
        return std::unexpected(valid_b.error());
    }
    const auto as = a.shape(), xs = x.shape(), bs = b.shape();
    if (as.columns != xs.rows || as.rows != bs.rows || xs.columns != bs.columns) {
        return fail(FailureCode::InvalidInput);
    }
    auto output = DenseMatrix<Scalar>::copy(bs, b.values(), plan, operation);
    if (!output) {
        return std::unexpected(output.error());
    }
    const auto* trailing = low(x);
    for (Index column = 0; column < bs.columns; ++column) {
        for (Index row = 0; row < bs.rows; ++row) {
            CompensatedSum real, imag;
            const auto rhs =
                b.values()[static_cast<std::size_t>(column * bs.leading_dimension + row)];
            const auto seed = real.add(std::real(rhs));
            if (!seed) {
                return std::unexpected(seed.error());
            }
            if constexpr (!std::same_as<Scalar, double>) {
                const auto seed_imag = imag.add(rhs.imag());
                if (!seed_imag) {
                    return std::unexpected(seed_imag.error());
                }
            }
            for (int plane = 0; plane < (trailing == nullptr ? 1 : 2); ++plane) {
                const auto values = plane == 0 ? high(x).values() : trailing->values();
                for (Index inner = 0; inner < as.columns; ++inner) {
                    const auto product = residual_product(
                        real, imag,
                        a.values()[static_cast<std::size_t>(inner * as.leading_dimension + row)],
                        values[static_cast<std::size_t>(column * xs.leading_dimension + inner)]);
                    if (!product) {
                        return std::unexpected(product.error());
                    }
                }
            }
            const auto real_result = real.result();
            if (!real_result) {
                return std::unexpected(real_result.error());
            }
            const auto index = static_cast<std::size_t>(column * bs.leading_dimension + row);
            if constexpr (std::same_as<Scalar, double>) {
                output->values()[index] = real_result->value;
                if (!bounds.empty()) {
                    bounds[static_cast<std::size_t>(column * bs.rows + row)] =
                        outward_add(std::abs(real_result->value), real_result->error);
                }
            } else {
                const auto imag_result = imag.result();
                if (!imag_result) {
                    return std::unexpected(imag_result.error());
                }
                output->values()[index] = {real_result->value, imag_result->value};
                if (!bounds.empty()) {
                    const double ru = outward_add(std::abs(real_result->value), real_result->error);
                    const double iu = outward_add(std::abs(imag_result->value), imag_result->error);
                    bounds[static_cast<std::size_t>(column * bs.rows + row)] = outward_add(ru, iu);
                }
            }
        }
    }
    return output;
}
template <StorageScalar Scalar>
Result<DenseMatrix<Scalar>>
original_residual(const DenseMatrix<Scalar>& a, const PrecisionBuffer<Scalar>& x,
                  const DenseMatrix<Scalar>& b, core::ResourcePlan& plan) noexcept {
    return evaluate_residual(a, x, b, plan, {});
}
template <StorageScalar Scalar>
Result<QualifiedProjection<Scalar>>
project_binary64(const PrecisionBuffer<Scalar>& source, const DenseMatrix<Scalar>& a,
                 const DenseMatrix<Scalar>& b, std::span<const double> scales, double tolerance,
                 core::NumericalGeneration target_generation, core::ScopeId consumer,
                 const ProjectionVerifier<Scalar>& verifier, core::ResourcePlan& plan) noexcept {
    return core::guard_entry(
        operation, FailureCode::PrecisionUnresolved, {},
        [&]() -> Result<QualifiedProjection<Scalar>> {
            if (!finite(tolerance)) {
                return fail(FailureCode::Nonfinite);
            }
            const auto bs = b.shape();
            const auto count =
                core::checked_extent_product(std::array<Index, 2>{bs.rows, bs.columns}, operation);
            if (!count) {
                return std::unexpected(count.error());
            }
            if (std::cmp_not_equal(scales.size(), *count) || tolerance < 0 ||
                target_generation == source.generation()) {
                return fail(FailureCode::InvalidInput);
            }
            for (const double scale : scales) {
                if (!finite(scale)) {
                    return fail(FailureCode::Nonfinite);
                }
                if (scale <= 0) {
                    return fail(FailureCode::InvalidInput);
                }
            }
            // Reject malformed residual inputs before any projection storage is charged.
            const auto valid_a = a.validate(operation);
            if (!valid_a) {
                return std::unexpected(valid_a.error());
            }
            const auto valid_b = b.validate(operation);
            if (!valid_b) {
                return std::unexpected(valid_b.error());
            }
            const auto as = a.shape(), xs = source.shape();
            if (as.columns != xs.rows || as.rows != bs.rows || xs.columns != bs.columns) {
                return fail(FailureCode::InvalidInput);
            }
            auto values =
                DenseMatrix<Scalar>::copy(source.shape(), high(source).values(), plan, operation);
            if (!values) {
                return std::unexpected(values.error());
            }
            const auto* trailing = low(source);
            if (trailing != nullptr) {
                for (Index column = 0; column < source.shape().columns; ++column) {
                    for (Index row = 0; row < source.shape().rows; ++row) {
                        const auto i = static_cast<std::size_t>(
                            column * source.shape().leading_dimension + row);
                        values->values()[i] += trailing->values()[i];
                        if (!finite(values->values()[i])) {
                            return fail(FailureCode::Nonfinite);
                        }
                    }
                }
            }
            auto projected =
                PrecisionBuffer<Scalar>::binary64(std::move(*values), target_generation);
            if (!projected) {
                return std::unexpected(projected.error());
            }
            auto bounds = plan.allocate<double>(*count, operation);
            if (!bounds) {
                return std::unexpected(bounds.error());
            }
            const auto residual = evaluate_residual(a, *projected, b, plan, bounds->values());
            if (!residual) {
                return std::unexpected(residual.error());
            }
            for (Index column = 0; column < bs.columns; ++column) {
                for (Index row = 0; row < bs.rows; ++row) {
                    const auto scale_index = static_cast<std::size_t>(column * bs.rows + row);
                    if (!within_budget(bounds->values()[scale_index], scales[scale_index],
                                       tolerance)) {
                        return fail(FailureCode::PrecisionUnresolved);
                    }
                }
            }
            const ProjectionContext context{source.generation(), target_generation, consumer};
            const auto checked = verifier.verify(source, high(*projected), context);
            if (!checked) {
                return std::unexpected(checked.error());
            }
            if (*checked != ProjectionAssessment::Verified) {
                return fail(FailureCode::PrecisionUnresolved);
            }
            // Only the exact immutable candidate evaluated above can acquire this qualification.
            auto& owned = std::get<Binary64<Scalar>>(projected->storage_).values;
            return QualifiedProjection<Scalar>{std::move(owned), context};
        });
}

template class PrecisionBuffer<double>;
template class PrecisionBuffer<std::complex<double>>;
template Result<PrecisionBuffer<double>> add_correction(const PrecisionBuffer<double>&,
                                                        const DenseReal&, core::NumericalGeneration,
                                                        core::ResourcePlan&) noexcept;
template Result<PrecisionBuffer<std::complex<double>>>
add_correction(const PrecisionBuffer<std::complex<double>>&, const DenseComplex&,
               core::NumericalGeneration, core::ResourcePlan&) noexcept;
template Result<DenseReal> original_residual(const DenseReal&, const PrecisionBuffer<double>&,
                                             const DenseReal&, core::ResourcePlan&) noexcept;
template Result<DenseComplex> original_residual(const DenseComplex&,
                                                const PrecisionBuffer<std::complex<double>>&,
                                                const DenseComplex&, core::ResourcePlan&) noexcept;
template Result<QualifiedProjection<double>>
project_binary64(const PrecisionBuffer<double>&, const DenseReal&, const DenseReal&,
                 std::span<const double>, double, core::NumericalGeneration, core::ScopeId,
                 const ProjectionVerifier<double>&, core::ResourcePlan&) noexcept;
template Result<QualifiedProjection<std::complex<double>>>
project_binary64(const PrecisionBuffer<std::complex<double>>&, const DenseComplex&,
                 const DenseComplex&, std::span<const double>, double, core::NumericalGeneration,
                 core::ScopeId, const ProjectionVerifier<std::complex<double>>&,
                 core::ResourcePlan&) noexcept;
} // namespace uvlas::numerics
