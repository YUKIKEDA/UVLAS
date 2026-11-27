#include "pg_reference_projection.hpp"

#include <algorithm>
#include <cmath>
#include <uvlas/structure/pg_shell_material.hpp>

namespace uvlas::structure {
namespace {
using detail::partition_entry;
using detail::partition_fail;
using detail::partition_take;
using detail::ReferenceArithmetic;
using detail::ReferenceRadicals;
using R = ExactRational;
constexpr std::size_t order = 6;
using Matrix = std::array<R, order * order>;
constexpr auto invalid = core::FailureCode::InvalidInput;
constexpr auto unresolved = core::FailureCode::PrecisionUnresolved;

bool engineering(PgMaterialComponents tag) {
    switch (tag) {
    case PgMaterialComponents::Engineering:
        return true;
    case PgMaterialComponents::Mandel:
        return false;
    }
    partition_fail(invalid);
}

const R& entry(const Matrix& matrix, std::size_t i, std::size_t j) {
    return matrix[i * order + j];
}

// Exact LDL pivots. Scaling only the diagonal by 1-tau tests
// D^-1 C D^-1 - tau I by congruence, with D^2 = diag(C).
bool positive(const Matrix& value, const R& diagonal_scale, ReferenceArithmetic& math) {
    Matrix remaining;
    for (std::size_t i = 0; i < order; ++i) {
        for (std::size_t j = 0; j < i; ++j)
            remaining[i * order + j] = math.copy(entry(value, i, j));
        remaining[i * order + i] = math.mul(entry(value, i, i), diagonal_scale);
    }
    for (std::size_t k = 0; k < order; ++k) {
        const auto& pivot = entry(remaining, k, k);
        if (pivot.sign() <= 0)
            return false;
        for (std::size_t i = k + 1; i < order; ++i) {
            const auto multiplier = math.div(entry(remaining, i, k), pivot);
            for (std::size_t j = k + 1; j <= i; ++j) {
                auto& cell = remaining[i * order + j];
                cell = math.sub(cell, math.mul(multiplier, entry(remaining, j, k)));
            }
        }
    }
    return true;
}

void qualify(const Matrix& value, const R& one, const R& reduced_diagonal,
             ReferenceArithmetic& math) {
    if (!positive(value, one, math))
        partition_fail(invalid);
    if (!positive(value, reduced_diagonal, math))
        partition_fail(unresolved);
}

Matrix isotropic(const PgIsotropicElasticity& value, const R& one, const R& two,
                 ReferenceArithmetic& math) {
    const auto plus = math.add(one, value.poisson);
    const auto minus = math.sub(one, math.mul(two, value.poisson));
    if (value.young_pascal.sign() <= 0 || plus.sign() <= 0 || minus.sign() <= 0)
        partition_fail(invalid);
    const auto two_mu = math.div(value.young_pascal, plus);
    const auto lambda = math.div(math.mul(two_mu, value.poisson), minus);
    Matrix result;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j)
            result[i * order + j] = i == j ? math.add(lambda, two_mu) : math.copy(lambda);
        result[(i + 3) * order + i + 3] = math.copy(two_mu);
    }
    return result;
}

Matrix orthotropic(const PgOrthotropicElasticity& value, const R& one, const R& two,
                   const R& reduced_diagonal, ReferenceArithmetic& math) {
    Matrix compliance;
    for (std::size_t i = 0; i < 3; ++i) {
        if (value.young_pascal[i].sign() <= 0 || value.shear_pascal[i].sign() <= 0)
            partition_fail(invalid);
        compliance[i * order + i] = math.div(one, value.young_pascal[i]);
        compliance[(i + 3) * order + i + 3] = math.div(one, math.mul(two, value.shear_pascal[i]));
    }
    constexpr std::array<std::array<std::size_t, 2>, 3> pair{{{0, 1}, {0, 2}, {1, 2}}};
    for (std::size_t k = 0; k < pair.size(); ++k) {
        const auto [i, j] = pair[k];
        compliance[i * order + j] = math.neg(math.div(value.poisson[k], value.young_pascal[i]));
        compliance[j * order + i] = math.copy(compliance[i * order + j]);
    }
    qualify(compliance, one, reduced_diagonal, math);
    // The normal block is the sole coupled block of the Mandel compliance.
    Matrix adjugate;
    for (std::size_t i = 0; i < 3; ++i) {
        const auto i1 = (i + 1) % 3, i2 = (i + 2) % 3;
        for (std::size_t j = 0; j < 3; ++j) {
            const auto j1 = (j + 1) % 3, j2 = (j + 2) % 3;
            adjugate[i * order + j] =
                math.sub(math.mul(entry(compliance, i1, j1), entry(compliance, i2, j2)),
                         math.mul(entry(compliance, i1, j2), entry(compliance, i2, j1)));
        }
    }
    R determinant;
    for (std::size_t j = 0; j < 3; ++j)
        determinant =
            math.add(determinant, math.mul(entry(compliance, 0, j), entry(adjugate, j, 0)));
    Matrix result;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j)
            result[i * order + j] = math.div(entry(adjugate, i, j), determinant);
        result[(i + 3) * order + i + 3] = math.mul(two, value.shear_pascal[i]);
    }
    return result;
}

struct SymmetricLaw {
    Matrix coefficients;
    R correction_squared;
};

SymmetricLaw effective(const PgEffectiveElasticity& input, const R& one, const R& two,
                       const R& tolerance, ReferenceArithmetic& math) {
    const bool input_engineering = engineering(input.components);
    const auto four = math.mul(two, two);
    const auto& raw = input.stiffness_pascal;
    R norm_squared, asymmetry_squared;
    for (std::size_t i = 0; i < order; ++i) {
        for (std::size_t j = 0; j < order; ++j) {
            const auto count = static_cast<unsigned>(i >= 3) + static_cast<unsigned>(j >= 3);
            const auto& scale_squared =
                !input_engineering || count == 0 ? one : (count == 1 ? two : four);
            norm_squared =
                math.add(norm_squared, math.mul(math.square(entry(raw, i, j)), scale_squared));
            const auto delta = math.sub(entry(raw, i, j), entry(raw, j, i));
            asymmetry_squared =
                math.add(asymmetry_squared, math.mul(math.square(delta), scale_squared));
        }
    }
    if (norm_squared.sign() == 0 ||
        math.compare(asymmetry_squared, math.mul(math.square(tolerance), norm_squared)) > 0)
        partition_fail(invalid);
    Matrix symmetric;
    for (std::size_t i = 0; i < order; ++i) {
        for (std::size_t j = 0; j < order; ++j) {
            symmetric[i * order + j] = math.div(math.add(entry(raw, i, j), entry(raw, j, i)), two);
        }
    }
    return {std::move(symmetric), math.div(asymmetry_squared, math.mul(four, norm_squared))};
}

void qualify_natural_metric(const std::array<PgReferenceScalar, 6>& strain, const R& one,
                            const R& two, ReferenceArithmetic& math) {
    ReferenceRadicals exact{math, two};
    const auto diagonal = [&](std::size_t i) {
        return exact.add(PgReferenceScalar{math.copy(one), {}}, exact.scale(strain[i], two));
    };
    const auto off_diagonal = [&](std::size_t i) {
        return PgReferenceScalar{math.mul(two, strain[i].radical), math.copy(strain[i].rational)};
    };
    const auto a = diagonal(0), b = diagonal(1), c = diagonal(2);
    const auto yz = off_diagonal(3), xz = off_diagonal(4), xy = off_diagonal(5);
    const auto leading = exact.sub(exact.mul(a, b), exact.mul(xy, xy));
    if (exact.sign(a) <= 0 || exact.sign(leading) <= 0)
        partition_fail(invalid);
    const auto determinant = exact.sub(
        exact.add(exact.mul(leading, c), exact.scale(exact.mul(xy, exact.mul(xz, yz)), two)),
        exact.add(exact.mul(a, exact.mul(yz, yz)), exact.mul(b, exact.mul(xz, xz))));
    if (exact.sign(determinant) <= 0)
        partition_fail(invalid);
}

using Mandel = std::array<PgReferenceScalar, 6>;

Mandel normalize_initial(const PgInitialGreenStrain& initial, const R& two,
                         ReferenceArithmetic& math) {
    const bool initial_engineering = engineering(initial.convention);
    Mandel strain;
    for (std::size_t i = 0; i < order; ++i) {
        if (initial_engineering && i >= 3)
            strain[i].radical = math.div(initial.components[i], two);
        else
            strain[i].rational = math.copy(initial.components[i]);
    }
    qualify_natural_metric(strain, math.constant("1"), two, math);
    return strain;
}

PgShellMaterialResponse response(const PgShellElasticLaw& law, const Mandel& initial,
                                 const Mandel& green, const R& two, ReferenceArithmetic& math) {
    if (law.density_kg_m3().sign() <= 0)
        partition_fail(invalid);
    ReferenceRadicals exact{math, two};
    Mandel elastic, stress;
    for (std::size_t i = 0; i < order; ++i)
        elastic[i] = exact.sub(green[i], initial[i]);
    PgReferenceScalar doubled_energy;
    for (std::size_t i = 0; i < order; ++i) {
        for (std::size_t j = 0; j < order; ++j)
            stress[i] = exact.add(stress[i], exact.mul(law.stiffness()[i * order + j], elastic[j]));
        doubled_energy = exact.add(doubled_energy, exact.mul(elastic[i], stress[i]));
    }
    PgReferenceScalar energy{math.div(doubled_energy.rational, two),
                             math.div(doubled_energy.radical, two)};
    return {std::move(energy), std::move(stress), math.info()};
}

template <class F>
auto evaluate(F&& f) noexcept {
    auto result = partition_entry(std::forward<F>(f));
    if (!result)
        result.error().operation = core::Operation::EvaluateStructure;
    return result;
}

using Bound = PgMaterialEnclosure;

struct MaterialBounds {
    ReferenceArithmetic& math;

    Bound add(const Bound& a, const Bound& b) const {
        return {math.add(a.lower, b.lower), math.add(a.upper, b.upper)};
    }
    Bound sub(const Bound& a, const Bound& b) const {
        return {math.sub(a.lower, b.upper), math.sub(a.upper, b.lower)};
    }
    Bound scale(const Bound& a, const R& factor) const {
        if (factor.sign() < 0)
            return {math.mul(a.upper, factor), math.mul(a.lower, factor)};
        return {math.mul(a.lower, factor), math.mul(a.upper, factor)};
    }
    Bound multiply(const Bound& a, const Bound& b) const {
        if ((a.lower.sign() == 0 && a.upper.sign() == 0) ||
            (b.lower.sign() == 0 && b.upper.sign() == 0))
            return {};
        std::array products{math.mul(a.lower, b.lower), math.mul(a.lower, b.upper),
                            math.mul(a.upper, b.lower), math.mul(a.upper, b.upper)};
        const auto [low, high] = std::ranges::minmax_element(
            products, [&](const R& x, const R& y) { return math.compare(x, y) < 0; });
        return {math.copy(*low), math.copy(*high)};
    }
    Bound square(const Bound& a) const {
        auto left = math.square(a.lower), right = math.square(a.upper);
        if (math.compare(left, right) > 0)
            std::swap(left, right);
        if (a.lower.sign() <= 0 && a.upper.sign() >= 0)
            left = {};
        return {std::move(left), std::move(right)};
    }
    Bound divide_positive(const Bound& a, const Bound& b, const R& one) const {
        return multiply(a, {math.div(one, b.upper), math.div(one, b.lower)});
    }
    Bound root(const R& g) const {
        if (g.sign() <= 0)
            partition_fail(invalid);
        const R zero;
        const detail::PgProjectionComparator squared{math, zero, 0};
        const detail::PgPositiveRootComparator comparator{math, squared, g};
        const auto high = detail::reference_binary64(math, detail::pg_round_word(comparator, {}));
        if (comparator.compare({}, high) == 0)
            return {math.copy(high), math.copy(high)};
        const auto negative_high = math.neg(high);
        const auto low = detail::pg_round_word(comparator, negative_high);
        const auto low_exact = detail::reference_binary64(math, low);
        if (comparator.compare(negative_high, low_exact) == 0) {
            auto root_value = math.add(high, low_exact);
            return {math.copy(root_value), std::move(root_value)};
        }
        // The two-word residual is rounded to nearest. Its adjacent words
        // therefore enclose the exact residual even at subnormal boundaries.
        const auto inf = std::numeric_limits<double>::infinity();
        return {math.add(high, detail::reference_binary64(math, std::nextafter(low, -inf))),
                math.add(high, detail::reference_binary64(math, std::nextafter(low, inf)))};
    }
    Bound radical(const PgReferenceScalar& value, const Bound& root) const {
        if (value.radical.sign() == 0)
            return {math.copy(value.rational), math.copy(value.rational)};
        auto term = scale(root, value.radical);
        return {math.add(value.rational, term.lower), math.add(value.rational, term.upper)};
    }
};

void require_positive(const Bound& value, core::FailureCode nonpositive = invalid) {
    if (value.upper.sign() <= 0)
        partition_fail(nonpositive);
    if (value.lower.sign() <= 0)
        partition_fail(unresolved);
}

void qualify_natural_metric(const std::array<Bound, order>& initial, const Bound& sqrt_two,
                            const R& two, const MaterialBounds& bounds) {
    const auto one = bounds.math.constant("1");
    const Bound unit{bounds.math.copy(one), bounds.math.copy(one)};
    const auto a = bounds.add(unit, bounds.scale(initial[0], two));
    const auto b = bounds.add(unit, bounds.scale(initial[1], two));
    const auto c = bounds.add(unit, bounds.scale(initial[2], two));
    const auto yz = bounds.multiply(sqrt_two, initial[3]);
    const auto xz = bounds.multiply(sqrt_two, initial[4]);
    const auto xy = bounds.multiply(sqrt_two, initial[5]);
    require_positive(a);
    const auto leading = bounds.sub(bounds.multiply(a, b), bounds.square(xy));
    require_positive(leading);
    const auto triple = bounds.scale(bounds.multiply(xy, bounds.multiply(xz, yz)), two);
    const auto cross =
        bounds.add(bounds.multiply(a, bounds.square(yz)), bounds.multiply(b, bounds.square(xz)));
    require_positive(bounds.sub(bounds.add(bounds.multiply(leading, c), triple), cross));
}
PgShellMaterialEnclosure enclose_response(const PgShellElasticLaw& law,
                                          const PgInitialGreenEnclosure& input,
                                          const std::array<Bound, order>& green, const R& two,
                                          const Bound& sqrt_two, const MaterialBounds& bounds) {
    auto& math = bounds.math;
    if (law.density_kg_m3().sign() <= 0 || input.radicand.sign() <= 0)
        partition_fail(invalid);
    const auto distinct_root =
        input.radicand != two && std::ranges::any_of(input.mandel, [](const auto& value) {
            return value.radical.sign() != 0;
        });
    const auto separate_root = distinct_root ? bounds.root(input.radicand) : Bound{};
    const auto& field_root = distinct_root ? separate_root : sqrt_two;
    std::array<Bound, order> initial, elastic, stress;
    for (std::size_t i = 0; i < order; ++i) {
        const auto& error = input.absolute_error[i];
        if (error.sign() < 0 || math.compare(green[i].lower, green[i].upper) > 0)
            partition_fail(invalid);
        auto centre = bounds.radical(input.mandel[i], field_root);
        initial[i] = {math.sub(centre.lower, error), math.add(centre.upper, error)};
        elastic[i] = bounds.sub(green[i], initial[i]);
    }
    qualify_natural_metric(initial, sqrt_two, two, bounds);
    Bound doubled_energy;
    for (std::size_t i = 0; i < order; ++i) {
        for (std::size_t j = 0; j < order; ++j) {
            const auto stiffness = bounds.radical(law.stiffness()[i * order + j], sqrt_two);
            const auto linear = bounds.multiply(stiffness, elastic[j]);
            stress[i] = bounds.add(stress[i], linear);
            if (i == j) {
                doubled_energy = bounds.add(doubled_energy,
                                            bounds.multiply(stiffness, bounds.square(elastic[i])));
            } else if (i < j) {
                const auto cross = bounds.multiply(elastic[i], linear);
                doubled_energy = bounds.add(doubled_energy, bounds.scale(cross, two));
            }
        }
    }
    if (doubled_energy.lower.sign() < 0)
        doubled_energy.lower = {}; // Positive definite law bounds the energy from below.
    Bound energy{math.div(doubled_energy.lower, two), math.div(doubled_energy.upper, two)};
    return PgShellMaterialEnclosure{std::move(energy), std::move(stress), math.info()};
}

using Tensor = PgShellTensorEnclosure;
using Vector = std::array<Bound, 3>;

Bound point(const R& value, ReferenceArithmetic& math) {
    return {math.copy(value), math.copy(value)};
}

Vector column(const Tensor& a, std::size_t j, ReferenceArithmetic& math) {
    return {Bound{math.copy(a[j].lower), math.copy(a[j].upper)},
            Bound{math.copy(a[3 + j].lower), math.copy(a[3 + j].upper)},
            Bound{math.copy(a[6 + j].lower), math.copy(a[6 + j].upper)}};
}

Bound dot(const Vector& a, const Vector& b, const MaterialBounds& bounds) {
    Bound result;
    for (std::size_t k = 0; k < a.size(); ++k)
        result = bounds.add(result, bounds.multiply(a[k], b[k]));
    return result;
}

Vector cross(const Vector& a, const Vector& b, const MaterialBounds& bounds) {
    return {bounds.sub(bounds.multiply(a[1], b[2]), bounds.multiply(a[2], b[1])),
            bounds.sub(bounds.multiply(a[2], b[0]), bounds.multiply(a[0], b[2])),
            bounds.sub(bounds.multiply(a[0], b[1]), bounds.multiply(a[1], b[0]))};
}

Bound determinant(const Tensor& a, const MaterialBounds& bounds) {
    const auto first = column(a, 0, bounds.math);
    const auto second = column(a, 1, bounds.math);
    const auto third = column(a, 2, bounds.math);
    return dot(first, cross(second, third, bounds), bounds);
}

Vector normalize(Vector a, const R& one, const MaterialBounds& bounds) {
    Bound norm_squared;
    for (const auto& term : a)
        norm_squared = bounds.add(norm_squared, bounds.square(term));
    require_positive(norm_squared, core::FailureCode::GeometryInvalid);
    const auto lower = bounds.root(norm_squared.lower);
    auto upper = norm_squared.lower == norm_squared.upper ? bounds.math.copy(lower.upper)
                                                          : bounds.root(norm_squared.upper).upper;
    const Bound norm{bounds.math.copy(lower.lower), std::move(upper)};
    require_positive(norm);
    for (auto& value : a)
        value = bounds.divide_positive(value, norm, one);
    return a;
}

Tensor material_axes(const PgShellExplicitAxes& input, const R& one, const MaterialBounds& bounds) {
    const R zero;
    Tensor axes;
    for (std::size_t j = 0; j < axes.size(); ++j)
        axes[j] = point(input.columns[j], bounds.math);
    for (std::size_t i = 0; i < 3; ++i) {
        const auto a = column(axes, i, bounds.math);
        for (std::size_t j = i; j < 3; ++j) {
            const auto entry = dot(a, column(axes, j, bounds.math), bounds);
            if (entry.lower != (i == j ? one : zero))
                partition_fail(core::FailureCode::GeometryInvalid);
        }
    }
    require_positive(determinant(axes, bounds), core::FailureCode::GeometryInvalid);
    return axes;
}

Tensor material_axes(const PgShellProjectedAxes& input, const R& one,
                     const MaterialBounds& bounds) {
    auto& math = bounds.math;
    const auto& [world_fibre, cosine, sine] = input.fibre;
    if (math.add(math.square(cosine), math.square(sine)) != one)
        partition_fail(core::FailureCode::GeometryInvalid);
    const auto& ref = input.reference;
    const auto root = bounds.root(ref.g);
    ReferenceRadicals exact{math, ref.g};
    PgReferenceVector fibre;
    for (std::size_t i = 0; i < 3; ++i)
        fibre[i].rational = math.copy(world_fibre[i]);
    const auto area = exact.cross(ref.jacobian_columns_metres[0], ref.jacobian_columns_metres[1]);
    const auto oriented_volume = exact.dot(area, ref.jacobian_columns_metres[2]);
    if (exact.sign(oriented_volume) <= 0)
        partition_fail(core::FailureCode::GeometryInvalid);
    // n x (w x n) is |n|^2 times the projected fibre. Form it before
    // enclosing either direction so exact parallelism stays distinguishable.
    const auto projected = exact.cross(area, exact.cross(fibre, area));
    const auto enclose = [&](const PgReferenceVector& value) {
        Vector result;
        for (std::size_t i = 0; i < result.size(); ++i) {
            if (exact.sign(value[i]) != 0)
                result[i] = bounds.radical(value[i], root);
        }
        return result;
    };
    const auto normal = normalize(enclose(area), one, bounds);
    const auto along = normalize(enclose(projected), one, bounds);
    const auto across = cross(normal, along, bounds);
    Tensor axes;
    for (std::size_t i = 0; i < 3; ++i) {
        axes[3 * i] = bounds.add(bounds.scale(along[i], cosine), bounds.scale(across[i], sine));
        axes[3 * i + 1] = bounds.sub(bounds.scale(across[i], cosine), bounds.scale(along[i], sine));
        axes[3 * i + 2] = {math.copy(normal[i].lower), math.copy(normal[i].upper)};
    }
    return axes;
}

enum class RightTensor { AsStored, Transposed };

Tensor multiply(const Tensor& a, const Tensor& b, RightTensor right, const MaterialBounds& bounds) {
    Tensor result;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            auto& x = result[3 * i + j];
            for (std::size_t k = 0; k < 3; ++k) {
                const auto index = right == RightTensor::Transposed ? 3 * j + k : 3 * k + j;
                x = bounds.add(x, bounds.multiply(a[3 * i + k], b[index]));
            }
        }
    }
    return result;
}

constexpr std::array mandel_axes{std::array{0U, 0U}, std::array{1U, 1U}, std::array{2U, 2U},
                                 std::array{1U, 2U}, std::array{0U, 2U}, std::array{0U, 1U}};

std::array<Bound, order> green_material(const Tensor& fa, const R& one, const R& half,
                                        const Bound& sqrt_two, const MaterialBounds& bounds) {
    const std::array axes{column(fa, 0, bounds.math), column(fa, 1, bounds.math),
                          column(fa, 2, bounds.math)};
    const auto unit = point(one, bounds.math);
    const auto scale_shear = bounds.scale(sqrt_two, half);
    std::array<Bound, order> green;
    for (std::size_t k = 0; k < order; ++k) {
        const auto [i, j] = mandel_axes[k];
        if (i == j) {
            Bound squared;
            for (const auto& x : axes[i])
                squared = bounds.add(squared, bounds.square(x));
            green[k] = bounds.scale(bounds.sub(squared, unit), half);
        } else {
            green[k] = bounds.multiply(dot(axes[i], axes[j], bounds), scale_shear);
        }
    }
    return green;
}

Tensor stress_tensor(const std::array<Bound, order>& stress, const Bound& sqrt_two, const R& half,
                     const MaterialBounds& bounds) {
    Tensor tensor;
    const auto inverse_sqrt_two = bounds.scale(sqrt_two, half);
    for (std::size_t k = 0; k < order; ++k) {
        const auto [i, j] = mandel_axes[k];
        if (i == j) {
            tensor[3 * i + i] = {bounds.math.copy(stress[k].lower),
                                 bounds.math.copy(stress[k].upper)};
        } else {
            auto shear = bounds.multiply(stress[k], inverse_sqrt_two);
            tensor[3 * j + i] = {bounds.math.copy(shear.lower), bounds.math.copy(shear.upper)};
            tensor[3 * i + j] = std::move(shear);
        }
    }
    return tensor;
}

void add_gradient_column(std::array<PgReferenceScalar, 9>& gradient, std::size_t row,
                         const PgReferenceVector& shape, const R& coefficient,
                         ReferenceRadicals& exact) {
    for (std::size_t column = 0; column < shape.size(); ++column) {
        auto& entry = gradient[3 * row + column];
        entry = exact.add(entry, exact.scale(shape[column], coefficient));
    }
}

Tensor first_pk_derivative(const PgShellElasticLaw& law, const PgShellMaterialPointResponse& value,
                           const Tensor& f, const Tensor& df, const MaterialBounds& bounds) {
    auto& math = bounds.math;
    const auto one = math.constant("1"), two = math.constant("2");
    const auto half = math.div(one, two);
    const auto sqrt_two = bounds.root(two);
    const auto& axes = value.reference_axes;
    const auto fa = multiply(f, axes, RightTensor::AsStored, bounds);
    const auto dfa = multiply(df, axes, RightTensor::AsStored, bounds);
    const std::array base_columns{column(fa, 0, math), column(fa, 1, math), column(fa, 2, math)};
    const std::array delta_columns{column(dfa, 0, math), column(dfa, 1, math),
                                   column(dfa, 2, math)};
    const auto mandel_shear = bounds.scale(sqrt_two, half);
    std::array<Bound, order> green, stress;
    for (std::size_t k = 0; k < order; ++k) {
        const auto [i, j] = mandel_axes[k];
        auto derivative = dot(base_columns[i], delta_columns[j], bounds);
        if (i != j)
            derivative = bounds.multiply(
                bounds.add(derivative, dot(delta_columns[i], base_columns[j], bounds)),
                mandel_shear);
        green[k] = std::move(derivative);
    }
    for (std::size_t i = 0; i < order; ++i) {
        for (std::size_t j = 0; j < order; ++j) {
            const auto stiffness = bounds.radical(law.stiffness()[i * order + j], sqrt_two);
            stress[i] = bounds.add(stress[i], bounds.multiply(stiffness, green[j]));
        }
    }
    const auto local_stress = stress_tensor(stress, sqrt_two, half, bounds);
    const auto world_stress = multiply(multiply(axes, local_stress, RightTensor::AsStored, bounds),
                                       axes, RightTensor::Transposed, bounds);
    auto derivative = multiply(f, world_stress, RightTensor::AsStored, bounds);
    const auto geometric =
        multiply(df, value.second_pk_world_pascal, RightTensor::AsStored, bounds);
    for (std::size_t k = 0; k < derivative.size(); ++k)
        derivative[k] = bounds.add(derivative[k], geometric[k]);
    return derivative;
}
} // namespace

PgShellElasticLaw::PgShellElasticLaw(std::array<PgReferenceScalar, 36>&& stiffness, R&& density,
                                     R&& correction_squared, PgReferenceCallInfo work) noexcept
    : stiffness_(std::move(stiffness)),
      density_(std::move(density)),
      correction_squared_(std::move(correction_squared)),
      work_(work) {}

core::Result<PgShellElasticLaw> PgShellElasticLaw::prepare(const PgElasticity& input,
                                                           const R& density,
                                                           const PgReferenceLimits& limits,
                                                           core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgShellElasticLaw> {
        ReferenceArithmetic math{limits, plan};
        math.cancel();
        if (density.sign() <= 0)
            partition_fail(invalid);
        // Covers borrowed inputs, compliance/adjugate, normalization and LDL copies,
        // and the two coefficients of every normalized value.
        math.admit(384);
        const auto one = math.constant("1"), two = math.constant("2");
        math.charge();
        const auto tolerance =
            partition_take(R::from_decimal("1e-12", math.leaf_limits, math.plan));
        const auto reduced_diagonal = math.sub(one, tolerance);
        auto law = std::visit(
            [&](const auto& value) -> SymmetricLaw {
                using Kind = std::decay_t<decltype(value)>;
                if constexpr (std::is_same_v<Kind, PgIsotropicElasticity>) {
                    return {isotropic(value, one, two, math), {}};
                } else if constexpr (std::is_same_v<Kind, PgOrthotropicElasticity>) {
                    return {orthotropic(value, one, two, reduced_diagonal, math), {}};
                } else {
                    return effective(value, one, two, tolerance, math);
                }
            },
            input);
        // Positive diagonal congruence to Mandel leaves the scaled eigenvalues unchanged.
        qualify(law.coefficients, one, reduced_diagonal, math);
        const auto* full = std::get_if<PgEffectiveElasticity>(&input);
        const bool engineering_stiffness = full && engineering(full->components);
        std::array<PgReferenceScalar, order * order> stiffness;
        for (std::size_t i = 0; i < order; ++i) {
            for (std::size_t j = 0; j < order; ++j) {
                auto& input_value = law.coefficients[i * order + j];
                auto& output = stiffness[i * order + j];
                if (!engineering_stiffness || (i < 3 && j < 3))
                    output.rational = std::move(input_value);
                else if (i >= 3 && j >= 3)
                    output.rational = math.mul(two, input_value);
                else
                    output.radical = std::move(input_value);
            }
        }
        auto rho = math.copy(density);
        return PgShellElasticLaw{std::move(stiffness), std::move(rho),
                                 std::move(law.correction_squared), math.info()};
    });
}

PgShellMaterial::PgShellMaterial(PgShellElasticLaw&& law, Mandel&& initial,
                                 PgReferenceCallInfo work) noexcept
    : law_(std::move(law)), initial_(std::move(initial)), work_(work) {}

core::Result<PgShellMaterial> PgShellMaterial::prepare(const PgElasticity& input, const R& density,
                                                       const PgInitialGreenStrain& initial,
                                                       const PgReferenceLimits& limits,
                                                       core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgShellMaterial> {
        ReferenceArithmetic math{limits, plan};
        auto law = partition_take(PgShellElasticLaw::prepare(input, density, limits, plan));
        math.charge(law.call_info().scalar_operations);
        math.admit();
        auto strain = normalize_initial(initial, math.constant("2"), math);
        return PgShellMaterial{std::move(law), std::move(strain), math.info()};
    });
}

core::Result<PgShellMaterialResponse>
evaluate_pg_shell_material(const PgShellMaterial& material, const Mandel& green,
                           const PgReferenceLimits& limits, core::ResourcePlan& plan) noexcept {
    return evaluate([&]() -> core::Result<PgShellMaterialResponse> {
        ReferenceArithmetic math{limits, plan};
        math.admit();
        return response(material.law(), material.initial_strain(), green, math.constant("2"), math);
    });
}

core::Result<PgShellMaterialResponse>
evaluate_pg_shell_material(const PgShellElasticLaw& law, const PgInitialGreenStrain& initial,
                           const Mandel& green, const PgReferenceLimits& limits,
                           core::ResourcePlan& plan) noexcept {
    return evaluate([&]() -> core::Result<PgShellMaterialResponse> {
        ReferenceArithmetic math{limits, plan};
        math.admit();
        const auto two = math.constant("2");
        const auto strain = normalize_initial(initial, two, math);
        return response(law, strain, green, two, math);
    });
}

core::Result<PgShellMaterialEnclosure>
enclose_pg_shell_material(const PgShellElasticLaw& law, const PgInitialGreenEnclosure& input,
                          const std::array<PgMaterialEnclosure, 6>& green,
                          const PgReferenceLimits& limits, core::ResourcePlan& plan) noexcept {
    return evaluate([&]() -> core::Result<PgShellMaterialEnclosure> {
        ReferenceArithmetic math{limits, plan};
        math.admit();
        const MaterialBounds bounds{math};
        const auto two = math.constant("2");
        return enclose_response(law, input, green, two, bounds.root(two), bounds);
    });
}

core::Result<PgShellMaterialPointResponse>
enclose_pg_shell_material_point(const PgShellElasticLaw& law, const PgInitialGreenEnclosure& input,
                                const PgShellReferenceAxes& axes_source, const Tensor& f,
                                const PgReferenceLimits& limits,
                                core::ResourcePlan& plan) noexcept {
    return evaluate([&]() -> core::Result<PgShellMaterialPointResponse> {
        ReferenceArithmetic math{limits, plan};
        math.admit(384);
        const MaterialBounds bounds{math};
        for (const auto& entry : f) {
            if (math.compare(entry.lower, entry.upper) > 0)
                partition_fail(invalid);
        }
        auto det = determinant(f, bounds);
        require_positive(det, core::FailureCode::GeometryInvalid);
        const auto one = math.constant("1"), two = math.constant("2");
        const auto half = math.div(one, two);
        const auto sqrt_two = bounds.root(two);
        auto axes = std::visit(
            [&](const auto& source) { return material_axes(source, one, bounds); }, axes_source);
        const auto fa = multiply(f, axes, RightTensor::AsStored, bounds);
        auto green = green_material(fa, one, half, sqrt_two, bounds);
        auto local = enclose_response(law, input, green, two, sqrt_two, bounds);
        const auto material_stress =
            stress_tensor(local.second_pk_mandel_pascal, sqrt_two, half, bounds);
        auto second = multiply(multiply(axes, material_stress, RightTensor::AsStored, bounds), axes,
                               RightTensor::Transposed, bounds);
        auto first = multiply(f, second, RightTensor::AsStored, bounds);
        auto cauchy = multiply(first, f, RightTensor::Transposed, bounds);
        for (auto& entry : cauchy)
            entry = bounds.divide_positive(entry, det, one);
        local.work = math.info();
        return PgShellMaterialPointResponse{std::move(local), std::move(axes),   std::move(green),
                                            std::move(det),   std::move(second), std::move(first),
                                            std::move(cauchy)};
    });
}

core::Result<PgShellPatchPoint> enclose_pg_shell_patch_point(const PgPreparedPatch& patch,
                                                             const PgShellPatchState& state,
                                                             const PgReferencePoint& query,
                                                             const PgReferenceLimits& limits,
                                                             core::ResourcePlan& plan) noexcept {
    return evaluate([&]() -> core::Result<PgShellPatchPoint> {
        ReferenceArithmetic math{limits, plan};
        math.cancel();
        const auto info = partition_take(patch.info());
        if (info.identity != state.patch)
            partition_fail(core::FailureCode::StaleGeneration);
        const auto q0 = patch.reference_coordinates();
        const auto q = state.coordinates_metres;
        if (q.size() != q0.size() || q.size() != info.basis.rank)
            partition_fail(invalid);
        // Keep both Cartesian vectors alive during reference evaluation. Each
        // returned column then holds eight more rational coefficients.
        const auto coordinates = detail::partition_mul(6, info.basis.rank);
        auto reference_limits = limits;
        math.reserve_nested_rationals(reference_limits.algebra, coordinates);
        auto row =
            partition_take(patch.evaluate_reference(std::span{&query, 1}, reference_limits, plan));
        math.charge(row.call_info.scalar_operations);
        math.visit(row.call_info.partition_visits);
        math.admit(detail::partition_add(384, detail::partition_mul(14, info.basis.rank)));

        auto& reference = row.points.values().front();
        ReferenceRadicals exact{math, reference.g};
        PgReferenceVector position;
        std::array<PgReferenceScalar, 9> gradient;
        for (std::size_t a = 0; a < 3; ++a) {
            position[a] = {math.copy(reference.position_metres[a].rational),
                           math.copy(reference.position_metres[a].radical)};
            gradient[4 * a].rational = math.constant("1");
        }
        // x = X0 + sum N_j (q_j - q0_j), F = I + sum (q_j - q0_j) (G_j)^T.
        // Accumulate exactly before enclosing to retain affine and translation
        // modes of the conditioned basis, including when q0 is much larger.
        for (std::size_t j = 0; j < q.size(); ++j) {
            math.visit(1);
            const auto& shape = row.columns.values()[j];
            for (std::size_t a = 0; a < 3; ++a) {
                if (q[j][a] == q0[j][a])
                    continue;
                const auto delta = math.sub(q[j][a], q0[j][a]);
                position[a] = exact.add(position[a], exact.scale(shape.value, delta));
                add_gradient_column(gradient, a, shape.physical_gradient_per_metre, delta, exact);
            }
        }
        const MaterialBounds bounds{math};
        const auto root = bounds.root(reference.g);
        std::array<Bound, 3> x;
        Tensor f;
        for (std::size_t j = 0; j < x.size(); ++j)
            x[j] = bounds.radical(position[j], root);
        for (std::size_t j = 0; j < f.size(); ++j)
            f[j] = bounds.radical(gradient[j], root);
        return PgShellPatchPoint{info.identity, std::move(reference), std::move(row.columns),
                                 std::move(x),  std::move(f),         math.info()};
    });
}

core::Result<PgShellMaterialPointResponse>
enclose_pg_shell_material_point(const PgShellElasticLaw& law, const PgInitialGreenEnclosure& input,
                                const PgShellPointAxes& axes, const PgShellPatchPoint& point,
                                const PgReferenceLimits& limits,
                                core::ResourcePlan& plan) noexcept {
    return evaluate([&]() -> core::Result<PgShellMaterialPointResponse> {
        ReferenceArithmetic math{limits, plan};
        math.cancel();
        auto point_limits = limits;
        // The reference and shapes remain live for later physical assembly.
        const auto kept = detail::partition_add(88, detail::partition_mul(8, point.shapes.size()));
        math.reserve_nested_rationals(point_limits.algebra, kept);
        const auto source = std::visit(
            [&](const auto& sampled) -> PgShellReferenceAxes {
                using Sample = std::decay_t<decltype(sampled)>;
                if constexpr (std::is_same_v<Sample, PgShellProjectedFibre>)
                    return PgShellProjectedAxes{point.reference, sampled};
                else
                    return sampled;
            },
            axes);
        return enclose_pg_shell_material_point(law, input, source, point.deformation_gradient,
                                               point_limits, plan);
    });
}

core::Result<PgShellMaterialPointLinearization> linearize_pg_shell_material_point(
    const PgShellElasticLaw& law, const PgInitialGreenEnclosure& input,
    const PgShellPointAxes& axes, const PgShellPatchPoint& point,
    const PgShellPatchDirection& direction, const PgReferenceLimits& limits,
    core::ResourcePlan& plan) noexcept {
    return evaluate([&]() -> core::Result<PgShellMaterialPointLinearization> {
        ReferenceArithmetic math{limits, plan};
        math.cancel();
        if (direction.patch != point.patch)
            partition_fail(core::FailureCode::StaleGeneration);
        const auto shapes = point.shapes.values();
        if (direction.coefficients_metres.size() != shapes.size())
            partition_fail(invalid);
        auto value_limits = limits;
        math.reserve_nested_rationals(value_limits.algebra,
                                      detail::partition_mul(3, shapes.size()));
        auto value = partition_take(
            enclose_pg_shell_material_point(law, input, axes, point, value_limits, plan));
        math.charge(value.material.work.scalar_operations);
        math.visit(value.material.work.partition_visits);
        math.admit(detail::partition_add(512, detail::partition_mul(11, shapes.size())));
        std::array<PgReferenceScalar, 9> gradient;
        const auto& g = point.reference.g;
        ReferenceRadicals exact{math, g};
        for (std::size_t j = 0; j < shapes.size(); ++j) {
            math.visit(1);
            for (std::size_t a = 0; a < 3; ++a) {
                const auto& coefficient = direction.coefficients_metres[j][a];
                if (coefficient.sign() != 0)
                    add_gradient_column(gradient, a, shapes[j].physical_gradient_per_metre,
                                        coefficient, exact);
            }
        }
        const MaterialBounds bounds{math};
        const auto root = bounds.root(g);
        Tensor df;
        for (std::size_t j = 0; j < df.size(); ++j)
            df[j] = bounds.radical(gradient[j], root);
        auto derivative = first_pk_derivative(law, value, point.deformation_gradient, df, bounds);
        value.material.work = math.info();
        return PgShellMaterialPointLinearization{std::move(value), std::move(derivative)};
    });
}
} // namespace uvlas::structure
