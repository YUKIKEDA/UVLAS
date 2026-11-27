#include "pg_reference_arithmetic.hpp"

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
} // namespace

PgShellMaterial::PgShellMaterial(std::array<PgReferenceScalar, 36>&& stiffness,
                                 std::array<PgReferenceScalar, 6>&& initial, R&& density,
                                 R&& correction_squared, PgReferenceCallInfo work) noexcept
    : stiffness_(std::move(stiffness)),
      initial_(std::move(initial)),
      density_(std::move(density)),
      correction_squared_(std::move(correction_squared)),
      work_(work) {}

core::Result<PgShellMaterial> PgShellMaterial::prepare(const PgElasticity& input, const R& density,
                                                       const PgInitialGreenStrain& initial,
                                                       const PgReferenceLimits& limits,
                                                       core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgShellMaterial> {
        ReferenceArithmetic math{limits, plan};
        math.cancel();
        if (density.sign() <= 0)
            partition_fail(invalid);
        const bool initial_engineering = engineering(initial.convention);
        // Covers borrowed inputs, compliance/adjugate, normalization and LDL copies,
        // the two coefficients of every normalized value and the natural metric.
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
        std::array<PgReferenceScalar, order> strain;
        for (std::size_t i = 0; i < order; ++i) {
            if (initial_engineering && i >= 3)
                strain[i].radical = math.div(initial.components[i], two);
            else
                strain[i].rational = math.copy(initial.components[i]);
        }
        qualify_natural_metric(strain, one, two, math);
        auto rho = math.copy(density);
        return PgShellMaterial{std::move(stiffness), std::move(strain), std::move(rho),
                               std::move(law.correction_squared), math.info()};
    });
}

core::Result<PgShellMaterialResponse>
evaluate_pg_shell_material(const PgShellMaterial& material,
                           const std::array<PgReferenceScalar, 6>& green,
                           const PgReferenceLimits& limits, core::ResourcePlan& plan) noexcept {
    auto result = partition_entry([&]() -> core::Result<PgShellMaterialResponse> {
        ReferenceArithmetic math{limits, plan};
        math.cancel();
        if (material.density_kg_m3().sign() <= 0)
            partition_fail(invalid);
        math.admit();
        const auto two = math.constant("2");
        ReferenceRadicals exact{math, two};
        std::array<PgReferenceScalar, order> elastic, stress;
        for (std::size_t i = 0; i < order; ++i)
            elastic[i] = exact.sub(green[i], material.initial_strain()[i]);
        PgReferenceScalar doubled_energy;
        for (std::size_t i = 0; i < order; ++i) {
            for (std::size_t j = 0; j < order; ++j)
                stress[i] = exact.add(stress[i],
                                      exact.mul(material.stiffness()[i * order + j], elastic[j]));
            doubled_energy = exact.add(doubled_energy, exact.mul(elastic[i], stress[i]));
        }
        PgReferenceScalar energy{math.div(doubled_energy.rational, two),
                                 math.div(doubled_energy.radical, two)};
        return PgShellMaterialResponse{std::move(energy), std::move(stress), math.info()};
    });
    if (!result)
        result.error().operation = core::Operation::EvaluateStructure;
    return result;
}
} // namespace uvlas::structure
