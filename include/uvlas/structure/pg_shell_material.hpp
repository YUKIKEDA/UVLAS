#pragma once

#include <uvlas/structure/pg_reference_geometry.hpp>

namespace uvlas::structure {

/// Order: xx, yy, zz, yz, xz, xy. The tag supplies the shear scaling.
enum class PgMaterialComponents : std::uint8_t { Mandel, Engineering };

struct PgIsotropicElasticity {
    const ExactRational& young_pascal;
    const ExactRational& poisson;
};

struct PgOrthotropicElasticity {
    const std::array<ExactRational, 3>& young_pascal; ///< E1, E2, E3.
    const std::array<ExactRational, 3>& poisson;      ///< Nu12, nu13, nu23.
    const std::array<ExactRational, 3>& shear_pascal; ///< G23, G13, G12.
};

struct PgEffectiveElasticity {
    const std::array<ExactRational, 36>& stiffness_pascal; ///< Row-major.
    PgMaterialComponents components;
};

using PgElasticity =
    std::variant<PgIsotropicElasticity, PgOrthotropicElasticity, PgEffectiveElasticity>;

struct PgInitialGreenStrain {
    const std::array<ExactRational, 6>& components;
    PgMaterialComponents convention;
};

/// Reference-volume elastic law shared by all three material input forms.
/// Each stored Mandel value represents rational + radical * sqrt(2) exactly.
class PgShellMaterial {
public:
    PgShellMaterial(PgShellMaterial&&) noexcept = default;
    PgShellMaterial& operator=(PgShellMaterial&&) noexcept = default;
    PgShellMaterial(const PgShellMaterial&) = delete;
    PgShellMaterial& operator=(const PgShellMaterial&) = delete;

    /// Qualify the supplied SI law, density and natural metric before publishing it.
    /// The input spans and values are borrowed only during this call.
    [[nodiscard]] static core::Result<PgShellMaterial>
    prepare(const PgElasticity&, const ExactRational& density_kg_m3, const PgInitialGreenStrain&,
            const PgReferenceLimits&, core::ResourcePlan&) noexcept;

    [[nodiscard]] const std::array<PgReferenceScalar, 36>& stiffness() const& noexcept {
        return stiffness_;
    }
    [[nodiscard]] const std::array<PgReferenceScalar, 36>& stiffness() const&& = delete;
    [[nodiscard]] const std::array<PgReferenceScalar, 6>& initial_strain() const& noexcept {
        return initial_;
    }
    [[nodiscard]] const std::array<PgReferenceScalar, 6>& initial_strain() const&& = delete;
    [[nodiscard]] const ExactRational& density_kg_m3() const& noexcept {
        return density_;
    }
    [[nodiscard]] const ExactRational& density_kg_m3() const&& = delete;
    /// Squared relative Frobenius correction in the normalized Mandel components.
    [[nodiscard]] const ExactRational& symmetry_correction_squared() const& noexcept {
        return correction_squared_;
    }
    [[nodiscard]] const ExactRational& symmetry_correction_squared() const&& = delete;
    [[nodiscard]] PgReferenceCallInfo call_info() const noexcept {
        return work_;
    }

private:
    PgShellMaterial(std::array<PgReferenceScalar, 36>&&, std::array<PgReferenceScalar, 6>&&,
                    ExactRational&&, ExactRational&&, PgReferenceCallInfo) noexcept;
    std::array<PgReferenceScalar, 36> stiffness_;
    std::array<PgReferenceScalar, 6> initial_;
    ExactRational density_, correction_squared_;
    PgReferenceCallInfo work_;
};

struct PgShellMaterialResponse {
    PgReferenceScalar energy_joules_per_cubic_metre;
    std::array<PgReferenceScalar, 6> second_pk_mandel_pascal;
    PgReferenceCallInfo work;
};

/// Evaluate the reference-volume law in the fixed material axes. Each
/// Green strain component is an exact Mandel value a+b sqrt(2). The stress
/// tangent with respect to these components is the prepared stiffness.
[[nodiscard]] core::Result<PgShellMaterialResponse>
evaluate_pg_shell_material(const PgShellMaterial&,
                           const std::array<PgReferenceScalar, 6>& green_mandel,
                           const PgReferenceLimits&, core::ResourcePlan&) noexcept;

} // namespace uvlas::structure
