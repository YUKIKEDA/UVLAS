#pragma once

#include <uvlas/structure/pg_prepared_patch.hpp>

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
class PgShellElasticLaw {
public:
    PgShellElasticLaw(PgShellElasticLaw&&) noexcept = default;
    PgShellElasticLaw& operator=(PgShellElasticLaw&&) noexcept = default;
    PgShellElasticLaw(const PgShellElasticLaw&) = delete;
    PgShellElasticLaw& operator=(const PgShellElasticLaw&) = delete;

    /// Qualify the supplied SI law and density independently of any material point.
    [[nodiscard]] static core::Result<PgShellElasticLaw> prepare(const PgElasticity&,
                                                                 const ExactRational& density_kg_m3,
                                                                 const PgReferenceLimits&,
                                                                 core::ResourcePlan&) noexcept;

    [[nodiscard]] const std::array<PgReferenceScalar, 36>& stiffness() const& noexcept {
        return stiffness_;
    }
    [[nodiscard]] const std::array<PgReferenceScalar, 36>& stiffness() const&& = delete;
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
    PgShellElasticLaw(std::array<PgReferenceScalar, 36>&&, ExactRational&&, ExactRational&&,
                      PgReferenceCallInfo) noexcept;
    std::array<PgReferenceScalar, 36> stiffness_;
    ExactRational density_, correction_squared_;
    PgReferenceCallInfo work_;
};

/// An elastic law and fixed initial Green strain at one material point.
class PgShellMaterial {
public:
    PgShellMaterial(PgShellMaterial&&) noexcept = default;
    PgShellMaterial& operator=(PgShellMaterial&&) noexcept = default;
    PgShellMaterial(const PgShellMaterial&) = delete;
    PgShellMaterial& operator=(const PgShellMaterial&) = delete;

    [[nodiscard]] static core::Result<PgShellMaterial>
    prepare(const PgElasticity&, const ExactRational& density_kg_m3, const PgInitialGreenStrain&,
            const PgReferenceLimits&, core::ResourcePlan&) noexcept;
    [[nodiscard]] const PgShellElasticLaw& law() const& noexcept {
        return law_;
    }
    [[nodiscard]] const PgShellElasticLaw& law() const&& = delete;
    [[nodiscard]] const std::array<PgReferenceScalar, 6>& initial_strain() const& noexcept {
        return initial_;
    }
    [[nodiscard]] const std::array<PgReferenceScalar, 6>& initial_strain() const&& = delete;
    [[nodiscard]] PgReferenceCallInfo call_info() const noexcept {
        return work_;
    }

private:
    PgShellMaterial(PgShellElasticLaw&&, std::array<PgReferenceScalar, 6>&&,
                    PgReferenceCallInfo) noexcept;
    PgShellElasticLaw law_;
    std::array<PgReferenceScalar, 6> initial_;
    PgReferenceCallInfo work_;
};

struct PgShellMaterialResponse {
    PgReferenceScalar energy_joules_per_cubic_metre;
    std::array<PgReferenceScalar, 6> second_pk_mandel_pascal;
    PgReferenceCallInfo work;
};

struct PgMaterialEnclosure {
    ExactRational lower, upper;
};

struct PgInitialGreenEnclosure {
    const ExactRational& radicand;
    const std::array<PgReferenceScalar, 6>& mandel;
    const std::array<ExactRational, 6>& absolute_error;
};

struct PgShellMaterialEnclosure {
    PgMaterialEnclosure energy_joules_per_cubic_metre;
    std::array<PgMaterialEnclosure, 6> second_pk_mandel_pascal;
    PgReferenceCallInfo work;
};

/// Row-major Cartesian tensor. Axis matrices store the basis as columns.
using PgShellTensorEnclosure = std::array<PgMaterialEnclosure, 9>;

struct PgShellExplicitAxes {
    const std::array<ExactRational, 9>& columns;
};

struct PgShellProjectedFibre {
    const std::array<ExactRational, 3>& world_fibre;
    /// Exact unit pair for the right-handed in-plane rotation.
    const ExactRational& cosine;
    const ExactRational& sine;
};

struct PgShellProjectedAxes {
    const PgReferenceValue& reference;
    PgShellProjectedFibre fibre;
};

struct PgShellFieldFibre {
    const PgReferenceVector& world_fibre;
    const std::array<ExactRational, 3>& direction_error;
    const PgMaterialEnclosure& cosine;
    const PgMaterialEnclosure& sine;
};

struct PgShellFieldAxes {
    const std::array<PgReferenceScalar, 9>& columns;
};

using PgShellSampledAxes = std::variant<PgShellFieldFibre, PgShellFieldAxes>;

struct PgShellReferencedAxes {
    const PgReferenceValue& reference;
    PgShellSampledAxes field;
};

using PgShellReferenceAxes =
    std::variant<PgShellExplicitAxes, PgShellProjectedAxes, PgShellReferencedAxes>;
using PgShellPointAxes =
    std::variant<PgShellExplicitAxes, PgShellProjectedFibre, PgShellFieldFibre, PgShellFieldAxes>;

struct PgShellPatchState {
    PgPreparedPatchIdentity patch;
    std::span<const PgCartesianExactPosition> coordinates_metres;
};

struct PgShellPatchPoint {
    PgPreparedPatchIdentity patch;
    PgReferenceValue reference;
    core::PlannedBuffer<PgConditionedReferenceColumn> shapes;
    std::array<PgMaterialEnclosure, 3> current_position_metres;
    PgShellTensorEnclosure deformation_gradient;
    PgReferenceCallInfo work;
};

/// Assemble x and F at a volume or one-sided boundary query from the exact
/// conditioned columns and the current and reference Cartesian coefficients.
[[nodiscard]] core::Result<PgShellPatchPoint>
enclose_pg_shell_patch_point(const PgPreparedPatch&, const PgShellPatchState&,
                             const PgReferencePoint&, const PgReferenceLimits&,
                             core::ResourcePlan&) noexcept;

struct PgShellMaterialPointResponse {
    PgShellMaterialEnclosure material;
    PgShellTensorEnclosure reference_axes;
    std::array<PgMaterialEnclosure, 6> green_material_mandel;
    PgMaterialEnclosure determinant;
    PgShellTensorEnclosure second_pk_world_pascal;
    PgShellTensorEnclosure first_pk_world_pascal;
    PgShellTensorEnclosure cauchy_world_pascal;
};

struct PgShellPatchDirection {
    PgPreparedPatchIdentity patch;
    std::span<const PgCartesianExactPosition> coefficients_metres;
};

struct PgShellMaterialPointLinearization {
    PgShellMaterialPointResponse value;
    PgShellTensorEnclosure first_pk_direction_pascal;
};

/// Evaluate all stress measures at one point with positive det(F).
/// F maps reference world increments to current world increments.
[[nodiscard]] core::Result<PgShellMaterialPointResponse>
enclose_pg_shell_material_point(const PgShellElasticLaw&, const PgInitialGreenEnclosure&,
                                const PgShellReferenceAxes&,
                                const PgShellTensorEnclosure& deformation_gradient,
                                const PgReferenceLimits&, core::ResourcePlan&) noexcept;

/// Use the same physical point and its retained layer tangents for projection.
[[nodiscard]] core::Result<PgShellMaterialPointResponse>
enclose_pg_shell_material_point(const PgShellElasticLaw&, const PgInitialGreenEnclosure&,
                                const PgShellPointAxes&, const PgShellPatchPoint&,
                                const PgReferenceLimits&, core::ResourcePlan&) noexcept;

/// Differentiate P=F S along one coefficient direction at the same fixed
/// material point, including both geometric and constitutive terms.
[[nodiscard]] core::Result<PgShellMaterialPointLinearization>
linearize_pg_shell_material_point(const PgShellElasticLaw&, const PgInitialGreenEnclosure&,
                                  const PgShellPointAxes&, const PgShellPatchPoint&,
                                  const PgShellPatchDirection&, const PgReferenceLimits&,
                                  core::ResourcePlan&) noexcept;

struct PgShellVolumeSample {
    PgDiscreteVolumePoint quadrature;
    core::LayerId layer;
    const PgShellElasticLaw& law;
    PgInitialGreenEnclosure initial;
    PgShellPointAxes axes;
};

/// Exact reference evaluation used by elastic volume assembly.
enum class PgShellReferencePath : std::uint8_t { Direct, Composed };

struct PgShellVolumeResolution {
    std::uint32_t fractional_bits; ///< Outward dyadic rounding of each SI contribution.
    const ExactRational& energy_width_joules;
    const ExactRational& force_width_newtons;
    const ExactRational& tangent_width_newtons;
    /// Composed permits counted preparation and falls back at unsupported points.
    PgShellReferencePath reference_path = PgShellReferencePath::Direct;
};

using PgShellCartesianEnclosure = std::array<PgMaterialEnclosure, 3>;

struct PgShellElasticVolumeAction {
    PgPreparedPatchIdentity patch;
    PgMaterialEnclosure energy_joules;
    core::PlannedBuffer<PgShellCartesianEnclosure> internal_force_newtons;
    core::PlannedBuffer<PgShellCartesianEnclosure> tangent_action_newtons;
    PgReferenceCallInfo work;
};

/// Assemble the energy, its gradient and an optional tangent action for one
/// finite rule. Positive samples must cover each partition volume in order,
/// and each layer sample belongs to its volume. Width limits refer to the
/// finite-rule result, including field and arithmetic errors.
[[nodiscard]] core::Result<PgShellElasticVolumeAction> assemble_pg_shell_elastic_volume(
    const PgPreparedPatch&, const PgShellPatchState&, std::span<const PgShellVolumeSample>,
    const std::optional<PgShellPatchDirection>&, const PgShellVolumeResolution&,
    const PgReferenceLimits&, core::ResourcePlan&) noexcept;

/// Same assembly for budgeted descriptor storage. Every slot must be bound.
[[nodiscard]] core::Result<PgShellElasticVolumeAction>
assemble_pg_shell_elastic_volume(const PgPreparedPatch&, const PgShellPatchState&,
                                 std::span<const std::optional<PgShellVolumeSample>>,
                                 const std::optional<PgShellPatchDirection>&,
                                 const PgShellVolumeResolution&, const PgReferenceLimits&,
                                 core::ResourcePlan&) noexcept;

struct PgShellLayerDensity {
    core::LayerId layer;
    const ExactRational& kilograms_per_cubic_metre;
};

struct PgShellMassResolution {
    std::uint32_t fractional_bits;
    const ExactRational& mass_width_kg;
    const ExactRational& action_width_kg_metres;
};

struct PgShellMassVolumeAction {
    PgPreparedPatchIdentity patch;
    PgMaterialEnclosure mass_kg;
    core::PlannedBuffer<PgShellCartesianEnclosure> action_kg_metres;
    PgReferenceCallInfo work;
};

/// Apply the consistent reference-volume mass on a separately supplied rule.
/// Densities are in the partition's physical layer order. The trial positions
/// and elastic fields do not enter the action.
[[nodiscard]] core::Result<PgShellMassVolumeAction>
apply_pg_shell_consistent_mass(const PgPreparedPatch&, std::span<const PgDiscreteVolumePoint>,
                               std::span<const PgShellLayerDensity>, const PgShellPatchDirection&,
                               const PgShellMassResolution&, const PgReferenceLimits&,
                               core::ResourcePlan&) noexcept;

[[nodiscard]] core::Result<PgShellMassVolumeAction>
apply_pg_shell_consistent_mass(const PgPreparedPatch&, std::span<const PgOwnedDiscreteVolumePoint>,
                               std::span<const std::optional<PgShellLayerDensity>>,
                               const PgShellPatchDirection&, const PgShellMassResolution&,
                               const PgReferenceLimits&, core::ResourcePlan&) noexcept;

/// Evaluate the whole input uncertainty in the fixed reference material frame.
/// Require a positive natural metric throughout the initial strain enclosure.
[[nodiscard]] core::Result<PgShellMaterialEnclosure>
enclose_pg_shell_material(const PgShellElasticLaw&, const PgInitialGreenEnclosure&,
                          const std::array<PgMaterialEnclosure, 6>& green_mandel,
                          const PgReferenceLimits&, core::ResourcePlan&) noexcept;

/// Evaluate the reference-volume law in the fixed material axes. Each
/// Green strain component is an exact Mandel value a+b sqrt(2). The stress
/// tangent with respect to these components is the prepared stiffness.
[[nodiscard]] core::Result<PgShellMaterialResponse>
evaluate_pg_shell_material(const PgShellMaterial&,
                           const std::array<PgReferenceScalar, 6>& green_mandel,
                           const PgReferenceLimits&, core::ResourcePlan&) noexcept;

/// Evaluate an explicitly sampled initial field against a reusable elastic law.
[[nodiscard]] core::Result<PgShellMaterialResponse>
evaluate_pg_shell_material(const PgShellElasticLaw&, const PgInitialGreenStrain&,
                           const std::array<PgReferenceScalar, 6>& green_mandel,
                           const PgReferenceLimits&, core::ResourcePlan&) noexcept;

} // namespace uvlas::structure
