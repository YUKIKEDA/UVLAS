#pragma once

#include "nastran_shell_cards.hpp"

namespace uvlas::io::detail {

struct NastranMaterialTag;
using NastranMaterialId = core::Identifier<NastranMaterialTag>;
enum class NastranMaterialKind : std::uint8_t { Mat1, Mat8 };

struct NastranMaterialDeclaration {
    NastranMaterialId id{0};
    NastranMaterialKind kind = NastranMaterialKind::Mat1;
    NastranBulkLexeme id_field;
    std::size_t card_ordinal = 0;
    std::optional<std::size_t> mat1_ordinal;
};

struct NastranMaterialReference {
    NastranMaterialId id{0};
    NastranBulkLexeme field;
    std::size_t declaration_ordinal = 0;
};

struct NastranMat1 {
    std::size_t declaration_ordinal = 0;
    std::optional<NastranShellReal> young, shear, poisson, density;
    std::optional<NastranShellReal> expansion, reference_temperature, damping;
    std::optional<NastranShellReal> tension_limit, compression_limit, shear_limit;
    std::optional<NastranShellMaterialAxis> material_coordinate;
};

enum class NastranSpecialBending : std::uint8_t { Zero, PlaneStrain };

struct NastranBendingDesignation {
    NastranSpecialBending kind = NastranSpecialBending::Zero;
    NastranBulkLexeme field;
};

using NastranBendingMaterial =
    std::variant<std::monostate, NastranMaterialReference, NastranBendingDesignation>;

struct NastranPshell {
    std::size_t property_ordinal = 0;
    std::optional<NastranMaterialReference> membrane;
    NastranBendingMaterial bending;
    std::optional<NastranMaterialReference> transverse_shear, coupling;
    std::optional<NastranShellReal> thickness, bending_ratio, shear_ratio, nonstructural_mass;
    std::optional<NastranShellReal> first_fiber, second_fiber;
    bool distinct_materials = false;
    bool orthotropic_material = false;
};

enum class NastranIsotropicShellFailure : std::uint8_t {
    Input,
    MissingField,
    InvalidField,
    ExtraField,
    DuplicateMaterial,
    UndefinedMaterial,
    UndefinedCoordinate,
    IncompatibleMaterial,
    Limit
};

struct NastranIsotropicShellDiagnostic {
    NastranIsotropicShellFailure reason = NastranIsotropicShellFailure::Input;
    NastranBulkPosition position;
};

struct NastranIsotropicShellLimits {
    std::uint64_t material_ids; // MAT1 and MAT8 share the material namespace.
    std::uint64_t work;
};

struct NastranIsotropicShellView {
    NastranShellView input;
    std::span<const NastranMaterialDeclaration> materials;
    std::span<const NastranMat1> mat1;
    std::span<const NastranPshell> pshell;
    std::uint64_t work;
};

// Owns exact MAT1/PSHELL fields, including checked material references and their card positions.
// A retained source value or role gives no guarantee of equivalent three-dimensional elasticity.
class NastranIsotropicShellCards {
public:
    NastranIsotropicShellCards(NastranIsotropicShellCards&&) noexcept = default;
    NastranIsotropicShellCards& operator=(NastranIsotropicShellCards&&) = delete;

    [[nodiscard]] static core::Result<NastranIsotropicShellCards>
    prepare(NastranShellCards&&, NastranIsotropicShellLimits, core::ResourcePlan&,
            const std::stop_token&,
            std::optional<NastranIsotropicShellDiagnostic>* = nullptr) noexcept;
    [[nodiscard]] std::optional<NastranIsotropicShellView> view() const& noexcept;
    [[nodiscard]] std::optional<NastranIsotropicShellView> view() const&& = delete;

private:
    NastranIsotropicShellCards(NastranShellCards&&,
                               core::PlannedBuffer<NastranMaterialDeclaration>&&,
                               core::PlannedBuffer<NastranMat1>&&,
                               core::PlannedBuffer<NastranPshell>&&, std::uint64_t) noexcept;
    NastranShellCards input_;
    core::PlannedBuffer<NastranMaterialDeclaration> materials_;
    core::PlannedBuffer<NastranMat1> mat1_;
    core::PlannedBuffer<NastranPshell> pshell_;
    std::uint64_t work_;
};

} // namespace uvlas::io::detail
