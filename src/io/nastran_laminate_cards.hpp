#pragma once

#include "nastran_isotropic_shell_cards.hpp"

namespace uvlas::io::detail {

struct NastranMat8 {
    std::size_t declaration_ordinal = 0;
    std::optional<NastranShellReal> longitudinal_young, transverse_young, poisson_12;
    std::optional<NastranShellReal> shear_12, shear_1z, shear_2z, density;
    std::optional<NastranShellReal> expansion_1, expansion_2, reference_temperature;
    std::optional<NastranShellReal> longitudinal_tension, longitudinal_compression;
    std::optional<NastranShellReal> transverse_tension, transverse_compression, shear_limit;
    std::optional<NastranShellReal> damping, tsai_wu_interaction, strain_allowables;
};

enum class NastranLaminateMode : std::uint8_t {
    FullStack,
    SymmetricHalf,
    MembraneOnly,
    BendingOnly,
    Smeared,
    SmearedWithCore
};

enum class NastranPcompFailureTheory : std::uint8_t {
    Hill,
    Hoffman,
    TsaiWu,
    MaximumStrain,
    Hashin,
    HashinTape,
    HashinFabric
};

struct NastranPcompTheoryField {
    NastranPcompFailureTheory theory = NastranPcompFailureTheory::Hill;
    NastranBulkLexeme field;
};

struct NastranPcompStressOutput {
    bool requested = false;
    NastranBulkLexeme field;
};

struct NastranPcompPly {
    NastranBulkLexeme source; // First nonblank field of this ply.
    std::optional<NastranMaterialReference> material;
    std::optional<NastranShellReal> thickness, angle;
    std::optional<NastranPcompStressOutput> stress_output;
    // Indices within this PCOMP, pointing to a ply with the corresponding explicit value.
    std::size_t material_source_ply = 0;
    std::size_t thickness_source_ply = 0;
};

struct NastranPcomp {
    std::size_t property_ordinal = 0;
    std::size_t first_ply = 0;
    std::size_t ply_count = 0;
    NastranLaminateMode mode = NastranLaminateMode::FullStack;
    std::optional<NastranBulkLexeme> mode_field;
    std::optional<NastranPcompTheoryField> failure_theory;
    std::optional<NastranShellReal> bottom_offset, nonstructural_mass, bond_shear_limit;
    std::optional<NastranShellReal> reference_temperature, damping;
};

enum class NastranLaminateFailure : std::uint8_t {
    Input,
    MissingField,
    InvalidField,
    UndefinedMaterial,
    UnsupportedField,
    Limit
};

struct NastranLaminateDiagnostic {
    NastranLaminateFailure reason = NastranLaminateFailure::Input;
    NastranBulkPosition position;
};

struct NastranLaminateLimits {
    std::uint64_t plies;
    std::uint64_t work;
};

struct NastranLaminateView {
    NastranIsotropicShellView input;
    std::span<const NastranMat8> mat8;
    std::span<const NastranPcomp> pcomp;
    std::span<const NastranPcompPly> plies; // Stored in declared order, from the bottom up.
    std::uint64_t work;

    [[nodiscard]] std::span<const NastranPcompPly> layers(const NastranPcomp&) const noexcept;
};

// Owns the exact source and declared laminate records. Indices encode source
// inheritance; physical ply expansion, default thickness/offset and SI data
// are determined after this stage.
class NastranLaminateCards {
public:
    NastranLaminateCards(NastranLaminateCards&&) noexcept = default;
    NastranLaminateCards& operator=(NastranLaminateCards&&) = delete;

    [[nodiscard]] static core::Result<NastranLaminateCards>
    prepare(NastranIsotropicShellCards&&, NastranLaminateLimits, core::ResourcePlan&,
            const std::stop_token&, std::optional<NastranLaminateDiagnostic>* = nullptr) noexcept;
    [[nodiscard]] std::optional<NastranLaminateView> view() const& noexcept;
    [[nodiscard]] std::optional<NastranLaminateView> view() const&& = delete;
    [[nodiscard]] const NastranIsotropicShellCards& isotropic() const& noexcept {
        return input_;
    }
    [[nodiscard]] const NastranIsotropicShellCards& isotropic() const&& = delete;

private:
    NastranLaminateCards(NastranIsotropicShellCards&&, core::PlannedBuffer<NastranMat8>&&,
                         core::PlannedBuffer<NastranPcomp>&&,
                         core::PlannedBuffer<NastranPcompPly>&&, std::uint64_t) noexcept;
    NastranIsotropicShellCards input_;
    core::PlannedBuffer<NastranMat8> mat8_;
    core::PlannedBuffer<NastranPcomp> pcomp_;
    core::PlannedBuffer<NastranPcompPly> plies_;
    std::uint64_t work_;
};

} // namespace uvlas::io::detail
