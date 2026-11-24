#pragma once

#include "nastran_beam_cards.hpp"
#include "nastran_isotropic_shell_cards.hpp"

namespace uvlas::io::detail {

using NastranPbeamReal = std::optional<NastranShellReal>;

struct NastranPbeamSection {
    NastranPbeamReal area, i1, i2, i12, torsion, nonstructural_mass;
};

struct NastranPbeamRecovery {
    // Pairs of local y,z coordinates for C, D, E and F, relative to the shear center.
    std::array<NastranPbeamReal, 8> locations;
    NastranBulkPosition position;
};

enum class NastranPbeamOutput : std::uint8_t { No, SectionPoints, EndAPoints };

struct NastranPbeamStation {
    NastranPbeamOutput output = NastranPbeamOutput::No;
    NastranBulkLexeme output_field;
    NastranPbeamReal fraction; // X/XB, retained as the source decimal.
    NastranPbeamSection section;
    std::optional<NastranPbeamRecovery> recovery;
};

struct NastranPbeamModifiers {
    NastranPbeamReal shear_factor_1, shear_factor_2, taper_relief_1, taper_relief_2;
    std::array<NastranPbeamReal, 2> nonstructural_inertia, warping; // A, B.
};

struct NastranPbeamOffsets {
    // (y,z) for end A followed by (y,z) for end B, relative to the shear center.
    std::array<NastranPbeamReal, 4> nonstructural_mass_center, neutral_axis;
};

struct NastranPbeam {
    NastranPbeamDeclaration declaration;
    NastranMaterialReference material;
    NastranPbeamSection end_a;
    std::optional<NastranPbeamRecovery> recovery_a;
    std::size_t first_station = 0, station_count = 0;
    std::optional<std::size_t> end_b_station; // Ordinal within this property.
    std::optional<NastranPbeamModifiers> modifiers;
    std::optional<NastranPbeamOffsets> offsets;
};

enum class NastranPbeamFailure : std::uint8_t {
    Input,
    SourceMismatch,
    MissingField,
    InvalidField,
    UnsupportedField,
    ExtraField,
    DuplicateStation,
    MissingEndB,
    UndefinedMaterial,
    WrongMaterialKind,
    Limit
};

struct NastranPbeamDiagnostic {
    NastranPbeamFailure reason = NastranPbeamFailure::Input;
    NastranBulkPosition position;
};

struct NastranPbeamLimits {
    std::uint64_t properties;
    std::uint64_t stations;
    std::uint64_t work;
};

struct NastranPbeamView {
    NastranBeamOrigin source;
    // Same order as NastranBeamView::properties, indexed by a CBEAM's PBEAM declaration ordinal.
    std::span<const NastranPbeam> properties;
    std::span<const NastranPbeamStation> stations;
    std::uint64_t work;

    [[nodiscard]] std::span<const NastranPbeamStation> sections(const NastranPbeam&) const noexcept;
};

// Owns the declared sections and lexical positions. Defaults, interpolation, sectional
// admissibility, beam geometry, stress recovery, units and K6 are computed at later stages.
class NastranPbeamCards {
public:
    NastranPbeamCards(NastranPbeamCards&&) noexcept;
    NastranPbeamCards& operator=(NastranPbeamCards&&) = delete;

    [[nodiscard]] static core::Result<NastranPbeamCards>
    prepare(const NastranIsotropicShellCards&, const NastranBeamCards&, NastranPbeamLimits,
            core::ResourcePlan&, const std::stop_token&,
            std::optional<NastranPbeamDiagnostic>* = nullptr) noexcept;

    [[nodiscard]] std::optional<NastranPbeamView> view() const& noexcept;
    [[nodiscard]] std::optional<NastranPbeamView> view() const&& = delete;

private:
    NastranPbeamCards(NastranBeamOrigin, core::PlannedBuffer<NastranPbeam>&&,
                      core::PlannedBuffer<NastranPbeamStation>&&, std::uint64_t) noexcept;

    std::optional<NastranBeamOrigin> source_;
    core::PlannedBuffer<NastranPbeam> properties_;
    core::PlannedBuffer<NastranPbeamStation> stations_;
    std::uint64_t work_;
};

} // namespace uvlas::io::detail
