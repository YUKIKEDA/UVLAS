#pragma once

#include "nastran_shell_cards.hpp"

namespace uvlas::io::detail {

struct NastranPbeamDeclaration {
    NastranPropertyId id{0};
    NastranBulkLexeme id_field;
    std::size_t card_ordinal = 0;
};

struct NastranPbeamReference {
    NastranPropertyId id{0};
    std::optional<NastranBulkLexeme> field; // Blank PID denotes EID.
    std::size_t declaration_ordinal = 0;
};

struct NastranBeamGridReference {
    NastranGridId id{0};
    NastranBulkLexeme field;
    std::size_t grid_ordinal = 0;
};

enum class NastranBeamOrientationBasis : std::uint8_t { GridDisplacement, Basic };
enum class NastranBeamOffsetBasis : std::uint8_t { GridDisplacement, OffsetSystem };

struct NastranBeamInterpretation {
    NastranBeamOrientationBasis orientation = NastranBeamOrientationBasis::GridDisplacement;
    NastranBeamOffsetBasis a = NastranBeamOffsetBasis::GridDisplacement;
    NastranBeamOffsetBasis b = NastranBeamOffsetBasis::GridDisplacement;
    std::optional<NastranBulkLexeme> field;
};

using NastranBeamVector = std::array<std::optional<NastranShellReal>, 3>;
using NastranBeamOrientation = std::variant<NastranBeamVector, NastranBeamGridReference>;

struct NastranBeamReleases {
    std::uint8_t element_components = 0;
    std::optional<NastranBulkLexeme> field;
};

struct NastranBeamWarpingPoint {
    NastranGridId id{0};
    NastranBulkLexeme field;
    // The adopted GRID form uses its external T1. Zero has no referenced grid.
    std::optional<std::size_t> grid_ordinal;
};

struct NastranBeamEnd {
    NastranBeamGridReference point;
    NastranBeamReleases releases;
    NastranBeamVector offset; // Grid point to the end of the shear-center axis.
    std::optional<NastranBeamWarpingPoint> warping;
};

struct NastranCbeam {
    NastranElementId id{0};
    NastranBulkLexeme id_field;
    std::size_t card_ordinal = 0;
    NastranPbeamReference property;
    NastranBeamEnd a, b;
    NastranBeamOrientation orientation;
    NastranBeamInterpretation interpretation;
};

enum class NastranBeamFailure : std::uint8_t {
    Input,
    MissingField,
    InvalidField,
    ExtraField,
    UnsupportedField,
    DuplicateEndpoint,
    OrientationEndpoint,
    UndefinedGrid,
    UndefinedWarpingPoint,
    DuplicateProperty,
    UndefinedProperty,
    WrongPropertyKind,
    Limit
};

struct NastranBeamDiagnostic {
    NastranBeamFailure reason = NastranBeamFailure::Input;
    NastranBulkPosition position;
};

struct NastranBeamLimits {
    std::uint64_t beams;
    std::uint64_t properties;
    std::uint64_t work;
};

struct NastranBeamOrigin {
    RootSnapshotLocation location;
    SourceDigest digest;
};

struct NastranBeamView {
    NastranBeamOrigin source;
    std::span<const NastranCbeam> beams;
    std::span<const NastranPbeamDeclaration> properties;
    std::uint64_t work;
};

// Owns the CBEAM records and PBEAM identities; the source is read only during
// qualification. Beam frames, K6, units, releases and warping physics are separate.
class NastranBeamCards {
public:
    NastranBeamCards(NastranBeamCards&&) noexcept;
    NastranBeamCards& operator=(NastranBeamCards&&) = delete;

    [[nodiscard]] static core::Result<NastranBeamCards>
    prepare(const NastranShellCards&, NastranBeamLimits, core::ResourcePlan&,
            const std::stop_token&, std::optional<NastranBeamDiagnostic>* = nullptr) noexcept;
    [[nodiscard]] std::optional<NastranBeamView> view() const& noexcept;
    [[nodiscard]] std::optional<NastranBeamView> view() const&& = delete;

private:
    NastranBeamCards(NastranBeamOrigin, core::PlannedBuffer<NastranCbeam>&&,
                     core::PlannedBuffer<NastranPbeamDeclaration>&&, std::uint64_t) noexcept;

    std::optional<NastranBeamOrigin> source_;
    core::PlannedBuffer<NastranCbeam> beams_;
    core::PlannedBuffer<NastranPbeamDeclaration> properties_;
    std::uint64_t work_;
};

} // namespace uvlas::io::detail
