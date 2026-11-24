#pragma once

#include "nastran_coordinate_cards.hpp"

#include <variant>

namespace uvlas::io::detail {

struct NastranElementTag;
struct NastranPropertyTag;
using NastranElementId = core::Identifier<NastranElementTag>;
using NastranPropertyId = core::Identifier<NastranPropertyTag>;

enum class NastranShellPropertyKind : std::uint8_t { Pshell, Pcomp };
enum class NastranShellThicknessMode : std::uint8_t { Absolute, Fraction };

struct NastranShellProperty {
    NastranPropertyId id{0};
    NastranShellPropertyKind kind = NastranShellPropertyKind::Pshell;
    NastranBulkLexeme id_field;
    std::size_t card_ordinal = 0;
};

struct NastranShellPropertyReference {
    NastranPropertyId id{0};
    std::optional<NastranBulkLexeme> field;
    std::size_t property_ordinal = 0;
};

struct NastranShellReal {
    NastranBulkLexeme field;
    NastranRealField value;
};

struct NastranShellCorner {
    NastranGridId id{0};
    NastranBulkLexeme field;
    std::size_t grid_ordinal = 0;
    std::optional<NastranShellReal> thickness;
};

struct NastranShellMaterialAxis {
    NastranCoordinateId id{0};
    NastranBulkLexeme field;
    std::optional<std::size_t> coordinate_ordinal; // Absent for basic (MCID 0).
};

// Absent orientation denotes the implicit zero THETA. Explicit reals are degrees.
using NastranShellOrientation =
    std::variant<std::monostate, NastranShellReal, NastranShellMaterialAxis>;
using NastranShellConnectivity =
    std::variant<std::array<NastranShellCorner, 3>, std::array<NastranShellCorner, 4>>;

struct NastranShellElement {
    NastranElementId id{0};
    NastranBulkLexeme id_field;
    std::size_t card_ordinal = 0;
    NastranShellPropertyReference property;
    NastranShellConnectivity corners;
    NastranShellOrientation orientation;
    std::optional<NastranShellReal> offset;
    NastranShellThicknessMode thickness_mode = NastranShellThicknessMode::Absolute;
    std::optional<NastranBulkLexeme> thickness_mode_field;
};

enum class NastranShellFailure : std::uint8_t {
    Input,
    MissingField,
    InvalidField,
    ExtraField,
    DuplicateElement,
    DuplicateProperty,
    DuplicateGrid,
    UndefinedGrid,
    UndefinedCoordinate,
    UndefinedProperty,
    WrongPropertyKind,
    Limit
};

struct NastranShellDiagnostic {
    NastranShellFailure reason = NastranShellFailure::Input;
    NastranBulkPosition position;
};

struct NastranShellLimits {
    std::uint64_t element_ids;  // Includes CBEAM, CONM2 and RBE2 identity checks.
    std::uint64_t property_ids; // Includes PBEAM identity checks.
    std::uint64_t work;
};

struct NastranShellView {
    NastranCoordinateView input;
    std::span<const NastranShellElement> elements;
    std::span<const NastranShellProperty> properties;
    std::uint64_t work;
};

// Owns the authenticated source and its typed shell declarations. Elemental
// geometry, units, material properties and property-driven thickness are deferred.
class NastranShellCards {
public:
    NastranShellCards(NastranShellCards&&) noexcept = default;
    NastranShellCards& operator=(NastranShellCards&&) = delete;

    [[nodiscard]] static core::Result<NastranShellCards>
    prepare(NastranCoordinateCards&&, NastranShellLimits, core::ResourcePlan&,
            const std::stop_token&, std::optional<NastranShellDiagnostic>* = nullptr) noexcept;
    [[nodiscard]] std::optional<NastranShellView> view() const& noexcept;
    [[nodiscard]] std::optional<NastranShellView> view() const&& = delete;

private:
    NastranShellCards(NastranCoordinateCards&&, core::PlannedBuffer<NastranShellElement>&&,
                      core::PlannedBuffer<NastranShellProperty>&&, std::uint64_t) noexcept;

    NastranCoordinateCards input_;
    core::PlannedBuffer<NastranShellElement> elements_;
    core::PlannedBuffer<NastranShellProperty> properties_;
    std::uint64_t work_;
};

} // namespace uvlas::io::detail
