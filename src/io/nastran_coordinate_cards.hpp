#pragma once

#include "nastran_bulk.hpp"

namespace uvlas::io::detail {

struct NastranGridTag;
struct NastranCoordinateTag;
using NastranGridId = core::Identifier<NastranGridTag>;
using NastranCoordinateId = core::Identifier<NastranCoordinateTag>;

enum class NastranCoordinateKind : std::uint8_t { Rectangular, Cylindrical, Spherical };

struct NastranCoordinateReference {
    NastranCoordinateId id{0}; // Zero denotes the implicit basic rectangular system.
    std::optional<NastranBulkLexeme> field;
};

struct NastranCoordinateValue {
    std::optional<NastranRealField> supplied;
    std::optional<NastranBulkLexeme> field;

    [[nodiscard]] std::string_view decimal() const& noexcept {
        return supplied ? supplied->decimal() : std::string_view{"0.0"};
    }
    [[nodiscard]] std::string_view decimal() const&& = delete;
};

struct NastranCoordinateDefinition {
    NastranCoordinateId id{0};
    NastranCoordinateKind kind = NastranCoordinateKind::Rectangular;
    NastranBulkLexeme name, id_field;
    NastranCoordinateReference reference;
    // Raw chart coordinates of A, B and C, in the RID system.
    std::array<NastranCoordinateValue, 3> a, b, c;
};

struct NastranGridPoint {
    NastranGridId id{0};
    NastranBulkLexeme name, id_field;
    NastranCoordinateReference position_system, displacement_system;
    std::array<NastranCoordinateValue, 3> position;
    // Bit zero describes external component 1. No UVLAS DOF is assigned here.
    std::uint8_t permanent_components = 0;
    std::optional<NastranBulkLexeme> permanent_field, superelement_field;
};

enum class NastranCoordinateFailure : std::uint8_t {
    Input,
    UnsupportedCard,
    UnsupportedField,
    MissingField,
    InvalidField,
    ExtraField,
    DuplicateGrid,
    DuplicateCoordinate,
    UndefinedCoordinate,
    CoordinateCycle,
    Limit
};

struct NastranCoordinateDiagnostic {
    NastranCoordinateFailure reason = NastranCoordinateFailure::Input;
    NastranBulkPosition position;
};

struct NastranCoordinateLimits {
    std::uint64_t grids;
    std::uint64_t coordinate_systems;
    std::uint64_t depth; // Counts explicit CORD records along a path to basic.
    std::uint64_t work;
};

struct NastranCoordinateView {
    NastranBulkView bulk;
    std::span<const NastranGridPoint> grids;
    std::span<const NastranCoordinateDefinition> coordinates;
    std::span<const std::size_t> dependency_first; // Ordinals in coordinates.
    std::size_t deferred_cards;
    std::uint64_t work;
};

// Owns the authenticated source, GRID/CORD2 records and their checked reference order.
// Other adopted structural cards remain visible in the bulk owner for later consumers.
class NastranCoordinateCards {
public:
    NastranCoordinateCards(NastranCoordinateCards&&) noexcept = default;
    NastranCoordinateCards& operator=(NastranCoordinateCards&&) = delete;

    [[nodiscard]] static core::Result<NastranCoordinateCards>
    prepare(NastranBulkDeck&&, NastranCoordinateLimits, core::ResourcePlan&, const std::stop_token&,
            std::optional<NastranCoordinateDiagnostic>* = nullptr) noexcept;

    [[nodiscard]] std::optional<NastranCoordinateView> view() const& noexcept;
    [[nodiscard]] std::optional<NastranCoordinateView> view() const&& = delete;

private:
    NastranCoordinateCards(NastranBulkDeck&&, core::PlannedBuffer<NastranGridPoint>&&,
                           core::PlannedBuffer<NastranCoordinateDefinition>&&,
                           core::PlannedBuffer<std::size_t>&&, std::size_t, std::uint64_t) noexcept;

    NastranBulkDeck bulk_;
    core::PlannedBuffer<NastranGridPoint> grids_;
    core::PlannedBuffer<NastranCoordinateDefinition> coordinates_;
    core::PlannedBuffer<std::size_t> dependency_first_;
    std::size_t deferred_cards_;
    std::uint64_t work_;
};

} // namespace uvlas::io::detail
