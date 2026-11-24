#pragma once

#include "nastran_shell_cards.hpp"

#include <variant>

namespace uvlas::io::detail {

struct NastranBasicMassOffset {};
struct NastranBasicMassCentre {};
struct NastranFramedMassOffset {
    NastranCoordinateId frame{0};
    std::size_t frame_ordinal = 0;
};

using NastranMassFrame =
    std::variant<NastranBasicMassOffset, NastranBasicMassCentre, NastranFramedMassOffset>;

struct NastranMassReal {
    std::optional<NastranRealField> supplied;
    std::optional<NastranBulkLexeme> field;
};

struct NastranMassInertia {
    NastranMassReal i11, i21, i22, i31, i32, i33;
};

struct NastranLumpedMass {
    NastranElementId id{0};
    NastranGridId grid{0};
    std::size_t grid_ordinal = 0;
    NastranMassFrame frame;
    NastranBulkLexeme card, id_field, grid_field;
    std::optional<NastranBulkLexeme> frame_field;
    NastranMassReal mass;
    std::array<NastranMassReal, 3> position;
    NastranMassInertia inertia; // Input products retain their original signs about the centre.
};

enum class NastranLumpedMassFailure : std::uint8_t {
    Input,
    MissingField,
    InvalidField,
    ExtraField,
    UnknownGrid,
    UnknownFrame,
    Limit
};

struct NastranLumpedMassDiagnostic {
    NastranLumpedMassFailure reason;
    NastranBulkPosition position;
};

struct NastranLumpedMassLimits {
    std::uint64_t records;
    std::uint64_t work;
};

struct NastranLumpedMassSource {
    RootSnapshotLocation location;
    SourceDigest digest;
};

struct NastranLumpedMassView {
    NastranLumpedMassSource source;
    std::span<const NastranLumpedMass> masses;
    std::uint64_t work;
};

class NastranLumpedMassCards {
public:
    NastranLumpedMassCards(NastranLumpedMassCards&&) noexcept;
    NastranLumpedMassCards& operator=(NastranLumpedMassCards&&) = delete;

    // Reuse the checked element ID namespace without consuming the deck or transforming the values.
    [[nodiscard]] static core::Result<NastranLumpedMassCards>
    prepare(const NastranShellCards&, NastranLumpedMassLimits, core::ResourcePlan&,
            const std::stop_token&, std::optional<NastranLumpedMassDiagnostic>* = nullptr) noexcept;

    [[nodiscard]] std::optional<NastranLumpedMassView> view() const& noexcept;
    [[nodiscard]] std::optional<NastranLumpedMassView> view() const&& = delete;

private:
    NastranLumpedMassCards(NastranLumpedMassSource, core::PlannedBuffer<NastranLumpedMass>&&,
                           std::uint64_t) noexcept;

    std::optional<NastranLumpedMassSource> source_;
    core::PlannedBuffer<NastranLumpedMass> masses_;
    std::uint64_t work_;
};

} // namespace uvlas::io::detail
