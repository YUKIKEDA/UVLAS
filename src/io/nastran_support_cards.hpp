#pragma once

#include "nastran_coordinate_cards.hpp"

#include <variant>

namespace uvlas::io::detail {

struct NastranSupportSetTag;
using NastranSupportSetId = core::Identifier<NastranSupportSetTag>;

struct NastranSupportPoint {
    NastranGridId grid{0};
    NastranBulkLexeme field;
    std::size_t grid_ordinal = 0;
};

struct NastranSupportRange {
    NastranGridId first{0}, last{0};
    NastranBulkLexeme first_field, thru_field, last_field;
};

using NastranSupportTarget = std::variant<NastranSupportPoint, NastranSupportRange>;

struct NastranSupportTerm {
    NastranSupportSetId set{0};
    NastranSupportTarget target;
    std::uint8_t components = 0; // External GRID components 1 through 6.
    NastranBulkLexeme card, set_field, component_field;
    std::optional<NastranBulkLexeme> motion_field;
    std::optional<NastranRealField> prescribed_motion;
};

struct NastranSupportUnionMember {
    NastranSupportSetId combination{0}, member{0};
    NastranBulkLexeme card, combination_field, member_field;
};

enum class NastranSupportFailure : std::uint8_t {
    Input,
    MissingField,
    InvalidField,
    ExtraField,
    UnknownGrid,
    UnknownSet,
    DuplicateCombination,
    DuplicateMember,
    NestedCombination,
    Limit
};

struct NastranSupportDiagnostic {
    NastranSupportFailure reason = NastranSupportFailure::Input;
    NastranBulkPosition position;
};

struct NastranSupportLimits {
    std::uint64_t terms;
    std::uint64_t union_members;
    std::uint64_t work;
};

struct NastranSupportOrigin {
    RootSnapshotLocation location;
    SourceDigest digest;
};

struct NastranSupportView {
    NastranSupportView(NastranSupportOrigin owner, std::span<const NastranSupportTerm> term_records,
                       std::span<const NastranSupportUnionMember> member_records,
                       std::uint64_t executed_work) noexcept
        : source(owner), terms(term_records), union_members(member_records), work(executed_work) {}

    NastranSupportOrigin source;
    std::span<const NastranSupportTerm> terms;
    std::span<const NastranSupportUnionMember> union_members;
    std::uint64_t work;
};

// Owns the parsed records. The authenticated deck can still move on to another
// card reader. Its location and digest identify the source of each field.
class NastranSupportCards {
public:
    NastranSupportCards(NastranSupportCards&&) noexcept;
    NastranSupportCards& operator=(NastranSupportCards&&) = delete;

    [[nodiscard]] static core::Result<NastranSupportCards>
    prepare(const NastranCoordinateCards&, NastranSupportLimits, core::ResourcePlan&,
            const std::stop_token&, std::optional<NastranSupportDiagnostic>* = nullptr) noexcept;

    [[nodiscard]] std::optional<NastranSupportView> view() const& noexcept;
    [[nodiscard]] std::optional<NastranSupportView> view() const&& = delete;

private:
    NastranSupportCards(NastranSupportOrigin, core::PlannedBuffer<NastranSupportTerm>&&,
                        core::PlannedBuffer<NastranSupportUnionMember>&&, std::uint64_t) noexcept;

    std::optional<NastranSupportOrigin> source_;
    core::PlannedBuffer<NastranSupportTerm> terms_;
    core::PlannedBuffer<NastranSupportUnionMember> members_;
    std::uint64_t work_;
};

} // namespace uvlas::io::detail
