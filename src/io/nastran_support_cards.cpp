#include "nastran_support_cards.hpp"

#include <algorithm>
#include <limits>
#include <new>
#include <tuple>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
using Reason = NastranSupportFailure;
using Lexeme = NastranBulkLexeme;
using Position = NastranBulkPosition;
constexpr auto operation = core::Operation::PrepareModel;

struct Abort {
    core::Failure failure;
    NastranSupportDiagnostic diagnostic;
};

[[noreturn]] void reject(Reason reason, Position position = {}, Code code = Code::InvalidInput) {
    throw Abort{{code, operation}, {reason, position}};
}

struct Work {
    std::uint64_t limit;
    const std::stop_token& stop;
    std::uint64_t used = 0;

    void visit(Position position = {}) {
        if (stop.stop_requested())
            reject(Reason::Limit, position, Code::Cancelled);
        if (used == limit)
            reject(Reason::Limit, position, Code::ResourceExceeded);
        ++used;
    }
};

[[nodiscard]] bool same_name(std::string_view text, std::string_view name) noexcept {
    return std::ranges::equal(text, name, [](char actual, char expected) {
        return (actual >= 'a' && actual <= 'z' ? actual - ('a' - 'A') : actual) == expected;
    });
}

struct Fields {
    const NastranBulkView& bulk;
    std::span<const Lexeme> fields;
    std::size_t used;
    Lexeme card;
    NastranBulkFormat format;

    [[nodiscard]] std::optional<Lexeme> at(std::size_t index) const noexcept {
        return index < fields.size() ? std::optional{fields[index]} : std::nullopt;
    }

    [[nodiscard]] Lexeme required(std::size_t index) const {
        const auto field = at(index);
        if (!field || field->size == 0)
            reject(Reason::MissingField, field ? field->position : card.position);
        return *field;
    }

    [[nodiscard]] std::uint64_t positive(Lexeme field) const {
        if (format == NastranBulkFormat::Free && field.size > 8)
            reject(Reason::InvalidField, field.position);
        const auto value = bulk.integer(field);
        if (!value || *value < 1)
            reject(Reason::InvalidField, field.position);
        return static_cast<std::uint64_t>(*value);
    }

    [[nodiscard]] std::uint8_t components(Lexeme field) const {
        const auto text = bulk.text(field);
        if (text == "0")
            reject(Reason::InvalidField, field.position, Code::Unsupported);
        std::uint8_t mask = 0;
        for (const char digit : text) {
            if (digit < '1' || digit > '6')
                reject(Reason::InvalidField, field.position);
            const auto bit = static_cast<std::uint8_t>(1U << static_cast<unsigned>(digit - '1'));
            if ((mask & bit) != 0)
                reject(Reason::InvalidField, field.position);
            mask = static_cast<std::uint8_t>(mask | bit);
        }
        return mask;
    }
};

template <class Emit>
void read_spc(const Fields& input, NastranSupportSetId set, Lexeme set_field, Work& work,
              Emit emit) {
    if (input.used > 7)
        reject(Reason::ExtraField, input.fields[7].position);
    for (const std::size_t index : {std::size_t{1}, std::size_t{4}}) {
        work.visit(input.card.position);
        if (index == 4 && input.used <= index)
            break;
        const auto grid = input.required(index);
        const auto flags = input.required(index + 1);
        const auto motion = input.at(index + 2);
        std::optional<NastranRealField> prescribed;
        if (motion && motion->size != 0) {
            const auto value = input.bulk.real(*motion);
            if (!value)
                reject(Reason::InvalidField, motion->position);
            prescribed = *value;
        }
        emit(NastranSupportTerm{
            set, NastranSupportPoint{NastranGridId{input.positive(grid)}, grid, 0},
            input.components(flags), input.card, set_field, flags, motion, prescribed});
    }
}

template <class Emit>
void read_spc1(const Fields& input, NastranSupportSetId set, Lexeme set_field, Work& work,
               Emit emit) {
    const auto flags = input.required(1);
    const auto components = input.components(flags);
    const auto first = input.required(2);
    const auto thru = input.at(3);
    if (thru && same_name(input.bulk.text(*thru), "THRU")) {
        work.visit(thru->position);
        const auto last = input.required(4);
        if (input.used > 5)
            reject(Reason::ExtraField, input.fields[5].position);
        const NastranGridId low{input.positive(first)};
        const NastranGridId high{input.positive(last)};
        if (high <= low)
            reject(Reason::InvalidField, last.position);
        emit(NastranSupportTerm{set,
                                NastranSupportRange{low, high, first, *thru, last},
                                components,
                                input.card,
                                set_field,
                                flags,
                                {},
                                {}});
        return;
    }

    for (const auto grid : input.fields.first(input.used).subspan(2)) {
        work.visit(grid.position);
        emit(NastranSupportTerm{set,
                                NastranSupportPoint{NastranGridId{input.positive(grid)}, grid, 0},
                                components,
                                input.card,
                                set_field,
                                flags,
                                {},
                                {}});
    }
}

template <class Term, class Member>
void read_cards(const NastranBulkView& bulk, Work& work, Term term, Member member) {
    for (const auto& card : bulk.cards) {
        work.visit(card.name.position);
        const auto name = bulk.text(card.name);
        const bool spc = same_name(name, "SPC");
        const bool spc1 = same_name(name, "SPC1");
        if (!spc && !spc1 && !same_name(name, "SPCADD"))
            continue;
        const auto data = bulk.data_fields(card);
        std::size_t used = data.size();
        while (used != 0 && data[used - 1].size == 0) {
            work.visit(data[used - 1].position);
            --used;
        }
        const Fields input{bulk, data, used, card.name, card.format};
        const auto sid = input.required(0);
        const NastranSupportSetId set{input.positive(sid)};
        if (spc) {
            read_spc(input, set, sid, work, term);
        } else if (spc1) {
            read_spc1(input, set, sid, work, term);
        } else {
            static_cast<void>(input.required(1));
            for (const auto field : data.first(used).subspan(1)) {
                work.visit(field.position);
                member(NastranSupportUnionMember{set, NastranSupportSetId{input.positive(field)},
                                                 card.name, sid, field});
            }
        }
    }
}

template <class T>
[[nodiscard]] auto allocate(core::ResourcePlan& memory, std::size_t size) {
    if (!std::in_range<core::Index>(size))
        reject(Reason::Limit, {}, Code::ResourceExceeded);
    auto result = memory.allocate<T>(static_cast<core::Index>(size), operation);
    if (!result)
        throw Abort{result.error(), {Reason::Limit, {}}};
    return std::move(*result);
}

void count(std::size_t& size, std::uint64_t maximum, Position position) {
    if (size >= maximum || size == std::numeric_limits<std::size_t>::max())
        reject(Reason::Limit, position, Code::ResourceExceeded);
    ++size;
}

struct GridIndex {
    NastranGridId id{0};
    std::size_t ordinal = 0;
};

void resolve_points(std::span<NastranSupportTerm> terms, std::span<const NastranGridPoint> grids,
                    Work& work, core::ResourcePlan& memory) {
    if (terms.empty())
        return;
    auto index = allocate<GridIndex>(memory, grids.size());
    auto ordered = index.values();
    for (std::size_t i = 0; i < grids.size(); ++i) {
        work.visit(grids[i].id_field.position);
        ordered[i] = {grids[i].id, i};
    }
    const auto compare = [&](NastranGridId lhs, NastranGridId rhs) {
        work.visit();
        return lhs < rhs;
    };
    std::ranges::sort(ordered, compare, &GridIndex::id);
    for (auto& term : terms) {
        work.visit(term.card.position);
        auto* point = std::get_if<NastranSupportPoint>(&term.target);
        if (point == nullptr)
            continue; // THRU covers existing grids; its endpoints need not be declared.
        const auto found = std::ranges::lower_bound(ordered, point->grid, compare, &GridIndex::id);
        if (found == ordered.end() || found->id != point->grid)
            reject(Reason::UnknownGrid, point->field.position);
        point->grid_ordinal = found->ordinal;
    }
}

struct MemberIndex {
    NastranSupportSetId combination{0}, member{0};
    std::size_t ordinal = 0;
};

void resolve_unions(std::span<const NastranSupportTerm> terms,
                    std::span<const NastranSupportUnionMember> members, Work& work,
                    core::ResourcePlan& memory) {
    if (members.empty())
        return;
    auto member_index = allocate<MemberIndex>(memory, members.size());
    const auto ordered = member_index.values();
    auto direct_index = allocate<std::uint64_t>(memory, terms.size());
    const auto direct = direct_index.values();
    for (std::size_t i = 0; i < terms.size(); ++i) {
        work.visit(terms[i].card.position);
        direct[i] = terms[i].set.value();
    }
    std::ranges::sort(direct, [&](const auto& lhs, const auto& rhs) {
        work.visit();
        return lhs < rhs;
    });
    for (std::size_t i = 0; i < members.size(); ++i) {
        work.visit(members[i].member_field.position);
        ordered[i] = {members[i].combination, members[i].member, i};
    }
    std::ranges::sort(ordered, [&](const MemberIndex& lhs, const MemberIndex& rhs) {
        work.visit();
        return std::tuple{lhs.combination, members[lhs.ordinal].card.position.byte, lhs.member,
                          lhs.ordinal} < std::tuple{rhs.combination,
                                                    members[rhs.ordinal].card.position.byte,
                                                    rhs.member, rhs.ordinal};
    });
    for (std::size_t i = 1; i < ordered.size(); ++i) {
        const auto& previous = ordered[i - 1];
        const auto& current = ordered[i];
        const auto& prior_card = members[previous.ordinal];
        const auto& card = members[current.ordinal];
        work.visit(card.member_field.position);
        if (previous.combination != current.combination)
            continue;
        if (prior_card.card.position != card.card.position)
            reject(Reason::DuplicateCombination, card.combination_field.position);
        if (previous.member == current.member)
            reject(Reason::DuplicateMember, card.member_field.position);
    }
    for (const auto& card : members) {
        work.visit(card.member_field.position);
        const auto compare = [&](const auto& lhs, const auto& rhs) {
            work.visit(card.member_field.position);
            return lhs < rhs;
        };
        if (std::ranges::binary_search(ordered, card.member, compare, &MemberIndex::combination))
            reject(Reason::NestedCombination, card.member_field.position);
        if (!std::ranges::binary_search(direct, card.member.value(), compare))
            reject(Reason::UnknownSet, card.member_field.position);
    }
}

} // namespace

NastranSupportCards::NastranSupportCards(NastranSupportOrigin source,
                                         core::PlannedBuffer<NastranSupportTerm>&& terms,
                                         core::PlannedBuffer<NastranSupportUnionMember>&& members,
                                         std::uint64_t work) noexcept
    : source_(source), terms_(std::move(terms)), members_(std::move(members)), work_(work) {}

NastranSupportCards::NastranSupportCards(NastranSupportCards&& other) noexcept
    : source_(std::exchange(other.source_, std::nullopt)),
      terms_(std::move(other.terms_)),
      members_(std::move(other.members_)),
      work_(other.work_) {}

core::Result<NastranSupportCards>
NastranSupportCards::prepare(const NastranCoordinateCards& coordinates, NastranSupportLimits limits,
                             core::ResourcePlan& memory, const std::stop_token& stop,
                             std::optional<NastranSupportDiagnostic>* diagnostic) noexcept {
    if (diagnostic != nullptr)
        diagnostic->reset();
    try {
        Work work{limits.work, stop};
        work.visit();
        const auto source = coordinates.view();
        if (!source)
            reject(Reason::Input);
        std::size_t terms = 0;
        std::size_t members = 0;
        read_cards(
            source->bulk, work,
            [&](const NastranSupportTerm& term) { count(terms, limits.terms, term.card.position); },
            [&](const NastranSupportUnionMember& member) {
                count(members, limits.union_members, member.card.position);
            });
        auto recorded_terms = allocate<NastranSupportTerm>(memory, terms);
        auto recorded_members = allocate<NastranSupportUnionMember>(memory, members);
        std::size_t term_slot = 0;
        std::size_t member_slot = 0;
        read_cards(
            source->bulk, work,
            [&](const NastranSupportTerm& term) { recorded_terms.values()[term_slot++] = term; },
            [&](const NastranSupportUnionMember& member) {
                recorded_members.values()[member_slot++] = member;
            });
        resolve_points(recorded_terms.values(), source->grids, work, memory);
        resolve_unions(recorded_terms.values(), recorded_members.values(), work, memory);
        return NastranSupportCards{{source->bulk.location, source->bulk.digest},
                                   std::move(recorded_terms),
                                   std::move(recorded_members),
                                   work.used};
    } catch (const Abort& error) {
        if (diagnostic != nullptr && error.failure.code != Code::Cancelled)
            *diagnostic = error.diagnostic;
        return std::unexpected(error.failure);
    } catch (const std::bad_alloc&) {
        return std::unexpected(core::Failure{Code::ResourceExceeded, operation});
    } catch (...) {
        return std::unexpected(core::Failure{Code::InvalidInput, operation});
    }
}

std::optional<NastranSupportView> NastranSupportCards::view() const& noexcept {
    if (!source_)
        return std::nullopt;
    return NastranSupportView{*source_, terms_.values(), members_.values(), work_};
}

} // namespace uvlas::io::detail
