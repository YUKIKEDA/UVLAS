#include "nastran_rigid_link_cards.hpp"

#include <algorithm>
#include <limits>
#include <new>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
using Reason = NastranRigidLinkFailure;
using Lexeme = NastranBulkLexeme;
using Position = NastranBulkPosition;
constexpr auto operation = core::Operation::PrepareModel;

struct Abort {
    core::Failure failure;
    NastranRigidLinkDiagnostic diagnostic;
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

[[nodiscard]] bool is_rbe2(std::string_view text) noexcept {
    return std::ranges::equal(text, std::string_view{"RBE2"}, [](char actual, char expected) {
        return (actual >= 'a' && actual <= 'z' ? actual - ('a' - 'A') : actual) == expected;
    });
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

void count(std::size_t& current, std::size_t increment, std::uint64_t limit, Position position) {
    if (increment > limit - current ||
        increment > std::numeric_limits<std::size_t>::max() - current)
        reject(Reason::Limit, position, Code::ResourceExceeded);
    current += increment;
}

struct CardFields {
    const NastranBulkView& bulk;
    const NastranBulkCard& card;
    std::span<const Lexeme> fields;

    [[nodiscard]] Lexeme required(std::size_t index) const {
        if (index >= fields.size())
            reject(Reason::MissingField, card.name.position);
        const auto field = fields[index];
        if (field.size == 0)
            reject(Reason::MissingField, field.position);
        return field;
    }

    [[nodiscard]] std::uint64_t positive(Lexeme field) const {
        const auto value = bulk.integer(field);
        if (!value || *value <= 0 || (card.format == NastranBulkFormat::Free && field.size > 8))
            reject(Reason::InvalidField, field.position);
        return static_cast<std::uint64_t>(*value);
    }

    [[nodiscard]] std::uint8_t components(Lexeme field) const {
        std::uint8_t mask = 0;
        for (const char digit : bulk.text(field)) {
            if (digit < '1' || digit > '6')
                reject(Reason::InvalidField, field.position);
            const auto bit = static_cast<std::uint8_t>(1U << static_cast<unsigned>(digit - '1'));
            if ((mask & bit) != 0)
                reject(Reason::InvalidField, field.position);
            mask = static_cast<std::uint8_t>(mask | bit);
        }
        return mask;
    }

    [[nodiscard]] NastranRigidThermalField thermal(Lexeme field) const {
        if (field.size == 0)
            reject(Reason::UnsupportedThermalForm, field.position, Code::Unsupported);
        const auto value = bulk.real(field);
        if (!value)
            reject(Reason::InvalidField, field.position);
        return {*value, field};
    }
};

struct ParsedLink {
    NastranRigidLink link;
    std::span<const Lexeme> dependents;
};

[[nodiscard]] ParsedLink parse(const CardFields& input, Work& work) {
    const auto eid = input.required(0);
    const auto independent = input.required(1);
    const auto components = input.required(2);
    static_cast<void>(input.required(3));
    std::size_t end = input.fields.size();
    while (input.fields[end - 1].size == 0) {
        work.visit(input.fields[end - 1].position);
        --end;
    }
    std::size_t stop = 3;
    for (; stop < end; ++stop) {
        const auto field = input.fields[stop];
        work.visit(field.position);
        if (!input.bulk.integer(field))
            break;
        static_cast<void>(input.positive(field));
    }
    if (stop == 3)
        reject(Reason::MissingField, input.fields[3].position);
    if (end - stop > 2)
        reject(Reason::ExtraField, input.fields[stop + 2].position);
    ParsedLink result{{NastranElementId{input.positive(eid)},
                       {NastranGridId{input.positive(independent)}, independent, 0},
                       input.components(components),
                       input.card.name,
                       eid,
                       components,
                       0,
                       stop - 3,
                       {},
                       {}},
                      input.fields.subspan(3, stop - 3)};
    if (stop < end)
        result.link.alpha = input.thermal(input.fields[stop]);
    if (end - stop == 2) {
        work.visit(input.fields[stop + 1].position);
        result.link.reference_temperature = input.thermal(input.fields[stop + 1]);
    }
    return result;
}

struct IndexedGrid {
    NastranGridId id{0};
    std::size_t ordinal = 0;
};

[[nodiscard]] auto grid_index(std::span<const NastranGridPoint> grids, Work& work,
                              core::ResourcePlan& memory) {
    auto result = allocate<IndexedGrid>(memory, grids.size());
    const auto index = result.values();
    for (std::size_t i = 0; i < grids.size(); ++i) {
        work.visit(grids[i].id_field.position);
        index[i] = {grids[i].id, i};
    }
    std::ranges::sort(
        index,
        [&](NastranGridId lhs, NastranGridId rhs) {
            work.visit();
            return lhs < rhs;
        },
        &IndexedGrid::id);
    return result;
}

[[nodiscard]] std::size_t find_grid(std::span<const IndexedGrid> index,
                                    NastranRigidGridReference point, Work& work) {
    const auto found = std::ranges::lower_bound(
        index, point.id,
        [&](NastranGridId lhs, NastranGridId rhs) {
            work.visit(point.field.position);
            return lhs < rhs;
        },
        &IndexedGrid::id);
    if (found == index.end() || found->id != point.id)
        reject(Reason::UnknownGrid, point.field.position);
    return found->ordinal;
}

void resolve(ParsedLink& parsed, const CardFields& fields, std::span<const IndexedGrid> grids,
             std::span<NastranRigidGridReference> output, std::span<IndexedGrid> order,
             Work& work) {
    auto& link = parsed.link;
    link.independent.ordinal = find_grid(grids, link.independent, work);
    for (std::size_t i = 0; i < parsed.dependents.size(); ++i) {
        const auto field = parsed.dependents[i];
        work.visit(field.position);
        const NastranGridId id{fields.positive(field)};
        if (id == link.independent.id)
            reject(Reason::RepeatedGrid, field.position);
        NastranRigidGridReference point{id, field, 0};
        point.ordinal = find_grid(grids, point, work);
        output[i] = point;
        order[i] = {id, i};
    }
    std::ranges::sort(order, [&](const IndexedGrid& lhs, const IndexedGrid& rhs) {
        work.visit(link.card.position);
        return std::pair{lhs.id, lhs.ordinal} < std::pair{rhs.id, rhs.ordinal};
    });
    for (std::size_t i = 1; i < order.size(); ++i) {
        const auto field = output[order[i].ordinal].field;
        work.visit(field.position);
        if (order[i - 1].id == order[i].id)
            reject(Reason::RepeatedGrid, field.position);
    }
}

} // namespace

NastranRigidLinkCards::NastranRigidLinkCards(
    NastranRigidLinkSource source, core::PlannedBuffer<NastranRigidLink>&& links,
    core::PlannedBuffer<NastranRigidGridReference>&& dependents, std::uint64_t work) noexcept
    : source_(source), links_(std::move(links)), dependents_(std::move(dependents)), work_(work) {}

NastranRigidLinkCards::NastranRigidLinkCards(NastranRigidLinkCards&& other) noexcept
    : source_(std::exchange(other.source_, std::nullopt)),
      links_(std::move(other.links_)),
      dependents_(std::move(other.dependents_)),
      work_(other.work_) {}

core::Result<NastranRigidLinkCards>
NastranRigidLinkCards::prepare(const NastranShellCards& shells, NastranRigidLinkLimits limits,
                               core::ResourcePlan& memory, const std::stop_token& stop,
                               std::optional<NastranRigidLinkDiagnostic>* diagnostic) noexcept {
    if (diagnostic != nullptr)
        diagnostic->reset();
    try {
        Work work{limits.work, stop};
        work.visit();
        const auto input = shells.view();
        if (!input)
            reject(Reason::Input);
        const auto& bulk = input->input.bulk;
        std::size_t link_count = 0, total = 0, largest = 0;
        for (const auto& card : bulk.cards) {
            work.visit(card.name.position);
            if (!is_rbe2(bulk.text(card.name)))
                continue;
            const auto parsed = parse({bulk, card, bulk.data_fields(card)}, work);
            count(link_count, 1, limits.links, card.name.position);
            count(total, parsed.dependents.size(), limits.dependent_grids, card.name.position);
            largest = std::max(largest, parsed.dependents.size());
        }
        auto links = allocate<NastranRigidLink>(memory, link_count);
        auto dependents = allocate<NastranRigidGridReference>(memory, total);
        if (link_count != 0) {
            const auto grids = grid_index(input->input.grids, work, memory);
            auto order = allocate<IndexedGrid>(memory, largest);
            std::size_t link_slot = 0, member_slot = 0;
            for (const auto& card : bulk.cards) {
                work.visit(card.name.position);
                if (!is_rbe2(bulk.text(card.name)))
                    continue;
                const CardFields fields{bulk, card, bulk.data_fields(card)};
                auto parsed = parse(fields, work);
                const auto size = parsed.dependents.size();
                resolve(parsed, fields, grids.values(),
                        dependents.values().subspan(member_slot, size), order.values().first(size),
                        work);
                parsed.link.first_dependent = member_slot;
                member_slot += size;
                links.values()[link_slot++] = parsed.link;
            }
        }
        return NastranRigidLinkCards{
            {bulk.location, bulk.digest}, std::move(links), std::move(dependents), work.used};
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

std::optional<NastranRigidLinkView> NastranRigidLinkCards::view() const& noexcept {
    if (!source_)
        return std::nullopt;
    return NastranRigidLinkView{*source_, links_.values(), dependents_.values(), work_};
}

} // namespace uvlas::io::detail
