#include "nastran_coordinate_cards.hpp"

#include <algorithm>
#include <limits>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
using Reason = NastranCoordinateFailure;
using Position = NastranBulkPosition;
using Lexeme = NastranBulkLexeme;
constexpr auto operation = core::Operation::PrepareModel;
constexpr auto no_index = std::numeric_limits<std::size_t>::max();
constexpr std::int64_t grid_id_exclusive_limit = 100000000;

struct Abort {
    core::Failure failure;
    NastranCoordinateDiagnostic diagnostic;
};

[[noreturn]] void reject(Code code, Reason reason, Position position = {}) {
    throw Abort{{code, operation}, {reason, position}};
}

class Work {
public:
    Work(std::uint64_t limit, const std::stop_token& stop) noexcept : limit_(limit), stop_(stop) {}

    void charge(std::uint64_t amount, Position position = {}) {
        if (stop_.stop_requested())
            reject(Code::Cancelled, Reason::Limit, position);
        if (amount > limit_ - used_)
            reject(Code::ResourceExceeded, Reason::Limit, position);
        used_ += amount;
    }

    [[nodiscard]] std::uint64_t used() const noexcept {
        return used_;
    }

private:
    std::uint64_t limit_;
    const std::stop_token& stop_;
    std::uint64_t used_ = 0;
};

enum class Card : std::uint8_t { Grid, Rectangular, Cylindrical, Spherical, Deferred };

bool keyword(std::string_view actual, std::string_view expected) noexcept {
    return std::ranges::equal(actual, expected, [](char a, char b) {
        return (a >= 'a' && a <= 'z' ? a - ('a' - 'A') : a) == b;
    });
}

Card classify(std::string_view name, Position position) {
    constexpr std::array cards{
        std::pair{"GRID", Card::Grid}, std::pair{"CORD2R", Card::Rectangular},
        std::pair{"CORD2C", Card::Cylindrical}, std::pair{"CORD2S", Card::Spherical}};
    for (const auto& [label, card] : cards)
        if (keyword(name, label))
            return card;
    // These adopted entries neither define a coordinate system nor change GRID defaults.
    constexpr std::array<std::string_view, 13> deferred{
        "CQUAD4", "CTRIA3", "MAT1", "PSHELL", "MAT8", "PCOMP", "CBEAM",
        "PBEAM",  "CONM2",  "RBE2", "SPC",    "SPC1", "SPCADD"};
    if (std::ranges::any_of(deferred, [name](auto label) { return keyword(name, label); }))
        return Card::Deferred;
    reject(Code::Unsupported, Reason::UnsupportedCard, position);
}

class Fields {
public:
    Fields(const NastranBulkView& bulk, const NastranBulkCard& card, Work& work) noexcept
        : bulk_(bulk),
          values_(bulk.data_fields(card)),
          name_(card.name),
          format_(card.format),
          work_(work) {}

    [[nodiscard]] std::optional<Lexeme> field(std::size_t index) const noexcept {
        if (index >= values_.size())
            return std::nullopt;
        return values_[index];
    }

    [[nodiscard]] Lexeme required(std::size_t index) const {
        const auto slot = field(index);
        if (!slot || slot->size == 0)
            reject(Code::InvalidInput, Reason::MissingField,
                   slot ? slot->position : name_.position);
        return *slot;
    }

    [[nodiscard]] std::int64_t integer(Lexeme field) const {
        work_.charge(field.size + 1, field.position);
        if (format_ == NastranBulkFormat::Free && field.size > 8)
            reject(Code::InvalidInput, Reason::InvalidField, field.position);
        const auto value = bulk_.integer(field);
        if (!value)
            reject(value.error().code, Reason::InvalidField, field.position);
        return *value;
    }

    [[nodiscard]] std::int64_t integer_or_basic(std::size_t index) const {
        const auto slot = field(index);
        return slot && slot->size != 0 ? integer(*slot) : 0;
    }

    [[nodiscard]] NastranCoordinateReference reference(std::size_t index) const {
        const auto slot = field(index);
        const auto id = integer_or_basic(index);
        if (id < 0)
            reject(Code::InvalidInput, Reason::InvalidField,
                   slot ? slot->position : name_.position);
        return {NastranCoordinateId{static_cast<std::uint64_t>(id)}, slot};
    }

    [[nodiscard]] NastranCoordinateValue real(std::size_t index, bool mandatory) const {
        const auto source = mandatory ? std::optional{required(index)} : field(index);
        if (!source || source->size == 0)
            return {std::nullopt, source};
        work_.charge(source->size + 1, source->position);
        const auto value = bulk_.real(*source);
        if (!value)
            reject(value.error().code, Reason::InvalidField, source->position);
        return {*value, source};
    }

    [[nodiscard]] std::array<NastranCoordinateValue, 3> triple(std::size_t index,
                                                               bool mandatory) const {
        return {real(index, mandatory), real(index + 1, mandatory), real(index + 2, mandatory)};
    }

    [[nodiscard]] std::uint8_t permanent(std::size_t index) const {
        const auto slot = field(index);
        if (!slot || slot->size == 0)
            return 0;
        const auto spelling = bulk_.text(*slot);
        work_.charge(slot->size + 1, slot->position);
        if (spelling == "0")
            return 0;
        std::uint8_t mask = 0;
        for (const char digit : spelling) {
            if (digit < '1' || digit > '6')
                reject(Code::InvalidInput, Reason::InvalidField, slot->position);
            mask = static_cast<std::uint8_t>(mask | (1U << static_cast<unsigned>(digit - '1')));
        }
        return mask;
    }

    void check_tail(std::size_t fields) const {
        for (const auto& slot : values_.subspan(std::min(fields, values_.size()))) {
            work_.charge(1, slot.position);
            if (slot.size != 0)
                reject(Code::InvalidInput, Reason::ExtraField, slot.position);
        }
    }

private:
    const NastranBulkView& bulk_;
    std::span<const Lexeme> values_;
    Lexeme name_;
    NastranBulkFormat format_;
    Work& work_;
};

NastranGridPoint grid(const NastranBulkCard& card, const Fields& fields) {
    const auto id_token = fields.required(0);
    const auto id = fields.integer(id_token);
    if (id <= 0 || id >= grid_id_exclusive_limit)
        reject(Code::InvalidInput, Reason::InvalidField, id_token.position);
    const auto cp = fields.reference(1);
    const auto position = fields.triple(2, false);
    const auto cd_token = fields.field(5);
    const auto cd = fields.integer_or_basic(5);
    if (cd < 0) {
        const auto where = cd_token ? cd_token->position : card.name.position;
        reject(cd == -1 ? Code::Unsupported : Code::InvalidInput,
               cd == -1 ? Reason::UnsupportedField : Reason::InvalidField, where);
    }
    const auto permanent = fields.permanent(6);
    const auto sid_token = fields.field(7);
    const auto superelement = fields.integer_or_basic(7);
    if (superelement != 0) {
        const auto where = sid_token ? sid_token->position : card.name.position;
        reject(superelement > 0 ? Code::Unsupported : Code::InvalidInput,
               superelement > 0 ? Reason::UnsupportedField : Reason::InvalidField, where);
    }
    fields.check_tail(8);
    const NastranCoordinateReference displacement{
        NastranCoordinateId{static_cast<std::uint64_t>(cd)}, cd_token};
    return {NastranGridId{static_cast<std::uint64_t>(id)},
            card.name,
            id_token,
            cp,
            displacement,
            position,
            permanent,
            fields.field(6),
            sid_token};
}

NastranCoordinateDefinition coordinate(const NastranBulkCard& card, Card kind,
                                       const Fields& fields) {
    const auto id_token = fields.required(0);
    const auto id = fields.integer(id_token);
    if (id <= 0)
        reject(Code::InvalidInput, Reason::InvalidField, id_token.position);
    auto chart = NastranCoordinateKind::Rectangular;
    if (kind == Card::Cylindrical)
        chart = NastranCoordinateKind::Cylindrical;
    if (kind == Card::Spherical)
        chart = NastranCoordinateKind::Spherical;
    const auto rid = fields.reference(1);
    const auto a = fields.triple(2, true);
    const auto b = fields.triple(5, true);
    const auto c = fields.triple(8, true);
    fields.check_tail(11);
    return {NastranCoordinateId{static_cast<std::uint64_t>(id)},
            chart,
            card.name,
            id_token,
            rid,
            a,
            b,
            c};
}

template <class T>
void admit(std::size_t count, std::uint64_t& remaining, Work& work) {
    if (!std::in_range<core::Index>(count) || count > remaining / sizeof(T))
        reject(Code::ResourceExceeded, Reason::Limit);
    work.charge(count);
    remaining -= count * sizeof(T);
}

template <class T>
auto allocate(core::ResourcePlan& plan, std::size_t count) {
    auto result = plan.allocate<T>(static_cast<core::Index>(count), operation);
    if (!result)
        throw Abort{result.error(), {Reason::Limit, {}}};
    return std::move(*result);
}

struct IndexEntry {
    std::uint64_t id = 0;
    std::size_t index = 0;
};

template <class Row>
void build_index(std::span<IndexEntry> index, std::span<const Row> rows, Reason duplicate,
                 Work& work) {
    for (std::size_t i = 0; i < index.size(); ++i) {
        work.charge(1, rows[i].id_field.position);
        index[i] = {rows[i].id.value(), i};
    }
    std::ranges::sort(index, [&](const IndexEntry& a, const IndexEntry& b) {
        work.charge(1, rows[a.index].id_field.position);
        return a.id < b.id || (a.id == b.id && a.index < b.index);
    });
    std::optional<std::size_t> first_duplicate;
    for (std::size_t i = 1; i < index.size(); ++i) {
        const auto& entry = index[i];
        work.charge(1, rows[entry.index].id_field.position);
        if (index[i - 1].id == entry.id && (!first_duplicate || entry.index < *first_duplicate))
            first_duplicate = entry.index;
    }
    if (first_duplicate)
        reject(Code::InvalidInput, duplicate, rows[*first_duplicate].id_field.position);
}

Position reference_position(const NastranCoordinateReference& reference, Lexeme owner) noexcept {
    return reference.field ? reference.field->position : owner.position;
}

std::size_t resolve(std::span<const IndexEntry> index, const NastranCoordinateReference& reference,
                    Lexeme owner, Work& work) {
    if (reference.id.value() == 0)
        return no_index;
    const auto position = reference_position(reference, owner);
    const auto found = std::ranges::lower_bound(
        index, reference.id.value(),
        [&](std::uint64_t a, std::uint64_t b) {
            work.charge(1, position);
            return a < b;
        },
        &IndexEntry::id);
    if (found == index.end() || found->id != reference.id.value())
        reject(Code::InvalidInput, Reason::UndefinedCoordinate, position);
    return found->index;
}

struct GraphNode {
    std::size_t parent = no_index;
    std::size_t depth = 0;
    bool active = false;
};

void dependency_order(std::span<const NastranCoordinateDefinition> coordinates,
                      std::span<GraphNode> graph, std::span<std::size_t> stack,
                      std::span<std::size_t> order, std::uint64_t limit, Work& work) {
    std::size_t finished = 0;
    for (std::size_t i = 0; i < graph.size(); ++i) {
        work.charge(1, coordinates[i].name.position);
        if (graph[i].depth != 0)
            continue;
        std::size_t current = i;
        std::size_t used = 0;
        std::size_t depth = 0;
        auto where = coordinates[i].name.position;
        while (current != no_index) {
            work.charge(1, where);
            auto& node = graph[current];
            if (node.active)
                reject(Code::InvalidInput, Reason::CoordinateCycle, where);
            if (used >= limit || node.depth > limit - used)
                reject(Code::ResourceExceeded, Reason::Limit, where);
            if (node.depth != 0) {
                depth = node.depth;
                break;
            }
            node.active = true;
            stack[used++] = current;
            where = reference_position(coordinates[current].reference, coordinates[current].name);
            current = node.parent;
        }
        while (used != 0) {
            const auto row = stack[--used];
            work.charge(1, coordinates[row].name.position);
            graph[row].depth = ++depth;
            graph[row].active = false;
            order[finished++] = row;
        }
    }
}

struct Counts {
    std::size_t grids = 0;
    std::size_t coordinates = 0;
    std::size_t deferred = 0;
};

} // namespace

NastranCoordinateCards::NastranCoordinateCards(
    NastranBulkDeck&& bulk, core::PlannedBuffer<NastranGridPoint>&& grids,
    core::PlannedBuffer<NastranCoordinateDefinition>&& coordinates,
    core::PlannedBuffer<std::size_t>&& order, std::size_t deferred_cards,
    std::uint64_t work) noexcept
    : bulk_(std::move(bulk)),
      grids_(std::move(grids)),
      coordinates_(std::move(coordinates)),
      dependency_first_(std::move(order)),
      deferred_cards_(deferred_cards),
      work_(work) {}

std::optional<NastranCoordinateView> NastranCoordinateCards::view() const& noexcept {
    const auto bulk = bulk_.view();
    if (!bulk)
        return std::nullopt;
    return NastranCoordinateView{
        *bulk, grids_.values(), coordinates_.values(), dependency_first_.values(), deferred_cards_,
        work_};
}

core::Result<NastranCoordinateCards>
NastranCoordinateCards::prepare(NastranBulkDeck&& deck, NastranCoordinateLimits limits,
                                core::ResourcePlan& plan, const std::stop_token& stop,
                                std::optional<NastranCoordinateDiagnostic>* diagnostic) noexcept {
    if (diagnostic)
        diagnostic->reset();
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<NastranCoordinateCards> {
            try {
                Work work{limits.work, stop};
                work.charge(0);
                const auto bulk = deck.view();
                if (!bulk)
                    reject(Code::InvalidInput, Reason::Input);
                Counts counts;
                for (const auto& entry : bulk->cards) {
                    work.charge(1, entry.name.position);
                    const auto kind = classify(bulk->text(entry.name), entry.name.position);
                    if (kind == Card::Grid) {
                        if (counts.grids == limits.grids)
                            reject(Code::ResourceExceeded, Reason::Limit, entry.name.position);
                        ++counts.grids;
                    } else if (kind == Card::Deferred) {
                        ++counts.deferred;
                    } else {
                        if (counts.coordinates == limits.coordinate_systems)
                            reject(Code::ResourceExceeded, Reason::Limit, entry.name.position);
                        ++counts.coordinates;
                    }
                }

                auto remaining = plan.remaining_bytes();
                const auto stack_size = static_cast<std::size_t>(
                    std::min<std::uint64_t>(limits.depth, counts.coordinates));
                admit<NastranGridPoint>(counts.grids, remaining, work);
                admit<NastranCoordinateDefinition>(counts.coordinates, remaining, work);
                admit<IndexEntry>(counts.grids, remaining, work);
                admit<IndexEntry>(counts.coordinates, remaining, work);
                admit<GraphNode>(counts.coordinates, remaining, work);
                admit<std::size_t>(stack_size, remaining, work);
                admit<std::size_t>(counts.coordinates, remaining, work);

                auto grids = allocate<NastranGridPoint>(plan, counts.grids);
                auto coordinates = allocate<NastranCoordinateDefinition>(plan, counts.coordinates);
                auto grid_index = allocate<IndexEntry>(plan, counts.grids);
                auto coordinate_index = allocate<IndexEntry>(plan, counts.coordinates);
                auto graph = allocate<GraphNode>(plan, counts.coordinates);
                auto stack = allocate<std::size_t>(plan, stack_size);
                auto order = allocate<std::size_t>(plan, counts.coordinates);
                std::size_t grid_slot = 0;
                std::size_t coordinate_slot = 0;
                for (const auto& entry : bulk->cards) {
                    work.charge(1, entry.name.position);
                    const auto kind = classify(bulk->text(entry.name), entry.name.position);
                    const Fields fields{*bulk, entry, work};
                    if (kind == Card::Grid)
                        grids.values()[grid_slot++] = grid(entry, fields);
                    else if (kind != Card::Deferred)
                        coordinates.values()[coordinate_slot++] = coordinate(entry, kind, fields);
                }

                const auto grid_rows = std::as_const(grids).values();
                const auto coordinate_rows = std::as_const(coordinates).values();
                build_index(grid_index.values(), grid_rows, Reason::DuplicateGrid, work);
                build_index(coordinate_index.values(), coordinate_rows, Reason::DuplicateCoordinate,
                            work);
                for (std::size_t i = 0; i < coordinate_rows.size(); ++i) {
                    const auto& row = coordinate_rows[i];
                    work.charge(1, row.name.position);
                    graph.values()[i].parent =
                        resolve(coordinate_index.values(), row.reference, row.name, work);
                }
                for (const auto& row : grid_rows) {
                    work.charge(1, row.name.position);
                    static_cast<void>(
                        resolve(coordinate_index.values(), row.position_system, row.name, work));
                    static_cast<void>(resolve(coordinate_index.values(), row.displacement_system,
                                              row.name, work));
                }

                dependency_order(coordinate_rows, graph.values(), stack.values(), order.values(),
                                 limits.depth, work);
                work.charge(0);
                return NastranCoordinateCards{std::move(deck),        std::move(grids),
                                              std::move(coordinates), std::move(order),
                                              counts.deferred,        work.used()};
            } catch (const Abort& error) {
                if (diagnostic && error.failure.code != Code::Cancelled)
                    *diagnostic = error.diagnostic;
                return std::unexpected(error.failure);
            }
        });
}

} // namespace uvlas::io::detail
