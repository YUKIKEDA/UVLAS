#include "nastran_shell_cards.hpp"

#include <algorithm>
#include <limits>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
using Reason = NastranShellFailure;
using Position = NastranBulkPosition;
using Lexeme = NastranBulkLexeme;
constexpr auto operation = core::Operation::PrepareModel;
constexpr std::uint64_t element_id_limit = 100000000;
constexpr std::uint64_t pcomp_id_limit = 10000000;
constexpr std::size_t tflag_slot = 9;

struct Abort {
    core::Failure failure;
    NastranShellDiagnostic diagnostic;
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

enum class Card : std::uint8_t { Other, Tri, Quad, OtherElement, Pshell, Pcomp, Pbeam };

Card classify(std::string_view name) noexcept {
    constexpr std::array cards{
        std::pair{"CTRIA3", Card::Tri},         std::pair{"CQUAD4", Card::Quad},
        std::pair{"CBEAM", Card::OtherElement}, std::pair{"CONM2", Card::OtherElement},
        std::pair{"RBE2", Card::OtherElement},  std::pair{"PSHELL", Card::Pshell},
        std::pair{"PCOMP", Card::Pcomp},        std::pair{"PBEAM", Card::Pbeam}};
    for (const auto& [label, card] : cards)
        if (std::ranges::equal(name, std::string_view{label}, [](char a, char b) {
                return (a >= 'a' && a <= 'z' ? a - ('a' - 'A') : a) == b;
            }))
            return card;
    return Card::Other;
}

class Fields {
public:
    Fields(const NastranBulkView& bulk, const NastranBulkCard& card, Work& work) noexcept
        : bulk_(bulk), card_(card), slots_(bulk.data_fields(card)), work_(work) {}

    [[nodiscard]] std::optional<Lexeme> field(std::size_t i) const noexcept {
        return i < slots_.size() ? std::optional{slots_[i]} : std::nullopt;
    }

    [[nodiscard]] Lexeme required(std::size_t i) const {
        const auto slot = field(i);
        if (!slot || slot->size == 0)
            reject(Code::InvalidInput, Reason::MissingField,
                   slot ? slot->position : card_.name.position);
        return *slot;
    }

    [[nodiscard]] core::Result<std::int64_t> integer(Lexeme field) const {
        work_.charge(field.size + 1, field.position);
        return bulk_.integer(field);
    }

    [[nodiscard]] std::uint64_t checked_integer(Lexeme field, std::int64_t value,
                                                std::uint64_t minimum,
                                                std::uint64_t exclusive_max) const {
        if ((card_.format == NastranBulkFormat::Free && field.size > 8) || value < 0 ||
            static_cast<std::uint64_t>(value) < minimum ||
            static_cast<std::uint64_t>(value) >= exclusive_max)
            reject(Code::InvalidInput, Reason::InvalidField, field.position);
        return static_cast<std::uint64_t>(value);
    }

    [[nodiscard]] std::uint64_t
    natural(Lexeme field, std::uint64_t minimum = 1,
            std::uint64_t exclusive_max = std::numeric_limits<std::uint64_t>::max()) const {
        const auto value = integer(field);
        if (!value)
            reject(value.error().code, Reason::InvalidField, field.position);
        return checked_integer(field, *value, minimum, exclusive_max);
    }

    [[nodiscard]] NastranShellReal real(Lexeme field) const {
        work_.charge(field.size + 1, field.position);
        const auto value = bulk_.real(field);
        if (!value)
            reject(value.error().code, Reason::InvalidField, field.position);
        return {field, *value};
    }

    [[nodiscard]] std::optional<NastranShellReal> optional_real(std::size_t i) const {
        const auto slot = field(i);
        if (!slot || slot->size == 0)
            return std::nullopt;
        return real(*slot);
    }

    void require_blanks(std::size_t first, std::size_t end) const {
        for (std::size_t i = first; i < std::min(end, slots_.size()); ++i) {
            const auto& slot = slots_[i];
            work_.charge(1, slot.position);
            if (slot.size != 0)
                reject(Code::InvalidInput, Reason::ExtraField, slot.position);
        }
    }

private:
    const NastranBulkView& bulk_;
    const NastranBulkCard& card_;
    std::span<const Lexeme> slots_;
    Work& work_;
};

struct Key {
    std::uint64_t id = 0;
    std::size_t ordinal = 0;
    Lexeme source;
};

void order(std::span<Key> keys, Work& work) {
    std::ranges::sort(keys, [&work](const Key& a, const Key& b) {
        work.charge(1, a.source.position);
        return a.id < b.id || (a.id == b.id && a.ordinal < b.ordinal);
    });
}

void unique(std::span<const Key> keys, Reason reason, Work& work) {
    std::optional<Position> duplicate;
    for (std::size_t i = 1; i < keys.size(); ++i) {
        const auto pos = keys[i].source.position;
        work.charge(1, pos);
        if (keys[i].id == keys[i - 1].id && (!duplicate || pos.byte < duplicate->byte))
            duplicate = pos;
    }
    if (duplicate)
        reject(Code::InvalidInput, reason, *duplicate);
}

std::optional<std::size_t> lookup(std::span<const Key> keys, std::uint64_t id, Position pos,
                                  Work& work) {
    const auto found = std::ranges::lower_bound(
        keys, id,
        [&work, pos](std::uint64_t a, std::uint64_t b) {
            work.charge(1, pos);
            return a < b;
        },
        &Key::id);
    if (found != keys.end() && found->id == id)
        return found->ordinal;
    return std::nullopt;
}

template <class Record>
void index_checked_rows(std::span<Key> keys, std::span<const Record> rows, Work& work) {
    for (std::size_t i = 0; i < rows.size(); ++i) {
        work.charge(1, rows[i].id_field.position);
        keys[i] = {rows[i].id.value(), i, rows[i].id_field};
    }
    order(keys, work);
}

template <class T>
void admit(std::size_t count, std::uint64_t& remaining, Work& work) {
    if (!std::in_range<core::Index>(count) || count > remaining / sizeof(T))
        reject(Code::ResourceExceeded, Reason::Limit);
    remaining -= count * sizeof(T);
    work.charge(count);
}

template <class T>
core::PlannedBuffer<T> allocate(core::ResourcePlan& plan, std::size_t count) {
    auto memory = plan.allocate<T>(static_cast<core::Index>(count), operation);
    if (!memory)
        throw Abort{memory.error(), {Reason::Limit, {}}};
    return std::move(*memory);
}

struct Counts {
    std::size_t elements = 0, shells = 0, properties = 0, beams = 0;
};

Counts count_cards(const NastranBulkView& bulk, NastranShellLimits limits, Work& work) {
    Counts counts;
    for (const auto& card : bulk.cards) {
        const auto pos = card.name.position;
        work.charge(1, pos);
        const auto kind = classify(bulk.text(card.name));
        if (kind == Card::Tri || kind == Card::Quad || kind == Card::OtherElement) {
            if (counts.elements == limits.element_ids)
                reject(Code::ResourceExceeded, Reason::Limit, pos);
            ++counts.elements;
            if (kind != Card::OtherElement)
                ++counts.shells;
        } else if (kind == Card::Pshell || kind == Card::Pcomp || kind == Card::Pbeam) {
            if (counts.properties + counts.beams == limits.property_ids)
                reject(Code::ResourceExceeded, Reason::Limit, pos);
            if (kind == Card::Pbeam)
                ++counts.beams;
            else
                ++counts.properties;
        }
    }
    return counts;
}

struct Indexes {
    std::span<Key> grids, coordinates, elements, properties, beam_properties;
};

void read_declarations(const NastranBulkView& bulk, Indexes keys,
                       std::span<NastranShellElement> elements,
                       std::span<NastranShellProperty> properties, Work& work) {
    std::size_t element = 0, shell = 0, property = 0, beam = 0;
    for (std::size_t i = 0; i < bulk.cards.size(); ++i) {
        const auto& card = bulk.cards[i];
        work.charge(1, card.name.position);
        const auto kind = classify(bulk.text(card.name));
        if (kind == Card::Other)
            continue;
        const Fields fields{bulk, card, work};
        const auto field = fields.required(0);
        if (kind == Card::Tri || kind == Card::Quad || kind == Card::OtherElement) {
            const auto id = fields.natural(field, 1, element_id_limit);
            keys.elements[element++] = {id, i, field};
            if (kind != Card::OtherElement) {
                auto& row = elements[shell++];
                row.id = NastranElementId{id};
                row.id_field = field;
                row.card_ordinal = i;
            }
        } else if (kind == Card::Pbeam) {
            keys.beam_properties[beam++] = {fields.natural(field), i, field};
        } else {
            const auto id = kind == Card::Pcomp ? fields.natural(field, 1, pcomp_id_limit)
                                                : fields.natural(field);
            properties[property] = {NastranPropertyId{id},
                                    kind == Card::Pcomp ? NastranShellPropertyKind::Pcomp
                                                        : NastranShellPropertyKind::Pshell,
                                    field, i};
            keys.properties[property] = {id, property, field};
            ++property;
        }
    }
    order(keys.elements, work);
    unique(keys.elements, Reason::DuplicateElement, work);
    order(keys.properties, work);
    unique(keys.properties, Reason::DuplicateProperty, work);
    order(keys.beam_properties, work);
}

NastranShellPropertyReference property_reference(const NastranShellElement& element,
                                                 const Fields& fields, Indexes keys, Work& work) {
    const auto field = fields.field(1);
    const auto id = field && field->size != 0 ? fields.natural(*field) : element.id.value();
    const auto pos = field ? field->position : element.id_field.position;
    if (const auto slot = lookup(keys.properties, id, pos, work))
        return {NastranPropertyId{id}, field, *slot};
    const auto reason = lookup(keys.beam_properties, id, pos, work) ? Reason::WrongPropertyKind
                                                                    : Reason::UndefinedProperty;
    reject(Code::InvalidInput, reason, pos);
}

NastranShellOrientation orientation(const Fields& fields, std::size_t slot,
                                    std::span<const Key> coordinates, Work& work) {
    const auto field = fields.field(slot);
    if (!field || field->size == 0)
        return std::monostate{};
    if (const auto integer = fields.integer(*field)) {
        const auto id =
            fields.checked_integer(*field, *integer, 0, std::numeric_limits<std::uint64_t>::max());
        std::optional<std::size_t> coordinate;
        if (id != 0) {
            coordinate = lookup(coordinates, id, field->position, work);
            if (!coordinate)
                reject(Code::InvalidInput, Reason::UndefinedCoordinate, field->position);
        }
        return NastranShellMaterialAxis{NastranCoordinateId{id}, *field, coordinate};
    }
    return fields.real(*field);
}

template <std::size_t N>
std::array<NastranShellCorner, N> corners(const Fields& fields, std::span<const Key> grids,
                                          NastranShellThicknessMode mode, Work& work) {
    std::array<NastranShellCorner, N> result;
    std::size_t zeroes = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto field = fields.required(2 + i);
        const auto id = fields.natural(field);
        const auto prior = std::span<const NastranShellCorner>{result}.first(i);
        if (std::ranges::any_of(prior,
                                [id](const auto& corner) { return corner.id.value() == id; }))
            reject(Code::InvalidInput, Reason::DuplicateGrid, field.position);
        const auto grid = lookup(grids, id, field.position, work);
        if (!grid)
            reject(Code::InvalidInput, Reason::UndefinedGrid, field.position);
        result[i] = {NastranGridId{id}, field, *grid, fields.optional_real(tflag_slot + 1 + i)};
        if (const auto& thickness = result[i].thickness) {
            const auto decimal = thickness->value.decimal();
            const auto mantissa = decimal.substr(0, decimal.find('e'));
            const bool zero = mantissa.find_first_of("123456789") == std::string_view::npos;
            if ((zero && mode == NastranShellThicknessMode::Fraction) ||
                (!zero && mantissa.starts_with('-')))
                reject(Code::InvalidInput, Reason::InvalidField, thickness->field.position);
            if (zero)
                ++zeroes;
            if (zeroes == N)
                reject(Code::InvalidInput, Reason::InvalidField, thickness->field.position);
        }
    }
    return result;
}

template <std::size_t N>
void read_connectivity(NastranShellElement& row, const Fields& fields, Indexes keys, Work& work) {
    constexpr std::size_t orientation_slot = 2 + N;
    row.property = property_reference(row, fields, keys, work);
    row.orientation = orientation(fields, orientation_slot, keys.coordinates, work);
    row.offset = fields.optional_real(orientation_slot + 1);
    fields.require_blanks(orientation_slot + 2, tflag_slot);
    row.thickness_mode_field = fields.field(tflag_slot);
    if (const auto field = row.thickness_mode_field; field && field->size != 0) {
        const auto flag = fields.natural(*field, 0, 2);
        if (flag == 1)
            row.thickness_mode = NastranShellThicknessMode::Fraction;
    }
    row.corners = corners<N>(fields, keys.grids, row.thickness_mode, work);
    fields.require_blanks(tflag_slot + 1 + N, std::numeric_limits<std::size_t>::max());
}

} // namespace

NastranShellCards::NastranShellCards(NastranCoordinateCards&& input,
                                     core::PlannedBuffer<NastranShellElement>&& elements,
                                     core::PlannedBuffer<NastranShellProperty>&& properties,
                                     std::uint64_t work) noexcept
    : input_(std::move(input)),
      elements_(std::move(elements)),
      properties_(std::move(properties)),
      work_(work) {}

std::optional<NastranShellView> NastranShellCards::view() const& noexcept {
    const auto input = input_.view();
    if (!input)
        return std::nullopt;
    return NastranShellView{*input, elements_.values(), properties_.values(), work_};
}

core::Result<NastranShellCards>
NastranShellCards::prepare(NastranCoordinateCards&& source, NastranShellLimits limits,
                           core::ResourcePlan& plan, const std::stop_token& stop,
                           std::optional<NastranShellDiagnostic>* diagnostic) noexcept {
    if (diagnostic)
        diagnostic->reset();
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<NastranShellCards> {
            try {
                Work work{limits.work, stop};
                work.charge(0);
                const auto input = source.view();
                if (!input)
                    reject(Code::InvalidInput, Reason::Input);
                const auto counts = count_cards(input->bulk, limits, work);
                auto remaining = plan.remaining_bytes();
                admit<NastranShellElement>(counts.shells, remaining, work);
                admit<NastranShellProperty>(counts.properties, remaining, work);
                admit<Key>(input->grids.size(), remaining, work);
                admit<Key>(input->coordinates.size(), remaining, work);
                admit<Key>(counts.elements, remaining, work);
                admit<Key>(counts.properties, remaining, work);
                admit<Key>(counts.beams, remaining, work);

                auto elements = allocate<NastranShellElement>(plan, counts.shells);
                auto properties = allocate<NastranShellProperty>(plan, counts.properties);
                auto grid_index = allocate<Key>(plan, input->grids.size());
                auto coord_index = allocate<Key>(plan, input->coordinates.size());
                auto element_index = allocate<Key>(plan, counts.elements);
                auto property_index = allocate<Key>(plan, counts.properties);
                auto beam_index = allocate<Key>(plan, counts.beams);
                const Indexes keys{grid_index.values(), coord_index.values(),
                                   element_index.values(), property_index.values(),
                                   beam_index.values()};
                index_checked_rows(keys.grids, input->grids, work);
                index_checked_rows(keys.coordinates, input->coordinates, work);
                read_declarations(input->bulk, keys, elements.values(), properties.values(), work);
                for (auto& row : elements.values()) {
                    const auto& card = input->bulk.cards[row.card_ordinal];
                    work.charge(1, card.name.position);
                    const Fields fields{input->bulk, card, work};
                    if (classify(input->bulk.text(card.name)) == Card::Tri)
                        read_connectivity<3>(row, fields, keys, work);
                    else
                        read_connectivity<4>(row, fields, keys, work);
                }
                work.charge(0);
                return NastranShellCards{std::move(source), std::move(elements),
                                         std::move(properties), work.used()};
            } catch (const Abort& error) {
                if (diagnostic && error.failure.code != Code::Cancelled)
                    *diagnostic = error.diagnostic;
                return std::unexpected(error.failure);
            }
        });
}

} // namespace uvlas::io::detail
