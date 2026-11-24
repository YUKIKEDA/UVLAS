#include "nastran_beam_cards.hpp"

#include <algorithm>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
using Reason = NastranBeamFailure;
using Position = NastranBulkPosition;
using Lexeme = NastranBulkLexeme;
constexpr auto operation = core::Operation::PrepareModel;

struct Abort {
    core::Failure failure;
    NastranBeamDiagnostic diagnostic;
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

constexpr char upper(char c) noexcept {
    return c >= 'a' && c <= 'z' ? static_cast<char>(c - ('a' - 'A')) : c;
}

enum class BeamCard : std::uint8_t { Other, Beam, Property };

BeamCard classify(std::string_view name) noexcept {
    const auto equals = [name](std::string_view expected) {
        return std::ranges::equal(name, expected, [](char a, char b) { return upper(a) == b; });
    };
    if (equals("CBEAM"))
        return BeamCard::Beam;
    if (equals("PBEAM"))
        return BeamCard::Property;
    return BeamCard::Other;
}

class Fields {
public:
    Fields(const NastranBulkView& bulk, const NastranBulkCard& card, Work& work) noexcept
        : bulk_(bulk), card_(card), slots_(bulk.data_fields(card)), work_(work) {}

    [[nodiscard]] std::optional<Lexeme> at(std::size_t i) const noexcept {
        return i < slots_.size() ? std::optional{slots_[i]} : std::nullopt;
    }

    [[nodiscard]] std::optional<Lexeme> supplied(std::size_t i) const noexcept {
        const auto field = at(i);
        return field && field->size != 0 ? field : std::nullopt;
    }

    [[nodiscard]] Position position(std::size_t i) const noexcept {
        const auto field = at(i);
        return field ? field->position : card_.name.position;
    }

    [[nodiscard]] Lexeme required(std::size_t i) const {
        const auto field = supplied(i);
        if (!field)
            reject(Code::InvalidInput, Reason::MissingField, position(i));
        return *field;
    }

    [[nodiscard]] core::Result<std::int64_t> integer(Lexeme field) const {
        work_.charge(field.size + 1, field.position);
        return bulk_.integer(field);
    }

    [[nodiscard]] std::uint64_t checked_id(Lexeme field, std::int64_t value,
                                           std::int64_t minimum = 1) const {
        if (value < minimum || (card_.format == NastranBulkFormat::Free && field.size > 8))
            reject(Code::InvalidInput, Reason::InvalidField, field.position);
        return static_cast<std::uint64_t>(value);
    }

    [[nodiscard]] std::uint64_t id(Lexeme field, std::int64_t minimum = 1) const {
        const auto value = integer(field);
        if (!value)
            reject(value.error().code, Reason::InvalidField, field.position);
        return checked_id(field, *value, minimum);
    }

    [[nodiscard]] std::optional<NastranShellReal> real(std::size_t i) const {
        const auto field = supplied(i);
        if (!field)
            return std::nullopt;
        work_.charge(field->size + 1, field->position);
        const auto result = bulk_.real(*field);
        if (!result)
            reject(result.error().code, Reason::InvalidField, field->position);
        return NastranShellReal{*field, *result};
    }

    [[nodiscard]] NastranBeamVector vector(std::size_t start) const {
        return {real(start), real(start + 1), real(start + 2)};
    }

    [[nodiscard]] std::string_view text(Lexeme field) const noexcept {
        return bulk_.text(field);
    }

    void require_blank(std::size_t i) const {
        if (const auto field = supplied(i))
            reject(Code::InvalidInput, Reason::ExtraField, field->position);
    }

    void check_tail(std::size_t first) const {
        for (const auto& field : slots_.subspan(std::min(first, slots_.size()))) {
            work_.charge(1, field.position);
            if (field.size != 0)
                reject(Code::InvalidInput, Reason::ExtraField, field.position);
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
};

template <class Row>
void index_rows(std::span<Key> keys, std::span<const Row> rows, Work& work) {
    for (std::size_t i = 0; i < rows.size(); ++i) {
        work.charge(1, rows[i].id_field.position);
        keys[i] = {rows[i].id.value(), i};
    }
    std::ranges::sort(keys, [&work](const Key& a, const Key& b) {
        work.charge(1);
        return a.id < b.id || (a.id == b.id && a.ordinal < b.ordinal);
    });
}

void check_properties(std::span<const Key> keys, std::span<const NastranPbeamDeclaration> rows,
                      Work& work) {
    std::optional<std::size_t> duplicate;
    for (std::size_t i = 1; i < keys.size(); ++i) {
        const auto& current = keys[i];
        work.charge(1, rows[current.ordinal].id_field.position);
        if (current.id == keys[i - 1].id && (!duplicate || current.ordinal < *duplicate))
            duplicate = current.ordinal;
    }
    if (duplicate)
        reject(Code::InvalidInput, Reason::DuplicateProperty, rows[*duplicate].id_field.position);
}

std::optional<std::size_t> lookup(std::span<const Key> keys, std::uint64_t id, Position position,
                                  Work& work) {
    const auto found = std::ranges::lower_bound(
        keys, id,
        [&work, position](std::uint64_t a, std::uint64_t b) {
            work.charge(1, position);
            return a < b;
        },
        &Key::id);
    if (found == keys.end() || found->id != id)
        return std::nullopt;
    return found->ordinal;
}

class References {
public:
    References(std::span<const Key> grid_keys, std::span<const Key> property_keys,
               std::span<const NastranShellProperty> shell_properties, Work& work) noexcept
        : grids_(grid_keys),
          properties_(property_keys),
          shell_properties_(shell_properties),
          work_(work) {}

    [[nodiscard]] NastranBeamGridReference grid(std::uint64_t id, Lexeme field) const {
        const auto ordinal = lookup(grids_, id, field.position, work_);
        if (!ordinal)
            reject(Code::InvalidInput, Reason::UndefinedGrid, field.position);
        return {NastranGridId{id}, field, *ordinal};
    }

    [[nodiscard]] NastranBeamWarpingPoint warping(std::uint64_t id, Lexeme field) const {
        if (id == 0)
            return {NastranGridId{id}, field, std::nullopt};
        const auto ordinal = lookup(grids_, id, field.position, work_);
        if (!ordinal)
            reject(Code::InvalidInput, Reason::UndefinedWarpingPoint, field.position);
        return {NastranGridId{id}, field, ordinal};
    }

    [[nodiscard]] NastranPbeamReference property(std::uint64_t id, std::optional<Lexeme> field,
                                                 Position position) const {
        if (const auto ordinal = lookup(properties_, id, position, work_))
            return {NastranPropertyId{id}, field, *ordinal};
        for (const auto& other : shell_properties_) {
            work_.charge(1, position);
            if (other.id.value() == id)
                reject(Code::InvalidInput, Reason::WrongPropertyKind, position);
        }
        reject(Code::InvalidInput, Reason::UndefinedProperty, position);
    }

private:
    std::span<const Key> grids_, properties_;
    std::span<const NastranShellProperty> shell_properties_;
    Work& work_;
};

NastranBeamReleases releases(const Fields& fields, std::size_t slot) {
    NastranBeamReleases result{0, fields.at(slot)};
    if (!result.field)
        return result;
    const auto text = fields.text(*result.field);
    if (text.size() > 5)
        reject(Code::InvalidInput, Reason::InvalidField, result.field->position);
    for (const char digit : text) {
        if (digit < '1' || digit > '6')
            reject(Code::InvalidInput, Reason::InvalidField, result.field->position);
        const auto flag = static_cast<std::uint8_t>(1U << static_cast<unsigned>(digit - '1'));
        if ((result.element_components & flag) != 0)
            reject(Code::InvalidInput, Reason::InvalidField, result.field->position);
        result.element_components = static_cast<std::uint8_t>(result.element_components | flag);
    }
    return result;
}

NastranBeamInterpretation interpretation(const Fields& fields) {
    constexpr std::size_t slot = 7;
    NastranBeamInterpretation result;
    result.field = fields.at(slot);
    const auto supplied = fields.supplied(slot);
    if (!supplied)
        return result;
    const auto option = fields.text(*supplied);
    if (option.size() != 3 || (upper(option[0]) != 'G' && upper(option[0]) != 'B') ||
        (upper(option[1]) != 'G' && upper(option[1]) != 'O') ||
        (upper(option[2]) != 'G' && upper(option[2]) != 'O'))
        reject(Code::Unsupported, Reason::UnsupportedField, supplied->position);
    if (upper(option[0]) == 'B')
        result.orientation = NastranBeamOrientationBasis::Basic;
    if (upper(option[1]) == 'O')
        result.a = NastranBeamOffsetBasis::OffsetSystem;
    if (upper(option[2]) == 'O')
        result.b = NastranBeamOffsetBasis::OffsetSystem;
    return result;
}

NastranBeamOrientation orientation(const Fields& fields, NastranGridId a, NastranGridId b,
                                   const References& references) {
    constexpr std::size_t first = 4;
    if (const auto field = fields.supplied(first)) {
        if (const auto number = fields.integer(*field)) {
            const auto id = fields.checked_id(*field, *number);
            if (id == a.value() || id == b.value())
                reject(Code::InvalidInput, Reason::OrientationEndpoint, field->position);
            fields.require_blank(first + 1);
            fields.require_blank(first + 2);
            return references.grid(id, *field);
        }
    }
    const auto vector = fields.vector(first);
    if (std::ranges::none_of(vector, [](const auto& field) { return field.has_value(); }))
        reject(Code::InvalidInput, Reason::MissingField, fields.position(first));
    return vector;
}

void read_end_fields(NastranBeamEnd& end, const Fields& fields, std::size_t release_slot,
                     std::size_t offset_slot, std::size_t warping_slot,
                     const References& references) {
    end.releases = releases(fields, release_slot);
    end.offset = fields.vector(offset_slot);
    if (const auto field = fields.supplied(warping_slot))
        end.warping = references.warping(fields.id(*field, 0), *field);
}

NastranCbeam read_cbeam(const Fields& fields, std::size_t card, const References& references) {
    NastranCbeam beam;
    beam.card_ordinal = card;
    beam.id_field = fields.required(0);
    beam.id = NastranElementId{fields.id(beam.id_field)};
    const auto pid = fields.supplied(1);
    beam.property = references.property(pid ? fields.id(*pid) : beam.id.value(), fields.at(1),
                                        fields.position(1));
    const auto ga = fields.required(2);
    const auto gb = fields.required(3);
    const auto a = fields.id(ga);
    const auto b = fields.id(gb);
    if (a == b)
        reject(Code::InvalidInput, Reason::DuplicateEndpoint, gb.position);
    beam.a.point = references.grid(a, ga);
    beam.b.point = references.grid(b, gb);
    beam.orientation = orientation(fields, beam.a.point.id, beam.b.point.id, references);
    beam.interpretation = interpretation(fields);
    read_end_fields(beam.a, fields, 8, 10, 16, references);
    read_end_fields(beam.b, fields, 9, 13, 17, references);
    fields.check_tail(18);
    return beam;
}

template <class T>
void admit(std::size_t count, std::uint64_t& remaining, Work& work) {
    if (!std::in_range<core::Index>(count) || count > remaining / sizeof(T))
        reject(Code::ResourceExceeded, Reason::Limit);
    remaining -= count * sizeof(T);
    work.charge(count);
}

template <class T>
auto allocate(core::ResourcePlan& plan, std::size_t count) {
    auto buffer = plan.allocate<T>(static_cast<core::Index>(count), operation);
    if (!buffer)
        throw Abort{buffer.error(), {Reason::Limit, {}}};
    return std::move(*buffer);
}

struct Counts {
    std::size_t beams = 0, properties = 0;
};

Counts count_cards(const NastranBulkView& bulk, NastranBeamLimits limits, Work& work) {
    Counts counts;
    for (const auto& card : bulk.cards) {
        work.charge(1, card.name.position);
        const auto kind = classify(bulk.text(card.name));
        if (kind == BeamCard::Beam) {
            if (counts.beams == limits.beams)
                reject(Code::ResourceExceeded, Reason::Limit, card.name.position);
            ++counts.beams;
        } else if (kind == BeamCard::Property) {
            if (counts.properties == limits.properties)
                reject(Code::ResourceExceeded, Reason::Limit, card.name.position);
            ++counts.properties;
        }
    }
    return counts;
}

void read_properties(const NastranBulkView& bulk, std::span<NastranPbeamDeclaration> rows,
                     Work& work) {
    std::size_t row = 0;
    for (std::size_t i = 0; i < bulk.cards.size(); ++i) {
        const auto& card = bulk.cards[i];
        work.charge(1, card.name.position);
        if (classify(bulk.text(card.name)) == BeamCard::Property) {
            const Fields fields{bulk, card, work};
            const auto id = fields.required(0);
            rows[row++] = {NastranPropertyId{fields.id(id)}, id, i};
        }
    }
}

} // namespace

NastranBeamCards::NastranBeamCards(NastranBeamOrigin source,
                                   core::PlannedBuffer<NastranCbeam>&& beams,
                                   core::PlannedBuffer<NastranPbeamDeclaration>&& properties,
                                   std::uint64_t work) noexcept
    : source_(source), beams_(std::move(beams)), properties_(std::move(properties)), work_(work) {}

NastranBeamCards::NastranBeamCards(NastranBeamCards&& other) noexcept
    : source_(std::exchange(other.source_, std::nullopt)),
      beams_(std::move(other.beams_)),
      properties_(std::move(other.properties_)),
      work_(other.work_) {}

std::optional<NastranBeamView> NastranBeamCards::view() const& noexcept {
    if (!source_)
        return std::nullopt;
    return NastranBeamView{*source_, beams_.values(), properties_.values(), work_};
}

core::Result<NastranBeamCards>
NastranBeamCards::prepare(const NastranShellCards& parent, NastranBeamLimits limits,
                          core::ResourcePlan& plan, const std::stop_token& stop,
                          std::optional<NastranBeamDiagnostic>* diagnostic) noexcept {
    if (diagnostic)
        diagnostic->reset();
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<NastranBeamCards> {
            try {
                Work work{limits.work, stop};
                work.charge(0);
                const auto source = parent.view();
                if (!source)
                    reject(Code::InvalidInput, Reason::Input);
                const auto& bulk = source->input.bulk;
                const auto counts = count_cards(bulk, limits, work);
                const auto grids =
                    counts.beams == 0 ? source->input.grids.first(0) : source->input.grids;
                auto remaining = plan.remaining_bytes();
                admit<NastranCbeam>(counts.beams, remaining, work);
                admit<NastranPbeamDeclaration>(counts.properties, remaining, work);
                admit<Key>(grids.size(), remaining, work);
                admit<Key>(counts.properties, remaining, work);

                auto beams = allocate<NastranCbeam>(plan, counts.beams);
                auto properties = allocate<NastranPbeamDeclaration>(plan, counts.properties);
                auto grid_keys = allocate<Key>(plan, grids.size());
                auto property_keys = allocate<Key>(plan, counts.properties);
                read_properties(bulk, properties.values(), work);
                const auto property_rows = std::as_const(properties).values();
                index_rows(grid_keys.values(), grids, work);
                index_rows(property_keys.values(), property_rows, work);
                check_properties(property_keys.values(), property_rows, work);
                const References references{grid_keys.values(), property_keys.values(),
                                            source->properties, work};
                std::size_t beam = 0;
                for (std::size_t i = 0; i < bulk.cards.size(); ++i) {
                    const auto& card = bulk.cards[i];
                    work.charge(1, card.name.position);
                    if (classify(bulk.text(card.name)) == BeamCard::Beam)
                        beams.values()[beam++] =
                            read_cbeam(Fields{bulk, card, work}, i, references);
                }
                work.charge(0);
                return NastranBeamCards{{bulk.location, bulk.digest},
                                        std::move(beams),
                                        std::move(properties),
                                        work.used()};
            } catch (const Abort& error) {
                if (diagnostic && error.failure.code != Code::Cancelled)
                    *diagnostic = error.diagnostic;
                return std::unexpected(error.failure);
            }
        });
}

} // namespace uvlas::io::detail
