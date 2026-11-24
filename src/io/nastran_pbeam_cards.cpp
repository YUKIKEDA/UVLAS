#include "nastran_pbeam_cards.hpp"

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
using Reason = NastranPbeamFailure;
using Position = NastranBulkPosition;
using Lexeme = NastranBulkLexeme;
constexpr auto operation = core::Operation::PrepareModel;
constexpr std::size_t group_width = 8;
constexpr std::size_t maximum_stations = 10;
constexpr std::string_view nonzero = "123456789";

struct Abort {
    core::Failure failure;
    NastranPbeamDiagnostic diagnostic;
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

bool same_snapshot(const NastranBulkView& bulk, const NastranBeamOrigin& origin) noexcept {
    const auto& a = bulk.location;
    const auto& b = origin.location;
    return bulk.digest == origin.digest && a.source == b.source && a.role == b.role &&
           a.declaration_node == b.declaration_node && a.parent.source == b.parent.source &&
           a.parent.contract_revision == b.parent.contract_revision &&
           a.parent.representation_revision == b.parent.representation_revision;
}

constexpr char upper(char c) noexcept {
    return c >= 'a' && c <= 'z' ? static_cast<char>(c - ('a' - 'A')) : c;
}

class Fields {
public:
    Fields(const NastranBulkView& bulk, const NastranBulkCard& card, Work& work) noexcept
        : bulk_(bulk), card_(card), slots_(bulk.data_fields(card)), work_(work) {}

    [[nodiscard]] std::optional<Lexeme> supplied(std::size_t i) const noexcept {
        if (i < slots_.size() && slots_[i].size != 0)
            return slots_[i];
        return std::nullopt;
    }

    [[nodiscard]] Position position(std::size_t i) const noexcept {
        return i < slots_.size() ? slots_[i].position : card_.name.position;
    }

    [[nodiscard]] bool has_group(std::size_t i) const noexcept {
        return i < slots_.size();
    }

    [[nodiscard]] Lexeme required(std::size_t i) const {
        if (const auto field = supplied(i))
            return *field;
        reject(Code::InvalidInput, Reason::MissingField, position(i));
    }

    [[nodiscard]] NastranMaterialReference material() const {
        const auto field = required(1);
        work_.charge(field.size + 1, field.position);
        const auto value = bulk_.integer(field);
        if (!value)
            reject(value.error().code, Reason::InvalidField, field.position);
        if (*value <= 0 || (card_.format == NastranBulkFormat::Free && field.size > 8))
            reject(Code::InvalidInput, Reason::InvalidField, field.position);
        return {NastranMaterialId{static_cast<std::uint64_t>(*value)}, field, 0};
    }

    [[nodiscard]] NastranPbeamReal real(std::size_t i) const {
        const auto field = supplied(i);
        if (!field)
            return std::nullopt;
        work_.charge(field->size + 1, field->position);
        const auto value = bulk_.real(*field);
        if (!value)
            reject(value.error().code, Reason::InvalidField, field->position);
        return NastranShellReal{*field, *value};
    }

    [[nodiscard]] NastranShellReal required_real(std::size_t i) const {
        if (const auto value = real(i))
            return *value;
        reject(Code::InvalidInput, Reason::MissingField, position(i));
    }

    [[nodiscard]] NastranShellReal positive_real(std::size_t i) const {
        const auto value = required_real(i);
        const auto decimal = value.value.decimal();
        const auto mantissa = decimal.substr(0, decimal.find('e'));
        if (mantissa.starts_with('-') || mantissa.find_first_of(nonzero) == std::string_view::npos)
            reject(Code::InvalidInput, Reason::InvalidField, value.field.position);
        return value;
    }

    [[nodiscard]] NastranPbeamReal nonnegative_real(std::size_t i) const {
        const auto value = real(i);
        if (value) {
            const auto decimal = value->value.decimal();
            const auto mantissa = decimal.substr(0, decimal.find('e'));
            if (mantissa.starts_with('-') &&
                mantissa.find_first_of(nonzero) != std::string_view::npos)
                reject(Code::InvalidInput, Reason::InvalidField, value->field.position);
        }
        return value;
    }

    [[nodiscard]] bool station_token(std::size_t i) const {
        if (const auto value = supplied(i)) {
            const auto initial = upper(bulk_.text(*value).front());
            return initial >= 'A' && initial <= 'Z';
        }
        return false;
    }

    [[nodiscard]] NastranPbeamOutput output(Lexeme field) const {
        const auto actual = bulk_.text(field);
        work_.charge(field.size + 1, field.position);
        const auto matches = [&](std::string_view option) {
            return std::ranges::equal(actual, option, [](char a, char b) { return upper(a) == b; });
        };
        if (matches("YES"))
            return NastranPbeamOutput::SectionPoints;
        if (matches("YESA"))
            return NastranPbeamOutput::EndAPoints;
        if (matches("NO"))
            return NastranPbeamOutput::No;
        reject(Code::Unsupported, Reason::UnsupportedField, field.position);
    }

    void check_tail(std::size_t i) const {
        for (const auto& field : slots_.subspan(std::min(i, slots_.size()))) {
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

// At most 80 source bytes bound the base-10 exponent and the significand.
// The exponent is held in a fixed 512-bit integer with no allocation.
struct FractionKey {
    boost::multiprecision::int512_t power{};
    std::array<char, 82> significand{};
};

FractionKey station_key(const NastranShellReal& value) {
    const auto decimal = value.value.decimal();
    const auto marker = decimal.find('e');
    const auto mantissa = decimal.substr(0, marker);
    const auto first = mantissa.find_first_of(nonzero);
    if (mantissa.starts_with('-') || first == std::string_view::npos)
        reject(Code::InvalidInput, Reason::InvalidField, value.field.position);
    FractionKey key;
    std::size_t digit_count = 0;
    for (const char digit : mantissa.substr(first))
        if (digit != '.')
            key.significand[digit_count++] = digit;
    while (key.significand[digit_count - 1] == '0')
        key.significand[--digit_count] = '\0';

    if (marker != std::string_view::npos) {
        auto digits = decimal.substr(marker + 1);
        const bool negative = digits.starts_with('-');
        if (negative || digits.starts_with('+'))
            digits.remove_prefix(1);
        for (const char digit : digits) {
            key.power *= 10;
            key.power += digit - '0';
        }
        if (negative)
            key.power = -key.power;
    }
    const auto distance = static_cast<int>(mantissa.find('.')) - static_cast<int>(first);
    key.power += distance > 0 ? distance - 1 : distance;
    if (key.power > 0 || (key.power == 0 && (digit_count != 1 || key.significand[0] != '1')))
        reject(Code::InvalidInput, Reason::InvalidField, value.field.position);
    return key;
}

NastranPbeamSection end_a(const Fields& fields) {
    return {fields.positive_real(2), fields.positive_real(3),    fields.positive_real(4),
            fields.real(5),          fields.nonnegative_real(6), fields.real(7)};
}

NastranPbeamSection station_section(const Fields& fields, std::size_t start) {
    return {fields.real(start),     fields.real(start + 1), fields.real(start + 2),
            fields.real(start + 3), fields.real(start + 4), fields.real(start + 5)};
}

NastranPbeamRecovery recovery(const Fields& fields, std::size_t start) {
    NastranPbeamRecovery result;
    result.position = fields.position(start);
    for (std::size_t i = 0; i < result.locations.size(); ++i)
        result.locations[i] = fields.real(start + i);
    return result;
}

NastranPbeamModifiers modifiers(const Fields& fields, std::size_t start) {
    return {fields.real(start),
            fields.real(start + 1),
            fields.real(start + 2),
            fields.real(start + 3),
            {fields.real(start + 4), fields.real(start + 5)},
            {fields.real(start + 6), fields.real(start + 7)}};
}

NastranPbeamOffsets offsets(const Fields& fields, std::size_t start) {
    return {{fields.real(start), fields.real(start + 1), fields.real(start + 2),
             fields.real(start + 3)},
            {fields.real(start + 4), fields.real(start + 5), fields.real(start + 6),
             fields.real(start + 7)}};
}

template <class Emit>
NastranPbeam read_property(const Fields& fields, NastranPbeamDeclaration declaration, Work& work,
                           Emit emit) {
    NastranPbeam property;
    property.declaration = declaration;
    property.material = fields.material();
    property.end_a = end_a(fields);
    std::size_t cursor = group_width;
    if (fields.has_group(cursor) && !fields.station_token(cursor)) {
        work.charge(1, fields.position(cursor));
        property.recovery_a = recovery(fields, cursor);
        cursor += group_width;
    }

    std::array<FractionKey, maximum_stations> keys;
    while (fields.station_token(cursor)) {
        const auto output_field = fields.required(cursor);
        work.charge(1, output_field.position);
        if (property.station_count == keys.size())
            reject(Code::InvalidInput, Reason::ExtraField, output_field.position);
        const auto output = fields.output(output_field);
        const auto fraction = fields.required_real(cursor + 1);
        const auto key = station_key(fraction);
        for (std::size_t i = 0; i < property.station_count; ++i) {
            work.charge(1, fraction.field.position);
            if (keys[i].power == key.power && keys[i].significand == key.significand)
                reject(Code::InvalidInput, Reason::DuplicateStation, fraction.field.position);
        }
        keys[property.station_count] = key;
        if (key.power == 0)
            property.end_b_station = property.station_count;
        NastranPbeamStation station{output, output_field, fraction,
                                    station_section(fields, cursor + 2), std::nullopt};
        cursor += group_width;
        if (output == NastranPbeamOutput::SectionPoints) {
            work.charge(1, fields.position(cursor));
            if (!fields.has_group(cursor) || fields.station_token(cursor))
                reject(Code::InvalidInput, Reason::MissingField, fields.position(cursor));
            station.recovery = recovery(fields, cursor);
            cursor += group_width;
        }
        emit(station);
        ++property.station_count;
    }
    if (property.station_count != 0 && !property.end_b_station)
        reject(Code::InvalidInput, Reason::MissingEndB, declaration.id_field.position);
    if (fields.has_group(cursor)) {
        work.charge(1, fields.position(cursor));
        property.modifiers = modifiers(fields, cursor);
        cursor += group_width;
    }
    if (fields.has_group(cursor)) {
        work.charge(1, fields.position(cursor));
        property.offsets = offsets(fields, cursor);
        cursor += group_width;
    }
    fields.check_tail(cursor);
    return property;
}

struct MaterialKey {
    std::uint64_t id = 0;
    std::size_t ordinal = 0;
};

void index_materials(std::span<MaterialKey> keys, std::span<const NastranMaterialDeclaration> rows,
                     Work& work) {
    for (std::size_t i = 0; i < rows.size(); ++i) {
        work.charge(1, rows[i].id_field.position);
        keys[i] = {rows[i].id.value(), i};
    }
    std::ranges::sort(
        keys,
        [&work](std::uint64_t a, std::uint64_t b) {
            work.charge(1);
            return a < b;
        },
        &MaterialKey::id);
}

void resolve_material(NastranMaterialReference& reference, std::span<const MaterialKey> keys,
                      std::span<const NastranMaterialDeclaration> rows, Work& work) {
    const auto pos = reference.field.position;
    const auto found = std::ranges::lower_bound(
        keys, reference.id.value(),
        [&work, pos](std::uint64_t a, std::uint64_t b) {
            work.charge(1, pos);
            return a < b;
        },
        &MaterialKey::id);
    if (found == keys.end() || found->id != reference.id.value())
        reject(Code::InvalidInput, Reason::UndefinedMaterial, pos);
    if (rows[found->ordinal].kind != NastranMaterialKind::Mat1)
        reject(Code::InvalidInput, Reason::WrongMaterialKind, pos);
    reference.declaration_ordinal = found->ordinal;
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

std::size_t preflight(const NastranBeamView& beams, const NastranBulkView& bulk,
                      NastranPbeamLimits limits, Work& work) {
    std::size_t stations = 0;
    for (std::size_t i = 0; i < beams.properties.size(); ++i) {
        const auto& declaration = beams.properties[i];
        const auto& card = bulk.cards[declaration.card_ordinal];
        work.charge(1, card.name.position);
        if (i == limits.properties)
            reject(Code::ResourceExceeded, Reason::Limit, declaration.id_field.position);
        const auto emit = [&](const NastranPbeamStation& station) {
            if (stations == limits.stations)
                reject(Code::ResourceExceeded, Reason::Limit, station.output_field.position);
            ++stations;
        };
        static_cast<void>(read_property(Fields{bulk, card, work}, declaration, work, emit));
    }
    return stations;
}

} // namespace

std::span<const NastranPbeamStation>
NastranPbeamView::sections(const NastranPbeam& key) const noexcept {
    if (key.first_station > stations.size() ||
        key.station_count > stations.size() - key.first_station)
        return {};
    return stations.subspan(key.first_station, key.station_count);
}

NastranPbeamCards::NastranPbeamCards(NastranBeamOrigin source,
                                     core::PlannedBuffer<NastranPbeam>&& properties,
                                     core::PlannedBuffer<NastranPbeamStation>&& stations,
                                     std::uint64_t work) noexcept
    : source_(source),
      properties_(std::move(properties)),
      stations_(std::move(stations)),
      work_(work) {}

NastranPbeamCards::NastranPbeamCards(NastranPbeamCards&& other) noexcept
    : source_(std::exchange(other.source_, std::nullopt)),
      properties_(std::move(other.properties_)),
      stations_(std::move(other.stations_)),
      work_(other.work_) {}

std::optional<NastranPbeamView> NastranPbeamCards::view() const& noexcept {
    if (!source_)
        return std::nullopt;
    return NastranPbeamView{*source_, properties_.values(), stations_.values(), work_};
}

core::Result<NastranPbeamCards>
NastranPbeamCards::prepare(const NastranIsotropicShellCards& source,
                           const NastranBeamCards& beam_cards, NastranPbeamLimits limits,
                           core::ResourcePlan& plan, const std::stop_token& stop,
                           std::optional<NastranPbeamDiagnostic>* diagnostic) noexcept {
    if (diagnostic)
        diagnostic->reset();
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<NastranPbeamCards> {
            try {
                Work work{limits.work, stop};
                work.charge(0);
                const auto material_view = source.view();
                const auto beams = beam_cards.view();
                if (!material_view || !beams)
                    reject(Code::InvalidInput, Reason::Input);
                const auto& bulk = material_view->input.input.bulk;
                if (!same_snapshot(bulk, beams->source))
                    reject(Code::InvalidInput, Reason::SourceMismatch);
                const auto count = preflight(*beams, bulk, limits, work);
                const auto material_rows = beams->properties.empty()
                                               ? material_view->materials.first(0)
                                               : material_view->materials;
                auto remaining = plan.remaining_bytes();
                admit<NastranPbeam>(beams->properties.size(), remaining, work);
                admit<NastranPbeamStation>(count, remaining, work);
                admit<MaterialKey>(material_rows.size(), remaining, work);

                auto properties = allocate<NastranPbeam>(plan, beams->properties.size());
                auto stations = allocate<NastranPbeamStation>(plan, count);
                auto keys = allocate<MaterialKey>(plan, material_rows.size());
                index_materials(keys.values(), material_rows, work);
                std::size_t station_slot = 0;
                for (std::size_t i = 0; i < beams->properties.size(); ++i) {
                    const auto& declaration = beams->properties[i];
                    const auto& card = bulk.cards[declaration.card_ordinal];
                    work.charge(1, card.name.position);
                    const auto first = station_slot;
                    const auto emit = [&](const NastranPbeamStation& station) {
                        stations.values()[station_slot++] = station;
                    };
                    auto row = read_property(Fields{bulk, card, work}, declaration, work, emit);
                    resolve_material(row.material, keys.values(), material_rows, work);
                    row.first_station = first;
                    properties.values()[i] = row;
                }
                work.charge(0);
                return NastranPbeamCards{beams->source, std::move(properties), std::move(stations),
                                         work.used()};
            } catch (const Abort& error) {
                if (diagnostic && error.failure.code != Code::Cancelled)
                    *diagnostic = error.diagnostic;
                return std::unexpected(error.failure);
            }
        });
}

} // namespace uvlas::io::detail
