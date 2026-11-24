#include "nastran_isotropic_shell_cards.hpp"

#include <algorithm>
#include <limits>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
using Reason = NastranIsotropicShellFailure;
using Position = NastranBulkPosition;
using Lexeme = NastranBulkLexeme;
constexpr auto operation = core::Operation::PrepareModel;
constexpr std::uint64_t mat8_id_limit = 100000000;

struct Abort {
    core::Failure failure;
    NastranIsotropicShellDiagnostic diagnostic;
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

std::optional<NastranMaterialKind> material_kind(std::string_view name) noexcept {
    using enum NastranMaterialKind;
    if (name.size() != 4 || (name[0] != 'M' && name[0] != 'm') ||
        (name[1] != 'A' && name[1] != 'a') || (name[2] != 'T' && name[2] != 't'))
        return std::nullopt;
    if (name[3] == '1')
        return Mat1;
    if (name[3] == '8')
        return Mat8;
    return std::nullopt;
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

    [[nodiscard]] std::int64_t integer(Lexeme field) const {
        work_.charge(field.size + 1, field.position);
        if (card_.format == NastranBulkFormat::Free && field.size > 8)
            reject(Code::InvalidInput, Reason::InvalidField, field.position);
        const auto result = bulk_.integer(field);
        if (!result)
            reject(result.error().code, Reason::InvalidField, field.position);
        return *result;
    }

    [[nodiscard]] std::uint64_t positive_id(Lexeme field) const {
        const auto id = integer(field);
        if (id <= 0)
            reject(Code::InvalidInput, Reason::InvalidField, field.position);
        return static_cast<std::uint64_t>(id);
    }

    [[nodiscard]] std::optional<NastranShellReal> real(std::size_t slot) const {
        const auto field = supplied(slot);
        if (!field)
            return std::nullopt;
        work_.charge(field->size + 1, field->position);
        const auto result = bulk_.real(*field);
        if (!result)
            reject(result.error().code, Reason::InvalidField, field->position);
        return NastranShellReal{*field, *result};
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

std::string_view mantissa(std::string_view decimal) noexcept {
    return decimal.substr(0, decimal.find('e'));
}

constexpr std::string_view nonzero_digits = "123456789";

void require_positive(const std::optional<NastranShellReal>& value) {
    if (!value)
        return;
    const auto digits = mantissa(value->value.decimal());
    if (digits.starts_with('-') || digits.find_first_of(nonzero_digits) == std::string_view::npos)
        reject(Code::InvalidInput, Reason::InvalidField, value->field.position);
}

int exponent(std::string_view digits) noexcept {
    const bool negative = digits.starts_with('-');
    if (negative || digits.starts_with('+'))
        digits.remove_prefix(1);
    // The mantissa occupies at most one physical line. Once the exponent is outside
    // this range, the sign alone determines the range of any nonzero mantissa here.
    constexpr int cap = 256;
    int value = 0;
    for (const char digit : digits)
        value = std::min(cap, 10 * value + digit - '0');
    return negative ? -value : value;
}

void check_poisson(const std::optional<NastranShellReal>& value) {
    if (!value)
        return;
    const auto decimal = value->value.decimal();
    const auto significand = mantissa(decimal);
    const auto first = significand.find_first_of(nonzero_digits);
    if (first == std::string_view::npos)
        return;
    const auto point = significand.find('.');
    const auto distance = static_cast<int>(point) - static_cast<int>(first);
    const auto leading_power =
        (distance > 0 ? distance - 1 : distance) +
        (significand.size() == decimal.size() ? 0
                                              : exponent(decimal.substr(significand.size() + 1)));
    // MAT1 requires -1 < NU < 0.5, including exact bounds written with exponents.
    const bool valid = significand.starts_with('-') ? leading_power < 0
                       : leading_power == -1        ? significand[first] < '5'
                                                    : leading_power < -1;
    if (!valid)
        reject(Code::InvalidInput, Reason::InvalidField, value->field.position);
}

struct Key {
    std::uint64_t id = 0;
    std::size_t ordinal = 0;
    Position position;
};

template <class Row>
void index_rows(std::span<Key> keys, std::span<const Row> rows, Work& work) {
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto position = rows[i].id_field.position;
        work.charge(1, position);
        keys[i] = {rows[i].id.value(), i, position};
    }
    std::ranges::sort(keys, [&work](const Key& a, const Key& b) {
        work.charge(1, a.position);
        return a.id < b.id || (a.id == b.id && a.ordinal < b.ordinal);
    });
}

void require_unique_materials(std::span<const Key> keys, Work& work) {
    std::optional<Position> duplicate;
    for (std::size_t i = 1; i < keys.size(); ++i) {
        const auto position = keys[i].position;
        work.charge(1, position);
        if (keys[i].id == keys[i - 1].id && (!duplicate || position.byte < duplicate->byte))
            duplicate = position;
    }
    if (duplicate)
        reject(Code::InvalidInput, Reason::DuplicateMaterial, *duplicate);
}

std::size_t resolve(std::span<const Key> keys, std::uint64_t id, Position position, Reason reason,
                    Work& work) {
    const auto found = std::ranges::lower_bound(
        keys, id,
        [&work, position](std::uint64_t a, std::uint64_t b) {
            work.charge(1, position);
            return a < b;
        },
        &Key::id);
    if (found == keys.end() || found->id != id)
        reject(Code::InvalidInput, reason, position);
    return found->ordinal;
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
    std::size_t materials = 0, mat1 = 0, pshell = 0;
};

Counts counts(const NastranShellView& input, NastranIsotropicShellLimits limits, Work& work) {
    Counts result;
    for (const auto& card : input.input.bulk.cards) {
        work.charge(1, card.name.position);
        const auto kind = material_kind(input.input.bulk.text(card.name));
        if (!kind)
            continue;
        if (result.materials == limits.material_ids)
            reject(Code::ResourceExceeded, Reason::Limit, card.name.position);
        ++result.materials;
        if (*kind == NastranMaterialKind::Mat1)
            ++result.mat1;
    }
    for (const auto& property : input.properties) {
        work.charge(1, property.id_field.position);
        if (property.kind == NastranShellPropertyKind::Pshell)
            ++result.pshell;
    }
    return result;
}

void read_material_declarations(const NastranBulkView& bulk,
                                std::span<NastranMaterialDeclaration> declarations, Work& work) {
    std::size_t material = 0, mat1 = 0;
    for (std::size_t i = 0; i < bulk.cards.size(); ++i) {
        const auto& card = bulk.cards[i];
        work.charge(1, card.name.position);
        const auto kind = material_kind(bulk.text(card.name));
        if (!kind)
            continue;
        const Fields fields{bulk, card, work};
        const auto field = fields.supplied(0);
        if (!field)
            reject(Code::InvalidInput, Reason::MissingField, fields.position(0));
        const auto id = fields.positive_id(*field);
        if (*kind == NastranMaterialKind::Mat8 && id >= mat8_id_limit)
            reject(Code::InvalidInput, Reason::InvalidField, field->position);
        declarations[material++] = {NastranMaterialId{id}, *kind, *field, i,
                                    *kind == NastranMaterialKind::Mat1 ? std::optional{mat1++}
                                                                       : std::nullopt};
    }
}

NastranMat1 read_mat1(const Fields& fields, std::size_t declaration,
                      std::span<const Key> coordinates, Work& work) {
    NastranMat1 result;
    result.declaration_ordinal = declaration;
    result.young = fields.real(1);
    result.shear = fields.real(2);
    result.poisson = fields.real(3);
    if (!result.young && !result.shear)
        reject(Code::InvalidInput, Reason::MissingField, fields.position(1));
    require_positive(result.young);
    require_positive(result.shear);
    check_poisson(result.poisson);
    result.density = fields.real(4);
    result.expansion = fields.real(5);
    result.reference_temperature = fields.real(6);
    result.damping = fields.real(7);
    result.tension_limit = fields.real(8);
    result.compression_limit = fields.real(9);
    result.shear_limit = fields.real(10);
    if (const auto field = fields.supplied(11)) {
        const auto id = fields.positive_id(*field);
        result.material_coordinate = NastranShellMaterialAxis{
            NastranCoordinateId{id}, *field,
            resolve(coordinates, id, field->position, Reason::UndefinedCoordinate, work)};
    }
    fields.check_tail(12);
    return result;
}

class Materials {
public:
    Materials(std::span<const NastranMaterialDeclaration> rows, std::span<const Key> keys,
              Work& work) noexcept
        : rows_(rows), keys_(keys), work_(work) {}

    [[nodiscard]] NastranMaterialReference reference(std::uint64_t id, Lexeme field) const {
        return {NastranMaterialId{id}, field,
                resolve(keys_, id, field.position, Reason::UndefinedMaterial, work_)};
    }

    [[nodiscard]] std::optional<NastranMaterialReference> reference(const Fields& fields,
                                                                    std::size_t slot) const {
        const auto field = fields.supplied(slot);
        return field ? std::optional{reference(fields.positive_id(*field), *field)} : std::nullopt;
    }

    [[nodiscard]] NastranMaterialKind kind(const NastranMaterialReference& reference) const {
        return rows_[reference.declaration_ordinal].kind;
    }

private:
    std::span<const NastranMaterialDeclaration> rows_;
    std::span<const Key> keys_;
    Work& work_;
};

NastranBendingMaterial read_bending(const Fields& fields, const Materials& materials) {
    constexpr std::size_t field_index = 3;
    const auto field = fields.supplied(field_index);
    if (!field)
        return std::monostate{};
    const auto id = fields.integer(*field);
    if (id < -1)
        reject(Code::InvalidInput, Reason::InvalidField, field->position);
    if (id <= 0) {
        return NastranBendingDesignation{
            id == -1 ? NastranSpecialBending::PlaneStrain : NastranSpecialBending::Zero, *field};
    }
    return materials.reference(static_cast<std::uint64_t>(id), *field);
}

void classify_pshell(NastranPshell& result, const Materials& materials) {
    std::optional<NastranMaterialId> first;
    for (const auto* reference : {result.membrane ? &*result.membrane : nullptr,
                                  std::get_if<NastranMaterialReference>(&result.bending),
                                  result.transverse_shear ? &*result.transverse_shear : nullptr,
                                  result.coupling ? &*result.coupling : nullptr}) {
        if (!reference)
            continue;
        if (!first)
            first = reference->id;
        else if (*first != reference->id)
            result.distinct_materials = true;
        if (materials.kind(*reference) != NastranMaterialKind::Mat1)
            result.orthotropic_material = true;
    }
}

NastranPshell read_pshell(const Fields& fields, std::size_t property, const Materials& materials) {
    NastranPshell result;
    result.property_ordinal = property;
    result.membrane = materials.reference(fields, 1);
    result.thickness = fields.real(2);
    result.bending = read_bending(fields, materials);
    result.bending_ratio = fields.real(4);
    require_positive(result.bending_ratio);
    result.transverse_shear = materials.reference(fields, 5);
    result.shear_ratio = fields.real(6);
    require_positive(result.shear_ratio);
    result.nonstructural_mass = fields.real(7);
    result.first_fiber = fields.real(8);
    result.second_fiber = fields.real(9);
    result.coupling = materials.reference(fields, 10);
    fields.check_tail(11);

    const auto* special = std::get_if<NastranBendingDesignation>(&result.bending);
    const bool plane_strain = special && special->kind == NastranSpecialBending::PlaneStrain;
    if (plane_strain &&
        (!result.membrane || materials.kind(*result.membrane) != NastranMaterialKind::Mat1))
        reject(Code::InvalidInput, Reason::IncompatibleMaterial, fields.position(1));
    if (result.transverse_shear &&
        (plane_strain || std::holds_alternative<std::monostate>(result.bending)))
        reject(Code::InvalidInput, Reason::IncompatibleMaterial,
               result.transverse_shear->field.position);
    if (result.coupling) {
        const auto* bending = std::get_if<NastranMaterialReference>(&result.bending);
        if (!result.membrane || !bending || result.coupling->id == result.membrane->id ||
            result.coupling->id == bending->id)
            reject(Code::InvalidInput, Reason::IncompatibleMaterial,
                   result.coupling->field.position);
    }
    classify_pshell(result, materials);
    return result;
}

} // namespace

NastranIsotropicShellCards::NastranIsotropicShellCards(
    NastranShellCards&& input, core::PlannedBuffer<NastranMaterialDeclaration>&& materials,
    core::PlannedBuffer<NastranMat1>&& mat1, core::PlannedBuffer<NastranPshell>&& pshell,
    std::uint64_t work) noexcept
    : input_(std::move(input)),
      materials_(std::move(materials)),
      mat1_(std::move(mat1)),
      pshell_(std::move(pshell)),
      work_(work) {}

std::optional<NastranIsotropicShellView> NastranIsotropicShellCards::view() const& noexcept {
    const auto input = input_.view();
    if (!input)
        return std::nullopt;
    return NastranIsotropicShellView{*input, materials_.values(), mat1_.values(), pshell_.values(),
                                     work_};
}

core::Result<NastranIsotropicShellCards> NastranIsotropicShellCards::prepare(
    NastranShellCards&& source, NastranIsotropicShellLimits limits, core::ResourcePlan& plan,
    const std::stop_token& stop,
    std::optional<NastranIsotropicShellDiagnostic>* diagnostic) noexcept {
    if (diagnostic)
        diagnostic->reset();
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<NastranIsotropicShellCards> {
            try {
                Work work{limits.work, stop};
                work.charge(0);
                const auto input = source.view();
                if (!input)
                    reject(Code::InvalidInput, Reason::Input);
                const auto& bulk = input->input.bulk;
                const auto count = counts(*input, limits, work);
                auto remaining = plan.remaining_bytes();
                admit<NastranMaterialDeclaration>(count.materials, remaining, work);
                admit<NastranMat1>(count.mat1, remaining, work);
                admit<NastranPshell>(count.pshell, remaining, work);
                admit<Key>(count.materials, remaining, work);
                admit<Key>(input->input.coordinates.size(), remaining, work);

                auto declarations = allocate<NastranMaterialDeclaration>(plan, count.materials);
                auto mat1 = allocate<NastranMat1>(plan, count.mat1);
                auto pshell = allocate<NastranPshell>(plan, count.pshell);
                auto material_index = allocate<Key>(plan, count.materials);
                auto coordinate_index = allocate<Key>(plan, input->input.coordinates.size());
                read_material_declarations(bulk, declarations.values(), work);
                const auto rows = std::as_const(declarations).values();
                index_rows(material_index.values(), rows, work);
                require_unique_materials(material_index.values(), work);
                index_rows(coordinate_index.values(), input->input.coordinates, work);
                const Materials materials{rows, material_index.values(), work};
                for (std::size_t i = 0; i < rows.size(); ++i) {
                    const auto& row = rows[i];
                    work.charge(1, row.id_field.position);
                    if (const auto slot = row.mat1_ordinal) {
                        const Fields fields{bulk, bulk.cards[row.card_ordinal], work};
                        mat1.values()[*slot] =
                            read_mat1(fields, i, coordinate_index.values(), work);
                    }
                }
                std::size_t pshell_slot = 0;
                for (std::size_t i = 0; i < input->properties.size(); ++i) {
                    const auto& property = input->properties[i];
                    work.charge(1, property.id_field.position);
                    if (property.kind == NastranShellPropertyKind::Pshell) {
                        const Fields fields{bulk, bulk.cards[property.card_ordinal], work};
                        pshell.values()[pshell_slot++] = read_pshell(fields, i, materials);
                    }
                }
                work.charge(0);
                return NastranIsotropicShellCards{std::move(source), std::move(declarations),
                                                  std::move(mat1), std::move(pshell), work.used()};
            } catch (const Abort& error) {
                if (diagnostic && error.failure.code != Code::Cancelled)
                    *diagnostic = error.diagnostic;
                return std::unexpected(error.failure);
            }
        });
}

} // namespace uvlas::io::detail
