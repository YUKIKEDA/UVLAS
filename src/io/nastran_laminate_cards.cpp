#include "nastran_laminate_cards.hpp"

#include <algorithm>
#include <charconv>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
using Reason = NastranLaminateFailure;
using Position = NastranBulkPosition;
using Lexeme = NastranBulkLexeme;
constexpr auto operation = core::Operation::PrepareModel;
constexpr std::string_view nonzero_digits = "123456789";

struct Abort {
    core::Failure failure;
    NastranLaminateDiagnostic diagnostic;
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

bool keyword(std::string_view actual, std::string_view expected) noexcept {
    return std::ranges::equal(actual, expected, [](char a, char b) {
        return (a >= 'a' && a <= 'z' ? a - ('a' - 'A') : a) == b;
    });
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

    [[nodiscard]] std::string_view text(Lexeme field) const noexcept {
        return bulk_.text(field);
    }

    [[nodiscard]] std::uint64_t positive_id(Lexeme field) const {
        work_.charge(field.size + 1, field.position);
        const auto result = bulk_.integer(field);
        if (!result)
            reject(result.error().code, Reason::InvalidField, field.position);
        if (*result <= 0 || (card_.format == NastranBulkFormat::Free && field.size > 8))
            reject(Code::InvalidInput, Reason::InvalidField, field.position);
        return static_cast<std::uint64_t>(*result);
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

    [[nodiscard]] std::size_t size() const noexcept {
        return slots_.size();
    }

private:
    const NastranBulkView& bulk_;
    const NastranBulkCard& card_;
    std::span<const Lexeme> slots_;
    Work& work_;
};

template <class Enum, std::size_t N>
Enum choice(const Fields& fields, Lexeme field,
            const std::array<std::pair<std::string_view, Enum>, N>& known, Work& work) {
    const auto input = fields.text(field);
    work.charge(field.size + 1, field.position);
    for (const auto& [word, value] : known)
        if (keyword(input, word))
            return value;
    reject(Code::Unsupported, Reason::UnsupportedField, field.position);
}

void require_sign(const std::optional<NastranShellReal>& value, bool allow_zero = false) {
    if (!value)
        return;
    const auto decimal = value->value.decimal();
    const auto significand = decimal.substr(0, decimal.find('e'));
    const bool zero = significand.find_first_of(nonzero_digits) == std::string_view::npos;
    if (zero ? !allow_zero : significand.starts_with('-'))
        reject(Code::InvalidInput, Reason::InvalidField, value->field.position);
}

bool is_one(std::string_view decimal) noexcept {
    const auto marker = decimal.find('e');
    const auto significand = decimal.substr(0, marker);
    const auto first = significand.find_first_of(nonzero_digits);
    if (significand.starts_with('-') || first == std::string_view::npos ||
        significand[first] != '1' ||
        significand.find_first_of(nonzero_digits, first + 1) != std::string_view::npos)
        return false;
    const auto point = static_cast<int>(significand.find('.'));
    const auto distance = point - static_cast<int>(first);
    const auto mantissa_power = distance > 0 ? distance - 1 : distance;
    if (marker == std::string_view::npos)
        return mantissa_power == 0;
    auto text = decimal.substr(marker + 1);
    if (text.starts_with('+'))
        text.remove_prefix(1);
    int exponent = 0;
    const auto end = text.data() + text.size();
    const auto parsed = std::from_chars(text.data(), end, exponent);
    return parsed.ec == std::errc{} && parsed.ptr == end && exponent == -mantissa_power;
}

NastranMat8 read_mat8(const Fields& fields, std::size_t material, Work& work) {
    NastranMat8 result;
    result.declaration_ordinal = material;
    for (const auto i : {1U, 2U, 3U})
        if (!fields.supplied(i))
            reject(Code::InvalidInput, Reason::MissingField, fields.position(i));
    result.longitudinal_young = fields.real(1);
    result.transverse_young = fields.real(2);
    result.poisson_12 = fields.real(3);
    result.shear_12 = fields.real(4);
    result.shear_1z = fields.real(5);
    result.shear_2z = fields.real(6);
    result.density = fields.real(7);
    result.expansion_1 = fields.real(8);
    result.expansion_2 = fields.real(9);
    result.reference_temperature = fields.real(10);
    result.longitudinal_tension = fields.real(11);
    result.longitudinal_compression = fields.real(12);
    result.transverse_tension = fields.real(13);
    result.transverse_compression = fields.real(14);
    result.shear_limit = fields.real(15);
    result.damping = fields.real(16);
    result.tsai_wu_interaction = fields.real(17);
    result.strain_allowables = fields.real(18);
    require_sign(result.longitudinal_young);
    require_sign(result.transverse_young);
    require_sign(result.shear_12);
    require_sign(result.shear_1z, true);
    require_sign(result.shear_2z, true);
    require_sign(result.longitudinal_tension);
    require_sign(result.longitudinal_compression);
    require_sign(result.transverse_tension);
    require_sign(result.transverse_compression);
    require_sign(result.shear_limit);
    if (result.strain_allowables && !is_one(result.strain_allowables->value.decimal()))
        reject(Code::InvalidInput, Reason::InvalidField, result.strain_allowables->field.position);
    for (std::size_t i = 19; i < fields.size(); ++i) {
        work.charge(1, fields.position(i));
        if (const auto field = fields.supplied(i))
            reject(Code::Unsupported, Reason::UnsupportedField, field->position);
    }
    return result;
}

NastranPcomp read_pcomp_header(const Fields& fields, std::size_t property, Work& work) {
    using Lam = NastranLaminateMode;
    constexpr std::array<std::pair<std::string_view, Lam>, 5> modes{
        {{"SYM", Lam::SymmetricHalf},
         {"MEM", Lam::MembraneOnly},
         {"BEND", Lam::BendingOnly},
         {"SMEAR", Lam::Smeared},
         {"SMCORE", Lam::SmearedWithCore}}};
    using Ft = NastranPcompFailureTheory;
    constexpr std::array<std::pair<std::string_view, Ft>, 7> theories{
        {{"HILL", Ft::Hill},
         {"HOFF", Ft::Hoffman},
         {"TSAI", Ft::TsaiWu},
         {"STRN", Ft::MaximumStrain},
         {"HFAIL", Ft::Hashin},
         {"HTAPE", Ft::HashinTape},
         {"HFABR", Ft::HashinFabric}}};
    NastranPcomp result;
    result.property_ordinal = property;
    result.bottom_offset = fields.real(1);
    result.nonstructural_mass = fields.real(2);
    result.bond_shear_limit = fields.real(3);
    require_sign(result.bond_shear_limit);
    if (const auto theory = fields.supplied(4)) {
        result.failure_theory =
            NastranPcompTheoryField{choice(fields, *theory, theories, work), *theory};
        if (!result.bond_shear_limit)
            reject(Code::InvalidInput, Reason::MissingField, fields.position(3));
    }
    result.reference_temperature = fields.real(5);
    result.damping = fields.real(6);
    if (const auto mode = fields.supplied(7)) {
        result.mode = choice(fields, *mode, modes, work);
        result.mode_field = mode;
    }
    return result;
}

struct Key {
    std::uint64_t id = 0;
    std::size_t ordinal = 0;
};

class MaterialIndex {
public:
    MaterialIndex(std::span<Key> keys, std::span<const NastranMaterialDeclaration> materials,
                  Work& work)
        : keys_(keys), work_(work) {
        for (std::size_t i = 0; i < materials.size(); ++i) {
            work_.charge(1, materials[i].id_field.position);
            keys[i] = {materials[i].id.value(), i};
        }
        std::ranges::sort(
            keys,
            [this](std::uint64_t lhs, std::uint64_t rhs) {
                work_.charge(1);
                return lhs < rhs;
            },
            &Key::id);
    }

    [[nodiscard]] NastranMaterialReference reference(const Fields& fields, Lexeme field) const {
        const auto id = fields.positive_id(field);
        const auto found = std::ranges::lower_bound(
            keys_, id,
            [this, field](std::uint64_t lhs, std::uint64_t rhs) {
                work_.charge(1, field.position);
                return lhs < rhs;
            },
            &Key::id);
        if (found == keys_.end() || found->id != id)
            reject(Code::InvalidInput, Reason::UndefinedMaterial, field.position);
        return {NastranMaterialId{id}, field, found->ordinal};
    }

private:
    std::span<const Key> keys_;
    Work& work_;
};

constexpr std::size_t first_ply_field = 8;
constexpr std::size_t fields_per_ply = 4;

template <class Emit>
void visit_plies(const Fields& fields, Work& work, Emit emit) {
    for (std::size_t i = first_ply_field; i < fields.size(); i += fields_per_ply) {
        work.charge(1, fields.position(i));
        for (std::size_t j = 0; j < fields_per_ply; ++j) {
            if (const auto source = fields.supplied(i + j)) {
                emit(i, *source);
                break;
            }
        }
    }
}

void read_plies(NastranPcomp& pcomp, const Fields& fields, std::span<NastranPcompPly> plies,
                const MaterialIndex& materials, Work& work) {
    std::size_t material_source = 0, thickness_source = 0;
    constexpr std::array<std::pair<std::string_view, bool>, 2> stress_modes{
        {{"YES", true}, {"NO", false}}};
    visit_plies(fields, work, [&](std::size_t start, Lexeme source) {
        NastranPcompPly ply;
        ply.source = source;
        const auto material = fields.supplied(start);
        if (material) {
            ply.material = materials.reference(fields, *material);
            material_source = pcomp.ply_count;
        } else if (pcomp.ply_count == 0) {
            reject(Code::InvalidInput, Reason::MissingField, fields.position(start));
        }
        ply.material_source_ply = material_source;
        ply.thickness = fields.real(start + 1);
        if (ply.thickness) {
            thickness_source = pcomp.ply_count;
        } else if (pcomp.ply_count == 0) {
            reject(Code::InvalidInput, Reason::MissingField, fields.position(start + 1));
        }
        ply.thickness_source_ply = thickness_source;
        ply.angle = fields.real(start + 2);
        if (const auto stress = fields.supplied(start + 3))
            ply.stress_output =
                NastranPcompStressOutput{choice(fields, *stress, stress_modes, work), *stress};

        plies[pcomp.first_ply + pcomp.ply_count++] = ply;
    });
    if (pcomp.ply_count == 0)
        reject(Code::InvalidInput, Reason::MissingField, fields.position(first_ply_field));
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
    std::size_t mat8 = 0, pcomp = 0, plies = 0;
};

Counts count_records(const NastranIsotropicShellView& input, NastranLaminateLimits limits,
                     Work& work) {
    const auto& bulk = input.input.input.bulk;
    Counts counts;
    for (const auto& material : input.materials) {
        work.charge(1, material.id_field.position);
        if (material.kind == NastranMaterialKind::Mat8)
            ++counts.mat8;
    }
    for (const auto& property : input.input.properties) {
        work.charge(1, property.id_field.position);
        if (property.kind != NastranShellPropertyKind::Pcomp)
            continue;
        ++counts.pcomp;
        const Fields fields{bulk, bulk.cards[property.card_ordinal], work};
        visit_plies(fields, work, [&](std::size_t, Lexeme source) {
            if (counts.plies == limits.plies)
                reject(Code::ResourceExceeded, Reason::Limit, source.position);
            ++counts.plies;
        });
    }
    return counts;
}

} // namespace

std::span<const NastranPcompPly>
NastranLaminateView::layers(const NastranPcomp& property) const noexcept {
    if (property.first_ply > plies.size() || property.ply_count > plies.size() - property.first_ply)
        return {};
    return plies.subspan(property.first_ply, property.ply_count);
}

NastranLaminateCards::NastranLaminateCards(NastranIsotropicShellCards&& input,
                                           core::PlannedBuffer<NastranMat8>&& mat8,
                                           core::PlannedBuffer<NastranPcomp>&& pcomp,
                                           core::PlannedBuffer<NastranPcompPly>&& plies,
                                           std::uint64_t work) noexcept
    : input_(std::move(input)),
      mat8_(std::move(mat8)),
      pcomp_(std::move(pcomp)),
      plies_(std::move(plies)),
      work_(work) {}

std::optional<NastranLaminateView> NastranLaminateCards::view() const& noexcept {
    const auto input = input_.view();
    if (!input)
        return std::nullopt;
    return NastranLaminateView{*input, mat8_.values(), pcomp_.values(), plies_.values(), work_};
}

core::Result<NastranLaminateCards>
NastranLaminateCards::prepare(NastranIsotropicShellCards&& source, NastranLaminateLimits limits,
                              core::ResourcePlan& plan, const std::stop_token& stop,
                              std::optional<NastranLaminateDiagnostic>* diagnostic) noexcept {
    if (diagnostic)
        diagnostic->reset();
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<NastranLaminateCards> {
            try {
                Work work{limits.work, stop};
                work.charge(0);
                const auto input = source.view();
                if (!input)
                    reject(Code::InvalidInput, Reason::Input);
                const auto& bulk = input->input.input.bulk;
                const auto counts = count_records(*input, limits, work);
                const auto index_size = counts.pcomp == 0 ? 0 : input->materials.size();
                auto remaining = plan.remaining_bytes();
                admit<NastranMat8>(counts.mat8, remaining, work);
                admit<NastranPcomp>(counts.pcomp, remaining, work);
                admit<NastranPcompPly>(counts.plies, remaining, work);
                admit<Key>(index_size, remaining, work);

                auto mat8 = allocate<NastranMat8>(plan, counts.mat8);
                auto pcomp = allocate<NastranPcomp>(plan, counts.pcomp);
                auto plies = allocate<NastranPcompPly>(plan, counts.plies);
                auto keys = allocate<Key>(plan, index_size);
                const MaterialIndex materials{keys.values(), input->materials.first(index_size),
                                              work};

                std::size_t mat8_slot = 0;
                for (std::size_t i = 0; i < input->materials.size(); ++i) {
                    const auto& material = input->materials[i];
                    work.charge(1, material.id_field.position);
                    if (material.kind == NastranMaterialKind::Mat8) {
                        const Fields fields{bulk, bulk.cards[material.card_ordinal], work};
                        mat8.values()[mat8_slot++] = read_mat8(fields, i, work);
                    }
                }
                std::size_t pcomp_slot = 0, ply_slot = 0;
                for (std::size_t i = 0; i < input->input.properties.size(); ++i) {
                    const auto& property = input->input.properties[i];
                    work.charge(1, property.id_field.position);
                    if (property.kind == NastranShellPropertyKind::Pcomp) {
                        const Fields fields{bulk, bulk.cards[property.card_ordinal], work};
                        auto description = read_pcomp_header(fields, i, work);
                        description.first_ply = ply_slot;
                        read_plies(description, fields, plies.values(), materials, work);
                        ply_slot += description.ply_count;
                        pcomp.values()[pcomp_slot++] = description;
                    }
                }
                work.charge(0);
                return NastranLaminateCards{std::move(source), std::move(mat8), std::move(pcomp),
                                            std::move(plies), work.used()};
            } catch (const Abort& error) {
                if (diagnostic && error.failure.code != Code::Cancelled)
                    *diagnostic = error.diagnostic;
                return std::unexpected(error.failure);
            }
        });
}

} // namespace uvlas::io::detail
