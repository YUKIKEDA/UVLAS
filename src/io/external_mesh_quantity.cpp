#include "external_mesh_quantity.hpp"

#include <algorithm>
#include <array>
#include <tuple>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Rational = structure::ExactRational;
using Code = core::FailureCode;
using Field = ImportedQuantityField;
constexpr auto operation = core::Operation::PrepareModel;
constexpr SchemaDimension thickness_dimension{{1, 0, 0, 0}};

[[nodiscard]] std::unexpected<core::Failure> fail(Code code) noexcept {
    return std::unexpected(core::Failure{code, operation});
}

[[nodiscard]] core::Result<Rational> imported_si(std::string_view value, const SiUnits& units,
                                                 SchemaDimension dimensions,
                                                 const structure::ExactLimits& limits,
                                                 core::ResourcePlan& plan) noexcept {
    auto input = Rational::from_decimal(value, limits, plan);
    if (!input) {
        return std::unexpected(input.error());
    }
    const SiRecipeView multiplier{SiDecimalView{1, "1", 0}, {}, units, dimensions, 0};
    const auto unit = materialize_rational(multiplier, limits, plan);
    if (!unit) {
        return std::unexpected(unit.error());
    }
    return input->multiply(*unit, limits, plan);
}

[[nodiscard]] core::Result<Rational> case_si(SchemaNumberView value, const SiUnits& units,
                                             const structure::ExactLimits& limits,
                                             core::ResourcePlan& plan) noexcept {
    const SiRecipeView recipe{
        SiDecimalView{value.sign, value.digits, value.exponent10}, {}, units, value.dimensions, 0};
    return materialize_rational(recipe, limits, plan);
}

struct MaterialCheck {
    const SiUnits& external_units;
    const SiUnits& case_units;
    const structure::ExactLimits& limits;
    core::ResourcePlan& plan;
    const std::stop_token& stop;
    std::optional<ImportedQuantityMismatch>* diagnostic;

    [[nodiscard]] core::Result<void> match(const NastranShellReal& external,
                                           const Rational& expected, SchemaDimension dimension,
                                           Field field) const noexcept {
        if (stop.stop_requested()) {
            return fail(Code::Cancelled);
        }
        const auto actual =
            imported_si(external.value.decimal(), external_units, dimension, limits, plan);
        if (!actual) {
            return std::unexpected(actual.error());
        }
        if (stop.stop_requested()) {
            return fail(Code::Cancelled);
        }
        if (*actual != expected) {
            if (diagnostic) {
                *diagnostic = ImportedQuantityMismatch{field, external.field.position, {}};
            }
            return fail(Code::InvalidInput);
        }
        return {};
    }

    [[nodiscard]] core::Result<void> decimal(const std::optional<NastranShellReal>& external,
                                             SchemaNumberView expected,
                                             Field field) const noexcept {
        if (stop.stop_requested()) {
            return fail(Code::Cancelled);
        }
        if (!external) {
            return {};
        }
        const auto exact = case_si(expected, case_units, limits, plan);
        if (!exact) {
            return std::unexpected(exact.error());
        }
        return match(*external, *exact, expected.dimensions, field);
    }

    [[nodiscard]] core::Result<void>
    isotropic_shear(const std::optional<NastranShellReal>& external,
                    const schema::MaterialIsotropic& expected) const noexcept {
        if (!external) {
            return {};
        }
        const auto young = case_si(expected.young.number(), case_units, limits, plan);
        const auto nu = case_si(expected.poisson.number(), case_units, limits, plan);
        const auto two = Rational::from_decimal("2", limits, plan);
        const auto one = Rational::from_decimal("1", limits, plan);
        for (const auto* result : {&young, &nu, &two, &one}) {
            if (!*result) {
                return std::unexpected(result->error());
            }
        }
        const auto plus = one->add(*nu, limits, plan);
        if (!plus) {
            return std::unexpected(plus.error());
        }
        const auto denominator = two->multiply(*plus, limits, plan);
        if (!denominator) {
            return std::unexpected(denominator.error());
        }
        const auto shear = young->divide(*denominator, limits, plan);
        if (!shear) {
            return std::unexpected(shear.error());
        }
        return match(*external, *shear, expected.young.number().dimensions, Field::Shear12);
    }
};

[[nodiscard]] SiUnits length_units(SiLengthUnit length) noexcept {
    return {length, SiMassUnit::Kilogram, SiTimeUnit::Second, SiAngleUnit::Radian};
}

} // namespace

core::Result<void> check_imported_uniform_thickness(
    ImportedUniformThickness external, const schema::Thickness& field, SiLengthUnit case_length,
    std::uint64_t maximum_vertices, const structure::ExactLimits& limits, core::ResourcePlan& plan,
    const std::stop_token& stop, std::optional<ImportedQuantityMismatch>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    if (field.values.size() > maximum_vertices) {
        return fail(Code::ResourceExceeded);
    }
    const auto actual = imported_si(external.decimal, length_units(external.length),
                                    thickness_dimension, limits, plan);
    if (!actual) {
        return std::unexpected(actual.error());
    }
    if (actual->sign() <= 0) {
        return fail(Code::InvalidInput);
    }
    for (core::Index i = 0; std::cmp_less(i, field.values.size()); ++i) {
        if (stop.stop_requested()) {
            return fail(Code::Cancelled);
        }
        const auto vertex = field.values.at(i);
        if (!vertex) {
            return std::unexpected(vertex.error());
        }
        const auto expected = case_si(vertex->number(), length_units(case_length), limits, plan);
        if (!expected) {
            return std::unexpected(expected.error());
        }
        if (*expected != *actual) {
            if (diagnostic) {
                *diagnostic = ImportedQuantityMismatch{Field::FullThickness, external.position, i};
            }
            return fail(Code::InvalidInput);
        }
    }
    return stop.stop_requested() ? core::Result<void>{fail(Code::Cancelled)} : core::Result<void>{};
}

core::Result<void> check_imported_pshell_thickness(
    const NastranShellElement& element, const NastranPshell& property, SiLengthUnit deck_length,
    const schema::Thickness& field, SiLengthUnit case_length, std::uint64_t maximum_vertices,
    const structure::ExactLimits& limits, core::ResourcePlan& plan, const std::stop_token& stop,
    std::optional<ImportedQuantityMismatch>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    if (element.property.property_ordinal != property.property_ordinal) {
        return fail(Code::InvalidInput);
    }
    const auto has_corner = [](const auto* corners) noexcept {
        return corners && std::ranges::any_of(*corners, [](const NastranShellCorner& corner) {
                   return corner.thickness.has_value();
               });
    };
    if (has_corner(std::get_if<0>(&element.corners)) ||
        has_corner(std::get_if<1>(&element.corners))) {
        return fail(Code::Unsupported);
    }
    if (!property.thickness) {
        return {};
    }
    return check_imported_uniform_thickness(
        {property.thickness->value.decimal(), deck_length, property.thickness->field.position},
        field, case_length, maximum_vertices, limits, plan, stop, diagnostic);
}

core::Result<void> check_imported_material(
    const NastranMat1& source, const ImportedCaseMaterial& target, const SiUnits& external_units,
    const SiUnits& case_units, const structure::ExactLimits& limits, core::ResourcePlan& plan,
    const std::stop_token& stop, std::optional<ImportedQuantityMismatch>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    const auto* expected = std::get_if<schema::MaterialIsotropic>(&target);
    if (!expected) {
        return fail(Code::Unsupported);
    }
    const MaterialCheck check{external_units, case_units, limits, plan, stop, diagnostic};
    if (const auto young = check.decimal(source.young, expected->young.number(), Field::Young1);
        !young) {
        return young;
    }
    if (const auto shear = check.isotropic_shear(source.shear, *expected); !shear) {
        return shear;
    }
    if (const auto nu = check.decimal(source.poisson, expected->poisson.number(), Field::Poisson12);
        !nu) {
        return nu;
    }
    return check.decimal(source.density, expected->density.number(), Field::Density);
}

core::Result<void> check_imported_material(
    const NastranMat8& source, const ImportedCaseMaterial& target, const SiUnits& external_units,
    const SiUnits& case_units, const structure::ExactLimits& limits, core::ResourcePlan& plan,
    const std::stop_token& stop, std::optional<ImportedQuantityMismatch>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    if (stop.stop_requested()) {
        return fail(Code::Cancelled);
    }
    const auto* expected = std::get_if<schema::MaterialOrthotropic>(&target);
    if (!expected) {
        return fail(Code::Unsupported);
    }
    const MaterialCheck check{external_units, case_units, limits, plan, stop, diagnostic};
    const auto array_field = [&](const std::optional<NastranShellReal>& card, auto values,
                                 core::Index ordinal, Field field) -> core::Result<void> {
        const auto value = values.at(ordinal);
        if (!value) {
            return std::unexpected(value.error());
        }
        return check.decimal(card, value->number(), field);
    };
    for (const auto& [card, ordinal, field] :
         std::array{std::tuple{&source.longitudinal_young, core::Index{0}, Field::Young1},
                    std::tuple{&source.transverse_young, core::Index{1}, Field::Young2}}) {
        if (const auto result = array_field(*card, expected->young, ordinal, field); !result) {
            return result;
        }
    }
    if (const auto result = array_field(source.poisson_12, expected->poisson, 0, Field::Poisson12);
        !result) {
        return result;
    }
    // FD-04 stores G23, G13, G12, while MAT8 lists G12, G1Z, G2Z.
    for (const auto& [card, ordinal, field] :
         std::array{std::tuple{&source.shear_12, core::Index{2}, Field::Shear12},
                    std::tuple{&source.shear_1z, core::Index{1}, Field::Shear13},
                    std::tuple{&source.shear_2z, core::Index{0}, Field::Shear23}}) {
        if (const auto result = array_field(*card, expected->shear, ordinal, field); !result) {
            return result;
        }
    }
    return check.decimal(source.density, expected->density.number(), Field::Density);
}

} // namespace uvlas::io::detail
