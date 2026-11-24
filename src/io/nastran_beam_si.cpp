#include "nastran_beam_si.hpp"

#include "../structure/pg_reference_arithmetic.hpp"
#include "nastran_position_arithmetic.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace uvlas::io::detail {
namespace {

namespace g = structure::detail;
namespace d = nastran_directed;
namespace m = d::m;
using g::partition_fail;
using g::partition_take;
using g::ReferenceArithmetic;
using Rational = structure::ExactRational;
using Code = core::FailureCode;
using Reason = NastranBeamSiFailure;
using Axes = std::array<d::Position, 3>;
using Section = std::array<Rational, 6>;
constexpr auto operation = core::Operation::PrepareModel;

[[nodiscard]] std::optional<NastranBulkPosition>
position(const std::optional<NastranBulkLexeme>& field) noexcept {
    return field ? std::optional{field->position} : std::nullopt;
}

struct BeamContext {
    ReferenceArithmetic& math;
    m::Work& work;
    NastranBeamSiDiagnostic& active;

    void stage(Reason reason, std::optional<NastranBulkPosition> field) const noexcept {
        active.reason = reason;
        active.position = field;
        active.grid.reset();
        active.coordinate.reset();
    }

    [[noreturn]] void
    forward_position_failure(const core::Failure& failure,
                             const std::optional<NastranPositionDiagnostic>& details) const {
        if (details) {
            if (details->position) {
                active.position = details->position;
            }
            active.grid = details->grid;
            active.coordinate = details->coordinate;
        }
        throw g::PartitionAbort{failure};
    }
};

[[nodiscard]] Rational scale(ReferenceArithmetic& math, const SiUnits& units,
                             SchemaDimension dimension) {
    math.charge();
    const SiRecipeView recipe{SiDecimalView{1, "1", 0}, {}, units, dimension, 0};
    return partition_take(materialize_rational(recipe, math.leaf_limits, math.plan));
}

[[nodiscard]] Rational real(ReferenceArithmetic& math, const NastranShellReal& field) {
    math.charge();
    return partition_take(
        Rational::from_decimal(field.value.decimal(), math.leaf_limits, math.plan));
}

[[nodiscard]] Rational dimensional(ReferenceArithmetic& math, const NastranPbeamReal& field,
                                   const Rational& unit) {
    return field ? math.mul(real(math, *field), unit) : Rational{};
}

void reject_end_options(const BeamContext& context, const NastranBeamEnd& end) {
    if (end.releases.element_components != 0) {
        context.stage(Reason::Release, position(end.releases.field));
        partition_fail(Code::Unsupported);
    }
    if (end.warping && end.warping->grid_ordinal) {
        context.stage(Reason::Warping, end.warping->field.position);
        partition_fail(Code::Unsupported);
    }
}

[[nodiscard]] m::Interval interval(const NastranCartesianInterval& input, m::Work& work) {
    return {d::number(input.lower, work).lower, d::number(input.upper, work).upper};
}

[[nodiscard]] d::Position intervals(const NastranSiBeamVector& point, m::Work& work) {
    return {interval(point[0], work), interval(point[1], work), interval(point[2], work)};
}

[[nodiscard]] NastranCartesianInterval retain(const m::Interval& value, const BeamContext& c) {
    return {d::retain(value.lower, c.math.leaf_limits, c.math.plan, c.work),
            d::retain(value.upper, c.math.leaf_limits, c.math.plan, c.work)};
}

[[nodiscard]] NastranSiBeamVector retain(const d::Position& value, const BeamContext& c) {
    return {retain(value[0], c), retain(value[1], c), retain(value[2], c)};
}

[[nodiscard]] d::Position add(const d::Position& a, const d::Position& b, m::Work& work) {
    return {m::add(a[0], b[0], work), m::add(a[1], b[1], work), m::add(a[2], b[2], work)};
}

[[nodiscard]] Axes basis(const d::Position& tangent, const d::Position& vector, m::Work& work) {
    const auto x = d::normalized(tangent, work);
    const auto z = d::normalized(d::cross(x, vector, work), work);
    return {x, d::cross(z, x, work), z};
}

struct BeamVector {
    d::Position values;
    std::uint8_t active = 0;
};

[[nodiscard]] BeamVector read_vector(const NastranBeamVector& values, const Rational& unit,
                                     const BeamContext& c) {
    BeamVector result;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (const auto& field = values[i]) {
            c.active.position = field->field.position;
            const auto value = c.math.mul(real(c.math, *field), unit);
            result.values[i] = d::number(value, c.work);
            if (value.sign() != 0) {
                result.active = static_cast<std::uint8_t>(result.active | (1U << i));
            }
        }
    }
    return result;
}

void add_basis_vector(d::Position& result, const d::Position& axis, const m::Interval& coefficient,
                      m::Work& work) {
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = m::add(result[i], m::multiply(axis[i], coefficient, work), work);
    }
}

[[nodiscard]] d::Position to_basic(const BeamVector& vector, NastranGridId grid,
                                   const NastranIsotropicShellCards& input, const BeamContext& c) {
    if (vector.active == 0) {
        return {};
    }
    std::optional<NastranPositionDiagnostic> diagnostic;
    const auto axes = prepare_nastran_displacement_axes(
        input, grid, vector.active, {c.math.leaf_limits, c.work.limit - c.work.used, c.work.stop},
        c.math.plan, &diagnostic);
    if (!axes) {
        c.forward_position_failure(axes.error(), diagnostic);
    }
    c.work.used += axes->work;
    d::Position result;
    for (std::size_t i = 0; i < axes->directions.size(); ++i) {
        if ((vector.active & (1U << i)) != 0) {
            add_basis_vector(result, intervals(axes->directions[i], c.work), vector.values[i],
                             c.work);
        }
    }
    return result;
}

[[nodiscard]] d::Position resolve_offset(const BeamVector& raw, const NastranBeamEnd& end,
                                         NastranBeamOffsetBasis interpretation,
                                         const std::optional<Axes>& offset_axes,
                                         const NastranIsotropicShellCards& input,
                                         const BeamContext& c) {
    if (raw.active == 0) {
        return {};
    }
    if (interpretation == NastranBeamOffsetBasis::GridDisplacement) {
        return to_basic(raw, end.point.id, input, c);
    }
    if (!offset_axes) {
        partition_fail(Code::InvalidInput);
    }
    d::Position result;
    for (std::size_t i = 0; i < offset_axes->size(); ++i) {
        add_basis_vector(result, (*offset_axes)[i], raw.values[i], c.work);
    }
    return result;
}

struct PhysicalAxis {
    NastranSiBeamEnd a, b;
    std::array<NastranSiBeamVector, 3> axes;
    NastranCartesianInterval length;
};

[[nodiscard]] PhysicalAxis physical_axis(const NastranIsotropicShellCards& input,
                                         const NastranCbeam& beam, SiLengthUnit units,
                                         const Rational& length_scale, const BeamContext& c) {
    std::array points{beam.a.point.id, beam.b.point.id, beam.a.point.id};
    const auto* g0 = std::get_if<NastranBeamGridReference>(&beam.orientation);
    if (g0) {
        points[2] = g0->id;
    }
    c.stage(Reason::Axis, beam.a.point.field.position);
    std::optional<NastranPositionDiagnostic> diagnostic;
    const auto positions = NastranPositionEnclosures::prepare(
        input, units, std::span{points}.first(g0 ? 3 : 2),
        {c.math.leaf_limits, 3, c.work.limit - c.work.used, c.work.stop}, c.math.plan, &diagnostic);
    if (!positions) {
        c.forward_position_failure(positions.error(), diagnostic);
    }
    const auto view = positions->view();
    if (!view || view->grids.size() != (g0 ? 3U : 2U)) {
        partition_fail(Code::InvalidInput);
    }
    c.work.used += view->work;
    const auto grid_a = intervals(view->grids[0].metres, c.work);
    const auto grid_b = intervals(view->grids[1].metres, c.work);
    auto angle_field = beam.id_field.position;
    if (g0) {
        angle_field = g0->field.position;
    } else {
        const auto& raw = std::get<NastranBeamVector>(beam.orientation);
        const auto source =
            std::ranges::find_if(raw, [](const auto& item) { return item.has_value(); });
        if (source != raw.end() && *source) {
            angle_field = (*source)->field.position;
        }
    }
    c.stage(Reason::Orientation, angle_field);
    d::Position direction;
    if (g0) {
        direction = d::subtract(intervals(view->grids[2].metres, c.work), grid_a, c.work);
    } else {
        const auto vector =
            read_vector(std::get<NastranBeamVector>(beam.orientation), c.math.constant("1"), c);
        direction = beam.interpretation.orientation == NastranBeamOrientationBasis::Basic
                        ? vector.values
                        : to_basic(vector, beam.a.point.id, input, c);
    }
    const auto a_offset = read_vector(beam.a.offset, length_scale, c);
    const auto b_offset = read_vector(beam.b.offset, length_scale, c);
    std::optional<Axes> offset_axes;
    if ((a_offset.active != 0 && beam.interpretation.a == NastranBeamOffsetBasis::OffsetSystem) ||
        (b_offset.active != 0 && beam.interpretation.b == NastranBeamOffsetBasis::OffsetSystem)) {
        c.stage(Reason::Orientation, angle_field);
        offset_axes = basis(d::subtract(grid_b, grid_a, c.work), direction, c.work);
    }
    c.stage(Reason::EndOffset, beam.a.point.field.position);
    const auto da = resolve_offset(a_offset, beam.a, beam.interpretation.a, offset_axes, input, c);
    c.stage(Reason::EndOffset, beam.b.point.field.position);
    const auto db = resolve_offset(b_offset, beam.b, beam.interpretation.b, offset_axes, input, c);
    const auto end_a = add(grid_a, da, c.work);
    const auto end_b = add(grid_b, db, c.work);
    c.stage(Reason::Axis, beam.b.point.field.position);
    const auto tangent = d::subtract(end_b, end_a, c.work);
    const auto x = d::normalized(tangent, c.work);
    const auto length = m::square_root(d::norm_squared(tangent, c.work), c.work);
    c.stage(Reason::Orientation, angle_field);
    const auto z = d::normalized(d::cross(x, direction, c.work), c.work);
    return {{retain(grid_a, c), retain(da, c), retain(end_a, c)},
            {retain(grid_b, c), retain(db, c), retain(end_b, c)},
            {retain(x, c), retain(d::cross(z, x, c.work), c), retain(z, c)},
            retain(length, c)};
}

[[nodiscard]] NastranSiBeamMaterial material(const NastranMat1& row, NastranBulkPosition source,
                                             const SiUnits& units, const BeamContext& c) {
    c.stage(Reason::Material, source);
    const unsigned count = static_cast<unsigned>(row.young.has_value()) +
                           static_cast<unsigned>(row.shear.has_value()) +
                           static_cast<unsigned>(row.poisson.has_value());
    if (count < 2) {
        partition_fail(Code::Unsupported);
    }
    auto& math = c.math;
    const auto pressure = scale(math, units, {{-1, 1, -2, 0}});
    auto e = dimensional(math, row.young, pressure);
    auto g_modulus = dimensional(math, row.shear, pressure);
    auto nu = row.poisson ? real(math, *row.poisson) : Rational{};
    if ((row.young && e.sign() <= 0) || (row.shear && g_modulus.sign() <= 0)) {
        partition_fail(Code::InvalidInput);
    }
    const auto one = math.constant("1"), two = math.constant("2");
    if (row.poisson) {
        c.active.position = row.poisson->field.position;
        const auto factor = math.mul(two, math.add(one, nu));
        if (factor.sign() <= 0 || math.compare(factor, math.constant("3")) >= 0) {
            partition_fail(Code::InvalidInput);
        }
        if (!row.young) {
            e = math.mul(g_modulus, factor);
        } else if (!row.shear) {
            g_modulus = math.div(e, factor);
        } else if (math.compare(e, math.mul(g_modulus, factor)) != 0) {
            partition_fail(Code::InvalidInput);
        }
    } else {
        const auto ratio = math.div(e, math.mul(two, g_modulus));
        if (ratio.sign() <= 0 || math.compare(ratio, math.div(math.constant("3"), two)) >= 0) {
            partition_fail(Code::InvalidInput);
        }
        nu = math.sub(ratio, one);
    }
    std::optional<Rational> density;
    if (const auto& field = row.density) {
        c.active.position = field->field.position;
        density = dimensional(math, row.density, scale(math, units, {{-3, 1, 0, 0}}));
        if (density->sign() < 0) {
            partition_fail(Code::InvalidInput);
        }
    }
    return {row, std::move(e), std::move(g_modulus), std::move(nu), std::move(density)};
}

struct PropertyScales {
    Rational area, fourth_moment, line_mass, mass_moment;

    PropertyScales(const SiUnits& units, ReferenceArithmetic& math)
        : area(scale(math, units, {{2, 0, 0, 0}})),
          fourth_moment(scale(math, units, {{4, 0, 0, 0}})),
          line_mass(scale(math, units, {{-1, 1, 0, 0}})),
          mass_moment(scale(math, units, {{1, 1, 0, 0}})) {}
};

struct SectionSource {
    std::array<const NastranPbeamReal*, 6> fields;
    std::array<const Rational*, 6> units;

    SectionSource(const NastranPbeamSection& raw, const PropertyScales& scale) noexcept
        : fields{&raw.area, &raw.i1, &raw.i2, &raw.i12, &raw.torsion, &raw.nonstructural_mass},
          units{&scale.area,          &scale.fourth_moment, &scale.fourth_moment,
                &scale.fourth_moment, &scale.fourth_moment, &scale.line_mass} {}
};

[[nodiscard]] Section section_values(const NastranPbeamSection& raw, const PropertyScales& scale,
                                     const Section* defaults, const BeamContext& c) {
    const SectionSource source{raw, scale};
    Section result;
    for (std::size_t i = 0; i < result.size(); ++i) {
        const auto& field = *source.fields[i];
        if (field) {
            c.active.position = field->field.position;
            result[i] = dimensional(c.math, field, *source.units[i]);
        } else if (defaults) {
            result[i] = c.math.copy((*defaults)[i]);
        }
    }
    return result;
}

[[nodiscard]] Section interpolate(const Section& a, const Section& b, const Rational& s,
                                  ReferenceArithmetic& math) {
    Section result;
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = math.add(a[i], math.mul(s, math.sub(b[i], a[i])));
    }
    return result;
}

void reject_nonzero(const NastranPbeamReal& field, Reason reason, const BeamContext& c) {
    if (field) {
        c.stage(reason, field->field.position);
        if (real(c.math, *field).sign() != 0) {
            partition_fail(Code::Unsupported);
        }
    }
}

[[nodiscard]] Rational shear_factor(const NastranPbeamReal& field, const BeamContext& c) {
    if (!field) {
        return c.math.constant("1");
    }
    c.stage(Reason::Modifier, field->field.position);
    auto factor = real(c.math, *field);
    if (factor.sign() <= 0) {
        // Zero has a specific MSC meaning: zero flexibility, not zero stiffness.
        partition_fail(factor.sign() == 0 ? Code::Unsupported : Code::InvalidInput);
    }
    return factor;
}

void reject_recovery(const NastranPbeam& property, std::span<const NastranPbeamStation> stations,
                     const BeamContext& c) {
    if (property.recovery_a) {
        for (const auto& field : property.recovery_a->locations) {
            if (field) {
                c.stage(Reason::Recovery, field->field.position);
                partition_fail(Code::Unsupported);
            }
        }
    }
    for (const auto& station : stations) {
        c.math.visit(1);
        if (station.output != NastranPbeamOutput::No) {
            c.stage(Reason::Recovery, station.output_field.position);
            partition_fail(Code::Unsupported);
        }
    }
}

[[nodiscard]] NastranSiBeamSectionOffsets read_section_offsets(const NastranPbeam& property,
                                                               const Rational& length,
                                                               const PropertyScales& scale,
                                                               const BeamContext& c) {
    NastranSiBeamSectionOffsets result;
    if (const auto& modifier = property.modifiers) {
        reject_nonzero(modifier->warping[0], Reason::Warping, c);
        reject_nonzero(modifier->warping[1], Reason::Warping, c);
        reject_nonzero(modifier->taper_relief_1, Reason::Modifier, c);
        reject_nonzero(modifier->taper_relief_2, Reason::Modifier, c);
        for (std::size_t i = 0; i < modifier->nonstructural_inertia.size(); ++i) {
            const auto& field = modifier->nonstructural_inertia[i];
            if (field) {
                c.stage(Reason::Modifier, field->field.position);
                result.nonstructural_rotary_kg_m[i] = dimensional(c.math, field, scale.mass_moment);
                if (result.nonstructural_rotary_kg_m[i].sign() < 0) {
                    partition_fail(Code::InvalidInput);
                }
            } else if (i == 1) {
                result.nonstructural_rotary_kg_m[1] =
                    c.math.copy(result.nonstructural_rotary_kg_m[0]);
            }
        }
    }
    if (const auto& offsets = property.offsets) {
        const auto pair = [&](const std::array<NastranPbeamReal, 4>& fields,
                              std::array<Rational, 4>& values) {
            for (std::size_t i = 0; i < fields.size(); ++i) {
                if (fields[i]) {
                    c.stage(Reason::Modifier, fields[i]->field.position);
                    values[i] = dimensional(c.math, fields[i], length);
                } else if (i >= 2) {
                    values[i] = c.math.copy(values[i - 2]);
                }
            }
        };
        pair(offsets->nonstructural_mass_center, result.nonstructural_mass_centre_metres);
        pair(offsets->neutral_axis, result.neutral_axis_metres);
    }
    return result;
}

[[nodiscard]] NastranSiBeamStation
projected_station(Rational&& fraction, const NastranPbeamSection& raw,
                  std::optional<NastranBulkPosition> field, Section&& section,
                  const NastranSiBeamMaterial& elastic, const std::array<Rational, 2>& factors,
                  const BeamContext& c) {
    auto& math = c.math;
    const auto reject = [&](bool invalid, const NastranPbeamReal& entry) {
        if (invalid) {
            c.stage(Reason::Section, entry ? std::optional{entry->field.position} : field);
            partition_fail(Code::InvalidInput);
        }
    };
    reject(section[0].sign() <= 0, raw.area);
    reject(section[1].sign() <= 0, raw.i1);
    reject(section[2].sign() <= 0, raw.i2);
    reject(section[4].sign() < 0, raw.torsion);
    reject(section[5].sign() < 0, raw.nonstructural_mass);
    reject(math.compare(math.mul(section[1], section[2]), math.square(section[3])) <= 0, raw.i12);
    std::optional<Rational> mass;
    if (const auto& density = elastic.density_kg_per_m3) {
        mass = math.add(math.mul(*density, section[0]), section[5]);
    }
    const auto ga = math.mul(elastic.shear_pa, section[0]);
    NastranSiBeamStiffnessScalars stiffness{math.mul(elastic.young_pa, section[0]),
                                            math.mul(ga, factors[0]),
                                            math.mul(ga, factors[1]),
                                            math.mul(elastic.shear_pa, section[4]),
                                            math.mul(elastic.young_pa, section[1]),
                                            math.mul(elastic.young_pa, section[2]),
                                            math.mul(elastic.young_pa, section[3])};
    return {std::move(fraction),
            raw,
            field,
            {std::move(section[0]), std::move(section[1]), std::move(section[2]),
             std::move(section[3]), std::move(section[4]), std::move(section[5])},
            std::move(stiffness),
            std::move(mass)};
}

[[nodiscard]] core::PlannedBuffer<NastranSiBeamStation>
sections(const NastranPbeam& property, std::span<const NastranPbeamStation> raw,
         const PropertyScales& scale, const NastranSiBeamMaterial& elastic, const BeamContext& c) {
    auto& math = c.math;
    const auto* modifiers = property.modifiers ? &*property.modifiers : nullptr;
    const std::array factors{
        shear_factor(modifiers ? modifiers->shear_factor_1 : NastranPbeamReal{}, c),
        shear_factor(modifiers ? modifiers->shear_factor_2 : NastranPbeamReal{}, c)};
    auto rows = partition_take(
        math.plan.allocate<NastranSiBeamStation>(g::partition_index(raw.size() + 1), operation));
    c.stage(Reason::Section, property.declaration.id_field.position);
    const auto a = section_values(property.end_a, scale, nullptr, c);
    Section end_a;
    std::ranges::transform(a, end_a.begin(),
                           [&](const Rational& value) { return math.copy(value); });
    rows.values()[0] =
        projected_station({}, property.end_a, {}, std::move(end_a), elastic, factors, c);
    if (raw.empty()) {
        return rows;
    }
    if (!property.end_b_station || *property.end_b_station >= raw.size()) {
        partition_fail(Code::InvalidInput);
    }
    const auto end_b = *property.end_b_station;
    const auto b = section_values(raw[end_b].section, scale, &a, c);
    for (std::size_t i = 0; i < raw.size(); ++i) {
        math.visit(1);
        const auto& station = raw[i];
        c.stage(Reason::Section, station.output_field.position);
        if (!station.fraction) {
            partition_fail(Code::InvalidInput);
        }
        const auto field = station.fraction->field.position;
        c.active.position = field;
        auto fraction = real(math, *station.fraction);
        const auto inherited = interpolate(a, b, fraction, math);
        auto values = section_values(station.section, scale, &inherited, c);
        rows.values()[i + 1] = projected_station(std::move(fraction), station.section, field,
                                                 std::move(values), elastic, factors, c);
    }
    return rows;
}

} // namespace

core::Result<NastranSiBeamInput>
prepare_nastran_si_beam(const NastranIsotropicShellCards& input, const NastranBeamCards& beams,
                        const NastranPbeamCards& props, NastranElementId id, const SiUnits& units,
                        const NastranBeamSiLimits& limits, core::ResourcePlan& plan,
                        std::optional<NastranBeamSiDiagnostic>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    NastranBeamSiDiagnostic active;
    auto result = g::partition_entry([&]() -> core::Result<NastranSiBeamInput> {
        try {
            ReferenceArithmetic math{limits.arithmetic, plan};
            m::Work work{0, limits.maximum_directed_work, limits.arithmetic.cancellation};
            const BeamContext context{math, work, active};
            math.cancel();
            const auto model = input.view();
            const auto elements = beams.view();
            const auto properties = props.view();
            if (!model || !elements || !properties) {
                partition_fail(Code::InvalidInput);
            }
            const auto& deck = model->input.input.bulk;
            if (!deck.same_snapshot(elements->source.location, elements->source.digest) ||
                !deck.same_snapshot(properties->source.location, properties->source.digest)) {
                partition_fail(Code::StaleGeneration);
            }

            active.beam = id;
            context.stage(Reason::Selection, {});
            const auto row = std::ranges::find_if(elements->beams, [&](const auto& beam) {
                math.visit(1);
                return beam.id == id;
            });
            if (row == elements->beams.end()) {
                partition_fail(Code::InvalidInput);
            }
            const auto& beam = *row;
            const auto ordinal = beam.property.declaration_ordinal;
            if (ordinal >= properties->properties.size() ||
                beam.property.id != properties->properties[ordinal].declaration.id) {
                partition_fail(Code::InvalidInput);
            }
            const auto& property = properties->properties[ordinal];
            const auto stations = properties->sections(property);
            if (stations.size() != property.station_count) {
                partition_fail(Code::InvalidInput);
            }
            if (stations.size() + 1 > limits.maximum_stations) {
                partition_fail(Code::ResourceExceeded);
            }
            reject_end_options(context, beam.a);
            reject_end_options(context, beam.b);
            reject_recovery(property, stations, context);

            const auto material_id = property.material.declaration_ordinal;
            if (material_id >= model->materials.size()) {
                partition_fail(Code::InvalidInput);
            }
            const auto& declaration = model->materials[material_id];
            if (!declaration.mat1_ordinal || *declaration.mat1_ordinal >= model->mat1.size() ||
                declaration.id != property.material.id) {
                partition_fail(Code::InvalidInput);
            }
            const auto slots = g::partition_add(128, g::partition_mul(16, stations.size() + 1));
            math.admit(slots);
            const auto length_scale = scale(math, units, {{1, 0, 0, 0}});
            auto elastic = material(model->mat1[*declaration.mat1_ordinal],
                                    declaration.id_field.position, units, context);
            const PropertyScales scales{units, math};
            auto offsets = read_section_offsets(property, length_scale, scales, context);
            auto samples = sections(property, stations, scales, elastic, context);
            auto geometry = physical_axis(input, beam, units.length, length_scale, context);
            math.cancel();
            return NastranSiBeamInput{elements->source,
                                      beam,
                                      property,
                                      units,
                                      std::move(elastic),
                                      std::move(offsets),
                                      std::move(geometry.a),
                                      std::move(geometry.b),
                                      std::move(geometry.axes),
                                      std::move(geometry.length),
                                      std::move(samples),
                                      math.info(),
                                      work.used};
        } catch (const m::ArithmeticFailure& failure) {
            partition_fail(failure.code);
        } catch (const std::overflow_error&) {
            partition_fail(Code::PrecisionUnresolved);
        } catch (const std::range_error&) {
            partition_fail(Code::PrecisionUnresolved);
        }
    });
    if (!result && diagnostic && result.error().code != Code::Cancelled) {
        *diagnostic = active;
    }
    return result;
}

} // namespace uvlas::io::detail
