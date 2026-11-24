#include "nastran_lumped_mass_geometry.hpp"

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
using Tensor = std::array<Rational, 6>;
using Code = core::FailureCode;
using Reason = NastranMassGeometryFailure;

[[nodiscard]] std::optional<NastranBulkPosition>
position(const std::optional<NastranBulkLexeme>& field) noexcept {
    return field ? std::optional{field->position} : std::nullopt;
}

[[nodiscard]] Rational scale(ReferenceArithmetic& math, const SiUnits& units,
                             SchemaDimension dimension) {
    math.charge();
    const SiRecipeView recipe{SiDecimalView{1, "1", 0}, {}, units, dimension, 0};
    return partition_take(materialize_rational(recipe, math.leaf_limits, math.plan));
}

[[nodiscard]] Rational to_si(ReferenceArithmetic& math, const NastranMassReal& field,
                             const Rational& factor) {
    if (!field.supplied) {
        return {};
    }
    math.charge();
    const auto value = partition_take(
        Rational::from_decimal(field.supplied->decimal(), math.leaf_limits, math.plan));
    return math.mul(value, factor);
}

[[nodiscard]] Tensor centroidal(ReferenceArithmetic& math, const NastranMassInertia& field,
                                const Rational& factor) {
    return {to_si(math, field.i11, factor),           math.neg(to_si(math, field.i21, factor)),
            to_si(math, field.i22, factor),           math.neg(to_si(math, field.i31, factor)),
            math.neg(to_si(math, field.i32, factor)), to_si(math, field.i33, factor)};
}

[[nodiscard]] Rational minor2(ReferenceArithmetic& math, const Rational& a, const Rational& b,
                              const Rational& cross) {
    return math.sub(math.mul(a, b), math.square(cross));
}

[[nodiscard]] Rational determinant(ReferenceArithmetic& math, const Tensor& v) {
    const auto a = math.mul(v[0], minor2(math, v[2], v[5], v[4]));
    const auto b = math.mul(v[1], math.sub(math.mul(v[1], v[5]), math.mul(v[4], v[3])));
    const auto c = math.mul(v[3], math.sub(math.mul(v[1], v[4]), math.mul(v[2], v[3])));
    return math.add(math.sub(a, b), c);
}

[[nodiscard]] NastranBasicMassKind qualify(ReferenceArithmetic& math, const Tensor& value) {
    if (std::ranges::all_of(value, [](const Rational& field) { return field.sign() == 0; })) {
        return NastranBasicMassKind::Particle;
    }
    if (value[0].sign() <= 0 || minor2(math, value[0], value[2], value[1]).sign() <= 0 ||
        determinant(math, value).sign() <= 0) {
        partition_fail(Code::InvalidInput);
    }

    // The eigenvalues of tr(J)/2 I - J are half the three principal
    // triangle slacks. Test PSD by principal minors using the exact inputs.
    const auto half_trace =
        math.div(math.add(math.add(value[0], value[2]), value[5]), math.constant("2"));
    const Tensor second_moment{math.sub(half_trace, value[0]),
                               math.neg(value[1]),
                               math.sub(half_trace, value[2]),
                               math.neg(value[3]),
                               math.neg(value[4]),
                               math.sub(half_trace, value[5])};
    if (second_moment[0].sign() < 0 || second_moment[2].sign() < 0 || second_moment[5].sign() < 0 ||
        minor2(math, second_moment[0], second_moment[2], second_moment[1]).sign() < 0 ||
        minor2(math, second_moment[0], second_moment[5], second_moment[3]).sign() < 0 ||
        minor2(math, second_moment[2], second_moment[5], second_moment[4]).sign() < 0 ||
        determinant(math, second_moment).sign() < 0) {
        partition_fail(Code::InvalidInput);
    }
    return NastranBasicMassKind::Rigid;
}

[[nodiscard]] NastranCartesianInterval retain(const m::Interval& interval,
                                              const ReferenceArithmetic& math, m::Work& work) {
    return {d::retain(interval.lower, math.leaf_limits, math.plan, work),
            d::retain(interval.upper, math.leaf_limits, math.plan, work)};
}

[[nodiscard]] d::Position grid_position(const NastranShellCards& shells, NastranGridId id,
                                        SiLengthUnit length, const ReferenceArithmetic& math,
                                        m::Work& work, NastranMassGeometryDiagnostic& active) {
    const std::array requested{id};
    std::optional<NastranPositionDiagnostic> diagnostic;
    const auto resolved = NastranPositionEnclosures::prepare(
        shells, length, requested, {math.leaf_limits, 1, work.limit - work.used, work.stop},
        math.plan, &diagnostic);
    if (!resolved) {
        if (diagnostic) {
            active.position = diagnostic->position;
            active.coordinate = diagnostic->coordinate;
        }
        throw g::PartitionAbort{resolved.error()};
    }
    const auto owner = resolved->view();
    if (!owner || owner->grids.size() != 1) {
        partition_fail(Code::InvalidInput);
    }
    work.used += owner->work;
    const auto& grid = owner->grids.front();
    d::Position result;
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = {d::number(grid.metres[i].lower, work).lower,
                     d::number(grid.metres[i].upper, work).upper};
    }
    return result;
}

[[nodiscard]] std::array<NastranCartesianInterval, 6>
shift_to_grid(const Tensor& inertia, const Rational& mass, const d::Position& lever,
              const ReferenceArithmetic& math, m::Work& work) {
    const auto m_kg = d::number(mass, work);
    std::array<NastranCartesianInterval, 6> value;
    std::size_t entry = 0;
    for (std::size_t i = 0; i < lever.size(); ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            m::Interval shifted;
            if (i == j) {
                const auto a = (i + 1) % lever.size(), b = (i + 2) % lever.size();
                shifted = m::add(m::multiply(lever[a], lever[a], work),
                                 m::multiply(lever[b], lever[b], work), work);
            } else {
                shifted = m::subtract({}, m::multiply(lever[i], lever[j], work), work);
            }
            const auto quantity =
                m::add(d::number(inertia[entry], work), m::multiply(m_kg, shifted, work), work);
            value[entry++] = retain(quantity, math, work);
        }
    }
    return value;
}

} // namespace

core::Result<NastranBasicMassGeometry>
prepare_nastran_basic_mass(const NastranShellCards& shells, const NastranLumpedMassCards& masses,
                           NastranElementId id, const SiUnits& units,
                           const NastranMassGeometryLimits& limits, core::ResourcePlan& plan,
                           std::optional<NastranMassGeometryDiagnostic>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    NastranMassGeometryDiagnostic active;
    auto result = g::partition_entry([&]() -> core::Result<NastranBasicMassGeometry> {
        try {
            ReferenceArithmetic math{limits.arithmetic, plan};
            math.cancel();
            const auto shells_view = shells.view();
            const auto mass_view = masses.view();
            if (!shells_view || !mass_view) {
                partition_fail(Code::InvalidInput);
            }
            const auto& deck = shells_view->input.bulk;
            if (!deck.same_snapshot(mass_view->source.location, mass_view->source.digest)) {
                partition_fail(Code::StaleGeneration);
            }
            active = {Reason::Selection, id, {}, {}};
            const auto found = std::ranges::find_if(mass_view->masses, [&](const auto& row) {
                math.visit(1);
                return row.id == id;
            });
            if (found == mass_view->masses.end()) {
                partition_fail(Code::InvalidInput);
            }
            const auto& card = *found;
            if (const auto* framed = std::get_if<NastranFramedMassOffset>(&card.frame)) {
                active = {Reason::Frame, id, position(card.frame_field), framed->frame};
                partition_fail(Code::Unsupported);
            }

            const bool absolute = std::holds_alternative<NastranBasicMassCentre>(card.frame);
            if (absolute) {
                for (const auto& field : card.position) {
                    if (!field.supplied) {
                        active = {Reason::Position,
                                  id,
                                  position(field.field).value_or(card.card.position),
                                  {}};
                        partition_fail(Code::Unsupported);
                    }
                }
            }

            // Live SI quantities, exact tensor minors and the result endpoint rationals.
            // Nested GRID preparation receives the remaining external estimate.
            math.admit(72);
            active = {Reason::Mass, id, position(card.mass.field), {}};
            const auto mass_factor = scale(math, units, {{0, 1, 0, 0}});
            auto mass_kg = to_si(math, card.mass, mass_factor);
            if (mass_kg.sign() <= 0) {
                partition_fail(Code::InvalidInput);
            }
            active = {Reason::Inertia,
                      id,
                      position(card.inertia.i11.field).value_or(card.card.position),
                      {}};
            const auto inertia_factor = scale(math, units, {{2, 1, 0, 0}});
            auto inertia = centroidal(math, card.inertia, inertia_factor);
            const auto kind = qualify(math, inertia);

            m::Work work{0, limits.maximum_directed_work, limits.arithmetic.cancellation};
            active = {Reason::Grid, id, card.grid_field.position, {}};
            const auto grid = grid_position(shells, card.grid, units.length, math, work, active);
            const auto length = scale(math, units, {{1, 0, 0, 0}});
            std::array<NastranCartesianInterval, 3> grid_metres, centre_metres,
                grid_to_centre_metres;
            d::Position lever;
            for (std::size_t i = 0; i < grid.size(); ++i) {
                active = {Reason::Arithmetic,
                          id,
                          position(card.position[i].field).value_or(card.card.position),
                          {}};
                const auto input = d::number(to_si(math, card.position[i], length), work);
                const auto centre = absolute ? input : m::add(grid[i], input, work);
                lever[i] = absolute ? m::subtract(input, grid[i], work) : input;
                grid_metres[i] = retain(grid[i], math, work);
                centre_metres[i] = retain(centre, math, work);
                grid_to_centre_metres[i] = retain(lever[i], math, work);
            }
            active = {Reason::Arithmetic, id, card.card.position, {}};
            auto shifted = shift_to_grid(inertia, mass_kg, lever, math, work);
            math.cancel();
            return NastranBasicMassGeometry{mass_view->source,
                                            card,
                                            units,
                                            kind,
                                            std::move(mass_kg),
                                            std::move(grid_metres),
                                            std::move(centre_metres),
                                            std::move(grid_to_centre_metres),
                                            std::move(inertia),
                                            std::move(shifted),
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
