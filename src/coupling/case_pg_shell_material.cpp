#include "case_pg_shell_material.hpp"

#include "aero/certified_trigonometry_private.hpp"
#include "case_pg_chart.hpp"
#include "io/si_case_detail.hpp"
#include "structure/pg_rank_storage.hpp"
#include "structure/pg_reference_arithmetic.hpp"

#include <tuple>

namespace uvlas::coupling::detail {
namespace {
namespace pg = structure;
namespace phase = aero::detail::p1;
using core::FailureCode;
using pg::detail::partition_entry;
using pg::detail::partition_fail;
using pg::detail::partition_index;
using pg::detail::partition_take;
using pg::detail::ReferenceArithmetic;
constexpr auto operation = core::Operation::PrepareModel;
constexpr io::SchemaDimension pressure{{-1, 1, -2, 0}};
constexpr io::SchemaDimension density{{-3, 1, 0, 0}};
constexpr io::SchemaDimension dimensionless{{0, 0, 0, 0}};
constexpr io::SchemaDimension length{{1, 0, 0, 0}};
constexpr io::SchemaDimension angle{{0, 0, 0, 1}};

template <class Field>
pg::ExactRational quantity(const io::SiCase& source, Field field, io::SchemaDimension dimension,
                           ReferenceArithmetic& math, bool pi_multiple = false) {
    math.charge();
    auto recipe = partition_take(source.quantity(field)).recipe;
    if (recipe.dimensions.exponents != dimension.exponents ||
        recipe.pi_power != static_cast<int>(pi_multiple))
        partition_fail(FailureCode::InvalidInput);
    if (pi_multiple) {
        recipe.units.angle = io::SiAngleUnit::Radian;
        recipe.pi_power = 0;
        const auto degrees =
            partition_take(io::materialize_rational(recipe, math.leaf_limits, math.plan));
        return math.div(degrees, math.constant("180"));
    }
    return partition_take(io::materialize_rational(recipe, math.leaf_limits, math.plan));
}

template <class Field>
std::array<pg::ExactRational, 3> triple(const io::SiCase& source, Field field,
                                        io::SchemaDimension dimension, ReferenceArithmetic& math) {
    std::array<pg::ExactRational, 3> values;
    for (std::size_t i = 0; i < values.size(); ++i)
        values[i] = quantity(source, partition_take(field.at(partition_index(i))), dimension, math);
    return values;
}

// Worst-case owned rational slots per reusable law (including both Mandel coefficients).
constexpr std::uint64_t law_rationals = 2 * 36 + 2;

pg::PgShellElasticLaw prepare_law(const io::SiCase& source, const auto& record,
                                  std::uint64_t other_rationals, ReferenceArithmetic& math) {
    using Kind = std::decay_t<decltype(record)>;
    const auto rho = quantity(source, record.density, density, math);
    const auto qualify = [&](const pg::PgElasticity& input) {
        auto remaining = math.limits;
        remaining.maximum_scalar_operations -= math.scalar_operations;
        math.reserve_nested_rationals(remaining.algebra, other_rationals);
        auto law = partition_take(pg::PgShellElasticLaw::prepare(input, rho, remaining, math.plan));
        math.charge(law.call_info().scalar_operations);
        return law;
    };
    if constexpr (std::is_same_v<Kind, io::schema::MaterialIsotropic>) {
        const auto young = quantity(source, record.young, pressure, math);
        const auto poisson = quantity(source, record.poisson, dimensionless, math);
        return qualify(pg::PgIsotropicElasticity{young, poisson});
    } else if constexpr (std::is_same_v<Kind, io::schema::MaterialOrthotropic>) {
        const auto young = triple(source, record.young, pressure, math);
        const auto poisson = triple(source, record.poisson, dimensionless, math);
        const auto shear = triple(source, record.shear, pressure, math);
        return qualify(pg::PgOrthotropicElasticity{young, poisson, shear});
    } else {
        std::array<pg::ExactRational, 36> stiffness;
        for (core::Index i = 0; i < 6; ++i) {
            const auto row = partition_take(record.stiffness.at(i));
            for (core::Index j = 0; j < 6; ++j) {
                stiffness[static_cast<std::size_t>(6 * i + j)] =
                    quantity(source, partition_take(row.at(j)), pressure, math);
            }
        }
        const auto tag = record.components.selected();
        const auto convention = tag == decltype(tag)::ItemEngineering
                                    ? pg::PgMaterialComponents::Engineering
                                    : pg::PgMaterialComponents::Mandel;
        return qualify(pg::PgEffectiveElasticity{stiffness, convention});
    }
}

pg::PgReferenceScalar integer_power(pg::detail::ReferenceRadicals& field,
                                    const pg::PgReferenceScalar& base, std::uint64_t exponent,
                                    const pg::ExactRational& one) {
    pg::PgReferenceScalar result{field.math.copy(one), {}};
    if (exponent == 0)
        return result;
    auto factor = field.scale(base, one);
    while (exponent != 0) {
        if (exponent % 2 != 0)
            result = field.mul(result, factor);
        exponent /= 2;
        if (exponent != 0)
            factor = field.mul(factor, factor);
    }
    return result;
}

template <std::size_t N>
void sample_polynomial(const io::SiCase& source, const io::schema::ReferenceFieldPolynomial& input,
                       const CasePgMaterialCoordinates& point,
                       std::array<pg::PgReferenceScalar, N>& values,
                       std::array<pg::ExactRational, N>& errors, io::SchemaDimension dimension,
                       bool pi_multiple, ReferenceArithmetic& math) {
    const auto coordinates = input.coordinates.selected();
    using Coordinates = decltype(input.coordinates.selected());
    if (coordinates == Coordinates::ItemMaterial && !point.has_material_metres)
        partition_fail(FailureCode::InvalidInput);
    const auto measure = coordinates == Coordinates::ItemChart ? dimensionless : length;
    const auto origin = triple(source, input.origin, measure, math);
    const auto lengths = triple(source, input.lengths, measure, math);
    const auto one = math.constant("1");
    pg::detail::ReferenceRadicals field{math, point.g};
    pg::PgReferenceVector normalized;
    for (std::size_t i = 0; i < normalized.size(); ++i) {
        if (lengths[i].sign() <= 0)
            partition_fail(FailureCode::InvalidInput);
        const auto scale = math.div(one, lengths[i]);
        if (coordinates == Coordinates::ItemChart) {
            normalized[i].rational = math.mul(math.sub(point.chart[i], origin[i]), scale);
        } else {
            const auto& position = coordinates == Coordinates::ItemMaterial
                                       ? point.material_metres[i]
                                       : point.world_reference_metres[i];
            normalized[i].rational = math.mul(math.sub(position.rational, origin[i]), scale);
            normalized[i].radical = math.mul(position.radical, scale);
        }
    }
    for (std::size_t i = 0; i < N; ++i) {
        const auto ordinal = partition_index(i);
        errors[i] =
            quantity(source, partition_take(input.error.at(ordinal)), dimension, math, pi_multiple);
        if (errors[i].sign() < 0)
            partition_fail(FailureCode::InvalidInput);
        const auto terms = partition_take(input.components.at(ordinal));
        for (core::Index t = 0; t < partition_index(terms.size()); ++t) {
            math.visit(1);
            const auto term = partition_take(terms.at(t)).record();
            pg::PgReferenceScalar value{
                quantity(source, term.coefficient, dimension, math, pi_multiple), {}};
            for (std::size_t axis = 0; axis < normalized.size(); ++axis) {
                const auto power = partition_take(term.powers.at(partition_index(axis)));
                value = field.mul(
                    value, integer_power(field, normalized[axis], power.unsigned_value(), one));
            }
            values[i] = field.add(values[i], value);
        }
    }
}

void sample_expression(const io::SiCase&, const io::schema::ReferenceFieldExpression&,
                       const CasePgMaterialCoordinates&, std::span<pg::PgReferenceScalar>,
                       std::span<pg::ExactRational>, io::SchemaDimension, bool,
                       ReferenceArithmetic&);

template <std::size_t N>
void sample_field(const io::SiCase& source, const auto& field,
                  const CasePgMaterialCoordinates& point,
                  std::array<pg::PgReferenceScalar, N>& values,
                  std::array<pg::ExactRational, N>& errors, io::SchemaDimension dimension,
                  bool pi_multiple, ReferenceArithmetic& math) {
    if (const auto* constant = std::get_if<io::schema::ReferenceFieldConstant>(&field)) {
        for (std::size_t i = 0; i < N; ++i)
            values[i].rational =
                quantity(source, partition_take(constant->value.at(partition_index(i))), dimension,
                         math, pi_multiple);
    } else if (const auto* polynomial = std::get_if<io::schema::ReferenceFieldPolynomial>(&field)) {
        sample_polynomial(source, *polynomial, point, values, errors, dimension, pi_multiple, math);
    } else {
        sample_expression(source, std::get<io::schema::ReferenceFieldExpression>(field), point,
                          values, errors, dimension, pi_multiple, math);
    }
}

phase::Interval phase_rational(const pg::ExactRational& value, phase::Work& work) {
    work.step();
    constexpr auto bits = static_cast<std::uint64_t>(std::numeric_limits<phase::Integer>::digits);
    if (value.numerator_bits() > bits || value.denominator_bits() > bits)
        partition_fail(FailureCode::ResourceExceeded);
    const auto terms = pg::detail::rank_fraction<phase::Integer>(value);
    const auto scalar = [](const phase::Integer& x) {
        return phase::Interval{phase::rounded(x, 0, false), phase::rounded(x, 0, true)};
    };
    return phase::divide(scalar(terms.numerator), scalar(terms.denominator), work);
}

pg::ExactRational phase_endpoint(const phase::Bound& value, ReferenceArithmetic& math) {
    if (value.mantissa == 0)
        return {};
    math.charge();
    const auto power = static_cast<std::int64_t>(value.power);
    constexpr std::int64_t wide_bits = 8192;
    const auto magnitude_bits = static_cast<std::int64_t>(phase::bits(value.mantissa));
    if (power <= -wide_bits || power > wide_bits - magnitude_bits)
        partition_fail(FailureCode::PrecisionUnresolved);
    using Wide = pg::detail::RankWide;
    Wide units{value.mantissa};
    if (power < 0)
        return pg::detail::rank_materialize_dyadic(units, static_cast<std::uint32_t>(-power),
                                                   math.leaf_limits, math.plan);
    units = pg::detail::rank_abs(units) << static_cast<std::uint32_t>(power);
    if (value.mantissa < 0)
        units = -units;
    return pg::detail::rank_exact(pg::detail::RankFraction<Wide>{std::move(units), Wide{1}},
                                  math.leaf_limits, math.plan);
}

phase::Interval phase_radical(const pg::PgReferenceScalar& value, const pg::ExactRational& g,
                              phase::Work& work) {
    auto base = phase_rational(value.rational, work);
    if (value.radical.sign() == 0)
        return base;
    const auto root = phase::square_root(phase_rational(g, work), work);
    return phase::add(base, phase::multiply(phase_rational(value.radical, work), root, work), work);
}

phase::Interval interval_power(phase::Interval base, std::int64_t exponent, phase::Work& work) {
    // The schema admits int32, including its minimum. Widen before taking abs.
    const auto power = static_cast<std::uint64_t>(exponent < 0 ? -exponent : exponent);
    if ((power % 2) == 0)
        base = phase::absolute(base);
    auto result = phase::integer(1);
    for (auto bits = power; bits != 0; bits >>= 1U) {
        if ((bits % 2) != 0)
            result = phase::multiply(result, base, work);
        if (bits > 1) {
            const auto absolute = phase::absolute(base);
            base = phase::multiply(absolute, absolute, work);
        }
    }
    return exponent < 0 ? phase::divide(phase::integer(1), result, work) : result;
}

phase::Interval exp_endpoint(const phase::Bound& x, phase::Work& work) {
    if (x.mantissa == 0)
        return phase::integer(1);
    const auto magnitude = static_cast<std::int64_t>(phase::bits(x.mantissa)) + x.power;
    // The returned SI enclosure has a bounded exponent representation.
    if (magnitude > 14)
        partition_fail(FailureCode::PrecisionUnresolved);
    const auto halves = static_cast<int>(std::max<std::int64_t>(0, magnitude + 1));
    const auto reduced = phase::scale({x, x}, -halves);
    auto term = phase::integer(1);
    auto total = term;
    for (std::int64_t i = 1; i <= 64; ++i) {
        term = phase::divide(phase::multiply(term, reduced, work), phase::integer(i), work);
        total = phase::add(total, term, work);
    }
    // |reduced| <= 1/2. Twice the first omitted term bounds the entire tail.
    const auto next = phase::divide(phase::multiply(term, reduced, work), phase::integer(65), work);
    const phase::Interval radius{{}, phase::scale(phase::absolute(next), 1).upper};
    auto result = phase::Interval{phase::subtract(total, radius, work).lower,
                                  phase::add(total, radius, work).upper};
    for (int i = 0; i < halves; ++i)
        result = phase::multiply(result, result, work);
    return result;
}

struct ExpressionMath {
    explicit ExpressionMath(ReferenceArithmetic& reference)
        : math(reference),
          work{0, reference.limits.maximum_scalar_operations - reference.scalar_operations,
               reference.limits.cancellation} {}

    const phase::Interval& pi() {
        if (!pi_value)
            pi_value = phase::pi_constant(work);
        return *pi_value;
    }

    phase::Interval constant(const io::SiCase& source, const io::schema::ExprNodeConstant& node) {
        work.step();
        auto recipe = partition_take(source.quantity(node.value)).recipe;
        const auto power = recipe.pi_power;
        recipe.pi_power = 0;
        recipe.units.angle = io::SiAngleUnit::Radian;
        const auto exact =
            partition_take(io::materialize_rational(recipe, math.leaf_limits, math.plan));
        auto result = phase_rational(exact, work);
        if (power != 0 && exact.sign() != 0) {
            const auto radians_per_degree = phase::divide(pi(), phase::integer(180), work);
            result = phase::multiply(result, interval_power(radians_per_degree, power, work), work);
        }
        return result;
    }

    phase::Interval analytic(const io::schema::ExprNodeAnalytic& node, const phase::Interval& x) {
        using Op = decltype(node.op.selected());
        switch (node.op.selected()) {
        case Op::ItemSin:
        case Op::ItemCos:
            try {
                const auto rotation = phase::sine_cosine(x, phase::PhaseUnits::Radians, pi(), work);
                return node.op.selected() == Op::ItemSin ? rotation.sine : rotation.cosine;
            } catch (const phase::ArithmeticFailure& failure) {
                if (failure.code != FailureCode::PrecisionUnresolved)
                    throw;
                // Every finite phase is in the function's domain, including
                // boxes crossing several quadrants.
                return {phase::integer(-1).lower, phase::integer(1).upper};
            }
        case Op::ItemExp:
            return {exp_endpoint(x.lower, work).lower, exp_endpoint(x.upper, work).upper};
        case Op::ItemLog:
            if (!ln2)
                ln2 = phase::logarithm_series(
                    phase::divide(phase::integer(1), phase::integer(3), work), work);
            return {phase::log_endpoint(x.lower, *ln2, work).lower,
                    phase::log_endpoint(x.upper, *ln2, work).upper};
        case Op::ItemSqrt:
            return phase::square_root(x, work);
        }
        partition_fail(FailureCode::InvalidInput);
    }

    ReferenceArithmetic& math;
    phase::Work work;
    std::optional<phase::Interval> pi_value;
    std::optional<phase::Interval> ln2;
};

using CoordinatesBox = std::array<phase::Interval, 3>;

struct ExpressionCoordinates {
    CoordinatesBox domain, point;
};

ExpressionCoordinates reference_coordinates(const io::SiCase& source,
                                            const io::schema::ExpressionDomainReference& domain,
                                            const CasePgMaterialCoordinates& point,
                                            ReferenceArithmetic& math) {
    using Coordinates = decltype(domain.coordinates.selected());
    const auto kind = domain.coordinates.selected();
    if (kind == Coordinates::ItemMaterial && !point.has_material_metres)
        partition_fail(FailureCode::InvalidInput);
    const auto dimension = kind == Coordinates::ItemChart ? dimensionless : length;
    const auto lower = triple(source, domain.lower, dimension, math);
    const auto upper = triple(source, domain.upper, dimension, math);
    pg::detail::ReferenceRadicals exact{math, point.g};
    ExpressionCoordinates coordinates{};
    phase::Work work{0, math.limits.maximum_scalar_operations - math.scalar_operations,
                     math.limits.cancellation};
    for (std::size_t i = 0; i < 3; ++i) {
        if (kind == Coordinates::ItemChart) {
            if (math.compare(point.chart[i], lower[i]) < 0 ||
                math.compare(point.chart[i], upper[i]) > 0)
                partition_fail(FailureCode::InvalidInput);
            coordinates.point[i] = phase_rational(point.chart[i], work);
        } else {
            const auto& value = kind == Coordinates::ItemMaterial ? point.material_metres[i]
                                                                  : point.world_reference_metres[i];
            const pg::PgReferenceScalar a{math.copy(lower[i]), {}}, b{math.copy(upper[i]), {}};
            if (exact.sign(exact.sub(value, a)) < 0 || exact.sign(exact.sub(value, b)) > 0)
                partition_fail(FailureCode::InvalidInput);
            coordinates.point[i] = phase_radical(value, point.g, work);
        }
        coordinates.domain[i] = {phase_rational(lower[i], work).lower,
                                 phase_rational(upper[i], work).upper};
    }
    math.charge(work.used);
    return coordinates;
}

void evaluate_expression(const io::SiCase& source, const io::schema::Expression& expression,
                         const CoordinatesBox& input, std::span<phase::Interval> rows,
                         ExpressionMath& arithmetic) {
    auto& work = arithmetic.work;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        arithmetic.math.visit(1);
        work.step();
        const auto node = partition_take(expression.nodes.at(partition_index(i))).record();
        const auto operand = [&](const auto& field) -> const phase::Interval& {
            const auto ordinal = field.unsigned_value();
            if (ordinal >= i)
                partition_fail(FailureCode::InvalidInput);
            return rows[static_cast<std::size_t>(ordinal)];
        };
        rows[i] = std::visit(
            [&](const auto& typed) -> phase::Interval {
                using Node = std::decay_t<decltype(typed)>;
                if constexpr (std::is_same_v<Node, io::schema::ExprNodeConstant>) {
                    return arithmetic.constant(source, typed);
                } else if constexpr (std::is_same_v<Node, io::schema::ExprNodeVariable>) {
                    const auto axis = typed.axis.unsigned_value();
                    if (axis >= input.size())
                        partition_fail(FailureCode::InvalidInput);
                    return input[static_cast<std::size_t>(axis)];
                } else if constexpr (std::is_same_v<Node, io::schema::ExprNodeBinary>) {
                    const auto& a = operand(typed.first);
                    const auto& b = operand(typed.second);
                    using Op = decltype(typed.op.selected());
                    switch (typed.op.selected()) {
                    case Op::ItemAdd:
                        return phase::add(a, b, work);
                    case Op::ItemSubtract:
                        return phase::subtract(a, b, work);
                    case Op::ItemMultiply:
                        return phase::multiply(a, b, work);
                    case Op::ItemDivide:
                        return phase::divide(a, b, work);
                    }
                    partition_fail(FailureCode::InvalidInput);
                } else if constexpr (std::is_same_v<Node, io::schema::ExprNodeNegate>) {
                    return phase::negate(operand(typed.operand), work);
                } else if constexpr (std::is_same_v<Node, io::schema::ExprNodePower>) {
                    return interval_power(operand(typed.operand), typed.exponent.signed_value(),
                                          work);
                } else {
                    return arithmetic.analytic(typed, operand(typed.operand));
                }
            },
            node);
    }
}

void sample_expression(const io::SiCase& source, const io::schema::ReferenceFieldExpression& field,
                       const CasePgMaterialCoordinates& point,
                       std::span<pg::PgReferenceScalar> values, std::span<pg::ExactRational> errors,
                       io::SchemaDimension dimension, bool pi_multiple, ReferenceArithmetic& math) {
    const auto expr = field.expression.record();
    const auto domain = expr.domain.record();
    const auto* reference = std::get_if<io::schema::ExpressionDomainReference>(&domain);
    if (!reference || expr.outputs.size() != values.size() || expr.nodes.size() == 0)
        partition_fail(FailureCode::InvalidInput);
    const auto& typed = partition_take(source.source_case()).get();
    const auto input_quotas = partition_take(typed.root()).resources.record().input.record();
    const auto nodes = expr.nodes.size();
    // The typed schema has already certified edge count, depth, exponent
    // bounds and operand order. Bound both passes before allocating their row.
    if (nodes > input_quotas.maximum_expression_nodes.unsigned_value() ||
        input_quotas.maximum_expression_points.unsigned_value() == 0 ||
        nodes > (math.limits.maximum_partition_visits - math.partition_visits) / 2)
        partition_fail(FailureCode::ResourceExceeded);
    try {
        const auto coordinates = reference_coordinates(source, *reference, point, math);
        auto storage =
            partition_take(math.plan.allocate<phase::Interval>(partition_index(nodes), operation));
        ExpressionMath arithmetic{math};
        evaluate_expression(source, expr, coordinates.domain, storage.values(), arithmetic);
        evaluate_expression(source, expr, coordinates.point, storage.values(), arithmetic);
        const auto two = math.constant("2");
        for (std::size_t i = 0; i < values.size(); ++i) {
            const auto ordinal =
                partition_take(expr.outputs.at(partition_index(i))).unsigned_value();
            if (ordinal >= nodes)
                partition_fail(FailureCode::InvalidInput);
            auto output = storage.values()[static_cast<std::size_t>(ordinal)];
            if (pi_multiple)
                output = phase::divide(output, arithmetic.pi(), arithmetic.work);
            const auto lo = phase_endpoint(output.lower, math);
            const auto hi = phase_endpoint(output.upper, math);
            values[i].rational = math.div(math.add(lo, hi), two);
            auto error = quantity(source, partition_take(field.error.at(partition_index(i))),
                                  dimension, math, pi_multiple);
            if (error.sign() < 0)
                partition_fail(FailureCode::InvalidInput);
            errors[i] = math.add(error, math.div(math.sub(hi, lo), two));
        }
        math.charge(arithmetic.work.used);
    } catch (const phase::ArithmeticFailure& failure) {
        partition_fail(failure.code);
    } catch (const std::overflow_error&) {
        partition_fail(FailureCode::PrecisionUnresolved);
    } catch (const std::range_error&) {
        partition_fail(FailureCode::PrecisionUnresolved);
    }
}

void enclose_rotation(const pg::PgReferenceScalar& centre, const pg::ExactRational& error,
                      const pg::ExactRational& g, bool pi_multiple, CasePgProjectedAxisValue& out,
                      ReferenceArithmetic& math) {
    phase::Work work{0, math.limits.maximum_scalar_operations - math.scalar_operations,
                     math.limits.cancellation};
    try {
        auto phase_value = phase_radical(centre, g, work);
        if (error.sign() != 0) {
            const auto radius = phase_rational(error, work);
            phase_value = {phase::subtract(phase_value, radius, work).lower,
                           phase::add(phase_value, radius, work).upper};
        }
        const auto pi = phase::pi_constant(work);
        const auto units = pi_multiple ? phase::PhaseUnits::PiMultiple : phase::PhaseUnits::Radians;
        const auto rotation = phase::sine_cosine(phase_value, units, pi, work);
        math.charge(work.used);
        out.cosine = {phase_endpoint(rotation.cosine.lower, math),
                      phase_endpoint(rotation.cosine.upper, math)};
        out.sine = {phase_endpoint(rotation.sine.lower, math),
                    phase_endpoint(rotation.sine.upper, math)};
    } catch (const phase::ArithmeticFailure& failure) {
        partition_fail(failure.code);
    }
}

void validate_layer_point(const CasePgShellLayer& interval, const CasePgMaterialCoordinates& point,
                          ReferenceArithmetic& math) {
    if (point.g.sign() <= 0)
        partition_fail(FailureCode::GeometryInvalid);
    if (math.compare(point.chart[2], interval.lower) < 0 ||
        math.compare(point.chart[2], interval.upper) > 0)
        partition_fail(FailureCode::InvalidInput);
}
} // namespace

CasePgShellMaterials::CasePgShellMaterials(io::SchemaIdentity source, core::Index shell,
                                           core::PlannedBuffer<CasePgShellLayer>&& layers,
                                           Laws&& laws, std::size_t law_count,
                                           pg::PgReferenceCallInfo work) noexcept
    : source_(source),
      shell_(shell),
      layers_(std::move(layers)),
      laws_(std::move(laws)),
      law_count_(law_count),
      work_(work) {}

core::Result<CasePgShellMaterials>
CasePgShellMaterials::prepare(const io::SiCase& source, std::uint64_t shell_id,
                              const pg::PgReferenceLimits& limits,
                              core::ResourcePlan& memory) noexcept {
    return partition_entry([&]() -> core::Result<CasePgShellMaterials> {
        ReferenceArithmetic math{limits, memory};
        math.cancel();
        const auto& typed = partition_take(source.source_case()).get();
        const auto identity = partition_take(typed.identity());
        const auto model = partition_take(typed.root()).model.record();
        core::Index shell_ordinal = 0;
        for (; shell_ordinal < partition_index(model.shells.size()); ++shell_ordinal) {
            math.visit(1);
            if (partition_take(model.shells.at(shell_ordinal)).record().id.unsigned_value() ==
                shell_id)
                break;
        }
        const auto shell = partition_take(model.shells.at(shell_ordinal)).record();
        const auto count = shell.layers.size();
        if (count == 0)
            partition_fail(FailureCode::InvalidInput);
        const auto bounds_rationals = pg::detail::partition_mul(count, 2);
        const auto all_laws =
            pg::detail::partition_mul(law_rationals, std::min(count, model.materials.size()));
        math.admit(
            pg::detail::partition_add(pg::detail::reference_live_rationals,
                                      pg::detail::partition_add(bounds_rationals, all_laws)));
        auto layers =
            partition_take(memory.allocate<CasePgShellLayer>(partition_index(count), operation));
        auto laws = partition_take(memory.allocate<std::optional<pg::PgShellElasticLaw>>(
            partition_index(std::min(count, model.materials.size())), operation));
        struct Index {
            std::uint64_t id{};
            core::Index ordinal{};
        };
        auto by_material = partition_take(
            memory.allocate<Index>(partition_index(model.materials.size()), operation));
        auto by_layer = partition_take(memory.allocate<Index>(partition_index(count), operation));
        const auto less = [&](std::uint64_t a, std::uint64_t b) {
            math.visit(1);
            return a < b;
        };
        for (std::size_t i = 0; i < by_material.size(); ++i) {
            math.visit(1);
            const auto ordinal = partition_index(i);
            const auto material = partition_take(model.materials.at(ordinal)).record();
            by_material.values()[i] = {
                std::visit([](const auto& value) { return value.id.unsigned_value(); }, material),
                ordinal};
        }
        auto lower = math.constant("-1");
        for (std::size_t i = 0; i < layers.size(); ++i) {
            math.visit(1);
            const auto ordinal = partition_index(i);
            const auto layer = partition_take(shell.layers.at(ordinal)).record();
            auto lo = quantity(source, partition_take(layer.interval.at(0)), dimensionless, math);
            auto hi = quantity(source, partition_take(layer.interval.at(1)), dimensionless, math);
            if (!(lo == lower) || math.compare(hi, lo) <= 0)
                partition_fail(FailureCode::InvalidInput);
            lower = math.copy(hi);
            const auto material_id = layer.material.unsigned_value();
            layers.values()[i] = {core::LayerId{layer.id.unsigned_value()}, material_id,
                                  std::move(lo), std::move(hi), 0};
            by_layer.values()[i] = {material_id, ordinal};
        }
        if (!(lower == math.constant("1")))
            partition_fail(FailureCode::InvalidInput);

        const auto index = by_material.values();
        std::ranges::sort(index, less, &Index::id);
        std::ranges::sort(by_layer.values(), less, &Index::id);
        std::optional<std::uint64_t> previous;
        std::size_t prepared = 0;
        for (const auto& entry : by_layer.values()) {
            math.visit(1);
            if (entry.id != previous) {
                const auto found = std::ranges::lower_bound(index, entry.id, less, &Index::id);
                if (found == index.end() || found->id != entry.id || prepared >= laws.size())
                    partition_fail(FailureCode::InvalidInput);
                const auto record = partition_take(model.materials.at(found->ordinal)).record();
                const auto held = pg::detail::partition_add(
                    bounds_rationals, pg::detail::partition_mul(prepared, law_rationals));
                laws.values()[prepared++].emplace(std::visit(
                    [&](const auto& value) { return prepare_law(source, value, held, math); },
                    record));
                previous = entry.id;
            }
            layers.values()[static_cast<std::size_t>(entry.ordinal)].law = prepared - 1;
        }
        math.cancel();
        return CasePgShellMaterials{identity,        shell_ordinal, std::move(layers),
                                    std::move(laws), prepared,      math.info()};
    });
}

core::Result<std::reference_wrapper<const pg::PgShellElasticLaw>>
CasePgShellMaterials::law_for(core::Index layer) const& noexcept {
    if (layer < 0 || static_cast<std::uint64_t>(layer) >= layers_.size())
        return std::unexpected(core::Failure{FailureCode::InvalidInput, operation});
    const auto index = layers_.values()[static_cast<std::size_t>(layer)].law;
    return std::cref(laws_.values()[index].value());
}

core::Result<io::schema::Layer>
CasePgShellMaterials::source_layer(const io::SiCase& source, core::Index layer) const noexcept {
    return partition_entry([&]() -> core::Result<io::schema::Layer> {
        static_cast<void>(partition_take(law_for(layer)));
        const auto& typed = partition_take(source.source_case()).get();
        if (!io::detail::same_identity(partition_take(typed.identity()), source_))
            partition_fail(FailureCode::StaleGeneration);
        const auto model = partition_take(typed.root()).model.record();
        const auto shell = partition_take(model.shells.at(shell_)).record();
        return partition_take(shell.layers.at(layer)).record();
    });
}

core::Result<CasePgInitialGreenValue> CasePgShellMaterials::sample_initial_green(
    const io::SiCase& source, core::Index layer, const CasePgMaterialCoordinates& point,
    const pg::PgReferenceLimits& limits, core::ResourcePlan& memory) const noexcept {
    auto result = partition_entry([&]() -> core::Result<CasePgInitialGreenValue> {
        ReferenceArithmetic math{limits, memory};
        math.cancel();
        const auto selected = partition_take(source_layer(source, layer));
        const auto& interval = layers_.values()[static_cast<std::size_t>(layer)];
        math.admit();
        validate_layer_point(interval, point, math);
        CasePgInitialGreenValue out{interval.id, math.copy(point.g), {}, {}, {}};
        sample_field(source, selected.initial_green.record(), point, out.mandel, out.declared_error,
                     dimensionless, false, math);
        out.work = math.info();
        return out;
    });
    if (!result)
        result.error().operation = core::Operation::EvaluateStructure;
    return result;
}

pg::PgShellPointAxes CasePgAxisValue::axes() const& noexcept {
    if (const auto* projected = std::get_if<CasePgProjectedAxisValue>(&value))
        return pg::PgShellFieldFibre{projected->direction, projected->direction_error,
                                     projected->cosine, projected->sine};
    return pg::PgShellFieldAxes{std::get<CasePgExplicitAxisValue>(value).columns};
}

core::Result<CasePgAxisValue> CasePgShellMaterials::sample_axes(
    const io::SiCase& source, core::Index layer, const CasePgMaterialCoordinates& point,
    const pg::PgReferenceLimits& limits, core::ResourcePlan& memory) const noexcept {
    auto result = partition_entry([&]() -> core::Result<CasePgAxisValue> {
        ReferenceArithmetic math{limits, memory};
        math.cancel();
        const auto selected = partition_take(source_layer(source, layer));
        const auto& interval = layers_.values()[static_cast<std::size_t>(layer)];
        math.admit(384);
        validate_layer_point(interval, point, math);
        const auto axis = selected.axis.record();
        if (const auto* projected = std::get_if<io::schema::AxisProjected>(&axis)) {
            CasePgProjectedAxisValue value;
            sample_field(source, projected->direction.record(), point, value.direction,
                         value.direction_error, dimensionless, false, math);
            const auto& typed = partition_take(source.source_case()).get();
            const auto units = partition_take(typed.root()).units.record();
            const auto unit = units.angle.selected();
            const bool pi_multiple = unit == decltype(unit)::ItemDeg;
            std::array<pg::PgReferenceScalar, 1> radians;
            std::array<pg::ExactRational, 1> error;
            sample_field(source, projected->angle.record(), point, radians, error, angle,
                         pi_multiple, math);
            enclose_rotation(radians.front(), error.front(), point.g, pi_multiple, value, math);
            return CasePgAxisValue{interval.id, std::move(value), math.info()};
        }
        CasePgExplicitAxisValue value;
        std::array<pg::ExactRational, 9> error;
        const auto& explicit_axis = std::get<io::schema::AxisExplicit>(axis);
        sample_field(source, explicit_axis.columns.record(), point, value.columns, error,
                     dimensionless, false, math);
        // A box of independently uncertain explicit entries cannot certify
        // A^T A=I. A projected axis retains that identity by construction.
        if (std::ranges::any_of(error, [](const auto& x) { return x.sign() != 0; }))
            partition_fail(FailureCode::PrecisionUnresolved);
        return CasePgAxisValue{interval.id, std::move(value), math.info()};
    });
    if (!result)
        result.error().operation = core::Operation::EvaluateStructure;
    return result;
}
namespace {
constexpr std::uint64_t volume_sample_rationals = 48;
constexpr std::uint64_t recovered_side_rationals = 128;
constexpr std::uint64_t live_reference_rationals = 110;

pg::PgReferenceLimits remaining(ReferenceArithmetic& math, std::uint64_t kept = 0) {
    auto limits = math.limits;
    limits.maximum_scalar_operations -= math.scalar_operations;
    limits.maximum_partition_visits -= math.partition_visits;
    math.reserve_nested_rationals(limits.algebra, kept);
    return limits;
}

void record_work(ReferenceArithmetic& math, pg::PgReferenceCallInfo work) {
    math.charge(work.scalar_operations);
    math.visit(work.partition_visits);
}

std::uint64_t material_slots(const CasePgShellMaterials& materials) {
    return pg::detail::partition_add(
        pg::detail::partition_mul(materials.law_count(), law_rationals),
        pg::detail::partition_mul(materials.layers().size(), 2));
}

std::uint64_t expression_fields(const io::schema::Layer& layer) {
    const auto expression = [](const auto& field) -> std::uint64_t {
        return std::holds_alternative<io::schema::ReferenceFieldExpression>(field.record()) ? 1 : 0;
    };
    const auto axes = std::visit(
        [&](const auto& axis) -> std::uint64_t {
            using Axis = std::decay_t<decltype(axis)>;
            if constexpr (std::is_same_v<Axis, io::schema::AxisProjected>)
                return expression(axis.direction) + expression(axis.angle);
            else
                return expression(axis.columns);
        },
        layer.axis.record());
    return expression(layer.initial_green) + axes;
}

struct ShellFields {
    CasePgInitialGreenValue initial;
    CasePgAxisValue axes;
};

ShellFields sample_shell_fields(const CasePgShellMaterials& materials, const io::SiCase& source,
                                core::Index layer, const pg::PgReferenceValue& reference,
                                const pg::PgReferenceVector* material_metres,
                                std::uint64_t retained, ReferenceArithmetic& math) {
    const std::array chart{math.copy(reference.s), math.copy(reference.t),
                           math.copy(reference.zeta)};
    const pg::PgReferenceVector absent;
    const CasePgMaterialCoordinates coordinates{chart, material_metres ? *material_metres : absent,
                                                reference.position_metres, reference.g,
                                                material_metres != nullptr};
    const auto held = pg::detail::partition_add(retained, live_reference_rationals);
    auto initial = partition_take(materials.sample_initial_green(source, layer, coordinates,
                                                                 remaining(math, held), math.plan));
    record_work(math, initial.work);
    // The initial strain keeps its radicand and six value and error entries.
    auto axes = partition_take(
        materials.sample_axes(source, layer, coordinates,
                              remaining(math, pg::detail::partition_add(held, 22)), math.plan));
    record_work(math, axes.work);
    return {std::move(initial), std::move(axes)};
}

core::Index find_horizontal_face(std::span<const pg::PgPartitionHorizontalFace> faces,
                                 core::Index triangle, core::Index layer,
                                 pg::PgPartitionFaceSide outward, ReferenceArithmetic& math) {
    using Key = std::tuple<core::Index, core::Index, pg::PgPartitionFaceSide>;
    const auto key = [](const pg::PgPartitionHorizontalFace& face) {
        return Key{face.triangle, face.layer, face.side};
    };
    const Key wanted{triangle, layer, outward};
    // The partition publishes horizontal faces by triangle, layer, then side.
    const auto found = std::ranges::lower_bound(
        faces, wanted,
        [&](const Key& a, const Key& b) {
            math.visit(1);
            return a < b;
        },
        key);
    if (found == faces.end() || key(*found) != wanted)
        partition_fail(FailureCode::InvalidInput);
    return static_cast<core::Index>(found - faces.begin());
}
} // namespace

CasePgShellVolume::CasePgShellVolume(io::SchemaIdentity source, pg::PgPreparedPatchIdentity patch,
                                     std::uint64_t shell, CasePgShellMaterials&& materials,
                                     core::PlannedBuffer<Sample>&& samples) noexcept
    : source_(source),
      patch_(patch),
      shell_(shell),
      materials_(std::move(materials)),
      samples_(std::move(samples)) {}

core::Result<CasePgShellVolume>
CasePgShellVolume::prepare(const io::SiCase& source, const CaseStaticPgPointPatch& bound,
                           std::span<const pg::PgOwnedDiscreteVolumePoint> rule,
                           std::span<const pg::PgReferenceVector> material_metres,
                           const pg::PgReferenceLimits& limits,
                           core::ResourcePlan& memory) noexcept {
    return partition_entry([&]() -> core::Result<CasePgShellVolume> {
        ReferenceArithmetic math{limits, memory};
        math.cancel();
        const auto& typed = partition_take(source.source_case()).get();
        const auto identity = partition_take(typed.identity());
        if (!io::detail::same_identity(identity, bound.source()))
            partition_fail(FailureCode::StaleGeneration);
        if (!material_metres.empty() && material_metres.size() != rule.size())
            partition_fail(FailureCode::InvalidInput);
        const auto& patch = bound.patch();
        const auto info = partition_take(patch.info());
        const auto* basis = partition_take(patch.conditioning().reference_basis());
        const auto geometry = partition_take(basis->source());
        const auto volumes = partition_take(geometry.partition.volumes());
        const auto layers = partition_take(geometry.partition.layers());
        const auto count = partition_index(volumes.size());
        core::Index previous = -1;
        for (const auto& sample : rule) {
            math.visit(1);
            if (sample.volume < 0 || sample.volume >= count || sample.volume < previous ||
                sample.volume - previous > 1 || sample.weight.sign() <= 0)
                partition_fail(FailureCode::InvalidInput);
            previous = sample.volume;
        }
        if (count == 0 || previous != count - 1)
            partition_fail(FailureCode::InvalidInput);
        const auto shell = bound.shell().value();
        auto materials =
            partition_take(CasePgShellMaterials::prepare(source, shell, remaining(math), memory));
        record_work(math, materials.call_info());
        if (materials.layers().size() != layers.size())
            partition_fail(FailureCode::StaleGeneration);
        for (std::size_t j = 0; j < layers.size(); ++j) {
            math.visit(1);
            const auto& case_layer = materials.layers()[j];
            if (case_layer.id != layers[j].id || case_layer.lower != layers[j].lower ||
                case_layer.upper != layers[j].upper)
                partition_fail(FailureCode::StaleGeneration);
        }
        const auto quotas = partition_take(typed.root()).resources.record().input.record();
        const auto point_cap = quotas.maximum_expression_points.unsigned_value();
        core::Index previous_layer = -1;
        std::uint64_t expressions = 0, evaluated_points = 0;
        for (const auto& sample : rule) {
            math.visit(1);
            const auto layer = volumes[static_cast<std::size_t>(sample.volume)].layer;
            if (layer != previous_layer) {
                expressions =
                    expression_fields(partition_take(materials.source_layer(source, layer)));
                previous_layer = layer;
            }
            if (expressions > point_cap - evaluated_points)
                partition_fail(FailureCode::ResourceExceeded);
            evaluated_points += expressions;
        }
        const auto kept_rule =
            pg::detail::partition_add(pg::detail::partition_mul(4, rule.size()),
                                      pg::detail::partition_mul(6, material_metres.size()));
        const auto persistent = pg::detail::partition_add(material_slots(materials), kept_rule);
        math.admit(pg::detail::partition_add(
            pg::detail::partition_add(persistent, 384),
            pg::detail::partition_mul(rule.size(), volume_sample_rationals)));
        auto samples =
            partition_take(memory.allocate<Sample>(partition_index(rule.size()), operation));
        for (std::size_t i = 0; i < rule.size(); ++i) {
            math.visit(1);
            const auto& point = rule[i];
            const auto layer = volumes[static_cast<std::size_t>(point.volume)].layer;
            const auto retained = pg::detail::partition_add(
                persistent, pg::detail::partition_mul(i, volume_sample_rationals));
            auto reference = partition_take(basis->evaluate_geometry(
                pg::PgReferenceVolumePoint{point.volume, point.u, point.v, point.depth},
                remaining(math, retained), memory));
            record_work(math, reference.call_info);
            auto fields = sample_shell_fields(
                materials, source, layer, reference,
                material_metres.empty() ? nullptr : &material_metres[i], retained, math);
            auto& slot = samples.values()[i];
            slot.quadrature = {point.volume, math.copy(point.u), math.copy(point.v),
                               math.copy(point.depth), math.copy(point.weight)};
            slot.layer = layer;
            slot.initial.emplace(std::move(fields.initial));
            slot.axes.emplace(std::move(fields.axes));
        }
        return CasePgShellVolume{identity, info.identity, shell, std::move(materials),
                                 std::move(samples)};
    });
}

const pg::PgPreparedPatch&
CasePgShellVolume::checked_patch(const CaseStaticPgPointPatch& bound) const {
    const auto& patch = bound.patch();
    if (bound.shell().value() != shell_ || !io::detail::same_identity(bound.source(), source_) ||
        partition_take(patch.info()).identity != patch_)
        partition_fail(FailureCode::StaleGeneration);
    if (samples_.size() == 0)
        partition_fail(FailureCode::InvalidInput);
    return patch;
}

std::uint64_t CasePgShellVolume::retained_slots() const {
    return pg::detail::partition_add(
        material_slots(materials_),
        pg::detail::partition_mul(samples_.size(), volume_sample_rationals));
}

core::Result<pg::PgShellElasticVolumeAction>
CasePgShellVolume::elastic(const CaseStaticPgPointPatch& bound, const pg::PgShellPatchState& state,
                           const std::optional<pg::PgShellPatchDirection>& direction,
                           const pg::PgShellVolumeResolution& resolution,
                           const pg::PgReferenceLimits& limits,
                           core::ResourcePlan& memory) const noexcept {
    auto out = partition_entry([&]() -> core::Result<pg::PgShellElasticVolumeAction> {
        ReferenceArithmetic math{limits, memory};
        math.cancel();
        const auto& patch = checked_patch(bound);
        if (state.patch != patch_ || (direction && direction->patch != patch_))
            partition_fail(FailureCode::StaleGeneration);
        const auto retained = retained_slots();
        auto rows = partition_take(memory.allocate<std::optional<pg::PgShellVolumeSample>>(
            partition_index(samples_.size()), core::Operation::EvaluateStructure));
        for (std::size_t i = 0; i < samples_.size(); ++i) {
            math.visit(1);
            const auto& sample = samples_.values()[i];
            const auto& q = sample.quadrature;
            const auto& law = partition_take(materials_.law_for(sample.layer)).get();
            rows.values()[i].emplace(
                pg::PgShellVolumeSample{{{q.volume, q.u, q.v, q.depth}, q.weight},
                                        sample.initial->layer,
                                        law,
                                        sample.initial->enclosure(),
                                        sample.axes->axes()});
        }
        auto action = partition_take(pg::assemble_pg_shell_elastic_volume(
            patch, state, rows.values(), direction, resolution, remaining(math, retained), memory));
        record_work(math, action.work);
        action.work = math.info();
        return action;
    });
    if (!out)
        out.error().operation = core::Operation::EvaluateStructure;
    return out;
}
core::Result<pg::PgShellMassVolumeAction> CasePgShellVolume::mass(
    const CaseStaticPgPointPatch& bound, std::span<const pg::PgOwnedDiscreteVolumePoint> rule,
    const pg::PgShellPatchDirection& direction, const pg::PgShellMassResolution& resolution,
    const pg::PgReferenceLimits& limits, core::ResourcePlan& memory) const noexcept {
    auto out = partition_entry([&]() -> core::Result<pg::PgShellMassVolumeAction> {
        ReferenceArithmetic math{limits, memory};
        math.cancel();
        const auto& patch = checked_patch(bound);
        if (direction.patch != patch_)
            partition_fail(FailureCode::StaleGeneration);
        if (direction.coefficients_metres.size() != partition_take(patch.info()).basis.rank)
            partition_fail(FailureCode::InvalidInput);
        const auto layers = materials_.layers();
        const auto retained =
            pg::detail::partition_add(retained_slots(), pg::detail::partition_mul(rule.size(), 4));
        auto densities = partition_take(memory.allocate<std::optional<pg::PgShellLayerDensity>>(
            partition_index(layers.size()), core::Operation::EvaluateStructure));
        for (std::size_t i = 0; i < layers.size(); ++i) {
            math.visit(1);
            const auto& law = partition_take(materials_.law_for(partition_index(i))).get();
            densities.values()[i].emplace(
                pg::PgShellLayerDensity{layers[i].id, law.density_kg_m3()});
        }
        auto action = partition_take(
            pg::apply_pg_shell_consistent_mass(patch, rule, densities.values(), direction,
                                               resolution, remaining(math, retained), memory));
        record_work(math, action.work);
        action.work = math.info();
        return action;
    });
    if (!out)
        out.error().operation = core::Operation::EvaluateStructure;
    return out;
}
core::Result<CasePgInterfaceStress> CasePgShellVolume::recover_interface(
    const io::SiCase& source, const CaseStaticPgPointPatch& bound,
    const pg::PgShellPatchState& state, const CasePgInterfacePoint& query,
    const pg::PgReferenceLimits& limits, core::ResourcePlan& memory) const noexcept {
    auto result = partition_entry([&]() -> core::Result<CasePgInterfaceStress> {
        using Face = pg::PgPartitionFaceSide;
        ReferenceArithmetic math{limits, memory};
        math.cancel();
        const auto& patch = checked_patch(bound);
        if (state.patch != patch_)
            partition_fail(FailureCode::StaleGeneration);
        const auto& typed = partition_take(source.source_case()).get();
        if (!io::detail::same_identity(partition_take(typed.identity()), source_))
            partition_fail(FailureCode::StaleGeneration);
        const auto layers = materials_.layers();
        const auto lower = std::ranges::find_if(layers, [&](const CasePgShellLayer& layer) {
            math.visit(1);
            return layer.id == query.lower_layer;
        });
        if (lower == layers.end() || lower + 1 == layers.end() ||
            (lower + 1)->id != query.upper_layer)
            partition_fail(FailureCode::InvalidInput);
        const auto lower_layer = static_cast<core::Index>(lower - layers.begin());
        const auto upper_layer = lower_layer + 1;
        const auto* basis = partition_take(patch.conditioning().reference_basis());
        const auto reference = partition_take(basis->source());
        const auto faces = partition_take(reference.partition.horizontal_faces());
        const auto lower_face =
            find_horizontal_face(faces, query.triangle, lower_layer, Face::Upper, math);
        const auto upper_face =
            find_horizontal_face(faces, query.triangle, upper_layer, Face::Lower, math);
        const auto quotas = partition_take(typed.root()).resources.record().input.record();
        const auto used =
            expression_fields(partition_take(materials_.source_layer(source, lower_layer))) +
            expression_fields(partition_take(materials_.source_layer(source, upper_layer)));
        if (used > quotas.maximum_expression_points.unsigned_value())
            partition_fail(FailureCode::ResourceExceeded);
        const auto rank = partition_take(patch.info()).basis.rank;
        if (state.coordinates_metres.size() != rank)
            partition_fail(FailureCode::InvalidInput);
        const auto retained = retained_slots();
        const auto live_columns = pg::detail::partition_mul(14, rank);
        math.admit(
            pg::detail::partition_add(512, pg::detail::partition_add(retained, live_columns)));
        const auto read = [&](core::Index layer, core::Index face, Face outward,
                              const pg::PgReferenceVector* material_metres, std::uint64_t keep) {
            auto point = partition_take(pg::enclose_pg_shell_patch_point(
                patch, state, pg::PgReferenceHorizontalPoint{face, query.u, query.v},
                remaining(math, keep), memory));
            record_work(math, point.work);
            auto& ref = point.reference;
            auto& surface = std::get<pg::PgReferenceHorizontalMeasure>(ref.measure);
            const auto& model_layer = layers[static_cast<std::size_t>(layer)];
            if (ref.layer != model_layer.id || surface.side != outward)
                partition_fail(FailureCode::StaleGeneration);
            auto fields = sample_shell_fields(materials_, source, layer, ref, material_metres,
                                              pg::detail::partition_add(keep, live_columns), math);
            const auto& law = partition_take(materials_.law_for(layer)).get();
            const auto response_limits =
                remaining(math, pg::detail::partition_add(pg::detail::partition_add(keep, 48),
                                                          pg::detail::partition_mul(6, rank)));
            auto stress = partition_take(pg::enclose_pg_shell_material_point(
                law, fields.initial.enclosure(), fields.axes.axes(), point, response_limits,
                memory));
            record_work(math, stress.material.work);
            return CasePgLayerStress{ref.layer,
                                     face,
                                     outward,
                                     std::move(ref.zeta),
                                     std::move(ref.g),
                                     std::move(ref.position_metres),
                                     std::move(surface.area.outward),
                                     std::move(point.current_position_metres),
                                     std::move(stress)};
        };
        auto below =
            read(lower_layer, lower_face, Face::Upper, query.lower_material_metres, retained);
        auto above = read(upper_layer, upper_face, Face::Lower, query.upper_material_metres,
                          pg::detail::partition_add(retained, recovered_side_rationals));
        return CasePgInterfaceStress{patch_, std::move(below), std::move(above), math.info()};
    });
    if (!result)
        result.error().operation = core::Operation::EvaluateStructure;
    return result;
}
} // namespace uvlas::coupling::detail
