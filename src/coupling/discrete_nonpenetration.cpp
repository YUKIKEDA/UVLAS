#include "discrete_nonpenetration.hpp"

#include "aero/local_geometry_bounds_private.hpp"
#include "aero/p1_interval_private.hpp"
#include "aero/p1_shared_kernel.hpp"
#include "aero/query_coordinates_private.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <uvlas/numerics/dense_solve.hpp>

namespace uvlas::coupling::detail {
namespace a = aero::detail;
namespace p = a::p1;
namespace linear = numerics::linear;
namespace {
constexpr auto evaluate = core::Operation::EvaluateAero;
using Vector = std::array<p::Interval, 3>;

[[noreturn]] void reject(core::FailureCode code) {
    throw p::ArithmeticFailure{code};
}
template <class Result>
auto require(Result result) {
    if (!result) {
        const auto& error = result.error();
        if constexpr (std::is_same_v<typename Result::error_type, core::Failure>)
            throw error;
        else if constexpr (std::is_same_v<typename Result::error_type, a::GlobalGeometryFailure>)
            throw error.origin.cause;
        else
            throw error.cause;
    }
    return std::move(*result);
}
std::uint64_t remaining(const p::Work& work) noexcept {
    return work.limit - work.used;
}
std::uint64_t sum_budget(std::uint64_t a, std::uint64_t b) {
    if (b > std::numeric_limits<std::uint64_t>::max() - a)
        reject(core::FailureCode::ResourceExceeded);
    return a + b;
}
void account(p::Work& work, std::uint64_t amount) {
    if (amount > remaining(work))
        reject(core::FailureCode::ResourceExceeded);
    work.used += amount;
}
template <class T>
core::PlannedBuffer<T> allocate(core::ResourcePlan& resources, std::uint64_t count) {
    const auto size = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
    if (count > size)
        reject(core::FailureCode::ResourceExceeded);
    return require(resources.allocate<T>(static_cast<core::Index>(count), evaluate));
}

p::Interval checked_word(double word) {
    if (!std::isfinite(word))
        reject(core::FailureCode::Nonfinite);
    return a::decode_p1_word(word);
}
p::Interval import(const a::ScalarWords& words, p::Work& work) {
    work.step();
    if (const auto* single = std::get_if<a::SingleWord>(&words))
        return checked_word(single->value);
    const auto* two = std::get_if<a::DoubleWord>(&words);
    if (two == nullptr)
        reject(core::FailureCode::InvalidInput);
    return p::add(checked_word(two->high), checked_word(two->low), work);
}
Vector import(const NonpenetrationVector& value, p::Work& work) {
    return {import(value[0], work), import(value[1], work), import(value[2], work)};
}
p::Interval enclose(const p::Interval& center, double radius, int scale, p::Work& work) {
    const auto error = checked_word(radius);
    return p::scale(
        p::Interval{p::subtract(center, error, work).lower, p::add(center, error, work).upper},
        scale);
}
p::Interval import(const a::ScaledEnclosure& value, p::Work& work) {
    return enclose(import(value.center(), work), value.radius(), value.exponent(), work);
}
p::Interval import(const io::SiEnclosure& value, p::Work& work) {
    return std::visit(
        [&](const auto& word) {
            using S = std::decay_t<decltype(word)>;
            const auto center = [&]() -> a::ScalarWords {
                if constexpr (std::is_same_v<S, io::ScaledBinary64>)
                    return a::SingleWord{word.center};
                else
                    return a::DoubleWord{word.high, word.low};
            }();
            return enclose(import(center, work), word.radius, p::exponent(word.exponent2), work);
        },
        value);
}
Vector cross(const Vector& a, const Vector& b, p::Work& work) {
    Vector result{};
    for (std::size_t axis = 0; axis < result.size(); ++axis) {
        const auto next = (axis + 1) % 3;
        const auto last = (axis + 2) % 3;
        result[axis] = p::subtract(p::multiply(a[next], b[last], work),
                                   p::multiply(a[last], b[next], work), work);
    }
    return result;
}
p::Interval dot(const Vector& a, const Vector& b, p::Work& work) {
    p::Interval result{};
    for (std::size_t axis = 0; axis < a.size(); ++axis)
        result = p::add(result, p::multiply(a[axis], b[axis], work), work);
    return result;
}
Vector difference(const Vector& a, const Vector& b, p::Work& work) {
    return {p::subtract(a[0], b[0], work), p::subtract(a[1], b[1], work),
            p::subtract(a[2], b[2], work)};
}
void require_positive(const p::Interval& interval) {
    if (p::compare(interval.upper, {}) <= 0)
        reject(core::FailureCode::GeometryInvalid);
    if (p::compare(interval.lower, {}) <= 0)
        reject(core::FailureCode::PrecisionUnresolved);
}
Vector qualify_normal(const a::OwnedGeometry& geometry, const a::UvlmQuad& panel,
                      const DiscreteCollocation& input, p::Work& work) {
    const auto tx = import(input.tangents[0], work);
    const auto ty = import(input.tangents[1], work);
    auto normal = cross(tx, ty, work);
    const auto length2 = dot(normal, normal, work);
    require_positive(length2);
    for (auto facet : panel.facets) {
        const auto ids = require(geometry.facet_points(facet));
        const auto x0 =
            import(require(geometry.point_words(ids[0], a::Configuration::Current)), work);
        const auto x1 =
            import(require(geometry.point_words(ids[1], a::Configuration::Current)), work);
        const auto x2 =
            import(require(geometry.point_words(ids[2], a::Configuration::Current)), work);
        const auto edge1 = difference(x1, x0, work);
        const auto edge2 = difference(x2, x0, work);
        const auto outward = cross(edge1, edge2, work);
        require_positive(dot(normal, outward, work));
    }
    const auto length = p::square_root(length2, work);
    for (auto& coordinate : normal)
        coordinate = p::divide(coordinate, length, work);
    return normal;
}

// The geometry context presently holds a binary64 time. Compare the exact
// dyadic, using the canonical limbs of the independently qualified flow jet.
bool same_legacy_time(double word, a::pchip::RationalView time, p::Work& work) {
    work.step();
    const auto value = checked_word(word).lower;
    if (value.mantissa == 0)
        return time.sign == 0;
    const bool negative = value.mantissa < 0;
    if (time.sign != (negative ? -1 : 1))
        return false;
    const auto magnitude = p::magnitude(value.mantissa);
    const auto trailing = boost::multiprecision::lsb(magnitude);
    const p::Integer normalized = magnitude >> trailing;
    const auto mantissa = normalized.convert_to<std::uint64_t>();
    const auto power =
        p::exponent(static_cast<std::int64_t>(value.power) + static_cast<std::int64_t>(trailing));
    const auto equal_words = [&](std::span<const std::uint32_t> limbs, std::uint64_t digits,
                                 unsigned shift) {
        const auto width = static_cast<unsigned>(std::bit_width(digits)) + shift;
        if (limbs.size() != (width + 31U) / 32U)
            return false;
        const auto first = shift / 32U;
        const auto rem = shift % 32U;
        const std::array<std::uint32_t, 3> words{
            static_cast<std::uint32_t>(digits << rem),
            static_cast<std::uint32_t>((digits << rem) >> 32U),
            rem == 0 ? 0U : static_cast<std::uint32_t>(digits >> (64U - rem))};
        for (std::size_t i = 0; i < limbs.size(); ++i) {
            work.step();
            if (limbs[i] != (i < first ? 0U : words[i - first]))
                return false;
        }
        return true;
    };
    return equal_words(time.numerator, mantissa, static_cast<unsigned>(std::max(power, 0))) &&
           equal_words(time.denominator, 1, static_cast<unsigned>(std::max(-power, 0)));
}
bool same_time(const a::OwnedGeometry& geometry, a::pchip::RationalView supplied, p::Work& work) {
    const auto exact = require(geometry.exact_time());
    if (!exact)
        return same_legacy_time(geometry.context().time, supplied, work);
    work.step();
    const auto equal = [&](std::span<const std::uint32_t> l, std::span<const std::uint32_t> r) {
        return std::equal(l.begin(), l.end(), r.begin(), r.end(),
                          [&](std::uint32_t a, std::uint32_t b) {
                              work.step();
                              return a == b;
                          });
    };
    return exact->sign == supplied.sign && equal(exact->numerator, supplied.numerator) &&
           equal(exact->denominator, supplied.denominator);
}

struct PanelRow {
    a::UvlmQuad panel{
        a::PhysicalFaceId{0}, a::AeroFaceRole::Bound, {a::FacetId{0}, a::FacetId{0}}, {}, 0};
    std::size_t caller = 0;
};
struct Transfer {
    a::PhysicalFaceId face;
    bool reverse;
};
struct Trace {
    a::P1StrengthId id{0};
    a::ScalarWords value{a::SingleWord{0}};
    a::ScalarWords direction{a::SingleWord{0}};
    std::optional<Transfer> transfer{};
};
struct Assembly {
    core::PlannedBuffer<PanelRow> bound;
    core::PlannedBuffer<a::P1FaceStrengths> maps;
    core::PlannedBuffer<Trace> traces;
    core::PlannedBuffer<a::P1StrengthId> ids;
    core::PlannedBuffer<a::P1ComplexWords> values, directions;
};
a::ScalarWords negative(a::ScalarWords words) {
    if (const auto* one = std::get_if<a::SingleWord>(&words))
        return a::SingleWord{-one->value};
    const auto* two = std::get_if<a::DoubleWord>(&words);
    if (two == nullptr)
        reject(core::FailureCode::InvalidInput);
    return a::DoubleWord{-two->high, -two->low};
}
template <class T, class Key, class Projection>
const T& find(std::span<const T> items, Key key, Projection project, p::Work& work) {
    const auto it =
        std::lower_bound(items.begin(), items.end(), key, [&](const T& value, Key query) {
            work.step();
            return std::invoke(project, value) < query;
        });
    if (it == items.end() || std::invoke(project, *it) != key)
        reject(core::FailureCode::InvalidInput);
    return *it;
}

struct CirculationRowData {
    Vector normal{}, relative{};
    p::Interval rhs;
};

void match_circulation_geometry(const DiscreteCirculationRequest& request,
                                const a::GeometryContext& geo) {
    const auto& context = request.scales.context;
    const auto& key = context.cache;
    if ((geo.diagnostics.trial && *geo.diagnostics.trial != context.trial) ||
        key.fingerprint != geo.fingerprint || key.model != geo.model ||
        key.geometry != geo.current || key.precision != geo.precision || context.time != geo.time)
        reject(core::FailureCode::StaleGeneration);
}

void validate_circulation(const DiscreteCirculationRequest& request, std::size_t order,
                          p::Work& work, bool transposed = false) {
    if (request.scales.context.cache.representation != core::Precision::Binary64)
        reject(core::FailureCode::PrecisionUnsupported);
    if (request.solver_limits.maximum_order < 0 || request.solver_limits.maximum_rhs_columns < 0)
        reject(core::FailureCode::InvalidInput);
    if (std::cmp_greater(order, request.solver_limits.maximum_order) ||
        request.solver_limits.maximum_rhs_columns < 1)
        reject(core::FailureCode::ResourceExceeded);
    if (request.scales.unknown.size() != order || request.scales.residual.size() != order ||
        request.scales.unknown_units.size() != order ||
        request.scales.residual_units.size() != order)
        reject(core::FailureCode::InvalidInput);
    if (request.algebraic.norm != linear::Norm::Infinity &&
        request.algebraic.norm != linear::Norm::L2)
        reject(core::FailureCode::InvalidInput);
    for (const auto budget : {request.matrix_absolute, request.rhs_absolute,
                              request.original_absolute, request.original_relative}) {
        if (!std::isfinite(budget))
            reject(core::FailureCode::Nonfinite);
        if (budget < 0)
            reject(core::FailureCode::InvalidInput);
    }
    constexpr linear::UnitDimensions circulation_units{2, 0, -1, 0, 0, 0, 0};
    constexpr linear::UnitDimensions speed_units{1, 0, -1, 0, 0, 0, 0};
    constexpr linear::UnitDimensions dual_circulation{-2, 0, 1, 0, 0, 0, 0};
    constexpr linear::UnitDimensions dual_speed{-1, 0, 1, 0, 0, 0, 0};
    const auto& unknown_units = transposed ? dual_speed : circulation_units;
    const auto& residual_units = transposed ? dual_circulation : speed_units;
    for (std::size_t i = 0; i < order; ++i) {
        work.step();
        const auto& scales = request.scales;
        if (!std::isfinite(scales.unknown[i]) || !std::isfinite(scales.residual[i]))
            reject(core::FailureCode::Nonfinite);
        if (scales.unknown[i] <= 0 || scales.residual[i] <= 0 ||
            scales.unknown_units[i] != unknown_units || scales.residual_units[i] != residual_units)
            reject(core::FailureCode::InvalidInput);
    }
}

class NoDenseConsumerChecks final : public linear::ConsumerEvaluator<double> {
public:
    [[nodiscard]] numerics::SupportedPrecision supported_precision() const noexcept override {
        return numerics::SupportedPrecision::Binary64AndTwoWord;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const numerics::PrecisionBuffer<double>&, const numerics::PrecisionBuffer<double>&,
           const linear::Context&, const linear::ConsumerCheck& request,
           core::CertificateContext context) const override {
        return core::ErrorCertificate::unevaluated(context, request.budget);
    }
};

std::pair<numerics::DenseReal, numerics::DenseReal> dense_inputs(std::size_t order,
                                                                 core::ResourcePlan& memory) {
    const auto maximum = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
    if (order != 0 && order > maximum / order)
        reject(core::FailureCode::ResourceExceeded);
    const auto n = static_cast<core::Index>(order);
    const auto zeros = allocate<double>(memory, sum_budget(order * order, order));
    return {require(numerics::DenseReal::copy({n, n, n}, zeros.values().first(order * order),
                                              memory, core::Operation::SolveLinear)),
            require(numerics::DenseReal::copy({n, 1, n}, zeros.values().subspan(order * order),
                                              memory, core::Operation::SolveLinear))};
}

linear::Solution<double> solve_projected(const numerics::DenseReal& matrix, numerics::DenseReal rhs,
                                         const DiscreteCirculationRequest& request,
                                         core::ResourcePlan& memory, p::Work& work) {
    work.step();
    const auto& context = request.scales.context;
    const linear::DenseOperator op{matrix, context};
    const auto rhs_values = require(numerics::PrecisionBuffer<double>::binary64(
        std::move(rhs), core::NumericalGeneration{request.rhs_generation.value()}));
    const linear::Rhs<double> right{rhs_values, context, request.rhs_generation};
    auto factor =
        require(linear::DensePrepared::prepare(op, request.scales, request.solver_limits, memory));
    work.step();
    const NoDenseConsumerChecks no_checks;
    const linear::Precision<double> precision{
        linear::RetainedPrecision{request.retained_generation}, {}, no_checks};
    return require(linear::solve_linear(op, right, request.scales, precision, request.algebraic,
                                        request.solver_limits, factor, memory));
}

void match_linearization_context(linear::Context actual, const linear::Context& primal,
                                 bool transposed) {
    const bool is_distinct_matrix = actual.cache.matrix != primal.cache.matrix;
    if (transposed != is_distinct_matrix)
        reject(core::FailureCode::StaleGeneration);
    if (transposed) {
        if (actual.cache.scales == primal.cache.scales)
            reject(core::FailureCode::StaleGeneration);
        const auto is_primal_layout = [&](core::ScopeId layout) {
            return layout == primal.unknown_layout || layout == primal.residual_layout;
        };
        if (is_primal_layout(actual.unknown_layout) || is_primal_layout(actual.residual_layout) ||
            actual.unknown_layout == actual.residual_layout)
            reject(core::FailureCode::StaleGeneration);
        actual.cache.matrix = primal.cache.matrix;
        actual.unknown_layout = primal.unknown_layout;
        actual.residual_layout = primal.residual_layout;
    }
    actual.cache.scales = primal.cache.scales;
    if (actual != primal)
        reject(core::FailureCode::StaleGeneration);
}

void match_reused_scales(const linear::Scales& actual, const linear::Certificate& primal,
                         p::Work& work) {
    if (actual.context.cache.scales != primal.context().cache.scales)
        return;
    const auto equal = [&](auto supplied, auto certified) {
        if (supplied.size() != certified.size())
            return false;
        for (std::size_t i = 0; i < supplied.size(); ++i) {
            work.step();
            if (supplied[i] != certified[i])
                return false;
        }
        return true;
    };
    if (!equal(actual.unknown, primal.unknown_scales()) ||
        !equal(actual.residual, primal.residual_scales()) ||
        !equal(actual.unknown_units, primal.unknown_units()) ||
        !equal(actual.residual_units, primal.residual_units()))
        reject(core::FailureCode::StaleGeneration);
}

double project_binary64(const p::Interval& source, double budget, a::LocalRootScratch& scratch,
                        p::Work& work) {
    work.step();
    const auto scaled = a::enclose_p1_value(source, core::Precision::Binary64, scratch);
    const auto center = scaled.center();
    const auto* word = std::get_if<a::SingleWord>(&center);
    if (word == nullptr)
        reject(core::FailureCode::PrecisionUnsupported);
    const double value = std::scalbn(word->value, scaled.exponent());
    if (!std::isfinite(value))
        reject(core::FailureCode::PrecisionUnresolved);
    const auto error = p::absolute(p::subtract(source, checked_word(value), work));
    if (p::compare(error.upper, checked_word(budget).lower) > 0)
        reject(core::FailureCode::PrecisionUnresolved);
    return value;
}

Vector native_velocity(const Vector& relative, const a::P1ActionValue& induced, p::Work& work) {
    Vector total{};
    for (std::size_t axis = 0; axis < total.size(); ++axis)
        total[axis] = p::add(relative[axis], import(induced.velocity()[axis].real, work), work);
    return total;
}

p::Interval native_residual(const Vector& normal, const Vector& relative,
                            const a::P1ActionValue& induced, p::Work& work) {
    return dot(normal, native_velocity(relative, induced, work), work);
}

Vector normal_projection(const Vector& normal, const Vector& value, const p::Interval& length,
                         p::Work& work) {
    const auto projection = dot(normal, value, work);
    Vector output{};
    for (std::size_t axis = 0; axis < output.size(); ++axis)
        output[axis] =
            p::divide(p::subtract(value[axis], p::multiply(normal[axis], projection, work), work),
                      length, work);
    return output;
}

p::Interval tangent_length(const Vector& x, const Vector& y, p::Work& work) {
    const auto area = cross(x, y, work);
    return p::square_root(dot(area, area, work), work);
}

Vector scale_vector(const Vector& value, const p::Interval& factor, p::Work& work) {
    return {p::multiply(value[0], factor, work), p::multiply(value[1], factor, work),
            p::multiply(value[2], factor, work)};
}

CirculationVectorCovector enclose_vector(const Vector& value, core::Precision representation,
                                         a::LocalRootScratch& rounding, p::Work& work) {
    const auto enclose_axis = [&](std::size_t axis) {
        work.step();
        return a::enclose_p1_value(value[axis], representation, rounding);
    };
    return {enclose_axis(0), enclose_axis(1), enclose_axis(2)};
}

p::Interval fixed_circulation_direction(const DiscreteCollocation& row,
                                        const CirculationCollocationDirection& direction,
                                        const Vector& normal, const Vector& relative,
                                        const Vector& ambient_direction,
                                        const a::P1ActionValue& induced, p::Work& work) {
    const auto x = import(row.tangents[0], work);
    const auto y = import(row.tangents[1], work);
    const auto dx = import(direction.tangents[0], work);
    const auto dy = import(direction.tangents[1], work);
    const auto length = tangent_length(x, y, work);
    const auto d_first = cross(dx, y, work);
    const auto d_second = cross(x, dy, work);
    Vector d_area{}, d_total{};
    for (std::size_t axis = 0; axis < d_area.size(); ++axis)
        d_area[axis] = p::add(d_first[axis], d_second[axis], work);
    const auto d_normal = normal_projection(normal, d_area, length, work);
    const auto d_relative = difference(ambient_direction, import(direction.velocity, work), work);
    const auto& induction_direction = induced.jvp();
    if (!induction_direction)
        reject(core::FailureCode::InvalidInput);
    for (std::size_t axis = 0; axis < d_area.size(); ++axis)
        d_total[axis] =
            p::add(d_relative[axis], import((*induction_direction)[axis].real, work), work);
    return p::add(dot(d_normal, native_velocity(relative, induced, work), work),
                  dot(normal, d_total, work), work);
}

struct DyadicWords {
    std::array<std::uint32_t, 128> numerator{}, denominator{};
    std::size_t numerator_size = 0, denominator_size = 1;
    int sign = 0;

    DyadicWords(a::ExactInteger number, unsigned shift, p::Work& work) {
        work.step();
        denominator[0] = 1;
        if (number == 0)
            return;
        sign = number < 0 ? -1 : 1;
        if (sign < 0)
            number = -number;
        const auto trim =
            std::min(shift, static_cast<unsigned>(boost::multiprecision::lsb(number)));
        number >>= trim;
        shift -= trim;
        const auto word = shift / 32U;
        if (word >= denominator.size())
            reject(core::FailureCode::ResourceExceeded);
        denominator[0] = 0;
        denominator[word] = std::uint32_t{1} << (shift % 32U);
        denominator_size = static_cast<std::size_t>(word) + 1;
        while (number != 0) {
            work.step();
            if (numerator_size == numerator.size())
                reject(core::FailureCode::ResourceExceeded);
            numerator[numerator_size++] = (number & 0xffffffffU).convert_to<std::uint32_t>();
            number >>= 32U;
        }
    }

    [[nodiscard]] a::pchip::RationalView view() const noexcept {
        return {sign, std::span{numerator}.first(numerator_size),
                std::span{denominator}.first(denominator_size)};
    }
};

struct KnownRowData {
    std::optional<CaseKnownFlowJet> jet;
    Vector velocity{}, time{};
    std::array<Vector, 3> spatial{};
    std::array<std::uint64_t, 3> points{};
    Vector weights{};
    std::size_t point_count = 0;
};

struct KnownBoundary {
    const DiscreteKnownFlowInput& source;
    std::span<KnownRowData> rows;
    p::Interval time_direction;
};

KnownRowData evaluate_known_row(const a::OwnedGeometry& geometry, const a::ResolvedQuery& query,
                                const DiscreteKnownFlowInput& input, core::ResourcePlan& memory,
                                p::Work& work) {
    auto coordinate_limits = input.coordinate_limits;
    coordinate_limits.max_steps = std::min(coordinate_limits.max_steps, remaining(work));
    const auto position = require(a::exact_query_coordinates(geometry, query, coordinate_limits,
                                                             input.coordinates, work.stop));
    account(work, position.used_steps);
    const std::array xyz{
        DyadicWords{position.numerator[0], a::ExactQueryCoordinates::denominator_shift, work},
        DyadicWords{position.numerator[1], a::ExactQueryCoordinates::denominator_shift, work},
        DyadicWords{position.numerator[2], a::ExactQueryCoordinates::denominator_shift, work}};
    auto sampled = require(evaluate_case_known_flow(
        input.source, input.current, {input.time, {xyz[0].view(), xyz[1].view(), xyz[2].view()}},
        input.projection, input.workspace, input.flow_limits, memory, work.stop,
        geometry.context().diagnostics));
    KnownRowData result{};
    for (std::size_t i = 0; i < result.velocity.size(); ++i) {
        result.velocity[i] = import(require(sampled.velocity(i)), work);
        result.time[i] = import(require(sampled.time_partial(i)), work);
        for (std::size_t j = 0; j < result.velocity.size(); ++j)
            result.spatial[i][j] = import(require(sampled.spatial_partial(i, j)), work);
    }
    for (std::size_t i = 0; i < position.point_count; ++i) {
        work.step();
        result.points[i] = position.point_ids[i];
        result.weights[i] = p::exact(position.weights[i],
                                     -static_cast<int>(a::ExactQueryCoordinates::weight_shift));
    }
    result.point_count = position.point_count;
    result.jet.emplace(std::move(sampled));
    return result;
}

template <class Id>
std::size_t material_column(std::span<const Id> points, Id id, p::Work& work) {
    const auto found =
        std::lower_bound(points.begin(), points.end(), id, [&](auto left, auto right) {
            work.step();
            return left < right;
        });
    if (found == points.end() || *found != id)
        reject(core::FailureCode::InvalidInput);
    return static_cast<std::size_t>(found - points.begin());
}

Vector ambient_direction_at(const KnownRowData& known, const DiscreteCirculationDirection& dir,
                            std::size_t row, const Vector& uniform, const p::Interval& time,
                            p::Work& work) {
    Vector dx{};
    if (const auto& free = dir.rows[row].free_query) {
        dx = import(*free, work);
    } else {
        for (std::size_t i = 0; i < known.point_count; ++i) {
            const auto column = material_column(dir.points, a::PointId{known.points[i]}, work);
            const auto point = import(dir.displacements[column], work);
            for (std::size_t axis = 0; axis < dx.size(); ++axis)
                dx[axis] = p::add(dx[axis], p::multiply(known.weights[i], point[axis], work), work);
        }
    }
    Vector result{};
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = p::add(
            uniform[i],
            p::add(dot(known.spatial[i], dx, work), p::multiply(known.time[i], time, work), work),
            work);
    return result;
}

struct PrimalCirculationData {
    linear::Solution<double> solution;
    core::PlannedBuffer<std::optional<CirculationResidualRow>> rows;
    std::optional<a::P1ActionValue> geometry_owner;
    core::PlannedBuffer<p::Interval> influence;
    std::uint64_t steps;
};

struct PullbackData {
    PullbackData(const DiscreteCirculationRequest& settings,
                 std::span<const CirculationLinearRhs> objective) noexcept
        : request(settings), seed(objective) {}
    const DiscreteCirculationRequest& request;
    std::span<const CirculationLinearRhs> seed;
    std::optional<DiscreteCirculationLinearization> transpose{};
    core::PlannedBuffer<std::uint64_t> point_ids;
    core::PlannedBuffer<Vector> point_covectors;
    core::PlannedBuffer<std::optional<CirculationCollocationCovector>> rows;
    Vector flow{};
    p::Interval flow_time{};
};

struct LoadEdgeOccurrence {
    a::EdgeId edge{0};
    std::array<a::PointId, 2> endpoints{a::PointId{0}, a::PointId{0}};
    std::size_t panel = 0;
    int vortex_sign = 0;
    bool trailing = false;
};
struct LoadEdge {
    a::EdgeId edge{0};
    std::array<a::PointId, 2> endpoints{a::PointId{0}, a::PointId{0}};
    p::Interval strength{};
};
struct NodalLoad {
    std::uint64_t point = 0;
    Vector force{};
};

Vector point_position(const a::OwnedGeometry& geometry, a::PointId point, p::Work& work) {
    work.step();
    return import(require(geometry.point_words(point, a::Configuration::Current)), work);
}

void accumulate_vector(Vector& destination, const Vector& source, p::Work& work) {
    for (std::size_t axis = 0; axis < destination.size(); ++axis)
        destination[axis] = p::add(destination[axis], source[axis], work);
}

struct LoadTotals {
    Vector force{}, moment{};
    void add(const Vector& position, const Vector& value, const Vector& origin, p::Work& work) {
        accumulate_vector(force, value, work);
        accumulate_vector(moment, cross(difference(position, origin, work), value, work), work);
    }

    void add_direction(const Vector& position, const Vector& value, const Vector& dx,
                       const Vector& df, const Vector& origin, p::Work& work) {
        add(position, df, origin, work);
        accumulate_vector(moment, cross(dx, value, work), work);
    }
};

std::size_t load_node(std::span<const NodalLoad> nodes, a::PointId point, p::Work& work) {
    work.step();
    const auto it = std::ranges::lower_bound(
        nodes, point.value(),
        [&](std::uint64_t left, std::uint64_t right) {
            work.step();
            return left < right;
        },
        &NodalLoad::point);
    if (it == nodes.end() || it->point != point.value())
        reject(core::FailureCode::InvalidInput);
    return static_cast<std::size_t>(it - nodes.begin());
}

void distribute_load(std::span<NodalLoad> nodes, std::span<const a::PointId> vertices,
                     const Vector& force, const p::Interval& weight, p::Work& work) {
    const auto piece = scale_vector(force, weight, work);
    for (const auto point : vertices) {
        accumulate_vector(nodes[load_node(nodes, point, work)].force, piece, work);
    }
}

struct LoadTangentData {
    const DiscreteImplicitVectorLoadDirection& input;
    const DiscreteCirculationRequest& request;
    std::optional<DiscreteCirculationLinearization> result{};
};

struct LoadDirections {
    explicit LoadDirections(const DiscreteVectorLoadDirection& direction) noexcept
        : input(direction) {}
    const DiscreteVectorLoadDirection& input;
    core::PlannedBuffer<p::Interval> gamma, rate, strength;
    core::PlannedBuffer<NodalLoad> nodal;
    core::PlannedBuffer<std::optional<DiscreteLoadSampleDirection>> edges, facets;
    core::PlannedBuffer<std::optional<DiscreteNodalLoad>> points;
    std::optional<DiscreteLoadVector> force{}, moment{};
    Vector ambient{};
    LoadTotals totals{};
};

Vector add_vectors(Vector first, const Vector& second, p::Work& work) {
    accumulate_vector(first, second, work);
    return first;
}

Vector point_direction(const LoadDirections& direction, a::PointId id, p::Work& work) {
    work.step();
    const auto ids = direction.input.points;
    const auto it =
        std::lower_bound(ids.begin(), ids.end(), id, [&](a::PointId first, a::PointId second) {
            work.step();
            return first < second;
        });
    if (it == ids.end() || *it != id)
        reject(core::FailureCode::InvalidInput);
    return import(direction.input.displacements[static_cast<std::size_t>(it - ids.begin())], work);
}

Vector induction_direction(const a::P1ActionValue& induced, p::Work& work) {
    const auto& jvp = induced.jvp();
    if (!jvp)
        reject(core::FailureCode::InvalidInput);
    return {import((*jvp)[0].real, work), import((*jvp)[1].real, work),
            import((*jvp)[2].real, work)};
}

struct LoadPullbackData {
    explicit LoadPullbackData(const DiscreteVectorLoadWeights& weights) noexcept : input(weights) {}
    const DiscreteVectorLoadWeights& input;
    Vector force{}, moment{}, ambient{};
    core::PlannedBuffer<Vector> node_weights, position;
    core::PlannedBuffer<std::uint64_t> point_ids;
    core::PlannedBuffer<p::Interval> gamma, rate, strength;
    core::PlannedBuffer<std::optional<CirculationPointCovector>> points;
    core::PlannedBuffer<std::optional<LoadEdgeVelocityCovector>> edges;
    core::PlannedBuffer<std::optional<LoadCirculationCovector>> bound;
    std::optional<DiscreteLoadVector> flow{};
};

struct LoadSampleCovector {
    Vector position{}, force{};
};

void add_position_dual(LoadPullbackData& target, a::PointId id, const Vector& value,
                       p::Work& work) {
    const auto points = target.point_ids.values();
    work.step();
    const auto it = std::lower_bound(points.begin(), points.end(), id.value(),
                                     [&](std::uint64_t left, std::uint64_t right) {
                                         work.step();
                                         return left < right;
                                     });
    if (it == points.end() || *it != id.value())
        reject(core::FailureCode::InvalidInput);
    accumulate_vector(target.position.values()[static_cast<std::size_t>(it - points.begin())],
                      value, work);
}

LoadSampleCovector
load_sample_covector(const LoadPullbackData& target, const DiscreteLoadSampleWeights& seed,
                     std::span<const NodalLoad> nodes, std::span<const a::PointId> vertices,
                     const Vector& position, const Vector& value, const Vector& origin,
                     const p::Interval& weight, p::Work& work) {
    LoadSampleCovector result{
        add_vectors(import(seed.position, work), cross(value, target.moment, work), work),
        add_vectors(import(seed.force, work), target.force, work)};
    accumulate_vector(result.force, cross(target.moment, difference(position, origin, work), work),
                      work);
    for (const auto point : vertices) {
        const auto i = load_node(nodes, point, work);
        accumulate_vector(result.force, scale_vector(target.node_weights.values()[i], weight, work),
                          work);
    }
    return result;
}

a::ScalarWords retained_word(const numerics::PrecisionBuffer<double>& retained, std::size_t i) {
    if (const auto* one = std::get_if<numerics::Binary64<double>>(&retained.storage()))
        return a::SingleWord{one->values.values()[i]};
    const auto* two = std::get_if<numerics::TwoWord<double>>(&retained.storage());
    if (two == nullptr)
        reject(core::FailureCode::PrecisionUnsupported);
    return a::DoubleWord{two->hi.values()[i], two->lo.values()[i]};
}

class NativeResidualCheck {
public:
    explicit NativeResidualCheck(const DiscreteCirculationRequest& request) noexcept
        : request_(request) {}

    void add(const p::Interval& original, const p::Interval& rhs, std::size_t row, p::Work& work) {
        const auto tolerance = p::add(
            p::multiply(checked_word(request_.original_absolute),
                        checked_word(request_.scales.residual[row]), work),
            p::multiply(checked_word(request_.original_relative), p::absolute(rhs), work), work);
        const auto magnitude = p::absolute(original);
        if (p::compare(magnitude.upper, {}) == 0)
            return;
        if (p::compare(magnitude.lower, tolerance.upper) > 0)
            reject(core::FailureCode::LinearUnconverged);
        if (p::compare(tolerance.lower, {}) <= 0)
            reject(core::FailureCode::PrecisionUnresolved);
        const auto ratio = p::divide(magnitude, tolerance, work);
        if (request_.algebraic.norm == linear::Norm::L2)
            sum_squared_ = p::add(sum_squared_, p::multiply(ratio, ratio, work), work);
        else
            check_at_most_one(ratio);
    }

    void finish() const {
        if (request_.algebraic.norm == linear::Norm::L2)
            check_at_most_one(sum_squared_);
    }

private:
    static void check_at_most_one(const p::Interval& value) {
        const auto one = p::integer(1).lower;
        if (p::compare(value.lower, one) > 0)
            reject(core::FailureCode::LinearUnconverged);
        if (p::compare(value.upper, one) > 0)
            reject(core::FailureCode::PrecisionUnresolved);
    }
    const DiscreteCirculationRequest& request_;
    p::Interval sum_squared_{};
};

} // namespace

struct KnownCirculationRun {
    std::optional<PrimalCirculationData> primal;
    std::optional<DiscreteCirculationLinearization> linearization;
    core::PlannedBuffer<KnownRowData> known;
    core::PlannedBuffer<std::optional<CirculationForcingRow>> forcing;
    core::PlannedBuffer<std::optional<CirculationPointCovector>> points;
    core::PlannedBuffer<std::optional<CirculationCollocationCovector>> dual_rows;
    std::optional<CirculationVectorCovector> uniform;
    std::optional<a::ScaledEnclosure> time;
};

struct NonpenetrationStorage {
    std::optional<DiscreteNonpenetrationRow> value;
};
DiscreteNonpenetration::DiscreteNonpenetration(core::PlannedBuffer<NonpenetrationStorage> storage,
                                               ProjectedCaseUniformFlow&& flow,
                                               std::optional<a::P1ActionValue> geometry_owner,
                                               a::GeometryContext geometry_context,
                                               a::P1ActionStamp trial, std::uint64_t steps) noexcept
    : rows_(std::move(storage)),
      flow_(std::move(flow)),
      geometry_owner_(std::move(geometry_owner)),
      geometry_context_(geometry_context),
      trial_stamp_(trial),
      used_steps_(steps) {}
DiscreteNonpenetration::DiscreteNonpenetration(DiscreteNonpenetration&&) noexcept = default;
DiscreteNonpenetration&
DiscreteNonpenetration::operator=(DiscreteNonpenetration&&) noexcept = default;
DiscreteNonpenetration::~DiscreteNonpenetration() = default;
std::size_t DiscreteNonpenetration::size() const noexcept {
    return rows_.size();
}
core::Result<DiscreteNonpenetrationRow> DiscreteNonpenetration::row(std::size_t i) const noexcept {
    if (i >= rows_.size())
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    const auto& result = rows_.values()[i].value;
    if (!result)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    return *result;
}
const ProjectedCaseUniformFlow& DiscreteNonpenetration::flow() const& noexcept {
    return flow_;
}
std::uint64_t DiscreteNonpenetration::used_steps() const noexcept {
    return flow_.source().jet.active() && (rows_.size() == 0 || geometry_owner_) ? used_steps_ : 0;
}

a::GeometryContext DiscreteNonpenetration::geometry_context() const noexcept {
    return geometry_context_;
}
a::P1ActionStamp DiscreteNonpenetration::trial_stamp() const noexcept {
    return trial_stamp_;
}

struct NonpenetrationBuilder {
    static Assembly assemble(const DiscreteNonpenetrationInput& input,
                             const DiscreteNonpenetrationLimits& limits, core::ResourcePlan& memory,
                             p::Work& work, bool solve) {
        const auto& stop = work.stop;
        work.step();
        if (!a::p1_environment_supported())
            reject(core::FailureCode::PrecisionUnsupported);
        if (const auto invalid = limits.geometry_accuracy.validation_error())
            reject(*invalid);
        if (const auto invalid = limits.induction_accuracy.validation_error())
            reject(*invalid);
        const auto geometry_owner = require(input.contact.geometry());
        const auto& geometry = geometry_owner.get();
        if (!input.wake.empty() && geometry.context().diagnostics.parent != input.history.parent)
            reject(core::FailureCode::StaleGeneration);

        const auto counts = require(geometry.counts());
        const auto face_count = counts.faces;
        if (counts.quads != face_count || counts.triangles != 0)
            reject(core::FailureCode::InvalidInput);
        const auto bound_count = input.collocations.size();
        if (bound_count > face_count || input.wake.size() != face_count - bound_count ||
            input.contact.trailing_edge_count() > input.wake.size())
            reject(core::FailureCode::InvalidInput);
        const auto maximum = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
        if (face_count > maximum || input.history.columns.size() > maximum - face_count)
            reject(core::FailureCode::ResourceExceeded);
        auto bound = allocate<PanelRow>(memory, bound_count);
        auto maps = allocate<a::P1FaceStrengths>(memory, face_count);
        auto traces = allocate<Trace>(memory, face_count + input.history.columns.size());
        auto attachments = allocate<std::uint64_t>(memory, input.contact.trailing_edge_count());
        const auto trace = traces.values();
        std::size_t num_traces = 0;
        for (std::size_t i = 0; i < bound_count; ++i) {
            const auto& row = input.collocations[i];
            const auto panel = require(input.contact.uvlm_quad(row.facet, remaining(work), stop));
            account(work, panel.used_steps);
            if (panel.role != a::AeroFaceRole::Bound)
                reject(core::FailureCode::InvalidInput);
            bound.values()[i] = PanelRow{panel, i};
            maps.values()[i] = a::P1FaceStrengths{panel.face, a::P0BoundColumn{}};
            const auto zero = a::ScalarWords{a::SingleWord{0}};
            trace[num_traces++] = {
                a::P1StrengthId::bound(panel.face), solve ? zero : row.circulation,
                solve ? zero : row.circulation_direction, Transfer{panel.face, false}};
        }
        std::ranges::sort(bound.values(), [&](const PanelRow& l, const PanelRow& r) {
            work.step();
            return l.panel.face < r.panel.face;
        });
        auto prior = input.history.columns.end();
        for (auto it = input.history.columns.begin(); it != input.history.columns.end(); ++it) {
            work.step();
            if (it->id.bound_face() ||
                (prior != input.history.columns.end() && prior->id >= it->id))
                reject(core::FailureCode::InvalidInput);
            prior = it;
            trace[num_traces++] = {it->id, it->strength, a::SingleWord{0}};
        }
        const auto order = [&](const Trace& l, const Trace& r) {
            work.step();
            return l.id < r.id;
        };
        std::ranges::sort(trace.first(num_traces), order);
        for (std::size_t i = 1; i < num_traces; ++i) {
            work.step();
            if (trace[i - 1].id == trace[i].id)
                reject(core::FailureCode::InvalidInput);
        }
        const auto existing = num_traces;
        std::size_t visit = 0, attached = 0;
        std::uint64_t spare = 0;
        const auto new_column = [&] {
            for (; visit < existing; ++visit) {
                work.step();
                const auto old = trace[visit].id.wake_value();
                // Independent wake IDs precede the bound identity domain.
                if (!old || *old > spare)
                    break;
                if (*old == spare)
                    ++spare;
            }
            return a::P1StrengthId{spare++};
        };
        const auto history = [&](a::P1StrengthId id) {
            if (id.bound_face())
                reject(core::FailureCode::InvalidInput);
            return find(input.history.columns, id, &FrozenWakeTrace::id, work).id;
        };
        for (std::size_t i = 0; i < input.wake.size(); ++i) {
            const auto& row = input.wake[i];
            const auto panel = require(input.contact.uvlm_quad(row.facet, remaining(work), stop));
            account(work, panel.used_steps);
            const auto* previous = std::get_if<a::P1StrengthId>(&row.leading);
            if (panel.role != a::AeroFaceRole::Wake ||
                panel.attachment.has_value() == (previous != nullptr))
                reject(core::FailureCode::InvalidInput);
            auto leading = a::P1StrengthId{0};
            if (previous != nullptr) {
                leading = history(*previous);
            } else {
                if (!panel.attachment || attached >= attachments.size())
                    reject(core::FailureCode::InvalidInput);
                const auto& entry = *panel.attachment;
                attachments.values()[attached++] = entry.edge.value();
                const auto& source = find(
                    std::span<const PanelRow>{bound.values()}, entry.bound_face,
                    [](const PanelRow& b) { return b.panel.face; }, work);
                if (entry.coefficient == -1) {
                    const auto& collocation = input.collocations[source.caller];
                    leading = new_column();
                    const auto transfer = Transfer{entry.bound_face, true};
                    if (solve) {
                        trace[num_traces++] = {leading, a::DoubleWord{0, 0}, a::DoubleWord{0, 0},
                                               transfer};
                    } else {
                        trace[num_traces++] = {leading, negative(collocation.circulation),
                                               negative(collocation.circulation_direction),
                                               transfer};
                    }
                } else {
                    leading = a::P1StrengthId::bound(entry.bound_face);
                }
            }
            maps.values()[bound_count + i] =
                a::P1FaceStrengths{panel.face, std::array{leading, history(row.trailing)}};
        }
        if (attached != attachments.size())
            reject(core::FailureCode::InvalidInput);
        std::ranges::sort(attachments.values(), [&](std::uint64_t l, std::uint64_t r) {
            work.step();
            return l < r;
        });
        for (std::size_t i = 1; i < attachments.size(); ++i) {
            work.step();
            if (attachments.values()[i - 1] == attachments.values()[i])
                reject(core::FailureCode::InvalidInput);
        }
        std::ranges::sort(maps.values(),
                          [&](const a::P1FaceStrengths& l, const a::P1FaceStrengths& r) {
                              work.step();
                              return l.face < r.face;
                          });
        for (std::size_t i = 1; i < maps.size(); ++i) {
            work.step();
            if (maps.values()[i - 1].face == maps.values()[i].face)
                reject(core::FailureCode::InvalidInput);
        }
        std::ranges::sort(trace.first(num_traces), order);
        auto ids = allocate<a::P1StrengthId>(memory, num_traces);
        auto values = allocate<a::P1ComplexWords>(memory, num_traces);
        auto directions = allocate<a::P1ComplexWords>(memory, num_traces);
        for (std::size_t i = 0; i < num_traces; ++i) {
            work.step();
            (void)import(trace[i].value, work);
            (void)import(trace[i].direction, work);
            ids.values()[i] = trace[i].id;
            values.values()[i].real = trace[i].value;
            directions.values()[i].real = trace[i].direction;
        }
        return {std::move(bound), std::move(maps),   std::move(traces),
                std::move(ids),   std::move(values), std::move(directions)};
    }

    static void admit(const Assembly& assembly, const DiscreteNonpenetrationInput& input,
                      const DiscreteNonpenetrationLimits& limits, const core::ResourcePlan& memory,
                      std::uint64_t preparations, std::uint64_t applications, p::Work& work) {
        std::uint64_t payload = 0, reserved_steps = 0;
        const auto multiply_budget = [](std::uint64_t value, std::uint64_t count) {
            if (count != 0 && value > std::numeric_limits<std::uint64_t>::max() / count)
                reject(core::FailureCode::ResourceExceeded);
            return value * count;
        };
        const auto application_reserve = multiply_budget(limits.action_steps, applications);
        for (const auto& item : assembly.bound.values()) {
            const auto& row = input.collocations[item.caller];
            const a::P1ActionSource source{input.contact, row.query, input.trial_stamp};
            auto preflight_limits = limits.induction;
            preflight_limits.steps = std::min(preflight_limits.steps, remaining(work));
            const auto plan = require(a::forecast_bound_wake_induction(
                source, assembly.maps.values(), assembly.ids.values(), preflight_limits,
                work.stop));
            account(work, plan.preflight_steps);
            payload = sum_budget(payload, multiply_budget(plan.payload_bytes, preparations));
            const auto prepare_reserve = multiply_budget(plan.reserved_steps, preparations);
            reserved_steps =
                sum_budget(reserved_steps, sum_budget(prepare_reserve, application_reserve));
        }
        if (payload > memory.remaining_bytes() || reserved_steps > remaining(work))
            reject(core::FailureCode::ResourceExceeded);
    }

    static a::SharedP1Action prepare_action(const Assembly& assembly,
                                            const a::P1ActionSource& source,
                                            const DiscreteNonpenetrationLimits& limits,
                                            a::SharedP1Workspace& shared,
                                            a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
                                            core::ResourcePlan& memory, p::Work& work) {
        auto induction_limits = limits.induction;
        induction_limits.steps = std::min(induction_limits.steps, remaining(work));
        auto action = require(a::prepare_bound_wake_induction(
            source, assembly.maps.values(), assembly.ids.values(), limits.geometry_accuracy,
            limits.induction_accuracy, induction_limits, shared, local, kernel, memory, work.stop));
        account(work, require(action.report()).used_steps);
        return action;
    }

    static void match_flow(const DiscreteNonpenetrationInput& input,
                           a::pchip::PchipInputStamp stamp, a::pchip::RequestedSide side,
                           a::pchip::RationalView time, p::Work& work) {
        const auto owner = require(input.contact.geometry());
        if (stamp != input.flow_stamp || side != input.flow_side ||
            !same_time(owner.get(), time, work))
            reject(core::FailureCode::StaleGeneration);
    }

    static void match_flow(const DiscreteNonpenetrationInput& input,
                           const ProjectedCaseUniformFlow& flow, p::Work& work) {
        const auto& jet = flow.source().jet;
        match_flow(input, require(jet.stamp()), flow.source().requested_side, require(jet.time()),
                   work);
    }

    static Vector ambient_velocity(const ProjectedCaseUniformFlow& flow, p::Work& work) {
        return {import(require(flow.velocity(0)), work), import(require(flow.velocity(1)), work),
                import(require(flow.velocity(2)), work)};
    }

    static DiscreteNonpenetration
    evaluate_rows(const DiscreteNonpenetrationInput& input, ProjectedCaseUniformFlow&& flow,
                  const DiscreteNonpenetrationLimits& limits, a::SharedP1Workspace& shared,
                  a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
                  core::ResourcePlan& memory, const std::stop_token& stop) {
        p::Work work{0, limits.steps, stop};
        const auto geometry_owner = require(input.contact.geometry());
        const auto& geometry = geometry_owner.get();
        match_flow(input, flow, work);
        const auto assembly = assemble(input, limits, memory, work, false);
        const auto bound = assembly.bound.values();
        auto rows = allocate<NonpenetrationStorage>(memory, bound.size());
        admit(assembly, input, limits, memory, 1, 2, work);
        const auto ambient = ambient_velocity(flow, work);
        std::optional<a::P1ActionValue> geometry_retainer;
        for (std::size_t i = 0; i < bound.size(); ++i) {
            const auto& item = bound[i];
            const auto& row = input.collocations[item.caller];
            const auto normal = qualify_normal(geometry, item.panel, row, work);
            const auto relative = difference(ambient, import(row.velocity, work), work);
            const a::P1ActionSource source{input.contact, row.query, input.trial_stamp};
            const auto action =
                prepare_action(assembly, source, limits, shared, local, kernel, memory, work);
            const auto report = require(action.report());
            auto full = require(
                a::apply_shared_p1(action, source, assembly.values.values(), std::nullopt,
                                   std::min(limits.action_steps, remaining(work)), shared, stop));
            account(work, full.used_steps());
            const auto linear = require(
                a::apply_shared_p1(action, source, assembly.directions.values(), std::nullopt,
                                   std::min(limits.action_steps, remaining(work)), shared, stop));
            account(work, linear.used_steps());
            Vector velocity{}, delta{};
            for (std::size_t axis = 0; axis < normal.size(); ++axis) {
                velocity[axis] =
                    p::add(import(full.velocity()[axis].real, work), relative[axis], work);
                delta[axis] = import(linear.velocity()[axis].real, work);
            }
            a::LocalRootScratch rounding;
            work.step();
            const auto value = a::enclose_p1_value(
                dot(normal, velocity, work), limits.induction_accuracy.representation, rounding);
            work.step();
            const auto direction = a::enclose_p1_value(
                dot(normal, delta, work), limits.induction_accuracy.representation, rounding);
            rows.values()[i].value.emplace(DiscreteNonpenetrationRow{
                item.panel.face, value, direction, report.far_facets != 0});
            if (i == 0)
                geometry_retainer.emplace(std::move(full));
        }
        work.step();
        return {std::move(rows),    std::move(flow),   std::move(geometry_retainer),
                geometry.context(), input.trial_stamp, work.used};
    }

    static std::size_t transfer_column(std::span<const PanelRow> bound, const Transfer& transfer,
                                       p::Work& work) {
        const auto& row =
            find(bound, transfer.face, [](const PanelRow& b) { return b.panel.face; }, work);
        return static_cast<std::size_t>(&row - bound.data());
    }

    static a::P1ActionAdjoint
    adjoint_axis(const a::SharedP1Action& action, const a::P1ActionSource& source,
                 std::span<const a::P1ComplexWords> strengths, std::size_t axis,
                 const DiscreteNonpenetrationLimits& limits, a::SharedP1Workspace& shared,
                 core::ResourcePlan& memory, p::Work& work) {
        std::array<a::P1ComplexWords, 3> dual{};
        dual[axis].real = a::SingleWord{1};
        auto result = require(a::adjoint_shared_p1(action, source, strengths, dual,
                                                   std::min(limits.action_steps, remaining(work)),
                                                   shared, memory, work.stop));
        account(work, result.used_steps());
        return result;
    }

    static void scatter_geometry(const a::P1ActionAdjoint& adjoint, const p::Interval& weight,
                                 std::span<Vector> points, std::optional<Vector>& query,
                                 p::Work& work) {
        constexpr std::size_t axes = std::tuple_size_v<Vector>;
        for (std::size_t i = 0; i < points.size(); ++i)
            for (std::size_t axis = 0; axis < axes; ++axis) {
                auto& entry = points[i][axis];
                const auto value = import(require(adjoint.point(i, axis)), work);
                entry = p::add(entry, p::multiply(weight, value, work), work);
            }
        for (std::size_t axis = 0; axis < axes; ++axis) {
            const auto coordinate = require(adjoint.free_query(axis));
            if (coordinate.has_value() != query.has_value())
                reject(core::FailureCode::InvalidInput);
            if (coordinate && query) {
                auto& entry = (*query)[axis];
                entry = p::add(entry, p::multiply(weight, import(*coordinate, work), work), work);
            }
        }
    }

    static CirculationCollocationCovector
    pullback_row(PullbackData& target, std::size_t i, const PanelRow& panel,
                 const DiscreteCollocation& row, const CirculationRowData& data,
                 const a::P1ActionValue& primal, const a::SharedP1Action& action,
                 const a::P1ActionSource& source, const Assembly& assembly,
                 const DiscreteNonpenetrationLimits& limits, a::SharedP1Workspace& shared,
                 core::ResourcePlan& memory, a::LocalRootScratch& rounding, p::Work& work,
                 const KnownRowData* known = nullptr) {
        const auto& transpose = target.transpose;
        if (!transpose)
            reject(core::FailureCode::InvalidInput);
        const auto weight =
            p::negate(import(retained_word(transpose->solution().values(), i), work), work);
        const auto velocity_dual = scale_vector(data.normal, weight, work);
        const auto x = import(row.tangents[0], work);
        const auto y = import(row.tangents[1], work);
        const auto direction = native_velocity(data.relative, primal, work);
        const auto area_dual = scale_vector(
            normal_projection(data.normal, direction, tangent_length(x, y, work), work), weight,
            work);
        std::optional<Vector> query;
        if (std::holds_alternative<a::FreePoint>(require(source.query.definition()).definition))
            query.emplace();
        if (i == 0) {
            target.point_ids = allocate<std::uint64_t>(memory, action.point_count());
            target.point_covectors = allocate<Vector>(memory, action.point_count());
            for (std::size_t point = 0; point < action.point_count(); ++point) {
                work.step();
                target.point_ids.values()[point] = require(action.point(point)).value();
            }
        }
        for (std::size_t axis = 0; axis < velocity_dual.size(); ++axis) {
            target.flow[axis] = p::add(target.flow[axis], velocity_dual[axis], work);
            const auto dual = adjoint_axis(action, source, assembly.values.values(), axis, limits,
                                           shared, memory, work);
            scatter_geometry(dual, velocity_dual[axis], target.point_covectors.values(), query,
                             work);
        }
        if (known != nullptr) {
            target.flow_time =
                p::add(target.flow_time, dot(velocity_dual, known->time, work), work);
            Vector spatial_dual{};
            for (std::size_t j = 0; j < spatial_dual.size(); ++j)
                for (std::size_t d = 0; d < velocity_dual.size(); ++d)
                    spatial_dual[j] =
                        p::add(spatial_dual[j],
                               p::multiply(known->spatial[d][j], velocity_dual[d], work), work);
            if (query) {
                for (std::size_t j = 0; j < spatial_dual.size(); ++j)
                    (*query)[j] = p::add((*query)[j], spatial_dual[j], work);
            } else {
                const std::span<const std::uint64_t> ids = target.point_ids.values();
                for (std::size_t j = 0; j < known->point_count; ++j) {
                    const auto column = material_column(ids, known->points[j], work);
                    auto& point = target.point_covectors.values()[column];
                    for (std::size_t axis = 0; axis < point.size(); ++axis)
                        point[axis] =
                            p::add(point[axis],
                                   p::multiply(known->weights[j], spatial_dual[axis], work), work);
                }
            }
        }
        const auto representation = limits.induction_accuracy.representation;
        std::optional<CirculationVectorCovector> enclosed_query;
        if (query)
            enclosed_query.emplace(enclose_vector(*query, representation, rounding, work));
        const auto x_dual = cross(y, area_dual, work);
        const auto y_dual = cross(area_dual, x, work);
        return {panel.panel.face,
                enclosed_query,
                {enclose_vector(x_dual, representation, rounding, work),
                 enclose_vector(y_dual, representation, rounding, work)},
                enclose_vector(scale_vector(velocity_dual, p::integer(-1), work), representation,
                               rounding, work)};
    }

    static void assemble_influence(const Assembly& assembly, const a::SharedP1Action& action,
                                   const a::P1ActionSource& source, const Vector& normal,
                                   a::SharedP1Workspace& shared, core::ResourcePlan& memory,
                                   const DiscreteNonpenetrationLimits& limits,
                                   std::span<p::Interval> influence, p::Work& work) {
        const auto transpose = [&](std::size_t axis) {
            return adjoint_axis(action, source, assembly.directions.values(), axis, limits, shared,
                                memory, work);
        };
        const std::array coefficients{transpose(0), transpose(1), transpose(2)};
        for (auto& entry : influence) {
            work.step();
            entry = {};
        }
        for (std::size_t i = 0; i < assembly.ids.size(); ++i) {
            work.step();
            const auto& transfer = assembly.traces.values()[i].transfer;
            if (!transfer)
                continue;
            const Vector spatial{import(require(coefficients[0].strength(i)).real, work),
                                 import(require(coefficients[1].strength(i)).real, work),
                                 import(require(coefficients[2].strength(i)).real, work)};
            const auto projected = dot(normal, spatial, work);
            auto& column = influence[transfer_column(assembly.bound.values(), *transfer, work)];
            column = transfer->reverse ? p::subtract(column, projected, work)
                                       : p::add(column, projected, work);
        }
    }

    static PrimalCirculationData solve_data(
        const DiscreteNonpenetrationInput& input, const ProjectedCaseUniformFlow* flow,
        const DiscreteNonpenetrationLimits& limits, const DiscreteCirculationRequest& request,
        a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
        core::ResourcePlan& memory, p::Work& work, const DiscreteCirculationDirection* direction,
        std::span<p::Interval> forcing, PullbackData* pullback = nullptr,
        std::optional<Assembly>* retained_assembly = nullptr, KnownBoundary* known = nullptr) {
        const auto& stop = work.stop;
        work.step();
        const auto geometry_owner = require(input.contact.geometry());
        const auto& geometry = geometry_owner.get();
        const auto order = input.collocations.size();
        match_circulation_geometry(request, geometry.context());
        validate_circulation(request, order, work);
        auto [matrix, rhs] = dense_inputs(order, memory);
        if (known != nullptr) {
            const auto& time = known->source.time;
            match_flow(input, time.stamp, time.side, time.time, work);
        } else {
            if (flow == nullptr)
                reject(core::FailureCode::InvalidInput);
            match_flow(input, *flow, work);
        }
        auto assembly = assemble(input, limits, memory, work, true);
        const auto bound = std::span<const PanelRow>{assembly.bound.values()};
        auto data = allocate<CirculationRowData>(memory, order);
        auto results = allocate<std::optional<CirculationResidualRow>>(memory, order);
        auto influence = allocate<p::Interval>(memory, order * order);
        Vector ambient_direction{};
        if (direction != nullptr) {
            const auto counts = require(geometry.counts());
            if (direction->rows.size() != order ||
                direction->displacements.size() != direction->points.size() ||
                direction->points.size() > counts.points ||
                (order == 0 && !direction->points.empty()))
                reject(core::FailureCode::InvalidInput);
            ambient_direction = import(direction->uniform_flow, work);
            for (std::size_t i = 0; i < order; ++i) {
                work.step();
                const auto& row = direction->rows[i];
                const auto query = require(input.collocations[bound[i].caller].query.definition());
                if (row.face != bound[i].panel.face ||
                    row.free_query.has_value() !=
                        std::holds_alternative<a::FreePoint>(query.definition))
                    reject(core::FailureCode::InvalidInput);
                if (row.free_query)
                    (void)import(*row.free_query, work);
                (void)import(row.tangents[0], work);
                (void)import(row.tangents[1], work);
                (void)import(row.velocity, work);
            }
        }
        core::PlannedBuffer<p::Interval> objective;
        if (pullback != nullptr) {
            objective = allocate<p::Interval>(memory, order);
            for (std::size_t i = 0; i < order; ++i) {
                work.step();
                const auto& seed = pullback->seed[i];
                if (seed.face != bound[i].panel.face)
                    reject(core::FailureCode::InvalidInput);
                objective.values()[i] = import(seed.value, work);
            }
        }
        admit(assembly, input, limits, memory, 2, pullback != nullptr ? 8U : 5U, work);

        const auto ambient = flow != nullptr ? ambient_velocity(*flow, work) : Vector{};
        a::LocalRootScratch rounding;
        const auto apply = [&](const a::SharedP1Action& action, const a::P1ActionSource& source,
                               const std::optional<a::P1ActionDirection>& variation) {
            auto value = require(
                a::apply_shared_p1(action, source, assembly.values.values(), variation,
                                   std::min(limits.action_steps, remaining(work)), shared, stop));
            account(work, value.used_steps());
            return value;
        };
        for (std::size_t i = 0; i < order; ++i) {
            const auto& item = bound[i];
            const auto& row = input.collocations[item.caller];
            const a::P1ActionSource source{input.contact, row.query, input.trial_stamp};
            auto& record = data.values()[i];
            record.normal = qualify_normal(geometry, item.panel, row, work);
            const auto& known_velocity = [&]() -> const Vector& {
                if (known == nullptr)
                    return ambient;
                auto& sample = known->rows[i];
                sample = evaluate_known_row(geometry, row.query, known->source, memory, work);
                return sample.velocity;
            }();
            record.relative = difference(known_velocity, import(row.velocity, work), work);
            const auto action =
                prepare_action(assembly, source, limits, shared, local, kernel, memory, work);
            if (direction != nullptr && i == 0) {
                if (action.point_count() != direction->points.size())
                    reject(core::FailureCode::InvalidInput);
                for (std::size_t j = 0; j < action.point_count(); ++j) {
                    work.step();
                    if (require(action.point(j)) != direction->points[j])
                        reject(core::FailureCode::InvalidInput);
                    (void)import(direction->displacements[j], work);
                }
            }
            const auto old_wake = apply(action, source, std::nullopt);
            record.rhs =
                p::negate(native_residual(record.normal, record.relative, old_wake, work), work);
            rhs.values()[i] = project_binary64(record.rhs, request.rhs_absolute, rounding, work);
            auto native_row = influence.values().subspan(i * order, order);
            assemble_influence(assembly, action, source, record.normal, shared, memory, limits,
                               native_row, work);
            for (std::size_t col = 0; col < order; ++col)
                matrix.values()[col * order + i] =
                    project_binary64(native_row[col], request.matrix_absolute, rounding, work);
        }

        auto solution = solve_projected(matrix, std::move(rhs), request, memory, work);
        if (pullback != nullptr) {
            const auto before = work.used;
            const auto face = [&](std::size_t i) { return bound[i].panel.face; };
            auto transpose = linearize_native(CirculationSense::Transpose, face, influence.values(),
                                              objective.values(), pullback->request, memory, work);
            transpose.steps_ = work.used - before;
            pullback->transpose.emplace(std::move(transpose));
        }

        for (std::size_t i = 0; i < assembly.ids.size(); ++i) {
            work.step();
            const auto& transfer = assembly.traces.values()[i].transfer;
            if (!transfer)
                continue;
            const auto coordinate =
                retained_word(solution.values(), transfer_column(bound, *transfer, work));
            assembly.values.values()[i].real =
                transfer->reverse ? negative(coordinate) : coordinate;
        }

        NativeResidualCheck check{request};
        std::optional<a::P1ActionValue> geometry_parent;
        for (std::size_t i = 0; i < order; ++i) {
            const auto& item = bound[i];
            const auto& row = input.collocations[item.caller];
            const auto& record = data.values()[i];
            const a::P1ActionSource source{input.contact, row.query, input.trial_stamp};
            const auto action =
                prepare_action(assembly, source, limits, shared, local, kernel, memory, work);
            std::optional<a::P1ActionDirection> variation;
            if (direction != nullptr)
                variation.emplace(
                    a::P1ActionDirection{direction->displacements, assembly.directions.values(),
                                         direction->rows[i].free_query, a::SingleWord{0}});
            auto full = apply(action, source, variation);
            const auto* sample = known != nullptr ? &known->rows[i] : nullptr;
            if (direction != nullptr) {
                const auto flow_direction =
                    known != nullptr
                        ? ambient_direction_at(known->rows[i], *direction, i, ambient_direction,
                                               known->time_direction, work)
                        : ambient_direction;
                forcing[i] =
                    fixed_circulation_direction(row, direction->rows[i], record.normal,
                                                record.relative, flow_direction, full, work);
            }
            const auto original = native_residual(record.normal, record.relative, full, work);
            check.add(original, record.rhs, i, work);
            if (pullback != nullptr)
                pullback->rows.values()[i].emplace(
                    pullback_row(*pullback, i, item, row, record, full, action, source, assembly,
                                 limits, shared, memory, rounding, work, sample));
            const auto enclosed =
                a::enclose_p1_value(original, limits.induction_accuracy.representation, rounding);
            results.values()[i].emplace(CirculationResidualRow{
                item.panel.face, enclosed, require(action.report()).far_facets != 0});
            if (i == 0)
                geometry_parent.emplace(std::move(full));
        }
        check.finish();
        work.step();
        if (retained_assembly != nullptr)
            retained_assembly->emplace(std::move(assembly));
        return {std::move(solution), std::move(results), std::move(geometry_parent),
                std::move(influence), work.used};
    }

    static DiscreteCirculation finish_primal(PrimalCirculationData&& data,
                                             ProjectedCaseUniformFlow&& flow) noexcept {
        return {std::move(data.solution),       std::move(flow),           std::move(data.rows),
                std::move(data.geometry_owner), std::move(data.influence), data.steps};
    }

    static DiscreteCirculation
    solve(const DiscreteNonpenetrationInput& input, ProjectedCaseUniformFlow&& flow,
          const DiscreteNonpenetrationLimits& limits, const DiscreteCirculationRequest& request,
          a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
          core::ResourcePlan& memory, const std::stop_token& stop) {
        p::Work work{0, limits.steps, stop};
        auto data = solve_data(input, &flow, limits, request, shared, local, kernel, memory, work,
                               nullptr, {});
        return finish_primal(std::move(data), std::move(flow));
    }

    static void validate_tangent(const DiscreteCirculationRequest& tangent,
                                 const DiscreteCirculationRequest& primal, std::size_t order,
                                 p::Work& work) {
        match_linearization_context(tangent.scales.context, primal.scales.context, false);
        if (tangent.rhs_generation == primal.rhs_generation ||
            tangent.retained_generation == primal.retained_generation)
            reject(core::FailureCode::StaleGeneration);
        validate_circulation(tangent, order, work);
    }

    static DiscreteCirculationLinearization
    solve_retained_tangent(const PrimalCirculationData& primal, std::span<p::Interval> forcing,
                           const DiscreteCirculationRequest& tangent_request,
                           core::ResourcePlan& memory, p::Work& work) {
        for (auto& row : forcing)
            row = p::negate(row, work);
        const auto before = work.used;
        const auto face = [&](std::size_t i) { return face_at(primal.rows.values(), i); };
        auto tangent = linearize_native(CirculationSense::Tangent, face, primal.influence.values(),
                                        forcing, tangent_request, memory, work);
        work.step();
        tangent.steps_ = work.used - before;
        return tangent;
    }

    static DiscreteCirculationDirectionalSolve
    solve_direction(const DiscreteNonpenetrationInput& input,
                    const DiscreteCirculationDirection& direction, ProjectedCaseUniformFlow&& flow,
                    const DiscreteNonpenetrationLimits& limits,
                    const DiscreteCirculationRequest& primal_request,
                    const DiscreteCirculationRequest& tangent_request, a::SharedP1Workspace& shared,
                    a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
                    core::ResourcePlan& memory, const std::stop_token& stop) {
        p::Work work{0, limits.steps, stop};
        work.step();
        const auto order = input.collocations.size();
        validate_tangent(tangent_request, primal_request, order, work);
        auto native_rhs = allocate<p::Interval>(memory, order);
        auto forcing = allocate<std::optional<CirculationForcingRow>>(memory, order);
        auto primal = solve_data(input, &flow, limits, primal_request, shared, local, kernel,
                                 memory, work, &direction, native_rhs.values());
        match_reused_scales(tangent_request.scales, primal.solution.certificate(), work);
        a::LocalRootScratch rounding;
        for (std::size_t i = 0; i < order; ++i) {
            const auto& row = primal.rows.values()[i];
            if (!row)
                reject(core::FailureCode::InvalidInput);
            forcing.values()[i].emplace(CirculationForcingRow{
                row->face,
                a::enclose_p1_value(native_rhs.values()[i],
                                    limits.induction_accuracy.representation, rounding)});
        }
        auto tangent =
            solve_retained_tangent(primal, native_rhs.values(), tangent_request, memory, work);
        return {finish_primal(std::move(primal), std::move(flow)), std::move(tangent),
                std::move(forcing), work.used};
    }

    static DiscreteCirculationPullback solve_pullback(
        const DiscreteNonpenetrationInput& input, std::span<const CirculationLinearRhs> seed,
        ProjectedCaseUniformFlow&& flow, const DiscreteNonpenetrationLimits& limits,
        const DiscreteCirculationRequest& primal_request,
        const DiscreteCirculationRequest& transpose_request, a::SharedP1Workspace& shared,
        a::LocalPacketWorkspace& local, a::P1Workspace& kernel, core::ResourcePlan& memory,
        const std::stop_token& stop) {
        p::Work work{0, limits.steps, stop};
        work.step();
        const auto order = input.collocations.size();
        if (seed.size() != order)
            reject(core::FailureCode::InvalidInput);
        match_linearization_context(transpose_request.scales.context, primal_request.scales.context,
                                    true);
        validate_circulation(transpose_request, order, work, true);
        PullbackData target{transpose_request, seed};
        target.rows = allocate<std::optional<CirculationCollocationCovector>>(memory, order);
        auto primal = solve_data(input, &flow, limits, primal_request, shared, local, kernel,
                                 memory, work, nullptr, {}, &target);
        auto& transpose = target.transpose;
        if (!transpose)
            reject(core::FailureCode::InvalidInput);
        auto points =
            allocate<std::optional<CirculationPointCovector>>(memory, target.point_ids.size());
        const auto representation = limits.induction_accuracy.representation;
        a::LocalRootScratch rounding;
        for (std::size_t i = 0; i < points.size(); ++i) {
            work.step();
            points.values()[i].emplace(
                CirculationPointCovector{a::PointId{target.point_ids.values()[i]},
                                         enclose_vector(target.point_covectors.values()[i],
                                                        representation, rounding, work)});
        }
        const auto uniform = enclose_vector(target.flow, representation, rounding, work);
        work.step();
        return {finish_primal(std::move(primal), std::move(flow)),
                std::move(*transpose),
                std::move(points),
                std::move(target.rows),
                uniform,
                work.used};
    }

    static DiscreteKnownCirculation
    solve_known(const DiscreteNonpenetrationInput& input, const DiscreteKnownFlowInput& flow,
                const KnownCirculationDerivative& derivative,
                const DiscreteNonpenetrationLimits& limits,
                const DiscreteCirculationRequest& request, a::SharedP1Workspace& shared,
                a::LocalPacketWorkspace& local, a::P1Workspace& kernel, core::ResourcePlan& memory,
                const std::stop_token& stop) {
        p::Work work{0, limits.steps, stop};
        work.step();
        const auto order = input.collocations.size();
        if (order == 0)
            reject(core::FailureCode::InvalidInput);
        const auto* tangent = std::get_if<KnownCirculationTangent>(&derivative);
        const auto* adjoint = std::get_if<KnownCirculationAdjoint>(&derivative);
        if (tangent != nullptr)
            validate_tangent(tangent->request, request, order, work);
        if (adjoint != nullptr) {
            if (adjoint->seed.size() != order)
                reject(core::FailureCode::InvalidInput);
            match_linearization_context(adjoint->request.scales.context, request.scales.context,
                                        true);
            validate_circulation(adjoint->request, order, work, true);
        }

        auto output = allocate<KnownCirculationRun>(memory, 1);
        auto& result = output.values().front();
        result.known = allocate<KnownRowData>(memory, order);
        KnownBoundary known{flow, result.known.values(),
                            tangent == nullptr ? p::Interval{} : import(tangent->flow_time, work)};
        core::PlannedBuffer<p::Interval> native_rhs;
        if (tangent != nullptr) {
            native_rhs = allocate<p::Interval>(memory, order);
            result.forcing = allocate<std::optional<CirculationForcingRow>>(memory, order);
        }
        std::optional<PullbackData> pullback;
        if (adjoint != nullptr) {
            pullback.emplace(adjoint->request, adjoint->seed);
            pullback->rows = allocate<std::optional<CirculationCollocationCovector>>(memory, order);
        }
        auto primal =
            solve_data(input, nullptr, limits, request, shared, local, kernel, memory, work,
                       tangent == nullptr ? nullptr : &tangent->state, native_rhs.values(),
                       pullback ? &*pullback : nullptr, nullptr, &known);
        a::LocalRootScratch rounding;
        const auto representation = limits.induction_accuracy.representation;
        if (tangent != nullptr) {
            match_reused_scales(tangent->request.scales, primal.solution.certificate(), work);
            for (std::size_t i = 0; i < order; ++i) {
                work.step();
                result.forcing.values()[i].emplace(CirculationForcingRow{
                    face_at(primal.rows.values(), i),
                    a::enclose_p1_value(native_rhs.values()[i], representation, rounding)});
            }
            result.linearization.emplace(solve_retained_tangent(primal, native_rhs.values(),
                                                                tangent->request, memory, work));
        }
        if (pullback) {
            auto& transpose = pullback->transpose;
            if (!transpose)
                reject(core::FailureCode::InvalidInput);
            result.linearization.emplace(std::move(*transpose));
            result.dual_rows = std::move(pullback->rows);
            result.points = allocate<std::optional<CirculationPointCovector>>(
                memory, pullback->point_ids.size());
            for (std::size_t i = 0; i < result.points.size(); ++i) {
                work.step();
                result.points.values()[i].emplace(
                    CirculationPointCovector{a::PointId{pullback->point_ids.values()[i]},
                                             enclose_vector(pullback->point_covectors.values()[i],
                                                            representation, rounding, work)});
            }
            result.uniform.emplace(enclose_vector(pullback->flow, representation, rounding, work));
            result.time.emplace(a::enclose_p1_value(pullback->flow_time, representation, rounding));
        }
        result.primal.emplace(std::move(primal));
        work.step();
        return DiscreteKnownCirculation{std::move(output)};
    }

    static std::pair<std::size_t, int>
    take_edge_panel(std::span<const LoadEdgeOccurrence> occurrences, std::size_t& cursor,
                    p::Work& work) {
        const auto key = occurrences[cursor].edge;
        const auto panel = occurrences[cursor].panel;
        int orientation = 0;
        while (cursor < occurrences.size() && occurrences[cursor].edge == key &&
               occurrences[cursor].panel == panel) {
            work.step();
            orientation += occurrences[cursor].vortex_sign;
            ++cursor;
        }
        if (orientation < -1 || orientation > 1)
            reject(core::FailureCode::GeometryInvalid);
        return {panel, orientation};
    }

    static void reduce_load_edges(std::span<LoadEdgeOccurrence> occurrences,
                                  std::span<const EdgeMidpointVelocity> velocities,
                                  std::span<LoadEdge> edges,
                                  const linear::Solution<double>& solution, p::Work& work,
                                  LoadDirections* direction) {
        std::ranges::sort(
            occurrences, [&](const LoadEdgeOccurrence& left, const LoadEdgeOccurrence& right) {
                work.step();
                return std::pair{left.edge, left.panel} < std::pair{right.edge, right.panel};
            });
        std::size_t emitted = 0;
        for (std::size_t i = 0; i < occurrences.size();) {
            work.step();
            const auto owner = occurrences[i];
            p::Interval strength{}, d_strength{};
            bool exterior = false;
            while (i < occurrences.size() && occurrences[i].edge == owner.edge) {
                const auto [panel, orientation] = take_edge_panel(occurrences, i, work);
                if (orientation == 0)
                    continue; // The two occurrences of the same quad diagonal cancel exactly.
                exterior = true;
                const auto gamma = import(retained_word(solution.values(), panel), work);
                strength = orientation == 1 ? p::add(strength, gamma, work)
                                            : p::subtract(strength, gamma, work);
                if (direction != nullptr) {
                    const auto& d_gamma = direction->gamma.values()[panel];
                    d_strength = orientation == 1 ? p::add(d_strength, d_gamma, work)
                                                  : p::subtract(d_strength, d_gamma, work);
                }
            }
            if (exterior && !owner.trailing) {
                if (emitted == velocities.size() || velocities[emitted].edge != owner.edge)
                    reject(core::FailureCode::InvalidInput);
                (void)import(velocities[emitted].velocity, work);
                edges[emitted] = {owner.edge, owner.endpoints, strength};
                if (direction != nullptr)
                    direction->strength.values()[emitted] = d_strength;
                ++emitted;
            }
        }
        if (emitted != velocities.size())
            reject(core::FailureCode::InvalidInput);
    }

    static std::span<NodalLoad> prepare_load_topology(std::span<const PanelRow> bound,
                                                      const a::OwnedGeometry& geometry,
                                                      const a::GlobalContactState& contact,
                                                      std::span<LoadEdgeOccurrence> occurrences,
                                                      std::span<NodalLoad> nodes, p::Work& work) {
        std::size_t k = 0;
        for (std::size_t panel = 0; panel < bound.size(); ++panel)
            for (const auto facet : bound[panel].panel.facets) {
                work.step();
                const auto points = require(geometry.facet_points(facet));
                const auto edges = require(geometry.facet_edges(facet));
                const auto incidence = require(geometry.incidence(facet));
                for (std::size_t side = 0; side < edges.size(); ++side) {
                    const auto edge =
                        require(contact.uvlm_load_edge(edges[side], remaining(work), work.stop));
                    account(work, edge.used_steps);
                    occurrences[k] = {edge.edge, edge.endpoints, panel,
                                      -static_cast<int>(incidence[side]), edge.trailing};
                    nodes[k].point = points[side].value();
                    ++k;
                }
            }
        std::ranges::sort(nodes, [&](const NodalLoad& left, const NodalLoad& right) {
            work.step();
            return left.point < right.point;
        });
        const auto duplicates =
            std::ranges::unique(nodes, [&](const NodalLoad& left, const NodalLoad& right) {
                work.step();
                return left.point == right.point;
            });
        return nodes.first(static_cast<std::size_t>(duplicates.begin() - nodes.begin()));
    }

    static void prepare_load_directions(LoadDirections& direction, Assembly& assembly,
                                        const DiscreteVectorLoadInput& loads,
                                        const a::OwnedGeometry& geometry,
                                        std::span<const NodalLoad> vertices,
                                        const LoadTangentData* tangent, core::ResourcePlan& memory,
                                        p::Work& work) {
        const auto bound = std::span<const PanelRow>{assembly.bound.values()};
        direction.gamma = allocate<p::Interval>(memory, bound.size());
        direction.rate = allocate<p::Interval>(memory, bound.size());
        direction.strength = allocate<p::Interval>(memory, loads.edges.size());
        direction.nodal = allocate<NodalLoad>(memory, vertices.size());
        direction.edges =
            allocate<std::optional<DiscreteLoadSampleDirection>>(memory, loads.edges.size());
        direction.facets =
            allocate<std::optional<DiscreteLoadSampleDirection>>(memory, 2 * bound.size());
        direction.points = allocate<std::optional<DiscreteNodalLoad>>(memory, vertices.size());
        direction.ambient = import(direction.input.uniform_flow, work);
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            work.step();
            direction.nodal.values()[i].point = vertices[i].point;
        }
        for (std::size_t i = 0; i < direction.input.points.size(); ++i) {
            work.step();
            const auto point = direction.input.points[i];
            if ((i != 0 && direction.input.points[i - 1] >= point) ||
                require(geometry.owner(point)) != point)
                reject(core::FailureCode::InvalidInput);
            (void)import(direction.input.displacements[i], work);
        }
        for (std::size_t i = 0; i < loads.edges.size(); ++i) {
            work.step();
            if (direction.input.edges[i].edge != loads.edges[i].edge)
                reject(core::FailureCode::InvalidInput);
            (void)import(direction.input.edges[i].velocity, work);
        }
        const auto gamma = [&](std::size_t i) {
            if (tangent == nullptr)
                return direction.input.bound[i].circulation;
            if (!tangent->result)
                reject(core::FailureCode::InvalidInput);
            return retained_word(tangent->result->solution().values(), i);
        };
        for (std::size_t i = 0; i < bound.size(); ++i) {
            work.step();
            const auto rate = [&] {
                if (tangent != nullptr)
                    return tangent->input.rates[i];
                const auto& row = direction.input.bound[i];
                return BoundCirculationRate{row.face, row.circulation_rate};
            }();
            if (rate.face != bound[i].panel.face)
                reject(core::FailureCode::InvalidInput);
            direction.gamma.values()[i] = import(gamma(i), work);
            direction.rate.values()[i] = import(rate.rate, work);
        }
        for (std::size_t i = 0; i < assembly.ids.size(); ++i) {
            work.step();
            const auto& transfer = assembly.traces.values()[i].transfer;
            if (!transfer)
                continue;
            const auto delta = gamma(transfer_column(bound, *transfer, work));
            assembly.directions.values()[i].real = transfer->reverse ? negative(delta) : delta;
        }
    }

    static void prepare_load_pullback(LoadPullbackData& pullback, std::span<const NodalLoad> nodes,
                                      std::size_t edges, std::size_t order,
                                      core::ResourcePlan& memory, p::Work& work) {
        if (pullback.input.points.size() != nodes.size())
            reject(core::FailureCode::InvalidInput);
        pullback.force = import(pullback.input.force, work);
        pullback.moment = import(pullback.input.moment, work);
        pullback.node_weights = allocate<Vector>(memory, nodes.size());
        pullback.strength = allocate<p::Interval>(memory, edges);
        pullback.gamma = allocate<p::Interval>(memory, order);
        pullback.rate = allocate<p::Interval>(memory, order);
        pullback.edges = allocate<std::optional<LoadEdgeVelocityCovector>>(memory, edges);
        pullback.bound = allocate<std::optional<LoadCirculationCovector>>(memory, order);
        for (std::size_t i = 0; i < nodes.size(); ++i)
            pullback.node_weights.values()[i] = import(pullback.input.points[i], work);
        if (edges == 0) {
            pullback.position = allocate<Vector>(memory, nodes.size());
            pullback.point_ids = allocate<std::uint64_t>(memory, nodes.size());
            for (std::size_t i = 0; i < nodes.size(); ++i) {
                work.step();
                pullback.point_ids.values()[i] = nodes[i].point;
            }
        }
    }

    static void init_midpoint_pullback(LoadPullbackData& pullback, const a::SharedP1Action& action,
                                       core::ResourcePlan& memory, p::Work& work) {
        pullback.position = allocate<Vector>(memory, action.point_count());
        pullback.point_ids = allocate<std::uint64_t>(memory, action.point_count());
        for (std::size_t i = 0; i < action.point_count(); ++i) {
            work.step();
            pullback.point_ids.values()[i] = require(action.point(i)).value();
        }
    }

    static void scatter_midpoint_induction(LoadPullbackData& pullback,
                                           const a::SharedP1Action& action,
                                           const a::P1ActionSource& source,
                                           const Assembly& assembly, const Vector& velocity_dual,
                                           const DiscreteNonpenetrationLimits& limits,
                                           a::SharedP1Workspace& shared, core::ResourcePlan& memory,
                                           p::Work& work) {
        const auto bound = assembly.bound.values();
        std::optional<Vector> free_query;
        for (std::size_t axis = 0; axis < velocity_dual.size(); ++axis) {
            const auto dual = adjoint_axis(action, source, assembly.values.values(), axis, limits,
                                           shared, memory, work);
            scatter_geometry(dual, velocity_dual[axis], pullback.position.values(), free_query,
                             work);
            for (std::size_t i = 0; i < assembly.ids.size(); ++i) {
                work.step();
                const auto& transfer = assembly.traces.values()[i].transfer;
                if (!transfer)
                    continue;
                const auto coefficient = import(require(dual.strength(i)).real, work);
                const auto contribution = p::multiply(velocity_dual[axis], coefficient, work);
                auto& gamma = pullback.gamma.values()[transfer_column(bound, *transfer, work)];
                gamma = transfer->reverse ? p::subtract(gamma, contribution, work)
                                          : p::add(gamma, contribution, work);
            }
        }
    }

    static void scatter_force_edge_strength(LoadPullbackData& pullback,
                                            std::span<const LoadEdgeOccurrence> occurrences,
                                            std::span<const LoadEdge> edges, p::Work& work) {
        std::size_t at = 0;
        for (std::size_t i = 0; i < edges.size(); ++i) {
            const auto key = edges[i].edge;
            while (at < occurrences.size() && occurrences[at].edge < key) {
                work.step();
                ++at;
            }
            const auto& dual = pullback.strength.values()[i];
            while (at < occurrences.size() && occurrences[at].edge == key) {
                const auto [panel, orientation] = take_edge_panel(occurrences, at, work);
                auto& gamma = pullback.gamma.values()[panel];
                if (orientation == 1)
                    gamma = p::add(gamma, dual, work);
                else if (orientation == -1)
                    gamma = p::subtract(gamma, dual, work);
            }
        }
    }

    static void finish_load_pullback(LoadPullbackData& pullback, std::span<const PanelRow> bound,
                                     core::Precision representation, a::LocalRootScratch& scratch,
                                     core::ResourcePlan& memory, p::Work& work) {
        pullback.points =
            allocate<std::optional<CirculationPointCovector>>(memory, pullback.point_ids.size());
        for (std::size_t i = 0; i < pullback.points.size(); ++i) {
            work.step();
            pullback.points.values()[i].emplace(CirculationPointCovector{
                a::PointId{pullback.point_ids.values()[i]},
                enclose_vector(pullback.position.values()[i], representation, scratch, work)});
        }
        for (std::size_t i = 0; i < bound.size(); ++i) {
            work.step();
            pullback.bound.values()[i].emplace(LoadCirculationCovector{
                bound[i].panel.face,
                a::enclose_p1_value(pullback.gamma.values()[i], representation, scratch),
                a::enclose_p1_value(pullback.rate.values()[i], representation, scratch)});
        }
        pullback.flow.emplace(enclose_vector(pullback.ambient, representation, scratch, work));
    }

    static void propagate_load_circulation(
        LoadPullbackData& loads, PullbackData& implicit, const PrimalCirculationData& primal,
        const Assembly& assembly, const DiscreteNonpenetrationInput& input, const Vector& ambient,
        const DiscreteNonpenetrationLimits& limits, a::SharedP1Workspace& shared,
        a::LocalPacketWorkspace& local, a::P1Workspace& kernel, core::ResourcePlan& memory,
        a::LocalRootScratch& rounding, p::Work& work) {
        const auto before = work.used;
        const auto face = [&](std::size_t i) { return face_at(primal.rows.values(), i); };
        auto transpose =
            linearize_native(CirculationSense::Transpose, face, primal.influence.values(),
                             loads.gamma.values(), implicit.request, memory, work);
        transpose.steps_ = work.used - before;
        implicit.transpose.emplace(std::move(transpose));
        implicit.rows =
            allocate<std::optional<CirculationCollocationCovector>>(memory, assembly.bound.size());
        admit(assembly, input, limits, memory, 1, 4, work);
        const auto geometry = require(input.contact.geometry());

        for (std::size_t i = 0; i < assembly.bound.size(); ++i) {
            work.step();
            const auto& item = assembly.bound.values()[i];
            const auto& row = input.collocations[item.caller];
            const a::P1ActionSource source{input.contact, row.query, input.trial_stamp};
            CirculationRowData record;
            record.normal = qualify_normal(geometry.get(), item.panel, row, work);
            record.relative = difference(ambient, import(row.velocity, work), work);
            const auto action =
                prepare_action(assembly, source, limits, shared, local, kernel, memory, work);
            const auto full = require(a::apply_shared_p1(
                action, source, assembly.values.values(), std::nullopt,
                std::min(limits.action_steps, remaining(work)), shared, work.stop));
            account(work, full.used_steps());
            implicit.rows.values()[i].emplace(pullback_row(implicit, i, item, row, record, full,
                                                           action, source, assembly, limits, shared,
                                                           memory, rounding, work));
        }

        const auto direct_points = loads.point_ids.values();
        const auto implicit_points = implicit.point_ids.values();
        std::size_t next = 0;
        for (std::size_t i = 0; i < implicit_points.size(); ++i) {
            work.step();
            if (next < direct_points.size() && direct_points[next] == implicit_points[i]) {
                accumulate_vector(implicit.point_covectors.values()[i],
                                  loads.position.values()[next], work);
                ++next;
            }
        }
        if (next != direct_points.size())
            reject(core::FailureCode::InvalidInput);
        loads.point_ids = std::move(implicit.point_ids);
        loads.position = std::move(implicit.point_covectors);
        accumulate_vector(loads.ambient, implicit.flow, work);
    }

    static DiscreteCirculationLoads
    solve_loads(const DiscreteNonpenetrationInput& input, const DiscreteVectorLoadInput& loads,
                ProjectedCaseUniformFlow&& flow, const DiscreteNonpenetrationLimits& limits,
                const DiscreteCirculationRequest& request, a::GeometryWorkspace& geometry_work,
                a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local,
                a::P1Workspace& kernel, core::ResourcePlan& memory, const std::stop_token& stop,
                LoadDirections* direction = nullptr, LoadPullbackData* pullback = nullptr,
                PullbackData* implicit = nullptr, LoadTangentData* tangent = nullptr) {
        p::Work work{0, limits.steps, stop};
        work.step();
        const auto order = input.collocations.size();
        if (implicit != nullptr) {
            match_linearization_context(implicit->request.scales.context, request.scales.context,
                                        true);
            if (implicit->request.rhs_generation == request.rhs_generation ||
                implicit->request.retained_generation == request.retained_generation)
                reject(core::FailureCode::StaleGeneration);
            validate_circulation(implicit->request, order, work, true);
        }
        if (tangent != nullptr)
            validate_tangent(tangent->request, request, order, work);
        if (loads.stamp != input.trial_stamp || loads.circulation != request.retained_generation)
            reject(core::FailureCode::StaleGeneration);
        if (loads.rates.size() != order)
            reject(core::FailureCode::InvalidInput);
        const auto maximum = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
        if (order > maximum / 6)
            reject(core::FailureCode::ResourceExceeded);
        if (loads.edges.size() > 6 * order)
            reject(core::FailureCode::InvalidInput);
        if (pullback != nullptr) {
            const auto& seed = pullback->input;
            if (seed.stamp != loads.stamp || seed.circulation != loads.circulation)
                reject(core::FailureCode::StaleGeneration);
            if (seed.edges.size() != loads.edges.size() || seed.facets.size() != 2 * order)
                reject(core::FailureCode::InvalidInput);
        }
        const auto density = import(loads.density, work);
        if (p::compare(density.lower, {}) <= 0)
            reject(core::FailureCode::InvalidInput);
        const auto origin = import(loads.moment_origin, work);
        const auto snapshot = require(input.contact.geometry());
        const auto& geometry = snapshot.get();
        const auto radius = import(require(geometry.core_radius()), work);
        if (!loads.edges.empty() && p::compare(radius.lower, {}) <= 0)
            reject(core::FailureCode::GeometryInvalid);
        if (direction != nullptr) {
            const auto& variation = direction->input;
            if (variation.stamp != loads.stamp || variation.circulation != loads.circulation)
                reject(core::FailureCode::StaleGeneration);
            const auto counts = require(geometry.counts());
            const auto rate_count =
                tangent == nullptr ? variation.bound.size() : tangent->input.rates.size();
            if (variation.points.size() != variation.displacements.size() ||
                variation.points.size() > counts.points || rate_count != order ||
                variation.edges.size() != loads.edges.size() ||
                (order == 0 && !variation.points.empty()))
                reject(core::FailureCode::InvalidInput);
        }

        std::optional<Assembly> retained;
        core::PlannedBuffer<p::Interval> forcing;
        if (tangent != nullptr)
            forcing = allocate<p::Interval>(memory, order);
        const auto* state = tangent == nullptr ? nullptr : &tangent->input.state;
        auto primal = solve_data(input, &flow, limits, request, shared, local, kernel, memory, work,
                                 state, forcing.values(), nullptr, &retained);
        if (tangent != nullptr) {
            match_reused_scales(tangent->request.scales, primal.solution.certificate(), work);
            tangent->result.emplace(
                solve_retained_tangent(primal, forcing.values(), tangent->request, memory, work));
        }
        if (!retained)
            reject(core::FailureCode::InvalidInput);
        auto& assembly = *retained;
        const auto bound = std::span<const PanelRow>{assembly.bound.values()};
        auto rates = allocate<p::Interval>(memory, order);
        for (std::size_t i = 0; i < order; ++i) {
            work.step();
            if (loads.rates[i].face != bound[i].panel.face)
                reject(core::FailureCode::InvalidInput);
            rates.values()[i] = import(loads.rates[i].rate, work);
        }

        auto occurrences = allocate<LoadEdgeOccurrence>(memory, 6 * order);
        auto nodes = allocate<NodalLoad>(memory, 6 * order);
        auto edge_loads =
            allocate<std::optional<DiscreteCirculatoryLoad>>(memory, loads.edges.size());
        auto facet_loads = allocate<std::optional<DiscreteNoncirculatoryLoad>>(memory, 2 * order);
        const auto vertex_loads = prepare_load_topology(bound, geometry, input.contact,
                                                        occurrences.values(), nodes.values(), work);
        if (pullback != nullptr)
            prepare_load_pullback(*pullback, vertex_loads, loads.edges.size(), order, memory, work);
        auto point_loads = allocate<std::optional<DiscreteNodalLoad>>(memory, vertex_loads.size());
        auto edges = allocate<LoadEdge>(memory, loads.edges.size());
        if (direction != nullptr)
            prepare_load_directions(*direction, assembly, loads, geometry, vertex_loads, tangent,
                                    memory, work);
        reduce_load_edges(occurrences.values(), loads.edges, edges.values(), primal.solution, work,
                          direction);

        a::LocalRootScratch rounding;
        const auto representation = limits.induction_accuracy.representation;
        const auto finish = [&](const Vector& value) {
            return enclose_vector(value, representation, rounding, work);
        };
        const auto half = p::scale(p::integer(1), -1);
        const auto third = p::divide(p::integer(1), p::integer(3), work);
        LoadTotals totals;
        const auto ambient = ambient_velocity(flow, work);
        for (std::size_t i = 0; i < edges.size(); ++i) {
            work.step();
            const auto& edge = edges.values()[i];
            const auto a_pos = point_position(geometry, edge.endpoints[0], work);
            const auto b_pos = point_position(geometry, edge.endpoints[1], work);
            auto center = a_pos;
            accumulate_vector(center, b_pos, work);
            center = scale_vector(center, half, work);
            // The query API has no work readback. Charge its full explicit allowance.
            auto query_limits = limits.induction.geometry;
            query_limits.max_steps = std::min(query_limits.max_steps, remaining(work));
            const auto query = require(a::resolve_query(
                geometry,
                {a::QueryId{edge.edge.value()}, a::SourceEdgePoint{edge.edge, a::SingleWord{.5}}},
                geometry.context(), query_limits, geometry_work, memory, stop));
            account(work, query_limits.max_steps);
            const a::P1ActionSource source{input.contact, query, input.trial_stamp};
            auto action_limits = limits.induction;
            action_limits.steps = std::min(action_limits.steps, remaining(work));
            const auto action = require(a::prepare_bound_wake_midpoint_induction(
                source, edge.edge, assembly.maps.values(), assembly.ids.values(),
                limits.geometry_accuracy, limits.induction_accuracy, action_limits, shared, local,
                kernel, memory, stop));
            const auto report = require(action.report());
            account(work, report.used_steps);
            std::optional<a::P1ActionDirection> variation;
            if (direction != nullptr) {
                if (i == 0) {
                    if (action.point_count() != direction->input.points.size())
                        reject(core::FailureCode::InvalidInput);
                    for (std::size_t point = 0; point < action.point_count(); ++point) {
                        work.step();
                        if (require(action.point(point)) != direction->input.points[point])
                            reject(core::FailureCode::InvalidInput);
                    }
                }
                variation.emplace(a::P1ActionDirection{direction->input.displacements,
                                                       assembly.directions.values(), std::nullopt,
                                                       a::SingleWord{0}});
            }
            const auto induced = require(
                a::apply_shared_p1(action, source, assembly.values.values(), variation,
                                   std::min(limits.action_steps, remaining(work)), shared, stop));
            account(work, induced.used_steps());
            const auto relative = difference(ambient, import(loads.edges[i].velocity, work), work);
            const auto effective = native_velocity(relative, induced, work);
            const auto edge_vector = difference(b_pos, a_pos, work);
            const auto velocity_cross_edge = cross(effective, edge_vector, work);
            const auto factor = p::multiply(density, edge.strength, work);
            const auto force = scale_vector(velocity_cross_edge, factor, work);
            if (direction != nullptr) {
                const auto da = point_direction(*direction, edge.endpoints[0], work);
                const auto db = point_direction(*direction, edge.endpoints[1], work);
                const auto dx = scale_vector(add_vectors(da, db, work), half, work);
                const auto d_edge = difference(db, da, work);
                const auto d_relative = difference(
                    direction->ambient, import(direction->input.edges[i].velocity, work), work);
                const auto d_effective =
                    add_vectors(d_relative, induction_direction(induced, work), work);
                const auto d_cross = add_vectors(cross(d_effective, edge_vector, work),
                                                 cross(effective, d_edge, work), work);
                const auto d_factor = p::multiply(density, direction->strength.values()[i], work);
                const auto df = add_vectors(scale_vector(velocity_cross_edge, d_factor, work),
                                            scale_vector(d_cross, factor, work), work);
                direction->totals.add_direction(center, force, dx, df, origin, work);
                distribute_load(direction->nodal.values(), edge.endpoints, df, half, work);
                direction->edges.values()[i].emplace(
                    DiscreteLoadSampleDirection{finish(dx), finish(df)});
            }
            if (pullback != nullptr) {
                if (i == 0)
                    init_midpoint_pullback(*pullback, action, memory, work);
                const auto dual =
                    load_sample_covector(*pullback, pullback->input.edges[i], vertex_loads,
                                         edge.endpoints, center, force, origin, half, work);
                pullback->strength.values()[i] =
                    p::multiply(density, dot(dual.force, velocity_cross_edge, work), work);
                const auto velocity_dual =
                    scale_vector(cross(edge_vector, dual.force, work), factor, work);
                const auto edge_dual =
                    scale_vector(cross(dual.force, effective, work), factor, work);
                const auto midpoint_dual = scale_vector(dual.position, half, work);
                add_position_dual(*pullback, edge.endpoints[0],
                                  difference(midpoint_dual, edge_dual, work), work);
                add_position_dual(*pullback, edge.endpoints[1],
                                  add_vectors(midpoint_dual, edge_dual, work), work);
                accumulate_vector(pullback->ambient, velocity_dual, work);
                pullback->edges.values()[i].emplace(LoadEdgeVelocityCovector{
                    edge.edge, finish(scale_vector(velocity_dual, p::integer(-1), work))});
                scatter_midpoint_induction(*pullback, action, source, assembly, velocity_dual,
                                           limits, shared, memory, work);
            }
            totals.add(center, force, origin, work);
            distribute_load(vertex_loads, edge.endpoints, force, half, work);
            edge_loads.values()[i].emplace(
                DiscreteCirculatoryLoad{edge.edge, edge.endpoints, finish(center), finish(force),
                                        report.excluded_bound_occurrences, report.far_facets != 0});
        }
        for (std::size_t i = 0; i < order; ++i) {
            const auto factor =
                p::multiply(half, p::multiply(density, rates.values()[i], work), work);
            for (std::size_t child = 0; child < bound[i].panel.facets.size(); ++child) {
                work.step();
                const auto facet = bound[i].panel.facets[child];
                const auto points = require(geometry.facet_points(facet));
                const auto x = point_position(geometry, points[0], work);
                const auto y = point_position(geometry, points[1], work);
                const auto z = point_position(geometry, points[2], work);
                const auto xy = difference(y, x, work);
                const auto xz = difference(z, x, work);
                const auto doubled_area = cross(xy, xz, work);
                const auto force = scale_vector(doubled_area, factor, work);
                auto center = x;
                accumulate_vector(center, y, work);
                accumulate_vector(center, z, work);
                center = scale_vector(center, third, work);
                if (direction != nullptr) {
                    const auto dx = point_direction(*direction, points[0], work);
                    const auto dy = point_direction(*direction, points[1], work);
                    const auto dz = point_direction(*direction, points[2], work);
                    const auto dcenter =
                        scale_vector(add_vectors(add_vectors(dx, dy, work), dz, work), third, work);
                    const auto darea = add_vectors(cross(difference(dy, dx, work), xz, work),
                                                   cross(xy, difference(dz, dx, work), work), work);
                    const auto d_factor = p::multiply(
                        half, p::multiply(density, direction->rate.values()[i], work), work);
                    const auto df = add_vectors(scale_vector(doubled_area, d_factor, work),
                                                scale_vector(darea, factor, work), work);
                    direction->totals.add_direction(center, force, dcenter, df, origin, work);
                    distribute_load(direction->nodal.values(), points, df, third, work);
                    direction->facets.values()[2 * i + child].emplace(
                        DiscreteLoadSampleDirection{finish(dcenter), finish(df)});
                }
                if (pullback != nullptr) {
                    const auto dual = load_sample_covector(
                        *pullback, pullback->input.facets[2 * i + child], vertex_loads, points,
                        center, force, origin, third, work);
                    auto& rate = pullback->rate.values()[i];
                    rate = p::add(rate,
                                  p::multiply(p::multiply(density, half, work),
                                              dot(dual.force, doubled_area, work), work),
                                  work);
                    const auto y_dual = scale_vector(cross(xz, dual.force, work), factor, work);
                    const auto z_dual = scale_vector(cross(dual.force, xy, work), factor, work);
                    const auto centroid_dual = scale_vector(dual.position, third, work);
                    add_position_dual(
                        *pullback, points[0],
                        difference(centroid_dual, add_vectors(y_dual, z_dual, work), work), work);
                    add_position_dual(*pullback, points[1],
                                      add_vectors(centroid_dual, y_dual, work), work);
                    add_position_dual(*pullback, points[2],
                                      add_vectors(centroid_dual, z_dual, work), work);
                }
                totals.add(center, force, origin, work);
                distribute_load(vertex_loads, points, force, third, work);
                facet_loads.values()[2 * i + child].emplace(DiscreteNoncirculatoryLoad{
                    bound[i].panel.face, facet, points, finish(center), finish(force)});
            }
        }
        for (std::size_t i = 0; i < vertex_loads.size(); ++i) {
            work.step();
            const auto point = a::PointId{vertex_loads[i].point};
            point_loads.values()[i].emplace(
                DiscreteNodalLoad{point, finish(vertex_loads[i].force)});
            if (direction != nullptr)
                direction->points.values()[i].emplace(
                    DiscreteNodalLoad{point, finish(direction->nodal.values()[i].force)});
        }
        if (direction != nullptr) {
            direction->force.emplace(finish(direction->totals.force));
            direction->moment.emplace(finish(direction->totals.moment));
        }
        if (pullback != nullptr) {
            scatter_force_edge_strength(*pullback, occurrences.values(), edges.values(), work);
            if (implicit != nullptr)
                propagate_load_circulation(*pullback, *implicit, primal, assembly, input, ambient,
                                           limits, shared, local, kernel, memory, rounding, work);
            finish_load_pullback(*pullback, bound, representation, rounding, memory, work);
        }
        const auto force = finish(totals.force);
        const auto moment = finish(totals.moment);
        const auto reference = finish(origin);
        work.step();
        return {finish_primal(std::move(primal), std::move(flow)),
                std::move(edge_loads),
                std::move(facet_loads),
                std::move(point_loads),
                force,
                moment,
                reference,
                work.used};
    }

    static DiscreteCirculationLoadJvp solve_load_direction(
        const DiscreteNonpenetrationInput& input, const DiscreteVectorLoadInput& loads,
        const DiscreteVectorLoadDirection& variation, ProjectedCaseUniformFlow&& flow,
        const DiscreteNonpenetrationLimits& limits, const DiscreteCirculationRequest& request,
        a::GeometryWorkspace& geometry_work, a::SharedP1Workspace& shared,
        a::LocalPacketWorkspace& local, a::P1Workspace& kernel, core::ResourcePlan& memory,
        const std::stop_token& stop) {
        LoadDirections direction{variation};
        auto value = solve_loads(input, loads, std::move(flow), limits, request, geometry_work,
                                 shared, local, kernel, memory, stop, &direction);
        return finish_load_direction(std::move(value), direction);
    }

    static DiscreteCirculationLoadJvp finish_load_direction(DiscreteCirculationLoads&& value,
                                                            LoadDirections& direction) {
        if (!direction.force || !direction.moment)
            reject(core::FailureCode::InvalidInput);
        return {std::move(value),
                std::move(direction.edges),
                std::move(direction.facets),
                std::move(direction.points),
                *direction.force,
                *direction.moment};
    }

    static DiscreteCirculationImplicitLoadJvp solve_implicit_load_direction(
        const DiscreteNonpenetrationInput& input, const DiscreteVectorLoadInput& loads,
        const DiscreteImplicitVectorLoadDirection& variation, ProjectedCaseUniformFlow&& flow,
        const DiscreteNonpenetrationLimits& limits, const DiscreteCirculationRequest& request,
        const DiscreteCirculationRequest& tangent, a::GeometryWorkspace& geometry_work,
        a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
        core::ResourcePlan& memory, const std::stop_token& stop) {
        LoadTangentData stage{variation, tangent};
        const DiscreteVectorLoadDirection load_stage{variation.stamp,
                                                     variation.circulation,
                                                     variation.state.points,
                                                     variation.state.displacements,
                                                     variation.edges,
                                                     {},
                                                     variation.state.uniform_flow};
        LoadDirections direction{load_stage};
        auto value =
            solve_loads(input, loads, std::move(flow), limits, request, geometry_work, shared,
                        local, kernel, memory, stop, &direction, nullptr, nullptr, &stage);
        if (!stage.result)
            reject(core::FailureCode::InvalidInput);
        return {finish_load_direction(std::move(value), direction), std::move(*stage.result)};
    }

    static DiscreteCirculationLoadVjp solve_load_pullback(
        const DiscreteNonpenetrationInput& input, const DiscreteVectorLoadInput& loads,
        const DiscreteVectorLoadWeights& weights, ProjectedCaseUniformFlow&& flow,
        const DiscreteNonpenetrationLimits& limits, const DiscreteCirculationRequest& request,
        a::GeometryWorkspace& geometry_work, a::SharedP1Workspace& shared,
        a::LocalPacketWorkspace& local, a::P1Workspace& kernel, core::ResourcePlan& memory,
        const std::stop_token& stop) {
        LoadPullbackData target{weights};
        auto value = solve_loads(input, loads, std::move(flow), limits, request, geometry_work,
                                 shared, local, kernel, memory, stop, nullptr, &target);
        if (!target.flow)
            reject(core::FailureCode::InvalidInput);
        return {std::move(value), std::move(target.points), std::move(target.edges),
                std::move(target.bound), *target.flow};
    }

    static DiscreteCirculationImplicitLoadVjp solve_implicit_load_pullback(
        const DiscreteNonpenetrationInput& input, const DiscreteVectorLoadInput& loads,
        const DiscreteVectorLoadWeights& weights, ProjectedCaseUniformFlow&& flow,
        const DiscreteNonpenetrationLimits& limits, const DiscreteCirculationRequest& request,
        const DiscreteCirculationRequest& transpose, a::GeometryWorkspace& geometry_work,
        a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
        core::ResourcePlan& memory, const std::stop_token& stop) {
        LoadPullbackData target{weights};
        PullbackData implicit{transpose, {}};
        auto value = solve_loads(input, loads, std::move(flow), limits, request, geometry_work,
                                 shared, local, kernel, memory, stop, nullptr, &target, &implicit);
        if (!target.flow || !implicit.transpose)
            reject(core::FailureCode::InvalidInput);
        return {DiscreteCirculationLoadVjp{std::move(value), std::move(target.points),
                                           std::move(target.edges), std::move(target.bound),
                                           *target.flow},
                std::move(*implicit.transpose), std::move(implicit.rows)};
    }

    static a::PhysicalFaceId face_at(std::span<const std::optional<CirculationResidualRow>> rows,
                                     std::size_t i) {
        const auto& row = rows[i];
        if (!row)
            reject(core::FailureCode::InvalidInput);
        return row->face;
    }

    static DiscreteCirculationLinearization
    linearize(const DiscreteCirculation& primal, CirculationSense sense,
              std::span<const CirculationLinearRhs> input,
              const DiscreteCirculationRequest& request, std::uint64_t max_steps,
              core::ResourcePlan& memory, const std::stop_token& stop) {
        p::Work work{0, max_steps, stop};
        work.step();
        if (sense != CirculationSense::Tangent && sense != CirculationSense::Transpose)
            reject(core::FailureCode::InvalidInput);
        const bool transposed = sense == CirculationSense::Transpose;
        const auto order = primal.rows_.size();
        const auto shape = primal.solution_.values().shape();
        if (input.size() != order || shape.columns != 1 || !std::cmp_equal(shape.rows, order))
            reject(core::FailureCode::InvalidInput);
        const auto& previous = primal.solution_.certificate();
        match_linearization_context(request.scales.context, previous.context(), transposed);
        match_reused_scales(request.scales, previous, work);
        validate_circulation(request, order, work, transposed);
        auto right = allocate<p::Interval>(memory, order);
        for (std::size_t i = 0; i < order; ++i) {
            work.step();
            const auto source = require(primal.row(i));
            if (input[i].face != source.face)
                reject(core::FailureCode::InvalidInput);
            const auto value = import(input[i].value, work);
            right.values()[i] = transposed ? value : p::negate(value, work);
        }
        const auto face = [&](std::size_t i) { return face_at(primal.rows_.values(), i); };
        return linearize_native(sense, face, primal.influence_.values(), right.values(), request,
                                memory, work);
    }

    template <class FaceAt>
    static DiscreteCirculationLinearization
    linearize_native(CirculationSense sense, FaceAt face_at, std::span<const p::Interval> native,
                     std::span<const p::Interval> right, const DiscreteCirculationRequest& request,
                     core::ResourcePlan& memory, p::Work& work) {
        const auto order = right.size();
        const bool transposed = sense == CirculationSense::Transpose;
        const auto coefficient = [&](std::size_t row, std::size_t column) -> const p::Interval& {
            return transposed ? native[column * order + row] : native[row * order + column];
        };
        auto [matrix, rhs] = dense_inputs(order, memory);
        auto rows = allocate<std::optional<CirculationLinearRow>>(memory, order);
        auto retained = allocate<p::Interval>(memory, order);
        a::LocalRootScratch rounding;
        for (std::size_t i = 0; i < order; ++i) {
            rhs.values()[i] = project_binary64(right[i], request.rhs_absolute, rounding, work);
            for (std::size_t j = 0; j < order; ++j)
                matrix.values()[j * order + i] =
                    project_binary64(coefficient(i, j), request.matrix_absolute, rounding, work);
        }
        auto solution = solve_projected(matrix, std::move(rhs), request, memory, work);
        for (std::size_t i = 0; i < order; ++i)
            retained.values()[i] = import(retained_word(solution.values(), i), work);

        NativeResidualCheck check{request};
        for (std::size_t i = 0; i < order; ++i) {
            const auto& original_rhs = right[i];
            auto original = p::negate(original_rhs, work);
            for (std::size_t j = 0; j < order; ++j)
                original = p::add(original,
                                  p::multiply(coefficient(i, j), retained.values()[j], work), work);
            check.add(original, original_rhs, i, work);
            rows.values()[i].emplace(CirculationLinearRow{
                face_at(i), a::enclose_p1_value(original, core::Precision::TwoWord, rounding)});
        }
        check.finish();
        work.step();
        return {std::move(solution), std::move(rows), sense, work.used};
    }
};

core::Result<DiscreteNonpenetration> evaluate_discrete_nonpenetration(
    const DiscreteNonpenetrationInput& input, ProjectedCaseUniformFlow&& flow,
    const DiscreteNonpenetrationLimits& limits, a::SharedP1Workspace& shared,
    a::LocalPacketWorkspace& local, a::P1Workspace& kernel, core::ResourcePlan& memory,
    const std::stop_token& stop) noexcept {
    try {
        return NonpenetrationBuilder::evaluate_rows(input, std::move(flow), limits, shared, local,
                                                    kernel, memory, stop);
    } catch (const p::ArithmeticFailure& error) {
        return std::unexpected(core::Failure{error.code, evaluate});
    } catch (const core::Failure& error) {
        return std::unexpected(error);
    } catch (const std::bad_alloc&) {
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, evaluate});
    } catch (...) {
        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved, evaluate});
    }
}

DiscreteCirculation::DiscreteCirculation(
    linear::Solution<double> solution, ProjectedCaseUniformFlow&& flow,
    core::PlannedBuffer<std::optional<CirculationResidualRow>> rows,
    std::optional<a::P1ActionValue> geometry_owner, core::PlannedBuffer<p::Interval> influence,
    std::uint64_t steps) noexcept
    : solution_(std::move(solution)),
      flow_(std::move(flow)),
      rows_(std::move(rows)),
      geometry_owner_(std::move(geometry_owner)),
      influence_(std::move(influence)),
      steps_(steps) {}

DiscreteCirculation::DiscreteCirculation(DiscreteCirculation&&) noexcept = default;
DiscreteCirculation& DiscreteCirculation::operator=(DiscreteCirculation&&) noexcept = default;
DiscreteCirculation::~DiscreteCirculation() = default;

core::Result<CirculationResidualRow> DiscreteCirculation::row(std::size_t i) const noexcept {
    if (i >= rows_.size() || !geometry_owner_)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    const auto& residual = rows_.values()[i];
    if (!residual)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    return *residual;
}

core::Result<DiscreteCirculation> solve_discrete_circulation(
    const DiscreteNonpenetrationInput& input, ProjectedCaseUniformFlow&& flow,
    const DiscreteNonpenetrationLimits& limits, const DiscreteCirculationRequest& request,
    a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
    core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    try {
        return NonpenetrationBuilder::solve(input, std::move(flow), limits, request, shared, local,
                                            kernel, memory, stop);
    } catch (const p::ArithmeticFailure& error) {
        return std::unexpected(core::Failure{error.code, evaluate});
    } catch (const core::Failure& error) {
        return std::unexpected(error);
    } catch (const std::bad_alloc&) {
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, evaluate});
    } catch (...) {
        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved, evaluate});
    }
}

DiscreteCirculationLinearization::DiscreteCirculationLinearization(
    linear::Solution<double> solution,
    core::PlannedBuffer<std::optional<CirculationLinearRow>> rows, CirculationSense sense,
    std::uint64_t steps) noexcept
    : solution_(std::move(solution)), rows_(std::move(rows)), sense_(sense), steps_(steps) {}

core::Result<CirculationLinearRow>
DiscreteCirculationLinearization::row(std::size_t i) const noexcept {
    if (i >= rows_.size())
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    const auto& result = rows_.values()[i];
    if (!result)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    return *result;
}

core::Result<DiscreteCirculationLinearization> solve_discrete_circulation_linearization(
    const DiscreteCirculation& primal, CirculationSense sense,
    std::span<const CirculationLinearRhs> rhs, const DiscreteCirculationRequest& request,
    std::uint64_t max_steps, core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    try {
        return NonpenetrationBuilder::linearize(primal, sense, rhs, request, max_steps, memory,
                                                stop);
    } catch (const p::ArithmeticFailure& error) {
        return std::unexpected(core::Failure{error.code, evaluate});
    } catch (const core::Failure& error) {
        return std::unexpected(error);
    } catch (const std::bad_alloc&) {
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, evaluate});
    } catch (...) {
        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved, evaluate});
    }
}

DiscreteCirculationDirectionalSolve::DiscreteCirculationDirectionalSolve(
    DiscreteCirculation&& primal, DiscreteCirculationLinearization&& tangent,
    core::PlannedBuffer<std::optional<CirculationForcingRow>> forcing, std::uint64_t steps) noexcept
    : primal_(std::move(primal)),
      tangent_(std::move(tangent)),
      forcing_(std::move(forcing)),
      steps_(steps) {}

core::Result<CirculationForcingRow>
DiscreteCirculationDirectionalSolve::forcing(std::size_t i) const noexcept {
    if (i >= forcing_.size())
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    const auto& row = forcing_.values()[i];
    if (!row)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    return *row;
}

core::Result<DiscreteCirculationDirectionalSolve> solve_discrete_circulation_direction(
    const DiscreteNonpenetrationInput& input, const DiscreteCirculationDirection& direction,
    ProjectedCaseUniformFlow&& flow, const DiscreteNonpenetrationLimits& limits,
    const DiscreteCirculationRequest& primal, const DiscreteCirculationRequest& tangent,
    a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
    core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    try {
        return NonpenetrationBuilder::solve_direction(input, direction, std::move(flow), limits,
                                                      primal, tangent, shared, local, kernel,
                                                      memory, stop);
    } catch (const p::ArithmeticFailure& error) {
        return std::unexpected(core::Failure{error.code, evaluate});
    } catch (const core::Failure& error) {
        return std::unexpected(error);
    } catch (const std::bad_alloc&) {
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, evaluate});
    } catch (...) {
        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved, evaluate});
    }
}

DiscreteCirculationPullback::DiscreteCirculationPullback(
    DiscreteCirculation&& primal, DiscreteCirculationLinearization&& transpose,
    core::PlannedBuffer<std::optional<CirculationPointCovector>> points,
    core::PlannedBuffer<std::optional<CirculationCollocationCovector>> rows,
    const CirculationVectorCovector& uniform_flow, std::uint64_t steps) noexcept
    : primal_(std::move(primal)),
      transpose_(std::move(transpose)),
      points_(std::move(points)),
      rows_(std::move(rows)),
      uniform_flow_(uniform_flow),
      steps_(steps) {}

core::Result<CirculationPointCovector>
DiscreteCirculationPullback::point(std::size_t i) const noexcept {
    if (i >= points_.size())
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    const auto& point = points_.values()[i];
    if (!point)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    return *point;
}

core::Result<CirculationCollocationCovector>
DiscreteCirculationPullback::row(std::size_t i) const noexcept {
    if (i >= rows_.size())
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    const auto& row = rows_.values()[i];
    if (!row)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    return *row;
}

core::Result<DiscreteCirculationPullback> solve_discrete_circulation_pullback(
    const DiscreteNonpenetrationInput& input, std::span<const CirculationLinearRhs> seed,
    ProjectedCaseUniformFlow&& flow, const DiscreteNonpenetrationLimits& limits,
    const DiscreteCirculationRequest& primal, const DiscreteCirculationRequest& transpose,
    a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
    core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    try {
        return NonpenetrationBuilder::solve_pullback(input, seed, std::move(flow), limits, primal,
                                                     transpose, shared, local, kernel, memory,
                                                     stop);
    } catch (const p::ArithmeticFailure& error) {
        return std::unexpected(core::Failure{error.code, evaluate});
    } catch (const core::Failure& error) {
        return std::unexpected(error);
    } catch (const std::bad_alloc&) {
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, evaluate});
    } catch (...) {
        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved, evaluate});
    }
}

DiscreteCirculationLoads::DiscreteCirculationLoads(
    DiscreteCirculation&& primal, core::PlannedBuffer<std::optional<DiscreteCirculatoryLoad>> edges,
    core::PlannedBuffer<std::optional<DiscreteNoncirculatoryLoad>> facets,
    core::PlannedBuffer<std::optional<DiscreteNodalLoad>> points, const DiscreteLoadVector& force,
    const DiscreteLoadVector& moment, const DiscreteLoadVector& origin,
    std::uint64_t steps) noexcept
    : primal_(std::move(primal)),
      edges_(std::move(edges)),
      facets_(std::move(facets)),
      points_(std::move(points)),
      force_(force),
      moment_(moment),
      origin_(origin),
      steps_(steps) {}

namespace {
template <class T>
core::Result<T> load_row(const core::PlannedBuffer<std::optional<T>>& rows,
                         std::size_t i) noexcept {
    if (i >= rows.size())
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    const auto& row = rows.values()[i];
    if (!row)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, evaluate});
    return *row;
}
} // namespace

core::Result<DiscreteCirculatoryLoad> DiscreteCirculationLoads::edge(std::size_t i) const noexcept {
    return load_row(edges_, i);
}
core::Result<DiscreteNoncirculatoryLoad>
DiscreteCirculationLoads::facet(std::size_t i) const noexcept {
    return load_row(facets_, i);
}
core::Result<DiscreteNodalLoad> DiscreteCirculationLoads::point(std::size_t i) const noexcept {
    return load_row(points_, i);
}

core::Result<DiscreteCirculationLoads> solve_discrete_circulation_loads(
    const DiscreteNonpenetrationInput& input, const DiscreteVectorLoadInput& loads,
    ProjectedCaseUniformFlow&& flow, const DiscreteNonpenetrationLimits& limits,
    const DiscreteCirculationRequest& request, a::GeometryWorkspace& geometry_work,
    a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
    core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    try {
        return NonpenetrationBuilder::solve_loads(input, loads, std::move(flow), limits, request,
                                                  geometry_work, shared, local, kernel, memory,
                                                  stop);
    } catch (const p::ArithmeticFailure& error) {
        return std::unexpected(core::Failure{error.code, evaluate});
    } catch (const core::Failure& error) {
        return std::unexpected(error);
    } catch (const std::bad_alloc&) {
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, evaluate});
    } catch (...) {
        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved, evaluate});
    }
}

DiscreteCirculationLoadJvp::DiscreteCirculationLoadJvp(
    DiscreteCirculationLoads&& loads,
    core::PlannedBuffer<std::optional<DiscreteLoadSampleDirection>> edges,
    core::PlannedBuffer<std::optional<DiscreteLoadSampleDirection>> facets,
    core::PlannedBuffer<std::optional<DiscreteNodalLoad>> points, const DiscreteLoadVector& force,
    const DiscreteLoadVector& moment) noexcept
    : loads_(std::move(loads)),
      edges_(std::move(edges)),
      facets_(std::move(facets)),
      points_(std::move(points)),
      force_(force),
      moment_(moment) {}

core::Result<DiscreteLoadSampleDirection>
DiscreteCirculationLoadJvp::edge(std::size_t i) const noexcept {
    return load_row(edges_, i);
}
core::Result<DiscreteLoadSampleDirection>
DiscreteCirculationLoadJvp::facet(std::size_t i) const noexcept {
    return load_row(facets_, i);
}
core::Result<DiscreteNodalLoad> DiscreteCirculationLoadJvp::point(std::size_t i) const noexcept {
    return load_row(points_, i);
}

core::Result<DiscreteCirculationLoadJvp> solve_discrete_circulation_load_direction(
    const DiscreteNonpenetrationInput& input, const DiscreteVectorLoadInput& loads,
    const DiscreteVectorLoadDirection& direction, ProjectedCaseUniformFlow&& flow,
    const DiscreteNonpenetrationLimits& limits, const DiscreteCirculationRequest& request,
    a::GeometryWorkspace& geometry_work, a::SharedP1Workspace& shared,
    a::LocalPacketWorkspace& local, a::P1Workspace& kernel, core::ResourcePlan& memory,
    const std::stop_token& stop) noexcept {
    try {
        return NonpenetrationBuilder::solve_load_direction(input, loads, direction, std::move(flow),
                                                           limits, request, geometry_work, shared,
                                                           local, kernel, memory, stop);
    } catch (const p::ArithmeticFailure& error) {
        return std::unexpected(core::Failure{error.code, evaluate});
    } catch (const core::Failure& error) {
        return std::unexpected(error);
    } catch (const std::bad_alloc&) {
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, evaluate});
    } catch (...) {
        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved, evaluate});
    }
}

DiscreteCirculationLoadVjp::DiscreteCirculationLoadVjp(
    DiscreteCirculationLoads&& loads,
    core::PlannedBuffer<std::optional<CirculationPointCovector>> points,
    core::PlannedBuffer<std::optional<LoadEdgeVelocityCovector>> edges,
    core::PlannedBuffer<std::optional<LoadCirculationCovector>> bound,
    const DiscreteLoadVector& flow) noexcept
    : loads_(std::move(loads)),
      points_(std::move(points)),
      edges_(std::move(edges)),
      bound_(std::move(bound)),
      flow_(flow) {}

core::Result<CirculationPointCovector>
DiscreteCirculationLoadVjp::point(std::size_t i) const noexcept {
    return load_row(points_, i);
}
core::Result<LoadEdgeVelocityCovector>
DiscreteCirculationLoadVjp::midpoint_velocity(std::size_t i) const noexcept {
    return load_row(edges_, i);
}
core::Result<LoadCirculationCovector>
DiscreteCirculationLoadVjp::bound(std::size_t i) const noexcept {
    return load_row(bound_, i);
}

core::Result<DiscreteCirculationLoadVjp> solve_discrete_circulation_load_pullback(
    const DiscreteNonpenetrationInput& input, const DiscreteVectorLoadInput& loads,
    const DiscreteVectorLoadWeights& weights, ProjectedCaseUniformFlow&& flow,
    const DiscreteNonpenetrationLimits& limits, const DiscreteCirculationRequest& request,
    a::GeometryWorkspace& geometry_work, a::SharedP1Workspace& shared,
    a::LocalPacketWorkspace& local, a::P1Workspace& kernel, core::ResourcePlan& memory,
    const std::stop_token& stop) noexcept {
    try {
        return NonpenetrationBuilder::solve_load_pullback(input, loads, weights, std::move(flow),
                                                          limits, request, geometry_work, shared,
                                                          local, kernel, memory, stop);
    } catch (const p::ArithmeticFailure& error) {
        return std::unexpected(core::Failure{error.code, evaluate});
    } catch (const core::Failure& error) {
        return std::unexpected(error);
    } catch (const std::bad_alloc&) {
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, evaluate});
    } catch (...) {
        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved, evaluate});
    }
}

core::Result<LoadCirculationRateCovector>
DiscreteCirculationImplicitLoadVjp::circulation_rate(std::size_t i) const noexcept {
    const auto dual = values_.bound(i);
    if (!dual)
        return std::unexpected(dual.error());
    return LoadCirculationRateCovector{dual->face, dual->circulation_rate};
}

core::Result<CirculationCollocationCovector>
DiscreteCirculationImplicitLoadVjp::collocation(std::size_t i) const noexcept {
    return load_row(rows_, i);
}

core::Result<DiscreteCirculationImplicitLoadVjp> solve_discrete_circulation_implicit_load_pullback(
    const DiscreteNonpenetrationInput& input, const DiscreteVectorLoadInput& loads,
    const DiscreteVectorLoadWeights& weights, ProjectedCaseUniformFlow&& flow,
    const DiscreteNonpenetrationLimits& limits, const DiscreteCirculationRequest& primal,
    const DiscreteCirculationRequest& transpose, a::GeometryWorkspace& geometry_work,
    a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
    core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    try {
        return NonpenetrationBuilder::solve_implicit_load_pullback(
            input, loads, weights, std::move(flow), limits, primal, transpose, geometry_work,
            shared, local, kernel, memory, stop);
    } catch (const p::ArithmeticFailure& error) {
        return std::unexpected(core::Failure{error.code, evaluate});
    } catch (const core::Failure& error) {
        return std::unexpected(error);
    } catch (const std::bad_alloc&) {
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, evaluate});
    } catch (...) {
        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved, evaluate});
    }
}

core::Result<DiscreteCirculationImplicitLoadJvp> solve_discrete_circulation_implicit_load_direction(
    const DiscreteNonpenetrationInput& input, const DiscreteVectorLoadInput& loads,
    const DiscreteImplicitVectorLoadDirection& variation, ProjectedCaseUniformFlow&& flow,
    const DiscreteNonpenetrationLimits& limits, const DiscreteCirculationRequest& primal,
    const DiscreteCirculationRequest& tangent, a::GeometryWorkspace& geometry_work,
    a::SharedP1Workspace& shared, a::LocalPacketWorkspace& local, a::P1Workspace& kernel,
    core::ResourcePlan& memory, const std::stop_token& stop) noexcept {
    try {
        return NonpenetrationBuilder::solve_implicit_load_direction(
            input, loads, variation, std::move(flow), limits, primal, tangent, geometry_work,
            shared, local, kernel, memory, stop);
    } catch (const p::ArithmeticFailure& error) {
        return std::unexpected(core::Failure{error.code, evaluate});
    } catch (const core::Failure& error) {
        return std::unexpected(error);
    } catch (const std::bad_alloc&) {
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, evaluate});
    } catch (...) {
        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved, evaluate});
    }
}

DiscreteKnownCirculation::DiscreteKnownCirculation(
    core::PlannedBuffer<KnownCirculationRun>&& data) noexcept
    : data_(std::move(data)) {}
DiscreteKnownCirculation::DiscreteKnownCirculation(DiscreteKnownCirculation&&) noexcept = default;
DiscreteKnownCirculation&
DiscreteKnownCirculation::operator=(DiscreteKnownCirculation&&) noexcept = default;
DiscreteKnownCirculation::~DiscreteKnownCirculation() = default;

namespace {
core::Failure absent() noexcept {
    return {core::FailureCode::InvalidInput, evaluate};
}
} // namespace
core::Result<std::reference_wrapper<const linear::Solution<double>>>
DiscreteKnownCirculation::solution() const& noexcept {
    if (data_.values().empty())
        return std::unexpected(absent());
    const auto& primal = data_.values().front().primal;
    if (!primal)
        return std::unexpected(absent());
    return std::cref(primal->solution);
}
core::Result<std::reference_wrapper<const DiscreteCirculationLinearization>>
DiscreteKnownCirculation::linearization() const& noexcept {
    if (data_.values().empty())
        return std::unexpected(absent());
    const auto& result = data_.values().front().linearization;
    if (!result)
        return std::unexpected(absent());
    return std::cref(*result);
}
core::Result<CirculationResidualRow> DiscreteKnownCirculation::row(std::size_t i) const noexcept {
    if (data_.values().empty())
        return std::unexpected(absent());
    const auto& primal = data_.values().front().primal;
    if (!primal || i >= primal->rows.size())
        return std::unexpected(absent());
    const auto& row = primal->rows.values()[i];
    if (!row)
        return std::unexpected(absent());
    return *row;
}
core::Result<std::reference_wrapper<const CaseKnownFlowJet>>
DiscreteKnownCirculation::flow(std::size_t i) const& noexcept {
    if (data_.values().empty() || i >= data_.values().front().known.size())
        return std::unexpected(absent());
    const auto& jet = data_.values().front().known.values()[i].jet;
    if (!jet)
        return std::unexpected(absent());
    return std::cref(*jet);
}
core::Result<CirculationForcingRow>
DiscreteKnownCirculation::forcing(std::size_t i) const noexcept {
    if (data_.values().empty() || i >= data_.values().front().forcing.size())
        return std::unexpected(absent());
    const auto& result = data_.values().front().forcing.values()[i];
    if (!result)
        return std::unexpected(absent());
    return *result;
}
core::Result<CirculationPointCovector>
DiscreteKnownCirculation::point(std::size_t i) const noexcept {
    if (data_.values().empty() || i >= data_.values().front().points.size())
        return std::unexpected(absent());
    const auto& result = data_.values().front().points.values()[i];
    if (!result)
        return std::unexpected(absent());
    return *result;
}
core::Result<CirculationCollocationCovector>
DiscreteKnownCirculation::collocation(std::size_t i) const noexcept {
    if (data_.values().empty() || i >= data_.values().front().dual_rows.size())
        return std::unexpected(absent());
    const auto& result = data_.values().front().dual_rows.values()[i];
    if (!result)
        return std::unexpected(absent());
    return *result;
}
core::Result<CirculationVectorCovector>
DiscreteKnownCirculation::uniform_component_dual() const noexcept {
    if (data_.values().empty())
        return std::unexpected(absent());
    const auto& result = data_.values().front().uniform;
    if (!result)
        return std::unexpected(absent());
    return *result;
}
core::Result<a::ScaledEnclosure> DiscreteKnownCirculation::flow_time_dual() const noexcept {
    if (data_.values().empty())
        return std::unexpected(absent());
    const auto& result = data_.values().front().time;
    if (!result)
        return std::unexpected(absent());
    return *result;
}

core::Result<DiscreteKnownCirculation> solve_discrete_known_circulation(
    const DiscreteNonpenetrationInput& input, const DiscreteKnownFlowInput& flow,
    const KnownCirculationDerivative& derivative, const DiscreteNonpenetrationLimits& limits,
    const DiscreteCirculationRequest& primal, a::SharedP1Workspace& shared,
    a::LocalPacketWorkspace& local, a::P1Workspace& kernel, core::ResourcePlan& memory,
    const std::stop_token& stop) noexcept {
    try {
        return NonpenetrationBuilder::solve_known(input, flow, derivative, limits, primal, shared,
                                                  local, kernel, memory, stop);
    } catch (const p::ArithmeticFailure& error) {
        return std::unexpected(core::Failure{error.code, evaluate});
    } catch (const core::Failure& error) {
        return std::unexpected(error);
    } catch (const std::bad_alloc&) {
        return std::unexpected(core::Failure{core::FailureCode::ResourceExceeded, evaluate});
    } catch (...) {
        return std::unexpected(core::Failure{core::FailureCode::PrecisionUnresolved, evaluate});
    }
}
} // namespace uvlas::coupling::detail
