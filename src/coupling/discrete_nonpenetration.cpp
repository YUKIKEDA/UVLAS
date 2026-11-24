#include "discrete_nonpenetration.hpp"

#include "aero/local_geometry_bounds_private.hpp"
#include "aero/p1_interval_private.hpp"
#include "aero/p1_shared_kernel.hpp"

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
    actual.cache.scales = primal.cache.scales;
    if (transposed) {
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
    if (actual != primal)
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
};

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
                             const ProjectedCaseUniformFlow& flow,
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
        const auto& jet = flow.source().jet;
        if (require(jet.stamp()) != input.flow_stamp ||
            flow.source().requested_side != input.flow_side ||
            !same_time(geometry, require(jet.time()), work))
            reject(core::FailureCode::StaleGeneration);
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
        const auto assembly = assemble(input, flow, limits, memory, work, false);
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
                 core::ResourcePlan& memory, a::LocalRootScratch& rounding, p::Work& work) {
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
        const auto representation = limits.induction_accuracy.representation;
        std::optional<CirculationVectorCovector> enclosed_query;
        if (query)
            enclosed_query.emplace(enclose_vector(*query, representation, rounding, work));
        const auto x_dual = cross(y, area_dual, work);
        const auto y_dual = cross(area_dual, x, work);
        return {panel.panel.face,
                std::move(enclosed_query),
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

    static PrimalCirculationData
    solve_data(const DiscreteNonpenetrationInput& input, const ProjectedCaseUniformFlow& flow,
               const DiscreteNonpenetrationLimits& limits,
               const DiscreteCirculationRequest& request, a::SharedP1Workspace& shared,
               a::LocalPacketWorkspace& local, a::P1Workspace& kernel, core::ResourcePlan& memory,
               p::Work& work, const DiscreteCirculationDirection* direction,
               std::span<p::Interval> forcing, PullbackData* pullback = nullptr) {
        const auto& stop = work.stop;
        work.step();
        const auto geometry_owner = require(input.contact.geometry());
        const auto& geometry = geometry_owner.get();
        const auto order = input.collocations.size();
        match_circulation_geometry(request, geometry.context());
        validate_circulation(request, order, work);
        auto [matrix, rhs] = dense_inputs(order, memory);
        auto assembly = assemble(input, flow, limits, memory, work, true);
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

        const auto ambient = ambient_velocity(flow, work);
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
            record.relative = difference(ambient, import(row.velocity, work), work);
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
            if (direction != nullptr)
                forcing[i] =
                    fixed_circulation_direction(row, direction->rows[i], record.normal,
                                                record.relative, ambient_direction, full, work);
            const auto original = native_residual(record.normal, record.relative, full, work);
            check.add(original, record.rhs, i, work);
            if (pullback != nullptr)
                pullback->rows.values()[i].emplace(
                    pullback_row(*pullback, i, item, row, record, full, action, source, assembly,
                                 limits, shared, memory, rounding, work));
            const auto enclosed =
                a::enclose_p1_value(original, limits.induction_accuracy.representation, rounding);
            results.values()[i].emplace(CirculationResidualRow{
                item.panel.face, enclosed, require(action.report()).far_facets != 0});
            if (i == 0)
                geometry_parent.emplace(std::move(full));
        }
        check.finish();
        work.step();
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
        auto data = solve_data(input, flow, limits, request, shared, local, kernel, memory, work,
                               nullptr, {});
        return finish_primal(std::move(data), std::move(flow));
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
        match_linearization_context(tangent_request.scales.context, primal_request.scales.context,
                                    false);
        validate_circulation(tangent_request, order, work);
        auto native_rhs = allocate<p::Interval>(memory, order);
        auto forcing = allocate<std::optional<CirculationForcingRow>>(memory, order);
        auto primal = solve_data(input, flow, limits, primal_request, shared, local, kernel, memory,
                                 work, &direction, native_rhs.values());
        a::LocalRootScratch rounding;
        for (std::size_t i = 0; i < order; ++i) {
            const auto& row = primal.rows.values()[i];
            if (!row)
                reject(core::FailureCode::InvalidInput);
            forcing.values()[i].emplace(CirculationForcingRow{
                row->face,
                a::enclose_p1_value(native_rhs.values()[i],
                                    limits.induction_accuracy.representation, rounding)});
            native_rhs.values()[i] = p::negate(native_rhs.values()[i], work);
        }
        const auto before_tangent = work.used;
        const auto face = [&](std::size_t i) { return face_at(primal.rows.values(), i); };
        auto tangent = linearize_native(CirculationSense::Tangent, face, primal.influence.values(),
                                        native_rhs.values(), tangent_request, memory, work);
        work.step();
        tangent.steps_ = work.used - before_tangent;
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
        auto primal = solve_data(input, flow, limits, primal_request, shared, local, kernel, memory,
                                 work, nullptr, {}, &target);
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
        const auto previous = primal.solution_.certificate().context();
        match_linearization_context(request.scales.context, previous, transposed);
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

} // namespace uvlas::coupling::detail
