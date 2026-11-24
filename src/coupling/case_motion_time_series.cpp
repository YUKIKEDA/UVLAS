#include "case_motion_time_series.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <type_traits>
#include <utility>

namespace uvlas::coupling::detail {
namespace p = aero::detail::pchip;
namespace s = io::schema;
namespace motion_time_detail {
namespace {

using Rational = structure::ExactRational;
using Words = io::detail::SiRationalWords;
using View = p::RationalView;
using Code = core::FailureCode;
using Stamp = p::PchipInputStamp;
constexpr auto preparation = core::Operation::PrepareModel;
constexpr auto evaluation = core::Operation::EvaluateStep;
constexpr auto index_max = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
constexpr io::SchemaDimension seconds{{0, 0, 1, 0}};
constexpr io::SchemaDimension metres{{1, 0, 0, 0}};
constexpr io::SchemaDimension radians{{0, 0, 0, 1}};

template <class T>
T take(core::Result<T> result) {
    if (!result)
        throw result.error();
    if constexpr (!std::is_void_v<T>)
        return std::move(*result);
}

template <class T, class Function>
core::Result<T> guarded(core::Operation op, core::FailureDetails details,
                        const Function& call) noexcept {
    return core::guard_entry(op, Code::InvalidInput, details, [&]() -> core::Result<T> {
        try {
            return call();
        } catch (const core::Failure& error) {
            return std::unexpected(error);
        }
    });
}

View rational_view(const Words& value) noexcept {
    const auto view = value.view();
    return {view.sign, view.numerator, view.denominator};
}

bool same_source(const io::SchemaIdentity& left, const io::SchemaIdentity& right) noexcept {
    return left.source == right.source && left.contract_revision == right.contract_revision &&
           left.representation_revision == right.representation_revision;
}

bool valid_order(p::JetOrder order) noexcept {
    return order == p::JetOrder::Value || order == p::JetOrder::First ||
           order == p::JetOrder::Second;
}

bool valid_side(p::RequestedSide side) noexcept {
    return side == p::RequestedSide::Automatic || side == p::RequestedSide::Left ||
           side == p::RequestedSide::Right;
}

struct SeriesDomain {
    Rational origin, scale, begin, end;
};

struct FormulaSeries {
    core::PlannedBuffer<Rational> coefficients;
    std::optional<SeriesDomain> polynomial;
};

using Series = std::variant<FormulaSeries, p::PchipCurveReady>;

struct WordsAndViews {
    core::PlannedBuffer<Words> words;
    core::PlannedBuffer<View> views;
};

class Call {
public:
    Call(const CaseMotionTimeLimits& limits, core::ResourcePlan& plan, const std::stop_token& stop,
         core::Operation operation, core::FailureDetails context) noexcept
        : exact(limits.conversion),
          curve(limits.curve),
          memory(plan),
          cancellation(stop),
          op(operation),
          details(context),
          operations_(limits.maximum_rational_operations) {}

    [[noreturn]] void fail(Code code) const {
        throw core::Failure{code, op, details};
    }

    void stop() const {
        if (cancellation.stop_requested())
            fail(Code::Cancelled);
    }

    [[nodiscard]] std::uint64_t sum(std::uint64_t left, std::uint64_t right) const {
        if (left > index_max || right > index_max - left)
            fail(Code::ResourceExceeded);
        return left + right;
    }

    [[nodiscard]] std::uint64_t times_three(std::uint64_t count) const {
        if (count > index_max / 3)
            fail(Code::ResourceExceeded);
        return 3 * count;
    }

    void operation() {
        stop();
        if (operations_ == 0)
            fail(Code::ResourceExceeded);
        --operations_;
    }

    void steps(std::uint64_t count) {
        stop();
        if (count > curve.arithmetic.steps)
            fail(Code::ResourceExceeded);
        curve.arithmetic.steps -= count;
    }

    p::CurveLimits reserve(const p::CurveAdmission& admission) {
        if (admission.payload_bytes > memory.remaining_bytes() ||
            admission.reserved_steps > curve.arithmetic.steps / 2)
            fail(Code::ResourceExceeded);
        auto execution = curve;
        execution.arithmetic.steps -= admission.reserved_steps;
        steps(2 * admission.reserved_steps);
        const auto selection = curve.metadata.selection();
        if (const auto* estimate = std::get_if<numerics::linear::EstimateAdmission>(&selection)) {
            if (admission.metadata_estimate_bytes > estimate->bytes)
                fail(Code::ResourceExceeded);
            curve.metadata = numerics::linear::EstimateAdmission{estimate->bytes -
                                                                 admission.metadata_estimate_bytes};
        }
        return execution;
    }

    std::strong_ordering compare(const Rational& left, const Rational& right) {
        operation();
        return take(left.compare(right, exact, memory));
    }

    Rational clone(const Rational& value) {
        operation();
        return take(value.clone(exact, memory));
    }

    Rational add(const Rational& left, const Rational& right) {
        operation();
        return take(left.add(right, exact, memory));
    }

    Rational subtract(const Rational& left, const Rational& right) {
        operation();
        return take(left.subtract(right, exact, memory));
    }

    Rational multiply(const Rational& left, const Rational& right) {
        operation();
        return take(left.multiply(right, exact, memory));
    }

    Rational divide(const Rational& left, const Rational& right) {
        operation();
        return take(left.divide(right, exact, memory));
    }

    [[nodiscard]] Words copy(const Rational& value) const {
        return take(io::detail::copy_rational_words(value, memory, cancellation));
    }

    template <class Field>
    Rational read(const io::SiCase& source, Field field, io::SchemaDimension dimensions,
                  std::int16_t pi_power) {
        operation();
        auto recipe = take(source.quantity(field)).recipe;
        if (recipe.dimensions.exponents != dimensions.exponents || recipe.pi_power != pi_power)
            fail(Code::InvalidInput);
        if (pi_power == 0)
            return take(io::materialize_rational(recipe, exact, memory));
        // Materialize the coefficient of pi with the same validated recipe path as scalar gust.
        recipe.units.angle = io::SiAngleUnit::Radian;
        recipe.pi_power = 0;
        if (!degree_denominator_) {
            operation();
            degree_denominator_.emplace(take(Rational::from_decimal("180", exact, memory)));
        }
        const auto raw = take(io::materialize_rational(recipe, exact, memory));
        return divide(raw, *degree_denominator_);
    }

    template <class T>
    [[nodiscard]] core::PlannedBuffer<T> allocate(std::uint64_t count) const {
        stop();
        if (count > index_max)
            fail(Code::ResourceExceeded);
        return take(memory.allocate<T>(static_cast<core::Index>(count), op));
    }

    [[nodiscard]] WordsAndViews allocate_words(std::uint64_t count) const {
        constexpr auto width = sizeof(Words) + sizeof(View);
        if (count > memory.remaining_bytes() / width)
            fail(Code::ResourceExceeded);
        return {allocate<Words>(count), allocate<View>(count)};
    }

    const structure::ExactLimits& exact;
    p::CurveLimits curve;
    core::ResourcePlan& memory;
    const std::stop_token& cancellation;
    core::Operation op;
    core::FailureDetails details;

private:
    std::uint64_t operations_;
    std::optional<Rational> degree_denominator_;
};

void store(WordsAndViews& storage, std::size_t index, const Rational& value, const Call& call) {
    auto& words = storage.words.values()[index];
    words = call.copy(value);
    storage.views.values()[index] = rational_view(words);
}

template <class Array>
void read_row(const io::SiCase& source, Array row, io::SchemaDimension dimension,
              std::int16_t pi_power, std::span<Rational> out, Call& call) {
    for (core::Index c = 0; c < 3; ++c) {
        out[static_cast<std::size_t>(c)] = call.read(source, take(row.at(c)), dimension, pi_power);
    }
}

using TypedSeries = decltype(std::declval<s::Motion>().translation.record());
using BoundaryPolicy = decltype(p::BoundaryRequest::policy);

std::uint64_t recipe_count(const TypedSeries& series, const Call& call) {
    return std::visit(
        [&](const auto& typed) -> std::uint64_t {
            using T = std::remove_cvref_t<decltype(typed)>;
            if constexpr (std::is_same_v<T, s::TimeSeriesConstant>) {
                return 3;
            } else if constexpr (std::is_same_v<T, s::TimeSeriesPolynomial>) {
                return call.sum(4, call.times_three(typed.coefficients.size()));
            } else if constexpr (std::is_same_v<T, s::TimeSeriesPchip>) {
                const auto period =
                    std::holds_alternative<s::OutsidePeriodic>(typed.outside.record());
                return call.sum(call.sum(typed.times.size(), call.times_three(typed.values.size())),
                                period ? 1U : 0U);
            } else {
                call.fail(Code::Unsupported);
            }
        },
        series);
}

BoundaryPolicy read_outside(const io::SiCase& source, const s::TimeSeriesPchip& typed,
                            WordsAndViews& storage, std::size_t last, Call& call) {
    const auto outside = typed.outside.record();
    if (const auto* periodic = std::get_if<s::OutsidePeriodic>(&outside)) {
        const auto period = call.read(source, periodic->period, seconds, 0);
        store(storage, last, period, call);
        return p::PeriodicBoundary{storage.views.values()[last]};
    }
    const auto tag = std::get<s::OutsideRejectHold>(outside).kind.selected();
    if (tag == decltype(tag)::ItemHold)
        return p::HoldBoundary{};
    return p::RejectBoundary{};
}

p::BitBound bit_shape(View value, const Call& call) {
    const auto cap = call.curve.arithmetic.input_bits;
    const auto bits = [&](std::span<const std::uint32_t> words) -> std::uint64_t {
        if (words.empty())
            return 0;
        if (words.size() - 1 > cap / 32)
            call.fail(Code::ResourceExceeded);
        const auto preceding = std::uint64_t{32} * (words.size() - 1);
        const auto leading = static_cast<std::uint64_t>(std::bit_width(words.back()));
        if (leading > cap - preceding)
            call.fail(Code::ResourceExceeded);
        return preceding + leading;
    };
    return {bits(value.numerator), bits(value.denominator)};
}

// For each source phase, the first occurrence strictly after begin is enough to
// decide whether any occurrence lies in the open interval, even over many periods.
template <class Arithmetic, class Value>
Value next_periodic_knot(Arithmetic& math, Value begin, Value knot, Value period) {
    using enum p::RationalOp;
    const auto offset = take(math.binary(Subtract, begin, knot));
    const auto quotient = take(math.binary(Divide, offset, period));
    const auto whole = take(math.floor(quotient));
    const auto one = take(math.integer(1));
    const auto next = take(math.binary(Add, whole, one));
    const auto shift = take(math.binary(Multiply, next, period));
    return take(math.binary(Add, knot, shift));
}

bool interior_motion_knot(const p::BoundaryRequest& window, View knot, Call& call) {
    call.operation();
    const auto* periodic = std::get_if<p::PeriodicBoundary>(&window.policy);
    p::Counter forecast{call.curve.arithmetic, call.details};
    const auto a = take(forecast.input(bit_shape(window.begin, call)));
    const auto b = take(forecast.input(bit_shape(window.end, call)));
    const auto t = take(forecast.input(bit_shape(knot, call)));
    if (periodic) {
        const auto period = take(forecast.input(bit_shape(periodic->period, call)));
        take(forecast.compare(next_periodic_knot(forecast, a, t, period), b));
    } else {
        take(forecast.compare(a, t));
        take(forecast.compare(t, b));
    }
    const auto admission = take(forecast.finish());
    call.steps(admission.reserved_steps());
    auto arena =
        take(p::Workspace::create(admission, call.memory, call.details, call.cancellation));
    auto math = take(arena.begin(call.cancellation, call.details));
    const auto begin = take(math.input(window.begin));
    const auto end = take(math.input(window.end));
    const auto time = take(math.input(knot));
    if (periodic) {
        const auto period = take(math.input(periodic->period));
        return take(math.compare(next_periodic_knot(math, begin, time, period), end)) < 0;
    }
    return take(math.compare(begin, time)) < 0 && take(math.compare(time, end)) < 0;
}

void check_c2(const p::PchipCurveReady& curve, p::CurveWorkspace& workspace, Call& call) {
    try {
        const auto input = take(curve.input());
        const auto expected = take(curve.retain_input());
        auto knots = input.samples.times;
        if (std::holds_alternative<p::PeriodicBoundary>(input.boundary.policy))
            knots = knots.first(knots.size() - 1); // The last knot repeats the first seam.
        else if (std::holds_alternative<p::RejectBoundary>(input.boundary.policy))
            knots = knots.subspan(1, knots.size() - 2);
        for (const auto knot : knots) {
            if (!interior_motion_knot(input.boundary, knot, call))
                continue;
            const p::PchipQuery query{knot, p::JetOrder::Second, p::RequestedSide::Automatic};
            const auto admission = take(p::forecast_pchip_evaluation(
                curve, expected, query, call.curve, call.cancellation, call.details));
            static_cast<void>(take(p::evaluate_pchip_curve(curve, expected, query, workspace,
                                                           call.reserve(admission), call.memory,
                                                           call.cancellation, call.details)));
        }
    } catch (const core::Failure& failure) {
        throw core::Failure{failure.code, preparation, failure.details};
    }
}

Series prepare_series(const io::SiCase& source, const TypedSeries& series,
                      io::SchemaDimension dimension, std::int16_t pi_power,
                      CaseMotionTimeWindow window, std::span<const View> events, Stamp stamp,
                      p::CurveWorkspace& workspace, Call& call) {
    return std::visit(
        [&](const auto& typed) -> Series {
            using T = std::remove_cvref_t<decltype(typed)>;
            if constexpr (std::is_same_v<T, s::TimeSeriesConstant>) {
                auto coefficients = call.allocate<Rational>(3);
                read_row(source, typed.value, dimension, pi_power, coefficients.values(), call);
                return FormulaSeries{std::move(coefficients), {}};
            } else if constexpr (std::is_same_v<T, s::TimeSeriesPolynomial>) {
                SeriesDomain domain{call.read(source, typed.origin, seconds, 0),
                                    call.read(source, typed.scale, seconds, 0),
                                    call.read(source, take(typed.domain.at(0)), seconds, 0),
                                    call.read(source, take(typed.domain.at(1)), seconds, 0)};
                if (call.compare(window.begin, domain.begin) < 0 ||
                    call.compare(window.end, domain.end) > 0)
                    call.fail(Code::InvalidInput);
                const auto count = typed.coefficients.size();
                auto coefficients = call.allocate<Rational>(call.times_three(count));
                for (core::Index i = 0; std::cmp_less(i, count); ++i)
                    read_row(source, take(typed.coefficients.at(i)), dimension, pi_power,
                             coefficients.values().subspan(static_cast<std::size_t>(i) * 3, 3),
                             call);
                return FormulaSeries{std::move(coefficients), std::move(domain)};
            } else if constexpr (std::is_same_v<T, s::TimeSeriesPchip>) {
                const auto count = typed.times.size();
                const auto sample_count = call.times_three(count);
                const auto samples_end = call.sum(count, sample_count);
                auto storage = call.allocate_words(recipe_count(series, call));
                for (core::Index i = 0; std::cmp_less(i, count); ++i) {
                    const auto index = static_cast<std::size_t>(i);
                    const auto time = call.read(source, take(typed.times.at(i)), seconds, 0);
                    store(storage, index, time, call);
                    const auto row = take(typed.values.at(i));
                    std::array<Rational, 3> values;
                    read_row(source, row, dimension, pi_power, values, call);
                    for (std::size_t c = 0; c < values.size(); ++c)
                        store(storage, static_cast<std::size_t>(count) + index * 3 + c, values[c],
                              call);
                }
                const auto boundary = read_outside(source, typed, storage,
                                                   static_cast<std::size_t>(samples_end), call);
                const auto begin = call.copy(window.begin);
                const auto end = call.copy(window.end);
                // Radian is dimensionless in the SI base-unit descriptor. The Motion wrapper
                // separately carries the angle role, exponent and source units.
                const auto length = static_cast<std::int8_t>(dimension.exponents[0]);
                const numerics::SiDimension scalar{{length, 0, 0, 0, 0, 0, 0}};
                const std::array dimensions{scalar, scalar, scalar};
                const auto views = storage.views.values();
                const p::PreparationRequest request{
                    {views.first(static_cast<std::size_t>(count)),
                     views.subspan(static_cast<std::size_t>(count),
                                   static_cast<std::size_t>(sample_count)),
                     3},
                    dimensions,
                    {rational_view(begin), rational_view(end), boundary},
                    events,
                    stamp};
                const auto admission = take(p::forecast_pchip_prepare(
                    request, call.curve, call.cancellation, call.details));
                auto curve =
                    take(p::prepare_pchip_curve(request, workspace, call.reserve(admission),
                                                call.memory, call.cancellation, call.details));
                if (window.continuity == CaseMotionContinuity::ContinuousC2)
                    check_c2(curve, workspace, call);
                return curve;
            } else {
                call.fail(Code::Unsupported);
            }
        },
        series);
}

void match(const FormulaSeries& actual, const FormulaSeries& current, Call& call) {
    if (actual.coefficients.size() != current.coefficients.size() ||
        actual.polynomial.has_value() != current.polynomial.has_value())
        call.fail(Code::StaleGeneration);
    const auto require_equal = [&](const Rational& left, const Rational& right) {
        if (call.compare(left, right) != 0)
            call.fail(Code::StaleGeneration);
    };
    if (actual.polynomial && current.polynomial) {
        require_equal(actual.polynomial->origin, current.polynomial->origin);
        require_equal(actual.polynomial->scale, current.polynomial->scale);
        require_equal(actual.polynomial->begin, current.polynomial->begin);
        require_equal(actual.polynomial->end, current.polynomial->end);
    }
    for (std::size_t i = 0; i < actual.coefficients.size(); ++i)
        require_equal(actual.coefficients.values()[i], current.coefficients.values()[i]);
}

std::array<Words, 9> formula_jet(const FormulaSeries& formula, const Rational& time,
                                 p::JetOrder order, Call& call) {
    std::array<Words, 9> result;
    const auto requested = static_cast<std::size_t>(order);
    if (!formula.polynomial) {
        const Rational zero;
        for (std::size_t k = 0; k <= requested; ++k)
            for (std::size_t c = 0; c < 3; ++c)
                result[3 * k + c] = call.copy(k == 0 ? formula.coefficients.values()[c] : zero);
        return result;
    }

    const auto& domain = *formula.polynomial;
    const auto shifted = call.subtract(time, domain.origin);
    const auto u = call.divide(shifted, domain.scale);
    const auto scale_squared =
        requested == 2 ? call.multiply(domain.scale, domain.scale) : Rational{};
    const auto coefficients = formula.coefficients.values();
    for (std::size_t c = 0; c < 3; ++c) {
        Rational value, first, second;
        for (std::size_t i = coefficients.size() / 3; i != 0; --i) {
            if (requested == 2) {
                const auto twice_first = call.add(first, first);
                second = call.add(call.multiply(u, second), twice_first);
            }
            if (requested >= 1)
                first = call.add(call.multiply(u, first), value);
            value = call.add(call.multiply(u, value), coefficients[3 * (i - 1) + c]);
        }
        result[c] = call.copy(value);
        if (requested >= 1)
            result[3 + c] = call.copy(call.divide(first, domain.scale));
        if (requested == 2)
            result[6 + c] = call.copy(call.divide(second, scale_squared));
    }
    return result;
}

} // namespace

struct Owner {
    std::optional<io::SiQuantityKey> declaration;
    std::uint64_t declaration_ordinal{};
    Stamp stamp{p::PchipIssuerId{0}, p::PchipInputGeneration{0}};
    std::int16_t pi_power{};
    Rational begin, end;
    CaseMotionContinuity continuity = CaseMotionContinuity::PiecewiseC2;
    WordsAndViews events;
    Series translation, rotation;
};

struct Builder {
    static CaseMotionTimeSeries prepare(const io::SiCase& source, CaseMotionField field,
                                        CaseMotionTimeWindow window, Stamp stamp,
                                        p::CurveWorkspace& workspace,
                                        const CaseMotionTimeLimits& limits, Call& call) {
        call.stop();
        const auto motion = field.record();
        const auto anchor = take(source.quantity(take(motion.origin.record().origin.at(0))));
        const auto pi_power =
            static_cast<std::int16_t>(anchor.recipe.units.angle == io::SiAngleUnit::Degree);
        const auto translation = motion.translation.record();
        const auto rotation = motion.rotation_vector.record();
        auto count = call.sum(recipe_count(translation, call), recipe_count(rotation, call));
        count = call.sum(count, motion.events.size());
        if (count > limits.maximum_exact_numbers)
            call.fail(Code::ResourceExceeded);
        if ((window.continuity != CaseMotionContinuity::PiecewiseC2 &&
             window.continuity != CaseMotionContinuity::ContinuousC2) ||
            call.compare(window.begin, window.end) > 0)
            call.fail(Code::InvalidInput);
        auto storage = call.allocate<Owner>(1);
        auto& owner = storage.values().front();
        owner.declaration.emplace(anchor.key);
        owner.declaration_ordinal = anchor.ledger_origin;
        owner.stamp = stamp;
        owner.pi_power = pi_power;
        owner.begin = call.clone(window.begin);
        owner.end = call.clone(window.end);
        owner.continuity = window.continuity;
        owner.events = call.allocate_words(motion.events.size());
        for (core::Index i = 0; std::cmp_less(i, motion.events.size()); ++i) {
            const auto event = call.read(source, take(motion.events.at(i)), seconds, 0);
            store(owner.events, static_cast<std::size_t>(i), event, call);
        }
        const auto events = owner.events.views.values();
        owner.translation =
            prepare_series(source, translation, metres, 0, window, events, stamp, workspace, call);
        owner.rotation = prepare_series(source, rotation, radians, pi_power, window, events, stamp,
                                        workspace, call);
        call.stop();
        return CaseMotionTimeSeries{std::move(storage)};
    }

    static CaseMotionTimeJet evaluate(const CaseMotionTimeSeries& prepared,
                                      const CaseMotionTimeSeries& current,
                                      CaseMotionTimeQuery request, p::CurveWorkspace& workspace,
                                      Call& call) {
        call.stop();
        if (!valid_order(request.order) || !valid_side(request.side) ||
            prepared.owner_.size() == 0 || current.owner_.size() == 0)
            call.fail(Code::InvalidInput);
        const auto& actual = prepared.owner_.values().front();
        const auto& expected = current.owner_.values().front();
        const auto& declaration = bind(actual, expected, request.stamp, call);
        if (call.compare(request.time, actual.begin) < 0 ||
            call.compare(request.time, actual.end) > 0)
            call.fail(Code::InvalidInput);
        auto instant = call.copy(request.time);
        auto translation = evaluate_series(actual.translation, expected.translation, request,
                                           rational_view(instant), workspace, call);
        auto rotation = evaluate_series(actual.rotation, expected.rotation, request,
                                        rational_view(instant), workspace, call);
        call.stop();
        return CaseMotionTimeJet{
            std::move(translation), std::move(rotation), std::move(instant), declaration,
            request.stamp,          actual.pi_power,     request.order,      request.side};
    }

private:
    static const io::SiQuantityKey& bind(const Owner& actual, const Owner& expected, Stamp stamp,
                                         Call& call) {
        if (stamp != actual.stamp || stamp != expected.stamp ||
            actual.pi_power != expected.pi_power || actual.continuity != expected.continuity ||
            !actual.declaration || !expected.declaration ||
            actual.declaration_ordinal != expected.declaration_ordinal ||
            !same_source(actual.declaration->source(), expected.declaration->source()))
            call.fail(Code::StaleGeneration);
        const auto& declaration = *actual.declaration;
        if (&actual == &expected)
            return declaration;
        if (call.compare(actual.begin, expected.begin) != 0 ||
            call.compare(actual.end, expected.end) != 0 ||
            actual.events.views.size() != expected.events.views.size())
            call.fail(Code::StaleGeneration);
        for (std::size_t i = 0; i < actual.events.views.size(); ++i) {
            const auto left = actual.events.views.values()[i];
            const auto right = expected.events.views.values()[i];
            call.steps(1);
            for (const auto count : {left.numerator.size(), left.denominator.size(),
                                     right.numerator.size(), right.denominator.size()})
                call.steps(count);
            if (left.sign != right.sign || !std::ranges::equal(left.numerator, right.numerator) ||
                !std::ranges::equal(left.denominator, right.denominator))
                call.fail(Code::StaleGeneration);
        }
        return declaration;
    }

    static CaseMotionTimeJet::SeriesJet evaluate_series(const Series& actual, const Series& current,
                                                        CaseMotionTimeQuery request, View time,
                                                        p::CurveWorkspace& workspace, Call& call) {
        if (actual.index() != current.index())
            call.fail(Code::StaleGeneration);
        if (const auto* formula = std::get_if<FormulaSeries>(&actual)) {
            if (&actual != &current)
                match(*formula, std::get<FormulaSeries>(current), call);
            return formula_jet(*formula, request.time, request.order, call);
        }
        const auto& curve = std::get<p::PchipCurveReady>(actual);
        const auto expected = take(std::get<p::PchipCurveReady>(current).retain_input());
        const p::PchipQuery query{time, request.order, request.side};
        const auto admission = take(p::forecast_pchip_evaluation(curve, expected, query, call.curve,
                                                                 call.cancellation, call.details));
        return take(p::evaluate_pchip_curve(curve, expected, query, workspace,
                                            call.reserve(admission), call.memory, call.cancellation,
                                            call.details));
    }
};

} // namespace motion_time_detail

CaseMotionTimeSeries::CaseMotionTimeSeries(
    core::PlannedBuffer<motion_time_detail::Owner> owner) noexcept
    : owner_(std::move(owner)) {}
CaseMotionTimeSeries::CaseMotionTimeSeries(CaseMotionTimeSeries&&) noexcept = default;
CaseMotionTimeSeries& CaseMotionTimeSeries::operator=(CaseMotionTimeSeries&&) noexcept = default;
CaseMotionTimeSeries::~CaseMotionTimeSeries() = default;

core::Result<io::SiQuantityKey> CaseMotionTimeSeries::declaration() const noexcept {
    const auto owners = owner_.values();
    if (!owners.empty()) {
        const auto& declaration = owners.front().declaration;
        if (declaration)
            return *declaration;
    }
    return std::unexpected(
        core::Failure{core::FailureCode::InvalidInput, core::Operation::PrepareModel});
}

std::span<const p::RationalView> CaseMotionTimeSeries::events() const& noexcept {
    return owner_.size() != 0 ? owner_.values().front().events.views.values()
                              : std::span<const p::RationalView>{};
}

CaseMotionTimeJet::CaseMotionTimeJet(SeriesJet translation, SeriesJet rotation,
                                     io::detail::SiRationalWords time,
                                     io::SiQuantityKey declaration, p::PchipInputStamp stamp,
                                     std::int16_t pi_power, p::JetOrder order,
                                     p::RequestedSide side) noexcept
    : translation_(std::move(translation)),
      rotation_(std::move(rotation)),
      time_(std::move(time)),
      declaration_(std::move(declaration)),
      stamp_(stamp),
      pi_power_(pi_power),
      order_(order),
      side_(side) {}

p::RationalView CaseMotionTimeJet::time() const& noexcept {
    return motion_time_detail::rational_view(time_);
}

core::Result<p::RationalView> CaseMotionTimeJet::component(const SeriesJet& series,
                                                           std::uint64_t axis,
                                                           p::JetOrder derivative) const noexcept {
    if (!motion_time_detail::valid_order(derivative) || derivative > order_ || axis >= 3 ||
        time().denominator.empty())
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::EvaluateStep});
    if (const auto* exact = std::get_if<p::ExactPchipJet>(&series))
        return exact->derivative(axis, derivative);
    const auto index = 3 * static_cast<std::size_t>(derivative) + static_cast<std::size_t>(axis);
    return motion_time_detail::rational_view(std::get<FormulaValues>(series)[index]);
}

core::Result<p::RationalView>
CaseMotionTimeJet::translation(std::uint64_t axis, p::JetOrder derivative) const& noexcept {
    return component(translation_, axis, derivative);
}

core::Result<p::RationalView>
CaseMotionTimeJet::rotation_coefficient(std::uint64_t axis,
                                        p::JetOrder derivative) const& noexcept {
    return component(rotation_, axis, derivative);
}

core::Result<p::PchipJetLocation>
CaseMotionTimeJet::location(const SeriesJet& series) const noexcept {
    if (time().denominator.empty())
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::EvaluateStep});
    if (const auto* exact = std::get_if<p::ExactPchipJet>(&series))
        return exact->location();
    return p::PchipJetLocation{std::nullopt, side_, p::JetCandidate::Only};
}

core::Result<p::PchipJetLocation> CaseMotionTimeJet::translation_location() const noexcept {
    return location(translation_);
}

core::Result<p::PchipJetLocation> CaseMotionTimeJet::rotation_location() const noexcept {
    return location(rotation_);
}

core::Result<CaseMotionTimeSeries> prepare_case_motion_time_series(
    const io::SiCase& source, CaseMotionField motion, CaseMotionTimeWindow window,
    p::PchipInputStamp stamp, p::CurveWorkspace& workspace, const CaseMotionTimeLimits& limits,
    core::ResourcePlan& plan, const std::stop_token& stop, core::FailureDetails details) noexcept {
    using namespace motion_time_detail;
    return guarded<CaseMotionTimeSeries>(preparation, details, [&] {
        Call call{limits, plan, stop, preparation, details};
        return Builder::prepare(source, motion, window, stamp, workspace, limits, call);
    });
}

core::Result<CaseMotionTimeJet> evaluate_case_motion_time_series(
    const CaseMotionTimeSeries& source, const CaseMotionTimeSeries& current,
    CaseMotionTimeQuery query, p::CurveWorkspace& workspace, const CaseMotionTimeLimits& limits,
    core::ResourcePlan& plan, const std::stop_token& stop, core::FailureDetails details) noexcept {
    using namespace motion_time_detail;
    return guarded<CaseMotionTimeJet>(evaluation, details, [&] {
        Call call{limits, plan, stop, evaluation, details};
        return Builder::evaluate(source, current, query, workspace, call);
    });
}

} // namespace uvlas::coupling::detail
