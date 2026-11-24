#include "case_flow_clock.hpp"

#include <array>
#include <charconv>
#include <limits>
#include <type_traits>

namespace uvlas::coupling::detail {
namespace a = aero::detail;
namespace p = a::pchip;
namespace s = io::schema;
namespace {
constexpr auto prepare = core::Operation::PrepareModel;
constexpr auto evaluate = core::Operation::EvaluateAero;
constexpr io::SchemaDimension seconds{{0, 0, 1, 0}};
using enum core::FailureCode;

template <class T>
T take(core::Result<T> result) {
    if (!result)
        throw result.error();
    return std::move(*result);
}

[[noreturn]] void fail(core::FailureCode code, core::Operation op, core::FailureDetails details) {
    throw core::Failure{code, op, details};
}

void check_stop(const std::stop_token& stop, core::Operation op, core::FailureDetails details) {
    if (stop.stop_requested())
        fail(Cancelled, op, details);
}

template <class T, class F>
core::Result<T> guarded(core::Operation op, core::FailureDetails details, F call) noexcept {
    return core::guard_entry(op, InvalidInput, details, [&]() -> core::Result<T> {
        try {
            return call();
        } catch (const core::Failure& error) {
            return std::unexpected(error);
        }
    });
}

p::RationalView rational_view(const io::detail::SiRationalWords& owned) noexcept {
    const auto words = owned.view();
    return {words.sign, words.numerator, words.denominator};
}

structure::ExactRational read_time(const io::SiCase& source, decltype(s::TimeInterval::start) field,
                                   const structure::ExactLimits& limits, core::ResourcePlan& plan,
                                   const std::stop_token& stop, core::FailureDetails details) {
    check_stop(stop, prepare, details);
    const auto quantity = take(source.quantity(field));
    if (quantity.recipe.dimensions.exponents != seconds.exponents)
        fail(InvalidInput, prepare, details);
    return take(io::materialize_rational(quantity.recipe, limits, plan));
}

template <class Arithmetic, class Value>
Value next_occurrence(Arithmetic& arithmetic, Value after, Value knot, Value period) {
    using enum p::RationalOp;
    const auto phase = take(arithmetic.binary(Subtract, after, knot));
    const auto cycles = take(arithmetic.binary(Divide, phase, period));
    const auto last = take(arithmetic.floor(cycles));
    const auto one = take(arithmetic.integer(1));
    const auto next = take(arithmetic.binary(Add, last, one));
    const auto offset = take(arithmetic.binary(Multiply, next, period));
    return take(arithmetic.binary(Add, knot, offset));
}

std::uint64_t step_count(const structure::ExactRational& count, std::uint64_t cap,
                         core::ResourcePlan& plan, const std::stop_token& stop,
                         core::FailureDetails details) {
    if (count.sign() != 1 || count.denominator_bits() != 1)
        fail(InvalidInput, prepare, details);
    if (count.numerator_bits() > std::numeric_limits<std::uint64_t>::digits)
        fail(ResourceExceeded, prepare, details);
    const auto words = take(io::detail::copy_rational_words(count, plan, stop));
    const auto value = words.view().numerator;
    auto steps = static_cast<std::uint64_t>(value[0]);
    if (value.size() == 2)
        steps |= static_cast<std::uint64_t>(value[1]) << 32U;
    if (steps > cap)
        fail(ResourceExceeded, prepare, details);
    return steps;
}

structure::ExactRational exact_index(std::uint64_t index, const structure::ExactLimits& limits,
                                     core::ResourcePlan& plan, core::FailureDetails details) {
    std::array<char, std::numeric_limits<std::uint64_t>::digits10 + 1> buffer{};
    const auto [end, error] = std::to_chars(buffer.data(), buffer.data() + buffer.size(), index);
    if (error != std::errc{})
        fail(ResourceExceeded, evaluate, details);
    return take(structure::ExactRational::from_decimal(
        {buffer.data(), static_cast<std::size_t>(end - buffer.data())}, limits, plan));
}
} // namespace

namespace case_clock_detail {
class TimeSeriesEventCheck {
public:
    TimeSeriesEventCheck(const io::SiCase& source, std::span<const Grid::Interval> intervals,
                         const structure::ExactLimits& conversion, p::RationalLimits arithmetic,
                         std::uint64_t work, core::ResourcePlan& memory,
                         const std::stop_token& stop, core::FailureDetails details) noexcept
        : source_(source),
          intervals_(intervals),
          conversion_(conversion),
          arithmetic_(arithmetic),
          work_(work),
          memory_(memory),
          stop_(stop),
          details_(details) {}

    void check_motions(const s::Model& model) {
        for (core::Index i = 0; std::cmp_less(i, model.aero_surfaces.size()); ++i) {
            charge(1);
            const auto declared = take(model.aero_surfaces.at(i)).record().attachment.record();
            if (const auto* prescribed = std::get_if<s::AeroAttachmentPrescribed>(&declared))
                check_motion(prescribed->motion.record());
        }
        for (core::Index i = 0; std::cmp_less(i, model.supports.size()); ++i) {
            charge(1);
            std::visit(
                [&](const auto& support) {
                    using Support = std::remove_cvref_t<decltype(support)>;
                    if constexpr (std::is_same_v<Support, s::SupportElastic>)
                        check_motion(support.reference.record());
                    else
                        check_motion(support.motion.record());
                },
                take(model.supports.at(i)).record());
        }
    }

    void check_model(const s::Model& model, const s::AnalysisTransient& analysis) {
        check_motions(model);
        check_loads(analysis.loads);
        check_member_loads(model.shells);
        check_member_loads(model.beams);
    }

private:
    template <class Members>
    void check_member_loads(const Members& members) {
        for (core::Index i = 0; std::cmp_less(i, members.size()); ++i) {
            charge(1);
            check_loads(take(members.at(i)).record().loads);
        }
    }

    template <class Loads>
    void check_loads(const Loads& loads) {
        for (core::Index i = 0; std::cmp_less(i, loads.size()); ++i) {
            charge(1);
            std::visit(
                [&](const auto& load) {
                    using Load = std::remove_cvref_t<decltype(load)>;
                    if constexpr (std::is_same_v<Load, s::LoadPointForce> ||
                                  std::is_same_v<Load, s::LoadCouple>)
                        check(load.series.record());
                    else if constexpr (std::is_same_v<Load, s::LoadGravity>)
                        check(load.acceleration.record());
                    else if constexpr (std::is_same_v<Load, s::LoadCurrentExpression>)
                        check_expression(load.expression.record());
                    else
                        check(load.multiplier.record());
                },
                take(loads.at(i)).record());
        }
    }

public:
    void check_events(decltype(s::Motion::events) events) {
        for (core::Index i = 0; std::cmp_less(i, events.size()); ++i) {
            charge(1);
            check_absolute(
                read_time(source_, take(events.at(i)), conversion_, memory_, stop_, details_));
        }
    }

    template <class Series>
    void check(const Series& series) {
        charge(1);
        if (const auto* expression = std::get_if<s::TimeSeriesExpression>(&series)) {
            check_expression(expression->expression.record());
            return;
        }
        const auto* pchip = std::get_if<s::TimeSeriesPchip>(&series);
        if (pchip == nullptr)
            return;
        const auto outside = pchip->outside.record();
        const auto* periodic = std::get_if<s::OutsidePeriodic>(&outside);
        std::optional<Period> period;
        if (periodic != nullptr) {
            auto value =
                read_time(source_, periodic->period, conversion_, memory_, stop_, details_);
            auto encoded = words(value);
            period.emplace(Period{std::move(value), std::move(encoded)});
        }
        // The last knot is the first knot one period later.
        const auto count = pchip->times.size() - (period ? 1U : 0U);
        for (core::Index i = 0; std::cmp_less(i, count); ++i) {
            charge(1);
            const auto knot =
                read_time(source_, take(pchip->times.at(i)), conversion_, memory_, stop_, details_);
            if (period) {
                const auto knot_words = words(knot);
                for (const auto& interval : intervals_)
                    check_periodic(interval, knot, knot_words, *period);
                continue;
            }
            check_absolute(knot);
        }
    }

private:
    void check_absolute(const structure::ExactRational& instant) {
        for (const auto& interval : intervals_) {
            charge(1);
            if (take(interval.begin.compare(instant, conversion_, memory_)) < 0 &&
                take(instant.compare(interval.end, conversion_, memory_)) < 0)
                fail(InvalidInput, prepare, details_);
        }
    }

    void check_motion(const s::Motion& motion) {
        check_events(motion.events);
        check(motion.translation.record());
        check(motion.rotation_vector.record());
    }

    void check_expression(const s::Expression& expression) {
        const auto domain = expression.domain.record();
        if (const auto* time = std::get_if<s::ExpressionDomainTime>(&domain))
            check_events(time->events);
        else if (const auto* current = std::get_if<s::ExpressionDomainCurrentPositionTime>(&domain))
            check_events(current->events);
    }

    struct Period {
        structure::ExactRational value;
        io::detail::SiRationalWords encoded;
    };

    void charge(std::uint64_t work) {
        check_stop(stop_, prepare, details_);
        if (work > work_)
            fail(ResourceExceeded, prepare, details_);
        work_ -= work;
    }

    static p::BitBound bound(const structure::ExactRational& value) noexcept {
        return {value.numerator_bits(), value.denominator_bits()};
    }

    [[nodiscard]] io::detail::SiRationalWords words(const structure::ExactRational& value) const {
        return take(io::detail::copy_rational_words(value, memory_, stop_));
    }

    void check_periodic(const Grid::Interval& interval, const structure::ExactRational& knot,
                        const io::detail::SiRationalWords& knot_words, const Period& period) {
        charge(1);
        auto limits = arithmetic_;
        limits.steps = std::min(limits.steps, work_);
        p::Counter counter{limits, details_};
        const auto begin_bound = take(counter.input(bound(interval.begin)));
        const auto end_bound = take(counter.input(bound(interval.end)));
        const auto knot_bound = take(counter.input(bound(knot)));
        const auto period_bound = take(counter.input(bound(period.value)));
        const auto next_bound = next_occurrence(counter, begin_bound, knot_bound, period_bound);
        if (auto compared = counter.compare(next_bound, end_bound); !compared)
            throw compared.error();
        const auto plan = take(counter.finish());
        charge(plan.reserved_steps());

        const auto begin_words = words(interval.begin);
        const auto end_words = words(interval.end);
        auto workspace = take(p::Workspace::create(plan, memory_, details_, stop_));
        auto session = take(workspace.begin(stop_, details_));
        const auto begin = take(session.input(rational_view(begin_words)));
        const auto end = take(session.input(rational_view(end_words)));
        const auto source_knot = take(session.input(rational_view(knot_words)));
        const auto source_period = take(session.input(rational_view(period.encoded)));
        const auto next = next_occurrence(session, begin, source_knot, source_period);
        if (take(session.compare(next, end)) < 0)
            fail(InvalidInput, prepare, details_);
    }

    const io::SiCase& source_;
    std::span<const Grid::Interval> intervals_;
    const structure::ExactLimits& conversion_;
    p::RationalLimits arithmetic_;
    std::uint64_t work_;
    core::ResourcePlan& memory_;
    const std::stop_token& stop_;
    core::FailureDetails details_;
};

struct Builder {
    static void check_model_events(const io::SiCase& source, const Grid& grid,
                                   const CaseFlowClockLimits& limits, core::ResourcePlan& memory,
                                   const std::stop_token& stop, core::FailureDetails details) {
        const auto root = take(take(source.source_case()).get().root());
        const auto transient = std::get<s::AnalysisTransient>(root.analysis.record());
        const auto arithmetic = limits.fields.flow.arithmetic;
        TimeSeriesEventCheck events{source,
                                    grid.intervals.values(),
                                    limits.fields.conversion,
                                    arithmetic,
                                    arithmetic.steps,
                                    memory,
                                    stop,
                                    details};
        events.check_model(root.model.record(), transient);
    }

    static Grid make_grid(const io::SiCase& source, p::PchipInputStamp stamp,
                          const CaseFlowClockLimits& limits, core::ResourcePlan& memory,
                          const std::stop_token& stop, core::FailureDetails details) {
        check_stop(stop, prepare, details);
        const auto& schema = take(source.source_case()).get();
        const auto origin = take(schema.identity());
        const auto root = take(schema.root());
        const auto analysis = root.analysis.record();
        const auto* transient = std::get_if<s::AnalysisTransient>(&analysis);
        if (transient == nullptr)
            fail(Unsupported, prepare, details);
        const auto size = transient->intervals.size();
        if (size == 0)
            fail(InvalidInput, prepare, details);
        if (size > limits.maximum_intervals || !std::in_range<core::Index>(size))
            fail(ResourceExceeded, prepare, details);
        auto intervals =
            take(memory.allocate<Grid::Interval>(static_cast<core::Index>(size), prepare));
        const auto& arithmetic = limits.fields.conversion;
        for (std::size_t i = 0; i < intervals.size(); ++i) {
            const auto input = take(transient->intervals.at(static_cast<core::Index>(i))).record();
            auto begin = read_time(source, input.start, arithmetic, memory, stop, details);
            auto end = read_time(source, input.end, arithmetic, memory, stop, details);
            auto step = read_time(source, input.step, arithmetic, memory, stop, details);
            const auto duration = take(end.subtract(begin, arithmetic, memory));
            const auto count = take(duration.divide(step, arithmetic, memory));
            const auto steps =
                step_count(count, limits.maximum_steps_per_interval, memory, stop, details);
            intervals.values()[i] = {std::move(begin), std::move(end), std::move(step), steps,
                                     input.restart.selected()};
        }
        check_stop(stop, prepare, details);
        return {origin, stamp, std::move(intervals)};
    }

    static CaseClockedFlows make(const io::SiCase& source, p::PchipInputStamp stamp,
                                 a::UniformFlowWorkspace& workspace,
                                 const CaseFlowClockLimits& limits, core::ResourcePlan& memory,
                                 const std::stop_token& stop, core::FailureDetails details) {
        auto grid = make_grid(source, stamp, limits, memory, stop, details);
        check_model_events(source, grid, limits, memory, stop, details);
        auto ambient = take(prepare_case_uniform_flow(source, stamp, {}, workspace, limits.fields,
                                                      memory, stop, details));
        auto wake = take(prepare_case_prescribed_wake(source, stamp, {}, workspace, limits.fields,
                                                      memory, stop, details));
        check_stop(stop, prepare, details);
        return {std::move(grid), std::move(ambient), std::move(wake)};
    }

    static CaseClockedKnownFlows make(const io::SiCase& source, p::PchipInputStamp stamp,
                                      a::UniformFlowWorkspace& workspace,
                                      const CaseClockedKnownFlowLimits& limits,
                                      core::ResourcePlan& memory, const std::stop_token& stop,
                                      core::FailureDetails details) {
        auto grid = make_grid(source, stamp, limits.clock_and_wake, memory, stop, details);
        check_model_events(source, grid, limits.clock_and_wake, memory, stop, details);
        auto ambient = take(prepare_case_known_flow(source, stamp, {}, workspace, limits.ambient,
                                                    memory, stop, details));
        auto wake = take(prepare_case_prescribed_wake(
            source, stamp, {}, workspace, limits.clock_and_wake.fields, memory, stop, details));
        check_stop(stop, prepare, details);
        return {std::move(grid), std::move(ambient), std::move(wake)};
    }

    static void bind(const Grid& source, const Grid& current, p::PchipInputStamp stamp,
                     const CaseFlowClockLimits& limits, const std::stop_token& stop,
                     core::FailureDetails details) {
        const auto prepared = source.intervals.values();
        const auto supplied = current.intervals.values();
        if (prepared.empty() || supplied.empty())
            fail(InvalidInput, evaluate, details);
        if (stamp != source.stamp || stamp != current.stamp || prepared.size() != supplied.size() ||
            source.source.contract_revision != current.source.contract_revision ||
            source.source.representation_revision != current.source.representation_revision)
            fail(StaleGeneration, evaluate, details);
        if (prepared.size() > limits.maximum_intervals)
            fail(ResourceExceeded, evaluate, details);
        const auto bits = limits.fields.conversion.maximum_value_bits;
        const auto fits = [bits](const structure::ExactRational& number) {
            return number.numerator_bits() <= bits && number.denominator_bits() <= bits;
        };
        for (std::size_t i = 0; i < prepared.size(); ++i) {
            check_stop(stop, evaluate, details);
            const auto& left = prepared[i];
            const auto& right = supplied[i];
            if (left.steps > limits.maximum_steps_per_interval ||
                right.steps > limits.maximum_steps_per_interval || !fits(left.begin) ||
                !fits(left.end) || !fits(left.step) || !fits(right.begin) || !fits(right.end) ||
                !fits(right.step))
                fail(ResourceExceeded, evaluate, details);
            if (left.restart != right.restart || left.steps != right.steps ||
                left.begin != right.begin || left.end != right.end || left.step != right.step)
                fail(StaleGeneration, evaluate, details);
        }
    }

    /// Both save and restore check bounds and qualifications with bind before copying.
    static Grid copy_grid(const Grid& source, const CaseFlowClockLimits& limits,
                          core::ResourcePlan& memory, const std::stop_token& stop,
                          core::FailureDetails details) {
        const auto intervals = source.intervals.values();
        const auto size = static_cast<core::Index>(intervals.size());
        auto copied = take(memory.allocate<Grid::Interval>(size, prepare));
        const auto& arithmetic = limits.fields.conversion;
        for (std::size_t i = 0; i < intervals.size(); ++i) {
            check_stop(stop, prepare, details);
            const auto& record = intervals[i];
            copied.values()[i] = {take(record.begin.clone(arithmetic, memory)),
                                  take(record.end.clone(arithmetic, memory)),
                                  take(record.step.clone(arithmetic, memory)), record.steps,
                                  record.restart};
        }
        return {source.source, source.stamp, std::move(copied)};
    }

    struct Node {
        const Grid::Interval& interval;
        std::optional<structure::ExactRational> interior;

        [[nodiscard]] const structure::ExactRational& time(std::uint64_t step) const noexcept {
            if (interior)
                return *interior;
            return step == 0 ? interval.begin : interval.end;
        }
    };

    static Node locate(const Grid& field, const Grid& current, CaseFlowClockQuery address,
                       const CaseFlowClockLimits& limits, core::ResourcePlan& memory,
                       const std::stop_token& stop, core::FailureDetails details) {
        check_stop(stop, evaluate, details);
        bind(field, current, address.stamp, limits, stop, details);
        if (address.interval >= field.intervals.size())
            fail(InvalidInput, evaluate, details);
        const auto& interval = field.intervals.values()[static_cast<std::size_t>(address.interval)];
        using enum p::RequestedSide;
        if (address.step > interval.steps ||
            (address.side != Left && address.side != Right && address.side != Automatic) ||
            (address.step == 0 && address.side != Right) ||
            (address.step == interval.steps && address.side != Left))
            fail(InvalidInput, evaluate, details);

        if (address.step == 0 || address.step == interval.steps)
            return {interval, std::nullopt};
        const auto& arithmetic = limits.fields.conversion;
        const auto k = exact_index(address.step, arithmetic, memory, details);
        const auto offset = take(interval.step.multiply(k, arithmetic, memory));
        return {interval, take(interval.begin.add(offset, arithmetic, memory))};
    }

    static p::CurveLimits execution_limits(p::CurveLimits limits, const p::CurveAdmission& forecast,
                                           core::FailureDetails details) {
        // The forecast does not report its actual scan work separately.
        // Reserve its full bound as well as the subsequent evaluation.
        if (forecast.reserved_steps > limits.arithmetic.steps / 2)
            fail(ResourceExceeded, evaluate, details);
        limits.arithmetic.steps -= forecast.reserved_steps;
        return limits;
    }

    static CaseClockedFlowSample
    sample(const CaseClockedFlows& field, const CaseClockedFlows& current,
           CaseFlowClockQuery address, const CaseFlowProjectionRequest& ambient,
           const CaseFlowProjectionRequest& wake, a::UniformFlowWorkspace& workspace,
           const CaseFlowClockLimits& limits, core::ResourcePlan& memory,
           const std::stop_token& stop, core::FailureDetails details) {
        const auto node =
            locate(field.grid_, current.grid_, address, limits, memory, stop, details);
        const auto& interval = node.interval;
        auto time = take(io::detail::copy_rational_words(node.time(address.step), memory, stop));
        const CaseUniformFlowQuery query{address.stamp, rational_view(time), address.side};
        const auto ambient_plan = take(forecast_case_uniform_flow(
            field.ambient_, current.ambient_, query, limits.fields.flow, stop, details));
        const auto wake_plan = take(forecast_case_prescribed_wake(
            field.wake_, current.wake_, query, limits.fields.flow, stop, details));
        const auto ambient_limits = execution_limits(limits.fields.flow, ambient_plan, details);
        const auto wake_limits = execution_limits(limits.fields.flow, wake_plan, details);
        auto step = take(io::detail::copy_rational_words(interval.step, memory, stop));
        const auto remaining = memory.remaining_bytes();
        if (ambient_plan.payload_bytes > remaining ||
            wake_plan.payload_bytes > remaining - ambient_plan.payload_bytes)
            fail(ResourceExceeded, evaluate, details);

        auto exact_ambient =
            take(evaluate_case_uniform_flow(field.ambient_, current.ambient_, query, workspace,
                                            ambient_limits, memory, stop, details));
        auto exact_wake = take(evaluate_case_prescribed_wake(
            field.wake_, current.wake_, query, workspace, wake_limits, memory, stop, details));
        auto projected_ambient =
            take(project_case_uniform_flow(std::move(exact_ambient), ambient, stop, details));
        auto projected_wake =
            take(project_case_prescribed_wake(std::move(exact_wake), wake, stop, details));
        check_stop(stop, evaluate, details);
        return {field.grid_.source,       address,         interval.restart,
                std::move(time),          std::move(step), std::move(projected_ambient),
                std::move(projected_wake)};
    }

    static CaseClockedKnownFlowSample
    sample(const CaseClockedKnownFlows& field, const CaseClockedKnownFlows& current,
           CaseFlowClockQuery address, const std::array<p::RationalView, 3>& position,
           const CaseKnownFlowProjection& ambient, const CaseFlowProjectionRequest& wake,
           a::UniformFlowWorkspace& workspace, const CaseClockedKnownFlowLimits& limits,
           core::ResourcePlan& memory, const std::stop_token& stop, core::FailureDetails details) {
        const auto node = locate(field.grid_, current.grid_, address, limits.clock_and_wake, memory,
                                 stop, details);
        const auto& interval = node.interval;
        const auto time =
            take(io::detail::copy_rational_words(node.time(address.step), memory, stop));
        const CaseUniformFlowQuery query{address.stamp, rational_view(time), address.side};
        const auto& wake_limits = limits.clock_and_wake.fields.flow;
        auto step = take(io::detail::copy_rational_words(interval.step, memory, stop));
        auto physical = take(evaluate_case_known_flow(field.ambient_, current.ambient_,
                                                      {query, position}, ambient, workspace,
                                                      limits.ambient, memory, stop, details));
        auto exact_wake = take(evaluate_case_prescribed_wake(
            field.wake_, current.wake_, query, workspace, wake_limits, memory, stop, details));
        auto projected_wake =
            take(project_case_prescribed_wake(std::move(exact_wake), wake, stop, details));
        check_stop(stop, evaluate, details);
        return {field.grid_.source,  address,
                interval.restart,    std::move(step),
                std::move(physical), std::move(projected_wake)};
    }

    static CaseClockedPrescribedInputs
    make(const io::SiCase& source, const s::AeroAttachmentPrescribed& attachment,
         p::PchipInputStamp stamp, a::UniformFlowWorkspace& flow_workspace,
         p::CurveWorkspace& motion_workspace, const CaseClockedPrescribedLimits& limits,
         core::ResourcePlan& memory, const std::stop_token& stop, core::FailureDetails details) {
        const auto& clock = limits.flow.clock_and_wake;
        auto grid = make_grid(source, stamp, clock, memory, stop, details);
        const auto intervals = grid.intervals.values();
        auto flow = take(prepare_case_known_flow(source, stamp, {}, flow_workspace,
                                                 limits.flow.ambient, memory, stop, details));
        auto moving = take(prepare_prescribed_aero_attachment(
            source, attachment, {intervals.front().begin, intervals.back().end}, stamp,
            motion_workspace, limits.motion, memory, stop, details));
        const auto root = take(take(source.source_case()).get().root());
        // Grid and ambient preparation already qualified these variants.
        const auto ambient = std::get<s::FlowKnown>(root.flow.record());
        const auto gust = ambient.gust.record();
        const auto transient = std::get<s::AnalysisTransient>(root.analysis.record());
        const auto wake = transient.wake.record();
        const auto* active = std::get_if<s::WakeActive>(&wake);
        if (active == nullptr)
            fail(Unsupported, prepare, details);
        const auto convection = active->convection.record();
        const auto* prescribed = std::get_if<s::ConvectionPrescribed>(&convection);
        TimeSeriesEventCheck events{source,
                                    intervals,
                                    clock.fields.conversion,
                                    clock.fields.flow.arithmetic,
                                    limits.maximum_event_work,
                                    memory,
                                    stop,
                                    details};
        events.check_model(root.model.record(), transient);
        events.check_events(ambient.events);
        events.check(ambient.uniform.record());
        if (const auto* sine = std::get_if<s::GustSine>(&gust))
            events.check(sine->amplitude.record());
        std::optional<CasePrescribedWakeReady> transport;
        if (prescribed != nullptr) {
            transport.emplace(take(prepare_case_prescribed_wake(
                source, stamp, {}, flow_workspace, clock.fields, memory, stop, details)));
            events.check(prescribed->velocity.record());
        }
        check_stop(stop, prepare, details);
        return {std::move(grid), std::move(flow), std::move(transport), std::move(moving)};
    }

    static CaseClockedPrescribedSnapshot
    save(const CaseClockedPrescribedInputs& source, a::UniformFlowWorkspace& field_workspace,
         p::CurveWorkspace& motion_workspace, const CaseClockedPrescribedLimits& limits,
         core::ResourcePlan& memory, const std::stop_token& stop, core::FailureDetails details) {
        check_stop(stop, prepare, details);
        const auto& clock = limits.flow.clock_and_wake;
        bind(source.grid_, source.grid_, source.grid_.stamp, clock, stop, details);
        auto grid = copy_grid(source.grid_, clock, memory, stop, details);
        auto ambient = take(snapshot_case_known_flow(source.ambient_, field_workspace,
                                                     limits.flow.ambient, memory, stop, details));
        std::optional<CasePrescribedWakeSnapshot> transport;
        if (const auto& prescribed = source.prescribed_)
            transport.emplace(take(snapshot_case_prescribed_wake(
                *prescribed, field_workspace, clock.fields.flow, memory, stop, details)));
        auto motion = take(snapshot_case_motion(source.attachment_.motion_, motion_workspace,
                                                limits.motion, memory, stop, details));
        check_stop(stop, prepare, details);
        return {std::move(grid), std::move(ambient), std::move(transport), std::move(motion)};
    }

    static CaseClockedPrescribedInputs
    restore(const CaseClockedPrescribedSnapshot& saved, const CaseClockedPrescribedInputs& current,
            a::UniformFlowWorkspace& field_workspace, p::CurveWorkspace& motion_workspace,
            const CaseClockedPrescribedLimits& limits, core::ResourcePlan& memory,
            const std::stop_token& stop, core::FailureDetails details) {
        check_stop(stop, prepare, details);
        const auto& clock = limits.flow.clock_and_wake;
        bind(saved.grid_, current.grid_, current.grid_.stamp, clock, stop, details);
        const auto& saved_wake = saved.prescribed_;
        const auto& live_wake = current.prescribed_;
        if (saved_wake.has_value() != live_wake.has_value())
            fail(StaleGeneration, prepare, details);
        auto ambient =
            take(restore_case_known_flow(saved.ambient_, current.ambient_, field_workspace,
                                         limits.flow.ambient, memory, stop, details));
        std::optional<CasePrescribedWakeReady> transport;
        if (saved_wake && live_wake)
            transport.emplace(
                take(restore_case_prescribed_wake(*saved_wake, *live_wake, field_workspace,
                                                  clock.fields.flow, memory, stop, details)));
        auto moving =
            take(restore_case_motion(take(saved.motion_.view()), current.attachment_.motion_,
                                     motion_workspace, limits.motion, memory, stop, details));
        auto grid = copy_grid(current.grid_, clock, memory, stop, details);
        check_stop(stop, prepare, details);
        return {std::move(grid), std::move(ambient), std::move(transport),
                PrescribedAeroAttachment{std::move(moving)}};
    }

    static CaseClockedPrescribedSample
    sample(const CaseClockedPrescribedInputs& field, const CaseClockedPrescribedInputs& current,
           CaseFlowClockQuery address, const std::array<p::RationalView, 3>& position,
           const CaseKnownFlowProjection& ambient,
           const std::optional<CaseFlowProjectionRequest>& wake, p::JetOrder motion_order,
           const CaseMotionProjection& motion_projection, a::UniformFlowWorkspace& flow_workspace,
           p::CurveWorkspace& motion_workspace, const CaseClockedPrescribedLimits& limits,
           core::ResourcePlan& memory, const std::stop_token& stop, core::FailureDetails details) {
        const auto node = locate(field.grid_, current.grid_, address, limits.flow.clock_and_wake,
                                 memory, stop, details);
        if (field.prescribed_.has_value() != current.prescribed_.has_value())
            fail(StaleGeneration, evaluate, details);
        if (field.prescribed_.has_value() != wake.has_value())
            fail(InvalidInput, evaluate, details);
        const auto& instant = node.time(address.step);
        auto motion = take(evaluate_case_motion(
            field.attachment_.motion_, current.attachment_.motion_,
            {address.stamp, instant, motion_order, address.side}, motion_projection,
            motion_workspace, limits.motion, memory, stop, details));
        const auto words = take(io::detail::copy_rational_words(instant, memory, stop));
        const CaseUniformFlowQuery time{address.stamp, rational_view(words), address.side};
        auto flow = take(evaluate_case_known_flow(field.ambient_, current.ambient_,
                                                  {time, position}, ambient, flow_workspace,
                                                  limits.flow.ambient, memory, stop, details));
        auto step = take(io::detail::copy_rational_words(node.interval.step, memory, stop));
        std::optional<ProjectedCasePrescribedWake> transport;
        if (wake) {
            auto exact = take(evaluate_case_prescribed_wake(
                *field.prescribed_, *current.prescribed_, time, flow_workspace,
                limits.flow.clock_and_wake.fields.flow, memory, stop, details));
            transport.emplace(
                take(project_case_prescribed_wake(std::move(exact), *wake, stop, details)));
        }
        check_stop(stop, evaluate, details);
        return {field.grid_.source, address,         node.interval.restart,
                std::move(step),    std::move(flow), std::move(transport),
                std::move(motion)};
    }

    static PrescribedAeroTrial
    attachment(const CaseClockedPrescribedInputs& field, const CaseClockedPrescribedInputs& current,
               CaseFlowClockQuery address, const PrescribedAeroReference& reference,
               const PrescribedAeroProjection& projection, p::CurveWorkspace& motion,
               a::GeometryWorkspace& geometry, a::GlobalWorkspace& contact,
               const CaseFlowClockLimits& clock_limits, const PrescribedAeroLimits& limits,
               core::ResourcePlan& memory, const std::stop_token& stop,
               core::FailureDetails details) {
        const auto node =
            locate(field.grid_, current.grid_, address, clock_limits, memory, stop, details);
        if (field.prescribed_.has_value() != current.prescribed_.has_value())
            fail(StaleGeneration, evaluate, details);
        return take(evaluate_prescribed_aero_attachment(
            field.attachment_, current.attachment_, reference,
            {address.stamp, node.time(address.step), p::JetOrder::First, address.side}, projection,
            motion, geometry, contact, limits, memory, stop, details));
    }
};
} // namespace case_clock_detail

CaseClockedFlows::CaseClockedFlows(case_clock_detail::Grid&& grid, CaseUniformFlowReady&& ambient,
                                   CasePrescribedWakeReady&& wake) noexcept
    : grid_(std::move(grid)), ambient_(std::move(ambient)), wake_(std::move(wake)) {}

core::Result<std::uint64_t>
case_clock_detail::Grid::interval_steps(std::uint64_t interval) const noexcept {
    if (interval >= intervals.size())
        return std::unexpected(core::Failure{InvalidInput, evaluate});
    return intervals.values()[static_cast<std::size_t>(interval)].steps;
}

core::Result<std::uint64_t>
CaseClockedFlows::interval_steps(std::uint64_t interval) const noexcept {
    return grid_.interval_steps(interval);
}

CaseClockedFlowSample::CaseClockedFlowSample(io::SchemaIdentity source, CaseFlowClockQuery address,
                                             CaseClockRestart restart,
                                             io::detail::SiRationalWords&& instant,
                                             io::detail::SiRationalWords&& interval_step,
                                             ProjectedCaseUniformFlow&& ambient,
                                             ProjectedCasePrescribedWake&& wake) noexcept
    : source_(source),
      address_(address),
      restart_(restart),
      instant_(std::move(instant)),
      interval_step_(std::move(interval_step)),
      ambient_(std::move(ambient)),
      wake_(std::move(wake)) {}

core::Result<p::RationalView> CaseClockedFlowSample::time() const& noexcept {
    if (!ambient_.source().jet.active())
        return std::unexpected(core::Failure{InvalidInput, evaluate});
    return rational_view(instant_);
}

core::Result<p::RationalView> CaseClockedFlowSample::step_seconds() const& noexcept {
    if (!ambient_.source().jet.active())
        return std::unexpected(core::Failure{InvalidInput, evaluate});
    return rational_view(interval_step_);
}

core::Result<CaseClockedFlows>
prepare_case_clocked_flows(const io::SiCase& source, p::PchipInputStamp stamp,
                           a::UniformFlowWorkspace& workspace, const CaseFlowClockLimits& limits,
                           core::ResourcePlan& memory, const std::stop_token& stop,
                           core::FailureDetails details) noexcept {
    return guarded<CaseClockedFlows>(prepare, details, [&] {
        return case_clock_detail::Builder::make(source, stamp, workspace, limits, memory, stop,
                                                details);
    });
}

core::Result<CaseClockedFlowSample>
evaluate_case_clocked_flows(const CaseClockedFlows& field, const CaseClockedFlows& current,
                            CaseFlowClockQuery address, const CaseFlowProjectionRequest& ambient,
                            const CaseFlowProjectionRequest& wake,
                            a::UniformFlowWorkspace& workspace, const CaseFlowClockLimits& limits,
                            core::ResourcePlan& memory, const std::stop_token& stop,
                            core::FailureDetails details) noexcept {
    return guarded<CaseClockedFlowSample>(evaluate, details, [&] {
        return case_clock_detail::Builder::sample(field, current, address, ambient, wake, workspace,
                                                  limits, memory, stop, details);
    });
}

core::Result<p::RationalView> CaseClockedKnownFlowSample::step_seconds() const& noexcept {
    if (!ambient_.time())
        return std::unexpected(core::Failure{InvalidInput, evaluate});
    return rational_view(step_);
}

core::Result<CaseClockedKnownFlows> prepare_case_clocked_known_flows(
    const io::SiCase& source, p::PchipInputStamp stamp, a::UniformFlowWorkspace& workspace,
    const CaseClockedKnownFlowLimits& limits, core::ResourcePlan& memory,
    const std::stop_token& stop, core::FailureDetails details) noexcept {
    return guarded<CaseClockedKnownFlows>(prepare, details, [&] {
        return case_clock_detail::Builder::make(source, stamp, workspace, limits, memory, stop,
                                                details);
    });
}

core::Result<CaseClockedKnownFlowSample> evaluate_case_clocked_known_flows(
    const CaseClockedKnownFlows& field, const CaseClockedKnownFlows& current,
    CaseFlowClockQuery address, const std::array<p::RationalView, 3>& position,
    const CaseKnownFlowProjection& ambient, const CaseFlowProjectionRequest& wake,
    a::UniformFlowWorkspace& workspace, const CaseClockedKnownFlowLimits& limits,
    core::ResourcePlan& memory, const std::stop_token& stop,
    core::FailureDetails details) noexcept {
    return guarded<CaseClockedKnownFlowSample>(evaluate, details, [&] {
        return case_clock_detail::Builder::sample(field, current, address, position, ambient, wake,
                                                  workspace, limits, memory, stop, details);
    });
}

core::Result<CaseClockedPrescribedInputs> prepare_case_clocked_prescribed_inputs(
    const io::SiCase& source, const s::AeroAttachmentPrescribed& attachment,
    p::PchipInputStamp stamp, a::UniformFlowWorkspace& flow_workspace,
    p::CurveWorkspace& motion_workspace, const CaseClockedPrescribedLimits& limits,
    core::ResourcePlan& memory, const std::stop_token& stop,
    core::FailureDetails details) noexcept {
    return guarded<CaseClockedPrescribedInputs>(prepare, details, [&] {
        return case_clock_detail::Builder::make(source, attachment, stamp, flow_workspace,
                                                motion_workspace, limits, memory, stop, details);
    });
}

core::Result<CaseClockedPrescribedSnapshot> snapshot_case_clocked_prescribed_inputs(
    const CaseClockedPrescribedInputs& source, a::UniformFlowWorkspace& field_workspace,
    p::CurveWorkspace& motion_workspace, const CaseClockedPrescribedLimits& limits,
    core::ResourcePlan& memory, const std::stop_token& stop,
    core::FailureDetails details) noexcept {
    return guarded<CaseClockedPrescribedSnapshot>(prepare, details, [&] {
        return case_clock_detail::Builder::save(source, field_workspace, motion_workspace, limits,
                                                memory, stop, details);
    });
}

core::Result<CaseClockedPrescribedInputs> restore_case_clocked_prescribed_inputs(
    const CaseClockedPrescribedSnapshot& saved, const CaseClockedPrescribedInputs& current,
    a::UniformFlowWorkspace& field_workspace, p::CurveWorkspace& motion_workspace,
    const CaseClockedPrescribedLimits& limits, core::ResourcePlan& memory,
    const std::stop_token& stop, core::FailureDetails details) noexcept {
    return guarded<CaseClockedPrescribedInputs>(prepare, details, [&] {
        return case_clock_detail::Builder::restore(saved, current, field_workspace,
                                                   motion_workspace, limits, memory, stop, details);
    });
}

core::Result<CaseClockedPrescribedSample> evaluate_case_clocked_prescribed_inputs(
    const CaseClockedPrescribedInputs& field, const CaseClockedPrescribedInputs& current,
    CaseFlowClockQuery address, const std::array<p::RationalView, 3>& position,
    const CaseKnownFlowProjection& ambient, const std::optional<CaseFlowProjectionRequest>& wake,
    p::JetOrder motion_order, const CaseMotionProjection& motion,
    a::UniformFlowWorkspace& flow_workspace, p::CurveWorkspace& motion_workspace,
    const CaseClockedPrescribedLimits& limits, core::ResourcePlan& memory,
    const std::stop_token& stop, core::FailureDetails details) noexcept {
    return guarded<CaseClockedPrescribedSample>(evaluate, details, [&] {
        return case_clock_detail::Builder::sample(field, current, address, position, ambient, wake,
                                                  motion_order, motion, flow_workspace,
                                                  motion_workspace, limits, memory, stop, details);
    });
}

core::Result<PrescribedAeroTrial> evaluate_case_clocked_prescribed_attachment(
    const CaseClockedPrescribedInputs& field, const CaseClockedPrescribedInputs& current,
    CaseFlowClockQuery address, const PrescribedAeroReference& reference,
    const PrescribedAeroProjection& projection, p::CurveWorkspace& motion,
    a::GeometryWorkspace& geometry, a::GlobalWorkspace& contact,
    const CaseFlowClockLimits& clock_limits, const PrescribedAeroLimits& limits,
    core::ResourcePlan& memory, const std::stop_token& stop,
    core::FailureDetails details) noexcept {
    return guarded<PrescribedAeroTrial>(evaluate, details, [&] {
        return case_clock_detail::Builder::attachment(field, current, address, reference,
                                                      projection, motion, geometry, contact,
                                                      clock_limits, limits, memory, stop, details);
    });
}

} // namespace uvlas::coupling::detail
