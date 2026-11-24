#include "case_flow_clock.hpp"

#include <array>
#include <charconv>
#include <limits>

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
struct Builder {
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

    struct Node {
        io::detail::SiRationalWords time;
        const Grid::Interval& interval;
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

        const auto& arithmetic = limits.fields.conversion;
        auto time = [&]() -> io::detail::SiRationalWords {
            if (address.step == 0)
                return take(io::detail::copy_rational_words(interval.begin, memory, stop));
            if (address.step == interval.steps)
                return take(io::detail::copy_rational_words(interval.end, memory, stop));
            const auto k = exact_index(address.step, arithmetic, memory, details);
            const auto offset = take(interval.step.multiply(k, arithmetic, memory));
            const auto instant = take(interval.begin.add(offset, arithmetic, memory));
            return take(io::detail::copy_rational_words(instant, memory, stop));
        }();
        return {std::move(time), interval};
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
        auto [time, interval] =
            locate(field.grid_, current.grid_, address, limits, memory, stop, details);
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
        auto [time, interval] = locate(field.grid_, current.grid_, address, limits.clock_and_wake,
                                       memory, stop, details);
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

} // namespace uvlas::coupling::detail
