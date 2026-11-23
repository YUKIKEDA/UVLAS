#include "case_uniform_flow.hpp"

#include "io/si_exact_words.hpp"

#include <limits>
#include <type_traits>

namespace uvlas::coupling::detail {
namespace a = aero::detail;
namespace p = a::pchip;
namespace s = io::schema;
using p::RationalView;

namespace {
constexpr auto preparation = core::Operation::PrepareModel;
constexpr auto evaluation = core::Operation::EvaluateAero;
constexpr auto index_max = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
constexpr std::size_t xyz = 3;
constexpr io::SchemaDimension seconds{{0, 0, 1, 0}};
constexpr io::SchemaDimension velocity{{1, 0, -1, 0}};
using enum core::FailureCode;

template <class T>
T take(core::Result<T> result) {
    if (!result)
        throw result.error();
    return std::move(*result);
}

[[noreturn]] void fail(core::FailureCode code, core::Operation operation,
                       core::FailureDetails details) {
    throw core::Failure{code, operation, details};
}

template <class T, class Call>
core::Result<T> guarded(core::Operation operation, core::FailureDetails details,
                        Call call) noexcept {
    return core::guard_entry(operation, InvalidInput, details, [&]() -> core::Result<T> {
        try {
            return call();
        } catch (const core::Failure& error) {
            return std::unexpected(error);
        }
    });
}

std::uint64_t add_count(std::uint64_t left, std::uint64_t right, core::FailureDetails details) {
    if (left > index_max || right > index_max - left)
        fail(ResourceExceeded, preparation, details);
    return left + right;
}

std::uint64_t components(std::uint64_t rows, core::FailureDetails details) {
    if (rows > index_max / xyz)
        fail(ResourceExceeded, preparation, details);
    return rows * xyz;
}

struct TimeWindowFields {
    decltype(s::TimeInterval::start) begin, end;
};

using TimeWindow = std::variant<TimeWindowFields, CaseUniformTimeWindow>;

TimeWindow window_for(const s::Case& root, std::optional<CaseUniformTimeWindow> supplied,
                      core::FailureDetails details) {
    const auto analysis = root.analysis.record();
    if (const auto* transient = std::get_if<s::AnalysisTransient>(&analysis)) {
        const auto count = transient->intervals.size();
        if (count == 0 || count > index_max)
            fail(InvalidInput, preparation, details);
        const auto begin = take(transient->intervals.at(0)).record().start;
        const auto end =
            take(transient->intervals.at(static_cast<core::Index>(count - 1))).record().end;
        return TimeWindowFields{begin, end};
    }
    if (!supplied)
        fail(Unsupported, preparation, details);
    return *supplied;
}

/// Holds every exact SI recipe result until the velocity provider copies or retains it.
class SiSeries {
public:
    SiSeries(const io::SiCase& source, std::uint64_t count, const CaseUniformFlowLimits& limits,
             core::ResourcePlan& memory, const std::stop_token& stop, core::FailureDetails details)
        : source_(source),
          limits_(limits.conversion),
          memory_(memory),
          stop_(stop),
          details_(details) {
        constexpr auto element_bytes = sizeof(io::detail::SiRationalWords) + sizeof(RationalView);
        check_stop();
        if (count > limits.maximum_exact_numbers || count > index_max ||
            count > memory.remaining_bytes() / element_bytes)
            fail(ResourceExceeded, preparation, details_);
        values_ = take(memory_.allocate<io::detail::SiRationalWords>(
            static_cast<core::Index>(count), preparation));
        check_stop();
        records_ =
            take(memory_.allocate<RationalView>(static_cast<core::Index>(count), preparation));
    }

    template <class Field>
    RationalView read(Field number, const io::SchemaDimension& expected) {
        check_stop();
        if (cursor_ == records_.values().size())
            fail(ResourceExceeded, preparation, details_);
        const auto quantity = take(source_.quantity(number));
        if (quantity.recipe.dimensions.exponents != expected.exponents)
            fail(InvalidInput, preparation, details_);
        auto& owned = values_.values()[cursor_];
        owned =
            take(io::detail::materialize_rational_words(quantity.recipe, limits_, memory_, stop_));
        const auto words = owned.view();
        const RationalView result{words.sign, words.numerator, words.denominator};
        records_.values()[cursor_++] = result;
        return result;
    }

    template <class Array>
    std::span<const RationalView> sequence(Array array, const io::SchemaDimension& expected) {
        const auto start = cursor_;
        for (std::uint64_t i = 0; i < array.size(); ++i)
            (void)read(take(array.at(static_cast<core::Index>(i))), expected);
        return records_.values().subspan(start, cursor_ - start);
    }

    template <class Array>
    std::array<RationalView, xyz> row(Array array) {
        if (array.size() != xyz)
            fail(InvalidInput, preparation, details_);
        std::array<RationalView, xyz> result{};
        for (std::size_t i = 0; i < xyz; ++i)
            result[i] = read(take(array.at(static_cast<core::Index>(i))), velocity);
        return result;
    }

    template <class Array>
    std::span<const RationalView> grid(Array array) {
        const auto start = cursor_;
        for (std::uint64_t i = 0; i < array.size(); ++i)
            (void)row(take(array.at(static_cast<core::Index>(i))));
        return records_.values().subspan(start, cursor_ - start);
    }

private:
    void check_stop() const {
        if (stop_.stop_requested())
            fail(Cancelled, preparation, details_);
    }

    const io::SiCase& source_;
    const structure::ExactLimits& limits_;
    core::ResourcePlan& memory_;
    const std::stop_token& stop_;
    core::FailureDetails details_;
    core::PlannedBuffer<io::detail::SiRationalWords> values_;
    core::PlannedBuffer<RationalView> records_;
    std::size_t cursor_ = 0;
};

using BoundaryPolicy = decltype(p::BoundaryRequest::policy);
using Outside = decltype(s::TimeSeriesPchip::outside.record());

BoundaryPolicy outside_policy(const Outside& outside, SiSeries& values) {
    if (const auto* period = std::get_if<s::OutsidePeriodic>(&outside))
        return p::PeriodicBoundary{values.read(period->period, seconds)};
    const auto tag = std::get<s::OutsideRejectHold>(outside).kind.selected();
    if (tag == decltype(tag)::ItemHold)
        return p::HoldBoundary{};
    return p::RejectBoundary{};
}

CaseUniformTimeWindow canonical_window(const TimeWindow& source, SiSeries& values) {
    if (const auto* typed = std::get_if<TimeWindowFields>(&source))
        return {values.read(typed->begin, seconds), values.read(typed->end, seconds)};
    return std::get<CaseUniformTimeWindow>(source);
}

struct ImportContext {
    const io::SiCase& source;
    const s::Case& root;
    const s::FlowKnown& flow;
    p::PchipInputStamp stamp;
    std::optional<CaseUniformTimeWindow> window;
    a::UniformFlowWorkspace& workspace;
    const CaseUniformFlowLimits& limits;
    core::ResourcePlan& memory;
    const std::stop_token& stop;
    core::FailureDetails details;

    a::UniformFlowReady prepare(const a::UniformFlowInput& request) const {
        return take(
            a::prepare_uniform_flow(request, workspace, limits.flow, memory, stop, details));
    }

    template <class Series>
    a::UniformFlowReady operator()(const Series& series) const {
        if constexpr (std::is_same_v<Series, s::TimeSeriesConstant>) {
            SiSeries values{
                source, add_count(flow.events.size(), xyz, details), limits, memory, stop, details};
            const auto events = values.sequence(flow.events, seconds);
            return prepare({stamp, events, a::ConstantUniformFlow{values.row(series.value)}});
        } else if constexpr (std::is_same_v<Series, s::TimeSeriesPolynomial>) {
            const auto count =
                add_count(4,
                          add_count(flow.events.size(),
                                    components(series.coefficients.size(), details), details),
                          details);
            SiSeries values{source, count, limits, memory, stop, details};
            const auto events = values.sequence(flow.events, seconds);
            return prepare(
                {stamp, events,
                 a::PolynomialUniformFlow{values.read(series.origin, seconds),
                                          values.read(series.scale, seconds),
                                          values.read(take(series.domain.at(0)), seconds),
                                          values.read(take(series.domain.at(1)), seconds),
                                          values.grid(series.coefficients)}});
        } else if constexpr (std::is_same_v<Series, s::TimeSeriesPchip>) {
            const auto requested = window_for(root, window, details);
            const auto outside = series.outside.record();
            auto count = add_count(flow.events.size(), series.times.size(), details);
            count = add_count(count, components(series.values.size(), details), details);
            if (std::holds_alternative<TimeWindowFields>(requested))
                count = add_count(count, 2, details);
            if (std::holds_alternative<s::OutsidePeriodic>(outside))
                count = add_count(count, 1, details);
            SiSeries values{source, count, limits, memory, stop, details};
            const auto events = values.sequence(flow.events, seconds);
            const auto times = values.sequence(series.times, seconds);
            const auto samples = values.grid(series.values);
            const auto interval = canonical_window(requested, values);
            const auto policy = outside_policy(outside, values);
            return prepare({stamp, events,
                            a::PchipUniformFlow{{times, samples, xyz},
                                                {interval.begin, interval.end, policy}}});
        } else {
            fail(Unsupported, preparation, details);
        }
    }
};

template <class Call>
auto with_current(const CaseUniformFlowReady& field, const CaseUniformFlowReady& current,
                  CaseUniformFlowQuery request, p::CurveLimits limits, const std::stop_token& stop,
                  core::FailureDetails details, Call call) {
    const auto source = take(field.uniform().input());
    const auto expected = take(current.uniform().input());
    if (request.stamp != source.stamp || request.stamp != expected.stamp ||
        source.series.index() != expected.series.index())
        fail(StaleGeneration, evaluation, details);
    a::UniformFlowQuery query{request.stamp, request.time, request.side};
    if (std::holds_alternative<a::PchipUniformFlow>(expected.series)) {
        const auto owned = take(current.uniform().retain_pchip_input());
        query.expected_pchip.emplace(owned);
        return call(query, limits, std::uint64_t{0});
    }
    const auto checked = take(a::match_uniform_formula_input(
        field.uniform(), current.uniform(), limits.arithmetic.steps, stop, details));
    limits.arithmetic.steps -= checked;
    return call(query, limits, checked);
}

} // namespace

namespace case_flow_detail {
struct Builder {
    static CaseUniformFlowReady prepare(const io::SiCase& source, p::PchipInputStamp stamp,
                                        std::optional<CaseUniformTimeWindow> window,
                                        a::UniformFlowWorkspace& workspace,
                                        const CaseUniformFlowLimits& limits,
                                        core::ResourcePlan& memory, const std::stop_token& stop,
                                        core::FailureDetails details) {
        if (stop.stop_requested())
            fail(Cancelled, preparation, details);
        const auto schema = take(source.source_case());
        const auto root = take(schema.get().root());
        const auto flow = root.flow.record();
        const auto* known = std::get_if<s::FlowKnown>(&flow);
        if (!known || !std::holds_alternative<s::GustNone>(known->gust.record()))
            fail(Unsupported, preparation, details);
        const ImportContext session{source,    root,   *known, stamp, window,
                                    workspace, limits, memory, stop,  details};
        return CaseUniformFlowReady{std::visit(session, known->uniform.record())};
    }

    static CaseUniformFlowReady restore(const a::UniformFlowSnapshotView& saved,
                                        const CaseUniformFlowReady& current,
                                        a::UniformFlowWorkspace& workspace, p::CurveLimits limits,
                                        core::ResourcePlan& memory, const std::stop_token& stop,
                                        core::FailureDetails details) {
        return CaseUniformFlowReady{take(a::restore_uniform_flow(
            saved, current.uniform(), workspace, limits, memory, stop, details))};
    }
};
} // namespace case_flow_detail

CaseUniformFlowReady::CaseUniformFlowReady(a::UniformFlowReady flow) noexcept
    : flow_(std::move(flow)) {}

const a::UniformFlowReady& CaseUniformFlowReady::uniform() const& noexcept {
    return flow_;
}

core::Result<CaseUniformFlowReady>
prepare_case_uniform_flow(const io::SiCase& source, p::PchipInputStamp stamp,
                          std::optional<CaseUniformTimeWindow> window,
                          a::UniformFlowWorkspace& workspace, const CaseUniformFlowLimits& limits,
                          core::ResourcePlan& memory, const std::stop_token& stop,
                          core::FailureDetails details) noexcept {
    return guarded<CaseUniformFlowReady>(preparation, details, [&] {
        return case_flow_detail::Builder::prepare(source, stamp, window, workspace, limits, memory,
                                                  stop, details);
    });
}

core::Result<p::CurveAdmission>
forecast_case_uniform_flow(const CaseUniformFlowReady& field, const CaseUniformFlowReady& current,
                           CaseUniformFlowQuery query, p::CurveLimits limits,
                           const std::stop_token& stop, core::FailureDetails details) noexcept {
    return guarded<p::CurveAdmission>(evaluation, details, [&] {
        return with_current(
            field, current, query, limits, stop, details,
            [&](a::UniformFlowQuery request, p::CurveLimits remaining, std::uint64_t compared) {
                auto report = take(a::forecast_uniform_flow_evaluation(field.uniform(), request,
                                                                       remaining, stop, details));
                report.reserved_steps += compared;
                return report;
            });
    });
}

core::Result<CaseUniformFlowEvaluation>
evaluate_case_uniform_flow(const CaseUniformFlowReady& field, const CaseUniformFlowReady& current,
                           CaseUniformFlowQuery query, a::UniformFlowWorkspace& workspace,
                           p::CurveLimits limits, core::ResourcePlan& memory,
                           const std::stop_token& stop, core::FailureDetails details) noexcept {
    return guarded<CaseUniformFlowEvaluation>(evaluation, details, [&] {
        return with_current(
            field, current, query, limits, stop, details,
            [&](a::UniformFlowQuery request, p::CurveLimits remaining, std::uint64_t compared) {
                auto jet = take(a::evaluate_uniform_flow(field.uniform(), request, workspace,
                                                         remaining, memory, stop, details));
                const auto used = compared + jet.used_steps();
                return CaseUniformFlowEvaluation{std::move(jet), used};
            });
    });
}

core::Result<a::UniformFlowSnapshot>
snapshot_case_uniform_flow(const CaseUniformFlowReady& field, a::UniformFlowWorkspace& workspace,
                           p::CurveLimits limits, core::ResourcePlan& memory,
                           const std::stop_token& stop, core::FailureDetails details) noexcept {
    return a::snapshot_uniform_flow(field.uniform(), workspace, limits, memory, stop, details);
}

core::Result<CaseUniformFlowReady>
restore_case_uniform_flow(const a::UniformFlowSnapshotView& saved,
                          const CaseUniformFlowReady& current, a::UniformFlowWorkspace& workspace,
                          p::CurveLimits limits, core::ResourcePlan& memory,
                          const std::stop_token& stop, core::FailureDetails details) noexcept {
    return guarded<CaseUniformFlowReady>(evaluation, details, [&] {
        return case_flow_detail::Builder::restore(saved, current, workspace, limits, memory, stop,
                                                  details);
    });
}

} // namespace uvlas::coupling::detail
