#include "uniform_flow.hpp"

#include "pchip_slot_private.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <type_traits>

namespace uvlas::aero::detail::uniform_flow_detail {
namespace {

using pchip::BitBound;
using pchip::Counter;
using pchip::CurveAdmission;
using pchip::CurveLimits;
using pchip::RationalOp;
using pchip::RationalPlan;
using pchip::RationalView;
using pchip::Session;
using pchip::Slot;
using pchip::slot_detail::CheckedCount;
using pchip::slot_detail::failure;
using pchip::slot_detail::operation;
using pchip::slot_detail::supplement;
using enum core::FailureCode;

constexpr std::size_t components = 3;
constexpr std::size_t jet_records = 1 + 2 * components;
constexpr std::size_t polynomial_header = 4;
constexpr std::uint64_t formula_version = 1;
constexpr auto max_index = static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max());
constexpr std::array<std::uint32_t, 1> one{1};
constexpr RationalView zero{0, {}, one};
constexpr numerics::SiDimension velocity_dimension{{1, 0, -1, 0, 0, 0, 0}};
constexpr std::array dimensions{velocity_dimension, velocity_dimension, velocity_dimension};
constexpr std::uint64_t allocator_estimate = 4 * sizeof(void*);

struct Refusal {
    core::Failure reason;
};

[[noreturn]] void refuse(core::FailureCode code, core::FailureDetails details) {
    throw Refusal{failure(code, details)};
}

template <class T>
T needed(core::Result<T> result) {
    if (!result)
        throw Refusal{result.error()};
    return std::move(*result);
}

template <class T, class Call>
core::Result<T> guarded(core::FailureDetails details, Call call) noexcept {
    return core::guard_entry(operation, InvalidInput, details, [&]() -> core::Result<T> {
        try {
            return call();
        } catch (const Refusal& refused) {
            return std::unexpected(supplement(refused.reason, details));
        }
    });
}

struct Meter {
    std::uint64_t limit;
    const std::stop_token& stop;
    core::FailureDetails details;
    std::uint64_t used = 0;

    void step() {
        if (stop.stop_requested())
            refuse(Cancelled, details);
        if (used == limit)
            refuse(ResourceExceeded, details);
        ++used;
    }
};

void admit(const CurveAdmission& admission, CurveLimits limits, CheckedCount checked,
           core::FailureDetails details) {
    if (!checked.ok || admission.reserved_steps > limits.arithmetic.steps)
        refuse(ResourceExceeded, details);
    const auto selection = limits.metadata.selection();
    const auto* estimate = std::get_if<numerics::linear::EstimateAdmission>(&selection);
    if (!estimate)
        refuse(Unsupported, details);
    if (admission.metadata_estimate_bytes > estimate->bytes)
        refuse(ResourceExceeded, details);
}

void admit_payload(const CurveAdmission& admission, const core::ResourcePlan& memory,
                   core::FailureDetails details) {
    if (memory.remaining_bytes() < admission.payload_bytes)
        refuse(ResourceExceeded, details);
}

template <class T>
core::PlannedBuffer<T> allocate(core::ResourcePlan& memory, std::uint64_t elements,
                                const std::stop_token& stop, core::FailureDetails details) {
    if (stop.stop_requested())
        refuse(Cancelled, details);
    if (elements > max_index)
        refuse(ResourceExceeded, details);
    return needed(memory.allocate<T>(static_cast<core::Index>(elements), operation));
}

BitBound shape(RationalView value, std::uint64_t bits, core::FailureDetails details) {
    const auto words = bits / 32 + (bits % 32 != 0 ? 1U : 0U);
    if (value.numerator.size() > words || value.denominator.size() > words)
        refuse(ResourceExceeded, details);
    if (value.denominator.empty() || value.denominator.back() == 0 ||
        (!value.numerator.empty() && value.numerator.back() == 0) ||
        (value.numerator.empty() ? value.sign != 0 : (value.sign != -1 && value.sign != 1)))
        refuse(InvalidInput, details);
    CheckedCount checked;
    const auto width = [&](std::span<const std::uint32_t> magnitude) {
        return magnitude.empty()
                   ? std::uint64_t{0}
                   : checked.add(checked.mul(magnitude.size() - 1, 32),
                                 static_cast<std::uint64_t>(std::bit_width(magnitude.back())));
    };
    const BitBound result{width(value.numerator), width(value.denominator)};
    if (!checked.ok || std::max(result.numerator, result.denominator) > bits)
        refuse(ResourceExceeded, details);
    return result;
}

RationalView pack(RationalView source, std::span<std::uint32_t>& available, Meter& meter) {
    meter.step();
    const auto copy = [&](std::span<const std::uint32_t> from) {
        if (from.size() > available.size())
            refuse(ResourceExceeded, meter.details);
        const auto result = available.first(from.size());
        for (std::size_t i = 0; i < from.size(); ++i) {
            meter.step();
            result[i] = from[i];
        }
        available = available.subspan(from.size());
        return std::span<const std::uint32_t>{result};
    };
    const auto numerator = copy(source.numerator);
    return {source.sign, numerator, copy(source.denominator)};
}

struct FormulaView {
    std::span<const RationalView> coefficients;
    std::optional<std::array<RationalView, polynomial_header>> time_map;
    std::span<const RationalView> events;
};

FormulaView formula(const UniformFlowInput& input, core::FailureDetails details) {
    if (const auto* constant = std::get_if<ConstantUniformFlow>(&input.series))
        return {constant->velocity, std::nullopt, input.events};
    const auto* polynomial = std::get_if<PolynomialUniformFlow>(&input.series);
    if (!polynomial || polynomial->coefficients.empty() ||
        polynomial->coefficients.size() % components != 0)
        refuse(InvalidInput, details);
    return {polynomial->coefficients,
            std::array{polynomial->origin, polynomial->scale, polynomial->begin, polynomial->end},
            input.events};
}

template <class Engine>
using Value = std::conditional_t<std::is_same_v<Engine, Counter>, BitBound, Slot>;

template <class Engine>
Value<Engine> read(Engine& engine, RationalView value, std::uint64_t bits,
                   core::FailureDetails details) {
    if constexpr (std::is_same_v<Engine, Counter>)
        return needed(engine.input(shape(value, bits, details)));
    else
        return needed(engine.input(value));
}

template <class Engine, class Predicate>
void compare(Engine& engine, const Value<Engine>& a, const Value<Engine>& b, Predicate predicate,
             core::FailureDetails details) {
    if constexpr (std::is_same_v<Engine, Counter>) {
        const auto compared = engine.compare(a, b);
        if (!compared)
            throw Refusal{compared.error()};
    } else if (!predicate(needed(engine.compare(a, b)))) {
        refuse(InvalidInput, details);
    }
}

template <class Engine, class Store>
void qualify_formula(Engine& engine, const FormulaView& source, std::uint64_t bits, Meter& meter,
                     Store store) {
    const auto import = [&](RationalView raw) {
        meter.step();
        const auto value = read(engine, raw, bits, meter.details);
        store(raw);
        return value;
    };
    const auto increasing = [](auto ordering) { return ordering < 0; };
    if (source.time_map) {
        const auto& [origin, scale, begin, end] = *source.time_map;
        (void)import(origin);
        (void)import(scale);
        if (scale.sign != 1)
            refuse(InvalidInput, meter.details);
        const auto first = import(begin);
        const auto last = import(end);
        compare(engine, first, last, increasing, meter.details);
    }
    for (const auto coefficient : source.coefficients)
        (void)import(coefficient);
    std::optional<Value<Engine>> previous;
    for (const auto event : source.events) {
        const auto value = import(event);
        if (previous)
            compare(engine, *previous, value, increasing, meter.details);
        previous = value;
    }
}

template <class Engine>
std::array<Value<Engine>, jet_records> evaluate_formula(Engine& engine, const FormulaView& source,
                                                        RationalView time, std::uint64_t bits,
                                                        Meter& meter) {
    const auto import = [&](RationalView raw) {
        meter.step();
        return read(engine, raw, bits, meter.details);
    };
    const auto t = import(time);
    const auto z = needed(engine.integer(0));
    const auto& coefficients = source.coefficients;
    if (!source.time_map)
        return {t, import(coefficients[0]), import(coefficients[1]), import(coefficients[2]), z, z,
                z};
    const auto& [origin, scale, begin, end] = *source.time_map;
    const auto start = import(begin);
    const auto finish = import(end);
    const auto closed = [](auto ordering) { return ordering <= 0; };
    compare(engine, start, t, closed, meter.details);
    compare(engine, t, finish, closed, meter.details);
    const auto call = [&](RationalOp op, const auto& a, const auto& b) {
        return needed(engine.binary(op, a, b));
    };
    const auto anchor = import(origin);
    const auto dtdu = import(scale);
    const auto u = call(RationalOp::Divide, call(RationalOp::Subtract, t, anchor), dtdu);
    std::array result{t, z, z, z, z, z, z};
    const auto terms = coefficients.size() / components;
    for (std::size_t c = 0; c < components; ++c) {
        auto value = import(coefficients[(terms - 1) * components + c]);
        auto derivative = z;
        for (std::size_t power = terms - 1; power > 0; --power) {
            derivative = call(RationalOp::Add, call(RationalOp::Multiply, derivative, u), value);
            const auto a = import(coefficients[(power - 1) * components + c]);
            value = call(RationalOp::Add, call(RationalOp::Multiply, value, u), a);
        }
        result[1 + c] = value;
        result[1 + components + c] = call(RationalOp::Divide, derivative, dtdu);
    }
    return result;
}

struct Plan {
    CurveAdmission admission;
    RationalPlan arena;
    std::uint64_t limbs = 0;
    std::uint64_t records = 0;
    std::uint64_t preflight = 0;
};

pchip::PreparationRequest pchip_input(const UniformFlowInput& input, const PchipUniformFlow& series,
                                      core::FailureDetails details) {
    if (series.samples.components != components)
        refuse(InvalidInput, details);
    return {series.samples, dimensions, series.boundary, input.events, input.stamp};
}

pchip::PchipQuery pchip_query(UniformFlowQuery query) {
    return {query.time, pchip::JetOrder::First, query.side};
}

} // namespace

struct Formula {
    pchip::PchipInputStamp stamp{pchip::PchipIssuerId{0}, pchip::PchipInputGeneration{0}};
    core::PlannedBuffer<RationalView> records;
    core::PlannedBuffer<std::uint32_t> limbs;
    std::size_t coefficient_count = 0;
    bool polynomial = false;
    [[nodiscard]] FormulaView view() const noexcept {
        const auto raw = records.values();
        if (!polynomial)
            return {raw.first(components), std::nullopt, raw.subspan(components)};
        return {raw.subspan(polynomial_header, coefficient_count),
                std::array{raw[0], raw[1], raw[2], raw[3]},
                raw.subspan(polynomial_header + coefficient_count)};
    }
};

struct Builder {
    struct Lease {
        std::atomic_flag& entered;
        ~Lease() {
            entered.clear(std::memory_order_release);
        }
    };

    static Plan plan_prepare(const FormulaView& source, CurveLimits limits,
                             const std::stop_token& stop, core::FailureDetails details) {
        CheckedCount checked;
        const auto count =
            checked.add(checked.add(source.coefficients.size(), source.events.size()),
                        source.time_map ? polynomial_header : 0);
        if (!checked.ok || count > max_index)
            refuse(ResourceExceeded, details);
        Meter scan{limits.arithmetic.steps, stop, details};
        Counter counter{limits.arithmetic, details};
        std::uint64_t limbs = 0;
        qualify_formula(counter, source, limits.arithmetic.input_bits, scan, [&](auto value) {
            limbs =
                checked.add(limbs, checked.add(value.numerator.size(), value.denominator.size()));
            if (!checked.ok)
                refuse(ResourceExceeded, details);
        });
        if (limbs > max_index)
            refuse(ResourceExceeded, details);
        const auto arena = needed(counter.finish());
        const auto storage = checked.add(checked.add(sizeof(Formula), arena.payload_bytes()),
                                         checked.add(checked.mul(count, sizeof(RationalView)),
                                                     checked.mul(limbs, sizeof(std::uint32_t))));
        const auto body = checked.add(arena.reserved_steps(), checked.add(limbs, 1));
        const auto total = checked.add(body, checked.add(checked.mul(count, 2), scan.used));
        const CurveAdmission admission{
            storage, 7 * allocator_estimate + sizeof(core::PlannedBuffer<Formula>), total,
            arena.integer_bits()};
        admit(admission, limits, checked, details);
        return {admission, arena, limbs, count, scan.used};
    }

    static Plan plan_evaluate(const Formula& source, UniformFlowQuery query, CurveLimits limits,
                              const std::stop_token& stop, core::FailureDetails details) {
        CheckedCount checked;
        Meter scan{limits.arithmetic.steps, stop, details};
        Counter counter{limits.arithmetic, details};
        (void)evaluate_formula(counter, source.view(), query.time, limits.arithmetic.input_bits,
                               scan);
        const auto arena = needed(counter.finish());
        const auto limbs = checked.mul(2 * jet_records, arena.words());
        if (!checked.ok || limbs > max_index)
            refuse(ResourceExceeded, details);
        const auto body = checked.add(arena.reserved_steps(), checked.add(limbs, jet_records + 1));
        const CurveAdmission admission{
            checked.add(arena.payload_bytes(), checked.mul(limbs, sizeof(std::uint32_t))),
            4 * allocator_estimate, checked.add(body, checked.mul(scan.used, 2)),
            arena.integer_bits()};
        admit(admission, limits, checked, details);
        return {admission, arena, limbs, jet_records, scan.used};
    }

    static UniformFlowReady prepare(const UniformFlowInput& input, UniformFlowWorkspace& workspace,
                                    CurveLimits limits, core::ResourcePlan& memory,
                                    const std::stop_token& stop, core::FailureDetails details) {
        if (workspace.entered_.test_and_set(std::memory_order_acquire))
            refuse(InvalidInput, details);
        const Lease lease{workspace.entered_};
        if (const auto* series = std::get_if<PchipUniformFlow>(&input.series))
            return UniformFlowReady{needed(
                pchip::prepare_pchip_curve(pchip_input(input, *series, details), workspace.pchip_,
                                           limits, memory, stop, details))};
        const auto source = formula(input, details);
        const auto plan = plan_prepare(source, limits, stop, details);
        admit_payload(plan.admission, memory, details);
        auto holder = std::make_shared<core::PlannedBuffer<Formula>>(
            allocate<Formula>(memory, 1, stop, details));
        auto& target = holder->values()[0];
        target.records = allocate<RationalView>(memory, plan.records, stop, details);
        target.limbs = allocate<std::uint32_t>(memory, plan.limbs, stop, details);
        auto arena = needed(pchip::Workspace::create(plan.arena, memory, details, stop));
        auto session = needed(arena.begin(stop, details));
        Meter meter{plan.admission.reserved_steps - plan.preflight - plan.arena.reserved_steps(),
                    stop, details};
        std::size_t next = 0;
        auto available = target.limbs.values();
        qualify_formula(session, source, limits.arithmetic.input_bits, meter, [&](auto value) {
            target.records.values()[next] = pack(value, available, meter);
            ++next;
        });
        meter.step();
        target.stamp = input.stamp;
        target.polynomial = source.time_map.has_value();
        target.coefficient_count = source.coefficients.size();
        const auto work = plan.preflight + meter.used + needed(session.usage()).steps;
        return UniformFlowReady{std::shared_ptr<const Formula>{std::move(holder), &target}, work};
    }

    static void match_stamp(const UniformFlowReady& field, UniformFlowQuery query,
                            core::FailureDetails details) {
        if (needed(field.input()).stamp != query.stamp)
            refuse(StaleGeneration, details);
        (void)needed(pchip::decide_side(pchip::SideLocation::Single, query.side,
                                        pchip::JetOrder::First, details));
    }

    static CurveAdmission forecast(const UniformFlowReady& field, UniformFlowQuery query,
                                   CurveLimits limits, const std::stop_token& stop,
                                   core::FailureDetails details) {
        match_stamp(field, query, details);
        if (const auto* curve = std::get_if<pchip::PchipCurveReady>(&field.data_)) {
            const auto source = needed(curve->retain_input());
            return needed(pchip::forecast_pchip_evaluation(*curve, source, pchip_query(query),
                                                           limits, stop, details));
        }
        return plan_evaluate(*std::get<0>(field.data_), query, limits, stop, details).admission;
    }

    static ExactUniformFlowJet evaluate(const UniformFlowReady& field, UniformFlowQuery query,
                                        UniformFlowWorkspace& workspace, CurveLimits limits,
                                        core::ResourcePlan& memory, const std::stop_token& stop,
                                        core::FailureDetails details) {
        if (workspace.entered_.test_and_set(std::memory_order_acquire))
            refuse(InvalidInput, details);
        const Lease lease{workspace.entered_};
        match_stamp(field, query, details);
        if (const auto* curve = std::get_if<pchip::PchipCurveReady>(&field.data_)) {
            const auto source = needed(curve->retain_input());
            return ExactUniformFlowJet{needed(
                pchip::evaluate_pchip_curve(*curve, source, pchip_query(query), workspace.pchip_,
                                            limits, memory, stop, details))};
        }
        const auto& parent = std::get<0>(field.data_);
        const auto plan = plan_evaluate(*parent, query, limits, stop, details);
        admit_payload(plan.admission, memory, details);
        ExactUniformFlowJet::FormulaJet result{
            parent, allocate<std::uint32_t>(memory, plan.limbs, stop, details), {}, 0};
        auto arena = needed(pchip::Workspace::create(plan.arena, memory, details, stop));
        auto session = needed(arena.begin(stop, details));
        Meter meter{plan.admission.reserved_steps - plan.preflight - plan.arena.reserved_steps(),
                    stop, details};
        const auto values = evaluate_formula(session, parent->view(), query.time,
                                             limits.arithmetic.input_bits, meter);
        auto available = result.words.values();
        for (std::size_t i = 0; i < result.records.size(); ++i)
            result.records[i] = pack(needed(session.view(values[i])), available, meter);
        meter.step();
        result.steps = plan.preflight + meter.used + needed(session.usage()).steps;
        return ExactUniformFlowJet{std::move(result)};
    }

    static CurveAdmission plan_snapshot(const UniformFlowReady& field, CurveLimits limits,
                                        const std::stop_token& stop, core::FailureDetails details) {
        if (!field.active())
            refuse(InvalidInput, details);
        if (const auto* curve = std::get_if<pchip::PchipCurveReady>(&field.data_))
            return needed(pchip::forecast_pchip_snapshot(*curve, limits, stop, details));
        Meter meter{limits.arithmetic.steps, stop, details};
        meter.step();
        return {0, 0, meter.used};
    }

    static UniformFlowSnapshot snapshot(const UniformFlowReady& field,
                                        UniformFlowWorkspace& workspace, CurveLimits limits,
                                        core::ResourcePlan& memory, const std::stop_token& stop,
                                        core::FailureDetails details) {
        if (workspace.entered_.test_and_set(std::memory_order_acquire))
            refuse(InvalidInput, details);
        const Lease lease{workspace.entered_};
        if (const auto* curve = std::get_if<pchip::PchipCurveReady>(&field.data_))
            return UniformFlowSnapshot{needed(pchip::snapshot_pchip_curve(
                *curve, workspace.pchip_, limits, memory, stop, details))};
        const auto plan = plan_snapshot(field, limits, stop, details);
        return UniformFlowSnapshot{UniformFlowReady{std::get<0>(field.data_), plan.reserved_steps}};
    }

    static CurveAdmission plan_formula_restore(const UniformFormulaSnapshotView& saved,
                                               const UniformFlowReady& current, CurveLimits limits,
                                               const std::stop_token& stop,
                                               core::FailureDetails details) {
        Meter meter{limits.arithmetic.steps, stop, details};
        meter.step();
        if (saved.formula_version != formula_version)
            refuse(InvalidInput, details);
        const auto expected = needed(current.input());
        if (saved.input.stamp != expected.stamp ||
            saved.input.series.index() != expected.series.index())
            refuse(StaleGeneration, details);
        const auto actual = formula(saved.input, details);
        const auto target = formula(expected, details);
        const auto equal_words = [&](auto left, auto right) {
            if (left.size() != right.size())
                refuse(StaleGeneration, details);
            for (std::size_t i = 0; i < left.size(); ++i) {
                meter.step();
                if (left[i] != right[i])
                    refuse(StaleGeneration, details);
            }
        };
        const auto equal_records = [&](std::span<const RationalView> left,
                                       std::span<const RationalView> right) {
            if (left.size() != right.size())
                refuse(StaleGeneration, details);
            for (std::size_t i = 0; i < left.size(); ++i) {
                meter.step();
                if (left[i].sign != right[i].sign)
                    refuse(StaleGeneration, details);
                equal_words(left[i].numerator, right[i].numerator);
                equal_words(left[i].denominator, right[i].denominator);
            }
        };
        if (actual.time_map && target.time_map)
            equal_records(*actual.time_map, *target.time_map);
        equal_records(actual.coefficients, target.coefficients);
        equal_records(actual.events, target.events);
        meter.step();
        return {0, 0, meter.used};
    }

    static CurveAdmission plan_restore(const UniformFlowSnapshotView& saved,
                                       const UniformFlowReady& current, CurveLimits limits,
                                       const std::stop_token& stop, core::FailureDetails details) {
        if (!current.active())
            refuse(InvalidInput, details);
        if (const auto* saved_pchip = std::get_if<pchip::PchipSnapshotView>(&saved)) {
            const auto* curve = std::get_if<pchip::PchipCurveReady>(&current.data_);
            if (!curve)
                refuse(StaleGeneration, details);
            const auto source = needed(curve->retain_input());
            return needed(
                pchip::forecast_pchip_restore(*saved_pchip, source, limits, stop, details));
        }
        return plan_formula_restore(std::get<0>(saved), current, limits, stop, details);
    }

    static UniformFlowReady restore(const UniformFlowSnapshotView& saved,
                                    const UniformFlowReady& current,
                                    UniformFlowWorkspace& workspace, CurveLimits limits,
                                    core::ResourcePlan& memory, const std::stop_token& stop,
                                    core::FailureDetails details) {
        if (workspace.entered_.test_and_set(std::memory_order_acquire))
            refuse(InvalidInput, details);
        const Lease lease{workspace.entered_};
        if (!current.active())
            refuse(InvalidInput, details);
        if (const auto* saved_pchip = std::get_if<pchip::PchipSnapshotView>(&saved)) {
            const auto* curve = std::get_if<pchip::PchipCurveReady>(&current.data_);
            if (!curve)
                refuse(StaleGeneration, details);
            const auto source = needed(curve->retain_input());
            return UniformFlowReady{needed(pchip::restore_pchip_curve(
                *saved_pchip, source, workspace.pchip_, limits, memory, stop, details))};
        }
        const auto plan = plan_formula_restore(std::get<0>(saved), current, limits, stop, details);
        return UniformFlowReady{std::get<0>(current.data_), plan.reserved_steps};
    }
};

} // namespace uvlas::aero::detail::uniform_flow_detail

namespace uvlas::aero::detail {
namespace {
using uniform_flow_detail::Builder;
using uniform_flow_detail::components;
using uniform_flow_detail::failure;
using uniform_flow_detail::Formula;
using uniform_flow_detail::formula_version;
using uniform_flow_detail::guarded;
using uniform_flow_detail::needed;
using uniform_flow_detail::zero;
constexpr auto invalid = core::FailureCode::InvalidInput;
} // namespace

UniformFlowReady::UniformFlowReady(std::shared_ptr<const Formula> owner,
                                   std::uint64_t steps) noexcept
    : data_(std::move(owner)), formula_work_(steps) {}
UniformFlowReady::UniformFlowReady(pchip::PchipCurveReady owner) noexcept
    : data_(std::move(owner)) {}

bool UniformFlowReady::active() const noexcept {
    if (const auto* curve = std::get_if<pchip::PchipCurveReady>(&data_))
        return curve->active();
    const auto* formula = std::get_if<std::shared_ptr<const Formula>>(&data_);
    return formula && *formula;
}

core::Result<UniformFlowInput> UniformFlowReady::input() const& noexcept {
    if (const auto* curve = std::get_if<pchip::PchipCurveReady>(&data_)) {
        const auto source = curve->input();
        if (!source)
            return std::unexpected(source.error());
        return UniformFlowInput{source->stamp, source->events,
                                PchipUniformFlow{source->samples, source->boundary}};
    }
    const auto* formula = std::get_if<std::shared_ptr<const Formula>>(&data_);
    if (!formula || !*formula)
        return std::unexpected(failure(invalid));
    const auto& owner = **formula;
    const auto source = owner.view();
    if (!source.time_map)
        return UniformFlowInput{owner.stamp, source.events,
                                ConstantUniformFlow{{source.coefficients[0], source.coefficients[1],
                                                     source.coefficients[2]}}};
    const auto& [origin, scale, begin, end] = *source.time_map;
    return UniformFlowInput{owner.stamp, source.events,
                            PolynomialUniformFlow{origin, scale, begin, end, source.coefficients}};
}

std::uint64_t UniformFlowReady::used_steps() const noexcept {
    if (!active())
        return 0;
    if (const auto* curve = std::get_if<pchip::PchipCurveReady>(&data_))
        return curve->used_steps();
    return formula_work_;
}

ExactUniformFlowJet::ExactUniformFlowJet(FormulaJet jet) noexcept : data_(std::move(jet)) {}
ExactUniformFlowJet::ExactUniformFlowJet(pchip::ExactPchipJet jet) noexcept
    : data_(std::move(jet)) {}

const ExactUniformFlowJet::FormulaJet* ExactUniformFlowJet::active_formula() const noexcept {
    const auto* formula = std::get_if<FormulaJet>(&data_);
    return formula && formula->parent ? formula : nullptr;
}

bool ExactUniformFlowJet::active() const noexcept {
    if (const auto* jet = std::get_if<pchip::ExactPchipJet>(&data_))
        return jet->active();
    return active_formula() != nullptr;
}
core::Result<pchip::PchipInputStamp> ExactUniformFlowJet::stamp() const noexcept {
    if (const auto* jet = std::get_if<pchip::ExactPchipJet>(&data_))
        return jet->stamp();
    if (const auto* formula = active_formula())
        return formula->parent->stamp;
    return std::unexpected(failure(invalid));
}
core::Result<pchip::RationalView> ExactUniformFlowJet::time() const& noexcept {
    if (const auto* jet = std::get_if<pchip::ExactPchipJet>(&data_))
        return jet->original_time();
    if (const auto* formula = active_formula())
        return formula->records[0];
    return std::unexpected(failure(invalid));
}
core::Result<pchip::RationalView>
ExactUniformFlowJet::velocity(std::size_t component) const& noexcept {
    if (component >= components)
        return std::unexpected(failure(invalid));
    if (const auto* jet = std::get_if<pchip::ExactPchipJet>(&data_))
        return jet->derivative(component, pchip::JetOrder::Value);
    if (const auto* formula = active_formula())
        return formula->records[1 + component];
    return std::unexpected(failure(invalid));
}
core::Result<pchip::RationalView>
ExactUniformFlowJet::time_partial(std::size_t component) const& noexcept {
    if (component >= components)
        return std::unexpected(failure(invalid));
    if (const auto* jet = std::get_if<pchip::ExactPchipJet>(&data_))
        return jet->derivative(component, pchip::JetOrder::First);
    if (const auto* formula = active_formula())
        return formula->records[1 + components + component];
    return std::unexpected(failure(invalid));
}
core::Result<pchip::RationalView>
ExactUniformFlowJet::spatial_partial(std::size_t component, std::size_t coordinate) const noexcept {
    if (!active() || component >= components || coordinate >= components)
        return std::unexpected(failure(invalid));
    return zero;
}
std::uint64_t ExactUniformFlowJet::used_steps() const noexcept {
    if (const auto* jet = std::get_if<pchip::ExactPchipJet>(&data_))
        return jet->used_steps();
    if (const auto* formula = active_formula())
        return formula->steps;
    return 0;
}

core::Result<pchip::CurveAdmission>
forecast_uniform_flow_prepare(const UniformFlowInput& input, pchip::CurveLimits limits,
                              const std::stop_token& stop, core::FailureDetails details) noexcept {
    return guarded<pchip::CurveAdmission>(details, [&] {
        if (const auto* curve = std::get_if<PchipUniformFlow>(&input.series))
            return needed(pchip::forecast_pchip_prepare(
                uniform_flow_detail::pchip_input(input, *curve, details), limits, stop, details));
        return Builder::plan_prepare(uniform_flow_detail::formula(input, details), limits, stop,
                                     details)
            .admission;
    });
}
core::Result<UniformFlowReady>
prepare_uniform_flow(const UniformFlowInput& input, UniformFlowWorkspace& workspace,
                     pchip::CurveLimits limits, core::ResourcePlan& memory,
                     const std::stop_token& stop, core::FailureDetails details) noexcept {
    return guarded<UniformFlowReady>(
        details, [&] { return Builder::prepare(input, workspace, limits, memory, stop, details); });
}
core::Result<pchip::CurveAdmission>
forecast_uniform_flow_evaluation(const UniformFlowReady& field, UniformFlowQuery query,
                                 pchip::CurveLimits limits, const std::stop_token& stop,
                                 core::FailureDetails details) noexcept {
    return guarded<pchip::CurveAdmission>(
        details, [&] { return Builder::forecast(field, query, limits, stop, details); });
}
core::Result<ExactUniformFlowJet>
evaluate_uniform_flow(const UniformFlowReady& field, UniformFlowQuery query,
                      UniformFlowWorkspace& workspace, pchip::CurveLimits limits,
                      core::ResourcePlan& memory, const std::stop_token& stop,
                      core::FailureDetails details) noexcept {
    return guarded<ExactUniformFlowJet>(details, [&] {
        return Builder::evaluate(field, query, workspace, limits, memory, stop, details);
    });
}

UniformFlowSnapshot::UniformFlowSnapshot(UniformFlowReady source) noexcept
    : data_(std::move(source)) {}
UniformFlowSnapshot::UniformFlowSnapshot(pchip::PchipSnapshot source) noexcept
    : data_(std::move(source)) {}

core::Result<UniformFlowSnapshotView> UniformFlowSnapshot::view() const& noexcept {
    if (const auto* curve = std::get_if<pchip::PchipSnapshot>(&data_))
        return curve->view().transform([](auto inputs) { return UniformFlowSnapshotView{inputs}; });
    if (const auto* formula = std::get_if<UniformFlowReady>(&data_))
        return formula->input().transform([](auto inputs) {
            return UniformFlowSnapshotView{UniformFormulaSnapshotView{inputs, formula_version}};
        });
    return std::unexpected(failure(invalid));
}

core::Result<pchip::CurveAdmission>
forecast_uniform_flow_snapshot(const UniformFlowReady& field, pchip::CurveLimits limits,
                               const std::stop_token& stop, core::FailureDetails details) noexcept {
    return guarded<pchip::CurveAdmission>(
        details, [&] { return Builder::plan_snapshot(field, limits, stop, details); });
}

core::Result<UniformFlowSnapshot>
snapshot_uniform_flow(const UniformFlowReady& field, UniformFlowWorkspace& workspace,
                      pchip::CurveLimits limits, core::ResourcePlan& memory,
                      const std::stop_token& stop, core::FailureDetails details) noexcept {
    return guarded<UniformFlowSnapshot>(details, [&] {
        return Builder::snapshot(field, workspace, limits, memory, stop, details);
    });
}

core::Result<pchip::CurveAdmission>
forecast_uniform_flow_restore(const UniformFlowSnapshotView& saved, const UniformFlowReady& current,
                              pchip::CurveLimits limits, const std::stop_token& stop,
                              core::FailureDetails details) noexcept {
    return guarded<pchip::CurveAdmission>(
        details, [&] { return Builder::plan_restore(saved, current, limits, stop, details); });
}

core::Result<UniformFlowReady>
restore_uniform_flow(const UniformFlowSnapshotView& saved, const UniformFlowReady& current,
                     UniformFlowWorkspace& workspace, pchip::CurveLimits limits,
                     core::ResourcePlan& memory, const std::stop_token& stop,
                     core::FailureDetails details) noexcept {
    return guarded<UniformFlowReady>(details, [&] {
        return Builder::restore(saved, current, workspace, limits, memory, stop, details);
    });
}

} // namespace uvlas::aero::detail
