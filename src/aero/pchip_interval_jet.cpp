#include "pchip_interval_jet.hpp"

#include "pchip_slot_private.hpp"

#include <array>
namespace uvlas::aero::detail::pchip {
namespace interval_detail {
using slot_detail::CheckedCount;
using slot_detail::failure;
using slot_detail::supplement;
struct Record {
    std::uint64_t numerator = 0, denominator = 0;
    int sign = 0;
};
struct Owner {
    std::uint64_t components = 0, interval = 0, fields = 0, words = 0, work = 0, slots = 0;
    JetOrder order = JetOrder::Value;
    core::FailureDetails details;
    core::PlannedBuffer<Record> records;
    core::PlannedBuffer<std::uint32_t> limbs;
};
struct Meter {
    std::uint64_t cap = 0, used = 0;
    const std::stop_token& stop;
    core::FailureDetails details;
    core::Result<void> visit(std::uint64_t count = 1) {
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));
        if (count > cap - used)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        used += count;
        return {};
    }
};
core::Result<RationalView> input(const CoefficientBlock& source, IntervalRequest query,
                                 std::uint64_t index) {
    if (index == 0)
        return query.time;
    if (index == 1)
        return source.time(query.interval);
    if (index == 2)
        return source.time(query.interval + 1);
    return source.coefficient(query.interval, (index - 3) / 4, (index - 3) % 4);
}
BitBound known_shape(RationalView value) {
    const auto bits = [](std::span<const std::uint32_t> words) {
        return words.empty() ? std::uint64_t{0}
                             : static_cast<std::uint64_t>(words.size() - 1) * 32 +
                                   static_cast<std::uint64_t>(std::bit_width(words.back()));
    };
    return {bits(value.numerator), bits(value.denominator)};
}
core::Result<BitBound> check_shape(RationalView value, RationalLimits limits,
                                   core::FailureDetails details) {
    const auto maximum_words = limits.input_bits / 32 + (limits.input_bits % 32 != 0 ? 1U : 0U);
    if (value.numerator.size() > maximum_words || value.denominator.size() > maximum_words)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    if (value.denominator.empty() || value.denominator.back() == 0 ||
        (!value.numerator.empty() && value.numerator.back() == 0) ||
        (value.numerator.empty() ? value.sign != 0 : (value.sign != 1 && value.sign != -1)))
        return std::unexpected(failure(core::FailureCode::InvalidInput, details));
    CheckedCount checked;
    const auto bits = [&](std::span<const std::uint32_t> words) {
        return words.empty()
                   ? std::uint64_t{0}
                   : checked.add(checked.mul(static_cast<std::uint64_t>(words.size() - 1), 32),
                                 static_cast<std::uint64_t>(std::bit_width(words.back())));
    };
    const BitBound bound{bits(value.numerator), bits(value.denominator)};
    if (!checked.ok || std::max(bound.numerator, bound.denominator) > limits.input_bits)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    return bound;
}
struct Value {
    BitBound bound;
    std::optional<Slot> slot;
};
struct Machine {
    const CoefficientBlock& source;
    IntervalRequest query;
    Counter counter;
    core::FailureDetails details;
    Session* session = nullptr;
    std::span<std::optional<Slot>> handles;
    std::optional<core::Failure> error;
    Machine(const CoefficientBlock& block, IntervalRequest request, RationalLimits limits,
            core::FailureDetails context)
        : source(block), query(request), counter(limits, context), details(context) {}
    Value capture(core::Result<BitBound> result) {
        if (!result) {
            error = result.error();
            return {};
        }
        return {*result, {}};
    }
    Value capture(core::Result<Slot> result) {
        if (!result) {
            error = result.error();
            return {};
        }
        return {{}, *result};
    }
    Value get(std::uint64_t i) {
        if (error)
            return {};
        if (session) {
            if (i >= handles.size()) {
                error = failure(core::FailureCode::InvalidInput, details);
                return {};
            }
            const auto& slot = handles[static_cast<std::size_t>(i)];
            if (!slot) {
                error = failure(core::FailureCode::InvalidInput, details);
                return {};
            }
            return {{}, slot};
        }
        const auto value = input(source, query, i);
        if (!value) {
            error = supplement(value.error(), details);
            return {};
        }
        return {known_shape(*value), {}};
    }
    Value coefficient(std::uint64_t component, std::uint64_t power) {
        return get(3 + 4 * component + power);
    }
    Value integer(std::int64_t value) {
        if (error)
            return {};
        return session ? capture(session->integer(value)) : capture(counter.integer(value));
    }
    Value binary(RationalOp op, Value left, Value right) {
        if (error)
            return {};
        if (session) {
            if (!left.slot || !right.slot) {
                error = failure(core::FailureCode::InvalidInput, details);
                return {};
            }
            return capture(session->binary(op, *left.slot, *right.slot));
        }
        return capture(counter.binary(op, left.bound, right.bound));
    }
    Value add(Value a, Value b) {
        return binary(RationalOp::Add, a, b);
    }
    Value subtract(Value a, Value b) {
        return binary(RationalOp::Subtract, a, b);
    }
    Value multiply(Value a, Value b) {
        return binary(RationalOp::Multiply, a, b);
    }
    Value divide(Value a, Value b) {
        return binary(RationalOp::Divide, a, b);
    }
    int compare(Value a, Value b) {
        if (error)
            return 0;
        if (session) {
            if (!a.slot || !b.slot) {
                error = failure(core::FailureCode::InvalidInput, details);
                return 0;
            }
            const auto result = session->compare(*a.slot, *b.slot);
            if (!result) {
                error = result.error();
                return 0;
            }
            return *result < 0 ? -1 : *result > 0 ? 1 : 0;
        }
        const auto result = counter.compare(a.bound, b.bound);
        if (!result)
            error = result.error();
        return 0;
    }
};
struct Arguments {
    Value h, u, h2, two, three, six;
};
Arguments coordinates(Machine& machine, JetOrder order) {
    const auto h = machine.subtract(machine.get(2), machine.get(1));
    Arguments a{h, machine.divide(machine.subtract(machine.get(0), machine.get(1)), h), {}, {}, {},
                {}};
    if (order >= JetOrder::First) {
        a.two = machine.integer(2);
        a.three = machine.integer(3);
    }
    if (order == JetOrder::Second) {
        a.six = machine.integer(6);
        a.h2 = machine.multiply(h, h);
    }
    return a;
}
// Identical exact expression graph is used by the structural counter and numerical executor.
std::array<Value, 3> component(Machine& machine, const Arguments& a, std::uint64_t c,
                               JetOrder order) {
    const auto c0 = machine.coefficient(c, 0), c1 = machine.coefficient(c, 1),
               c2 = machine.coefficient(c, 2), c3 = machine.coefficient(c, 3);
    std::array<Value, 3> result{
        machine.add(
            machine.multiply(
                machine.add(machine.multiply(machine.add(machine.multiply(c3, a.u), c2), a.u), c1),
                a.u),
            c0),
        {},
        {}};
    if (order >= JetOrder::First) {
        const auto twice_c2 = machine.multiply(a.two, c2);
        result[1] = machine.divide(
            machine.add(
                machine.multiply(
                    machine.add(machine.multiply(machine.multiply(a.three, c3), a.u), twice_c2),
                    a.u),
                c1),
            a.h);
        if (order == JetOrder::Second)
            result[2] = machine.divide(
                machine.add(machine.multiply(machine.multiply(a.six, c3), a.u), twice_c2), a.h2);
    }
    return result;
}
struct Builder {
    static core::Result<IntervalPlan> forecast(const CoefficientBlock& source,
                                               IntervalRequest query, RationalLimits limits,
                                               const std::stop_token& stop,
                                               core::FailureDetails details) noexcept {
        if (!source.active() || source.nodes() < 2 || query.interval >= source.nodes() - 1 ||
            static_cast<std::uint64_t>(query.order) > 2)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        IntervalPlan plan;
        CheckedCount checked;
        plan.m_ = source.components();
        plan.fields_ = checked.add(static_cast<std::uint64_t>(query.order), 1);
        plan.inputs_ = checked.add(checked.mul(4, plan.m_), 3);
        plan.records_ = checked.add(checked.mul(plan.m_, plan.fields_), 2);
        const auto max_index =
            static_cast<std::uint64_t>(std::numeric_limits<std::ptrdiff_t>::max());
        if (!checked.ok || plan.inputs_ > max_index || plan.records_ > max_index)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        Machine machine{source, query, limits, details};
        Meter scan{limits.steps, 0, stop, details};
        for (std::uint64_t i = 0; i < plan.inputs_; ++i) {
            if (const auto visited = scan.visit(); !visited)
                return std::unexpected(visited.error());
            const auto value = input(source, query, i);
            if (!value)
                return std::unexpected(supplement(value.error(), details));
            const auto bound = check_shape(*value, limits, details);
            if (!bound)
                return std::unexpected(bound.error());
            const auto counted = machine.counter.input(*bound);
            if (!counted)
                return std::unexpected(counted.error());
        }
        if (const auto visited = scan.visit(2); !visited)
            return std::unexpected(visited.error());
        machine.compare(machine.get(0), machine.get(1));
        machine.compare(machine.get(0), machine.get(2));
        const auto a = coordinates(machine, query.order);
        for (std::uint64_t c = 0; c < plan.m_; ++c) {
            if (const auto visited = scan.visit(); !visited)
                return std::unexpected(visited.error());
            component(machine, a, c, query.order);
            if (machine.error)
                return std::unexpected(*machine.error);
        }
        const auto arithmetic = machine.counter.finish();
        if (!arithmetic)
            return std::unexpected(arithmetic.error());
        plan.arithmetic_ = *arithmetic;
        plan.preflight_ = scan.used;
        const auto w = plan.arithmetic_.words();
        plan.limbs_ = checked.mul(checked.mul(2, plan.records_), w);
        plan.outer_ =
            checked.add(checked.mul(8, checked.add(checked.add(plan.inputs_, plan.records_), 1)),
                        checked.mul(4, checked.mul(plan.records_, w)));
        plan.work_ = checked.add(plan.arithmetic_.reserved_steps(), plan.outer_);
        plan.payload_ = checked.add(
            plan.arithmetic_.payload_bytes(),
            checked.add(checked.mul(plan.inputs_, sizeof(std::optional<Slot>)),
                        checked.add(sizeof(Owner),
                                    checked.add(checked.mul(plan.records_, sizeof(Record)),
                                                checked.mul(plan.limbs_, sizeof(std::uint32_t))))));
        if (!checked.ok || plan.work_ > limits.steps || plan.preflight_ > plan.outer_)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        return plan;
    }
    static RationalView view(const Owner& owner, std::uint64_t record) {
        const auto& meta = owner.records.values()[static_cast<std::size_t>(record)];
        const auto first = static_cast<std::size_t>(2 * record * owner.words);
        return {meta.sign,
                owner.limbs.values().subspan(first, static_cast<std::size_t>(meta.numerator)),
                owner.limbs.values().subspan(first + static_cast<std::size_t>(owner.words),
                                             static_cast<std::size_t>(meta.denominator))};
    }
    static core::Result<void> copy(Owner& owner, std::uint64_t record, RationalView value,
                                   Meter& meter) {
        if (const auto visited = meter.visit(); !visited)
            return std::unexpected(visited.error());
        auto& meta = owner.records.values()[static_cast<std::size_t>(record)];
        meta = {value.numerator.size(), value.denominator.size(), value.sign};
        auto words =
            owner.limbs.values().subspan(static_cast<std::size_t>(2 * record * owner.words),
                                         static_cast<std::size_t>(2 * owner.words));
        for (std::size_t i = 0; i < value.numerator.size(); ++i) {
            if (const auto visited = meter.visit(); !visited)
                return std::unexpected(visited.error());
            words[i] = value.numerator[i];
        }
        for (std::size_t i = 0; i < value.denominator.size(); ++i) {
            if (const auto visited = meter.visit(); !visited)
                return std::unexpected(visited.error());
            words[static_cast<std::size_t>(owner.words) + i] = value.denominator[i];
        }
        return {};
    }
    static core::Result<IntervalValues>
    evaluate(const CoefficientBlock& source, IntervalRequest query, IntervalWorkspace& workspace,
             RationalLimits limits, core::ResourcePlan& memory, const std::stop_token& stop,
             core::FailureDetails details) noexcept {
        if (workspace.busy_.test_and_set(std::memory_order_acquire))
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        struct Lease {
            IntervalWorkspace& workspace;
            ~Lease() {
                workspace.busy_.clear(std::memory_order_release);
            }
        };
        const Lease lease{workspace};
        const auto counted = forecast(source, query, limits, stop, details);
        if (!counted)
            return std::unexpected(counted.error());
        const auto& plan = *counted;
        if (plan.payload_ > memory.remaining_bytes())
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        auto arithmetic = Workspace::create(plan.arithmetic_, memory, details, stop);
        if (!arithmetic)
            return std::unexpected(arithmetic.error());
        auto handles = memory.allocate<std::optional<Slot>>(static_cast<core::Index>(plan.inputs_),
                                                            core::Operation::EvaluateAero);
        if (!handles)
            return std::unexpected(supplement(handles.error(), details));
        auto header = memory.allocate<Owner>(1, core::Operation::EvaluateAero);
        if (!header)
            return std::unexpected(supplement(header.error(), details));
        auto records = memory.allocate<Record>(static_cast<core::Index>(plan.records_),
                                               core::Operation::EvaluateAero);
        if (!records)
            return std::unexpected(supplement(records.error(), details));
        auto limbs = memory.allocate<std::uint32_t>(static_cast<core::Index>(plan.limbs_),
                                                    core::Operation::EvaluateAero);
        if (!limbs)
            return std::unexpected(supplement(limbs.error(), details));
        auto& owner = header->values()[0];
        owner.components = plan.m_;
        owner.interval = query.interval;
        owner.fields = plan.fields_;
        owner.words = plan.arithmetic_.words();
        owner.order = query.order;
        owner.details = details;
        owner.records = std::move(*records);
        owner.limbs = std::move(*limbs);
        Meter meter{plan.outer_, plan.preflight_, stop, details};
        auto numerical = arithmetic->begin(stop, details);
        if (!numerical)
            return std::unexpected(numerical.error());
        auto& session = *numerical;
        for (std::uint64_t i = 0; i < plan.inputs_; ++i) {
            if (const auto visited = meter.visit(); !visited)
                return std::unexpected(visited.error());
            const auto value = input(source, query, i);
            if (!value)
                return std::unexpected(supplement(value.error(), details));
            const auto imported = session.input(*value);
            if (!imported)
                return std::unexpected(imported.error());
            handles->values()[static_cast<std::size_t>(i)] = *imported;
        }
        Machine machine{source, query, limits, details};
        machine.session = &session;
        machine.handles = handles->values();
        if (const auto visited = meter.visit(2); !visited)
            return std::unexpected(visited.error());
        const auto lower = machine.compare(machine.get(0), machine.get(1)),
                   upper = machine.compare(machine.get(0), machine.get(2));
        if (machine.error)
            return std::unexpected(*machine.error);
        if (lower < 0 || upper > 0)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const auto a = coordinates(machine, query.order);
        if (machine.error)
            return std::unexpected(*machine.error);
        const auto save = [&](std::uint64_t record, const Value& value) -> core::Result<void> {
            const auto& slot = value.slot;
            if (!slot)
                return std::unexpected(failure(core::FailureCode::InvalidInput, details));
            const auto read = session.view(*slot);
            if (!read)
                return std::unexpected(read.error());
            return copy(owner, record, *read, meter);
        };
        if (const auto saved = save(0, machine.get(0)); !saved)
            return std::unexpected(saved.error());
        if (const auto saved = save(1, a.u); !saved)
            return std::unexpected(saved.error());
        for (std::uint64_t c = 0; c < plan.m_; ++c) {
            if (const auto visited = meter.visit(); !visited)
                return std::unexpected(visited.error());
            const auto values = component(machine, a, c, query.order);
            if (machine.error)
                return std::unexpected(*machine.error);
            for (std::uint64_t d = 0; d < plan.fields_; ++d)
                if (const auto saved =
                        save(2 + c * plan.fields_ + d, values[static_cast<std::size_t>(d)]);
                    !saved)
                    return std::unexpected(saved.error());
        }
        const auto used = session.usage();
        if (!used)
            return std::unexpected(used.error());
        if (const auto visited = meter.visit(0); !visited)
            return std::unexpected(visited.error());
        owner.work = used->steps + meter.used;
        owner.slots = used->slots;
        return IntervalValues{std::move(*header)};
    }
    static core::Failure invalid(const core::PlannedBuffer<Owner>& owner) {
        return failure(core::FailureCode::InvalidInput,
                       owner.values().empty() ? core::FailureDetails{} : owner.values()[0].details);
    }
};
} // namespace interval_detail
IntervalValues::IntervalValues(core::PlannedBuffer<interval_detail::Owner> owner) noexcept
    : owner_(std::move(owner)) {}
IntervalValues::IntervalValues(IntervalValues&&) noexcept = default;
IntervalValues& IntervalValues::operator=(IntervalValues&&) noexcept = default;
IntervalValues::~IntervalValues() = default;
bool IntervalValues::active() const noexcept {
    return !owner_.values().empty();
}
std::uint64_t IntervalValues::components() const noexcept {
    return active() ? owner_.values()[0].components : 0;
}
core::Result<JetOrder> IntervalValues::requested_order() const noexcept {
    if (!active())
        return std::unexpected(interval_detail::Builder::invalid(owner_));
    return owner_.values()[0].order;
}
core::Result<std::uint64_t> IntervalValues::interval() const noexcept {
    if (!active())
        return std::unexpected(interval_detail::Builder::invalid(owner_));
    return owner_.values()[0].interval;
}
core::Result<RationalView> IntervalValues::time() const& noexcept {
    if (!active())
        return std::unexpected(interval_detail::Builder::invalid(owner_));
    return interval_detail::Builder::view(owner_.values()[0], 0);
}
core::Result<RationalView> IntervalValues::coordinate_u() const& noexcept {
    if (!active())
        return std::unexpected(interval_detail::Builder::invalid(owner_));
    return interval_detail::Builder::view(owner_.values()[0], 1);
}
core::Result<RationalView> IntervalValues::derivative(std::uint64_t c, JetOrder d) const& noexcept {
    if (!active() || c >= components() ||
        static_cast<std::uint64_t>(d) >= owner_.values()[0].fields)
        return std::unexpected(interval_detail::Builder::invalid(owner_));
    return interval_detail::Builder::view(owner_.values()[0], 2 + c * owner_.values()[0].fields +
                                                                  static_cast<std::uint64_t>(d));
}
std::uint64_t IntervalValues::used_steps() const noexcept {
    return active() ? owner_.values()[0].work : 0;
}
std::uint64_t IntervalValues::used_slots() const noexcept {
    return active() ? owner_.values()[0].slots : 0;
}
core::Result<IntervalPlan> forecast_interval(const CoefficientBlock& source, IntervalRequest query,
                                             RationalLimits limits, const std::stop_token& stop,
                                             core::FailureDetails details) noexcept {
    return interval_detail::Builder::forecast(source, query, limits, stop, details);
}
core::Result<IntervalValues> evaluate_interval(const CoefficientBlock& source,
                                               IntervalRequest query, IntervalWorkspace& workspace,
                                               RationalLimits limits, core::ResourcePlan& memory,
                                               const std::stop_token& stop,
                                               core::FailureDetails details) noexcept {
    return interval_detail::Builder::evaluate(source, query, workspace, limits, memory, stop,
                                              details);
}
} // namespace uvlas::aero::detail::pchip
