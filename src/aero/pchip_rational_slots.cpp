#include "pchip_slot_private.hpp"

#include <algorithm>
#include <type_traits>

namespace uvlas::aero::detail::pchip::slot_detail {
core::Result<std::uint64_t> next_epoch(std::uint64_t previous) noexcept {
    if (previous == std::numeric_limits<std::uint64_t>::max())
        return std::unexpected(failure(core::FailureCode::ResourceExceeded));
    return previous + 1;
}
core::Result<std::uint64_t> workspace_bytes(std::uint64_t slots, std::uint64_t words) noexcept {
    CheckedCount c;
    const auto integers = c.add(c.mul(slots, 2), 4);
    const auto limbs = c.mul(integers, words);
    const auto result = c.add(sizeof(WorkspaceData), c.add(c.mul(slots, sizeof(SlotMetadata)),
                                                           c.mul(limbs, sizeof(std::uint32_t))));
    if (!c.ok || slots > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()) ||
        limbs > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()))
        return std::unexpected(failure(core::FailureCode::ResourceExceeded));
    return result;
}
core::Result<std::uint64_t> opcode_cost(RationalOp op, std::uint64_t bits,
                                        std::uint64_t words) noexcept {
    CheckedCount c;
    const auto M = c.add(c.mul(words, words), c.mul(4, words));
    const auto A = c.add(c.mul(3, words), 3);
    const auto D = c.add(c.add(c.mul(2, words), c.mul(bits, c.add(c.mul(3, words), 3))), 1);
    const auto N = c.add(c.mul(c.add(c.mul(2, bits), 3), D), c.add(c.mul(4, words), 4));
    std::uint64_t body = 0;
    switch (op) {
    case RationalOp::Input:
        body = c.add(N, c.add(c.mul(2, words), 3));
        break;
    case RationalOp::Clone:
    case RationalOp::Negate:
        body = c.add(c.mul(2, words), 1);
        break;
    case RationalOp::Add:
    case RationalOp::Subtract:
        body = c.add(c.add(c.mul(3, M), A), N);
        break;
    case RationalOp::Multiply:
    case RationalOp::Divide:
        body = c.add(c.mul(2, M), N);
        break;
    case RationalOp::Compare:
        body = c.add(c.mul(2, M), A);
        break;
    case RationalOp::Floor:
        body = c.add(c.add(D, A), c.add(c.mul(2, words), 3));
        break;
    default:
        return std::unexpected(failure(core::FailureCode::InvalidInput));
    }
    const auto result =
        c.add(body, 1); // one outer opcode dispatch; inner loop visits are in M/A/D/N.
    if (!c.ok)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded));
    return result;
}
namespace {
core::Failure session_error(WorkspaceData* data) noexcept {
    return data ? Engine{*data}.error() : failure(core::FailureCode::StaleGeneration);
}
bool ready(WorkspaceData* data) noexcept {
    return data && !data->error.has_value();
}
} // namespace
} // namespace uvlas::aero::detail::pchip::slot_detail
namespace uvlas::aero::detail::pchip {
using namespace slot_detail;
Counter::Counter(RationalLimits limits, core::FailureDetails details) noexcept : details_(details) {
    plan_.limits_ = limits;
}
core::Result<BitBound> Counter::count(RationalOp op, BitBound bound,
                                      std::uint64_t intermediate) noexcept {
    if (bound.denominator == 0 || static_cast<std::size_t>(op) >= plan_.counts_.size())
        return std::unexpected(failure(core::FailureCode::InvalidInput, details_));
    auto next = plan_;
    CheckedCount c;
    auto& count = next.counts_[static_cast<std::size_t>(op)];
    count = c.add(count, 1);
    if (op != RationalOp::Compare)
        next.slots_ = c.add(next.slots_, 1);
    next.bits_ = std::max({next.bits_, bound.numerator, bound.denominator, intermediate});
    if (!c.ok || next.bits_ > next.limits_.work_bits || next.slots_ > next.limits_.slots)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details_));
    plan_ = next;
    return bound;
}
core::Result<BitBound> Counter::input(BitBound shape) noexcept {
    if (shape.denominator == 0)
        return std::unexpected(failure(core::FailureCode::InvalidInput, details_));
    if (std::max(shape.numerator, shape.denominator) > plan_.limits_.input_bits)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details_));
    return count(RationalOp::Input, shape, 0);
}
core::Result<BitBound> Counter::integer(std::int64_t value) noexcept {
    const auto absolute = value < 0 ? static_cast<std::uint64_t>(-(value + 1)) + 1
                                    : static_cast<std::uint64_t>(value);
    return input({static_cast<std::uint64_t>(std::bit_width(absolute)), 1});
}
core::Result<BitBound> Counter::unary(RationalOp op, BitBound value) noexcept {
    if (op != RationalOp::Clone && op != RationalOp::Negate)
        return std::unexpected(failure(core::FailureCode::InvalidInput, details_));
    return count(op, value, 0);
}
core::Result<BitBound> Counter::binary(RationalOp op, BitBound a, BitBound b) noexcept {
    if (a.denominator == 0 || b.denominator == 0)
        return std::unexpected(failure(core::FailureCode::InvalidInput, details_));
    CheckedCount c;
    BitBound out;
    if (op == RationalOp::Add || op == RationalOp::Subtract) {
        out = {c.add(std::max(c.add(a.numerator, b.denominator), c.add(b.numerator, a.denominator)),
                     1),
               c.add(a.denominator, b.denominator)};
    } else if (op == RationalOp::Multiply)
        out = {c.add(a.numerator, b.numerator), c.add(a.denominator, b.denominator)};
    else if (op == RationalOp::Divide)
        out = {c.add(a.numerator, b.denominator), c.add(a.denominator, b.numerator)};
    else
        return std::unexpected(failure(core::FailureCode::InvalidInput, details_));
    // A zero numerator is a valid bound; the caller's value branch still decides divide-by-zero.
    if (op == RationalOp::Divide && b.numerator == 0)
        out.denominator = a.denominator;
    if (!c.ok)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details_));
    return count(op, out, 0);
}
core::Result<void> Counter::compare(BitBound a, BitBound b) noexcept {
    if (a.denominator == 0 || b.denominator == 0)
        return std::unexpected(failure(core::FailureCode::InvalidInput, details_));
    CheckedCount c;
    const auto bits =
        std::max(c.add(a.numerator, b.denominator), c.add(b.numerator, a.denominator));
    if (!c.ok)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details_));
    auto result = count(RationalOp::Compare, {0, 1}, bits);
    if (!result)
        return std::unexpected(result.error());
    return {};
}
core::Result<BitBound> Counter::floor(BitBound a) noexcept {
    if (a.denominator == 0)
        return std::unexpected(failure(core::FailureCode::InvalidInput, details_));
    return count(RationalOp::Floor, {a.numerator, 1}, std::max(a.numerator, a.denominator));
}
core::Result<RationalPlan> Counter::finish() const noexcept {
    auto result = plan_;
    CheckedCount c;
    const auto storage_bits = c.add(c.add(result.bits_, 1), 31);
    result.words_ = storage_bits / 32;
    result.work_ = 0;
    for (std::size_t i = 0; i < result.counts_.size(); ++i) {
        if (result.counts_[i] == 0)
            continue;
        const auto cost = opcode_cost(static_cast<RationalOp>(i), result.bits_, result.words_);
        if (!cost)
            return std::unexpected(supplement(cost.error(), details_));
        result.work_ = c.add(result.work_, c.mul(*cost, result.counts_[i]));
    }
    if (!c.ok || result.work_ > result.limits_.steps)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details_));
    const auto bytes = workspace_bytes(result.slots_, result.words_);
    if (!bytes)
        return std::unexpected(supplement(bytes.error(), details_));
    result.payload_ = *bytes;
    result.valid_ = true;
    return result;
}
Workspace::Workspace(core::PlannedBuffer<WorkspaceData> data) noexcept : data_(std::move(data)) {}
Workspace::Workspace(Workspace&&) noexcept = default;
Workspace& Workspace::operator=(Workspace&&) noexcept = default;
Workspace::~Workspace() = default;
bool Workspace::active() const noexcept {
    return !data_.values().empty();
}
core::Result<Workspace> Workspace::create(const RationalPlan& p, core::ResourcePlan& resources,
                                          core::FailureDetails details,
                                          const std::stop_token& stop) noexcept {
    if (!p.valid_)
        return std::unexpected(failure(core::FailureCode::InvalidInput, details));
    if (stop.stop_requested())
        return std::unexpected(failure(core::FailureCode::Cancelled, details));
    if (p.payload_ > resources.remaining_bytes())
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    // The full forecast is checked before even the workspace owner payload is allocated.
    auto state = resources.allocate<WorkspaceData>(1, operation);
    if (!state)
        return std::unexpected(supplement(state.error(), details));
    auto& d = state->values()[0];
    d.plan = p;
    const auto limbs =
        (p.slots_ * 2 + 4) * p.words_; // finish checked this exact count within Index.
    auto words = resources.allocate<std::uint32_t>(static_cast<core::Index>(limbs), operation);
    if (!words)
        return std::unexpected(supplement(words.error(), details));
    d.limbs = std::move(*words);
    auto values = resources.allocate<SlotMetadata>(static_cast<core::Index>(p.slots_), operation);
    if (!values)
        return std::unexpected(supplement(values.error(), details));
    d.values = std::move(*values);
    // An allocator may request cancellation while the last prepared array is being made.
    if (stop.stop_requested())
        return std::unexpected(failure(core::FailureCode::Cancelled, details));
    return Workspace{std::move(*state)};
}
core::Result<Session> Workspace::begin(std::stop_token stop,
                                       core::FailureDetails details) & noexcept {
    if (!active())
        return std::unexpected(failure(core::FailureCode::InvalidInput, details));
    auto& d = data_.values()[0];
    if (d.busy.test_and_set(std::memory_order_acquire))
        return std::unexpected(failure(core::FailureCode::InvalidInput, details));
    if (stop.stop_requested()) {
        d.busy.clear(std::memory_order_release);
        return std::unexpected(failure(core::FailureCode::Cancelled, details));
    }
    const auto epoch = next_epoch(d.epoch);
    if (!epoch) {
        d.busy.clear(std::memory_order_release);
        return std::unexpected(supplement(epoch.error(), details));
    }
    d.epoch = *epoch;
    d.next = 0;
    d.steps = 0;
    d.error.reset();
    d.details = details;
    d.stop = std::move(stop);
    d.remaining = d.plan.counts_;
    d.scratch_used = {};
    return Session{&d};
}
Session::Session(WorkspaceData* data) noexcept : data_(data) {}
Session::Session(Session&& source) noexcept : data_(std::exchange(source.data_, nullptr)) {}
Session& Session::operator=(Session&& source) noexcept {
    if (this != &source) {
        if (data_)
            data_->busy.clear(std::memory_order_release);
        data_ = std::exchange(source.data_, nullptr);
    }
    return *this;
}
Session::~Session() {
    if (data_)
        data_->busy.clear(std::memory_order_release);
}
core::Result<Slot> Session::input(RationalView input) noexcept {
    if (!ready(data_))
        return std::unexpected(session_error(data_));
    Engine e{*data_};
    if (!e.enter(RationalOp::Input))
        return std::unexpected(e.error());
    auto out = e.allocate();
    if (!out)
        return out;
    if (!e.import(e.value(*out), input))
        return std::unexpected(e.error());
    return out;
}
core::Result<Slot> Session::integer(std::int64_t input) noexcept {
    if (!ready(data_))
        return std::unexpected(session_error(data_));
    Engine e{*data_};
    const auto magnitude = input < 0 ? static_cast<std::uint64_t>(-(input + 1)) + 1
                                     : static_cast<std::uint64_t>(input);
    if (static_cast<std::uint64_t>(std::bit_width(magnitude)) > data_->plan.limits_.input_bits) {
        e.reject(core::FailureCode::ResourceExceeded);
        return std::unexpected(e.error());
    }
    if (!e.enter(RationalOp::Input))
        return std::unexpected(e.error());
    auto out = e.allocate();
    if (!out)
        return out;
    auto r = e.value(*out);
    if (!e.integer(r.numerator, magnitude) || !e.integer(r.denominator, 1))
        return std::unexpected(e.error());
    r.sign = input < 0 ? -1 : input > 0 ? 1 : 0;
    return out;
}
core::Result<Slot> Session::unary(RationalOp op, const Slot& a) noexcept {
    if (!ready(data_))
        return std::unexpected(session_error(data_));
    Engine e{*data_};
    if (op != RationalOp::Clone && op != RationalOp::Negate) {
        e.reject(core::FailureCode::InvalidInput);
        return std::unexpected(e.error());
    }
    if (!e.validate(a) || !e.enter(op))
        return std::unexpected(e.error());
    auto out = e.allocate();
    if (!out)
        return out;
    const auto source = e.view(a);
    auto r = e.value(*out);
    if (!e.copy(r.numerator, {source.numerator}) || !e.copy(r.denominator, {source.denominator}))
        return std::unexpected(e.error());
    r.sign = op == RationalOp::Negate ? -source.sign : source.sign;
    return out;
}
core::Result<Slot> Session::binary(RationalOp op, const Slot& a, const Slot& b) noexcept {
    if (!ready(data_))
        return std::unexpected(session_error(data_));
    Engine e{*data_};
    if (op != RationalOp::Add && op != RationalOp::Subtract && op != RationalOp::Multiply &&
        op != RationalOp::Divide) {
        e.reject(core::FailureCode::InvalidInput);
        return std::unexpected(e.error());
    }
    if (!e.validate(a) || !e.validate(b) || !e.enter(op))
        return std::unexpected(e.error());
    auto out = e.allocate();
    if (!out)
        return out;
    if (!e.arithmetic(op, e.value(*out), e.view(a), e.view(b)))
        return std::unexpected(e.error());
    return out;
}
core::Result<std::strong_ordering> Session::compare(const Slot& a, const Slot& b) noexcept {
    if (!ready(data_))
        return std::unexpected(session_error(data_));
    Engine e{*data_};
    if (!e.validate(a) || !e.validate(b) || !e.enter(RationalOp::Compare))
        return std::unexpected(e.error());
    auto order = std::strong_ordering::equal;
    if (!e.compare(e.view(a), e.view(b), order))
        return std::unexpected(e.error());
    return order;
}
core::Result<Slot> Session::floor(const Slot& a) noexcept {
    if (!ready(data_))
        return std::unexpected(session_error(data_));
    Engine e{*data_};
    if (!e.validate(a) || !e.enter(RationalOp::Floor))
        return std::unexpected(e.error());
    auto out = e.allocate();
    if (!out)
        return out;
    if (!e.floor(e.value(*out), e.view(a)))
        return std::unexpected(e.error());
    return out;
}
core::Result<RationalView> Session::view(const Slot& a) const& noexcept {
    if (!ready(data_))
        return std::unexpected(session_error(data_));
    Engine e{*data_};
    if (!e.validate(a) || !e.step(0))
        return std::unexpected(e.error());
    return e.view(a);
}
core::Result<SlotUsage> Session::usage() const noexcept {
    if (!data_)
        return std::unexpected(failure(core::FailureCode::StaleGeneration));
    return SlotUsage{data_->steps, data_->next};
}
} // namespace uvlas::aero::detail::pchip
