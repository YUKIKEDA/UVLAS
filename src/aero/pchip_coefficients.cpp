#include "pchip_coefficients.hpp"

#include "pchip_slot_private.hpp"

#include <algorithm>
#include <array>
namespace uvlas::aero::detail::pchip {
core::Result<void> Counter::cover_alternative(const Counter& other) noexcept {
    const auto& a = plan_.limits_;
    const auto& b = other.plan_.limits_;
    if (a.input_bits != b.input_bits || a.work_bits != b.work_bits || a.steps != b.steps ||
        a.slots != b.slots)
        return std::unexpected(slot_detail::failure(core::FailureCode::InvalidInput, details_));
    auto next = plan_;
    slot_detail::CheckedCount c;
    next.slots_ = 0;
    for (std::size_t i = 0; i < next.counts_.size(); ++i) {
        next.counts_[i] = std::max(next.counts_[i], other.plan_.counts_[i]);
        if (i != static_cast<std::size_t>(RationalOp::Compare))
            next.slots_ = c.add(next.slots_, next.counts_[i]);
    }
    next.bits_ = std::max(next.bits_, other.plan_.bits_);
    if (!c.ok || next.slots_ > next.limits_.slots || next.bits_ > next.limits_.work_bits)
        return std::unexpected(slot_detail::failure(core::FailureCode::ResourceExceeded, details_));
    plan_ = next;
    return {};
}
namespace coefficient_detail {
using slot_detail::CheckedCount;
using slot_detail::failure;
using slot_detail::supplement;
struct Record {
    std::uint64_t numerator = 0, denominator = 0;
    int sign = 0;
};
struct Owner {
    std::uint64_t n = 0, m = 0, w = 0, work = 0, slots = 0;
    core::FailureDetails details;
    core::PlannedBuffer<Record> records;
    core::PlannedBuffer<std::uint32_t> limbs;
};
struct Value {
    BitBound bound;
    std::optional<Slot> slot;
};
BitBound shape(RationalView v) noexcept {
    const auto bits = [](std::span<const std::uint32_t> words) {
        return words.empty() ? std::uint64_t{0}
                             : static_cast<std::uint64_t>(words.size() - 1) * 32 +
                                   static_cast<std::uint64_t>(std::bit_width(words.back()));
    };
    return {bits(v.numerator), bits(v.denominator)};
}
struct Machine {
    using V = Value;
    SampleGrid samples;
    Counter counter;
    Session* session = nullptr;
    std::span<std::optional<Slot>> handles;
    std::optional<core::Failure> error;
    void reject(core::FailureCode code) {
        if (!error)
            error = failure(code);
    }
    V accept(core::Result<BitBound> result) {
        if (!result) {
            if (!error)
                error = result.error();
            return {};
        }
        return {*result, {}};
    }
    V accept(core::Result<Slot> result) {
        if (!result) {
            if (!error)
                error = result.error();
            return {};
        }
        return {{}, *result};
    }
    V read(std::uint64_t index) {
        if (error)
            return {};
        if (session) {
            if (index >= handles.size() || !handles[static_cast<std::size_t>(index)]) {
                reject(core::FailureCode::InvalidInput);
                return {};
            }
            return {{}, handles[static_cast<std::size_t>(index)]};
        }
        return {shape(index < samples.times.size()
                          ? samples.times[static_cast<std::size_t>(index)]
                          : samples.values[static_cast<std::size_t>(index - samples.times.size())]),
                {}};
    }
    V time(std::uint64_t i) {
        return read(i);
    }
    V sample(std::uint64_t i, std::uint64_t component) {
        return read(static_cast<std::uint64_t>(samples.times.size()) + i * samples.components +
                    component);
    }
    V integer(std::int64_t i) {
        if (error)
            return {};
        return session ? accept(session->integer(i)) : accept(counter.integer(i));
    }
    V binary(RationalOp op, V a, V b) {
        if (error)
            return {};
        if (session) {
            if (!a.slot || !b.slot) {
                reject(core::FailureCode::InvalidInput);
                return {};
            }
            return accept(session->binary(op, *a.slot, *b.slot));
        }
        return accept(counter.binary(op, a.bound, b.bound));
    }
    V add(V a, V b) {
        return binary(RationalOp::Add, a, b);
    }
    V subtract(V a, V b) {
        return binary(RationalOp::Subtract, a, b);
    }
    V multiply(V a, V b) {
        return binary(RationalOp::Multiply, a, b);
    }
    V divide(V a, V b) {
        return binary(RationalOp::Divide, a, b);
    }
    V negate(V a) {
        if (error)
            return {};
        if (session) {
            if (!a.slot) {
                reject(core::FailureCode::InvalidInput);
                return {};
            }
            return accept(session->unary(RationalOp::Negate, *a.slot));
        }
        return accept(counter.unary(RationalOp::Negate, a.bound));
    }
    int compare(V a, V b) {
        if (error)
            return 0;
        if (session) {
            if (!a.slot || !b.slot) {
                reject(core::FailureCode::InvalidInput);
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
    int sign(V a) {
        return compare(a, integer(0));
    }
    template <class Left, class Right>
    V choose(bool condition, Left left, Right right) {
        if (error)
            return {};
        if (session)
            return condition ? left(*this) : right(*this);
        auto a = *this, b = *this;
        const auto x = left(a), y = right(b);
        if (a.error) {
            error = a.error;
            return {};
        }
        if (b.error) {
            error = b.error;
            return {};
        }
        const auto joined = a.counter.cover_alternative(b.counter);
        if (!joined) {
            error = joined.error();
            return {};
        }
        counter = a.counter;
        return {{std::max(x.bound.numerator, y.bound.numerator),
                 std::max(x.bound.denominator, y.bound.denominator)},
                {}};
    }
};
// One expression definition drives both the no-allocation branch forecast and exact arithmetic.
template <class M>
typename M::V width(M& m, std::uint64_t i) {
    return m.subtract(m.time(i + 1), m.time(i));
}
template <class M>
typename M::V secant(M& m, std::uint64_t i, std::uint64_t c) {
    return m.divide(m.subtract(m.sample(i + 1, c), m.sample(i, c)), width(m, i));
}
template <class M>
typename M::V absolute(M& m, typename M::V value, int sign) {
    return m.choose(sign < 0, [&](M& a) { return a.negate(value); }, [&](M&) { return value; });
}
template <class M>
typename M::V slope(M& m, std::uint64_t index, std::uint64_t c) {
    const auto n = static_cast<std::uint64_t>(m.samples.times.size());
    if (n == 2)
        return secant(m, 0, c);
    if (index == 0 || index == n - 1) {
        const auto first = index == 0 ? std::uint64_t{0} : n - 2,
                   other = index == 0 ? std::uint64_t{1} : n - 3;
        const auto h = width(m, first), h1 = width(m, other), d = secant(m, first, c),
                   d1 = secant(m, other, c);
        const auto z = m.divide(
            m.subtract(m.multiply(m.add(m.multiply(m.integer(2), h), h1), d), m.multiply(h, d1)),
            m.add(h, h1));
        const int sz = m.sign(z), sd = m.sign(d), sd1 = m.sign(d1);
        const auto three_d = m.multiply(m.integer(3), d);
        const auto magnitude = absolute(m, z, sz), limit = absolute(m, three_d, sd);
        const int clipped = m.compare(magnitude, limit);
        return m.choose(
            sz != sd, [&](M& a) { return a.integer(0); },
            [&](M& a) {
                return a.choose(
                    sd != sd1 && clipped > 0, [&](M&) { return three_d; }, [&](M&) { return z; });
            });
    }
    const auto left = secant(m, index - 1, c), right = secant(m, index, c);
    const int sl = m.sign(left), sr = m.sign(right);
    return m.choose(
        sl == 0 || sr == 0 || sl != sr, [&](M& a) { return a.integer(0); },
        [&](M& a) {
            const auto before = width(a, index - 1), after = width(a, index);
            const auto w1 = a.add(a.multiply(a.integer(2), after), before),
                       w2 = a.add(after, a.multiply(a.integer(2), before));
            return a.divide(a.add(w1, w2), a.add(a.divide(w1, left), a.divide(w2, right)));
        });
}
template <class M>
std::array<typename M::V, 4> interval(M& m, std::uint64_t i, std::uint64_t c) {
    const auto h = width(m, i), a = m.sample(i, c), b = m.sample(i + 1, c), left = slope(m, i, c),
               right = slope(m, i + 1, c);
    return {a, m.multiply(h, left),
            m.subtract(m.multiply(m.integer(3), m.subtract(b, a)),
                       m.multiply(h, m.add(m.multiply(m.integer(2), left), right))),
            m.add(m.multiply(m.integer(2), m.subtract(a, b)), m.multiply(h, m.add(left, right)))};
}
struct Meter {
    std::uint64_t cap = 0, used = 0;
    const std::stop_token& stop;
    core::FailureDetails details;
    std::optional<core::Failure> error;
    bool visit(std::uint64_t n = 1) {
        if (error)
            return false;
        if (stop.stop_requested()) {
            error = failure(core::FailureCode::Cancelled, details);
            return false;
        }
        if (n > cap - used) {
            error = failure(core::FailureCode::ResourceExceeded, details);
            return false;
        }
        used += n;
        return true;
    }
};
struct Builder {
    static core::Result<CoefficientPlan> forecast(SampleGrid samples, RationalLimits limits,
                                                  const std::stop_token& stop,
                                                  core::FailureDetails details) noexcept {
        CoefficientPlan plan;
        CheckedCount checked;
        plan.n_ = samples.times.size();
        plan.m_ = samples.components;
        const auto nm = checked.mul(plan.n_, plan.m_);
        plan.inputs_ = checked.add(plan.n_, nm);
        plan.records_ = checked.add(
            plan.inputs_, checked.mul(checked.mul(plan.n_ > 0 ? plan.n_ - 1 : 0, plan.m_), 4));
        if (!checked.ok ||
            plan.records_ > static_cast<std::uint64_t>(std::numeric_limits<std::ptrdiff_t>::max()))
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        if (plan.n_ < 2 || plan.m_ == 0 || nm != samples.values.size())
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        Meter scan{limits.steps, 0, stop, details, {}};
        Machine machine{samples, Counter{limits, details}, nullptr, {}, {}};
        for (std::uint64_t i = 0; i < plan.inputs_; ++i) {
            if (!scan.visit())
                return std::unexpected(*scan.error);
            const auto v = i < plan.n_ ? samples.times[static_cast<std::size_t>(i)]
                                       : samples.values[static_cast<std::size_t>(i - plan.n_)];
            // Check limb counts before multiplying them into a bit index or scanning any words.
            const auto max_words = limits.input_bits / 32 + (limits.input_bits % 32 != 0 ? 1U : 0U);
            if (v.numerator.size() > max_words || v.denominator.size() > max_words)
                return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
            if (v.denominator.empty() || v.denominator.back() == 0 ||
                (!v.numerator.empty() && v.numerator.back() == 0) ||
                (v.numerator.empty() ? v.sign != 0 : (v.sign != 1 && v.sign != -1)))
                return std::unexpected(failure(core::FailureCode::InvalidInput, details));
            CheckedCount width_check;
            const auto bit_count = [&](std::span<const std::uint32_t> words) {
                return words.empty()
                           ? std::uint64_t{0}
                           : width_check.add(
                                 width_check.mul(static_cast<std::uint64_t>(words.size() - 1), 32),
                                 static_cast<std::uint64_t>(std::bit_width(words.back())));
            };
            const BitBound bound{bit_count(v.numerator), bit_count(v.denominator)};
            if (!width_check.ok)
                return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
            const auto result = machine.counter.input(bound);
            if (!result)
                return std::unexpected(result.error());
        }
        for (std::uint64_t i = 1; i < plan.n_; ++i) {
            if (!scan.visit())
                return std::unexpected(*scan.error);
            machine.compare(machine.time(i - 1), machine.time(i));
            if (machine.error)
                return std::unexpected(supplement(*machine.error, details));
        }
        for (std::uint64_t i = 0; i < plan.n_ - 1; ++i)
            for (std::uint64_t c = 0; c < plan.m_; ++c) {
                if (!scan.visit())
                    return std::unexpected(*scan.error);
                interval(machine, i, c);
                if (machine.error)
                    return std::unexpected(supplement(*machine.error, details));
            }
        plan.preflight_ = scan.used;
        const auto arithmetic = machine.counter.finish();
        if (!arithmetic)
            return std::unexpected(arithmetic.error());
        plan.arithmetic_ = *arithmetic;
        const auto w = plan.arithmetic_.words();
        plan.limbs_ = checked.mul(checked.mul(2, plan.records_), w);
        // Each record copy visits its metadata plus both used limb arrays, at most 2w.
        // Scan, import, order, interval/component and four coefficient dispatches have slack
        // bounded by 8 visits per input/output record, rather than an uncounted second pass.
        plan.outer_ =
            checked.add(checked.mul(8, checked.add(checked.add(plan.inputs_, plan.records_), 1)),
                        checked.mul(4, checked.mul(plan.records_, w)));
        plan.work_ = checked.add(plan.arithmetic_.reserved_steps(), plan.outer_);
        const auto handles = checked.mul(plan.inputs_, sizeof(std::optional<Slot>));
        const auto record_bytes = checked.mul(plan.records_, sizeof(Record));
        const auto limbs =
            checked.mul(checked.mul(checked.mul(2, plan.records_), w), sizeof(std::uint32_t));
        plan.payload_ = checked.add(
            plan.arithmetic_.payload_bytes(),
            checked.add(handles, checked.add(sizeof(Owner), checked.add(record_bytes, limbs))));
        if (!checked.ok || plan.work_ > limits.steps || plan.preflight_ > plan.outer_)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        return plan;
    }
    static RationalView view(const Owner& owner, std::uint64_t i) {
        const auto& r = owner.records.values()[static_cast<std::size_t>(i)];
        const auto start = static_cast<std::size_t>(2 * i * owner.w);
        return {r.sign, owner.limbs.values().subspan(start, static_cast<std::size_t>(r.numerator)),
                owner.limbs.values().subspan(start + static_cast<std::size_t>(owner.w),
                                             static_cast<std::size_t>(r.denominator))};
    }
    static bool copy(Owner& owner, std::uint64_t i, RationalView value, Meter& meter) {
        if (!meter.visit())
            return false;
        auto& record = owner.records.values()[static_cast<std::size_t>(i)];
        record = {value.numerator.size(), value.denominator.size(), value.sign};
        auto words = owner.limbs.values().subspan(static_cast<std::size_t>(2 * i * owner.w),
                                                  static_cast<std::size_t>(2 * owner.w));
        for (std::size_t j = 0; j < value.numerator.size(); ++j) {
            if (!meter.visit())
                return false;
            words[j] = value.numerator[j];
        }
        for (std::size_t j = 0; j < value.denominator.size(); ++j) {
            if (!meter.visit())
                return false;
            words[static_cast<std::size_t>(owner.w) + j] = value.denominator[j];
        }
        return true;
    }
    static core::Result<CoefficientBlock> prepare(SampleGrid samples,
                                                  CoefficientWorkspace& workspace,
                                                  RationalLimits limits, core::ResourcePlan& memory,
                                                  const std::stop_token& stop,
                                                  core::FailureDetails details) noexcept {
        if (workspace.busy_.test_and_set(std::memory_order_acquire))
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        struct Lease {
            CoefficientWorkspace& workspace;
            ~Lease() {
                workspace.busy_.clear(std::memory_order_release);
            }
        };
        const Lease lease{workspace};
        const auto counted = forecast(samples, limits, stop, details);
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
        owner.n = plan.n_;
        owner.m = plan.m_;
        owner.w = plan.arithmetic_.words();
        owner.details = details;
        owner.records = std::move(*records);
        owner.limbs = std::move(*limbs);
        Meter meter{plan.outer_, plan.preflight_, stop, details, {}};
        // Snapshot every original exact word before canonical GCD checks and numerical work.
        for (std::uint64_t i = 0; i < plan.inputs_; ++i) {
            const auto v = i < plan.n_ ? samples.times[static_cast<std::size_t>(i)]
                                       : samples.values[static_cast<std::size_t>(i - plan.n_)];
            if (!copy(owner, i, v, meter))
                return std::unexpected(*meter.error);
        }
        auto lease_arithmetic = arithmetic->begin(stop, details);
        if (!lease_arithmetic)
            return std::unexpected(lease_arithmetic.error());
        auto& session = *lease_arithmetic;
        for (std::uint64_t i = 0; i < plan.inputs_; ++i) {
            if (!meter.visit())
                return std::unexpected(*meter.error);
            auto value = session.input(view(owner, i));
            if (!value)
                return std::unexpected(value.error());
            handles->values()[static_cast<std::size_t>(i)] = *value;
        }
        Machine machine{samples, Counter{limits, details}, &session, handles->values(), {}};
        for (std::uint64_t i = 1; i < plan.n_; ++i) {
            if (!meter.visit())
                return std::unexpected(*meter.error);
            const auto order = machine.compare(machine.time(i - 1), machine.time(i));
            if (machine.error)
                return std::unexpected(supplement(*machine.error, details));
            if (order >= 0)
                return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        }
        for (std::uint64_t i = 0; i < plan.n_ - 1; ++i)
            for (std::uint64_t c = 0; c < plan.m_; ++c) {
                if (!meter.visit())
                    return std::unexpected(*meter.error);
                const auto values = interval(machine, i, c);
                if (machine.error)
                    return std::unexpected(supplement(*machine.error, details));
                for (std::size_t p = 0; p < 4; ++p) {
                    if (!values[p].slot)
                        return std::unexpected(failure(core::FailureCode::InvalidInput, details));
                    const auto value = session.view(*values[p].slot);
                    if (!value)
                        return std::unexpected(value.error());
                    if (!copy(owner, plan.inputs_ + (i * plan.m_ + c) * 4 + p, *value, meter))
                        return std::unexpected(*meter.error);
                }
            }
        const auto used = session.usage();
        if (!used)
            return std::unexpected(used.error());
        if (!meter.visit(0))
            return std::unexpected(*meter.error);
        owner.work = used->steps + meter.used;
        owner.slots = used->slots;
        return CoefficientBlock{std::move(*header)};
    }
};
} // namespace coefficient_detail
CoefficientBlock::CoefficientBlock(core::PlannedBuffer<coefficient_detail::Owner> owner) noexcept
    : owner_(std::move(owner)) {}
CoefficientBlock::CoefficientBlock(CoefficientBlock&&) noexcept = default;
CoefficientBlock& CoefficientBlock::operator=(CoefficientBlock&&) noexcept = default;
CoefficientBlock::~CoefficientBlock() = default;
bool CoefficientBlock::active() const noexcept {
    return !owner_.values().empty();
}
std::uint64_t CoefficientBlock::nodes() const noexcept {
    return active() ? owner_.values()[0].n : 0;
}
std::uint64_t CoefficientBlock::components() const noexcept {
    return active() ? owner_.values()[0].m : 0;
}
std::uint64_t CoefficientBlock::used_steps() const noexcept {
    return active() ? owner_.values()[0].work : 0;
}
std::uint64_t CoefficientBlock::used_slots() const noexcept {
    return active() ? owner_.values()[0].slots : 0;
}
core::Result<RationalView> CoefficientBlock::time(std::uint64_t i) const& noexcept {
    if (!active() || i >= nodes())
        return std::unexpected(
            slot_detail::failure(core::FailureCode::InvalidInput,
                                 active() ? owner_.values()[0].details : core::FailureDetails{}));
    return coefficient_detail::Builder::view(owner_.values()[0], i);
}
core::Result<RationalView> CoefficientBlock::value(std::uint64_t i,
                                                   std::uint64_t c) const& noexcept {
    if (!active() || i >= nodes() || c >= components())
        return std::unexpected(
            slot_detail::failure(core::FailureCode::InvalidInput,
                                 active() ? owner_.values()[0].details : core::FailureDetails{}));
    return coefficient_detail::Builder::view(owner_.values()[0], nodes() + i * components() + c);
}
core::Result<RationalView> CoefficientBlock::coefficient(std::uint64_t i, std::uint64_t c,
                                                         std::uint64_t p) const& noexcept {
    if (!active() || i >= nodes() - 1 || c >= components() || p >= 4)
        return std::unexpected(
            slot_detail::failure(core::FailureCode::InvalidInput,
                                 active() ? owner_.values()[0].details : core::FailureDetails{}));
    return coefficient_detail::Builder::view(
        owner_.values()[0], nodes() + nodes() * components() + (i * components() + c) * 4 + p);
}
core::Result<CoefficientPlan> forecast_coefficients(SampleGrid samples, RationalLimits limits,
                                                    const std::stop_token& stop,
                                                    core::FailureDetails details) noexcept {
    return coefficient_detail::Builder::forecast(samples, limits, stop, details);
}
core::Result<CoefficientBlock>
prepare_coefficients(SampleGrid samples, CoefficientWorkspace& workspace, RationalLimits limits,
                     core::ResourcePlan& memory, const std::stop_token& stop,
                     core::FailureDetails details) noexcept {
    return coefficient_detail::Builder::prepare(samples, workspace, limits, memory, stop, details);
}
} // namespace uvlas::aero::detail::pchip
