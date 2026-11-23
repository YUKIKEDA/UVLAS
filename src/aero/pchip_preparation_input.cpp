#include "pchip_preparation_input.hpp"

#include "pchip_slot_private.hpp"

#include <algorithm>
#include <array>
#include <type_traits>

namespace uvlas::aero::detail::pchip {
namespace preparation_detail {
using slot_detail::CheckedCount;
using slot_detail::failure;
using slot_detail::supplement;

enum class Policy : std::uint8_t { Reject, Hold, Periodic };

struct Layout {
    std::size_t times, values, events;
    Policy policy;

    [[nodiscard]] std::size_t event_offset() const noexcept {
        return times + values;
    }
    [[nodiscard]] std::size_t demand() const noexcept {
        return event_offset() + events;
    }
    [[nodiscard]] std::size_t period() const noexcept {
        return demand() + 2;
    }
};

struct Owner {
    Layout layout{};
    PchipInputStamp stamp{PchipIssuerId{0}, PchipInputGeneration{0}};
    std::uint64_t formula = CoefficientBlock::formula_version();
    std::uint64_t work = 0;
    core::FailureDetails details;
    core::PlannedBuffer<RationalView> records;
    core::PlannedBuffer<std::uint32_t> limbs;
    core::PlannedBuffer<numerics::SiDimension> dimensions;
};

struct Visits {
    std::uint64_t cap = 0, used = 0;
    const std::stop_token& stop;
    core::FailureDetails details;

    [[nodiscard]] core::Result<void> next() noexcept {
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));
        if (used == cap)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        ++used;
        return {};
    }
};

[[nodiscard]] Policy which(const BoundaryRequest& request) noexcept {
    if (std::holds_alternative<PeriodicBoundary>(request.policy))
        return Policy::Periodic;
    return std::holds_alternative<HoldBoundary>(request.policy) ? Policy::Hold : Policy::Reject;
}

[[nodiscard]] RationalView raw_record(const PreparationRequest& request, const Layout& layout,
                                      std::size_t index) noexcept {
    if (index < layout.times)
        return request.samples.times[index];
    if (index < layout.event_offset())
        return request.samples.values[index - layout.times];
    if (index < layout.demand())
        return request.events[index - layout.event_offset()];
    if (index == layout.demand())
        return request.boundary.begin;
    if (index == layout.demand() + 1)
        return request.boundary.end;
    // Only a Periodic layout contains a record at this index.
    const auto* policy = std::get_if<PeriodicBoundary>(&request.boundary.policy);
    return policy == nullptr ? RationalView{} : policy->period;
}

[[nodiscard]] core::Result<BitBound> shape(RationalView value, std::uint64_t limit,
                                           core::FailureDetails details) noexcept {
    const auto words = limit / 32 + (limit % 32 == 0 ? 0U : 1U);
    if (value.numerator.size() > words || value.denominator.size() > words)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    if (value.denominator.empty() || value.denominator.back() == 0 ||
        (!value.numerator.empty() && value.numerator.back() == 0) ||
        (value.numerator.empty() ? value.sign != 0 : (value.sign != -1 && value.sign != 1)))
        return std::unexpected(failure(core::FailureCode::InvalidInput, details));
    CheckedCount checked;
    const auto bits = [&](std::span<const std::uint32_t> magnitude) {
        return magnitude.empty()
                   ? std::uint64_t{0}
                   : checked.add(checked.mul(magnitude.size() - 1, 32),
                                 static_cast<std::uint64_t>(std::bit_width(magnitude.back())));
    };
    const BitBound result{bits(value.numerator), bits(value.denominator)};
    if (!checked.ok || std::max(result.numerator, result.denominator) > limit)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    return result;
}

template <class Arithmetic, class Reader>
class Program {
    static constexpr bool counting = std::is_same_v<Arithmetic, Counter>;
    using Value = std::conditional_t<counting, BitBound, Slot>;

public:
    Program(Arithmetic& arithmetic, Reader read, Visits& visits, RationalLimits limits) noexcept
        : arithmetic_(arithmetic), read_(read), visits_(visits), limits_(limits) {}

    core::Result<void> check(const Layout& layout) noexcept {
        const auto knots = increasing(0, layout.times);
        if (!knots)
            return std::unexpected(knots.error());
        for (std::size_t i = layout.times; i < layout.event_offset(); ++i) {
            const auto value = input(i);
            if (!value)
                return std::unexpected(value.error());
        }
        if (layout.events > 0) {
            const auto events = increasing(layout.event_offset(), layout.events);
            if (!events)
                return std::unexpected(events.error());
        }
        const auto begin = input(layout.demand());
        if (!begin)
            return std::unexpected(begin.error());
        const auto end = input(layout.demand() + 1);
        if (!end)
            return std::unexpected(end.error());
        const auto ordered = compare(*begin, *end).and_then([&](const auto relation) {
            return relation <= 0 ? core::Result<void>{} : invalid();
        });
        if (!ordered || layout.policy != Policy::Periodic)
            return ordered;

        const auto period = input(layout.period());
        if (!period)
            return std::unexpected(period.error());
        const auto width = visits_.next().and_then(
            [&] { return arithmetic_.binary(RationalOp::Subtract, (*knots)[1], (*knots)[0]); });
        if (!width)
            return std::unexpected(width.error());
        // Strictly increasing knots prove that matching the width also proves positivity.
        return compare(*period, *width).and_then([&](const auto relation) {
            return relation == 0 ? core::Result<void>{} : invalid();
        });
    }

private:
    Arithmetic& arithmetic_;
    Reader read_;
    Visits& visits_;
    RationalLimits limits_;

    [[nodiscard]] core::Result<std::array<Value, 2>> increasing(std::size_t start,
                                                                std::size_t length) noexcept {
        const auto first = input(start);
        if (!first)
            return std::unexpected(first.error());
        auto previous = *first;
        for (std::size_t i = 1; i < length; ++i) {
            const auto next = input(start + i);
            if (!next)
                return std::unexpected(next.error());
            const auto relation = compare(previous, *next);
            if (!relation)
                return std::unexpected(relation.error());
            if constexpr (!counting) {
                if (*relation >= 0)
                    return std::unexpected(
                        failure(core::FailureCode::InvalidInput, visits_.details));
            }
            previous = *next;
        }
        return std::array{*first, previous};
    }

    [[nodiscard]] core::Result<Value> input(std::size_t index) noexcept {
        const auto visit = visits_.next();
        if (!visit)
            return std::unexpected(visit.error());
        const auto value = read_(index);
        if constexpr (counting) {
            return shape(value, limits_.input_bits, visits_.details)
                .and_then([&](const auto bound) { return arithmetic_.input(bound); });
        } else {
            return arithmetic_.input(value);
        }
    }

    [[nodiscard]] core::Result<std::strong_ordering> compare(const Value& left,
                                                             const Value& right) noexcept {
        const auto visit = visits_.next();
        if (!visit)
            return std::unexpected(visit.error());
        if constexpr (counting) {
            return arithmetic_.compare(left, right).transform([] {
                return std::strong_ordering::equal;
            });
        } else {
            return arithmetic_.compare(left, right);
        }
    }

    [[nodiscard]] core::Result<void> invalid() const noexcept {
        return std::unexpected(failure(core::FailureCode::InvalidInput, visits_.details));
    }
};

struct Builder {
    struct Lease {
        PreparationWorkspace& workspace;
        ~Lease() {
            workspace.busy_.clear(std::memory_order_release);
        }
    };

    static core::Result<PreparationPlan> forecast(const PreparationRequest& request,
                                                  RationalLimits limits,
                                                  const std::stop_token& stop,
                                                  core::FailureDetails details) noexcept {
        CheckedCount count;
        const auto& grid = request.samples;
        const auto values = count.mul(grid.times.size(), grid.components);
        if (!count.ok)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        if (grid.times.size() < 2 || grid.components == 0 || grid.values.size() != values ||
            request.dimensions.size() != grid.components)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const Layout layout{grid.times.size(), grid.values.size(), request.events.size(),
                            which(request.boundary)};
        const auto fixed = layout.policy == Policy::Periodic ? 3U : 2U;
        const auto records =
            count.add(count.add(count.add(layout.times, values), layout.events), fixed);
        if (!count.ok || !fits(records))
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));

        Visits visits{limits.steps, 0, stop, details};
        std::uint64_t limbs = 0;
        for (std::uint64_t i = 0; i < records; ++i) {
            const auto value = raw_record(request, layout, static_cast<std::size_t>(i));
            const auto bound = shape(value, limits.input_bits, details);
            if (!bound)
                return std::unexpected(bound.error());
            const auto visit = visits.next();
            if (!visit)
                return std::unexpected(visit.error());
            limbs = count.add(limbs, count.add(value.numerator.size(), value.denominator.size()));
            if (!count.ok || !fits(limbs))
                return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        }

        Counter counter{limits, details};
        const auto read = [&](std::size_t index) { return raw_record(request, layout, index); };
        const auto inventory_visits = visits.used;
        Program program{counter, read, visits, limits};
        const auto valid = program.check(layout);
        if (!valid)
            return std::unexpected(valid.error());
        const auto arithmetic = counter.finish();
        if (!arithmetic)
            return std::unexpected(arithmetic.error());

        const auto copy = count.add(count.add(records, limbs), grid.components);
        const auto runtime = count.add(count.add(copy, visits.used - inventory_visits), 1);
        const auto work = count.add(count.add(visits.used, runtime), arithmetic->reserved_steps());
        const auto result_bytes = count.add(
            sizeof(Owner),
            count.add(count.mul(records, sizeof(RationalView)),
                      count.add(count.mul(limbs, sizeof(std::uint32_t)),
                                count.mul(grid.components, sizeof(numerics::SiDimension)))));
        const auto bytes = count.add(result_bytes, arithmetic->payload_bytes());
        if (!count.ok || work > limits.steps)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        PreparationPlan result;
        result.arithmetic_ = *arithmetic;
        result.records_ = records;
        result.limbs_ = limbs;
        result.bytes_ = bytes;
        result.work_ = work;
        result.scan_ = visits.used;
        return result;
    }

    static core::Result<PchipPreparationInput>
    make(const PreparationRequest& request, PreparationWorkspace& workspace, RationalLimits limits,
         core::ResourcePlan& resources, const std::stop_token& stop,
         core::FailureDetails details) noexcept {
        if (workspace.busy_.test_and_set(std::memory_order_acquire))
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const Lease lease{workspace};
        const auto plan = forecast(request, limits, stop, details);
        if (!plan)
            return std::unexpected(plan.error());
        return run(request, *plan, limits, resources, stop, details, plan->scan_);
    }

    static core::Result<PchipPreparationInput>
    run(const PreparationRequest& request, const PreparationPlan& admission, RationalLimits limits,
        core::ResourcePlan& resources, const std::stop_token& stop, core::FailureDetails details,
        std::uint64_t preflight = 0) noexcept {
        const auto* plan = &admission;
        CheckedCount count;
        const auto& grid = request.samples;
        const auto values = count.mul(grid.times.size(), grid.components);
        const auto fixed = which(request.boundary) == Policy::Periodic ? 3U : 2U;
        const auto records_count = count.add(
            count.add(count.add(grid.times.size(), values), request.events.size()), fixed);
        if (!count.ok)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        if (grid.times.size() < 2 || grid.components == 0 || grid.values.size() != values ||
            request.dimensions.size() != grid.components || records_count != plan->records_)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        if (plan->payload_bytes() > resources.remaining_bytes())
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));

        auto header = allocate<Owner>(resources, 1, details);
        if (!header)
            return std::unexpected(header.error());
        auto records = allocate<RationalView>(resources, plan->records_, details);
        if (!records)
            return std::unexpected(records.error());
        auto limbs = allocate<std::uint32_t>(resources, plan->limbs_, details);
        if (!limbs)
            return std::unexpected(limbs.error());
        auto dimensions =
            allocate<numerics::SiDimension>(resources, request.dimensions.size(), details);
        if (!dimensions)
            return std::unexpected(dimensions.error());
        auto arena = Workspace::create(plan->arithmetic_, resources, details, stop);
        if (!arena)
            return std::unexpected(arena.error());
        auto session = arena->begin(stop, details);
        if (!session)
            return std::unexpected(session.error());

        auto& owner = header->values()[0];
        owner.layout = {request.samples.times.size(), request.samples.values.size(),
                        request.events.size(), which(request.boundary)};
        owner.stamp = request.stamp;
        owner.details = details;
        owner.records = std::move(*records);
        owner.limbs = std::move(*limbs);
        owner.dimensions = std::move(*dimensions);
        const auto cap = plan->work_ - plan->arithmetic_.reserved_steps();
        Visits visits{cap, preflight, stop, details};
        const auto copied = copy(request, owner, visits);
        if (!copied)
            return std::unexpected(copied.error());
        const auto read = [&](std::size_t index) { return owner.records.values()[index]; };
        Program program{*session, read, visits, limits};
        const auto validated = program.check(owner.layout);
        if (!validated)
            return std::unexpected(validated.error());
        const auto usage = session->usage();
        if (!usage)
            return std::unexpected(usage.error());
        const auto complete = visits.next();
        if (!complete)
            return std::unexpected(complete.error());
        owner.work = visits.used + usage->steps;
        return PchipPreparationInput{std::move(*header)};
    }

    [[nodiscard]] static BoundaryRequest boundary(const Owner& owner) noexcept {
        const auto data = owner.records.values();
        const auto demand = owner.layout.demand();
        BoundaryRequest request{data[demand], data[demand + 1], RejectBoundary{}};
        if (owner.layout.policy == Policy::Hold)
            request.policy = HoldBoundary{};
        if (owner.layout.policy == Policy::Periodic)
            request.policy = PeriodicBoundary{data[owner.layout.period()]};
        return request;
    }

    [[nodiscard]] static core::Result<std::uint64_t>
    match(const PchipPreparationInput& actual, const PchipPreparationInput& expected,
          std::uint64_t limit, const std::stop_token& stop, core::FailureDetails details) noexcept {
        if (!actual.active() || !expected.active())
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const auto& a = actual.owner_.values()[0];
        const auto& b = expected.owner_.values()[0];
        Visits visits{limit, 0, stop, details};
        const auto started = visits.next();
        if (!started)
            return std::unexpected(started.error());
        const auto stale = failure(core::FailureCode::StaleGeneration, details);
        if (a.stamp != b.stamp || a.formula != b.formula || a.layout.policy != b.layout.policy ||
            a.layout.times != b.layout.times || a.layout.events != b.layout.events ||
            a.dimensions.size() != b.dimensions.size())
            return std::unexpected(stale);

        // Length mismatch terminates comparison before reading either sequence.
        CheckedCount count;
        const auto remaining = count.add(count.add(a.dimensions.size(), a.records.size()),
                                         count.add(a.limbs.size(), 1));
        if (!count.ok || remaining > limit - visits.used)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        for (std::size_t i = 0; i < a.dimensions.size(); ++i) {
            const auto next = visits.next();
            if (!next)
                return std::unexpected(next.error());
            if (a.dimensions.values()[i].powers != b.dimensions.values()[i].powers)
                return std::unexpected(stale);
        }
        for (std::size_t i = 0; i < a.records.size(); ++i) {
            const auto next = visits.next();
            if (!next)
                return std::unexpected(next.error());
            const auto l = a.records.values()[i];
            const auto r = b.records.values()[i];
            if (l.sign != r.sign)
                return std::unexpected(stale);
            const auto numerator = equal_words(l.numerator, r.numerator, visits);
            if (!numerator)
                return std::unexpected(numerator.error());
            const auto denominator = equal_words(l.denominator, r.denominator, visits);
            if (!denominator)
                return std::unexpected(denominator.error());
        }
        const auto done = visits.next();
        if (!done)
            return std::unexpected(done.error());
        return visits.used;
    }

private:
    [[nodiscard]] static bool fits(std::uint64_t count) noexcept {
        return count <= static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()) &&
               count <= std::numeric_limits<std::size_t>::max();
    }

    template <class T>
    [[nodiscard]] static core::Result<core::PlannedBuffer<T>>
    allocate(core::ResourcePlan& memory, std::uint64_t count,
             core::FailureDetails details) noexcept {
        auto result = memory.allocate<T>(static_cast<core::Index>(count), slot_detail::operation);
        if (!result)
            return std::unexpected(supplement(result.error(), details));
        return result;
    }

    [[nodiscard]] static core::Result<void> copy(const PreparationRequest& input, Owner& owner,
                                                 Visits& visits) noexcept {
        std::size_t cursor = 0;
        auto data = owner.limbs.values();
        const auto words = [&](std::span<const std::uint32_t> source)
            -> core::Result<std::span<const std::uint32_t>> {
            if (source.size() > data.size() - cursor)
                return std::unexpected(
                    failure(core::FailureCode::ResourceExceeded, visits.details));
            const auto target = data.subspan(cursor, source.size());
            for (std::size_t j = 0; j < source.size(); ++j) {
                const auto next = visits.next();
                if (!next)
                    return std::unexpected(next.error());
                target[j] = source[j];
            }
            cursor += source.size();
            return target;
        };
        auto records = owner.records.values();
        for (std::size_t i = 0; i < records.size(); ++i) {
            const auto next = visits.next();
            if (!next)
                return std::unexpected(next.error());
            const auto raw = raw_record(input, owner.layout, i);
            const auto n = words(raw.numerator);
            if (!n)
                return std::unexpected(n.error());
            const auto d = words(raw.denominator);
            if (!d)
                return std::unexpected(d.error());
            records[i] = {raw.sign, *n, *d};
        }
        for (std::size_t i = 0; i < input.dimensions.size(); ++i) {
            const auto next = visits.next();
            if (!next)
                return std::unexpected(next.error());
            owner.dimensions.values()[i] = input.dimensions[i];
        }
        return {};
    }

    [[nodiscard]] static core::Result<void> equal_words(std::span<const std::uint32_t> a,
                                                        std::span<const std::uint32_t> b,
                                                        Visits& visits) noexcept {
        const auto stale = failure(core::FailureCode::StaleGeneration, visits.details);
        if (a.size() != b.size())
            return std::unexpected(stale);
        for (std::size_t i = 0; i < a.size(); ++i) {
            const auto next = visits.next();
            if (!next)
                return std::unexpected(next.error());
            if (a[i] != b[i])
                return std::unexpected(stale);
        }
        return {};
    }
};
} // namespace preparation_detail

PchipPreparationInput::PchipPreparationInput(
    core::PlannedBuffer<preparation_detail::Owner> owner) noexcept
    : owner_(std::move(owner)) {}
PchipPreparationInput::PchipPreparationInput(PchipPreparationInput&&) noexcept = default;
PchipPreparationInput& PchipPreparationInput::operator=(PchipPreparationInput&&) noexcept = default;
PchipPreparationInput::~PchipPreparationInput() = default;

bool PchipPreparationInput::active() const noexcept {
    return owner_.size() != 0;
}

core::Result<PreparationRequest> PchipPreparationInput::view() const& noexcept {
    if (!active())
        return std::unexpected(preparation_detail::failure(core::FailureCode::InvalidInput));
    const auto& owner = owner_.values()[0];
    const auto records = owner.records.values();
    const auto& layout = owner.layout;
    const SampleGrid samples{records.first(layout.times),
                             records.subspan(layout.times, layout.values), owner.dimensions.size()};
    return PreparationRequest{samples, owner.dimensions.values(),
                              preparation_detail::Builder::boundary(owner),
                              records.subspan(layout.event_offset(), layout.events), owner.stamp};
}

core::Result<std::uint64_t> PchipPreparationInput::formula_version() const noexcept {
    if (!active())
        return std::unexpected(preparation_detail::failure(core::FailureCode::InvalidInput));
    return owner_.values()[0].formula;
}

std::uint64_t PchipPreparationInput::used_steps() const noexcept {
    return active() ? owner_.values()[0].work : 0;
}

core::Result<PreparationPlan> forecast_preparation_input(const PreparationRequest& request,
                                                         RationalLimits limits,
                                                         const std::stop_token& stop,
                                                         core::FailureDetails details) noexcept {
    return preparation_detail::Builder::forecast(request, limits, stop, details);
}

core::Result<PchipPreparationInput>
make_preparation_input(const PreparationRequest& request, PreparationWorkspace& workspace,
                       RationalLimits limits, core::ResourcePlan& resources,
                       const std::stop_token& stop, core::FailureDetails details) noexcept {
    return preparation_detail::Builder::make(request, workspace, limits, resources, stop, details);
}

core::Result<PchipPreparationInput>
make_planned_preparation_input(const PreparationRequest& request, const PreparationPlan& plan,
                               RationalLimits limits, core::ResourcePlan& memory,
                               const std::stop_token& stop, core::FailureDetails details) noexcept {
    return preparation_detail::Builder::run(request, plan, limits, memory, stop, details);
}

core::Result<std::uint64_t> match_preparation_input(const PchipPreparationInput& actual,
                                                    const PchipPreparationInput& expected,
                                                    std::uint64_t work_limit,
                                                    const std::stop_token& stop,
                                                    core::FailureDetails details) noexcept {
    return preparation_detail::Builder::match(actual, expected, work_limit, stop, details);
}
} // namespace uvlas::aero::detail::pchip
