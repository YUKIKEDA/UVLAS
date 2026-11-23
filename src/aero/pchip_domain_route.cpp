#include "pchip_domain_route.hpp"

#include "pchip_slot_private.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <limits>
namespace uvlas::aero::detail::pchip {
namespace domain_detail {
using slot_detail::CheckedCount;
using slot_detail::failure;
using slot_detail::supplement;
struct Meter {
    std::uint64_t cap = 0, used = 0;
    const std::stop_token& stop;
    core::FailureDetails details;
    core::Result<void> visit() {
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));
        if (used >= cap)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        ++used;
        return {};
    }
};
core::Result<BitBound> shape(RationalView value, RationalLimits limits,
                             core::FailureDetails details) {
    const auto words = limits.input_bits / 32 + (limits.input_bits % 32 != 0 ? 1U : 0U);
    if (value.numerator.size() > words || value.denominator.size() > words)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    if (value.denominator.empty() || value.denominator.back() == 0 ||
        (!value.numerator.empty() && value.numerator.back() == 0) ||
        (value.numerator.empty() ? value.sign != 0 : (value.sign != -1 && value.sign != 1)))
        return std::unexpected(failure(core::FailureCode::InvalidInput, details));
    CheckedCount checked;
    const auto bits = [&](std::span<const std::uint32_t> v) {
        return v.empty() ? std::uint64_t{0}
                         : checked.add(checked.mul(static_cast<std::uint64_t>(v.size() - 1), 32),
                                       static_cast<std::uint64_t>(std::bit_width(v.back())));
    };
    const BitBound result{bits(value.numerator), bits(value.denominator)};
    if (!checked.ok || std::max(result.numerator, result.denominator) > limits.input_bits)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    return result;
}

template <class Arithmetic, class Value>
core::Result<Value> reduce_time(Arithmetic& arithmetic, const Value& query, const Value& period,
                                const Value& first) noexcept {
    // Share the arithmetic graph between admission and numerical evaluation.
    return arithmetic.binary(RationalOp::Subtract, query, first)
        .and_then([&](const Value& elapsed) {
            return arithmetic.binary(RationalOp::Divide, elapsed, period);
        })
        .and_then([&](const Value& quotient) { return arithmetic.floor(quotient); })
        .and_then([&](const Value& cycles) {
            return arithmetic.binary(RationalOp::Multiply, cycles, period);
        })
        .and_then([&](const Value& shift) {
            return arithmetic.binary(RationalOp::Subtract, query, shift);
        });
}

struct Builder {
    struct Lease {
        DomainWorkspace& workspace;
        ~Lease() {
            workspace.busy_.clear(std::memory_order_release);
        }
    };
    static core::Result<DomainPlan> forecast(const CoefficientBlock& source, RationalView query,
                                             DomainPolicy policy, RationalLimits limits,
                                             const std::stop_token& stop,
                                             core::FailureDetails details) noexcept {
        if (policy != DomainPolicy::Reject && policy != DomainPolicy::Hold)
            return std::unexpected(failure(core::FailureCode::Unsupported, details));
        if (!source.active() || source.nodes() < 2)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        Meter scan{limits.steps, 0, stop, details};
        const auto first = source.time(0), last = source.time(source.nodes() - 1);
        if (!first)
            return std::unexpected(supplement(first.error(), details));
        if (!last)
            return std::unexpected(supplement(last.error(), details));
        Counter counter{limits, details};
        BitBound bounds[3];
        const RationalView inputs[]{query, *first, *last};
        for (std::size_t i = 0; i < 3; ++i) {
            if (const auto visited = scan.visit(); !visited)
                return std::unexpected(visited.error());
            const auto bound = shape(inputs[i], limits, details);
            if (!bound)
                return std::unexpected(bound.error());
            bounds[i] = *bound;
            const auto imported = counter.input(*bound);
            if (!imported)
                return std::unexpected(imported.error());
        }
        for (std::size_t i = 1; i < 3; ++i) {
            if (const auto visited = scan.visit(); !visited)
                return std::unexpected(visited.error());
            const auto compared = counter.compare(bounds[0], bounds[i]);
            if (!compared)
                return std::unexpected(compared.error());
        }
        const auto arithmetic = counter.finish();
        if (!arithmetic)
            return std::unexpected(arithmetic.error());
        DomainPlan plan;
        plan.arithmetic_ = *arithmetic;
        plan.preflight_ = scan.used;
        CheckedCount checked;
        // Three imports, two comparisons, and completion follow the five count visits.
        plan.outer_ = checked.add(plan.preflight_, 6);
        plan.work_ = checked.add(plan.outer_, arithmetic->reserved_steps());
        if (!checked.ok || plan.work_ > limits.steps)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        return plan;
    }
    static core::Result<DomainRoute> route(const CoefficientBlock& source, RationalView query,
                                           DomainPolicy policy, DomainWorkspace& workspace,
                                           RationalLimits limits, core::ResourcePlan& memory,
                                           const std::stop_token& stop,
                                           core::FailureDetails details) noexcept {
        if (workspace.busy_.test_and_set(std::memory_order_acquire))
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const Lease lease{workspace};
        const auto counted = forecast(source, query, policy, limits, stop, details);
        if (!counted)
            return std::unexpected(counted.error());
        auto arithmetic = Workspace::create(counted->arithmetic_, memory, details, stop);
        if (!arithmetic)
            return std::unexpected(arithmetic.error());
        return route_prepared(source, query, policy, *counted, *arithmetic, stop, details,
                              counted->preflight_);
    }

    static core::Result<DomainRoute>
    route_prepared(const CoefficientBlock& source, RationalView query, DomainPolicy policy,
                   const DomainPlan& plan, Workspace& arithmetic, const std::stop_token& stop,
                   core::FailureDetails details, std::uint64_t preflight = 0) noexcept {
        if (policy != DomainPolicy::Reject && policy != DomainPolicy::Hold)
            return std::unexpected(failure(core::FailureCode::Unsupported, details));
        if (!source.active() || source.nodes() < 2)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        auto numerical = arithmetic.begin(stop, details);
        if (!numerical)
            return std::unexpected(numerical.error());
        auto& session = *numerical;
        Meter meter{plan.outer_ - plan.preflight_ + preflight, preflight, stop, details};
        const auto first = source.time(0), last = source.time(source.nodes() - 1);
        if (!first)
            return std::unexpected(supplement(first.error(), details));
        if (!last)
            return std::unexpected(supplement(last.error(), details));
        const auto import = [&](RationalView value) -> core::Result<Slot> {
            if (const auto visited = meter.visit(); !visited)
                return std::unexpected(visited.error());
            return session.input(value);
        };
        const auto query_slot = import(query);
        if (!query_slot)
            return std::unexpected(query_slot.error());
        const auto first_slot = import(*first);
        if (!first_slot)
            return std::unexpected(first_slot.error());
        const auto last_slot = import(*last);
        if (!last_slot)
            return std::unexpected(last_slot.error());
        if (const auto visited = meter.visit(); !visited)
            return std::unexpected(visited.error());
        const auto at_first = session.compare(*query_slot, *first_slot);
        if (!at_first)
            return std::unexpected(at_first.error());
        if (const auto visited = meter.visit(); !visited)
            return std::unexpected(visited.error());
        const auto at_last = session.compare(*query_slot, *last_slot);
        if (!at_last)
            return std::unexpected(at_last.error());
        const auto relation = *at_first < 0    ? DomainRelation::Before
                              : *at_first == 0 ? DomainRelation::AtFirst
                              : *at_last > 0   ? DomainRelation::After
                              : *at_last == 0  ? DomainRelation::AtLast
                                               : DomainRelation::Inside;
        if (policy == DomainPolicy::Reject &&
            (relation == DomainRelation::Before || relation == DomainRelation::After))
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        if (const auto visited = meter.visit(); !visited)
            return std::unexpected(visited.error());
        const auto used = session.usage();
        if (!used)
            return std::unexpected(used.error());
        return DomainRoute{relation,           policy == DomainPolicy::Hold,
                           source.nodes() - 1, used->steps + meter.used,
                           used->slots,        details};
    }

    static core::Result<std::array<RationalView, 4>>
    periodic_inputs(const CoefficientBlock& source, RationalView query, RationalView period,
                    core::FailureDetails details) noexcept {
        if (!source.active() || source.nodes() < 2)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const auto first = source.time(0);
        const auto last = source.time(source.nodes() - 1);
        if (!first)
            return std::unexpected(supplement(first.error(), details));
        if (!last)
            return std::unexpected(supplement(last.error(), details));
        return std::array{query, period, *first, *last};
    }

    static core::Result<PeriodicTimePlan> forecast_periodic(const CoefficientBlock& source,
                                                            RationalView query, RationalView period,
                                                            RationalLimits limits,
                                                            const std::stop_token& stop,
                                                            core::FailureDetails details) noexcept {
        const auto inputs = periodic_inputs(source, query, period, details);
        if (!inputs)
            return std::unexpected(inputs.error());
        Meter scan{limits.steps, 0, stop, details};
        Counter counter{limits, details};
        std::array<BitBound, 4> bounds{};
        for (std::size_t i = 0; i < bounds.size(); ++i) {
            if (const auto visited = scan.visit(); !visited)
                return std::unexpected(visited.error());
            const auto bound = shape((*inputs)[i], limits, details);
            if (!bound)
                return std::unexpected(bound.error());
            const auto imported = counter.input(*bound);
            if (!imported)
                return std::unexpected(imported.error());
            bounds[i] = *imported;
        }
        const auto& [time, cycle, first, last] = bounds;
        const auto width = counter.binary(RationalOp::Subtract, last, first);
        if (!width)
            return std::unexpected(width.error());
        if (const auto checked_width = counter.compare(*width, cycle); !checked_width)
            return std::unexpected(checked_width.error());
        const auto normalized = reduce_time(counter, time, cycle, first);
        if (!normalized)
            return std::unexpected(normalized.error());
        if (const auto compare_first = counter.compare(*normalized, first); !compare_first)
            return std::unexpected(compare_first.error());
        if (const auto compare_last = counter.compare(*normalized, last); !compare_last)
            return std::unexpected(compare_last.error());
        const auto arithmetic = counter.finish();
        if (!arithmetic)
            return std::unexpected(arithmetic.error());

        PeriodicTimePlan plan;
        plan.normalized_ = *normalized;
        plan.arithmetic_ = *arithmetic;
        plan.preflight_ = scan.used;
        CheckedCount checked;
        plan.limbs_ = checked.mul(2, arithmetic->words());
        plan.payload_ = checked.add(arithmetic->payload_bytes(),
                                    checked.mul(plan.limbs_, sizeof(std::uint32_t)));
        // Four runtime imports, two bounded limb copies and one final cancellation check.
        plan.outer_ = checked.add(checked.add(plan.preflight_, 5), plan.limbs_);
        plan.work_ = checked.add(plan.outer_, arithmetic->reserved_steps());
        if (!checked.ok || plan.work_ > limits.steps ||
            plan.limbs_ > static_cast<std::uint64_t>(std::numeric_limits<core::Index>::max()) ||
            plan.limbs_ > std::numeric_limits<std::size_t>::max() / sizeof(std::uint32_t))
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        return plan;
    }

    static core::Result<PeriodicTime> periodic(const CoefficientBlock& source, RationalView query,
                                               RationalView period, DomainWorkspace& workspace,
                                               RationalLimits limits, core::ResourcePlan& memory,
                                               const std::stop_token& stop,
                                               core::FailureDetails details) noexcept {
        if (workspace.busy_.test_and_set(std::memory_order_acquire))
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const Lease lease{workspace};
        const auto counted = forecast_periodic(source, query, period, limits, stop, details);
        if (!counted)
            return std::unexpected(counted.error());
        auto storage = periodic_storage(*counted, memory, stop, details);
        if (!storage)
            return std::unexpected(storage.error());
        return periodic_prepared(source, query, period, *counted, *storage, stop, details,
                                 counted->preflight_);
    }

    static core::Result<PeriodicStorage> periodic_storage(const PeriodicTimePlan& plan,
                                                          core::ResourcePlan& memory,
                                                          const std::stop_token& stop,
                                                          core::FailureDetails details) noexcept {
        if (plan.payload_bytes() > memory.remaining_bytes())
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        auto arithmetic = Workspace::create(plan.arithmetic_, memory, details, stop);
        if (!arithmetic)
            return std::unexpected(arithmetic.error());
        auto limbs = memory.allocate<std::uint32_t>(static_cast<core::Index>(plan.limbs_),
                                                    slot_detail::operation);
        if (!limbs)
            return std::unexpected(supplement(limbs.error(), details));
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));
        return PeriodicStorage{std::move(*arithmetic), std::move(*limbs)};
    }

    static core::Result<PeriodicTime>
    periodic_prepared(const CoefficientBlock& source, RationalView query, RationalView period,
                      const PeriodicTimePlan& plan, PeriodicStorage& storage,
                      const std::stop_token& stop, core::FailureDetails details,
                      std::uint64_t preflight = 0) noexcept {
        if (storage.limbs_.size() != plan.limbs_ || plan.limbs_ == 0)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        auto numerical = storage.arena_.begin(stop, details);
        if (!numerical)
            return std::unexpected(numerical.error());
        auto& session = *numerical;
        Meter meter{plan.outer_ - plan.preflight_ + preflight, preflight, stop, details};
        const auto inputs = periodic_inputs(source, query, period, details);
        if (!inputs)
            return std::unexpected(inputs.error());
        const auto import = [&](RationalView value) -> core::Result<Slot> {
            if (const auto visited = meter.visit(); !visited)
                return std::unexpected(visited.error());
            return session.input(value);
        };
        const auto time = import((*inputs)[0]);
        if (!time)
            return std::unexpected(time.error());
        const auto cycle = import((*inputs)[1]);
        if (!cycle)
            return std::unexpected(cycle.error());
        const auto first = import((*inputs)[2]);
        if (!first)
            return std::unexpected(first.error());
        const auto last = import((*inputs)[3]);
        if (!last)
            return std::unexpected(last.error());
        const auto width = session.binary(RationalOp::Subtract, *last, *first);
        if (!width)
            return std::unexpected(width.error());
        const auto same_width = session.compare(*width, *cycle);
        if (!same_width)
            return std::unexpected(same_width.error());
        // CoefficientBlock admits only strictly increasing times. Equality proves positivity.
        if (*same_width != 0)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));

        const auto normalized = reduce_time(session, *time, *cycle, *first);
        if (!normalized)
            return std::unexpected(normalized.error());
        const auto compare_first = session.compare(*normalized, *first);
        if (!compare_first)
            return std::unexpected(compare_first.error());
        const auto compare_last = session.compare(*normalized, *last);
        if (!compare_last)
            return std::unexpected(compare_last.error());
        if (*compare_first < 0 || *compare_last >= 0)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const auto value = session.view(*normalized);
        if (!value)
            return std::unexpected(value.error());
        const auto words = static_cast<std::size_t>(plan.arithmetic_.words());
        auto output = storage.limbs_.values();
        if (value->numerator.size() > words || value->denominator.size() > words)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        const auto copy = [&](std::span<const std::uint32_t> from,
                              std::span<std::uint32_t> to) -> core::Result<void> {
            for (std::size_t i = 0; i < from.size(); ++i) {
                if (const auto visited = meter.visit(); !visited)
                    return std::unexpected(visited.error());
                to[i] = from[i];
            }
            return {};
        };
        if (const auto copied = copy(value->numerator, output.first(words)); !copied)
            return std::unexpected(copied.error());
        if (const auto copied = copy(value->denominator, output.subspan(words)); !copied)
            return std::unexpected(copied.error());
        const auto used = session.usage();
        if (!used)
            return std::unexpected(used.error());
        if (const auto visited = meter.visit(); !visited)
            return std::unexpected(visited.error());
        return PeriodicTime{std::move(storage.limbs_), value->numerator.size(),
                            value->denominator.size(), value->sign,
                            *compare_first == 0,       used->steps + meter.used};
    }
};
} // namespace domain_detail
core::Result<std::uint64_t> DomainRoute::constant_node() const noexcept {
    if (!has_constant())
        return std::unexpected(slot_detail::failure(core::FailureCode::InvalidInput, details_));
    return relation_ == DomainRelation::Before || relation_ == DomainRelation::AtFirst
               ? std::uint64_t{0}
               : last_;
}
core::Result<DomainPlan> forecast_domain(const CoefficientBlock& source, RationalView query,
                                         DomainPolicy policy, RationalLimits limits,
                                         const std::stop_token& stop,
                                         core::FailureDetails details) noexcept {
    return domain_detail::Builder::forecast(source, query, policy, limits, stop, details);
}
core::Result<Workspace> prepare_domain_storage(const DomainPlan& plan, core::ResourcePlan& memory,
                                               const std::stop_token& stop,
                                               core::FailureDetails details) noexcept {
    return Workspace::create(plan.arithmetic(), memory, details, stop);
}

core::Result<DomainRoute> route_prepared_domain(const CoefficientBlock& source, RationalView query,
                                                DomainPolicy policy, const DomainPlan& plan,
                                                Workspace&& storage, const std::stop_token& stop,
                                                core::FailureDetails details) noexcept {
    return domain_detail::Builder::route_prepared(source, query, policy, plan, storage, stop,
                                                  details);
}

core::Result<DomainRoute> route_domain(const CoefficientBlock& source, RationalView query,
                                       DomainPolicy policy, DomainWorkspace& workspace,
                                       RationalLimits limits, core::ResourcePlan& memory,
                                       const std::stop_token& stop,
                                       core::FailureDetails details) noexcept {
    return domain_detail::Builder::route(source, query, policy, workspace, limits, memory, stop,
                                         details);
}

core::Result<RationalView> PeriodicTime::time() const& noexcept {
    if (limbs_.size() == 0)
        return std::unexpected(slot_detail::failure(core::FailureCode::InvalidInput));
    const auto data = limbs_.values();
    return RationalView{sign_, data.first(numerator_), data.subspan(data.size() / 2, denominator_)};
}
core::Result<PeriodicTimePlan> forecast_periodic_time(const CoefficientBlock& source,
                                                      RationalView query, RationalView period,
                                                      RationalLimits limits,
                                                      const std::stop_token& stop,
                                                      core::FailureDetails details) noexcept {
    return domain_detail::Builder::forecast_periodic(source, query, period, limits, stop, details);
}
core::Result<PeriodicStorage> prepare_periodic_storage(const PeriodicTimePlan& plan,
                                                       core::ResourcePlan& memory,
                                                       const std::stop_token& stop,
                                                       core::FailureDetails details) noexcept {
    return domain_detail::Builder::periodic_storage(plan, memory, stop, details);
}

core::Result<PeriodicTime> normalize_prepared_periodic_time(const CoefficientBlock& source,
                                                            RationalView query, RationalView period,
                                                            const PeriodicTimePlan& plan,
                                                            PeriodicStorage&& storage,
                                                            const std::stop_token& stop,
                                                            core::FailureDetails details) noexcept {
    return domain_detail::Builder::periodic_prepared(source, query, period, plan, storage, stop,
                                                     details);
}

core::Result<PeriodicTime>
normalize_periodic_time(const CoefficientBlock& source, RationalView query, RationalView period,
                        DomainWorkspace& workspace, RationalLimits limits,
                        core::ResourcePlan& memory, const std::stop_token& stop,
                        core::FailureDetails details) noexcept {
    return domain_detail::Builder::periodic(source, query, period, workspace, limits, memory, stop,
                                            details);
}
} // namespace uvlas::aero::detail::pchip
