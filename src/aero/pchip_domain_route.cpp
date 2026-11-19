#include "pchip_domain_route.hpp"

#include "pchip_slot_private.hpp"

#include <algorithm>
#include <bit>
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
struct Builder {
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
        struct Lease {
            DomainWorkspace& workspace;
            ~Lease() {
                workspace.busy_.clear(std::memory_order_release);
            }
        };
        const Lease lease{workspace};
        const auto counted = forecast(source, query, policy, limits, stop, details);
        if (!counted)
            return std::unexpected(counted.error());
        const auto& plan = *counted;
        if (plan.payload_bytes() > memory.remaining_bytes())
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        auto arithmetic = Workspace::create(plan.arithmetic_, memory, details, stop);
        if (!arithmetic)
            return std::unexpected(arithmetic.error());
        auto numerical = arithmetic->begin(stop, details);
        if (!numerical)
            return std::unexpected(numerical.error());
        auto& session = *numerical;
        Meter meter{plan.outer_, plan.preflight_, stop, details};
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
core::Result<DomainRoute> route_domain(const CoefficientBlock& source, RationalView query,
                                       DomainPolicy policy, DomainWorkspace& workspace,
                                       RationalLimits limits, core::ResourcePlan& memory,
                                       const std::stop_token& stop,
                                       core::FailureDetails details) noexcept {
    return domain_detail::Builder::route(source, query, policy, workspace, limits, memory, stop,
                                         details);
}
} // namespace uvlas::aero::detail::pchip
