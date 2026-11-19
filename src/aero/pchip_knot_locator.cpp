#include "pchip_knot_locator.hpp"

#include "pchip_slot_private.hpp"

#include <algorithm>
#include <bit>
namespace uvlas::aero::detail::pchip {
namespace locator_detail {
using slot_detail::CheckedCount;
using slot_detail::failure;
using slot_detail::supplement;
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
core::Result<BitBound> shape(RationalView value, RationalLimits limits,
                             core::FailureDetails details) {
    const auto maximum_words = limits.input_bits / 32 + (limits.input_bits % 32 != 0 ? 1U : 0U);
    if (value.numerator.size() > maximum_words || value.denominator.size() > maximum_words)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    if (value.denominator.empty() || value.denominator.back() == 0 ||
        (!value.numerator.empty() && value.numerator.back() == 0) ||
        (value.numerator.empty() ? value.sign != 0 : (value.sign != -1 && value.sign != 1)))
        return std::unexpected(failure(core::FailureCode::InvalidInput, details));
    CheckedCount checked;
    const auto bits = [&](std::span<const std::uint32_t> words) {
        return words.empty()
                   ? std::uint64_t{0}
                   : checked.add(checked.mul(static_cast<std::uint64_t>(words.size() - 1), 32),
                                 static_cast<std::uint64_t>(std::bit_width(words.back())));
    };
    const BitBound result{bits(value.numerator), bits(value.denominator)};
    if (!checked.ok || std::max(result.numerator, result.denominator) > limits.input_bits)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    return result;
}
struct Builder {
    static core::Result<LocatorPlan> forecast(const CoefficientBlock& source, RationalView query,
                                              RationalLimits limits, const std::stop_token& stop,
                                              core::FailureDetails details) noexcept {
        if (!source.active() || source.nodes() < 2)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        LocatorPlan plan;
        plan.nodes_ = source.nodes();
        CheckedCount checked;
        // For nonempty n, bit_width(n) equals ceil(log2(n+1)), without forming n+1.
        plan.comparisons_ = checked.add(static_cast<std::uint64_t>(std::bit_width(plan.nodes_)), 3);
        Meter scan{limits.steps, 0, stop, details};
        if (const auto visited = scan.visit(); !visited)
            return std::unexpected(visited.error());
        const auto query_shape = shape(query, limits, details);
        if (!query_shape)
            return std::unexpected(query_shape.error());
        BitBound maximum{0, 1};
        for (std::uint64_t i = 0; i < plan.nodes_; ++i) {
            if (const auto visited = scan.visit(); !visited)
                return std::unexpected(visited.error());
            const auto time = source.time(i);
            if (!time)
                return std::unexpected(supplement(time.error(), details));
            const auto bound = shape(*time, limits, details);
            if (!bound)
                return std::unexpected(bound.error());
            maximum.numerator = std::max(maximum.numerator, bound->numerator);
            maximum.denominator = std::max(maximum.denominator, bound->denominator);
        }
        Counter counter{limits, details};
        const auto imported_query = counter.input(*query_shape);
        if (!imported_query)
            return std::unexpected(imported_query.error());
        for (std::uint64_t i = 0; i < plan.comparisons_; ++i) {
            if (const auto visited = scan.visit(); !visited)
                return std::unexpected(visited.error());
            const auto imported = counter.input(maximum);
            if (!imported)
                return std::unexpected(imported.error());
            const auto comparison = counter.compare(*query_shape, maximum);
            if (!comparison)
                return std::unexpected(comparison.error());
        }
        const auto arithmetic = counter.finish();
        if (!arithmetic)
            return std::unexpected(arithmetic.error());
        plan.arithmetic_ = *arithmetic;
        plan.preflight_ = scan.used;
        // One query import, one outer visit per probe, and one scalar-result completion.
        plan.outer_ = checked.add(plan.preflight_, checked.add(plan.comparisons_, 2));
        plan.work_ = checked.add(plan.arithmetic_.reserved_steps(), plan.outer_);
        if (!checked.ok || plan.work_ > limits.steps)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        return plan;
    }
    static core::Result<KnotLocation> locate(const CoefficientBlock& source, RationalView query,
                                             LocatorWorkspace& workspace, RationalLimits limits,
                                             core::ResourcePlan& memory,
                                             const std::stop_token& stop,
                                             core::FailureDetails details) noexcept {
        if (workspace.busy_.test_and_set(std::memory_order_acquire))
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        struct Lease {
            LocatorWorkspace& workspace;
            ~Lease() {
                workspace.busy_.clear(std::memory_order_release);
            }
        };
        const Lease lease{workspace};
        const auto counted = forecast(source, query, limits, stop, details);
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
        if (const auto visited = meter.visit(); !visited)
            return std::unexpected(visited.error());
        const auto imported = session.input(query);
        if (!imported)
            return std::unexpected(imported.error());
        const auto query_slot = *imported;
        std::uint64_t comparisons = 0;
        const auto probe = [&](std::uint64_t index) -> core::Result<std::strong_ordering> {
            if (comparisons >= plan.comparisons_)
                return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
            if (const auto visited = meter.visit(); !visited)
                return std::unexpected(visited.error());
            const auto time = source.time(index);
            if (!time)
                return std::unexpected(supplement(time.error(), details));
            const auto node = session.input(*time);
            if (!node)
                return std::unexpected(node.error());
            const auto order = session.compare(query_slot, *node);
            if (!order)
                return std::unexpected(order.error());
            ++comparisons;
            return *order;
        };
        const auto first = probe(0);
        if (!first)
            return std::unexpected(first.error());
        const auto last = probe(plan.nodes_ - 1);
        if (!last)
            return std::unexpected(last.error());
        if (*first < 0 || *last > 0)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        std::uint64_t low = 0, high = plan.nodes_;
        while (low < high) {
            const auto mid = low + (high - low) / 2;
            const auto order = probe(mid);
            if (!order)
                return std::unexpected(order.error());
            if (*order <= 0)
                high = mid;
            else
                low = mid + 1;
        }
        if (low >= plan.nodes_)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const auto equality = probe(low);
        if (!equality)
            return std::unexpected(equality.error());
        const bool exact = *equality == 0;
        if (*equality > 0 || (!exact && low == 0))
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        if (const auto visited = meter.visit(); !visited)
            return std::unexpected(visited.error());
        const auto used = session.usage();
        if (!used)
            return std::unexpected(used.error());
        return KnotLocation{plan.nodes_, low,    exact, comparisons, used->steps + meter.used,
                            used->slots, details};
    }
};
} // namespace locator_detail
core::Result<std::uint64_t> KnotLocation::knot_index() const noexcept {
    if (!exact_)
        return std::unexpected(slot_detail::failure(core::FailureCode::InvalidInput, details_));
    return position_;
}
core::Result<std::uint64_t> KnotLocation::left_interval() const noexcept {
    if (!has_left())
        return std::unexpected(slot_detail::failure(core::FailureCode::InvalidInput, details_));
    return position_ - 1;
}
core::Result<std::uint64_t> KnotLocation::right_interval() const noexcept {
    if (!has_right())
        return std::unexpected(slot_detail::failure(core::FailureCode::InvalidInput, details_));
    return exact_ ? position_ : position_ - 1;
}
core::Result<std::uint64_t> KnotLocation::choose(KnotSide side) const noexcept {
    if (side == KnotSide::Left)
        return left_interval();
    if (side == KnotSide::Right)
        return right_interval();
    return std::unexpected(slot_detail::failure(core::FailureCode::InvalidInput, details_));
}
core::Result<std::uint64_t> KnotLocation::unique_interval() const noexcept {
    if (!exact_)
        return position_ - 1;
    if (position_ == 0)
        return std::uint64_t{0};
    if (position_ == nodes_ - 1)
        return position_ - 1;
    return std::unexpected(slot_detail::failure(core::FailureCode::InvalidInput, details_));
}
core::Result<LocatorPlan> forecast_location(const CoefficientBlock& source, RationalView query,
                                            RationalLimits limits, const std::stop_token& stop,
                                            core::FailureDetails details) noexcept {
    return locator_detail::Builder::forecast(source, query, limits, stop, details);
}
core::Result<KnotLocation> locate_time(const CoefficientBlock& source, RationalView query,
                                       LocatorWorkspace& workspace, RationalLimits limits,
                                       core::ResourcePlan& memory, const std::stop_token& stop,
                                       core::FailureDetails details) noexcept {
    return locator_detail::Builder::locate(source, query, workspace, limits, memory, stop, details);
}
} // namespace uvlas::aero::detail::pchip
