#include "pchip_side.hpp"

#include "pchip_slot_private.hpp"

#include <algorithm>
#include <array>

namespace uvlas::aero::detail::pchip {
namespace {

using slot_detail::CheckedCount;
using slot_detail::failure;
using slot_detail::supplement;

struct Meter {
    std::uint64_t limit = 0, used = 0;
    const std::stop_token& stop;
    core::FailureDetails details;

    [[nodiscard]] core::Result<void> visit() noexcept {
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));
        if (used == limit)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        ++used;
        return {};
    }
};

[[nodiscard]] core::Result<void> validate_jet(const IntervalValues& values, JetOrder expected,
                                              core::FailureDetails details) noexcept {
    const auto actual = values.requested_order();
    if (!actual)
        return std::unexpected(supplement(actual.error(), details));
    if (*actual != expected || values.components() == 0)
        return std::unexpected(failure(core::FailureCode::InvalidInput, details));
    return {};
}

[[nodiscard]] core::Result<std::array<RationalView, 2>>
second_derivatives(const IntervalValues& left, const IntervalValues& right, std::uint64_t component,
                   core::FailureDetails details) noexcept {
    const auto l = left.derivative(component, JetOrder::Second);
    if (!l)
        return std::unexpected(supplement(l.error(), details));
    const auto r = right.derivative(component, JetOrder::Second);
    if (!r)
        return std::unexpected(supplement(r.error(), details));
    return std::array{*l, *r};
}

[[nodiscard]] core::Result<void> equal_words(std::span<const std::uint32_t> left,
                                             std::span<const std::uint32_t> right,
                                             Meter& meter) noexcept {
    if (left.size() != right.size())
        return std::unexpected(failure(core::FailureCode::InvalidInput, meter.details));
    for (std::size_t i = 0; i < left.size(); ++i) {
        const auto visited = meter.visit();
        if (!visited)
            return std::unexpected(visited.error());
        if (left[i] != right[i])
            return std::unexpected(failure(core::FailureCode::InvalidInput, meter.details));
    }
    return {};
}

[[nodiscard]] core::Result<void> require_zero(const IntervalValues& values, std::uint64_t component,
                                              JetOrder order, Meter& meter) noexcept {
    const auto visited = meter.visit();
    if (!visited)
        return std::unexpected(visited.error());
    const auto derivative = values.derivative(component, order);
    if (!derivative)
        return std::unexpected(supplement(derivative.error(), meter.details));
    if (derivative->sign != 0)
        return std::unexpected(failure(core::FailureCode::InvalidInput, meter.details));
    return {};
}

} // namespace

core::Result<SideDecision> decide_side(SideLocation location, RequestedSide side, JetOrder order,
                                       core::FailureDetails details) noexcept {
    if ((side != RequestedSide::Automatic && side != RequestedSide::Left &&
         side != RequestedSide::Right) ||
        (order != JetOrder::Value && order != JetOrder::First && order != JetOrder::Second))
        return std::unexpected(failure(core::FailureCode::InvalidInput, details));
    switch (location) {
    case SideLocation::Single:
        return JetCandidate::Only;
    case SideLocation::C1Join:
    case SideLocation::HoldAtFirst:
    case SideLocation::HoldAtLast:
        break;
    default:
        return std::unexpected(failure(core::FailureCode::InvalidInput, details));
    }
    if (side == RequestedSide::Left)
        return JetCandidate::Left;
    if (side == RequestedSide::Right)
        return JetCandidate::Right;
    if (location == SideLocation::C1Join) {
        if (order == JetOrder::Second)
            return CompareSecondDerivatives{};
        return JetCandidate::Right;
    }
    const auto inside =
        location == SideLocation::HoldAtFirst ? JetCandidate::Right : JetCandidate::Left;
    if (order == JetOrder::Value)
        return inside;
    return CompareHeldDerivatives{inside, order};
}

core::Result<SideAgreement> agree_side(CompareSecondDerivatives, const IntervalValues& left,
                                       const IntervalValues& right, std::uint64_t work_limit,
                                       const std::stop_token& stop,
                                       core::FailureDetails details) noexcept {
    const auto checked_left = validate_jet(left, JetOrder::Second, details);
    if (!checked_left)
        return std::unexpected(checked_left.error());
    const auto checked_right = validate_jet(right, JetOrder::Second, details);
    if (!checked_right)
        return std::unexpected(checked_right.error());
    if (left.components() != right.components())
        return std::unexpected(failure(core::FailureCode::InvalidInput, details));

    Meter meter{work_limit, 0, stop, details};
    CheckedCount checked;
    std::uint64_t execution = 1; // Final cancellation and completion.
    for (std::uint64_t component = 0; component < left.components(); ++component) {
        const auto scanned = meter.visit();
        if (!scanned)
            return std::unexpected(scanned.error());
        const auto pair = second_derivatives(left, right, component, details);
        if (!pair)
            return std::unexpected(pair.error());
        const auto& [l, r] = *pair;
        const auto numerators = std::max(l.numerator.size(), r.numerator.size());
        const auto denominators = std::max(l.denominator.size(), r.denominator.size());
        // One record comparison and both complete limb sequences in the worst case.
        execution = checked.add(execution, checked.add(1, checked.add(numerators, denominators)));
        if (!checked.ok || execution > work_limit - meter.used)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
    }

    // Canonicality of IntervalValues makes equality a record and limb comparison.
    for (std::uint64_t component = 0; component < left.components(); ++component) {
        const auto visited = meter.visit();
        if (!visited)
            return std::unexpected(visited.error());
        const auto pair = second_derivatives(left, right, component, details);
        if (!pair)
            return std::unexpected(pair.error());
        const auto& [l, r] = *pair;
        if (l.sign != r.sign)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const auto equal_numerator = equal_words(l.numerator, r.numerator, meter);
        if (!equal_numerator)
            return std::unexpected(equal_numerator.error());
        const auto equal_denominator = equal_words(l.denominator, r.denominator, meter);
        if (!equal_denominator)
            return std::unexpected(equal_denominator.error());
    }
    const auto completed = meter.visit();
    if (!completed)
        return std::unexpected(completed.error());
    return SideAgreement{JetCandidate::Right, meter.used};
}

core::Result<SideAgreement> agree_side(CompareHeldDerivatives comparison,
                                       const IntervalValues& inside, std::uint64_t work_limit,
                                       const std::stop_token& stop,
                                       core::FailureDetails details) noexcept {
    if ((comparison.inside != JetCandidate::Left && comparison.inside != JetCandidate::Right) ||
        (comparison.through != JetOrder::First && comparison.through != JetOrder::Second))
        return std::unexpected(failure(core::FailureCode::InvalidInput, details));
    const auto valid = validate_jet(inside, comparison.through, details);
    if (!valid)
        return std::unexpected(valid.error());
    if (stop.stop_requested())
        return std::unexpected(failure(core::FailureCode::Cancelled, details));
    CheckedCount checked;
    const auto derivatives = static_cast<std::uint64_t>(comparison.through);
    const auto required = checked.add(checked.mul(inside.components(), derivatives), 1);
    if (!checked.ok || required > work_limit)
        return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));

    Meter meter{work_limit, 0, stop, details};
    for (std::uint64_t component = 0; component < inside.components(); ++component) {
        const auto first = require_zero(inside, component, JetOrder::First, meter);
        if (!first)
            return std::unexpected(first.error());
        if (comparison.through == JetOrder::Second) {
            const auto second = require_zero(inside, component, JetOrder::Second, meter);
            if (!second)
                return std::unexpected(second.error());
        }
    }
    const auto completed = meter.visit();
    if (!completed)
        return std::unexpected(completed.error());
    return SideAgreement{comparison.inside, meter.used};
}

} // namespace uvlas::aero::detail::pchip
