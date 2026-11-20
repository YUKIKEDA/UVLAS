#pragma once

#include "pchip_interval_jet.hpp"

#include <variant>

namespace uvlas::aero::detail::pchip {

enum class RequestedSide : std::uint8_t { Automatic, Left, Right };

/// Classification established by the calling route on one immutable curve.
/// Single covers interval interiors, Reject endpoints, and exterior Hold constants.
/// C1Join covers an interior knot or a qualified periodic seam.
enum class SideLocation : std::uint8_t { Single, C1Join, HoldAtFirst, HoldAtLast };

enum class JetCandidate : std::uint8_t { Only, Left, Right };

struct CompareSecondDerivatives {};
struct CompareHeldDerivatives {
    JetCandidate inside;
    JetOrder through;
};

using SideDecision = std::variant<JetCandidate, CompareSecondDerivatives, CompareHeldDerivatives>;

/// Determine which candidates a caller needs. A direct choice never requires the opposite jet.
[[nodiscard]] core::Result<SideDecision> decide_side(SideLocation, RequestedSide, JetOrder,
                                                     core::FailureDetails details = {}) noexcept;

struct SideAgreement {
    JetCandidate selected;
    std::uint64_t used_steps;
};

/// Compare all second derivatives in the caller's left and right canonical candidates.
/// The caller establishes their common curve and exact query (including normalization at a seam).
/// Agreement selects the right interval. This operation has no allocation or physical provider.
[[nodiscard]] core::Result<SideAgreement>
agree_side(CompareSecondDerivatives, const IntervalValues& left, const IntervalValues& right,
           std::uint64_t work_limit, const std::stop_token& stop = {},
           core::FailureDetails details = {}) noexcept;

/// Compare the requested derivatives in the inside jet against the canonical zero of the
/// exterior Hold constant. The caller establishes that its value comes from the same endpoint.
/// Agreement selects the inside interval.
[[nodiscard]] core::Result<SideAgreement>
agree_side(CompareHeldDerivatives, const IntervalValues& inside, std::uint64_t work_limit,
           const std::stop_token& stop = {}, core::FailureDetails details = {}) noexcept;

} // namespace uvlas::aero::detail::pchip
