#pragma once
#include "p1_induction.hpp"

namespace uvlas::aero::detail {
namespace p1 {
struct Interval;
}
[[nodiscard]] bool p1_environment_supported() noexcept;
[[nodiscard]] p1::Interval decode_p1_word(double);
/// Retain the outward interval plus the final word rounding.
[[nodiscard]] ScaledEnclosure enclose_p1_value(const p1::Interval&, core::Precision,
                                               LocalRootScratch&);
/// Internal partition evaluation for a qualified shared-edge assembly. A disabled edge
/// removes its analytic contour term before numerical evaluation, retaining its sheet
/// term. Far evaluation already represents the full surface and is never masked.
[[nodiscard]] AeroResult<P1Basis>
evaluate_p1_partition(const LocalGeometryPacket&, const std::array<bool, 3>& boundary, P1Accuracy,
                      P1Limits, P1Workspace&, core::ResourcePlan&, const std::stop_token&) noexcept;
} // namespace uvlas::aero::detail
