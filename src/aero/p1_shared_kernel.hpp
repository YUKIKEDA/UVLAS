#pragma once
#include "p1_induction.hpp"

namespace uvlas::aero::detail {
class SharedP1Action;
class SharedP1Workspace;
struct P1ActionSource;
template <class Result>
class SharedP1ReuseSlot;
namespace p1 {
struct Interval;
}
[[nodiscard]] bool p1_environment_supported() noexcept;
[[nodiscard]] p1::Interval decode_p1_word(double);
/// Retain the outward interval plus the final word rounding.
[[nodiscard]] ScaledEnclosure enclose_p1_value(const p1::Interval&, core::Precision,
                                               LocalRootScratch&);
/// Fill caller scratch with the real strength transpose for the three unit
/// velocity covectors, rounded exactly as the corresponding full adjoints.
/// A failed call leaves scratch unspecified. No result borrows this storage.
[[nodiscard]] AeroResult<std::uint64_t>
project_real_strength_transpose(const SharedP1Action&, const P1ActionSource&,
                                std::span<p1::Interval> strength_major_xyz, std::uint64_t steps,
                                SharedP1Workspace&, const std::stop_token&) noexcept;
/// Internal partition evaluation for a qualified shared-edge assembly. A disabled edge
/// removes its analytic contour term before numerical evaluation, retaining its sheet
/// term. Far evaluation already represents the full surface and is never masked.
[[nodiscard]] AeroResult<P1Basis>
evaluate_p1_partition(const LocalGeometryPacket&, const std::array<bool, 3>& boundary, P1Accuracy,
                      P1Limits, P1Workspace&, core::ResourcePlan&, const std::stop_token&,
                      SharedP1ReuseSlot<P1Basis>* = nullptr) noexcept;
/// Analytic P0 contour columns indexed by the three negative-oriented facet edges.
[[nodiscard]] AeroResult<P1Basis>
evaluate_bound_partition(const LocalGeometryPacket&, const std::array<bool, 3>& boundary,
                         P1Accuracy, P1Limits, P1Workspace&, core::ResourcePlan&,
                         const std::stop_token&, SharedP1ReuseSlot<P1Basis>* = nullptr) noexcept;
/// Near/bound payloads contain values; far payloads retain the fixed-rule qualification.
[[nodiscard]] AeroResult<P1Basis>
evaluate_partition_values(const LocalGeometryPacket&, const std::array<bool, 3>& boundary,
                          bool bound, P1Accuracy, P1Limits, P1Workspace&, core::ResourcePlan&,
                          const std::stop_token&, SharedP1ReuseSlot<P1Basis>*) noexcept;
} // namespace uvlas::aero::detail
