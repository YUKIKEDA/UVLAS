#pragma once
#include "p1_exact_product_private.hpp"
#include "p1_induction.hpp"
#include "p1_interval_private.hpp"

namespace uvlas::aero::detail {
struct LocalRootScratch;
namespace p1 {
/// Basis-major component-major jets: value followed by thirteen ambient derivatives.
using FarComponents = std::array<Interval, 126>;
/// Compare actual fixed-rule intervals. A rule difference is never a continuous error bound.
[[nodiscard]] P1QuadratureDifference qualify_far(std::span<const Interval, 126> low,
                                                 std::span<const Interval, 126> high,
                                                 const Interval& longest_edge, int local_scale,
                                                 Work&);
/// Exact physical full-width predicate; all supplied scalars are finite and nonnegative.
/// This internal predicate cannot construct a qualified result or replace its owned geometry.
[[nodiscard]] bool full_width_at_most(double radius, int output_scale, double tolerance,
                                      const ExactInteger& longest_squared, bool derivative,
                                      ExactProductScratch&, LocalRootScratch&);
} // namespace p1
} // namespace uvlas::aero::detail
