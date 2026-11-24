#pragma once

#include "geometry_snapshot.hpp"

namespace uvlas::aero::detail {

/// Exact dyadic point for a qualified query and the fixed-chart material pullback.
/// Free points have no material stencil. Stencil weights have denominator 2^1074.
struct ExactQueryCoordinates {
    std::array<ExactInteger, 3> numerator;
    std::array<ExactInteger, 3> weights;
    std::array<std::uint64_t, 3> point_ids{};
    std::size_t point_count = 0;
    static constexpr unsigned weight_shift = 1074;
    static constexpr unsigned denominator_shift = 2 * weight_shift;
    std::uint64_t used_steps = 0;
};

[[nodiscard]] AeroResult<ExactQueryCoordinates>
exact_query_coordinates(const OwnedGeometry&, const ResolvedQuery&, const GeometryLimits&,
                        GeometryWorkspace&, std::stop_token = {},
                        Configuration = Configuration::Current) noexcept;

} // namespace uvlas::aero::detail
