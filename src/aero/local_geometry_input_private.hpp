#pragma once

#include "exact_triangle_private.hpp"
#include "geometry_snapshot_storage.hpp"

namespace uvlas::aero::detail {
class LocalPacketEngine;
/// Borrowed only during the caller's genuine GeometryWorkspace lease.
struct LocalGeometryInput {
    std::shared_ptr<const GeometryData> owner;
    QueryStorage query;
    core::Index facet;
    std::reference_wrapper<ExactPoints> points;
    using Classification =
        std::variant<ExactTriangleResult, std::reference_wrapper<const TriangleOutcome>>;
    Classification classification;
    std::uint64_t preparation_steps;
    LocalGeometryInput(std::shared_ptr<const GeometryData> data, QueryStorage expression,
                       core::Index target, ExactPoints& coordinates, Classification outcome,
                       std::uint64_t steps) noexcept
        : owner(std::move(data)),
          query(expression),
          facet(target),
          points(coordinates),
          classification(std::move(outcome)),
          preparation_steps(steps) {}

    [[nodiscard]] const TriangleOutcome* outcome() const noexcept {
        if (const auto* owned = std::get_if<ExactTriangleResult>(&classification))
            return owned->outcome();
        const auto* borrowed =
            std::get_if<std::reference_wrapper<const TriangleOutcome>>(&classification);
        return borrowed != nullptr ? &borrowed->get() : nullptr;
    }
};
/// Sealed A1→A2 bridge: no input rounding, forged exact points, or second lease.
class GeometryPacketAccess {
    friend class LocalPacketEngine;
    [[nodiscard]] static OwnedGeometry retain(const LocalGeometryInput&) noexcept;
    [[nodiscard]] static AeroResult<LocalGeometryInput>
    load(const OwnedGeometry&, const ResolvedQuery&, FacetId, const GeometryLimits&,
         GeometryWorkspace&, const PredicateLease&, core::ResourcePlan&, std::stop_token,
         TriangleOutcome*) noexcept;
};
} // namespace uvlas::aero::detail
