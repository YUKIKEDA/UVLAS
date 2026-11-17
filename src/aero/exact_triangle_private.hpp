#pragma once

#include "exact_triangle_predicates.hpp"

#include <array>
#include <functional>

namespace uvlas::aero::detail {
class GeometryEngine;
/// Fixed q=2^-2148 coordinate slots, accessible only to the validated factories.
using ExactPoints = std::array<std::array<ExactInteger, 3>, 4>;
class ExactPredicateAccess {
    friend class GeometryEngine;
    friend core::Result<ExactTriangleResult> classify_triangle(CoordinateView, PredicateWorkspace&,
                                                               core::ResourcePlan&, PredicateLimits,
                                                               std::stop_token,
                                                               core::FailureDetails) noexcept;
    /// The existing lease must belong to this workspace; no nested acquire occurs.
    [[nodiscard]] static core::Result<std::reference_wrapper<ExactPoints>>
    coordinates(PredicateWorkspace&, const PredicateLease&, core::FailureDetails) noexcept;
    [[nodiscard]] static core::Result<bool>
    nondegenerate(PredicateWorkspace&, const PredicateLease&, core::FailureDetails) noexcept;
    /// All coordinate values use fixed q and the proved 3173-bit magnitude bound.
    [[nodiscard]] static core::Result<ExactTriangleResult>
    finish(PredicateWorkspace&, const PredicateLease&, core::PlannedBuffer<TriangleOutcome>,
           PredicateLimits, std::stop_token, core::FailureDetails, std::uint64_t) noexcept;
};
} // namespace uvlas::aero::detail
