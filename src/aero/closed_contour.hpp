#pragma once

#include "p1_shared_action.hpp"

#include <atomic>

namespace uvlas::aero::detail {

/// Ordinary SI endpoints for a certified Binary64 interval.
struct ContourInterval {
    double lower{};
    double upper{};
};
using ContourVector = std::array<ContourInterval, 3>;

struct ClosedContourData;
struct ClosedContourEndpoint;

/// Reusable per-thread endpoint storage. Failed evaluations leave it reusable.
class ClosedContourWorkspace {
public:
    ClosedContourWorkspace(ClosedContourWorkspace&&) noexcept;
    ClosedContourWorkspace& operator=(ClosedContourWorkspace&&) noexcept;
    ~ClosedContourWorkspace();

private:
    friend class ClosedContour;
    explicit ClosedContourWorkspace(core::PlannedBuffer<ClosedContourEndpoint>&&) noexcept;
    core::PlannedBuffer<ClosedContourEndpoint> endpoints_;
    std::atomic_flag entered_ = ATOMIC_FLAG_INIT;
};

/// Fixed, complete constant-strength contour on an owned current snapshot.
/// Use a separate workspace and output row for each concurrent caller.
class ClosedContour {
public:
    ClosedContour(ClosedContour&&) noexcept;
    ClosedContour& operator=(ClosedContour&&) noexcept;
    ~ClosedContour();

    [[nodiscard]] static AeroResult<ClosedContour>
    prepare(const GlobalContactState&, std::span<const P1FaceStrengths>,
            std::span<const P1StrengthId>, GlobalLimits, GlobalWorkspace&, core::ResourcePlan&,
            std::stop_token = {}) noexcept;

    [[nodiscard]] AeroResult<ClosedContourWorkspace> workspace(core::ResourcePlan&) const noexcept;
    [[nodiscard]] std::uint64_t preparation_steps() const noexcept;

    /// Project the induced velocity of each strength column onto a supplied normal enclosure.
    /// Entries have units of inverse length. Discard the row on failure.
    [[nodiscard]] AeroResult<std::uint64_t>
    project(const ResolvedQuery&, const ContourVector& normal, std::span<ContourInterval> row,
            const GeometryLimits&, GlobalLimits, GeometryWorkspace&, GlobalWorkspace&,
            ClosedContourWorkspace&, std::stop_token = {}) const noexcept;

    /// Enclose the induced velocity at a free point or a physical load-edge midpoint.
    /// Strengths are in native column order; midpoint evaluation omits its own bound edge.
    [[nodiscard]] AeroResult<std::pair<ContourVector, std::uint64_t>>
    apply(const ResolvedQuery&, std::optional<EdgeId> self, std::span<const ContourInterval>,
          const GeometryLimits&, GlobalLimits, GeometryWorkspace&, GlobalWorkspace&,
          ClosedContourWorkspace&, std::stop_token = {}) const noexcept;

private:
    ClosedContour(core::PlannedBuffer<ClosedContourData>&&, GeometryContext) noexcept;
    core::PlannedBuffer<ClosedContourData> data_;
    GeometryContext context_;
};

} // namespace uvlas::aero::detail
