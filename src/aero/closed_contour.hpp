#pragma once

#include "p1_shared_action.hpp"

#include <atomic>
#include <cstdint>

namespace uvlas::aero::detail {

/// Ordinary SI endpoints for a certified Binary64 interval.
struct ContourInterval {
    double lower{};
    double upper{};
};
using ContourVector = std::array<ContourInterval, 3>;

/// One already enclosed age-quadrature sample for the fixed straight wake.
struct StraightWakeAgeSample {
    ContourInterval age_seconds;
    ContourInterval real_weight_seconds;      ///< q (cos(omega age) - 1).
    ContourInterval imaginary_weight_seconds; ///< -q sin(omega age).
};

struct StraightWakeColumnUpdate {
    ContourInterval real;
    ContourInterval imaginary;
};

struct StraightWakeAgeQuadrature {
    core::PlannedBuffer<StraightWakeAgeSample> samples;
    std::uint64_t interval_steps{};
};

/// Own the finite Gauss16 rule on panels with strictly increasing starts from
/// zero. The final endpoint may be enclosed when the physical age is not dyadic.
/// Frequency is in rad/s. Error covers nodes, weights and phase of the finite rule.
[[nodiscard]] AeroResult<StraightWakeAgeQuadrature> prepare_straight_wake_age_quadrature(
    const GeometryContext&, std::span<const double> panel_starts_seconds,
    ContourInterval retained_end_seconds, double angular_frequency, std::uint64_t maximum_samples,
    std::uint64_t maximum_interval_steps, core::ResourcePlan&, std::stop_token = {}) noexcept;

/// Project the age derivative of a directed, constant-strength wake contour at the
/// supplied samples. The contour goes from first to second, downstream, then back
/// to first. Core radius is strictly positive. The caller qualifies the geometry
/// and quadrature; these intervals enclose the supplied finite sum in 1/m.
[[nodiscard]] AeroResult<StraightWakeColumnUpdate> project_straight_wake_column_update(
    const GeometryContext&, const ContourVector& first, const ContourVector& second,
    const ContourVector& convection_velocity, ContourInterval core_radius,
    const ContourVector& point, const ContourVector& normal, std::span<const StraightWakeAgeSample>,
    std::uint64_t maximum_samples, std::stop_token = {}) noexcept;

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

/// Different interval evaluations of the same finite-edge cross product.
enum class ContourCrossProduct : std::uint8_t { EndpointRelative, MidpointRelative };

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
            ClosedContourWorkspace&, std::stop_token = {},
            ContourCrossProduct = ContourCrossProduct::EndpointRelative) const noexcept;

    /// Enclose the induced velocity at a free point or a physical load-edge midpoint.
    /// Strengths are in native column order; midpoint evaluation omits its own bound edge.
    [[nodiscard]] AeroResult<std::pair<ContourVector, std::uint64_t>>
    apply(const ResolvedQuery&, std::optional<EdgeId> self, std::span<const ContourInterval>,
          const GeometryLimits&, GlobalLimits, GeometryWorkspace&, GlobalWorkspace&,
          ClosedContourWorkspace&, std::stop_token = {},
          ContourCrossProduct = ContourCrossProduct::EndpointRelative) const noexcept;

private:
    ClosedContour(core::PlannedBuffer<ClosedContourData>&&, GeometryContext) noexcept;
    core::PlannedBuffer<ClosedContourData> data_;
    GeometryContext context_;
};

} // namespace uvlas::aero::detail
