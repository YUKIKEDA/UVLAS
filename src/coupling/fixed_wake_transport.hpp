#pragma once

#include "case_flow_clock.hpp"

namespace uvlas::coupling::detail {

/// One declared material record in a contact-qualified wake snapshot.
/// The release key and the current mesh point are distinct roles.
struct WakeMaterialPoint {
    aero::detail::MaterialReleaseKey release;
    aero::detail::PointId point;
    aero::detail::FacetId wake_facet;
};

struct WakeMaterialSnapshot {
    const aero::detail::GlobalContactState& contact;
    std::span<const WakeMaterialPoint> points; ///< Strictly ascending release keys.
};

struct WakeTransportClock {
    WakeTransportClock(const case_clock_detail::Grid& grid,
                       const case_clock_detail::Grid& current_grid,
                       CaseFlowClockQuery address) noexcept
        : clock(grid), current(current_grid), destination(address) {}

    const case_clock_detail::Grid& clock;
    const case_clock_detail::Grid& current;
    CaseFlowClockQuery destination;
};

struct WakeTransportLimits {
    CaseFlowClockLimits clock;
    std::uint64_t maximum_particles{}, maximum_steps{}, maximum_integer_bits{};
    core::Precision representation{};
};

struct WakeTransportRow {
    aero::detail::MaterialReleaseKey release;
    DiscreteLoadVector
        residual; ///< m/s, discrete transport rate minus the corresponding fluid velocity.
    std::optional<DiscreteLoadVector> direction;
    std::optional<DiscreteLoadVector> prescribed_position; ///< Explicit fixed-step candidate, m.
    std::optional<DiscreteLoadVector> prescribed_position_direction;
};

struct WakeTransportHistoryDirection {
    std::span<const NonpenetrationVector> previous, older; ///< XYZ, in release key order.
};

struct WakeTransportPositionDual {
    aero::detail::MaterialReleaseKey release;
    // Current, previous, older terms of the time stencil. The free velocity
    // transpose on current source geometry belongs to the same returned action.
    std::array<DiscreteLoadVector, 3> stencil;
};

struct WakeTransportData;

class WakeTransportResult {
public:
    WakeTransportResult(WakeTransportResult&&) noexcept;
    WakeTransportResult& operator=(WakeTransportResult&&) noexcept;
    WakeTransportResult(const WakeTransportResult&) = delete;
    WakeTransportResult& operator=(const WakeTransportResult&) = delete;
    ~WakeTransportResult();

    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] core::Result<WakeTransportRow> row(std::size_t) const noexcept;
    [[nodiscard]] core::Result<WakeTransportPositionDual> position_dual(std::size_t) const noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const DiscreteFreeWakeVelocity>>
    parent_free_velocity() const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const DiscreteFreeWakeVelocity>>
    parent_free_velocity() const&& = delete;
    [[nodiscard]] core::Result<aero::detail::ScaledEnclosure> prescribed_time_dual() const noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const DiscreteFreeWakeVelocity>>
    free_velocity() const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const DiscreteFreeWakeVelocity>>
    free_velocity() const&& = delete;

private:
    friend struct WakeTransportBuilder;
    explicit WakeTransportResult(core::PlannedBuffer<WakeTransportData>&&) noexcept;
    core::PlannedBuffer<WakeTransportData> data_;
};

/// Evaluate the old-particle residual at the same current contact-qualified
/// state as the actual velocity. Both past samples name physical release IDs.
[[nodiscard]] core::Result<WakeTransportResult> evaluate_fixed_free_wake_transport(
    const WakeTransportClock&, const std::array<WakeMaterialSnapshot, 3>&,
    const DiscreteNonpenetrationInput&, const DiscreteKnownFlowInput&,
    std::span<const aero::detail::PointId> particles, const DiscreteFreeWakeDirection*,
    const WakeTransportHistoryDirection&, std::span<const NonpenetrationVector> seeds,
    const WakeTransportLimits&, const DiscreteFreeWakeLimits&, aero::detail::SharedP1Workspace&,
    aero::detail::LocalPacketWorkspace&, aero::detail::P1Workspace&, core::ResourcePlan&,
    const std::stop_token& = {}) noexcept;

struct FirstFreeWakeEndpoint {
    const DiscreteNonpenetrationInput& source;
    const DiscreteKnownFlowInput& flow;
    std::span<const aero::detail::PointId> particles; ///< In material row order.
    const DiscreteFreeWakeDirection* direction{};
};

struct FirstFreeWakePoint {
    WakeMaterialPoint current;
    /// Empty at birth: the parent endpoint is the canonical trailing-edge point.
    std::optional<WakeMaterialPoint> parent;
};

/// Trapezoidal transport for newly released points, and for retained free
/// points on the first accepted interval step. Each fluid action and its
/// tangent/pullback are evaluated at that endpoint's own geometry and clock.
[[nodiscard]] core::Result<WakeTransportResult> evaluate_first_free_wake_transport(
    const WakeTransportClock&, const FirstFreeWakeEndpoint& current,
    const FirstFreeWakeEndpoint& parent, std::optional<aero::detail::ReleaseEpoch> birth_epoch,
    std::span<const FirstFreeWakePoint>, std::span<const NonpenetrationVector> seeds,
    const WakeTransportLimits&, const DiscreteFreeWakeLimits&, aero::detail::SharedP1Workspace&,
    aero::detail::LocalPacketWorkspace&, aero::detail::P1Workspace&, core::ResourcePlan&,
    const std::stop_token& = {}) noexcept;

/// Directions of the residual, including a common change to the evaluated
/// prescribed-flow time while the three-node stencil and h stay fixed.
struct PrescribedWakeTransportDirection {
    std::span<const NonpenetrationVector> current;
    WakeTransportHistoryDirection history;
    aero::detail::ScalarWords flow_time;
};
using PrescribedWakeClockSample =
    std::variant<std::reference_wrapper<const CaseClockedFlowSample>,
                 std::reference_wrapper<const CaseClockedKnownFlowSample>>;

[[nodiscard]] core::Result<WakeTransportResult> evaluate_fixed_prescribed_wake_transport(
    const WakeTransportClock&, const std::array<WakeMaterialSnapshot, 3>&,
    const PrescribedWakeClockSample&, const PrescribedWakeTransportDirection*,
    std::span<const NonpenetrationVector> seeds, const WakeTransportLimits&, core::ResourcePlan&,
    const std::stop_token& = {}) noexcept;

} // namespace uvlas::coupling::detail
