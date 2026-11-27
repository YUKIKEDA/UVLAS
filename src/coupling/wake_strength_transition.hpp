#pragma once

#include "fixed_wake_transport.hpp"

namespace uvlas::coupling::detail {

struct WakeStrengthEndpoint {
    const DiscreteNonpenetrationInput& source;
    std::span<const WakeMaterialPoint> material;
};

/// Oriented span endpoints identify a physical band independent of numeric source columns.
struct WakeTraceMaterial {
    aero::detail::PointId from{0}, to{0};
    aero::detail::ReleaseEpoch epoch{0};
};

struct RetiredWakeTrace {
    WakeTraceMaterial material;
    FrozenWakeTrace trace;
};

/// One signed entry per trial history column in the material transfer.
struct WakeStrengthParent {
    aero::detail::P1StrengthId id;
    bool reversed{};
};

/// Own material correspondence and the exact frozen words qualified for the trial.
class WakeStrengthTransition {
public:
    WakeStrengthTransition(WakeStrengthTransition&&) noexcept = default;
    WakeStrengthTransition& operator=(WakeStrengthTransition&&) noexcept = default;

    [[nodiscard]] const WakeMaterialTransition& material() const& noexcept {
        return material_;
    }
    [[nodiscard]] const WakeMaterialTransition& material() const&& = delete;
    [[nodiscard]] std::span<const FrozenWakeTrace> history() const& noexcept {
        return history_.values();
    }
    [[nodiscard]] std::span<const FrozenWakeTrace> history() const&& = delete;
    /// Physical identities in history() order, even when the trial renumbers columns.
    [[nodiscard]] std::span<const WakeTraceMaterial> history_material() const& noexcept {
        return identities_.values();
    }
    [[nodiscard]] std::span<const WakeTraceMaterial> history_material() const&& = delete;
    [[nodiscard]] std::span<const RetiredWakeTrace> retired_history() const& noexcept {
        return retired_.values();
    }
    [[nodiscard]] std::span<const RetiredWakeTrace> retired_history() const&& = delete;
    [[nodiscard]] std::span<const WakeStrengthParent> parents() const& noexcept {
        return parents_.values();
    }
    [[nodiscard]] std::span<const WakeStrengthParent> parents() const&& = delete;

private:
    friend core::Result<WakeStrengthTransition>
    prepare_wake_strength_transition(const io::SiCase&, const WakeTransportClock&,
                                     const WakeStrengthEndpoint&, const WakeStrengthEndpoint&,
                                     aero::detail::ReleaseEpoch, const WakeTransportLimits&,
                                     core::ResourcePlan&, const std::stop_token&) noexcept;
    WakeStrengthTransition(WakeMaterialTransition&& material,
                           core::PlannedBuffer<FrozenWakeTrace>&& history,
                           core::PlannedBuffer<WakeTraceMaterial>&& identities,
                           core::PlannedBuffer<RetiredWakeTrace>&& retired,
                           core::PlannedBuffer<WakeStrengthParent>&& parents) noexcept
        : material_(std::move(material)),
          history_(std::move(history)),
          identities_(std::move(identities)),
          retired_(std::move(retired)),
          parents_(std::move(parents)) {}
    WakeMaterialTransition material_;
    core::PlannedBuffer<FrozenWakeTrace> history_;
    core::PlannedBuffer<WakeTraceMaterial> identities_;
    core::PlannedBuffer<RetiredWakeTrace> retired_;
    core::PlannedBuffer<WakeStrengthParent> parents_;
};

/// Derive each strength's material band and release epoch from its checked D2
/// quadrilateral. Newly detached strength comes from the parent's signed T Gamma.
/// Retained strength copies the parent words, and the retired material alone may
/// lose its trace. The current attached coefficient stays implicit.
[[nodiscard]] core::Result<WakeStrengthTransition> prepare_wake_strength_transition(
    const io::SiCase&, const WakeTransportClock&, const WakeStrengthEndpoint& parent,
    const WakeStrengthEndpoint& trial, aero::detail::ReleaseEpoch, const WakeTransportLimits&,
    core::ResourcePlan&, const std::stop_token& = {}) noexcept;

/// Missing parent keys have zero variation. Supplied IDs are in source order
/// and may include bound circulation IDs as well as accepted wake strengths.
[[nodiscard]] core::Result<core::PlannedBuffer<FrozenWakeTrace>> apply_wake_strength_transfer(
    const WakeStrengthTransition&, std::span<const FrozenWakeTrace> parent_direction,
    std::uint64_t maximum_steps, core::ResourcePlan&, const std::stop_token& = {}) noexcept;

class WakeStrengthPullback {
public:
    [[nodiscard]] std::size_t size() const noexcept {
        return rows_.size();
    }
    [[nodiscard]] core::Result<FreeWakeStrengthCovector> row(std::size_t) const noexcept;

private:
    friend core::Result<WakeStrengthPullback> pullback_wake_strength_transfer(
        const WakeStrengthTransition&, std::span<const FreeWakeStrengthCovector>, core::Precision,
        std::uint64_t, core::ResourcePlan&, const std::stop_token&) noexcept;
    explicit WakeStrengthPullback(
        core::PlannedBuffer<std::optional<FreeWakeStrengthCovector>>&& rows) noexcept
        : rows_(std::move(rows)) {}

    core::PlannedBuffer<std::optional<FreeWakeStrengthCovector>> rows_;
};

/// Seeds use the exact trial history IDs and order. Return one enclosure per
/// contributing parent strength, summing any repeated parent IDs before rounding.
[[nodiscard]] core::Result<WakeStrengthPullback> pullback_wake_strength_transfer(
    const WakeStrengthTransition&, std::span<const FreeWakeStrengthCovector> seeds, core::Precision,
    std::uint64_t maximum_steps, core::ResourcePlan&, const std::stop_token& = {}) noexcept;

struct FreeWakeTrialEndpoint {
    WakeStrengthEndpoint wake;
    const DiscreteKnownFlowInput& flow;
    // For the trial endpoint, history must be empty and is derived from the parent.
    const DiscreteFreeWakeDirection* direction{};
};

/// Only BDF2 rows use older positions or the explicit two history directions.
struct FreeWakeTrialHistory {
    const WakeMaterialSnapshot* older{};
    WakeTransportHistoryDirection retained_direction;
};

struct FreeWakeTrialLimits {
    WakeTransportLimits material;
    WakeTransportLimits transport; ///< Each independent stencil evaluation.
    DiscreteFreeWakeLimits velocity;
};

/// Both parts use the same checked trial geometry, circulation and released keys.
struct FreeWakeTrial {
    WakeStrengthTransition strength;
    std::optional<WakeTransportResult> trapezoid;
    std::optional<WakeTransportResult> bdf2;
};

/// New releases use the trapezoid; retained points do so only at an interval
/// restart. Seeds use trial material order. Each returned part retains key order.
[[nodiscard]] core::Result<FreeWakeTrial> evaluate_free_wake_trial(
    const io::SiCase&, const WakeTransportClock&, const FreeWakeTrialEndpoint& parent,
    const FreeWakeTrialEndpoint& trial, const FreeWakeTrialHistory&, aero::detail::ReleaseEpoch,
    std::span<const NonpenetrationVector> seeds, const FreeWakeTrialLimits&,
    aero::detail::SharedP1Workspace&, aero::detail::LocalPacketWorkspace&,
    aero::detail::P1Workspace&, core::ResourcePlan&, const std::stop_token& = {}) noexcept;

} // namespace uvlas::coupling::detail
