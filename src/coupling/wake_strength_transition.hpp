#pragma once

#include "fixed_wake_transport.hpp"

namespace uvlas::coupling::detail {

struct WakeStrengthEndpoint {
    const DiscreteNonpenetrationInput& source;
    std::span<const WakeMaterialPoint> material;
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

private:
    friend core::Result<WakeStrengthTransition>
    prepare_wake_strength_transition(const io::SiCase&, const WakeTransportClock&,
                                     const WakeStrengthEndpoint&, const WakeStrengthEndpoint&,
                                     aero::detail::ReleaseEpoch, const WakeTransportLimits&,
                                     core::ResourcePlan&, const std::stop_token&) noexcept;
    WakeStrengthTransition(WakeMaterialTransition&& material,
                           core::PlannedBuffer<FrozenWakeTrace>&& history) noexcept
        : material_(std::move(material)), history_(std::move(history)) {}
    WakeMaterialTransition material_;
    core::PlannedBuffer<FrozenWakeTrace> history_;
};

/// Derive each strength's material band and release epoch from its checked D2
/// quadrilateral. Newly detached strength comes from the parent's signed T Gamma.
/// Retained strength copies the parent words, and the retired material alone may
/// lose its trace. The current attached coefficient stays implicit.
[[nodiscard]] core::Result<WakeStrengthTransition> prepare_wake_strength_transition(
    const io::SiCase&, const WakeTransportClock&, const WakeStrengthEndpoint& parent,
    const WakeStrengthEndpoint& trial, aero::detail::ReleaseEpoch, const WakeTransportLimits&,
    core::ResourcePlan&, const std::stop_token& = {}) noexcept;

} // namespace uvlas::coupling::detail
