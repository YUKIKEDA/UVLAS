#pragma once
#include "pchip_coefficients.hpp"
namespace uvlas::aero::detail::pchip {
/// Closed request for the mathematical derivatives with respect to physical input time.
enum class JetOrder : std::uint8_t { Value = 0, First = 1, Second = 2 };
/// Explicitly selected interval, without any outside policy, search or automatic side selection.
struct IntervalRequest {
    std::uint64_t interval = 0; ///< Interval in the immutable coefficient block.
    RationalView
        time; ///< Canonical exact time, synchronously borrowed and checked against both ends.
    JetOrder order = JetOrder::Value; ///< Includes value and all derivatives through this order.
};
namespace interval_detail {
struct Builder;
struct Owner;
} // namespace interval_detail
/// Immutable payload/work forecast; no curve or physical qualification is conferred.
class IntervalPlan {
public:
    /// Planned scratch plus owned-result heap object/array payload; excludes allocator overhead.
    [[nodiscard]] std::uint64_t payload_bytes() const noexcept {
        return payload_;
    }
    /// Numeric primitive and outer scan/import/copy visits for the entire call.
    [[nodiscard]] std::uint64_t reserved_steps() const noexcept {
        return work_;
    }
    [[nodiscard]] std::uint64_t preflight_steps() const noexcept {
        return preflight_;
    }
    [[nodiscard]] std::uint64_t local_preflight_steps() const noexcept {
        return preflight_ - all_intervals_scan_;
    }
    /// Exact slot forecast used by the shared prepared rational executor.
    [[nodiscard]] const RationalPlan& arithmetic() const& noexcept {
        return arithmetic_;
    }
    [[nodiscard]] const RationalPlan& arithmetic() const&& = delete;

private:
    friend struct interval_detail::Builder;
    RationalPlan arithmetic_;
    std::uint64_t m_ = 0, inputs_ = 0, records_ = 0, limbs_ = 0, payload_ = 0, work_ = 0,
                  outer_ = 0, preflight_ = 0, fields_ = 0, all_intervals_scan_ = 0;
};
/// Scratch and result storage reserved before curve routing and exact arithmetic.
class IntervalStorage {
public:
    IntervalStorage(IntervalStorage&&) noexcept;
    IntervalStorage& operator=(IntervalStorage&&) noexcept;
    ~IntervalStorage();

private:
    friend struct interval_detail::Builder;
    IntervalStorage(Workspace, core::PlannedBuffer<std::optional<Slot>>,
                    core::PlannedBuffer<interval_detail::Owner>) noexcept;
    Workspace arithmetic_;
    core::PlannedBuffer<std::optional<Slot>> handles_;
    core::PlannedBuffer<interval_detail::Owner> owner_;
};

/// Exclusive synchronous call lease. All arena allocation precedes numeric execution.
class IntervalWorkspace {
public:
    /// Construct an idle lease without allocation; keep it alive until the call returns.
    IntervalWorkspace() noexcept = default;
    IntervalWorkspace(const IntervalWorkspace&) = delete;
    IntervalWorkspace& operator=(const IntervalWorkspace&) = delete;

private:
    friend struct interval_detail::Builder;
    std::atomic_flag busy_ = ATOMIC_FLAG_INIT;
};
/// Owned private interval mathematics only, never ExactPchipJet/CurveReady/KnownFlow.
/// These values do not certify input generation, C1/side/outside policy or a parent curve.
class IntervalValues {
public:
    /// Transfer the canonical owner; source becomes inactive.
    IntervalValues(IntervalValues&&) noexcept;
    /// Release earlier storage and transfer ownership.
    IntervalValues& operator=(IntervalValues&&) noexcept;
    /// Release values independently of input/coefficient/workspace/ResourcePlan lifetime.
    ~IntervalValues();
    IntervalValues(const IntervalValues&) = delete;
    IntervalValues& operator=(const IntervalValues&) = delete;
    /// Whether this wrapper still owns its values.
    [[nodiscard]] bool active() const noexcept;
    /// Number of ordered components, or zero for an inactive wrapper.
    [[nodiscard]] std::uint64_t components() const noexcept;
    /// Return the requested derivative order, or InvalidInput for an inactive wrapper.
    [[nodiscard]] core::Result<JetOrder> requested_order() const noexcept;
    /// Return the explicitly selected interval; no search or side decision is implied.
    [[nodiscard]] core::Result<std::uint64_t> interval() const noexcept;
    /// Borrow the original exact query time in physical time units.
    [[nodiscard]] core::Result<RationalView> time() const& noexcept;
    [[nodiscard]] core::Result<RationalView> time() const&& = delete;
    /// Borrow u=(time-t_i)/(t_{i+1}-t_i). This is not a policy-normalized physical time.
    [[nodiscard]] core::Result<RationalView> coordinate_u() const& noexcept;
    [[nodiscard]] core::Result<RationalView> coordinate_u() const&& = delete;
    /// Borrow value/order1/order2 for a component. Unrequested orders return InvalidInput.
    [[nodiscard]] core::Result<RationalView> derivative(std::uint64_t, JetOrder) const& noexcept;
    [[nodiscard]] core::Result<RationalView> derivative(std::uint64_t, JetOrder) const&& = delete;
    /// Recorded numeric and structural visits, including preflight scans.
    [[nodiscard]] std::uint64_t used_steps() const noexcept;
    /// Number of rational SSA values materialized during this call.
    [[nodiscard]] std::uint64_t used_slots() const noexcept;
    /// Fixed mathematical expression version, without a physical-provider qualification.
    [[nodiscard]] static constexpr std::uint64_t formula_version() noexcept {
        return 1;
    }

private:
    friend struct interval_detail::Builder;
    explicit IntervalValues(core::PlannedBuffer<interval_detail::Owner>) noexcept;
    core::PlannedBuffer<interval_detail::Owner> owner_;
};
/// Count the exact selected-interval expression without allocating. Canonical GCD follows
/// admission. The existing slot input cap covers imported coefficients/query and fixed integer
/// literals.
[[nodiscard]] core::Result<IntervalPlan>
forecast_interval(const CoefficientBlock&, IntervalRequest, RationalLimits,
                  const std::stop_token& stop = {}, core::FailureDetails details = {}) noexcept;
/// Count the same expression before the interval and its exact query have been selected.
[[nodiscard]] core::Result<IntervalPlan>
forecast_interval_bound(const CoefficientBlock&, BitBound query, JetOrder, RationalLimits,
                        const std::stop_token& stop = {},
                        core::FailureDetails details = {}) noexcept;
[[nodiscard]] core::Result<IntervalStorage>
prepare_interval_storage(const IntervalPlan&, core::ResourcePlan&, const std::stop_token& stop = {},
                         core::FailureDetails details = {}) noexcept;

/// The plan and storage cover the selected interval, query and order.
[[nodiscard]] core::Result<IntervalValues>
evaluate_staged_interval(const CoefficientBlock&, IntervalRequest, const IntervalPlan&,
                         IntervalStorage&&, RationalLimits, const std::stop_token& stop = {},
                         core::FailureDetails details = {}) noexcept;

/// Recount/admit all arrays, then evaluate on a prepared arena without allocation in numeric loops.
/// Reject queries outside the selected closed interval. Keep source/query borrows immutable to
/// return.
[[nodiscard]] core::Result<IntervalValues>
evaluate_interval(const CoefficientBlock&, IntervalRequest, IntervalWorkspace&, RationalLimits,
                  core::ResourcePlan&, const std::stop_token& stop = {},
                  core::FailureDetails details = {}) noexcept;
} // namespace uvlas::aero::detail::pchip
