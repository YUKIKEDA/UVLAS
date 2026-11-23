#pragma once
#include "pchip_rational_slots.hpp"

#include <atomic>
#include <optional>
namespace uvlas::aero::detail::pchip {
/// Borrow canonical row-major samples for one synchronous mathematical preparation.
class SharedPreparationInput;

struct SampleGrid {
    std::span<const RationalView> times;  ///< Exactly n>=2 strictly increasing exact times.
    std::span<const RationalView> values; ///< n*m exact values, node then component order.
    std::uint64_t components =
        0; ///< Positive component count; no physical dimension qualification.
};
namespace coefficient_detail {
struct Builder;
struct Owner;
} // namespace coefficient_detail
/// Checked forecast for one complete private preparation, not a qualified time curve.
class CoefficientPlan {
public:
    /// Exact sum of all planned heap object/array payload, excluding allocator/OS overhead.
    [[nodiscard]] std::uint64_t payload_bytes() const noexcept {
        return payload_;
    }
    /// Combined numeric and outer scan/copy/loop visit reservation.
    [[nodiscard]] std::uint64_t reserved_steps() const noexcept {
        return work_;
    }
    /// Structural visits already performed by a separate coefficient forecast.
    [[nodiscard]] std::uint64_t preflight_steps() const noexcept {
        return preflight_;
    }
    /// Rational slot/bit/work forecast used by the shared exact executor.
    [[nodiscard]] const RationalPlan& arithmetic() const& noexcept {
        return arithmetic_;
    }
    [[nodiscard]] const RationalPlan& arithmetic() const&& = delete;

private:
    friend struct coefficient_detail::Builder;
    RationalPlan arithmetic_;
    std::uint64_t n_ = 0, m_ = 0, inputs_ = 0, records_ = 0, payload_ = 0, work_ = 0, outer_ = 0,
                  limbs_ = 0, preflight_ = 0;
};
/// Storage reserved before an integrated preparation enters exact arithmetic.
class CoefficientStorage {
public:
    CoefficientStorage(CoefficientStorage&&) noexcept;
    CoefficientStorage& operator=(CoefficientStorage&&) noexcept;
    ~CoefficientStorage();

private:
    friend struct coefficient_detail::Builder;
    CoefficientStorage(Workspace, core::PlannedBuffer<std::optional<Slot>>,
                       core::PlannedBuffer<coefficient_detail::Owner>) noexcept;
    Workspace arithmetic_;
    core::PlannedBuffer<std::optional<Slot>> handles_;
    core::PlannedBuffer<coefficient_detail::Owner> owner_;
};

/// Exclusive call lease. Preparation allocates its entire planned arena once before arithmetic.
/// Destroy only after the synchronous call returns; distinct workspaces may run concurrently.
class CoefficientWorkspace {
public:
    CoefficientWorkspace() noexcept = default; ///< Construct an idle lease without allocating.
    CoefficientWorkspace(const CoefficientWorkspace&) = delete;
    CoefficientWorkspace& operator=(const CoefficientWorkspace&) = delete;

private:
    friend struct coefficient_detail::Builder;
    std::atomic_flag busy_ = ATOMIC_FLAG_INIT;
};
/// Private owned mathematical coefficients; no outside-policy/C1/issuer/physical qualification.
class CoefficientBlock {
public:
    /// Transfer the exact input and coefficient owner; source becomes inactive.
    CoefficientBlock(CoefficientBlock&&) noexcept;
    /// Release earlier owned storage and transfer the source.
    CoefficientBlock& operator=(CoefficientBlock&&) noexcept;
    /// Destroy the owned canonical records independently of scratch/input/ResourcePlan lifetime.
    ~CoefficientBlock();
    CoefficientBlock(const CoefficientBlock&) = delete;
    CoefficientBlock& operator=(const CoefficientBlock&) = delete;
    /// Whether this wrapper still owns its prepared records.
    [[nodiscard]] bool active() const noexcept;
    /// Number of original sample nodes, or zero for an inactive owner.
    [[nodiscard]] std::uint64_t nodes() const noexcept;
    /// Number of ordered mathematical components, or zero for an inactive owner.
    [[nodiscard]] std::uint64_t components() const noexcept;
    /// Fixed coefficient formula identifier; does not qualify a physical provider.
    [[nodiscard]] static constexpr std::uint64_t formula_version() noexcept {
        return 1;
    }
    /// Borrow one exact original time; invalid index/inactive owner returns InvalidInput.
    [[nodiscard]] core::Result<RationalView> time(std::uint64_t) const& noexcept;
    [[nodiscard]] core::Result<RationalView> time(std::uint64_t) const&& = delete;
    /// Borrow an original sample in node/component order.
    [[nodiscard]] core::Result<RationalView> value(std::uint64_t, std::uint64_t) const& noexcept;
    [[nodiscard]] core::Result<RationalView> value(std::uint64_t, std::uint64_t) const&& = delete;
    /// Borrow one ascending-power u coefficient in interval/component/power order.
    [[nodiscard]] core::Result<RationalView> coefficient(std::uint64_t, std::uint64_t,
                                                         std::uint64_t) const& noexcept;
    [[nodiscard]] core::Result<RationalView> coefficient(std::uint64_t, std::uint64_t,
                                                         std::uint64_t) const&& = delete;
    /// Read charged numeric and outer visits after successful preparation.
    [[nodiscard]] std::uint64_t used_steps() const noexcept;
    /// Read the number of live rational SSA values materialized during preparation.
    [[nodiscard]] std::uint64_t used_slots() const noexcept;

private:
    friend struct coefficient_detail::Builder;
    explicit CoefficientBlock(core::PlannedBuffer<coefficient_detail::Owner>) noexcept;
    core::PlannedBuffer<coefficient_detail::Owner> owner_;
};
/// Count all coefficient branches without allocating or accepting noncanonical math as qualified.
/// GCD canonicality is checked after admission by prepare_coefficients.
[[nodiscard]] core::Result<CoefficientPlan>
forecast_coefficients(SampleGrid, RationalLimits, const std::stop_token& stop = {},
                      core::FailureDetails details = {}, bool retain_original = false) noexcept;

[[nodiscard]] core::Result<CoefficientStorage>
prepare_coefficient_storage(const CoefficientPlan&, core::ResourcePlan&,
                            const std::stop_token& stop = {},
                            core::FailureDetails details = {}) noexcept;

/// Run a plan counted for these inputs in retained-sample mode.
[[nodiscard]] core::Result<CoefficientBlock>
prepare_staged_coefficients(const SharedPreparationInput&, const CoefficientPlan&,
                            CoefficientStorage&&, RationalLimits, const std::stop_token& stop = {},
                            core::FailureDetails details = {}) noexcept;

/// Admit all scratch plus owned-result payload before allocating. Every numeric loop is prepared.
/// Caller keeps the sample borrow immutable until return. Failure leaves it and old blocks intact.
[[nodiscard]] core::Result<CoefficientBlock>
prepare_coefficients(SampleGrid, CoefficientWorkspace&, RationalLimits, core::ResourcePlan&,
                     const std::stop_token& stop = {}, core::FailureDetails details = {}) noexcept;
} // namespace uvlas::aero::detail::pchip
