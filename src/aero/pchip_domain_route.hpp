#pragma once
#include "pchip_coefficients.hpp"
namespace uvlas::aero::detail::pchip {
/// Nonperiodic routing only; continuity and derivative-side selection remain separate.
enum class DomainPolicy : std::uint8_t { Reject, Hold };
/// Exact relation to the closed source domain.
enum class DomainRelation : std::uint8_t { Before, AtFirst, Inside, AtLast, After };
namespace domain_detail {
struct Builder;
}
/// Complete arithmetic and structural admission for two endpoint comparisons.
class DomainPlan {
public:
    /// Exact arena payload; excludes allocator, stack and process overhead.
    [[nodiscard]] std::uint64_t payload_bytes() const noexcept {
        return arithmetic_.payload_bytes();
    }
    /// Full call reservation, including endpoint shape checks.
    [[nodiscard]] std::uint64_t reserved_steps() const noexcept {
        return work_;
    }
    /// Arithmetic capacity for three imports and two comparisons.
    [[nodiscard]] const RationalPlan& arithmetic() const& noexcept {
        return arithmetic_;
    }
    [[nodiscard]] const RationalPlan& arithmetic() const&& = delete;

private:
    friend struct domain_detail::Builder;
    RationalPlan arithmetic_;
    std::uint64_t work_ = 0, outer_ = 0, preflight_ = 0;
};
/// Exclusive synchronous call lease; keep alive until the call returns.
class DomainWorkspace {
public:
    /// Construct an idle lease without allocating.
    DomainWorkspace() noexcept = default;
    DomainWorkspace(const DomainWorkspace&) = delete;
    DomainWorkspace& operator=(const DomainWorkspace&) = delete;

private:
    friend struct domain_detail::Builder;
    std::atomic_flag busy_ = ATOMIC_FLAG_INIT;
};
/// Owned scalar candidates without source references, identity or public curve qualification.
class DomainRoute {
public:
    /// Exact query relation; strict interior includes interior knots.
    [[nodiscard]] DomainRelation relation() const noexcept {
        return relation_;
    }
    /// Whether the original query may be passed to locate_time and IntervalRequest.
    [[nodiscard]] bool has_inside() const noexcept {
        return relation_ != DomainRelation::Before && relation_ != DomainRelation::After;
    }
    /// Whether an exterior constant candidate is available. Hold endpoints retain both candidates.
    [[nodiscard]] bool has_constant() const noexcept {
        return hold_ && relation_ != DomainRelation::Inside;
    }
    /// Endpoint index for CoefficientBlock::value; absent candidates return InvalidInput.
    [[nodiscard]] core::Result<std::uint64_t> constant_node() const noexcept;
    /// Actual arithmetic and structural visits for this call.
    [[nodiscard]] std::uint64_t used_steps() const noexcept {
        return work_;
    }
    /// Actual imported rational slots (query and both endpoints).
    [[nodiscard]] std::uint64_t used_slots() const noexcept {
        return slots_;
    }

private:
    friend struct domain_detail::Builder;
    DomainRoute(DomainRelation r, bool hold, std::uint64_t last, std::uint64_t work,
                std::uint64_t slots, core::FailureDetails details) noexcept
        : relation_(r), hold_(hold), last_(last), work_(work), slots_(slots), details_(details) {}
    DomainRelation relation_;
    bool hold_;
    std::uint64_t last_, work_, slots_;
    core::FailureDetails details_;
};
/// Count without allocating, inspecting only query and endpoint shapes. Unknown policy is
/// Unsupported.
[[nodiscard]] core::Result<DomainPlan> forecast_domain(const CoefficientBlock&, RationalView,
                                                       DomainPolicy, RationalLimits,
                                                       const std::stop_token& stop = {},
                                                       core::FailureDetails details = {}) noexcept;
/// Classify exact time; caller preserves immutable source/query borrows until return.
/// Reject outside is InvalidInput. No Automatic continuity or derivative choice is made here.
[[nodiscard]] core::Result<DomainRoute> route_domain(const CoefficientBlock&, RationalView,
                                                     DomainPolicy, DomainWorkspace&, RationalLimits,
                                                     core::ResourcePlan&,
                                                     const std::stop_token& stop = {},
                                                     core::FailureDetails details = {}) noexcept;

/// Admission for exact reduction by a supplied positive period.
class PeriodicTimePlan {
public:
    [[nodiscard]] std::uint64_t payload_bytes() const noexcept {
        return payload_;
    }
    [[nodiscard]] std::uint64_t reserved_steps() const noexcept {
        return work_;
    }
    [[nodiscard]] const RationalPlan& arithmetic() const& noexcept {
        return arithmetic_;
    }
    [[nodiscard]] const RationalPlan& arithmetic() const&& = delete;

private:
    friend struct domain_detail::Builder;
    RationalPlan arithmetic_;
    std::uint64_t payload_ = 0, work_ = 0, outer_ = 0, preflight_ = 0, limbs_ = 0;
};

/// Owned canonical physical time in [first,last), independent of the source and workspace.
class PeriodicTime {
public:
    [[nodiscard]] core::Result<RationalView> time() const& noexcept;
    [[nodiscard]] core::Result<RationalView> time() const&& = delete;
    /// At a seam, use the source's last time for a left interval and this time for a right
    /// interval.
    [[nodiscard]] bool is_seam() const noexcept {
        return limbs_.size() != 0 && seam_;
    }
    [[nodiscard]] std::uint64_t used_steps() const noexcept {
        return used_;
    }

private:
    friend struct domain_detail::Builder;
    PeriodicTime(core::PlannedBuffer<std::uint32_t> limbs, std::size_t numerator,
                 std::size_t denominator, int sign, bool seam, std::uint64_t used) noexcept
        : limbs_(std::move(limbs)),
          numerator_(numerator),
          denominator_(denominator),
          sign_(sign),
          seam_(seam),
          used_(used) {}

    core::PlannedBuffer<std::uint32_t> limbs_;
    std::size_t numerator_, denominator_;
    int sign_;
    bool seam_;
    std::uint64_t used_;
};

/// Require an explicit canonical period; an empty view is rejected as malformed input.
[[nodiscard]] core::Result<PeriodicTimePlan>
forecast_periodic_time(const CoefficientBlock&, RationalView query, RationalView period,
                       RationalLimits, const std::stop_token& stop = {},
                       core::FailureDetails details = {}) noexcept;

/// Check the period against the exact domain width and reduce the immutable query.
/// This private operation only handles time; seam value and slope qualification is separate.
[[nodiscard]] core::Result<PeriodicTime>
normalize_periodic_time(const CoefficientBlock&, RationalView query, RationalView period,
                        DomainWorkspace&, RationalLimits, core::ResourcePlan&,
                        const std::stop_token& stop = {},
                        core::FailureDetails details = {}) noexcept;

} // namespace uvlas::aero::detail::pchip
