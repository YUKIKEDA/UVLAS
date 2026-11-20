#pragma once
#include "pchip_coefficients.hpp"

#include <variant>

namespace uvlas::aero::detail::pchip {

struct RejectBoundary {};
struct HoldBoundary {};
struct PeriodicBoundary {
    RationalView period;
};

/// Closed policy and exact, inclusive physical interval requiring C1.
struct BoundaryRequest {
    RationalView begin, end;
    std::variant<RejectBoundary, HoldBoundary, PeriodicBoundary> policy;
};

namespace boundary_detail {
struct Builder;
}

/// Checked reservation for one complete private boundary check.
class BoundaryPlan {
public:
    [[nodiscard]] std::uint64_t payload_bytes() const noexcept {
        return arithmetic_.payload_bytes();
    }
    [[nodiscard]] std::uint64_t reserved_steps() const noexcept {
        return work_;
    }
    [[nodiscard]] const RationalPlan& arithmetic() const& noexcept {
        return arithmetic_;
    }
    [[nodiscard]] const RationalPlan& arithmetic() const&& = delete;

private:
    friend struct boundary_detail::Builder;
    RationalPlan arithmetic_;
    std::uint64_t work_ = 0, scan_ = 0;
};

/// Exclusive synchronous lease. No payload is acquired until the whole call is admitted.
class BoundaryWorkspace {
public:
    BoundaryWorkspace() noexcept = default;
    BoundaryWorkspace(const BoundaryWorkspace&) = delete;
    BoundaryWorkspace& operator=(const BoundaryWorkspace&) = delete;

private:
    friend struct boundary_detail::Builder;
    std::atomic_flag busy_ = ATOMIC_FLAG_INIT;
};

/// Scalar accounting for the boundary check. No curve, identity, or provider is constructed.
struct BoundaryUsage {
    std::uint64_t steps, slots;
};

/// Admit every possible boundary visited by this policy without allocating.
[[nodiscard]] core::Result<BoundaryPlan>
forecast_boundary(const CoefficientBlock&, BoundaryRequest, RationalLimits,
                  const std::stop_token& stop = {}, core::FailureDetails details = {}) noexcept;

/// Check Reject containment, Hold endpoint slopes on the demanded interval, or periodic seam C1.
/// Source and all argument borrows remain immutable until return.
[[nodiscard]] core::Result<BoundaryUsage>
check_boundary(const CoefficientBlock&, BoundaryRequest, BoundaryWorkspace&, RationalLimits,
               core::ResourcePlan&, const std::stop_token& stop = {},
               core::FailureDetails details = {}) noexcept;

} // namespace uvlas::aero::detail::pchip
