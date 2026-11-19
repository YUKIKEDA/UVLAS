#pragma once
#include "pchip_coefficients.hpp"
namespace uvlas::aero::detail::pchip {
/// Explicit in-domain candidate only; no outside policy or automatic continuity decision.
enum class KnotSide : std::uint8_t { Left, Right };
namespace locator_detail {
struct Builder;
}
/// Bounded search forecast; the returned scalar location requires no heap storage.
class LocatorPlan {
public:
    /// Prepared arithmetic arena payload, excluding allocator/stack/process overhead.
    [[nodiscard]] std::uint64_t payload_bytes() const noexcept {
        return arithmetic_.payload_bytes();
    }
    /// Full numeric and structural visit reservation.
    [[nodiscard]] std::uint64_t reserved_steps() const noexcept {
        return work_;
    }
    /// Range2 + lower_bound bit_width(n) + equality1 comparison bound.
    [[nodiscard]] std::uint64_t maximum_comparisons() const noexcept {
        return comparisons_;
    }
    /// Checked exact arithmetic capacity used for query and node probes.
    [[nodiscard]] const RationalPlan& arithmetic() const& noexcept {
        return arithmetic_;
    }
    [[nodiscard]] const RationalPlan& arithmetic() const&& = delete;

private:
    friend struct locator_detail::Builder;
    RationalPlan arithmetic_;
    std::uint64_t nodes_ = 0, comparisons_ = 0, work_ = 0, outer_ = 0, preflight_ = 0;
};
/// Exclusive synchronous search lease. Keep it alive until the search returns.
class LocatorWorkspace {
public:
    /// Construct an idle lease without allocating.
    LocatorWorkspace() noexcept = default;
    LocatorWorkspace(const LocatorWorkspace&) = delete;
    LocatorWorkspace& operator=(const LocatorWorkspace&) = delete;

private:
    friend struct locator_detail::Builder;
    std::atomic_flag busy_ = ATOMIC_FLAG_INIT;
};
/// Owned scalar arithmetic location; no borrowed data, parent identity or curve qualification.
/// Do not use this as a cache identity for another source or as a continuity certificate.
class KnotLocation {
public:
    /// Whether the exact query equals one source knot.
    [[nodiscard]] bool is_knot() const noexcept {
        return exact_;
    }
    /// Source node count recorded by this search, without retaining the source.
    [[nodiscard]] std::uint64_t nodes() const noexcept {
        return nodes_;
    }
    /// Exact knot index; an interval-interior query returns InvalidInput.
    [[nodiscard]] core::Result<std::uint64_t> knot_index() const noexcept;
    /// Whether an in-domain left candidate exists.
    [[nodiscard]] bool has_left() const noexcept {
        return !exact_ || position_ > 0;
    }
    /// Whether an in-domain right candidate exists.
    [[nodiscard]] bool has_right() const noexcept {
        return !exact_ || position_ < nodes_ - 1;
    }
    /// Left candidate interval; nonexistent endpoint sides return InvalidInput.
    [[nodiscard]] core::Result<std::uint64_t> left_interval() const noexcept;
    /// Right candidate interval; nonexistent endpoint sides return InvalidInput.
    [[nodiscard]] core::Result<std::uint64_t> right_interval() const noexcept;
    /// Select an explicit existing side. Interior queries select the same interval for both.
    [[nodiscard]] core::Result<std::uint64_t> choose(KnotSide) const noexcept;
    /// Select the only distinct candidate; an internal knot with two candidates is InvalidInput.
    [[nodiscard]] core::Result<std::uint64_t> unique_interval() const noexcept;
    /// Actual range/search/equality comparisons, including endpoint range checks.
    [[nodiscard]] std::uint64_t comparisons() const noexcept {
        return comparisons_;
    }
    /// Actual arithmetic and outer visits, including preflight source scans.
    [[nodiscard]] std::uint64_t used_steps() const noexcept {
        return work_;
    }
    /// Actual query and node import slots.
    [[nodiscard]] std::uint64_t used_slots() const noexcept {
        return slots_;
    }

private:
    friend struct locator_detail::Builder;
    KnotLocation(std::uint64_t n, std::uint64_t p, bool exact, std::uint64_t comparisons,
                 std::uint64_t work, std::uint64_t slots, core::FailureDetails details) noexcept
        : nodes_(n),
          position_(p),
          comparisons_(comparisons),
          work_(work),
          slots_(slots),
          exact_(exact),
          details_(details) {}
    std::uint64_t nodes_, position_, comparisons_, work_, slots_;
    bool exact_;
    core::FailureDetails details_;
};
/// Scan immutable source time widths and count a bounded exact search without allocating.
[[nodiscard]] core::Result<LocatorPlan>
forecast_location(const CoefficientBlock&, RationalView, RationalLimits,
                  const std::stop_token& stop = {}, core::FailureDetails details = {}) noexcept;
/// Search only the source's closed domain. Caller preserves source/query borrows to return.
/// The selected index can be passed directly to IntervalRequest; policy decisions remain separate.
[[nodiscard]] core::Result<KnotLocation> locate_time(const CoefficientBlock&, RationalView,
                                                     LocatorWorkspace&, RationalLimits,
                                                     core::ResourcePlan&,
                                                     const std::stop_token& stop = {},
                                                     core::FailureDetails details = {}) noexcept;
} // namespace uvlas::aero::detail::pchip
