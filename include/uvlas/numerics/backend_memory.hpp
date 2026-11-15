#pragma once

#include <cstdint>
#include <optional>
#include <uvlas/core/failure.hpp>
#include <variant>

namespace uvlas::numerics::linear {
namespace detail {
class BackendPlanFactory;
}

/// Factory-generated forecast for external numerical-backend payload, in bytes.
/// These values are neither measured allocations nor a process memory limit.
class BackendMemoryPlan {
public:
    /// Estimated simultaneous external payload, including declared preparation scratch.
    [[nodiscard]] std::uint64_t peak_payload_estimate_bytes() const noexcept {
        return estimate_;
    }
    /// A proven payload ceiling, when the specific backend/configuration has one.
    /// Missing proof must not be interpreted as an unlimited or zero ceiling.
    [[nodiscard]] std::optional<std::uint64_t>
    proven_peak_payload_upper_bound_bytes() const noexcept {
        return proven_;
    }

private:
    friend class detail::BackendPlanFactory;
    BackendMemoryPlan(std::uint64_t estimate, std::optional<std::uint64_t> proven) noexcept
        : estimate_(estimate), proven_(proven) {}
    std::uint64_t estimate_;
    std::optional<std::uint64_t> proven_;
};
/// Explicit permission to admit a declared estimate; does not enforce actual allocations.
struct EstimateAdmission {
    std::uint64_t bytes; ///< Maximum accepted forecast payload.
};
/// Require a genuine backend-specific payload upper bound within this limit.
struct RequireProvenBound {
    std::uint64_t bytes; ///< Maximum accepted proven payload ceiling.
};
/// Callers must explicitly choose forecast admission or a proven bound requirement.
class BackendMemoryLimit {
public:
    /// Select forecast admission explicitly; no default selection exists.
    BackendMemoryLimit(EstimateAdmission request) noexcept : selection_(request) {}
    /// Select a proven-bound requirement explicitly.
    BackendMemoryLimit(RequireProvenBound request) noexcept : selection_(request) {}
    /// Copy the selected tagged request without borrowing internal storage.
    [[nodiscard]] std::variant<EstimateAdmission, RequireProvenBound> selection() const noexcept {
        return selection_;
    }

private:
    std::variant<EstimateAdmission, RequireProvenBound> selection_;
};

/// Immutable successful comparison of a factory-generated plan and an explicit limit.
/// Passing forecast admission grants no allocation/RSS or whole-process guarantee.
class BackendAdmission {
public:
    /// Copy the factory's declared plan, preserving whether a proof exists.
    [[nodiscard]] BackendMemoryPlan plan() const noexcept {
        return plan_;
    }
    /// Explicit admission mode and limit requested by the caller.
    [[nodiscard]] BackendMemoryLimit limit() const noexcept {
        return limit_;
    }
    /// Planned amount compared to the limit, not an observed allocation count.
    [[nodiscard]] std::uint64_t admitted_payload_bytes() const noexcept {
        return selected_;
    }

private:
    friend core::Result<BackendAdmission> admit_backend_memory(BackendMemoryPlan,
                                                               BackendMemoryLimit) noexcept;
    BackendAdmission(BackendMemoryPlan plan, BackendMemoryLimit limit,
                     std::uint64_t selected) noexcept
        : plan_(plan), limit_(limit), selected_(selected) {}
    BackendMemoryPlan plan_;
    BackendMemoryLimit limit_;
    std::uint64_t selected_;
};
/// Compare a factory plan before constructing backend allocations, without allocating.
/// RequireProvenBound with no proof returns Unsupported. A selected planned amount
/// above the requested limit returns ResourceExceeded. No estimated number is reported
/// as measured required/observed resource usage in Failure. Equality passes, including
/// a zero-byte plan and limit. This does not modify ResourcePlan's owned-array ledger.
[[nodiscard]] core::Result<BackendAdmission>
admit_backend_memory(BackendMemoryPlan plan, BackendMemoryLimit limit) noexcept;
} // namespace uvlas::numerics::linear
