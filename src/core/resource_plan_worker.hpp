#pragma once

#include <uvlas/core/resource_plan.hpp>

namespace uvlas::core::detail {
// Only the parent thread touches this reservation. Settle after all child plans
// have stopped, including failure paths with already successful allocations.
class WorkerPayloadReservation {
public:
    explicit WorkerPayloadReservation(ResourcePlan& plan) noexcept : plan_(plan) {}
    WorkerPayloadReservation(const WorkerPayloadReservation&) = delete;
    WorkerPayloadReservation& operator=(const WorkerPayloadReservation&) = delete;
    ~WorkerPayloadReservation() {
        if (active_)
            plan_.reserved_bytes_ -= reserved_;
    }
    bool begin(std::uint64_t total) noexcept {
        if (active_ || total > plan_.remaining_bytes())
            return false;
        plan_.reserved_bytes_ += total;
        reserved_ = total;
        active_ = true;
        return true;
    }
    bool settle(std::uint64_t actual) noexcept {
        if (!active_ || actual > reserved_)
            return false;
        plan_.reserved_bytes_ -= reserved_;
        plan_.allocated_bytes_ += actual;
        active_ = false;
        return true;
    }

private:
    ResourcePlan& plan_;
    std::uint64_t reserved_ = 0;
    bool active_ = false;
};
} // namespace uvlas::core::detail
