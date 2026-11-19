#pragma once
#include "pchip_rational_slots.hpp"

#include <atomic>
#include <bit>
#include <limits>
#include <optional>
#include <utility>

namespace uvlas::aero::detail::pchip::slot_detail {
constexpr auto operation = core::Operation::EvaluateAero;
struct SlotMetadata {
    std::uint64_t numerator = 0, denominator = 1;
    int sign = 0;
};
struct WorkspaceData {
    RationalPlan plan;
    core::PlannedBuffer<std::uint32_t> limbs;
    core::PlannedBuffer<SlotMetadata> values;
    std::array<std::uint64_t, 4> scratch_used{};
    std::array<std::uint64_t, static_cast<std::size_t>(RationalOp::Count)> remaining{};
    std::atomic_flag busy = ATOMIC_FLAG_INIT;
    std::stop_token stop;
    core::FailureDetails details;
    std::optional<core::Failure> error;
    std::uint64_t epoch = 0, next = 0, steps = 0;
};
/// Nonthrowing preflight arithmetic: once overflowed, no forecast can be published.
struct CheckedCount {
    bool ok = true;
    std::uint64_t add(std::uint64_t a, std::uint64_t b) noexcept {
        if (b > std::numeric_limits<std::uint64_t>::max() - a) {
            ok = false;
            return 0;
        }
        return a + b;
    }
    std::uint64_t mul(std::uint64_t a, std::uint64_t b) noexcept {
        if (a != 0 && b > std::numeric_limits<std::uint64_t>::max() / a) {
            ok = false;
            return 0;
        }
        return a * b;
    }
};
inline core::Failure failure(core::FailureCode code, core::FailureDetails details = {}) noexcept {
    return {code, operation, details};
}
inline core::Failure supplement(core::Failure value, const core::FailureDetails& source) noexcept {
    auto& d = value.details;
    if (!d.model)
        d.model = source.model;
    if (!d.parent)
        d.parent = source.parent;
    if (!d.trial)
        d.trial = source.trial;
    if (!d.numerical_generation)
        d.numerical_generation = source.numerical_generation;
    if (!d.component)
        d.component = source.component;
    if (!d.problem)
        d.problem = source.problem;
    if (!d.resources)
        d.resources = source.resources;
    if (!d.last_accepted)
        d.last_accepted = source.last_accepted;
    return value;
}
struct IntRead {
    std::span<const std::uint32_t> words;
    [[nodiscard]] std::uint64_t bits() const noexcept {
        return words.empty() ? 0
                             : static_cast<std::uint64_t>(words.size() - 1) * 32U + 32U -
                                   static_cast<std::uint64_t>(std::countl_zero(words.back()));
    }
};
struct IntRef {
    std::span<std::uint32_t> words;
    IntRef(std::span<std::uint32_t> buffer, std::uint64_t& count) noexcept
        : words(buffer), count_(&count) {}
    [[nodiscard]] std::uint64_t& used() const noexcept {
        return *count_;
    }
    [[nodiscard]] IntRead read() const noexcept {
        return {words.first(static_cast<std::size_t>(used()))};
    }

private:
    std::uint64_t* count_;
};
struct RatRef {
    IntRef numerator, denominator;
    int& sign;
};
/// Four scratch registers: rotating GCD a/b/remainder and quotient. Add/sub reuse a/b.
/// All routine exits, including cancellation/refusal, use status rather than C++ exceptions.
struct Engine {
    WorkspaceData& data;
    bool reject(core::FailureCode) noexcept;
    [[nodiscard]] core::Failure error() const noexcept;
    bool step(std::uint64_t count = 1) noexcept;
    bool enter(RationalOp) noexcept;
    core::Result<Slot> allocate() noexcept;
    bool validate(const Slot&) noexcept;
    RatRef value(const Slot&) noexcept; // caller has already validated the live private handle.
    IntRef scratch(std::size_t) noexcept;
    RationalView view(const Slot&) noexcept;
    bool integer(IntRef, std::uint64_t) noexcept;
    bool copy(IntRef, IntRead) noexcept;
    bool compare(IntRead, IntRead, int&) noexcept;
    bool add(IntRef, IntRead, IntRead) noexcept;
    bool subtract(IntRef, IntRead, IntRead) noexcept;
    bool multiply(IntRef, IntRead, IntRead) noexcept;
    bool divide(IntRef quotient, IntRef remainder, IntRead, IntRead) noexcept;
    bool normalize(RatRef, bool require_canonical = false) noexcept;
    bool import(RatRef, RationalView) noexcept;
    bool arithmetic(RationalOp, RatRef, RationalView, RationalView) noexcept;
    bool compare(RationalView, RationalView, std::strong_ordering&) noexcept;
    bool floor(RatRef, RationalView) noexcept;

private:
    bool check_width(IntRead) noexcept;
};
core::Result<std::uint64_t> opcode_cost(RationalOp, std::uint64_t bits,
                                        std::uint64_t words) noexcept;
core::Result<std::uint64_t> workspace_bytes(std::uint64_t slots, std::uint64_t words) noexcept;
core::Result<std::uint64_t> next_epoch(std::uint64_t) noexcept;
} // namespace uvlas::aero::detail::pchip::slot_detail
