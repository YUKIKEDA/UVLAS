#pragma once

#include <array>
#include <compare>
#include <cstdint>
#include <memory>
#include <span>
#include <stop_token>
#include <uvlas/core/resource_plan.hpp>

namespace uvlas::aero::detail::pchip {
/// Private arithmetic only: none of these types qualifies a physical curve or state.
/// Numerator and denominator words are unsigned magnitudes, least-significant word first.
struct RationalView {
    int sign = 0; ///< -1 or +1 for a nonzero numerator; zero for canonical 0/1.
    std::span<const std::uint32_t> numerator;   ///< Empty only for exact zero.
    std::span<const std::uint32_t> denominator; ///< Positive, nonempty, no leading zero word.
};
/// Explicit caps for the private preparation/evaluation arena; zero is never unlimited.
struct RationalLimits {
    std::uint64_t input_bits = 0; ///< Maximum bits in either canonical input integer.
    std::uint64_t work_bits = 0;  ///< Maximum unreduced integer width before carry storage.
    std::uint64_t steps = 0;      ///< Maximum charged limb/bit/dispatch visits.
    std::uint64_t slots = 0;      ///< Maximum number of rational SSA values.
};
/// A conservative mathematical bit shape, not a qualified rational value.
struct BitBound {
    std::uint64_t numerator = 0;   ///< Absolute numerator bit bound; zero may have bound zero.
    std::uint64_t denominator = 1; ///< Positive denominator bit bound.
};
/// Closed set of arithmetic bodies counted before any storage allocation.
enum class RationalOp : std::uint8_t {
    Input,
    Clone,
    Negate,
    Add,
    Subtract,
    Multiply,
    Divide,
    Compare,
    Floor,
    Count
};
namespace slot_detail {
struct WorkspaceData;
struct Engine;
} // namespace slot_detail
/// Immutable capacity forecast; constructed by Counter without allocating a DAG.
class RationalPlan {
public:
    /// Construct an invalid/unprepared plan, rejected by Workspace::create.
    RationalPlan() noexcept = default;
    /// Whether nonallocating counting completed successfully.
    [[nodiscard]] bool valid() const noexcept {
        return valid_;
    }
    /// Exact total C++ object/array payload reserved for a new workspace.
    [[nodiscard]] std::uint64_t payload_bytes() const noexcept {
        return payload_;
    }
    /// Maximum integer mathematical width, before one carry bit and word rounding.
    [[nodiscard]] std::uint64_t integer_bits() const noexcept {
        return bits_;
    }
    /// Prepared 32-bit words per integer slot, including carry capacity.
    [[nodiscard]] std::uint64_t words() const noexcept {
        return words_;
    }
    /// Maximum rational SSA slots, including imported constants/values.
    [[nodiscard]] std::uint64_t slots() const noexcept {
        return slots_;
    }
    /// Complete reserved work, including dispatch and normalization visits.
    [[nodiscard]] std::uint64_t reserved_steps() const noexcept {
        return work_;
    }

private:
    friend class Counter;
    friend class Workspace;
    friend class Session;
    friend struct slot_detail::Engine;
    RationalLimits limits_{};
    std::array<std::uint64_t, static_cast<std::size_t>(RationalOp::Count)> counts_{};
    std::uint64_t bits_ = 0, words_ = 0, slots_ = 0, work_ = 0, payload_ = 0;
    bool valid_ = false;
};
/// Nonallocating capacity visitor. Bounds use unreduced cross products, never hoped-for GCDs.
class Counter {
public:
    /// Keep explicit limits/diagnostic values; no allocation or numeric backend is invoked.
    explicit Counter(RationalLimits, core::FailureDetails details = {}) noexcept;
    /// Count a canonical input shape; denominator_bits must be positive.
    [[nodiscard]] core::Result<BitBound> input(BitBound) noexcept;
    /// Count an exact signed 64-bit constant, including INT64_MIN without signed overflow.
    [[nodiscard]] core::Result<BitBound> integer(std::int64_t) noexcept;
    /// Count a copied or negated value with unchanged mathematical bit bounds.
    [[nodiscard]] core::Result<BitBound> unary(RationalOp, BitBound) noexcept;
    /// Count add/subtract/multiply/divide, rejecting other opcode values.
    [[nodiscard]] core::Result<BitBound> binary(RationalOp, BitBound, BitBound) noexcept;
    /// Count a signed exact comparison, including both cross products; allocates no value slot.
    [[nodiscard]] core::Result<void> compare(BitBound, BitBound) noexcept;
    /// Count mathematical floor, including the negative noninteger correction.
    [[nodiscard]] core::Result<BitBound> floor(BitBound) noexcept;
    /// Complete checked bit/slot/byte/work admission; no memory is acquired.
    [[nodiscard]] core::Result<RationalPlan> finish() const noexcept;

private:
    friend struct slot_detail::Engine;
    core::Result<BitBound> count(RationalOp, BitBound, std::uint64_t intermediate) noexcept;
    RationalPlan plan_;
    core::FailureDetails details_;
};
/// Ephemeral private reference. Valid only in its producing Session, while its Workspace lives.
/// It never owns a product result, and must not be cached outside that session.
/// Epoch checks cover re-leasing the same still-live Workspace only. A Slot from a destroyed
/// Workspace is outside its borrow lifetime and must never be submitted to a new Workspace.
class Slot {
public:
    /// No public construction can manufacture a reference to an uncomputed value.
    Slot(const Slot&) noexcept = default;
    /// Copy the private reference without copying arithmetic storage.
    Slot& operator=(const Slot&) noexcept = default;

private:
    friend class Session;
    friend struct slot_detail::Engine;
    Slot(const slot_detail::WorkspaceData* owner, std::uint64_t epoch, std::uint64_t index) noexcept
        : owner_(owner), epoch_(epoch), index_(index) {}
    const slot_detail::WorkspaceData* owner_;
    std::uint64_t epoch_, index_;
};
/// Actual usage of one leased numerical session, separate from its upper-bound reservation.
struct SlotUsage {
    std::uint64_t steps = 0; ///< Charged limb/bit/dispatch visits already executed.
    std::uint64_t slots = 0; ///< SSA slots already assigned, including a failed unpublished slot.
};
class Session;
/// Move-only owner of prepared integer slots. All payload is charged through ResourcePlan.
/// Allocator overhead/RSS are separate; no claim of an allocator metadata bound is made.
/// Owned storage must outlive every borrowing Session/view. Destroying or move-assigning over
/// the last owner while a Session still borrows that storage violates this private lifetime
/// contract.
class Workspace {
public:
    /// Transfer owned storage; existing private references follow storage, not wrapper addresses.
    Workspace(Workspace&&) noexcept;
    /// Release previous storage and transfer ownership; no active session may borrow the old owner.
    Workspace& operator=(Workspace&&) noexcept;
    /// After all borrowing Sessions have ended, release arrays without refunding payload.
    ~Workspace();
    Workspace(const Workspace&) = delete;
    Workspace& operator=(const Workspace&) = delete;
    /// Preflight the complete forecast, then allocate one state, limb array, and slot metadata
    /// array. Failure preserves supplied diagnostics and earlier successful cumulative charges.
    [[nodiscard]] static core::Result<Workspace> create(const RationalPlan&, core::ResourcePlan&,
                                                        core::FailureDetails details = {},
                                                        const std::stop_token& stop = {}) noexcept;
    /// Acquire exclusive use, advance the checked session epoch, and reset arithmetic usage.
    /// Reentry returns InvalidInput; cancellation returns Cancelled without arithmetic.
    [[nodiscard]] core::Result<Session> begin(std::stop_token stop = {},
                                              core::FailureDetails details = {}) & noexcept;
    /// Moving-from or default-inactive owners cannot supply private scratch.
    [[nodiscard]] bool active() const noexcept;

private:
    friend class Session;
    explicit Workspace(core::PlannedBuffer<slot_detail::WorkspaceData>) noexcept;
    core::PlannedBuffer<slot_detail::WorkspaceData> data_;
};
/// Exclusive lease over private scratch, never a solver state or a qualification packet.
/// A numeric failure latches the session; ending it releases the lease for another call.
class Session {
public:
    /// Transfer the lease without allocating.
    Session(Session&&) noexcept;
    /// Release an earlier lease and transfer the source lease.
    Session& operator=(Session&&) noexcept;
    /// Release only this session's lease, including after numeric failure.
    ~Session();
    Session(const Session&) = delete;
    Session& operator=(const Session&) = delete;
    /// Copy and check canonical input after admission. Nonreduced inputs are InvalidInput.
    [[nodiscard]] core::Result<Slot> input(RationalView) noexcept;
    /// Materialize one signed integer from prepared slots; INT64_MIN is accepted exactly.
    [[nodiscard]] core::Result<Slot> integer(std::int64_t) noexcept;
    /// Copy or negate one live slot; old slots are never modified.
    [[nodiscard]] core::Result<Slot> unary(RationalOp, const Slot&) noexcept;
    /// Add/subtract/multiply/divide live values; division by zero is InvalidInput.
    [[nodiscard]] core::Result<Slot> binary(RationalOp, const Slot&, const Slot&) noexcept;
    /// Compare canonical values exactly without materializing a rounded difference.
    [[nodiscard]] core::Result<std::strong_ordering> compare(const Slot&, const Slot&) noexcept;
    /// Return mathematical floor, including negative nonintegers.
    [[nodiscard]] core::Result<Slot> floor(const Slot&) noexcept;
    /// Borrow canonical output until this Session ends/resets or its Workspace is destroyed.
    /// Copy into the future product owner's planned storage before releasing this lease.
    [[nodiscard]] core::Result<RationalView> view(const Slot&) const& noexcept;
    [[nodiscard]] core::Result<RationalView> view(const Slot&) const&& = delete;
    /// Read actual usage even after a latched arithmetic failure; closed leases are stale.
    [[nodiscard]] core::Result<SlotUsage> usage() const noexcept;

private:
    friend class Workspace;
    explicit Session(slot_detail::WorkspaceData*) noexcept;
    slot_detail::WorkspaceData* data_ = nullptr;
};
} // namespace uvlas::aero::detail::pchip
