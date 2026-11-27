#pragma once

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <type_traits>
#include <utility>
#include <uvlas/core/checked_dimensions.hpp>

namespace uvlas::core {

namespace detail {
class WorkerPayloadReservation;
}
class ResourcePlan;

/// Owns a fixed-size, value-initialized array; cannot grow outside its allocation budget.
/// Moving transfers ownership and leaves an empty source. Borrowed views do not own storage.
template <class T>
class PlannedBuffer {
public:
    /// Construct an empty buffer without allocation.
    PlannedBuffer() noexcept = default;
    /// Transfer storage and invalidate views into the source's previous ownership.
    PlannedBuffer(PlannedBuffer&& other) noexcept
        : storage_(std::move(other.storage_)), size_(std::exchange(other.size_, 0)) {}
    /// Release this buffer's old storage and transfer ownership, leaving an empty source.
    PlannedBuffer& operator=(PlannedBuffer&& other) noexcept {
        if (this != &other) {
            storage_ = std::move(other.storage_);
            size_ = std::exchange(other.size_, 0);
        }
        return *this;
    }
    /// Copying is forbidden because storage has a single owner.
    PlannedBuffer(const PlannedBuffer&) = delete;
    /// Copy assignment is forbidden because storage has a single owner.
    PlannedBuffer& operator=(const PlannedBuffer&) = delete;
    /// Release owned storage; destruction does not refund the preparation budget.
    ~PlannedBuffer() = default;

    /// Number of elements, including zero for an empty or moved-from buffer.
    [[nodiscard]] std::size_t size() const noexcept {
        return size_;
    }
    /// Borrow writable elements until this buffer is moved from, assigned to or destroyed.
    [[nodiscard]] std::span<T> values() & noexcept {
        return {storage_.get(), size_};
    }
    /// Borrow read-only elements until this buffer is moved from, assigned to or destroyed.
    [[nodiscard]] std::span<const T> values() const& noexcept {
        return {storage_.get(), size_};
    }
    /// A temporary cannot supply a borrowed view.
    [[nodiscard]] std::span<T> values() && = delete;
    /// A const temporary cannot supply a borrowed view.
    [[nodiscard]] std::span<const T> values() const&& = delete;

private:
    friend class ResourcePlan;
    PlannedBuffer(std::unique_ptr<T[]> storage, std::size_t size) noexcept
        : storage_(std::move(storage)), size_(size) {}
    std::unique_ptr<T[]> storage_;
    std::size_t size_ = 0;
};

/// Enforces a cumulative payload-byte budget during buffer preparation.
/// Share one plan across arrays that must fit together. Successfully allocated payload
/// remains charged even after a buffer is destroyed, so release order cannot make an
/// incomplete plan appear to fit. Create a new plan for a new preparation, not per array.
/// Counts exclude allocator metadata, array cookies, process overhead and external-library
/// allocations; include those separately when setting an overall process memory budget.
/// Allocations made during element construction may use the same plan: in-progress
/// requests reserve their payload so nested requests cannot reuse that capacity.
/// Use from one thread at a time. This object cannot be copied or moved to duplicate credit.
class ResourcePlan {
public:
    /// Set an explicit payload limit in bytes; zero permits only empty buffers.
    explicit ResourcePlan(std::uint64_t byte_limit) noexcept : byte_limit_(byte_limit) {}
    /// Copying would duplicate the available budget and is forbidden.
    ResourcePlan(const ResourcePlan&) = delete;
    /// Copy assignment would duplicate the available budget and is forbidden.
    ResourcePlan& operator=(const ResourcePlan&) = delete;
    /// Moving a plan is forbidden; buffers own storage independently of the plan's lifetime.
    ResourcePlan(ResourcePlan&&) = delete;
    /// Move assignment is forbidden.
    ResourcePlan& operator=(ResourcePlan&&) = delete;
    /// Destroy accounting only; allocated buffers retain ownership of their storage.
    ~ResourcePlan() = default;

    /// Configured payload limit in bytes.
    [[nodiscard]] std::uint64_t byte_limit() const noexcept {
        return byte_limit_;
    }
    /// Cumulative payload bytes of successful allocations, excluding in-progress requests.
    /// A failed request is not charged, but successful nested allocations remain charged.
    [[nodiscard]] std::uint64_t allocated_bytes() const noexcept {
        return allocated_bytes_;
    }
    /// Bytes available after successful allocations and in-progress reservations.
    [[nodiscard]] std::uint64_t remaining_bytes() const noexcept {
        return byte_limit_ - allocated_bytes_ - reserved_bytes_;
    }

    /// Check dimensions and budget before allocating a value-initialized owned array.
    /// @param count Nonnegative element count; zero succeeds without allocation.
    /// @param operation Operation to identify in any returned Failure.
    /// @return Owned buffer, InvalidInput for a negative count or an element-construction
    /// exception, or ResourceExceeded for overflow, insufficient budget or bad_alloc.
    /// Budget rejection occurs before allocation and element construction. Capacity is
    /// reserved before either can reenter this plan, and that reservation is released
    /// on failure. Successful nested allocations retain their cumulative charges even
    /// if this request fails. No exception escapes; T needs a nonthrowing destructor.
    /// Resource diagnostics capture the required total and already charged or reserved
    /// payload at admission, when representable; byte_limit() is the configured cap.
    template <class T>
        requires std::default_initializable<T> && std::is_nothrow_destructible_v<T>
    [[nodiscard]] Result<PlannedBuffer<T>> allocate(Index count, Operation operation) noexcept {
        static_assert(sizeof(std::size_t) <= sizeof(std::uint64_t));
        const auto bytes = checked_byte_count(count, sizeof(T), operation);
        if (!bytes) {
            return std::unexpected(bytes.error());
        }
        FailureDetails details;
        const auto committed_or_reserved = byte_limit_ - remaining_bytes();
        if (*bytes <= std::numeric_limits<std::uint64_t>::max() - committed_or_reserved) {
            details.resources = ResourceUsage{ResourceKind::Bytes, committed_or_reserved + *bytes,
                                              committed_or_reserved};
        }
        // N4-02: compare by subtraction before allocating or adding to the ledger.
        // An unrepresentable requested total is rejected without inventing a saturated count.
        if (*bytes > remaining_bytes()) {
            return std::unexpected(Failure{FailureCode::ResourceExceeded, operation, details});
        }
        // Allocators and element constructors can reenter this plan on the same thread.
        // Reserve first, and release only this request's reservation after all exception
        // unwinding; restoring an old ledger snapshot would erase nested successes.
        reserved_bytes_ += *bytes;
        auto result = guard_entry(
            operation, FailureCode::InvalidInput, details, [&]() -> Result<PlannedBuffer<T>> {
                if (count == 0) {
                    return PlannedBuffer<T>{};
                }
                auto storage = std::make_unique<T[]>(static_cast<std::size_t>(count));
                return PlannedBuffer<T>{std::move(storage), static_cast<std::size_t>(count)};
            });
        reserved_bytes_ -= *bytes;
        if (result) {
            // The reservation guaranteed room. No user code runs between releasing
            // it and recording success; ownership transfers are nonthrowing.
            allocated_bytes_ += *bytes;
        }
        return result;
    }

private:
    friend class detail::WorkerPayloadReservation;
    std::uint64_t byte_limit_;
    std::uint64_t allocated_bytes_ = 0;
    std::uint64_t reserved_bytes_ = 0;
};

} // namespace uvlas::core
