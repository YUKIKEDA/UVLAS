#pragma once

#include "pchip_preparation_input.hpp"

#include <functional>
#include <memory>
#include <uvlas/numerics/backend_memory.hpp>

namespace uvlas::aero::detail::pchip {
namespace shared_preparation_detail {
struct Builder;
struct Holder;
} // namespace shared_preparation_detail

struct SharedPreparationLimits {
    numerics::linear::BackendMemoryLimit metadata;
    std::uint64_t steps = 0;
};

/// Additional charge for transferring an existing input into shared ownership.
class SharedPreparationPlan {
public:
    [[nodiscard]] std::uint64_t payload_bytes() const noexcept {
        return payload_;
    }
    [[nodiscard]] std::uint64_t metadata_estimate_bytes() const noexcept {
        return metadata_;
    }
    [[nodiscard]] std::uint64_t reserved_steps() const noexcept {
        return steps_;
    }

private:
    friend struct shared_preparation_detail::Builder;
    SharedPreparationPlan(std::uint64_t payload, std::uint64_t metadata,
                          std::uint64_t steps) noexcept
        : payload_(payload), metadata_(metadata), steps_(steps) {}
    std::uint64_t payload_, metadata_, steps_;
};

class SharedPreparationWorkspace {
public:
    SharedPreparationWorkspace() noexcept = default;
    SharedPreparationWorkspace(const SharedPreparationWorkspace&) = delete;
    SharedPreparationWorkspace& operator=(const SharedPreparationWorkspace&) = delete;

private:
    friend struct shared_preparation_detail::Builder;
    std::atomic_flag busy_ = ATOMIC_FLAG_INIT;
};

/// Unpublished wrapper and metadata, allocated before exact preparation.
class SharedPreparationStorage {
public:
    SharedPreparationStorage(SharedPreparationStorage&&) noexcept = default;
    SharedPreparationStorage& operator=(SharedPreparationStorage&&) noexcept = default;

private:
    friend struct shared_preparation_detail::Builder;
    explicit SharedPreparationStorage(
        std::shared_ptr<shared_preparation_detail::Holder> holder) noexcept
        : holder_(std::move(holder)) {}
    std::shared_ptr<shared_preparation_detail::Holder> holder_;
};

/// Read-only handle to the original owned input. Copies retain the same immutable input.
class SharedPreparationInput {
public:
    SharedPreparationInput(const SharedPreparationInput&) noexcept = default;
    SharedPreparationInput& operator=(const SharedPreparationInput&) noexcept = default;
    SharedPreparationInput(SharedPreparationInput&&) noexcept = default;
    SharedPreparationInput& operator=(SharedPreparationInput&&) noexcept = default;

    [[nodiscard]] bool active() const noexcept {
        return input_ != nullptr;
    }

    /// Compare the live shared ownership of the immutable input.
    [[nodiscard]] bool shares_input_with(const SharedPreparationInput& other) const noexcept {
        return active() && other.active() && !input_.owner_before(other.input_) &&
               !other.input_.owner_before(input_);
    }

    /// Borrow the exact input during this handle's lifetime.
    [[nodiscard]] core::Result<std::reference_wrapper<const PchipPreparationInput>>
    input() const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const PchipPreparationInput>>
    input() const&& = delete;

private:
    friend struct shared_preparation_detail::Builder;
    explicit SharedPreparationInput(std::shared_ptr<const PchipPreparationInput> input) noexcept
        : input_(std::move(input)) {}
    std::shared_ptr<const PchipPreparationInput> input_;
};

/// Nonallocating estimate before the original input has been copied.
[[nodiscard]] core::Result<SharedPreparationPlan>
forecast_shared_preparation(SharedPreparationLimits, const std::stop_token& stop = {},
                            core::FailureDetails details = {}) noexcept;

/// Nonallocating estimate for the additional wrapper and opaque shared metadata.
[[nodiscard]] core::Result<SharedPreparationPlan>
forecast_shared_preparation(const PchipPreparationInput&, SharedPreparationLimits,
                            const std::stop_token& stop = {},
                            core::FailureDetails details = {}) noexcept;

[[nodiscard]] core::Result<SharedPreparationStorage>
prepare_shared_preparation_storage(const SharedPreparationPlan&, core::ResourcePlan&,
                                   const std::stop_token& stop = {},
                                   core::FailureDetails details = {}) noexcept;

/// Publish into a prepared wrapper without allocating.
[[nodiscard]] core::Result<SharedPreparationInput>
publish_shared_preparation(PchipPreparationInput&&, SharedPreparationStorage&&,
                           const std::stop_token& stop = {},
                           core::FailureDetails details = {}) noexcept;

/// Consume the source only after all allocation, admission and cancellation checks succeed.
/// The caller keeps exclusive access to the source wrapper for the synchronous call.
[[nodiscard]] core::Result<SharedPreparationInput>
share_preparation_input(PchipPreparationInput&&, SharedPreparationWorkspace&,
                        SharedPreparationLimits, core::ResourcePlan&,
                        const std::stop_token& stop = {},
                        core::FailureDetails details = {}) noexcept;

} // namespace uvlas::aero::detail::pchip
