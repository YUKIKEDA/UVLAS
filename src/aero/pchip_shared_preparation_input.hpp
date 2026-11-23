#pragma once

#include "pchip_preparation_input.hpp"

#include <functional>
#include <memory>
#include <uvlas/numerics/backend_memory.hpp>

namespace uvlas::aero::detail::pchip {
namespace shared_preparation_detail {
struct Builder;
}

struct SharedPreparationLimits {
    numerics::linear::BackendMemoryLimit metadata;
    std::uint64_t steps;
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

/// Nonallocating estimate for the additional wrapper and opaque shared metadata.
[[nodiscard]] core::Result<SharedPreparationPlan>
forecast_shared_preparation(const PchipPreparationInput&, SharedPreparationLimits,
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
