#include "pchip_shared_preparation_input.hpp"

#include "pchip_slot_private.hpp"

#include <optional>

namespace uvlas::aero::detail::pchip {
namespace shared_preparation_detail {
using slot_detail::failure;
using slot_detail::supplement;

struct Payload {
    std::optional<PchipPreparationInput> input;
};

/// The payload has a separate, exactly counted allocation, as in GlobalTopologyOwner.
struct Holder {
    explicit Holder(core::PlannedBuffer<Payload> storage) noexcept : payload(std::move(storage)) {}
    core::PlannedBuffer<Payload> payload;
};

constexpr std::uint64_t metadata_estimate = sizeof(Holder) + 4 * sizeof(void*);
constexpr std::uint64_t wrapper_visits = 4; // Forecast, payload, holder, publication.

struct Builder {
    struct Lease {
        SharedPreparationWorkspace& workspace;
        ~Lease() {
            workspace.busy_.clear(std::memory_order_release);
        }
    };

    static core::Result<SharedPreparationPlan> forecast(SharedPreparationLimits limits,
                                                        const std::stop_token& stop,
                                                        core::FailureDetails details) noexcept {
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));
        const auto metadata = limits.metadata.selection();
        const auto* estimate = std::get_if<numerics::linear::EstimateAdmission>(&metadata);
        if (estimate == nullptr)
            return std::unexpected(failure(core::FailureCode::Unsupported, details));
        if (estimate->bytes < metadata_estimate || limits.steps < wrapper_visits)
            return std::unexpected(failure(core::FailureCode::ResourceExceeded, details));
        return SharedPreparationPlan{sizeof(Payload), metadata_estimate, wrapper_visits};
    }

    static core::Result<SharedPreparationInput>
    share(PchipPreparationInput& source, SharedPreparationWorkspace& workspace,
          SharedPreparationLimits limits, core::ResourcePlan& memory, const std::stop_token& stop,
          core::FailureDetails details) noexcept {
        if (workspace.busy_.test_and_set(std::memory_order_acquire))
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const Lease lease{workspace};
        if (!source.active())
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        const auto plan = forecast(limits, stop, details);
        if (!plan)
            return std::unexpected(plan.error());

        auto storage = provision(memory, stop, details);
        if (!storage)
            return std::unexpected(storage.error());
        return publish(source, *storage, stop, details);
    }

    static core::Result<SharedPreparationStorage> provision(core::ResourcePlan& memory,
                                                            const std::stop_token& stop,
                                                            core::FailureDetails details) noexcept {
        auto payload = memory.allocate<Payload>(1, slot_detail::operation);
        if (!payload)
            return std::unexpected(supplement(payload.error(), details));
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));

        return core::guard_entry(slot_detail::operation, core::FailureCode::InvalidInput, details,
                                 [&]() -> core::Result<SharedPreparationStorage> {
                                     auto holder = std::make_shared<Holder>(std::move(*payload));
                                     if (stop.stop_requested())
                                         return std::unexpected(
                                             failure(core::FailureCode::Cancelled, details));
                                     return SharedPreparationStorage{std::move(holder)};
                                 });
    }

    static core::Result<SharedPreparationInput> publish(PchipPreparationInput& source,
                                                        SharedPreparationStorage& storage,
                                                        const std::stop_token& stop,
                                                        core::FailureDetails details) noexcept {
        if (!source.active() || !storage.holder_)
            return std::unexpected(failure(core::FailureCode::InvalidInput, details));
        if (stop.stop_requested())
            return std::unexpected(failure(core::FailureCode::Cancelled, details));

        // Publication is now nonthrowing; no failure can consume the old owner.
        const auto& value = storage.holder_->payload.values()[0].input.emplace(std::move(source));
        return SharedPreparationInput{
            std::shared_ptr<const PchipPreparationInput>{std::move(storage.holder_), &value}};
    }
};
} // namespace shared_preparation_detail

core::Result<std::reference_wrapper<const PchipPreparationInput>>
SharedPreparationInput::input() const& noexcept {
    if (!active())
        return std::unexpected(slot_detail::failure(core::FailureCode::InvalidInput));
    return std::cref(*input_);
}

core::Result<SharedPreparationPlan>
forecast_shared_preparation(const PchipPreparationInput& input, SharedPreparationLimits limits,
                            const std::stop_token& stop, core::FailureDetails details) noexcept {
    if (!input.active())
        return std::unexpected(slot_detail::failure(core::FailureCode::InvalidInput, details));
    return forecast_shared_preparation(limits, stop, details);
}

core::Result<SharedPreparationPlan>
forecast_shared_preparation(SharedPreparationLimits limits, const std::stop_token& stop,
                            core::FailureDetails details) noexcept {
    return shared_preparation_detail::Builder::forecast(limits, stop, details);
}

core::Result<SharedPreparationStorage>
prepare_shared_preparation_storage(const SharedPreparationPlan&, core::ResourcePlan& memory,
                                   const std::stop_token& stop,
                                   core::FailureDetails details) noexcept {
    return shared_preparation_detail::Builder::provision(memory, stop, details);
}

core::Result<SharedPreparationInput>
publish_shared_preparation(PchipPreparationInput&& source, SharedPreparationStorage&& storage,
                           const std::stop_token& stop, core::FailureDetails details) noexcept {
    return shared_preparation_detail::Builder::publish(source, storage, stop, details);
}

core::Result<SharedPreparationInput>
share_preparation_input(PchipPreparationInput&& input, SharedPreparationWorkspace& workspace,
                        SharedPreparationLimits limits, core::ResourcePlan& memory,
                        const std::stop_token& stop, core::FailureDetails details) noexcept {
    return shared_preparation_detail::Builder::share(input, workspace, limits, memory, stop,
                                                     details);
}
} // namespace uvlas::aero::detail::pchip
