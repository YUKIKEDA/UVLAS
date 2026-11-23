#pragma once

#include "pchip_boundary.hpp"

#include <uvlas/numerics/storage.hpp>

namespace uvlas::aero::detail::pchip {

namespace preparation_detail {
struct IssuerTag;
struct InputGenerationTag;
struct Builder;
struct Owner;
} // namespace preparation_detail

using PchipIssuerId = core::Identifier<preparation_detail::IssuerTag>;
using PchipInputGeneration = core::Identifier<preparation_detail::InputGenerationTag>;

/// Supplied and advanced by the preparation host, independently of content.
struct PchipInputStamp {
    PchipIssuerId issuer;
    PchipInputGeneration generation;
    bool operator==(const PchipInputStamp&) const noexcept = default;
};

/// Borrows remain immutable for the synchronous count and owned copy.
struct PreparationRequest {
    PreparationRequest(SampleGrid source_samples,
                       std::span<const numerics::SiDimension> source_dimensions,
                       const BoundaryRequest& source_boundary,
                       std::span<const RationalView> source_events,
                       PchipInputStamp input_stamp) noexcept
        : samples(source_samples),
          dimensions(source_dimensions),
          boundary(source_boundary),
          events(source_events),
          stamp(input_stamp) {}

    SampleGrid samples;
    std::span<const numerics::SiDimension> dimensions;
    BoundaryRequest boundary;
    std::span<const RationalView> events;
    PchipInputStamp stamp;
};

class PreparationPlan {
public:
    [[nodiscard]] std::uint64_t payload_bytes() const noexcept {
        return bytes_;
    }
    [[nodiscard]] std::uint64_t reserved_steps() const noexcept {
        return work_;
    }
    [[nodiscard]] const RationalPlan& arithmetic() const& noexcept {
        return arithmetic_;
    }
    [[nodiscard]] const RationalPlan& arithmetic() const&& = delete;

private:
    friend struct preparation_detail::Builder;
    RationalPlan arithmetic_;
    std::uint64_t records_ = 0, limbs_ = 0, bytes_ = 0, work_ = 0, scan_ = 0;
};

class PreparationWorkspace {
public:
    PreparationWorkspace() noexcept = default;
    PreparationWorkspace(const PreparationWorkspace&) = delete;
    PreparationWorkspace& operator=(const PreparationWorkspace&) = delete;

private:
    friend struct preparation_detail::Builder;
    std::atomic_flag busy_ = ATOMIC_FLAG_INIT;
};

/// Owned canonical preparation input, independent of coefficient and curve qualification.
/// The view borrows this owner; moving it leaves the source wrapper inactive.
class PchipPreparationInput {
public:
    PchipPreparationInput(PchipPreparationInput&&) noexcept;
    PchipPreparationInput& operator=(PchipPreparationInput&&) noexcept;
    ~PchipPreparationInput();
    PchipPreparationInput(const PchipPreparationInput&) = delete;
    PchipPreparationInput& operator=(const PchipPreparationInput&) = delete;

    [[nodiscard]] bool active() const noexcept;
    [[nodiscard]] core::Result<PreparationRequest> view() const& noexcept;
    [[nodiscard]] core::Result<PreparationRequest> view() const&& = delete;
    [[nodiscard]] core::Result<std::uint64_t> formula_version() const noexcept;
    [[nodiscard]] std::uint64_t used_steps() const noexcept;

private:
    friend struct preparation_detail::Builder;
    explicit PchipPreparationInput(core::PlannedBuffer<preparation_detail::Owner>) noexcept;
    core::PlannedBuffer<preparation_detail::Owner> owner_;
};

[[nodiscard]] core::Result<PreparationPlan>
forecast_preparation_input(const PreparationRequest&, RationalLimits,
                           const std::stop_token& stop = {},
                           core::FailureDetails details = {}) noexcept;

[[nodiscard]] core::Result<PchipPreparationInput>
make_preparation_input(const PreparationRequest&, PreparationWorkspace&, RationalLimits,
                       core::ResourcePlan&, const std::stop_token& stop = {},
                       core::FailureDetails details = {}) noexcept;

/// No allocation or hash-based shortcut. The expected input comes from the current host.
/// Returns charged work on equality, StaleGeneration on any stamp/content mismatch.
[[nodiscard]] core::Result<std::uint64_t>
match_preparation_input(const PchipPreparationInput& actual, const PchipPreparationInput& expected,
                        std::uint64_t work_limit, const std::stop_token& stop = {},
                        core::FailureDetails details = {}) noexcept;

} // namespace uvlas::aero::detail::pchip
