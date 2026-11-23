#pragma once

#include "pchip_shared_preparation_input.hpp"
#include "pchip_side.hpp"

#include <memory>

namespace uvlas::aero::detail::pchip {
namespace curve_detail {
struct Builder;
struct CurveOwner;
struct JetOwner;
} // namespace curve_detail

struct CurveLimits {
    RationalLimits arithmetic;
    numerics::linear::BackendMemoryLimit metadata;
};

/// Upper bounds for one complete call. Payload is cumulative. Opaque memory is an estimate.
struct CurveAdmission {
    std::uint64_t payload_bytes = 0;
    std::uint64_t metadata_estimate_bytes = 0;
    std::uint64_t reserved_steps = 0;
    std::uint64_t maximum_integer_bits = 0;
};

class CurveWorkspace {
public:
    CurveWorkspace() noexcept = default;
    CurveWorkspace(const CurveWorkspace&) = delete;
    CurveWorkspace& operator=(const CurveWorkspace&) = delete;

private:
    friend struct curve_detail::Builder;
    std::atomic_flag busy_ = ATOMIC_FLAG_INIT;
};

/// Owns the exact input and the qualified (interval, component, power) coefficients.
class PchipCurveReady {
public:
    PchipCurveReady(PchipCurveReady&&) noexcept = default;
    PchipCurveReady& operator=(PchipCurveReady&&) noexcept = default;
    PchipCurveReady(const PchipCurveReady&) = delete;
    PchipCurveReady& operator=(const PchipCurveReady&) = delete;

    [[nodiscard]] bool active() const noexcept;
    [[nodiscard]] core::Result<PreparationRequest> input() const& noexcept;
    [[nodiscard]] core::Result<PreparationRequest> input() const&& = delete;
    [[nodiscard]] core::Result<std::uint64_t> formula_version() const noexcept;
    [[nodiscard]] core::Result<RationalView> coefficient(std::uint64_t interval,
                                                         std::uint64_t component,
                                                         std::uint64_t power) const& noexcept;
    [[nodiscard]] core::Result<RationalView> coefficient(std::uint64_t, std::uint64_t,
                                                         std::uint64_t) const&& = delete;
    [[nodiscard]] std::uint64_t used_steps() const noexcept;

private:
    friend struct curve_detail::Builder;
    explicit PchipCurveReady(std::shared_ptr<const curve_detail::CurveOwner> owner) noexcept
        : owner_(std::move(owner)) {}
    std::shared_ptr<const curve_detail::CurveOwner> owner_;
};

struct PchipQuery {
    RationalView time;
    JetOrder order = JetOrder::Value;
    RequestedSide side = RequestedSide::Automatic;
};

struct PchipJetLocation {
    std::optional<std::uint64_t> interval;
    RequestedSide requested = RequestedSide::Automatic;
    JetCandidate selected = JetCandidate::Only;
};

/// Owns exact query/result data and retains the immutable parent even if its wrapper is gone.
class ExactPchipJet {
public:
    ExactPchipJet(ExactPchipJet&&) noexcept;
    ExactPchipJet& operator=(ExactPchipJet&&) noexcept;
    ~ExactPchipJet();
    ExactPchipJet(const ExactPchipJet&) = delete;
    ExactPchipJet& operator=(const ExactPchipJet&) = delete;

    [[nodiscard]] bool active() const noexcept;
    [[nodiscard]] std::uint64_t components() const noexcept;
    [[nodiscard]] std::uint64_t used_steps() const noexcept;
    [[nodiscard]] core::Result<PchipInputStamp> stamp() const noexcept;
    [[nodiscard]] core::Result<std::uint64_t> formula_version() const noexcept;
    [[nodiscard]] core::Result<PchipJetLocation> location() const noexcept;
    [[nodiscard]] core::Result<JetOrder> order() const noexcept;
    [[nodiscard]] core::Result<RationalView> original_time() const& noexcept;
    [[nodiscard]] core::Result<RationalView> original_time() const&& = delete;
    [[nodiscard]] core::Result<RationalView> normalized_time() const& noexcept;
    [[nodiscard]] core::Result<RationalView> normalized_time() const&& = delete;
    [[nodiscard]] core::Result<RationalView> derivative(std::uint64_t component,
                                                        JetOrder) const& noexcept;
    [[nodiscard]] core::Result<RationalView> derivative(std::uint64_t, JetOrder) const&& = delete;

private:
    friend struct curve_detail::Builder;
    explicit ExactPchipJet(core::PlannedBuffer<curve_detail::JetOwner>) noexcept;
    core::PlannedBuffer<curve_detail::JetOwner> owner_;
};

/// Borrowed representation of an in-memory save. Coefficients use (interval, component, power).
struct PchipSnapshotView {
    PreparationRequest input;
    std::uint64_t formula_version;
    std::span<const RationalView> coefficients;
};

/// An owned in-memory save with immutable input and coefficient records.
class PchipSnapshot {
public:
    PchipSnapshot(PchipSnapshot&&) noexcept = default;
    PchipSnapshot& operator=(PchipSnapshot&&) noexcept = default;
    PchipSnapshot(const PchipSnapshot&) = delete;
    PchipSnapshot& operator=(const PchipSnapshot&) = delete;

    [[nodiscard]] core::Result<PchipSnapshotView> view() const& noexcept;
    [[nodiscard]] core::Result<PchipSnapshotView> view() const&& = delete;

private:
    friend struct curve_detail::Builder;
    PchipSnapshot(std::shared_ptr<const curve_detail::CurveOwner>,
                  core::PlannedBuffer<RationalView>) noexcept;
    std::shared_ptr<const curve_detail::CurveOwner> owner_;
    core::PlannedBuffer<RationalView> coefficients_;
};

[[nodiscard]] core::Result<CurveAdmission>
forecast_pchip_prepare(const PreparationRequest&, CurveLimits, const std::stop_token& stop = {},
                       core::FailureDetails details = {}) noexcept;
[[nodiscard]] core::Result<PchipCurveReady>
prepare_pchip_curve(const PreparationRequest&, CurveWorkspace&, CurveLimits, core::ResourcePlan&,
                    const std::stop_token& stop = {}, core::FailureDetails details = {}) noexcept;

/// The caller supplies its current immutable input independently of the curve.
[[nodiscard]] core::Result<CurveAdmission>
forecast_pchip_evaluation(const PchipCurveReady&, const SharedPreparationInput&, PchipQuery,
                          CurveLimits, const std::stop_token& stop = {},
                          core::FailureDetails details = {}) noexcept;
[[nodiscard]] core::Result<ExactPchipJet>
evaluate_pchip_curve(const PchipCurveReady&, const SharedPreparationInput&, PchipQuery,
                     CurveWorkspace&, CurveLimits, core::ResourcePlan&,
                     const std::stop_token& stop = {}, core::FailureDetails details = {}) noexcept;

[[nodiscard]] core::Result<CurveAdmission>
forecast_pchip_snapshot(const PchipCurveReady&, CurveLimits, const std::stop_token& stop = {},
                        core::FailureDetails details = {}) noexcept;
[[nodiscard]] core::Result<PchipSnapshot>
snapshot_pchip_curve(const PchipCurveReady&, CurveWorkspace&, CurveLimits, core::ResourcePlan&,
                     const std::stop_token& stop = {}, core::FailureDetails details = {}) noexcept;

[[nodiscard]] core::Result<CurveAdmission>
forecast_pchip_restore(const PchipSnapshotView&, const SharedPreparationInput&, CurveLimits,
                       const std::stop_token& stop = {},
                       core::FailureDetails details = {}) noexcept;
[[nodiscard]] core::Result<PchipCurveReady>
restore_pchip_curve(const PchipSnapshotView&, const SharedPreparationInput&, CurveWorkspace&,
                    CurveLimits, core::ResourcePlan&, const std::stop_token& stop = {},
                    core::FailureDetails details = {}) noexcept;

} // namespace uvlas::aero::detail::pchip
