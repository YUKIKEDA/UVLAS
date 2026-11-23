#pragma once

#include "pchip_curve.hpp"

#include <array>
#include <variant>

namespace uvlas::aero::detail {

/// Exact SI values in world X, Y, Z order. Time is in seconds.
struct ConstantUniformFlow {
    std::array<pchip::RationalView, 3> velocity;
};

/// Coefficients in ascending powers of (t - origin)/scale, with XYZ within each power.
struct PolynomialUniformFlow {
    pchip::RationalView origin, scale, begin, end;
    std::span<const pchip::RationalView> coefficients;
};

struct PchipUniformFlow {
    pchip::SampleGrid samples;
    pchip::BoundaryRequest boundary;
};

/// AE-07, FD-06. No-gust uniform velocity input. Borrows are synchronous.
struct UniformFlowInput {
    using Series = std::variant<ConstantUniformFlow, PolynomialUniformFlow, PchipUniformFlow>;

    UniformFlowInput(pchip::PchipInputStamp input_stamp,
                     std::span<const pchip::RationalView> input_events,
                     Series input_series) noexcept
        : stamp(input_stamp), events(input_events), series(input_series) {}

    pchip::PchipInputStamp stamp;
    std::span<const pchip::RationalView> events;
    Series series;
};

namespace uniform_flow_detail {
struct Formula;
struct Builder;
} // namespace uniform_flow_detail

class UniformFlowWorkspace {
public:
    UniformFlowWorkspace() noexcept = default;
    UniformFlowWorkspace(const UniformFlowWorkspace&) = delete;
    UniformFlowWorkspace& operator=(const UniformFlowWorkspace&) = delete;

private:
    friend struct uniform_flow_detail::Builder;
    std::atomic_flag entered_ = ATOMIC_FLAG_INIT;
    pchip::CurveWorkspace pchip_;
};

/// Immutable source used by the value and both physical derivatives.
class UniformFlowReady {
public:
    UniformFlowReady(UniformFlowReady&&) noexcept = default;
    UniformFlowReady& operator=(UniformFlowReady&&) noexcept = default;
    UniformFlowReady(const UniformFlowReady&) = delete;
    UniformFlowReady& operator=(const UniformFlowReady&) = delete;
    [[nodiscard]] bool active() const noexcept;
    [[nodiscard]] core::Result<UniformFlowInput> input() const& noexcept;
    [[nodiscard]] core::Result<UniformFlowInput> input() const&& = delete;
    [[nodiscard]] std::uint64_t used_steps() const noexcept;

private:
    friend struct uniform_flow_detail::Builder;
    explicit UniformFlowReady(std::shared_ptr<const uniform_flow_detail::Formula>,
                              std::uint64_t steps) noexcept;
    explicit UniformFlowReady(pchip::PchipCurveReady) noexcept;
    std::variant<std::shared_ptr<const uniform_flow_detail::Formula>, pchip::PchipCurveReady> data_;
    std::uint64_t formula_work_ = 0;
};

struct UniformFlowQuery {
    pchip::PchipInputStamp stamp;
    pchip::RationalView time;
    pchip::RequestedSide side = pchip::RequestedSide::Automatic;
};

/// Exact world velocity [m/s], time partial [m/s²] and spatial Jacobian [1/s].
/// The spatial partial is identically zero at every world position.
class ExactUniformFlowJet {
public:
    ExactUniformFlowJet(ExactUniformFlowJet&&) noexcept = default;
    ExactUniformFlowJet& operator=(ExactUniformFlowJet&&) noexcept = default;
    ExactUniformFlowJet(const ExactUniformFlowJet&) = delete;
    ExactUniformFlowJet& operator=(const ExactUniformFlowJet&) = delete;
    [[nodiscard]] bool active() const noexcept;
    [[nodiscard]] core::Result<pchip::PchipInputStamp> stamp() const noexcept;
    [[nodiscard]] core::Result<pchip::RationalView> time() const& noexcept;
    [[nodiscard]] core::Result<pchip::RationalView> time() const&& = delete;
    [[nodiscard]] core::Result<pchip::RationalView> velocity(std::size_t component) const& noexcept;
    [[nodiscard]] core::Result<pchip::RationalView> velocity(std::size_t) const&& = delete;
    [[nodiscard]] core::Result<pchip::RationalView>
    time_partial(std::size_t component) const& noexcept;
    [[nodiscard]] core::Result<pchip::RationalView> time_partial(std::size_t) const&& = delete;
    [[nodiscard]] core::Result<pchip::RationalView>
    spatial_partial(std::size_t component, std::size_t coordinate) const noexcept;
    [[nodiscard]] std::uint64_t used_steps() const noexcept;

private:
    friend struct uniform_flow_detail::Builder;
    struct FormulaJet {
        std::shared_ptr<const uniform_flow_detail::Formula> parent;
        core::PlannedBuffer<std::uint32_t> words;
        std::array<pchip::RationalView, 7> records;
        std::uint64_t steps;
    };
    [[nodiscard]] const FormulaJet* active_formula() const noexcept;
    explicit ExactUniformFlowJet(FormulaJet) noexcept;
    explicit ExactUniformFlowJet(pchip::ExactPchipJet) noexcept;
    std::variant<FormulaJet, pchip::ExactPchipJet> data_;
};

[[nodiscard]] core::Result<pchip::CurveAdmission>
forecast_uniform_flow_prepare(const UniformFlowInput&, pchip::CurveLimits,
                              const std::stop_token& = {}, core::FailureDetails = {}) noexcept;
[[nodiscard]] core::Result<UniformFlowReady>
prepare_uniform_flow(const UniformFlowInput&, UniformFlowWorkspace&, pchip::CurveLimits,
                     core::ResourcePlan&, const std::stop_token& = {},
                     core::FailureDetails = {}) noexcept;

[[nodiscard]] core::Result<pchip::CurveAdmission>
forecast_uniform_flow_evaluation(const UniformFlowReady&, UniformFlowQuery, pchip::CurveLimits,
                                 const std::stop_token& = {}, core::FailureDetails = {}) noexcept;
[[nodiscard]] core::Result<ExactUniformFlowJet>
evaluate_uniform_flow(const UniformFlowReady&, UniformFlowQuery, UniformFlowWorkspace&,
                      pchip::CurveLimits, core::ResourcePlan&, const std::stop_token& = {},
                      core::FailureDetails = {}) noexcept;

struct UniformFormulaSnapshotView {
    UniformFlowInput input;
    std::uint64_t formula_version = 0;
};

using UniformFlowSnapshotView = std::variant<UniformFormulaSnapshotView, pchip::PchipSnapshotView>;

/// Owns the source, and the actual coefficient block for the PCHIP case.
class UniformFlowSnapshot {
public:
    UniformFlowSnapshot(UniformFlowSnapshot&&) noexcept = default;
    UniformFlowSnapshot& operator=(UniformFlowSnapshot&&) noexcept = default;
    UniformFlowSnapshot(const UniformFlowSnapshot&) = delete;
    UniformFlowSnapshot& operator=(const UniformFlowSnapshot&) = delete;
    [[nodiscard]] core::Result<UniformFlowSnapshotView> view() const& noexcept;
    [[nodiscard]] core::Result<UniformFlowSnapshotView> view() const&& = delete;

private:
    friend struct uniform_flow_detail::Builder;
    explicit UniformFlowSnapshot(UniformFlowReady) noexcept;
    explicit UniformFlowSnapshot(pchip::PchipSnapshot) noexcept;
    std::variant<UniformFlowReady, pchip::PchipSnapshot> data_;
};

[[nodiscard]] core::Result<pchip::CurveAdmission>
forecast_uniform_flow_snapshot(const UniformFlowReady&, pchip::CurveLimits,
                               const std::stop_token& = {}, core::FailureDetails = {}) noexcept;
[[nodiscard]] core::Result<UniformFlowSnapshot>
snapshot_uniform_flow(const UniformFlowReady&, UniformFlowWorkspace&, pchip::CurveLimits,
                      core::ResourcePlan&, const std::stop_token& = {},
                      core::FailureDetails = {}) noexcept;

/// Restore against the independent, presently qualified uniform-flow source.
[[nodiscard]] core::Result<pchip::CurveAdmission>
forecast_uniform_flow_restore(const UniformFlowSnapshotView&, const UniformFlowReady& current,
                              pchip::CurveLimits, const std::stop_token& = {},
                              core::FailureDetails = {}) noexcept;
[[nodiscard]] core::Result<UniformFlowReady>
restore_uniform_flow(const UniformFlowSnapshotView&, const UniformFlowReady& current,
                     UniformFlowWorkspace&, pchip::CurveLimits, core::ResourcePlan&,
                     const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

} // namespace uvlas::aero::detail
