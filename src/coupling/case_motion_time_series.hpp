#pragma once

#include "aero/pchip_curve.hpp"
#include "io/si_exact_words.hpp"

#include <array>
#include <uvlas/io/si_case.hpp>

namespace uvlas::coupling::detail {
namespace motion_time_detail {
struct Owner;
struct SnapshotOwner;
struct Builder;
} // namespace motion_time_detail

using CaseMotionField = decltype(io::schema::AeroAttachmentPrescribed::motion);

enum class CaseMotionContinuity : std::uint8_t {
    PiecewiseC2, ///< At a join the caller can request either one-sided second derivative.
    ContinuousC2 ///< The second derivatives on both sides agree at every interior join.
};

struct CaseMotionTimeWindow {
    const structure::ExactRational& begin;
    const structure::ExactRational& end;
    CaseMotionContinuity continuity = CaseMotionContinuity::PiecewiseC2;
};

struct CaseMotionTimeLimits {
    structure::ExactLimits conversion;
    aero::detail::pchip::CurveLimits curve;
    std::uint64_t maximum_exact_numbers{};
    std::uint64_t maximum_rational_operations{};
};

struct CaseMotionTimeQuery {
    CaseMotionTimeQuery(
        aero::detail::pchip::PchipInputStamp identity, const structure::ExactRational& instant,
        aero::detail::pchip::JetOrder requested = aero::detail::pchip::JetOrder::Value,
        aero::detail::pchip::RequestedSide direction =
            aero::detail::pchip::RequestedSide::Automatic) noexcept
        : stamp(identity), time(instant), order(requested), side(direction) {}

    aero::detail::pchip::PchipInputStamp stamp;
    const structure::ExactRational& time;
    aero::detail::pchip::JetOrder order;
    aero::detail::pchip::RequestedSide side;
};

/// Original constant or polynomial in ascending coefficient order, xyz per term.
struct CaseMotionFormulaSnapshotView {
    std::span<const aero::detail::pchip::RationalView> coefficients;
    /// When present: time origin, time scale, domain begin, domain end, in SI seconds.
    std::optional<std::array<aero::detail::pchip::RationalView, 4>> domain;
};

using CaseMotionSeriesSnapshotView =
    std::variant<CaseMotionFormulaSnapshotView, aero::detail::pchip::PchipSnapshotView>;

/// Borrowed model-input save. Events and the demand window retain absolute SI time.
struct CaseMotionTimeSnapshotView {
    io::SchemaIdentity source{{io::SourceId{0}, 0, {}}, 0, 0};
    std::uint64_t declaration_ordinal{};
    aero::detail::pchip::PchipInputStamp stamp{aero::detail::pchip::PchipIssuerId{0},
                                               aero::detail::pchip::PchipInputGeneration{0}};
    std::int16_t pi_power{};
    CaseMotionContinuity continuity = CaseMotionContinuity::PiecewiseC2;
    aero::detail::pchip::RationalView begin, end;
    std::span<const aero::detail::pchip::RationalView> events;
    CaseMotionSeriesSnapshotView translation, rotation;
};

class CaseMotionTimeSnapshot {
public:
    CaseMotionTimeSnapshot(CaseMotionTimeSnapshot&&) noexcept;
    CaseMotionTimeSnapshot& operator=(CaseMotionTimeSnapshot&&) noexcept;
    CaseMotionTimeSnapshot(const CaseMotionTimeSnapshot&) = delete;
    CaseMotionTimeSnapshot& operator=(const CaseMotionTimeSnapshot&) = delete;
    ~CaseMotionTimeSnapshot();

    [[nodiscard]] core::Result<CaseMotionTimeSnapshotView> view() const& noexcept;
    [[nodiscard]] core::Result<CaseMotionTimeSnapshotView> view() const&& = delete;

private:
    friend struct motion_time_detail::Builder;
    explicit CaseMotionTimeSnapshot(
        core::PlannedBuffer<motion_time_detail::SnapshotOwner>) noexcept;
    core::PlannedBuffer<motion_time_detail::SnapshotOwner> owner_;
};

/// Owns both time curves and the original event times. Frame qualification is separate.
class CaseMotionTimeSeries {
public:
    CaseMotionTimeSeries(CaseMotionTimeSeries&&) noexcept;
    CaseMotionTimeSeries& operator=(CaseMotionTimeSeries&&) noexcept;
    ~CaseMotionTimeSeries();
    CaseMotionTimeSeries(const CaseMotionTimeSeries&) = delete;
    CaseMotionTimeSeries& operator=(const CaseMotionTimeSeries&) = delete;

    [[nodiscard]] core::Result<io::SiQuantityKey> declaration() const noexcept;
    [[nodiscard]] std::span<const aero::detail::pchip::RationalView> events() const& noexcept;
    [[nodiscard]] std::span<const aero::detail::pchip::RationalView> events() const&& = delete;

private:
    friend struct motion_time_detail::Builder;
    explicit CaseMotionTimeSeries(core::PlannedBuffer<motion_time_detail::Owner>) noexcept;
    core::PlannedBuffer<motion_time_detail::Owner> owner_;
};

/// All requested derivatives in SI seconds. Rotation words carry the same pi factor at each order.
class CaseMotionTimeJet {
public:
    CaseMotionTimeJet(CaseMotionTimeJet&&) noexcept = default;
    CaseMotionTimeJet& operator=(CaseMotionTimeJet&&) noexcept = default;
    CaseMotionTimeJet(const CaseMotionTimeJet&) = delete;
    CaseMotionTimeJet& operator=(const CaseMotionTimeJet&) = delete;

    [[nodiscard]] io::SiQuantityKey declaration() const noexcept {
        return declaration_;
    }
    [[nodiscard]] aero::detail::pchip::PchipInputStamp stamp() const noexcept {
        return stamp_;
    }
    [[nodiscard]] std::int16_t rotation_pi_power() const noexcept {
        return pi_power_;
    }
    [[nodiscard]] aero::detail::pchip::JetOrder order() const noexcept {
        return order_;
    }
    [[nodiscard]] aero::detail::pchip::RequestedSide requested_side() const noexcept {
        return side_;
    }
    [[nodiscard]] aero::detail::pchip::RationalView time() const& noexcept;
    [[nodiscard]] aero::detail::pchip::RationalView time() const&& = delete;
    [[nodiscard]] core::Result<aero::detail::pchip::RationalView>
    translation(std::uint64_t component, aero::detail::pchip::JetOrder) const& noexcept;
    [[nodiscard]] core::Result<aero::detail::pchip::RationalView>
    translation(std::uint64_t, aero::detail::pchip::JetOrder) const&& = delete;
    [[nodiscard]] core::Result<aero::detail::pchip::RationalView>
    rotation_coefficient(std::uint64_t component, aero::detail::pchip::JetOrder) const& noexcept;
    [[nodiscard]] core::Result<aero::detail::pchip::RationalView>
    rotation_coefficient(std::uint64_t, aero::detail::pchip::JetOrder) const&& = delete;
    [[nodiscard]] core::Result<aero::detail::pchip::PchipJetLocation>
    translation_location() const noexcept;
    [[nodiscard]] core::Result<aero::detail::pchip::PchipJetLocation>
    rotation_location() const noexcept;

private:
    friend struct motion_time_detail::Builder;
    using FormulaValues = std::array<io::detail::SiRationalWords, 9>;
    using SeriesJet = std::variant<FormulaValues, aero::detail::pchip::ExactPchipJet>;
    CaseMotionTimeJet(SeriesJet, SeriesJet, io::detail::SiRationalWords, io::SiQuantityKey,
                      aero::detail::pchip::PchipInputStamp, std::int16_t,
                      aero::detail::pchip::JetOrder, aero::detail::pchip::RequestedSide) noexcept;
    [[nodiscard]] core::Result<aero::detail::pchip::RationalView>
    component(const SeriesJet&, std::uint64_t, aero::detail::pchip::JetOrder) const noexcept;
    [[nodiscard]] core::Result<aero::detail::pchip::PchipJetLocation>
    location(const SeriesJet&) const noexcept;
    SeriesJet translation_;
    SeriesJet rotation_;
    io::detail::SiRationalWords time_;
    io::SiQuantityKey declaration_;
    aero::detail::pchip::PchipInputStamp stamp_;
    std::int16_t pi_power_;
    aero::detail::pchip::JetOrder order_;
    aero::detail::pchip::RequestedSide side_;
};

[[nodiscard]] core::Result<CaseMotionTimeSeries> prepare_case_motion_time_series(
    const io::SiCase&, CaseMotionField, CaseMotionTimeWindow, aero::detail::pchip::PchipInputStamp,
    aero::detail::pchip::CurveWorkspace&, const CaseMotionTimeLimits&, core::ResourcePlan&,
    const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

[[nodiscard]] core::Result<CaseMotionTimeSnapshot>
snapshot_case_motion_time_series(const CaseMotionTimeSeries&, aero::detail::pchip::CurveWorkspace&,
                                 const CaseMotionTimeLimits&, core::ResourcePlan&,
                                 const std::stop_token& = {}, core::FailureDetails = {}) noexcept;
[[nodiscard]] core::Result<CaseMotionTimeSeries> restore_case_motion_time_series(
    const CaseMotionTimeSnapshotView&, const CaseMotionTimeSeries& current,
    aero::detail::pchip::CurveWorkspace&, const CaseMotionTimeLimits&, core::ResourcePlan&,
    const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

[[nodiscard]] core::Result<CaseMotionTimeJet>
evaluate_case_motion_time_series(const CaseMotionTimeSeries&, const CaseMotionTimeSeries& current,
                                 CaseMotionTimeQuery, aero::detail::pchip::CurveWorkspace&,
                                 const CaseMotionTimeLimits&, core::ResourcePlan&,
                                 const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

} // namespace uvlas::coupling::detail
