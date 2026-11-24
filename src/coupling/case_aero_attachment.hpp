#pragma once

#include "case_motion.hpp"
#include "discrete_nonpenetration.hpp"

namespace uvlas::coupling::detail {

/// A fixed collocation expression on the reference facet. The reference chart
/// and circulation words are supplied by the discrete model.
struct PrescribedAeroRow {
    aero::detail::QueryId query{0};
    aero::detail::SourceFacePoint material_point;
    aero::detail::ScalarWords circulation, circulation_direction;
};

/// One prescribed attachment in a qualified discrete model. Only the listed
/// canonical points and their registered aliases move with this attachment.
/// The moving set is exactly the vertices of the selected bound quadrilaterals.
/// Other current points, such as released wake particles, keep their supplied words.
struct PrescribedAeroReference {
    const aero::detail::GlobalTopology& topology;
    const aero::detail::GeometryInput& geometry;
    std::span<const aero::detail::PointId> moving_points;
    std::span<const PrescribedAeroRow> rows;
};

struct PrescribedAeroLimits {
    CaseMotionLimits motion;
    aero::detail::GeometryLimits geometry;
    aero::detail::GlobalLimits contact;
    std::uint64_t work{};
};

struct PrescribedAeroProjection {
    PrescribedAeroProjection(core::Precision format,
                             std::span<const CaseFlowComponentBudget> pose_components,
                             std::span<const CaseFlowComponentBudget> point_components) noexcept
        : representation(format), pose(pose_components), point(point_components) {}
    core::Precision representation;
    std::span<const CaseFlowComponentBudget> pose;  ///< Value and first derivative.
    std::span<const CaseFlowComponentBudget> point; ///< Shared xyz position/velocity budgets.
};

namespace case_clock_detail {
struct Builder;
}

namespace aero_attachment_detail {
struct Builder;
struct Trial;
} // namespace aero_attachment_detail

class PrescribedAeroAttachment {
public:
    PrescribedAeroAttachment(PrescribedAeroAttachment&&) noexcept = default;
    PrescribedAeroAttachment& operator=(PrescribedAeroAttachment&&) noexcept = default;
    PrescribedAeroAttachment(const PrescribedAeroAttachment&) = delete;
    PrescribedAeroAttachment& operator=(const PrescribedAeroAttachment&) = delete;

private:
    friend struct aero_attachment_detail::Builder;
    friend struct case_clock_detail::Builder;
    explicit PrescribedAeroAttachment(CaseMotionReady&& motion) noexcept
        : motion_(std::move(motion)) {}
    CaseMotionReady motion_;
};

class PrescribedAeroTrial {
public:
    PrescribedAeroTrial(PrescribedAeroTrial&&) noexcept;
    PrescribedAeroTrial& operator=(PrescribedAeroTrial&&) noexcept;
    PrescribedAeroTrial(const PrescribedAeroTrial&) = delete;
    PrescribedAeroTrial& operator=(const PrescribedAeroTrial&) = delete;
    ~PrescribedAeroTrial();

    [[nodiscard]] core::Result<std::reference_wrapper<const CaseMotionEvaluation>>
    motion() const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const aero::detail::GlobalContactState>>
    contact() const& noexcept;
    [[nodiscard]] core::Result<io::SiEnclosure>
    material_point(aero::detail::PointId, std::size_t axis,
                   aero::detail::pchip::JetOrder) const noexcept;
    [[nodiscard]] std::size_t collocation_count() const noexcept;
    [[nodiscard]] core::Result<DiscreteCollocation> collocation(std::size_t) const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const CaseMotionEvaluation>>
    motion() const&& = delete;
    [[nodiscard]] core::Result<std::reference_wrapper<const aero::detail::GlobalContactState>>
    contact() const&& = delete;
    [[nodiscard]] core::Result<DiscreteCollocation> collocation(std::size_t) const&& = delete;

private:
    friend struct aero_attachment_detail::Builder;
    explicit PrescribedAeroTrial(core::PlannedBuffer<aero_attachment_detail::Trial>&&) noexcept;
    core::PlannedBuffer<aero_attachment_detail::Trial> trial_;
};

[[nodiscard]] core::Result<PrescribedAeroAttachment>
prepare_prescribed_aero_attachment(const io::SiCase&, const io::schema::AeroAttachmentPrescribed&,
                                   CaseMotionTimeWindow, aero::detail::pchip::PchipInputStamp,
                                   aero::detail::pchip::CurveWorkspace&, const CaseMotionLimits&,
                                   core::ResourcePlan&, const std::stop_token& = {},
                                   core::FailureDetails = {}) noexcept;

/// Produces material collocations for the checked rounded geometry and the
/// qualified pose's velocity, with the same original facet IDs and exact time.
[[nodiscard]] core::Result<PrescribedAeroTrial> evaluate_prescribed_aero_attachment(
    const PrescribedAeroAttachment&, const PrescribedAeroAttachment& current,
    const PrescribedAeroReference&, CaseMotionTimeQuery, const PrescribedAeroProjection&,
    aero::detail::pchip::CurveWorkspace&, aero::detail::GeometryWorkspace&,
    aero::detail::GlobalWorkspace&, const PrescribedAeroLimits&, core::ResourcePlan&,
    const std::stop_token& = {}, core::FailureDetails = {}) noexcept;

} // namespace uvlas::coupling::detail
