#pragma once

#include "nastran_position_enclosures.hpp"
#include "nastran_rigid_link_cards.hpp"

#include <uvlas/structure/pg_reference_geometry.hpp>

namespace uvlas::io::detail {

using NastranRigidIntervalVector = std::array<NastranCartesianInterval, 3>;
using NastranRigidExactVector = std::array<structure::ExactRational, 3>;

enum class NastranAttachmentFailure : std::uint8_t {
    Source,
    Selection,
    Thermal,
    ReferencePosition,
    DisplacementAxes,
    Pose,
    LoadSelection,
};

struct NastranAttachmentDiagnostic {
    NastranAttachmentFailure reason = NastranAttachmentFailure::Source;
    std::optional<NastranElementId> link;
    std::optional<NastranGridId> grid;
    std::optional<NastranBulkPosition> position;
    std::optional<NastranCoordinateId> coordinate;
};

struct NastranRigidAttachmentLimits {
    structure::PgReferenceLimits arithmetic;
    std::uint64_t maximum_directed_work = 0;
    std::uint64_t maximum_dependents = 0;
};

struct NastranRigidPointReference {
    NastranRigidGridReference source;
    NastranRigidIntervalVector point_metres, offset_from_master_metres;
    std::array<std::optional<NastranComponentAxis>, 3> component_axes;
};

struct NastranRigidReferenceView {
    NastranRigidLinkSource source;
    const NastranRigidLink& card;
    const NastranRigidIntervalVector& master_metres;
    std::span<const NastranRigidPointReference> dependents;
    structure::PgReferenceCallInfo exact_work;
    std::uint64_t directed_work = 0;
};

class NastranRigidAttachment {
public:
    NastranRigidAttachment(NastranRigidAttachment&&) noexcept = default;
    NastranRigidAttachment& operator=(NastranRigidAttachment&&) = delete;

    /// Own the source-ordered RBE2 points and their fixed reference CD axes.
    [[nodiscard]] static core::Result<NastranRigidAttachment>
    prepare(const NastranShellCards&, const NastranRigidLinkCards&, NastranElementId, SiLengthUnit,
            const NastranRigidAttachmentLimits&, core::ResourcePlan&,
            std::optional<NastranAttachmentDiagnostic>* = nullptr) noexcept;

    [[nodiscard]] std::optional<NastranRigidReferenceView> view() const& noexcept;
    [[nodiscard]] std::optional<NastranRigidReferenceView> view() const&& = delete;

private:
    NastranRigidAttachment(NastranRigidLinkSource, const NastranRigidLink&,
                           NastranRigidIntervalVector&&,
                           core::PlannedBuffer<NastranRigidPointReference>&&,
                           structure::PgReferenceCallInfo, std::uint64_t) noexcept;
    NastranRigidLinkSource source_;
    NastranRigidLink card_;
    NastranRigidIntervalVector master_;
    core::PlannedBuffer<NastranRigidPointReference> dependents_;
    structure::PgReferenceCallInfo exact_work_;
    std::uint64_t directed_work_;
};

struct NastranRigidMotion {
    const NastranRigidExactVector& origin_metres;
    // Row-major right-handed rotation from the reference world axes to current world axes.
    const std::array<structure::ExactRational, 9>& rotation;
    const NastranRigidExactVector& world_velocity_metres_per_second;
    const NastranRigidExactVector& world_spin_rad_per_second;
};

struct NastranRigidPointLoad {
    NastranGridId point{0};
    // Positive applied loads in selected reference CD directions. Unselected components are zero.
    NastranRigidExactVector force_n, couple_nm;
};

struct NastranRigidProjectedPoint {
    NastranGridId grid{0};
    NastranRigidIntervalVector rigid_target_metres;
    // Only selected components are exposed; a rigid target does not constrain the other lanes.
    std::array<std::optional<NastranCartesianInterval>, 3> translation_target,
        translational_velocity, spin;
    // Present when all three rotational components are selected.
    std::optional<std::array<NastranComponentAxis, 3>> rigid_target_orientation;
};

struct NastranRigidProjection {
    NastranRigidLinkSource source;
    NastranElementId link{0};
    core::PlannedBuffer<NastranRigidProjectedPoint> dependents;
    NastranRigidIntervalVector world_force_n, world_moment_about_master_nm;
    structure::PgReferenceCallInfo exact_work;
    std::uint64_t directed_work = 0;
};

// Evaluate the finite rigid point target and the CM-selected velocity directions.
// Return the transpose action of the supplied point loads on the master world twist.
[[nodiscard]] core::Result<NastranRigidProjection>
evaluate_nastran_rigid_attachment(const NastranRigidAttachment&, const NastranRigidMotion&,
                                  std::span<const NastranRigidPointLoad>,
                                  const NastranRigidAttachmentLimits&, core::ResourcePlan&,
                                  std::optional<NastranAttachmentDiagnostic>* = nullptr) noexcept;

} // namespace uvlas::io::detail
