#pragma once

#include "geometry_snapshot.hpp"

#include <functional>
#include <uvlas/core/error_certificate.hpp>
#include <uvlas/core/instrumentation.hpp>

namespace uvlas::aero::detail {
class LocalPacketEngine;
class LocalBounds;
struct LocalRootScratch;
struct LocalPacketStorage;
struct LocalPacketData;
/// Explicit output representation and maximum dimensionless enclosure width.
struct LocalAccuracy {
    core::Precision representation;
    double tau;
    bool collect_timings;
    /// Timing is optional instrumentation and never changes numerical qualification.
    LocalAccuracy(core::Precision precision, double width, bool timings = false) noexcept
        : representation(precision), tau(width), collect_timings(timings) {}

    /// Reject invalid request values before preparing any facet.
    [[nodiscard]] std::optional<core::FailureCode> validation_error() const noexcept;
};
/// Wall-clock preparation phases. Root and width are subsets of enclosure, not extra work.
struct LocalPacketTimings {
    core::StageTiming classification;
    core::StageTiming qualification;
    core::StageTiming enclosure;
    core::StageTiming root;
    core::StageTiming width;
};
/// A certified real value 2^exponent * (exact word sum ﾂｱ radius).
class ScaledEnclosure {
public:
    [[nodiscard]] ScalarWords center() const noexcept {
        return center_;
    }
    [[nodiscard]] double radius() const noexcept {
        return radius_;
    }
    [[nodiscard]] int exponent() const noexcept {
        return exponent_;
    }

private:
    friend class LocalBounds;
    friend class P1Engine;
    ScaledEnclosure(ScalarWords words, double error, int scale) noexcept
        : center_(words), radius_(error), exponent_(scale) {}
    ScalarWords center_;
    double radius_;
    int exponent_;
};
/// Membership is established from stable topology and checked physical position.
struct VerifiedIncidence {
    std::variant<PointId, EdgeId, FacetId> support;
    explicit VerifiedIncidence(std::variant<PointId, EdgeId, FacetId> entity) noexcept
        : support(entity) {}
};
/// The query is outside the finite closed triangle; no proximity self exclusion.
struct Detached {};
/// A point on an oriented bound quad away from its physical contour.
struct BoundInterior {};
using LocalQualification = std::variant<VerifiedIncidence, Detached, BoundInterior>;
using ClosestEntity = std::variant<PointId, EdgeId, FacetId>;
/// The order of the twelve owned physical output enclosures.
enum class LocalQuantity : std::uint8_t {
    PositionX,
    PositionY,
    PositionZ,
    NormalX,
    NormalY,
    NormalZ,
    AreaX,
    AreaY,
    AreaZ,
    Length,
    Distance,
    Radius
};
/// Ambient order is query xyz, triangle vertices xyz, radius; quad order is the original four.
struct LocalAmbientDescriptor {
    FacetId facet;
    PhysicalFaceId face;
    BandId band;
    std::array<PointId, 3> triangle_vertices;
    std::array<EdgeId, 3> triangle_edges;
    std::array<std::int8_t, 3> positive_boundary_incidence;
    std::variant<std::array<PointId, 3>, std::array<PointId, 4>> physical_vertices;
    FacetScatter scatter;
    /// All identities and orders are required; no placeholder geometry is published.
    LocalAmbientDescriptor(FacetId target, PhysicalFaceId physical_face, BandId span_band,
                           std::array<PointId, 3> vertices, std::array<EdgeId, 3> edges,
                           std::array<std::int8_t, 3> incidence,
                           std::variant<std::array<PointId, 3>, std::array<PointId, 4>> original,
                           FacetScatter mapping) noexcept
        : facet(target),
          face(physical_face),
          band(span_band),
          triangle_vertices(vertices),
          triangle_edges(edges),
          positive_boundary_incidence(incidence),
          physical_vertices(original),
          scatter(mapping) {}
};
/// Actual planned workspace, shared exclusively across all metadata and arithmetic.
class LocalPacketWorkspace {
    class Key {
        friend class LocalPacketWorkspace;
        Key() = default;
    };

public:
    /// Factory-only key prevents unplanned or empty arithmetic storage injection.
    LocalPacketWorkspace(Key, std::unique_ptr<GeometryWorkspace>,
                         core::PlannedBuffer<LocalRootScratch>) noexcept;
    [[nodiscard]] static core::Result<std::unique_ptr<LocalPacketWorkspace>>
    create(core::ResourcePlan&, core::FailureDetails = {}) noexcept;
    LocalPacketWorkspace(const LocalPacketWorkspace&) = delete;
    LocalPacketWorkspace& operator=(const LocalPacketWorkspace&) = delete;
    LocalPacketWorkspace(LocalPacketWorkspace&&) = delete;
    LocalPacketWorkspace& operator=(LocalPacketWorkspace&&) = delete;
    ~LocalPacketWorkspace();
    /// Acquire the same lease used by every packet preparation.
    [[nodiscard]] core::Result<PredicateLease> acquire(core::FailureDetails = {}) noexcept;
    [[nodiscard]] static std::uint64_t payload_bytes() noexcept;

private:
    friend class LocalPacketEngine;
    std::unique_ptr<GeometryWorkspace> geometry_;
    core::PlannedBuffer<LocalRootScratch> scratch_;
};
/// Independent immutable local result; it never confers global geometry qualification.
class LocalGeometryPacket {
public:
    LocalGeometryPacket(LocalGeometryPacket&&) noexcept;
    LocalGeometryPacket& operator=(LocalGeometryPacket&&) noexcept;
    LocalGeometryPacket(const LocalGeometryPacket&) = delete;
    LocalGeometryPacket& operator=(const LocalGeometryPacket&) = delete;
    ~LocalGeometryPacket();
    [[nodiscard]] bool active() const noexcept;
    [[nodiscard]] AeroResult<std::reference_wrapper<const ScaledEnclosure>>
        value(LocalQuantity) const& noexcept;
    [[nodiscard]] AeroResult<std::reference_wrapper<const ScaledEnclosure>>
        value(LocalQuantity) const&& = delete;
    [[nodiscard]] AeroResult<LocalQualification> qualification() const noexcept;
    [[nodiscard]] AeroResult<ClosestEntity> closest() const noexcept;
    [[nodiscard]] AeroResult<FarRelation> far() const noexcept;
    [[nodiscard]] AeroResult<QueryRequest> query() const noexcept;
    /// Own immutable original words, topology, and context; no workspace lifetime is required.
    [[nodiscard]] AeroResult<std::reference_wrapper<const OwnedGeometry>>
    geometry() const& noexcept;
    [[nodiscard]] AeroResult<std::reference_wrapper<const OwnedGeometry>>
    geometry() const&& = delete;
    /// Positive geometric boundary is distinct from the reversed induction contour.
    [[nodiscard]] AeroResult<LocalAmbientDescriptor> ambient() const noexcept;
    [[nodiscard]] AeroResult<FacetScatter> scatter() const noexcept;
    [[nodiscard]] AeroResult<std::array<std::byte, 32>> snapshot_identity() const noexcept;
    [[nodiscard]] AeroResult<LocalAccuracy> accuracy() const noexcept;
    /// Limits admitted for this call, separate from immutable snapshot identity.
    [[nodiscard]] AeroResult<GeometryLimits> admitted_limits() const noexcept;
    [[nodiscard]] AeroResult<std::uint64_t> steps() const noexcept;
    /// Absent when disabled; timings do not qualify whole-model or physical kernel performance.
    [[nodiscard]] AeroResult<std::optional<LocalPacketTimings>> timings() const noexcept;
    [[nodiscard]] GeometryContext context() const noexcept {
        return context_;
    }
    [[nodiscard]] GlobalGeometryNotChecked global_status() const noexcept {
        return {};
    }
    [[nodiscard]] static std::uint64_t payload_bytes() noexcept;

private:
    friend class LocalPacketEngine;
    [[nodiscard]] const LocalPacketData* data() const noexcept;
    LocalGeometryPacket(core::PlannedBuffer<LocalPacketStorage>, GeometryContext) noexcept;
    core::PlannedBuffer<LocalPacketStorage> storage_;
    GeometryContext context_;
};
/// Build a current-configuration packet from original exact words and material expression.
/// Failures preserve input owners and all existing results. The target remains local-only.
[[nodiscard]] AeroResult<LocalGeometryPacket>
prepare_local_packet(const OwnedGeometry&, const ResolvedQuery&, FacetId, LocalAccuracy,
                     const GeometryLimits&, LocalPacketWorkspace&, core::ResourcePlan&,
                     const std::stop_token& = {}) noexcept;
/// Evaluate bound-contour geometry. Interior faces and artificial diagonals
/// carry no vortex source, while unqualified contact with a physical edge fails.
[[nodiscard]] AeroResult<LocalGeometryPacket>
prepare_bound_local_packet(const OwnedGeometry&, const ResolvedQuery&, FacetId, LocalAccuracy,
                           const GeometryLimits&, LocalPacketWorkspace&, core::ResourcePlan&,
                           const std::stop_token&) noexcept;
} // namespace uvlas::aero::detail
