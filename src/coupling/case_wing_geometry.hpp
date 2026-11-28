#pragma once

#include "aero/global_geometry.hpp"
#include "case_flow_clock.hpp"
#include "case_motion.hpp"
#include "case_wing_chart.hpp"

namespace uvlas::coupling::detail {

struct CaseWingProjection {
    core::Precision representation;
    CaseFlowProjectionLimits exact;
    std::array<CaseFlowComponentBudget, 3> position, velocity;
    std::array<CaseFlowComponentBudget, 3> chord_tangent, span_tangent, normal;
    std::array<CaseFlowComponentBudget, 24> motion_pose;
    CaseFlowComponentBudget core_radius, time;
};

struct CaseWingLimits {
    structure::PgReferenceLimits chart;
    CaseMotionLimits motion;
    CaseFlowClockLimits clock;
    std::uint64_t clock_event_work;
    aero::detail::GeometryLimits geometry;
    aero::detail::GlobalLimits global;
};

struct CaseWingPoint {
    aero::detail::PointId original{0}, canonical{0};
    std::uint64_t surface{}, chord{}, span{}, owner_chord{}, owner_span{}, last_incident_span{},
        strip{};
    structure::PgPartitionPoint chart;
    std::array<io::SiEnclosure, 3> reference, position, velocity;
};

struct CaseWingCell {
    std::uint64_t surface{}, chord{}, span{}, strip{};
    bool reversed{};
    aero::detail::PhysicalFaceId face{0};
    std::array<aero::detail::FacetId, 2> facets{aero::detail::FacetId{0}, aero::detail::FacetId{0}};
    structure::PgPartitionPoint chart;
    std::array<io::SiEnclosure, 3> position, velocity;
    // Tangents are partial derivatives of the Case chart. Only the unit normal is oriented.
    std::array<io::SiEnclosure, 3> chord_tangent, span_tangent;
    std::array<io::SiEnclosure, 3> normal; ///< Includes the Case tangent projection errors.
    std::optional<aero::detail::ResolvedQuery> collocation;
};

struct CaseWingSurfaceView {
    std::uint64_t surface, geometry, source_leading, source_trailing;
    bool reversed;
    std::span<const wing_chart::Rational> chord, span;
    std::span<const wing_chart::Knot> knots;
    std::span<const wing_chart::Point> leading, trailing; ///< Original declaration order.
};

struct CaseWingSharing {
    std::uint64_t own_surface{}, own_boundary{}, other_surface{}, other_boundary{};
    bool reversed{};
};

struct CaseWingTrailingBand {
    std::uint64_t surface{}, source_boundary{}, span{};
    std::array<std::uint64_t, 2> point_strips;
    aero::detail::TrailingEdgeRow transfer;
    std::array<wing_chart::Rational, 2> original_parameters; ///< Boundary line ordinal + fraction.
};

/// Request a finite straight closure from the Case, optionally asserting its exact advance.
struct CaseWingStraightWake {
    std::optional<std::span<const wing_chart::Rational, 3>> expected_displacement_metres;
};

struct CaseWingStraightWakeView {
    const std::array<wing_chart::Rational, 3>& displacement_metres;
    const std::array<wing_chart::Rational, 3>& velocity_metres_per_second;
    const io::SiEnclosure& duration_seconds;
};

struct CaseWingWakeBand {
    aero::detail::EdgeId trailing_edge{0};
    aero::detail::PhysicalFaceId bound_face{0}, wake_face{0};
    aero::detail::FacetId facet{0};
    std::array<aero::detail::PointId, 2> front{aero::detail::PointId{0}, aero::detail::PointId{0}};
    std::array<aero::detail::PointId, 2> terminal{aero::detail::PointId{0},
                                                  aero::detail::PointId{0}};
};

namespace wing_geometry {
struct Owner;
struct Build;
} // namespace wing_geometry

inline constexpr std::uint64_t case_wing_pairing_version = 1;

/// A bound lattice and its global contact qualification, owned by one original SI Case and clock.
class CaseWingGeometry {
public:
    CaseWingGeometry(CaseWingGeometry&&) noexcept;
    CaseWingGeometry& operator=(CaseWingGeometry&&) noexcept;
    CaseWingGeometry(const CaseWingGeometry&) = delete;
    CaseWingGeometry& operator=(const CaseWingGeometry&) = delete;
    ~CaseWingGeometry();

    /// A transient Case must supply a qualified clock address; a static Case uses a zero instant.
    [[nodiscard]] static core::Result<CaseWingGeometry>
    prepare(const io::SiCase&, const aero::detail::GeometryContext&,
            std::optional<CaseFlowClockQuery>, const CaseWingProjection&, const CaseWingLimits&,
            aero::detail::GeometryWorkspace&, aero::detail::GlobalWorkspace&,
            aero::detail::pchip::CurveWorkspace&, core::ResourcePlan&,
            CaseWingDiagnostic* = nullptr, std::optional<CaseWingStraightWake> = {}) noexcept;

    /// Evaluate a transient Case stage against this Case's immutable reference
    /// topology. Requalify contacts, point kinematics and collocation, retaining
    /// the same topology owner and trailing-edge transfer across stages.
    /// A finite straight-wake anchor returns Unsupported.
    [[nodiscard]] core::Result<CaseWingGeometry>
    evaluate_stage(const io::SiCase&, const aero::detail::GeometryContext&, CaseFlowClockQuery,
                   const CaseWingProjection&, const CaseWingLimits&,
                   aero::detail::GeometryWorkspace&, aero::detail::GlobalWorkspace&,
                   aero::detail::pchip::CurveWorkspace&, core::ResourcePlan&,
                   CaseWingDiagnostic* = nullptr) const noexcept;

    [[nodiscard]] core::Result<void> validate_for(const io::SiCase&,
                                                  const aero::detail::GeometryContext&,
                                                  std::optional<CaseFlowClockQuery>) const noexcept;

    [[nodiscard]] core::Result<io::SchemaIdentity> source() const noexcept;
    [[nodiscard]] std::optional<CaseFlowClockQuery> clock_address() const noexcept;
    [[nodiscard]] std::size_t surface_count() const noexcept;
    [[nodiscard]] core::Result<CaseWingSurfaceView> surface(std::size_t) const& noexcept;
    [[nodiscard]] core::Result<CaseWingSurfaceView> surface(std::size_t) const&& = delete;
    [[nodiscard]] std::span<const CaseWingSharing> sharing() const& noexcept;
    [[nodiscard]] std::span<const CaseWingSharing> sharing() const&& = delete;
    [[nodiscard]] std::span<const CaseWingPoint> points() const& noexcept;
    [[nodiscard]] std::span<const CaseWingPoint> points() const&& = delete;
    [[nodiscard]] core::Result<DiscreteCollocation> collocation(std::size_t) const& noexcept;
    [[nodiscard]] core::Result<DiscreteCollocation> collocation(std::size_t) const&& = delete;
    [[nodiscard]] std::span<const CaseWingCell> cells() const& noexcept;
    [[nodiscard]] std::span<const CaseWingCell> cells() const&& = delete;
    [[nodiscard]] std::span<const CaseWingTrailingBand> trailing_bands() const& noexcept;
    [[nodiscard]] std::span<const CaseWingTrailingBand> trailing_bands() const&& = delete;
    [[nodiscard]] std::optional<CaseWingStraightWakeView> straight_wake() const& noexcept;
    [[nodiscard]] std::optional<CaseWingStraightWakeView> straight_wake() const&& = delete;
    [[nodiscard]] std::span<const CaseWingWakeBand> wake_bands() const& noexcept;
    [[nodiscard]] std::span<const CaseWingWakeBand> wake_bands() const&& = delete;
    [[nodiscard]] core::Result<std::reference_wrapper<const aero::detail::GlobalContactState>>
    global() const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const aero::detail::GlobalContactState>>
    global() const&& = delete;

    [[nodiscard]] core::Result<std::reference_wrapper<const aero::detail::AliasKinematicsState>>
    kinematics() const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const aero::detail::AliasKinematicsState>>
    kinematics() const&& = delete;

private:
    friend struct wing_geometry::Build;
    explicit CaseWingGeometry(core::PlannedBuffer<wing_geometry::Owner>) noexcept;
    core::PlannedBuffer<wing_geometry::Owner> owner_;
};

} // namespace uvlas::coupling::detail
