#pragma once

#include "discrete_nonpenetration.hpp"

#include <memory>
#include <uvlas/structure/pg_boundary_tangents.hpp>
#include <uvlas/structure/pg_material_tangents.hpp>

namespace uvlas::coupling::detail {

/// Physical face and its fixed structural collocation. A detached point uses
/// the actual PG coordinates. A source expression, when supplied, must agree
/// with the same PG position on the contact-qualified finite facet.
struct PgAeroCollocation {
    aero::detail::FacetId facet{0};
    aero::detail::QueryId query{0};
    std::optional<aero::detail::SourceFacePoint> incidence;
    bool reverse_orientation = false;
};

class PgAeroMap {
public:
    PgAeroMap(PgAeroMap&&) noexcept = default;
    PgAeroMap& operator=(PgAeroMap&&) noexcept = default;
    PgAeroMap(const PgAeroMap&) = delete;
    PgAeroMap& operator=(const PgAeroMap&) = delete;

    /// Material points are ordered as canonical vertices followed by collocations.
    /// An aerodynamic alias always uses its owner's row.
    [[nodiscard]] static core::Result<PgAeroMap>
    prepare(const structure::PgPreparedPatch&, std::span<const aero::detail::PointId>,
            std::span<const PgAeroCollocation>, std::span<const structure::PgBoundaryPositionPoint>,
            const structure::PgPreparedReferenceProjection&, double maximum_tangent_error,
            const structure::PgReferenceLimits&, core::ResourcePlan&) noexcept;

    /// The role order is the same. Material points may be within a layer.
    /// Directions apply only to collocations and are exact in (s,t,zeta).
    [[nodiscard]] static core::Result<PgAeroMap>
    prepare_material(const structure::PgPreparedPatch&, std::span<const aero::detail::PointId>,
                     std::span<const PgAeroCollocation>,
                     std::span<const structure::PgReferencePoint>,
                     std::span<const structure::PgMaterialDirections>,
                     const structure::PgPreparedReferenceProjection&, double maximum_tangent_error,
                     const structure::PgReferenceLimits&, core::ResourcePlan&) noexcept;

    [[nodiscard]] structure::PgPreparedPatchIdentity source() const noexcept;

private:
    friend struct PgAeroBuilder;
    friend struct PgAeroTrialData;
    struct Data;
    explicit PgAeroMap(std::shared_ptr<const Data>) noexcept;
    template <class Position, class Tangents>
    [[nodiscard]] static core::Result<PgAeroMap>
    prepare_impl(std::span<const aero::detail::PointId>, std::span<const PgAeroCollocation>,
                 std::size_t material_count, const structure::PgReferenceLimits&,
                 core::ResourcePlan&, const Position&, const Tangents&) noexcept;
    std::shared_ptr<const Data> data_;
};

struct PgAeroState {
    structure::PgPreparedPatchIdentity patch;
    aero::detail::GeometryContext context;
    structure::PgCartesianWords position, velocity;
    aero::detail::P1ActionStamp stamp;
};

struct PgAeroLimits {
    /// Separate dimensional projection limits, for metres and metres per second.
    structure::PgVolumePositionLimits position, velocity, tangents;
    aero::detail::GeometryLimits geometry;
    aero::detail::GlobalLimits contact;
    std::uint64_t maximum_steps{};
    double maximum_midpoint_velocity_error{};
    double maximum_reference_geometry_error{}; ///< Euclidean displacement, in metres.
};

struct PgAeroChart {
    const aero::detail::GlobalTopology& topology;
    const aero::detail::GeometryInput& geometry;
    std::span<const aero::detail::EdgeId> load_edges; ///< Ascending canonical non-TE edges.
};

struct PgAeroTrialData;

class PgAeroTrial {
public:
    PgAeroTrial(PgAeroTrial&&) noexcept;
    PgAeroTrial& operator=(PgAeroTrial&&) noexcept;
    PgAeroTrial(const PgAeroTrial&) = delete;
    PgAeroTrial& operator=(const PgAeroTrial&) = delete;
    ~PgAeroTrial();

    [[nodiscard]] core::Result<std::reference_wrapper<const aero::detail::GlobalContactState>>
    contact() const& noexcept;
    [[nodiscard]] core::Result<std::reference_wrapper<const aero::detail::GlobalContactState>>
    contact() const&& = delete;
    [[nodiscard]] std::size_t collocation_count() const noexcept;
    [[nodiscard]] core::Result<DiscreteCollocation> collocation(std::size_t) const& noexcept;
    [[nodiscard]] core::Result<DiscreteCollocation> collocation(std::size_t) const&& = delete;
    [[nodiscard]] std::size_t midpoint_count() const noexcept;
    [[nodiscard]] core::Result<EdgeMidpointVelocity> midpoint(std::size_t) const noexcept;

    /// Scatter canonical assembled point forces into the retained PG column order.
    /// Project the exact force enclosures without dropping either words or radii.
    [[nodiscard]] core::Result<core::PlannedBuffer<std::optional<DiscreteLoadVector>>>
    pullback_nodal_forces(const DiscreteCirculationLoads&, core::Precision,
                          double maximum_absolute_error, std::uint64_t maximum_steps,
                          core::ResourcePlan&, const std::stop_token& = {}) const noexcept;

private:
    friend struct PgAeroBuilder;
    explicit PgAeroTrial(core::PlannedBuffer<PgAeroTrialData>&&) noexcept;
    core::PlannedBuffer<PgAeroTrialData> data_;
};

[[nodiscard]] core::Result<PgAeroTrial>
evaluate_pg_aero_trial(const PgAeroMap&, const PgAeroState&, const PgAeroChart&,
                       const PgAeroLimits&, core::Workspace&, aero::detail::GeometryWorkspace&,
                       aero::detail::GlobalWorkspace&, core::ResourcePlan&,
                       const std::stop_token& = {}) noexcept;

} // namespace uvlas::coupling::detail
