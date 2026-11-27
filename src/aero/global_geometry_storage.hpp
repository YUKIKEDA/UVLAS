#pragma once
#include "geometry_snapshot_storage.hpp"
#include "global_geometry.hpp"

namespace uvlas::aero::detail {
namespace global_detail {
using Vector = std::array<ExactInteger, 3>;
/// Homogeneous point with a strictly positive denominator; never rounded to binary64.
struct RationalPoint {
    Vector numerator;
    ExactInteger denominator = 1;
};
/// Fixed original-point ordinals for the six closed bounds; no giant integers per node.
struct Node {
    std::array<core::Index, 6> extrema{};
    std::uint64_t left = 0, right = 0, count = 0, facet = 0;
};
struct NodePair {
    std::uint64_t first = 0, second = 0;
};
struct BoundEdge {
    std::uint64_t original = 0, owner = 0;
    std::int8_t raw_basis = 0;
    bool seen = false;
    bool canonical = false;
};
struct PackedExact {
    std::uint64_t offset = 0;
    std::uint16_t first = 0, count = 0;
    std::int8_t sign = 0;
};
struct StoredDof {
    std::uint64_t id = 0, basis = 0;
    numerics::linear::UnitDimensions units{};
    std::array<ParameterStorage, 3> basis_vector{};
    ParameterStorage velocity{};
};
struct StoredJacobian {
    std::uint64_t column = 0;
    ParameterStorage value{};
};
struct StoredPointKinematics {
    std::uint64_t id = 0, owner = 0;
    std::array<std::uint64_t, 3> start{}, length{};
    std::array<ParameterStorage, 3> time{};
    std::array<PackedExact, 3> velocity{}, assertion{};
    bool asserted = false;
};
} // namespace global_detail

/// Prepared candidate/product arena; scalar temporaries have bounded fixed stack size.
struct GlobalScratch {
    std::array<global_detail::RationalPoint, 24> candidates;
    std::array<global_detail::RationalPoint, 4> cuts;
    std::array<std::uint32_t, 1024> product_left{}, product_right{}, product_output{};
    std::array<std::uint32_t, 256> packed{};
};
struct GlobalTopologyData {
    static constexpr std::uint8_t bound_role = 1;
    std::optional<OwnedGeometry> geometry;
    core::PlannedBuffer<std::uint8_t> roles;
    core::PlannedBuffer<global_detail::BoundEdge> edges;
    core::PlannedBuffer<global_detail::Node> nodes;
    core::PlannedBuffer<std::uint64_t> order, trailing_points;
    std::uint64_t trailing_point_size = 0;
    core::PlannedBuffer<TrailingEdgeRow> trailing;
    GlobalWorkPlan plan;
};
/// Allocator control block owns a separately planned object payload.
struct GlobalTopologyOwner {
    core::PlannedBuffer<GlobalTopologyData> payload;
    explicit GlobalTopologyOwner(core::PlannedBuffer<GlobalTopologyData> p) noexcept
        : payload(std::move(p)) {}
};
struct GlobalContactData {
    std::optional<OwnedGeometry> geometry;
    std::shared_ptr<const GlobalTopologyData> topology;
    GlobalWorkPlan plan;
};
struct AliasKinematicsData {
    std::optional<OwnedGeometry> geometry;
    std::uint64_t layout_generation = 0;
    core::PlannedBuffer<global_detail::StoredDof> dofs;
    core::PlannedBuffer<global_detail::StoredPointKinematics> points;
    core::PlannedBuffer<global_detail::StoredJacobian> coefficients;
    core::PlannedBuffer<std::uint32_t> limbs;
    GlobalWorkPlan plan;
};
} // namespace uvlas::aero::detail
