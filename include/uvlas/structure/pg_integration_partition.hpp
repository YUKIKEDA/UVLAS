#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <stop_token>
#include <string_view>
#include <uvlas/core/identifier.hpp>
#include <uvlas/core/state.hpp>
#include <uvlas/structure/pg_graph.hpp>
#include <variant>

namespace uvlas::structure {
/// Owned exact dimensionless chart coordinates. Neither coordinate is physical xyz.
struct PgPartitionPoint {
    ExactRational s; ///< First chart coordinate.
    ExactRational t; ///< Second chart coordinate.
};
/// Borrowed continuous piecewise-linear physical full thickness, in metres.
struct PgPartitionThickness {
    std::span<const PgPartitionPoint> vertices;            ///< Chart points, in caller order.
    std::span<const std::array<core::Index, 3>> triangles; ///< Vertex indices, either orientation.
    std::span<const ExactRational> values; ///< Positive full thickness, one per vertex.
};
/// One material interval. Ordered intervals must exactly cover [-1,1].
struct PgPartitionLayer {
    /// Borrow both exact zeta endpoints and preserve the explicitly supplied ID.
    PgPartitionLayer(core::LayerId material, const ExactRational& bottom,
                     const ExactRational& top) noexcept
        : id(material), lower(bottom), upper(top) {}
    core::LayerId id;           ///< Explicit layer identity, independent of array index.
    const ExactRational& lower; ///< Dimensionless lower zeta coordinate.
    const ExactRational& upper; ///< Dimensionless upper zeta coordinate.
};
/// A named, ordered subset of the actual polygon boundary; no physical tolerance.
struct PgPartitionBoundary {
    std::string_view name;                      ///< Nonempty UTF-8 name copied into the result.
    std::span<const PgPartitionPoint> vertices; ///< At least two ordered chart points.
};
/// An additional partition cut. Extending its supporting line creates no outer face.
struct PgPartitionRefinement {
    /// Borrow two explicit endpoints and preserve cut identity and revision.
    PgPartitionRefinement(std::uint64_t cut_id, std::uint64_t revision, const PgPartitionPoint& a,
                          const PgPartitionPoint& b) noexcept
        : id(cut_id), version(revision), first(a), second(b) {}
    std::uint64_t id;               ///< Explicit caller cut identity; zero is valid.
    std::uint64_t version;          ///< Explicit version of that cut definition.
    const PgPartitionPoint& first;  ///< First chart endpoint.
    const PgPartitionPoint& second; ///< Distinct second chart endpoint.
};
/// Read-only inputs borrowed only during the synchronous preparation call.
struct PgPartitionRequest {
    const PgGraph& graph;                     ///< Existing immutable graph and n=1/2/4 cells.
    std::span<const PgPartitionPoint> domain; ///< Simple positive-area polygon in the unit square.
    PgPartitionThickness thickness; ///< Continuous positive full-thickness mesh over this domain.
    std::span<const PgPartitionLayer> layers; ///< Ordered, nonempty full-thickness partition.
    std::span<const PgPartitionBoundary> boundaries;    ///< Named physical boundary subsets.
    std::span<const PgPartitionRefinement> refinements; ///< Internal computational cuts.
};
/// Explicit whole-call limits; zero never means unlimited. Input counts and
/// conservative output capacities are checked before allocation, not clipped.
struct PgPartitionLimits {
    ExactLimits algebra; ///< Exact scalar value/work bits and external-memory policy.
    std::uint64_t maximum_input_points; ///< All point occurrences, including cut endpoints.
    std::uint64_t maximum_segments;     ///< All polygon, thickness, named and refinement segments.
    std::uint64_t maximum_lines;  ///< Supporting-line capacity including auxiliary/grid/box lines.
    std::uint64_t maximum_pieces; ///< Conservative positive-area-piece capacity.
    std::uint64_t maximum_vertices;           ///< Conservative line-pair intersection capacity.
    std::uint64_t maximum_edges;              ///< Conservative atomic edge-occurrence capacity.
    std::uint64_t maximum_layers;             ///< Material interval count.
    std::uint64_t maximum_numeric_operations; ///< Scalar opcodes, separate from topology visits.
    std::uint64_t maximum_topology_visits;    ///< Record/loop/copy work; text bytes count here.
    std::uint64_t maximum_text_bytes;         ///< Owned names, excluding fixed numeric IDs.
    std::stop_token cancellation; ///< Checked before publication and within bounded loops.
};
/// Chart partition content, not a complete physical reference qualification.
struct PgPartitionIdentity {
    std::array<std::byte, 32> bytes; ///< Canonical SHA-256, independent of addresses and limb ABI.
    bool operator==(const PgPartitionIdentity&) const = default;
};
/// Optional nonoverlapping preparation phases; whole includes plan/allocation overhead.
struct PgPartitionTimings {
    core::StageTiming whole;         ///< Entire successful synchronous preparation.
    core::StageTiming input;         ///< Owned source clone and original geometry validation.
    core::StageTiming intersections; ///< Canonical original lines and line-pair vertices.
    core::StageTiming planar;        ///< Convex common refinement and cell/domain classification.
    core::StageTiming traces; ///< Shared edges, H continuity, components, all faces and hash.
};
/// Whole-call mathematical and allocation forecasts, independent of live heap/RSS.
struct PgPartitionPlan {
    std::uint64_t input_points = 0;          ///< All borrowed point occurrences.
    std::uint64_t segments = 0;              ///< All original segments with provenance.
    std::uint64_t line_capacity = 0;         ///< Before duplicate-support removal.
    std::uint64_t intersection_capacity = 0; ///< All line pairs.
    std::uint64_t piece_capacity = 0;       ///< Full square arrangement, including inactive pieces.
    std::uint64_t clip_vertex_capacity = 0; ///< One intermediate convex-vertex buffer.
    std::uint64_t edge_capacity = 0;        ///< Both-side final vertex/edge occurrences.
    std::uint64_t numeric_operations = 0;   ///< Conservative opcode count, not elapsed time.
    std::uint64_t topology_visits = 0;      ///< Conservative topology/copy operation count.
    std::uint64_t value_bits = 0;           ///< Bound for retained exact values.
    std::uint64_t algebraic_work_bits = 0;  ///< Before guard/word rounding.
    std::uint64_t planned_work_bits = 0;    ///< Guarded temporary integer capacity.
    std::uint64_t owned_payload_bytes =
        0; ///< Cumulative planned allocations, including scalar storage.
    std::uint64_t estimated_external_bytes =
        0;                        ///< Fixed-Boost limb-payload forecast, not hard cap.
    std::uint64_t text_bytes = 0; ///< Exact UTF-8 name-copy bytes.
    std::uint64_t allocated_owned_payload_bytes =
        0; ///< Observed successful plan delta; zero in forecasts.
    std::uint64_t observed_numeric_operations = 0; ///< Actually executed scalar opcodes.
    std::uint64_t observed_topology_visits = 0;    ///< Actually visited/copy/hash work.
    std::uint64_t observed_value_bits = 0; ///< Actual maximum scalar numerator/denominator bits.
    std::optional<PgPartitionTimings>
        timings; ///< Present only when explicitly enabled; forecasts omit it.
};
/// Original polygon edge identity; no unused material/cut fields are present.
struct PgPartitionDomainSource {
    core::Index edge = 0; ///< Original edge ordinal.
};
/// Original thickness triangle edge identity.
struct PgPartitionThicknessSource {
    core::Index triangle = 0; ///< Original triangle ordinal.
    core::Index edge = 0;     ///< Local edge ordinal, in [0,3).
};
/// Original named subset segment identity.
struct PgPartitionNamedSource {
    core::Index boundary = 0; ///< Original named boundary ordinal.
    core::Index segment = 0;  ///< Original segment ordinal within that boundary.
};
/// Computational cut identity, retaining every explicitly supplied field.
struct PgPartitionRefinementSource {
    core::Index ordinal = 0;   ///< Original cut ordinal.
    std::uint64_t id = 0;      ///< Original cut ID; zero is valid.
    std::uint64_t version = 0; ///< Original cut revision.
};
/// Closed provenance alternatives; absence of a field is distinct from ID zero.
using PgPartitionSourceDescriptor =
    std::variant<PgPartitionDomainSource, PgPartitionThicknessSource, PgPartitionNamedSource,
                 PgPartitionRefinementSource>;
/// Owned source descriptor; endpoint indices address the original-point array.
struct PgPartitionSource {
    PgPartitionSourceDescriptor descriptor; ///< Exact source-specific information.
    std::array<core::Index, 2> endpoints{}; ///< Original directed chart endpoints.
};
/// Exact positive chart area, with cell-major association and no rank claim.
struct PgPartitionCell {
    PgCellIndex cell{}; ///< Every background cell appears exactly once in standard order.
    ExactRational area; ///< Area of actual domain intersection, zero allowed.
    std::uint64_t components =
        0;               ///< Interior-connected components; point contacts do not join them.
    bool active = false; ///< Derived only from exact positive area by the private factory.
};
/// One positive convex piece; offsets index the returned flattened vertex array.
struct PgPartitionPiece {
    PgCellIndex cell{};                 ///< Unique background cell containing this piece interior.
    core::Index component = 0;          ///< Cell-local component, ordered by least piece index.
    core::Index thickness_triangle = 0; ///< Unique full-thickness mesh triangle on its interior.
    core::Index first_vertex = 0;       ///< Offset into piece_vertex_indices().
    core::Index vertex_count = 0;       ///< CCW polygon, starting at least VertexIndex.
    ExactRational area;                 ///< Exact chart area, strictly positive.
};
/// Owned layer coordinates; the result retains both ID and exact endpoints.
struct PgPartitionOwnedLayer {
    core::LayerId id{0}; ///< Original material-layer identity.
    ExactRational lower; ///< Dimensionless lower zeta.
    ExactRational upper; ///< Dimensionless upper zeta.
};
/// Owned named-boundary metadata, with bytes stored in boundary_text().
struct PgPartitionBoundaryRecord {
    core::Index text_offset = 0; ///< UTF-8 byte offset, no hidden terminator.
    core::Index text_count = 0;  ///< Name length in bytes.
    core::Index first_point = 0; ///< Offset in source_points().
    core::Index point_count = 0; ///< Original ordered boundary points.
};
/// Borrowed snapshot of the original thickness triangulation, owned by the result.
struct PgPartitionThicknessView {
    std::span<const PgPartitionPoint> vertices;            ///< Original chart vertices.
    std::span<const std::array<core::Index, 3>> triangles; ///< Original connectivity and order.
    std::span<const ExactRational> values; ///< Original positive full thickness in metres.
};
/// One shared atomic segment. Its canonical direction is lower to higher vertex index.
struct PgPartitionEdge {
    std::array<core::Index, 2> vertices{}; ///< Indices into vertices().
    std::optional<core::Index> left;       ///< Active piece left of the canonical tangent.
    std::optional<core::Index> right;      ///< Active piece right of the canonical tangent.
    core::Index first_source = 0;          ///< Offset into edge_source_indices().
    core::Index source_count = 0; ///< All exact containing original segments, in input order.
    bool external = false; ///< Exactly one active side AND original-domain boundary provenance.
};
/// Positive triangle from a deterministic fan of the owning planar piece.
struct PgPartitionTriangle {
    core::Index piece = 0;                 ///< Owning planar piece.
    std::array<core::Index, 3> vertices{}; ///< Positive chart orientation; uses common vertices().
    ExactRational twice_area; ///< Exact positive chart Jacobian on the reference triangle.
};
/// Triangle x normalized-thickness interval, without physical quadrature weights.
struct PgPartitionVolume {
    core::Index triangle = 0; ///< Common planar triangle index.
    core::Index layer = 0;    ///< Owned interval index, preserving its LayerId.
};
/// Material side of a horizontal face; interface records remain distinct.
enum class PgPartitionFaceSide : std::uint8_t {
    Lower, ///< Outward direction is decreasing zeta.
    Upper  ///< Outward direction is increasing zeta.
};
/// One lower/upper material trace using the same planar triangle as its volume.
struct PgPartitionHorizontalFace {
    core::Index triangle = 0; ///< Common planar triangle index.
    core::Index layer = 0;    ///< Material interval, not a merged interface.
    PgPartitionFaceSide side = PgPartitionFaceSide::Lower; ///< Outward zeta sign: lower negative.
};
/// Physical outer segment x layer. Geometry consumers form the physical normal.
struct PgPartitionLateralFace {
    core::Index edge = 0;            ///< An external atomic edge.
    core::Index layer = 0;           ///< Same interval as volume/horizontal faces.
    bool canonical_direction = true; ///< True when lower-to-higher vertex tangent is CCW domain.
};
/// One named original boundary trace, including original orientation and adjacent cell.
struct PgPartitionEdgeTrace {
    core::Index edge = 0;   ///< Shared external edge geometry; not a duplicate physical boundary.
    core::Index source = 0; ///< NamedBoundary source descriptor.
    core::Index piece = 0;  ///< Interior-adjacent piece, including cell and thickness association.
    bool canonical_direction = true; ///< Canonical tangent agrees with the named segment direction.
};
/// A closed-cell zero-measure intersection with an original domain vertex/segment.
struct PgPartitionContact {
    PgCellIndex cell{};     ///< Background cell whose open-area intersection is zero.
    core::Index source = 0; ///< Original domain source, not a caller-provided contact flag.
    std::array<core::Index, 2> vertices{}; ///< Equal endpoints mean a point contact.
};
namespace detail {
struct PgPartitionData;
struct PgPartitionFactory;
} // namespace detail

/// Immutable owned exact planar/through-thickness partition. It qualifies
/// activity only, not rank, physical Jacobians, quadrature accuracy or a Model.
class PgIntegrationPartition {
public:
    /// Transfer ownership; queries on moved-from owners return InvalidInput.
    PgIntegrationPartition(PgIntegrationPartition&&) noexcept;
    /// Transfer storage, invalidating every prior borrowed view of this owner.
    PgIntegrationPartition& operator=(PgIntegrationPartition&&) noexcept;
    /// Copying cannot bypass the explicit preparation budget.
    PgIntegrationPartition(const PgIntegrationPartition&) = delete;
    /// Implicit copies of owned geometry and exact values are forbidden.
    PgIntegrationPartition& operator=(const PgIntegrationPartition&) = delete;
    /// Release all owned arrays; the original ResourcePlan is not required.
    ~PgIntegrationPartition();
    /// Make an allocation-free conservative whole-call forecast. Invalid shapes,
    /// unsupported hard-memory requests, overflow, cancellation and quota failures
    /// return Failure values. No input, workspace or ResourcePlan is changed.
    [[nodiscard]] static core::Result<PgPartitionPlan> plan(const PgPartitionRequest&,
                                                            const PgPartitionLimits&) noexcept;
    /// Own the complete exact common refinement and all material/boundary traces.
    /// Holds an exclusive workspace lease for the synchronous call; separate
    /// workspaces/plans can prepare the same immutable inputs in parallel.
    /// Invalid input/geometry, stale associations, resource limits and cancellation
    /// fail without changing input or an earlier result. Exceptions never escape.
    [[nodiscard]] static core::Result<PgIntegrationPartition> prepare(const PgPartitionRequest&,
                                                                      const PgPartitionLimits&,
                                                                      core::Workspace&,
                                                                      core::ResourcePlan&) noexcept;
    /// Compare actual source content, including graph, domain, thickness, layers,
    /// names and refinements; a digest alone never proves compatibility.
    [[nodiscard]] core::Result<void> validate_for(const PgPartitionRequest&) const noexcept;
    /// Borrow the independently owned graph snapshot; this is not a physical Model.
    [[nodiscard]] core::Result<PgGraphView> graph() const& noexcept;
    [[nodiscard]] core::Result<PgGraphView> graph() const&& = delete;
    /// Borrow immutable exact vertices. All views expire when this owner moves or dies.
    [[nodiscard]] core::Result<std::span<const PgPartitionPoint>> vertices() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgPartitionPoint>> vertices() const&& = delete;
    /// Borrow reporting metadata, including inactive/zero-area cells. This span
    /// cannot confer activity qualification; consumers must take the factory-only
    /// PgIntegrationPartition and validate its actual source association.
    [[nodiscard]] core::Result<std::span<const PgPartitionCell>> activity_map() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgPartitionCell>> activity_map() const&& = delete;
    /// Borrow deterministic positive pieces; all topology indices address this result.
    [[nodiscard]] core::Result<std::span<const PgPartitionPiece>> pieces() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgPartitionPiece>> pieces() const&& = delete;
    /// Borrow each piece's flattened CCW vertex list.
    [[nodiscard]] core::Result<std::span<const core::Index>> piece_vertex_indices() const& noexcept;
    [[nodiscard]] core::Result<std::span<const core::Index>>
    piece_vertex_indices() const&& = delete;
    /// Borrow atomic shared/external edge topology and physical provenance.
    [[nodiscard]] core::Result<std::span<const PgPartitionEdge>> edges() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgPartitionEdge>> edges() const&& = delete;
    /// Borrow flattened edge-to-original-source indices.
    [[nodiscard]] core::Result<std::span<const core::Index>> edge_source_indices() const& noexcept;
    [[nodiscard]] core::Result<std::span<const core::Index>> edge_source_indices() const&& = delete;
    /// Borrow the retained source segment descriptors, preserving original direction.
    [[nodiscard]] core::Result<std::span<const PgPartitionSource>> sources() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgPartitionSource>> sources() const&& = delete;
    /// Borrow the original owned chart endpoint array used by source descriptors.
    [[nodiscard]] core::Result<std::span<const PgPartitionPoint>> source_points() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgPartitionPoint>> source_points() const&& = delete;
    /// Borrow retained full-thickness values and original triangle connectivity.
    [[nodiscard]] core::Result<PgPartitionThicknessView> thickness() const& noexcept;
    [[nodiscard]] core::Result<PgPartitionThicknessView> thickness() const&& = delete;
    /// Borrow exact owned intervals and material IDs, including both endpoints.
    [[nodiscard]] core::Result<std::span<const PgPartitionOwnedLayer>> layers() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgPartitionOwnedLayer>> layers() const&& = delete;
    /// Borrow named-boundary metadata; offsets address this result's source/text arrays.
    [[nodiscard]] core::Result<std::span<const PgPartitionBoundaryRecord>>
    boundaries() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgPartitionBoundaryRecord>>
    boundaries() const&& = delete;
    /// Borrow original UTF-8 boundary names without null terminators.
    [[nodiscard]] core::Result<std::span<const char>> boundary_text() const& noexcept;
    [[nodiscard]] core::Result<std::span<const char>> boundary_text() const&& = delete;
    /// Borrow the positive common planar integration triangles.
    [[nodiscard]] core::Result<std::span<const PgPartitionTriangle>> triangles() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgPartitionTriangle>> triangles() const&& = delete;
    /// Borrow triangle x layer volume domains, without physical weights.
    [[nodiscard]] core::Result<std::span<const PgPartitionVolume>> volumes() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgPartitionVolume>> volumes() const&& = delete;
    /// Borrow all external and both material-interface horizontal traces.
    [[nodiscard]] core::Result<std::span<const PgPartitionHorizontalFace>>
    horizontal_faces() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgPartitionHorizontalFace>>
    horizontal_faces() const&& = delete;
    /// Borrow only physical outer lateral faces; auxiliary cuts have none.
    [[nodiscard]] core::Result<std::span<const PgPartitionLateralFace>>
    lateral_faces() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgPartitionLateralFace>>
    lateral_faces() const&& = delete;
    /// Borrow original named-boundary subsets, not additional outer faces.
    [[nodiscard]] core::Result<std::span<const PgPartitionEdgeTrace>> edge_traces() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgPartitionEdgeTrace>>
    edge_traces() const&& = delete;
    /// Borrow point/segment contacts in zero-area cells.
    [[nodiscard]] core::Result<std::span<const PgPartitionContact>> contacts() const& noexcept;
    [[nodiscard]] core::Result<std::span<const PgPartitionContact>> contacts() const&& = delete;
    /// Return content identity or InvalidInput for a moved-from result.
    [[nodiscard]] core::Result<PgPartitionIdentity> identity() const noexcept;
    /// Borrow the whole-call forecast, independent of the original plan lifetime.
    [[nodiscard]] const PgPartitionPlan& preparation() const& noexcept;
    [[nodiscard]] const PgPartitionPlan& preparation() const&& = delete;

private:
    friend struct detail::PgPartitionFactory;
    explicit PgIntegrationPartition(core::PlannedBuffer<detail::PgPartitionData>,
                                    PgPartitionPlan) noexcept;
    core::PlannedBuffer<detail::PgPartitionData> data_;
    PgPartitionPlan info_;
};
} // namespace uvlas::structure
