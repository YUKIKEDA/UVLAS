#pragma once

#include "geometry_snapshot.hpp"

#include <atomic>
#include <functional>
#include <uvlas/numerics/linear_solve.hpp>

namespace uvlas::aero::detail {
namespace global_tags {
struct Epoch;
struct Layout;
struct Basis;
} // namespace global_tags
/// Epoch of a material release; all uint64 values are valid.
using ReleaseEpoch = core::Identifier<global_tags::Epoch>;
/// Revision of the complete physical degree-of-freedom layout.
using DofLayoutGeneration = core::Identifier<global_tags::Layout>;
/// Model-defined identity of a physical component basis.
using DofBasisId = core::Identifier<global_tags::Basis>;
/// Role in the bound subcomplex; a wake face does not count as a second bound owner.
enum class AeroFaceRole : std::uint8_t { Bound, Wake };
/// Explicit role of one physical face. Every retained face needs one entry.
struct FaceRoleInput {
    PhysicalFaceId face; ///< Existing physical face, never a child-facet ordinal.
    AeroFaceRole role;   ///< Bound or wake classification.
};
/// This edge is the canonical owner of its explicitly bound aliases.
struct CanonicalEdge {};
/// An edge owns itself, or directly names a canonical edge in the same snapshot.
using EdgeOwnership = std::variant<CanonicalEdge, EdgeId>;
/// Edge aliasing is explicit and never inferred from coincident endpoints.
struct EdgeBindingInput {
    EdgeId edge;             ///< Original edge identity.
    EdgeOwnership ownership; ///< CanonicalEdge means canonical; EdgeId names its owner.
    /// Retain the explicit edge and ownership declaration.
    EdgeBindingInput(EdgeId e, EdgeOwnership o) noexcept : edge(e), ownership(o) {}
};
/// Wake-induction traversal along a declared trailing edge, using original point IDs.
struct TrailingEdgeInput {
    EdgeId edge;                   ///< Canonical or explicitly aliased edge.
    std::array<PointId, 2> points; ///< Ordered wake-induction endpoints, not geometry boundary.
};
/// Complete immutable-model inputs for global topology construction.
struct GlobalTopologyInput {
    std::span<const FaceRoleInput> roles;              ///< One role per physical face.
    std::span<const EdgeBindingInput> edges;           ///< One binding per original edge.
    std::span<const TrailingEdgeInput> trailing_edges; ///< Explicit directed TE set.
};
/// Caller-selected work and opaque ownership-metadata admission limits.
struct GlobalLimits {
    std::uint64_t max_steps =
        0; ///< Whole-call fixed-capacity primitive units; zero is a real limit.
    numerics::linear::BackendMemoryLimit metadata; ///< Explicit forecast/proven-bound choice.
};
/// A global failure retains its first entity and, for contacts, the second facet.
struct GlobalGeometryFailure {
    OwnedGeometryFailure origin;        ///< Original cause, details, entity and context.
    std::optional<FacetId> other_facet; ///< Second facet in a rejected unordered pair.
};
/// A complete owning qualification or its preserved diagnostic; no partial success.
template <class T>
using GlobalResult = std::expected<T, GlobalGeometryFailure>;
/// Static admitted and observed counters, excluding opaque allocator metadata.
struct GlobalWorkPlan {
    std::uint64_t payload_bytes = 0;  ///< All successful planned payload allocations in this call.
    std::uint64_t reserved_steps = 0; ///< Checked whole-call work admission.
    std::uint64_t executed_steps = 0; ///< Charged fixed-limb units, not processor instructions.
    std::uint64_t candidates = 0;     ///< Closed-box leaf pairs, including allowed contacts.
    std::uint64_t metadata_forecast = 0; ///< Explicit estimate, never a proven RSS limit.
    std::uint64_t traversal_reserved_steps =
        0; ///< Two complete BVH passes, excluding narrow pairs.
    std::uint64_t traversal_executed_steps =
        0; ///< Actual box/stack/decode work across both passes.
};
/// A sparse integer row of the TE transfer T; the unique bound coefficient is never averaged.
struct TrailingEdgeRow {
    EdgeId edge{0};               ///< Canonical global edge.
    PhysicalFaceId bound_face{0}; ///< Unique bound physical owner.
    std::int8_t coefficient = 0;  ///< T entry in {-1,+1} in a qualified result.
    std::int8_t wake_sign = 0;    ///< Wake induction sign in canonical endpoint order.
};
/// Contact-qualified fixed quadrilateral and its declared leading attachment.
struct UvlmQuad {
    PhysicalFaceId face;
    AeroFaceRole role;
    std::array<FacetId, 2> facets;
    std::optional<TrailingEdgeRow> attachment;
    std::uint64_t used_steps;
};
/// Stable key of a released material point; no lossy packing or hash identity.
struct MaterialReleaseKey {
    PointId point;      ///< Canonical physical TE point.
    ReleaseEpoch epoch; ///< Material release epoch.
    auto operator<=>(const MaterialReleaseKey&) const = default;
};
/// A point is attached or released; equal coordinates cannot convert its history.
using ReleaseIdentity = std::variant<PointId, MaterialReleaseKey>;
/// No result in this module certifies a physical mapping or a continuous deformation path.
struct PhysicalMappingNotChecked {};
class GlobalEngine;
struct GlobalTopologyData;
struct GlobalContactData;
struct AliasKinematicsData;
struct GlobalScratch;
/// Owned, immutable global connectivity, reference BVH and TE transfer.
class GlobalTopology {
public:
    /// Transfer the immutable owner; a moved-from object is inactive.
    GlobalTopology(GlobalTopology&&) noexcept = default;
    /// Transfer ownership without modifying the retained source snapshot.
    GlobalTopology& operator=(GlobalTopology&&) noexcept = default;
    GlobalTopology(const GlobalTopology&) = delete;
    GlobalTopology& operator=(const GlobalTopology&) = delete;
    /// Release this owner; dependent results retain their own ownership.
    ~GlobalTopology();
    /// Whether the object owns a successfully constructed topology.
    [[nodiscard]] bool active() const noexcept;
    /// Read the sparse TE rows in canonical edge-ID order.
    [[nodiscard]] std::span<const TrailingEdgeRow> trailing_edges() const& noexcept;
    [[nodiscard]] std::span<const TrailingEdgeRow> trailing_edges() const&& = delete;
    /// Read canonical TE point IDs in strictly increasing order.
    [[nodiscard]] std::size_t trailing_point_count() const noexcept;
    /// Copy a canonical TE point at its stable sorted ordinal.
    [[nodiscard]] GlobalResult<PointId> trailing_point(std::size_t) const noexcept;
    /// Copy a material key for a declared TE point, without advancing or committing a wake.
    [[nodiscard]] GlobalResult<MaterialReleaseKey> release_key(PointId,
                                                               ReleaseEpoch) const noexcept;
    /// Original snapshot kept alive by this result.
    [[nodiscard]] GlobalResult<std::reference_wrapper<const OwnedGeometry>>
    geometry() const& noexcept;
    [[nodiscard]] GlobalResult<std::reference_wrapper<const OwnedGeometry>>
    geometry() const&& = delete;
    /// Planned payload and work for topology construction.
    [[nodiscard]] GlobalResult<GlobalWorkPlan> work_plan() const noexcept;
    /// No physical-mapping qualification is granted by connectivity checks.
    [[nodiscard]] static constexpr PhysicalMappingNotChecked mapping_status() noexcept {
        return {};
    }

private:
    friend class GlobalEngine;
    explicit GlobalTopology(std::shared_ptr<const GlobalTopologyData>, GeometryContext) noexcept;
    std::shared_ptr<const GlobalTopologyData> data_;
    GeometryContext context_;
};
/// A complete discrete contact check tied to the retained current snapshot and topology.
class GlobalContactState {
public:
    /// Transfer the qualified immutable owner.
    GlobalContactState(GlobalContactState&&) noexcept;
    /// Transfer the qualified owner; the source becomes inactive.
    GlobalContactState& operator=(GlobalContactState&&) noexcept;
    GlobalContactState(const GlobalContactState&) = delete;
    GlobalContactState& operator=(const GlobalContactState&) = delete;
    /// Release this contact qualification only.
    ~GlobalContactState();
    /// Whether all finite facets passed the global contact check.
    [[nodiscard]] bool active() const noexcept;
    /// Current geometry actually checked, independent of the caller's lifetime.
    [[nodiscard]] GlobalResult<std::reference_wrapper<const OwnedGeometry>>
    geometry() const& noexcept;
    [[nodiscard]] GlobalResult<std::reference_wrapper<const OwnedGeometry>>
    geometry() const&& = delete;
    /// Inspect the fixed quad of a source facet. A wake's L edge uses quad vertices 0→1.
    [[nodiscard]] GlobalResult<UvlmQuad> uvlm_quad(FacetId, std::uint64_t max_steps,
                                                   const std::stop_token& = {}) const noexcept;
    [[nodiscard]] std::size_t trailing_edge_count() const noexcept;
    /// Two-pass traversal, narrow-pair and planned-payload counters.
    [[nodiscard]] GlobalResult<GlobalWorkPlan> work_plan() const noexcept;
    /// Contact clearance does not certify physical mapping or state paths.
    [[nodiscard]] static constexpr PhysicalMappingNotChecked mapping_status() noexcept {
        return {};
    }

private:
    friend class GlobalEngine;
    friend class SharedP1Engine;
    explicit GlobalContactState(core::PlannedBuffer<GlobalContactData>, GeometryContext) noexcept;
    core::PlannedBuffer<GlobalContactData> data_;
    GeometryContext context_;
};
/// A common physical scalar DOF and its retained SI unit/basis description.
struct GlobalDofInput {
    core::DofId id; ///< Strictly increasing identity in the common layout.
    numerics::linear::UnitDimensions units{}; ///< Seven SI exponents of q, in existing N4 order.
    DofBasisId basis;                        ///< Identity in the retained model's basis definition.
    std::array<ScalarWords, 3> basis_vector; ///< Exact components in the global Cartesian frame.
    ScalarWords velocity;                    ///< qdot in the selected physical DOF unit per second.
    /// Retain a complete physical DOF description and original velocity words.
    GlobalDofInput(core::DofId i, numerics::linear::UnitDimensions u, DofBasisId b,
                   std::array<ScalarWords, 3> axes, ScalarWords v) noexcept
        : id(i), units(u), basis(b), basis_vector(axes), velocity(v) {}
};
/// One sparse exact coefficient; duplicate or unsorted columns are invalid.
struct GlobalJacobianEntry {
    core::DofId column; ///< Column in the common global layout.
    ScalarWords value;  ///< Exact m/(DOF unit) coefficient; exact zero must be omitted.
    /// Retain an identified coefficient without dropping its low word.
    GlobalJacobianEntry(core::DofId c, ScalarWords v) noexcept : column(c), value(v) {}
};
/// Read-only exact velocity assertion in units 2^-4296 m/s.
struct ExactVelocityInput {
    std::int8_t sign = 0;                 ///< -1/+1 for nonzero; zero requires an empty magnitude.
    std::uint16_t first_limb = 0;         ///< Offset of the lowest retained 32-bit limb.
    std::span<const std::uint32_t> limbs; ///< Little-endian magnitude, nonzero endpoint limbs.
};
/// Full xyz rows of one original point, including explicit time terms.
struct PointKinematicsInput {
    PointId point; ///< Original point ID, including each alias.
    std::array<std::span<const GlobalJacobianEntry>, 3>
        jacobian;                             ///< All global columns per xyz row.
    std::array<ScalarWords, 3> explicit_time; ///< gt in m/s.
    std::optional<std::array<ExactVelocityInput, 3>>
        velocity; ///< Optional exact consistency assertions.
    /// Bind complete xyz rows and time terms to an original point.
    PointKinematicsInput(PointId p, std::array<std::span<const GlobalJacobianEntry>, 3> j,
                         std::array<ScalarWords, 3> t,
                         std::optional<std::array<ExactVelocityInput, 3>> v = {}) noexcept
        : point(p), jacobian(j), explicit_time(t), velocity(v) {}
};
/// Complete data for a full-column alias comparison at one geometry state.
struct AliasKinematicsInput {
    DofLayoutGeneration generation;               ///< Revision of the full common layout.
    std::span<const GlobalDofInput> dofs;         ///< Strictly increasing global DOFs.
    std::span<const PointKinematicsInput> points; ///< Exactly one entry per original point.
    /// Describe the complete rows in one common physical layout.
    AliasKinematicsInput(DofLayoutGeneration g, std::span<const GlobalDofInput> d,
                         std::span<const PointKinematicsInput> p) noexcept
        : generation(g), dofs(d), points(p) {}
};
/// Read-only packed exact result; the containing state owns every limb.
struct ExactVelocityView {
    std::int8_t sign{};                   ///< Zero iff the limb span is empty.
    std::uint16_t first_limb{};           ///< Offset in the fixed 2^-4296 integer grid.
    std::span<const std::uint32_t> limbs; ///< Exact little-endian magnitude, not rounded words.
};
/// Full alias-kinematics consistency; physical mapping validity remains separate.
class AliasKinematicsState {
public:
    /// Transfer the immutable qualified payload.
    AliasKinematicsState(AliasKinematicsState&&) noexcept;
    /// Transfer the immutable payload and deactivate the source.
    AliasKinematicsState& operator=(AliasKinematicsState&&) noexcept;
    AliasKinematicsState(const AliasKinematicsState&) = delete;
    AliasKinematicsState& operator=(const AliasKinematicsState&) = delete;
    /// Release the retained source geometry and packed rows.
    ~AliasKinematicsState();
    /// Whether all supplied alias rows and velocity assertions passed.
    [[nodiscard]] bool active() const noexcept;
    /// Exact derived velocity of a retained point in xyz component 0..2.
    [[nodiscard]] GlobalResult<ExactVelocityView> velocity(PointId, std::size_t) const& noexcept;
    [[nodiscard]] GlobalResult<ExactVelocityView> velocity(PointId, std::size_t) const&& = delete;
    /// Original geometry retained with the complete source layout and rows.
    [[nodiscard]] GlobalResult<std::reference_wrapper<const OwnedGeometry>>
    geometry() const& noexcept;
    [[nodiscard]] GlobalResult<std::reference_wrapper<const OwnedGeometry>>
    geometry() const&& = delete;
    /// Work and packed-payload admission for this call.
    [[nodiscard]] GlobalResult<GlobalWorkPlan> work_plan() const noexcept;
    /// Row equality never supplies a physical mapping certificate.
    [[nodiscard]] static constexpr PhysicalMappingNotChecked mapping_status() noexcept {
        return {};
    }

private:
    friend class GlobalEngine;
    explicit AliasKinematicsState(core::PlannedBuffer<AliasKinematicsData>,
                                  GeometryContext) noexcept;
    core::PlannedBuffer<AliasKinematicsData> data_;
    GeometryContext context_;
};
/// Exclusive prepared exact scratch; mutable arithmetic never lives in a result.
class GlobalWorkspace {
public:
    /// Allocate fixed arithmetic scratch from the cumulative plan before use.
    [[nodiscard]] static core::Result<std::unique_ptr<GlobalWorkspace>>
    create(core::ResourcePlan&, numerics::linear::BackendMemoryLimit,
           core::FailureDetails = {}) noexcept;
    GlobalWorkspace(const GlobalWorkspace&) = delete;
    GlobalWorkspace& operator=(const GlobalWorkspace&) = delete;
    GlobalWorkspace(GlobalWorkspace&&) = delete;
    GlobalWorkspace& operator=(GlobalWorkspace&&) = delete;
    /// Release the prepared scratch when no operation is active.
    ~GlobalWorkspace();
    /// Exact fixed payload charged at creation, excluding allocator metadata.
    [[nodiscard]] static std::uint64_t payload_bytes() noexcept;

private:
    friend class GlobalEngine;
    explicit GlobalWorkspace(core::PlannedBuffer<GlobalScratch>) noexcept;
    core::PlannedBuffer<GlobalScratch> scratch_;
    std::atomic_flag busy_ = ATOMIC_FLAG_INIT;
};
/// Own global topology and the explicit TE transfer, without changing the source snapshot.
[[nodiscard]] GlobalResult<GlobalTopology>
build_global_topology(const OwnedGeometry&, const GlobalTopologyInput&, GlobalLimits,
                      GlobalWorkspace&, core::ResourcePlan&, std::stop_token = {}) noexcept;
/// Check all current finite-facet contacts using a refit of the retained reference tree.
[[nodiscard]] GlobalResult<GlobalContactState>
check_global_geometry(const GlobalTopology&, const OwnedGeometry&, GlobalLimits, GlobalWorkspace&,
                      core::ResourcePlan&, std::stop_token = {}) noexcept;
/// Compare all physical alias columns and derive exact velocities; never qualify a provider.
[[nodiscard]] GlobalResult<AliasKinematicsState>
check_alias_kinematics(const OwnedGeometry&, const AliasKinematicsInput&, GlobalLimits,
                       GlobalWorkspace&, core::ResourcePlan&, std::stop_token = {}) noexcept;
} // namespace uvlas::aero::detail
