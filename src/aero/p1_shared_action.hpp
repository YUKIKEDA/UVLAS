#pragma once

#include "global_geometry.hpp"
#include "p1_induction.hpp"

namespace uvlas::aero::detail {
class SharedP1Engine;
struct SharedP1Data;
struct SharedP1AdjointData;
struct BoundWakeFacetStorageData;
struct P1StrengthTag;
/// Independent wake identifiers and bound circulation identities have disjoint domains.
class P1StrengthId {
    using WakeId = core::Identifier<P1StrengthTag>;

public:
    P1StrengthId() noexcept = default;
    explicit P1StrengthId(std::uint64_t wake) noexcept : identity_(WakeId{wake}) {}

    /// Explicitly tie a trace to a bound circulation column.
    [[nodiscard]] static P1StrengthId bound(PhysicalFaceId face) noexcept {
        return P1StrengthId{face};
    }

    [[nodiscard]] std::optional<PhysicalFaceId> bound_face() const noexcept {
        if (const auto* face = std::get_if<PhysicalFaceId>(&identity_))
            return *face;
        return std::nullopt;
    }

    [[nodiscard]] std::optional<std::uint64_t> wake_value() const noexcept {
        if (const auto* wake = std::get_if<WakeId>(&identity_))
            return wake->value();
        return std::nullopt;
    }

    /// Wake IDs sort first. Each domain uses its own stable numeric order.
    auto operator<=>(const P1StrengthId&) const = default;

private:
    explicit P1StrengthId(PhysicalFaceId face) noexcept : identity_(face) {}
    std::variant<WakeId, PhysicalFaceId> identity_{WakeId{0}};
};

/// One independent circulation column, keyed by its physical face identifier.
struct P0BoundColumn {};
/// Wake triangles hold three columns and wake quads hold their ordered L/R.
using P1FaceColumns =
    std::variant<std::array<P1StrengthId, 3>, std::array<P1StrengthId, 2>, P0BoundColumn>;
struct P1FaceStrengths {
    P1FaceStrengths() noexcept = default; ///< Empty staging slot filled before qualification.
    P1FaceStrengths(PhysicalFaceId physical, const P1FaceColumns& basis) noexcept
        : face(physical), columns(basis) {}
    PhysicalFaceId face{0};
    P1FaceColumns columns{P0BoundColumn{}};
};

struct P1ActionStamp {
    std::uint64_t issuer = 0, generation = 0;
    bool operator==(const P1ActionStamp&) const = default;
};

/// Current contact qualification and original query are both checked before reuse.
struct P1ActionSource {
    const GlobalContactState& contact;
    const ResolvedQuery& query;
    P1ActionStamp stamp;
};

struct P1AssemblyLimits {
    GeometryLimits geometry;
    P1Limits kernel;
    std::uint64_t steps = 0;
};

struct P1AssemblyPlan {
    std::uint64_t payload_bytes = 0;
    std::uint64_t estimated_metadata_bytes = 0;
    std::uint64_t reserved_steps = 0;
    std::uint64_t preflight_steps = 0; ///< Measured visits spent producing the forecast.
};

/// The first successful action populates storage for later actions of the same shape.
struct P1SerialAssemblyPlan {
    P1AssemblyPlan ordinary;
    std::uint64_t initial_payload_bytes{};
    std::uint64_t repeated_payload_bytes{};
    std::uint64_t repeated_reserved_steps{};
    bool reusable{};
};

struct P1AssemblyReport {
    P1AssemblyPlan plan;
    std::uint64_t facets = 0, far_facets = 0;
    std::uint64_t cancelled_analytic_pairs = 0, used_steps = 0;
    std::uint64_t excluded_bound_occurrences = 0; ///< Structural midpoint self-edge only.
};

class SharedP1Workspace {
public:
    SharedP1Workspace() noexcept = default;
    SharedP1Workspace(const SharedP1Workspace&) = delete;
    SharedP1Workspace& operator=(const SharedP1Workspace&) = delete;

private:
    friend class SharedP1Engine;
    std::atomic_flag entered_ = ATOMIC_FLAG_INIT;
};

/// Retain the immutable sparse induction action and all local qualifications.
class SharedP1Action {
public:
    SharedP1Action(SharedP1Action&&) noexcept = default;
    SharedP1Action& operator=(SharedP1Action&&) noexcept = default;
    SharedP1Action(const SharedP1Action&) = delete;
    SharedP1Action& operator=(const SharedP1Action&) = delete;

    [[nodiscard]] bool active() const noexcept;
    [[nodiscard]] std::size_t point_count() const noexcept;
    [[nodiscard]] std::size_t strength_count() const noexcept;
    [[nodiscard]] AeroResult<PointId> point(std::size_t) const noexcept;
    [[nodiscard]] AeroResult<P1StrengthId> strength(std::size_t) const noexcept;
    [[nodiscard]] AeroResult<P1AssemblyReport> report() const noexcept;
    [[nodiscard]] AeroResult<P1EvaluationDiagnostic> facet_diagnostic(std::size_t) const noexcept;

private:
    friend class SharedP1Engine;
    SharedP1Action(std::shared_ptr<const SharedP1Data> data, GeometryContext context) noexcept
        : data_(std::move(data)), context_(context) {}
    std::shared_ptr<const SharedP1Data> data_;
    GeometryContext context_;
};

/// Real input uses a zero imaginary word. Both words of each component are consumed.
struct P1ComplexWords {
    ScalarWords real{SingleWord{0}};
    ScalarWords imaginary{SingleWord{0}};
};
struct P1ComplexEnclosure {
    ScaledEnclosure real, imaginary;
};

/// Point directions follow the action's canonical point order; strength directions follow
/// its strength order. A free query has a direction; a material query uses point directions.
struct P1ActionDirection {
    std::span<const std::array<ScalarWords, 3>> points;
    std::span<const P1ComplexWords> strengths;
    std::optional<std::array<ScalarWords, 3>> free_query;
    ScalarWords radius{SingleWord{0}};
};

/// Owning value and optional directional derivative. Enclosures concern the retained
/// evaluation formula; a far quadrature diagnostic is still available on the action.
class P1ActionValue {
public:
    [[nodiscard]] const std::array<P1ComplexEnclosure, 3>& velocity() const& noexcept;
    [[nodiscard]] const std::array<P1ComplexEnclosure, 3>& velocity() const&& = delete;
    [[nodiscard]] const std::optional<std::array<P1ComplexEnclosure, 3>>& jvp() const& noexcept;
    [[nodiscard]] const std::optional<std::array<P1ComplexEnclosure, 3>>& jvp() const&& = delete;
    [[nodiscard]] std::uint64_t used_steps() const noexcept;

private:
    friend class SharedP1Engine;
    P1ActionValue(std::shared_ptr<const SharedP1Data>, const std::array<P1ComplexEnclosure, 3>&,
                  const std::optional<std::array<P1ComplexEnclosure, 3>>&, std::uint64_t) noexcept;
    std::shared_ptr<const SharedP1Data> parent_;
    std::array<P1ComplexEnclosure, 3> velocity_;
    std::optional<std::array<P1ComplexEnclosure, 3>> jvp_;
    std::uint64_t steps_;
};

/// Pullback for the real inner product Re(conj(left) right) on complex strengths
/// and velocity, and the ordinary real inner product on geometry and radius.
class P1ActionAdjoint {
public:
    P1ActionAdjoint(P1ActionAdjoint&&) noexcept;
    P1ActionAdjoint& operator=(P1ActionAdjoint&&) noexcept;
    ~P1ActionAdjoint();
    P1ActionAdjoint(const P1ActionAdjoint&) = delete;
    P1ActionAdjoint& operator=(const P1ActionAdjoint&) = delete;

    [[nodiscard]] AeroResult<ScaledEnclosure> point(std::size_t, std::size_t axis) const noexcept;
    [[nodiscard]] AeroResult<P1ComplexEnclosure> strength(std::size_t) const noexcept;
    [[nodiscard]] AeroResult<std::optional<ScaledEnclosure>>
    free_query(std::size_t axis) const noexcept;
    [[nodiscard]] AeroResult<ScaledEnclosure> radius() const noexcept;
    [[nodiscard]] std::uint64_t used_steps() const noexcept;

private:
    friend class SharedP1Engine;
    P1ActionAdjoint(core::PlannedBuffer<SharedP1AdjointData>, GeometryContext) noexcept;
    core::PlannedBuffer<SharedP1AdjointData> data_;
    GeometryContext context_;
};

/// Count the complete assembly under the retained geometry shape before allocation.
[[nodiscard]] AeroResult<P1AssemblyPlan> forecast_shared_p1(const P1ActionSource&,
                                                            std::span<const P1FaceStrengths>,
                                                            std::span<const P1StrengthId>,
                                                            const P1AssemblyLimits&,
                                                            const std::stop_token& = {}) noexcept;

[[nodiscard]] AeroResult<SharedP1Action>
prepare_shared_p1(const P1ActionSource&, std::span<const P1FaceStrengths>,
                  std::span<const P1StrengthId>, LocalAccuracy, P1Accuracy, const P1AssemblyLimits&,
                  SharedP1Workspace&, LocalPacketWorkspace&, P1Workspace&, core::ResourcePlan&,
                  const std::stop_token& = {}) noexcept;

/// Forecast or prepare the UW-D01 mixed action. Bound quadrilaterals require
/// P0BoundColumn, while wake faces retain their existing P1 scatter.
[[nodiscard]] AeroResult<P1AssemblyPlan>
forecast_bound_wake_induction(const P1ActionSource&, std::span<const P1FaceStrengths>,
                              std::span<const P1StrengthId>, const P1AssemblyLimits&,
                              const std::stop_token& = {}) noexcept;

[[nodiscard]] AeroResult<SharedP1Action> prepare_bound_wake_induction(
    const P1ActionSource&, std::span<const P1FaceStrengths>, std::span<const P1StrengthId>,
    LocalAccuracy, P1Accuracy, const P1AssemblyLimits&, SharedP1Workspace&, LocalPacketWorkspace&,
    P1Workspace&, core::ResourcePlan&, const std::stop_token& = {}) noexcept;

/// Optional facet, local packet and basis storage for a serial sequence of mixed actions.
/// The accounting plan must outlive this storage. A previous action and its owning
/// results must all be destroyed before preparing another action in this storage.
class BoundWakeFacetStorage {
public:
    explicit BoundWakeFacetStorage(core::ResourcePlan& memory) noexcept : memory_(memory) {}
    BoundWakeFacetStorage(const BoundWakeFacetStorage&) = delete;
    BoundWakeFacetStorage& operator=(const BoundWakeFacetStorage&) = delete;

private:
    friend class SharedP1Engine;
    core::ResourcePlan& memory_;
    std::atomic_flag preparing_ = ATOMIC_FLAG_INIT;
    std::shared_ptr<BoundWakeFacetStorageData> data_;
};

/// Forecast a serial sequence with the same source extents and scoped results.
[[nodiscard]] AeroResult<P1SerialAssemblyPlan>
forecast_bound_wake_serial_induction(const P1ActionSource&, std::span<const P1FaceStrengths>,
                                     std::span<const P1StrengthId>, const P1AssemblyLimits&,
                                     const std::stop_token& = {}) noexcept;

/// Use this storage's resource plan for a new qualified query of the same facet count.
[[nodiscard]] AeroResult<SharedP1Action> prepare_bound_wake_induction_reusing_facets(
    const P1ActionSource&, std::span<const P1FaceStrengths>, std::span<const P1StrengthId>,
    LocalAccuracy, P1Accuracy, const P1AssemblyLimits&, SharedP1Workspace&, LocalPacketWorkspace&,
    P1Workspace&, BoundWakeFacetStorage&, const std::stop_token& = {}) noexcept;

/// Matrix assembly may apply values and transpose strengths. Coordinate derivatives are excluded.
[[nodiscard]] AeroResult<SharedP1Action> prepare_bound_wake_values(
    const P1ActionSource&, std::span<const P1FaceStrengths>, std::span<const P1StrengthId>,
    LocalAccuracy, P1Accuracy, const P1AssemblyLimits&, SharedP1Workspace&, LocalPacketWorkspace&,
    P1Workspace&, core::ResourcePlan&, BoundWakeFacetStorage*, const std::stop_token&) noexcept;

/// Prepare mixed induction at exactly the midpoint of a non-TE canonical bound edge.
/// Exclude all occurrences of that bound contour only. Retain the full wake sheet.
[[nodiscard]] AeroResult<SharedP1Action> prepare_bound_wake_midpoint_induction(
    const P1ActionSource&, EdgeId self_edge, std::span<const P1FaceStrengths>,
    std::span<const P1StrengthId>, LocalAccuracy, P1Accuracy, const P1AssemblyLimits&,
    SharedP1Workspace&, LocalPacketWorkspace&, P1Workspace&, core::ResourcePlan&,
    const std::stop_token& = {}) noexcept;

/// Include canonical-edge lookup and contour exclusion in the serial forecast.
[[nodiscard]] AeroResult<P1SerialAssemblyPlan> forecast_bound_wake_midpoint_serial_induction(
    const P1ActionSource&, std::span<const P1FaceStrengths>, std::span<const P1StrengthId>,
    const P1AssemblyLimits&, const std::stop_token& = {}) noexcept;

[[nodiscard]] AeroResult<SharedP1Action> prepare_bound_wake_midpoint_induction(
    const P1ActionSource&, EdgeId self_edge, std::span<const P1FaceStrengths>,
    std::span<const P1StrengthId>, LocalAccuracy, P1Accuracy, const P1AssemblyLimits&,
    SharedP1Workspace&, LocalPacketWorkspace&, P1Workspace&, BoundWakeFacetStorage&,
    const std::stop_token& = {}) noexcept;

/// The input coefficient order matches the prepared IDs exactly. Geometry is unchanged.
[[nodiscard]] AeroResult<P1ActionValue>
apply_shared_p1(const SharedP1Action&, const P1ActionSource&, std::span<const P1ComplexWords>,
                const std::optional<P1ActionDirection>&, std::uint64_t max_steps,
                SharedP1Workspace&, const std::stop_token& = {}) noexcept;

[[nodiscard]] AeroResult<P1ActionAdjoint>
adjoint_shared_p1(const SharedP1Action&, const P1ActionSource&, std::span<const P1ComplexWords>,
                  const std::array<P1ComplexWords, 3>& dual, std::uint64_t max_steps,
                  SharedP1Workspace&, core::ResourcePlan&, const std::stop_token& = {}) noexcept;

} // namespace uvlas::aero::detail
