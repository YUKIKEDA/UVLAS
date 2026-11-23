#pragma once

#include "global_geometry.hpp"
#include "p1_induction.hpp"

namespace uvlas::aero::detail {
class SharedP1Engine;
struct SharedP1Data;
struct SharedP1AdjointData;
struct P1StrengthTag;
using P1StrengthId = core::Identifier<P1StrengthTag>;

/// One independent circulation column, keyed by its physical face identifier.
struct P0BoundColumn {};
/// Wake triangles hold three columns and wake quads hold their ordered L/R.
using P1FaceColumns =
    std::variant<std::array<P1StrengthId, 3>, std::array<P1StrengthId, 2>, P0BoundColumn>;
struct P1FaceStrengths {
    P1FaceStrengths(PhysicalFaceId physical, const P1FaceColumns& basis) noexcept
        : face(physical), columns(basis) {}
    PhysicalFaceId face;
    P1FaceColumns columns;
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
};

struct P1AssemblyReport {
    P1AssemblyPlan plan;
    std::uint64_t facets = 0, far_facets = 0;
    std::uint64_t cancelled_analytic_pairs = 0, used_steps = 0;
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
