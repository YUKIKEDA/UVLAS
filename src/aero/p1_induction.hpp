#pragma once

#include "local_geometry_packet.hpp"

namespace uvlas::aero::detail {
class P1Engine;
struct P1Storage;
struct P1Data;
struct P1Scratch;

/// Desired returned representation and dimensionless full widths of B and dB/dgeometry.
struct P1Accuracy {
    /// Returned centre representation; input words are never projected.
    core::Precision representation;
    /// Maximum full width after multiplying a basis component by the longest edge.
    double tau_basis;
    /// Maximum full width after multiplying a derivative by squared edge length.
    double tau_jacobian;
    /// Require all numerical choices explicitly; no uninitialized request exists.
    P1Accuracy(core::Precision precision, double basis_width, double derivative_width) noexcept
        : representation(precision), tau_basis(basis_width), tau_jacobian(derivative_width) {}

    /// Reject invalid request values before evaluating any basis.
    [[nodiscard]] std::optional<core::FailureCode> validation_error() const noexcept;
};
/// Explicit arithmetic work and external-memory admission; zero means no budget.
struct P1Limits {
    /// Explicit forecast acceptance or a requirement for a proven backend ceiling.
    numerics::linear::BackendMemoryLimit external_memory;
    /// Entire-call scalar-operation and root-trial limit; cancellation is checked inside.
    std::uint64_t max_steps;
    /// Require both backend admission mode and entire-call work budget.
    P1Limits(numerics::linear::BackendMemoryLimit memory, std::uint64_t steps) noexcept
        : external_memory(memory), max_steps(steps) {}
};
/// Formula used for the returned local basis, independent of global geometry status.
enum class P1Branch : std::uint8_t {
    Analytic,              ///< Regularized boundary and constant in-plane sheet formula.
    IdealCoplanarExterior, ///< Zero-radius plane exterior limit retaining normal derivatives.
    FarGauss12,            ///< Fixed tensor twelve-by-twelve area quadrature.
    BoundContour           ///< Three constant directed edge kernels, without a sheet.
};
/// Difference of two fixed quadrature rules; never a certified continuous-integral bound.
struct P1QuadratureDifference {
    /// Per-basis dimensionless norm difference including both rules roundoff radii.
    std::array<double, 3> basis;
    /// Per-basis/per-ambient-column dimensionless norm upper bound, basis-major.
    std::array<double, 39> jacobian_columns;
};
/// The analytic result has no fixed-rule quadrature-difference diagnostic.
struct P1AnalyticEvaluation {};
/// A quadrature diagnostic is explicitly distinct from analytic evaluation.
using P1EvaluationDiagnostic = std::variant<P1AnalyticEvaluation, P1QuadratureDifference>;
/// Checked static expression shape and cumulative arithmetic admission.
struct P1WorkPlan {
    /// Static maximum simultaneously retained expression and accumulator slots.
    std::uint64_t nodes;
    /// Whole-call admission ceiling from the scalar opcode costs.
    std::uint64_t maximum_steps;
    /// Actual sizeof-based workspace and result payload, excluding process overhead.
    std::uint64_t payload_bytes;
    /// Compare every reported admission count exactly.
    bool operator==(const P1WorkPlan&) const = default;
};
/// Exclusive reusable storage. Factories reserve actual node/jet payload before construction.
class P1Workspace {
    class Key {
        friend class P1Workspace;
        Key() = default;
    };

public:
    /// Factory-only key prevents construction with unplanned storage.
    P1Workspace(Key, core::PlannedBuffer<P1Scratch>) noexcept;
    /// Release workspace payload independently of previously returned results.
    ~P1Workspace();
    /// A workspace lease cannot be shared by copying.
    P1Workspace(const P1Workspace&) = delete;
    /// A workspace lease cannot be reassigned by copying.
    P1Workspace& operator=(const P1Workspace&) = delete;
    /// Keep the exclusive workspace address stable.
    P1Workspace(P1Workspace&&) = delete;
    /// Keep the exclusive workspace address stable.
    P1Workspace& operator=(P1Workspace&&) = delete;
    /// Count expressions without allocation, then reserve and construct all arithmetic storage.
    [[nodiscard]] static core::Result<std::unique_ptr<P1Workspace>>
    create(core::ResourcePlan&, core::FailureDetails = {}) noexcept;
    /// Exact cumulatively charged object and arithmetic payload, excluding allocator overhead.
    [[nodiscard]] static std::uint64_t payload_bytes() noexcept;
    /// Largest scalar work count of one local evaluation, including a cold workspace.
    [[nodiscard]] static std::uint64_t maximum_steps() noexcept;
    /// Cold-work ceiling for a bound contour with any active edge subset.
    [[nodiscard]] static std::uint64_t maximum_bound_steps() noexcept;

private:
    friend class P1Engine;
    std::atomic_flag entered_ = ATOMIC_FLAG_INIT;
    core::PlannedBuffer<P1Scratch> scratch_;
};
/// Own all three P1 basis vectors and their ambient thirteen-variable Jacobians.
/// The geometry owner, original query, context, and local-only status survive workspace reuse.
class P1Basis {
public:
    /// Transfer the owner and leave the source inactive.
    P1Basis(P1Basis&&) noexcept;
    /// Release this result and transfer another owner without touching either workspace.
    P1Basis& operator=(P1Basis&&) noexcept;
    /// Results transfer ownership instead of copying their charged payload.
    P1Basis(const P1Basis&) = delete;
    /// Results transfer ownership instead of copying their charged payload.
    P1Basis& operator=(const P1Basis&) = delete;
    /// Release this result and retained immutable geometry.
    ~P1Basis();
    /// Whether this object still owns a complete result.
    [[nodiscard]] bool active() const noexcept;
    /// B[k,xyz] in 1/metre. Invalid indices and moved-from objects fail by value.
    [[nodiscard]] AeroResult<std::reference_wrapper<const ScaledEnclosure>>
    basis(std::size_t k, std::size_t xyz) const& noexcept;
    [[nodiscard]] AeroResult<std::reference_wrapper<const ScaledEnclosure>>
    basis(std::size_t, std::size_t) const&& = delete;
    /// Ambient order: query xyz, vertex0 xyz, vertex1 xyz, vertex2 xyz, radius.
    [[nodiscard]] AeroResult<std::reference_wrapper<const ScaledEnclosure>>
    derivative(std::size_t k, std::size_t xyz, std::size_t variable) const& noexcept;
    [[nodiscard]] AeroResult<std::reference_wrapper<const ScaledEnclosure>>
    derivative(std::size_t, std::size_t, std::size_t) const&& = delete;
    /// Original immutable coordinates and connectivity, retained independently of the packet.
    [[nodiscard]] AeroResult<std::reference_wrapper<const OwnedGeometry>>
    geometry() const& noexcept;
    [[nodiscard]] AeroResult<std::reference_wrapper<const OwnedGeometry>>
    geometry() const&& = delete;
    /// Original unrounded query expression and stable query identity.
    [[nodiscard]] AeroResult<QueryRequest> query() const noexcept;
    /// Original facet, vertex, edge, and scatter descriptors.
    [[nodiscard]] AeroResult<LocalAmbientDescriptor> ambient() const noexcept;
    /// Selected analytic, coplanar exterior, or fixed far branch.
    [[nodiscard]] AeroResult<P1Branch> branch() const noexcept;
    /// Owned branch-specific diagnostic; a difference is never a certified error bound.
    [[nodiscard]] AeroResult<P1EvaluationDiagnostic> diagnostic() const noexcept;
    /// Explicit centre representation and requested normalized full widths.
    [[nodiscard]] AeroResult<P1Accuracy> accuracy() const noexcept;
    /// Scalar operations and root candidates actually executed.
    [[nodiscard]] AeroResult<std::uint64_t> steps() const noexcept;
    /// Static admitted shape, work ceiling, and payload for this evaluator.
    [[nodiscard]] AeroResult<P1WorkPlan> work_plan() const noexcept;
    /// Identity of the actual retained geometry owner.
    [[nodiscard]] AeroResult<std::array<std::byte, 32>> snapshot_identity() const noexcept;
    /// Immutable context captured with this evaluation.
    [[nodiscard]] GeometryContext context() const noexcept {
        return context_;
    }
    /// Local evaluation never qualifies intersections, trailing edges, or the whole model.
    [[nodiscard]] static constexpr GlobalGeometryNotChecked global_status() noexcept {
        return {};
    }
    /// Exact cumulatively charged object and arithmetic payload, excluding allocator overhead.
    [[nodiscard]] static std::uint64_t payload_bytes() noexcept;

private:
    friend class P1Engine;
    [[nodiscard]] const P1Data* data() const noexcept;
    P1Basis(core::PlannedBuffer<P1Storage>, GeometryContext) noexcept;
    core::PlannedBuffer<P1Storage> storage_;
    GeometryContext context_;
};
/// Evaluate all P1 basis columns without projecting original low words or material coordinates.
/// Failures preserve the input packet and every previously returned owning result.
[[nodiscard]] AeroResult<P1Basis> evaluate_p1_basis(const LocalGeometryPacket&, P1Accuracy,
                                                    P1Limits, P1Workspace&, core::ResourcePlan&,
                                                    const std::stop_token& = {}) noexcept;
} // namespace uvlas::aero::detail
