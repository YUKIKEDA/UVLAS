#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stop_token>
#include <utility>
#include <uvlas/structure/exact_polynomial.hpp>

namespace uvlas::structure {
/// Adopted graph and half-cell convention; this tag does not qualify a PG basis.
enum class PgGraphConvention : std::uint8_t {
    Version1 ///< Rational quartic graph, adopted half-cells and first-cell anchor.
};
/// Interior points have one cell. At a boundary, choose its lower or upper cell.
/// At chart zero only Upper is inward; at chart one only Lower is inward.
enum class PgSide : std::uint8_t {
    Interior, ///< A point strictly inside one cell in this axis.
    Lower,    ///< The cell immediately below this boundary.
    Upper     ///< The cell immediately above this boundary.
};
/// Differentiate in global chart coordinates s,t or selected local coordinates u,v.
enum class PgDerivativeFrame : std::uint8_t {
    ChartST, ///< Derivatives in dimensionless chart s,t.
    CellUV   ///< Derivatives in the selected cell's dimensionless u,v.
};
/// Requested cell indices, checked against the receiving graph before use.
struct PgCellIndex {
    core::Index first = 0;  ///< Cell index along chart s.
    core::Index second = 0; ///< Cell index along chart t.
    bool operator==(const PgCellIndex&) const = default;
};
/// Requested mixed derivative; zero in both coordinates returns the position.
struct PgDerivativeOrders {
    core::Index first = 0;  ///< Order along s or local u.
    core::Index second = 0; ///< Order along t or local v.
};
/// Canonical SHA-256 graph-content digest. No factory accepts this as proof input.
/// Physical orientation, origin, reference length and thickness are not covered.
struct PgGraphContentIdentity {
    /// Own an explicit digest; the value alone cannot construct qualified results.
    explicit PgGraphContentIdentity(std::array<std::byte, 32> digest) noexcept : bytes(digest) {}
    std::array<std::byte, 32> bytes; ///< SHA-256 output bytes, independent of host limb width.
    bool operator==(const PgGraphContentIdentity&) const = default;
};
/// Borrowed input for a synchronous factory; every referenced owner must outlive it.
struct PgGraphRequest {
    PgGraphRequest(const ExactRational& origin_x, const ExactRational& width_value,
                   const ExactRational& skew_value, const ExactChartPolynomial& height_value,
                   core::Index grid_size, PgGraphConvention revision) noexcept
        : x0(origin_x),
          width(width_value),
          skew(skew_value),
          height(height_value),
          n(grid_size),
          convention(revision) {}
    const ExactRational& x0;            ///< Dimensionless chart x offset.
    const ExactRational& width;         ///< Strictly positive dimensionless chart width.
    const ExactRational& skew;          ///< Dimensionless affine skew coefficient.
    const ExactChartPolynomial& height; ///< Dimensionless global graph f(s,t), degree at most four.
    core::Index n;                      ///< Adopted internal grid, one, two or four.
    PgGraphConvention convention;       ///< Explicit adopted cell/anchor convention.
};
/// One call's explicit algebra, cell-count and cancellation limits. The external
/// admission applies to all retained call objects plus its largest active scratch.
struct PgGraphLimits {
    PgGraphLimits(ExactLimits algebra_limits, std::uint64_t cell_limit,
                  std::stop_token cancellation_token) noexcept
        : algebra(algebra_limits),
          maximum_cells(cell_limit),
          cancellation(std::move(cancellation_token)) {}
    ExactLimits algebra;          ///< Exact limits and whole-call external admission.
    std::uint64_t maximum_cells;  ///< Maximum full background-cell count; no clipping.
    std::stop_token cancellation; ///< Checked at bounded preparation phases and before publication.
};
/// Owned diagnostics of one successful call. Planned bytes are an admission
/// envelope; actual bytes are the cumulative ResourcePlan delta, not live heap/RSS.
struct PgGraphCallInfo {
    ExactPreparationInfo
        arithmetic; ///< Whole-call maximum bit/raw and aggregate external forecast.
    std::uint64_t planned_owned_payload_upper_bytes =
        0;                                           ///< Cumulative allocation admission envelope.
    std::uint64_t allocated_owned_payload_bytes = 0; ///< Successful ResourcePlan delta.
};
/// Immutable borrowed graph definition. Valid until its graph owner moves or dies.
struct PgGraphView {
    /// Borrow every field explicitly; the source owner supplies their validity.
    PgGraphView(const ExactRational& origin_x, const ExactRational& width_value,
                const ExactRational& skew_value, const ExactChartPolynomial& height_value,
                core::Index grid, PgGraphConvention revision,
                std::span<const PgCellIndex> all_cells,
                PgGraphContentIdentity graph_identity) noexcept
        : x0(origin_x),
          width(width_value),
          skew(skew_value),
          height(height_value),
          n(grid),
          convention(revision),
          cells(all_cells),
          identity(graph_identity) {}
    const ExactRational& x0;            ///< Dimensionless chart x offset.
    const ExactRational& width;         ///< Strictly positive dimensionless chart width.
    const ExactRational& skew;          ///< Dimensionless affine skew coefficient.
    const ExactChartPolynomial& height; ///< Dimensionless global graph f(s,t), degree at most four.
    core::Index n;                      ///< Adopted internal grid, one, two or four.
    PgGraphConvention convention;       ///< Explicit adopted cell/anchor convention.
    std::span<const PgCellIndex> cells; ///< All cells in cs-major, ct-minor order.
    PgGraphContentIdentity
        identity; ///< Producing graph-content identity, not physical qualification.
};
/// Exact chart/local coordinates and an explicit one-sided cell association.
struct PgCellPointView {
    /// Borrow explicit coordinates; a view alone cannot create an owned point.
    PgCellPointView(const ExactRational& chart_s, const ExactRational& chart_t,
                    const ExactRational& local_u, const ExactRational& local_v, PgCellIndex index,
                    PgSide first_side, PgSide second_side,
                    PgGraphContentIdentity graph_identity) noexcept
        : s(chart_s),
          t(chart_t),
          u(local_u),
          v(local_v),
          cell(index),
          side_s(first_side),
          side_t(second_side),
          identity(graph_identity) {}
    const ExactRational& s; ///< Exact chart s in [0,1].
    const ExactRational& t; ///< Exact chart t in [0,1].
    const ExactRational& u; ///< Exact local u in [-1,1].
    const ExactRational& v; ///< Exact local v in [-1,1].
    PgCellIndex cell;       ///< Selected one-sided cell.
    PgSide side_s;          ///< Explicit s-boundary choice.
    PgSide side_t;          ///< Explicit t-boundary choice.
    PgGraphContentIdentity
        identity; ///< Producing graph-content identity, not physical qualification.
};
/// Cell polynomials in x,y,z order, with their producing graph/cell association.
struct PgCellPolynomialsView {
    /// Borrow explicit components and association without conferring qualification.
    PgCellPolynomialsView(std::span<const ExactCellPolynomial> values, PgCellIndex index,
                          PgGraphContentIdentity graph_identity) noexcept
        : components(values), cell(index), identity(graph_identity) {}
    std::span<const ExactCellPolynomial> components; ///< Exactly three polynomials in x,y,z order.
    PgCellIndex cell;                                ///< Selected one-sided cell.
    PgGraphContentIdentity
        identity; ///< Producing graph-content identity, not physical qualification.
};
/// Exact vector in x,y,z order; these are dimensionless graph-chart quantities.
struct PgGraphValueView {
    /// Borrow explicit values and their producing derivative association.
    PgGraphValueView(std::span<const ExactRational> values, PgCellIndex index,
                     PgDerivativeOrders derivative_orders, PgDerivativeFrame derivative_frame,
                     PgGraphContentIdentity graph_identity) noexcept
        : components(values),
          cell(index),
          orders(derivative_orders),
          frame(derivative_frame),
          identity(graph_identity) {}
    std::span<const ExactRational> components; ///< Exactly three values in x,y,z order.
    PgCellIndex cell;                          ///< Selected one-sided cell.
    PgDerivativeOrders orders;                 ///< Mixed derivative used for all three components.
    PgDerivativeFrame frame;                   ///< Coordinate frame of the derivative values.
    PgGraphContentIdentity
        identity; ///< Producing graph-content identity, not physical qualification.
};
namespace detail {
struct PgGraphData;
struct PgCellPointData;
struct PgCellPolynomialsData;
struct PgGraphValueData;
struct PgRectangleActivityData;
struct PgGraphFactory;
} // namespace detail
class PgCellPoint;
class PgCellPolynomials;
class PgGraphValue;
class QualifiedRectangularActivityMap;

/// Immutable owned rational graph and all half-cell definitions. No physical
/// reference, rank, conditioning, quadrature or material qualification is implied.
class PgGraph {
public:
    /// Transfer ownership. Operations on the moved-from object return InvalidInput.
    PgGraph(PgGraph&&) noexcept;
    PgGraph& operator=(PgGraph&&) noexcept;
    PgGraph(const PgGraph&) = delete;
    PgGraph& operator=(const PgGraph&) = delete;
    ~PgGraph();
    /// Validate width>0, degree<=4, n=1/2/4 and convention, then own exact clones.
    /// Invalid definitions return InvalidInput; budget/overflow/bad_alloc return
    /// ResourceExceeded. Cancellation is Cancelled; unproven hard bounds are Unsupported.
    [[nodiscard]] static core::Result<PgGraph> prepare(const PgGraphRequest&, const PgGraphLimits&,
                                                       core::ResourcePlan&) noexcept;
    /// Make an independent budgeted copy; no mutable cache or workspace is shared.
    [[nodiscard]] core::Result<PgGraph> clone(const PgGraphLimits&,
                                              core::ResourcePlan&) const noexcept;
    /// Borrow definition and all cells, or InvalidInput for a moved-from owner.
    [[nodiscard]] core::Result<PgGraphView> view() const& noexcept;
    [[nodiscard]] core::Result<PgGraphView> view() const&& = delete;
    /// Produce an exact location. Interior forbids boundaries; Lower/Upper are
    /// required on boundaries and must select an inward, existing half-cell.
    [[nodiscard]] core::Result<PgCellPoint> locate(const ExactRational& s, const ExactRational& t,
                                                   PgSide side_s, PgSide side_t,
                                                   const PgGraphLimits&,
                                                   core::ResourcePlan&) const noexcept;
    /// Invert a specified cell's local coordinates. Endpoint sides must agree
    /// with that cell; out-of-range coordinates and outward choices are InvalidInput.
    [[nodiscard]] core::Result<PgCellPoint> from_local(PgCellIndex, const ExactRational& u,
                                                       const ExactRational& v, PgSide side_s,
                                                       PgSide side_t, const PgGraphLimits&,
                                                       core::ResourcePlan&) const noexcept;
    /// Return owned local x,y,z polynomials. No output borrows the graph storage.
    /// Invalid cell indices are InvalidInput; limits and allocations can fail.
    [[nodiscard]] core::Result<PgCellPolynomials>
    cell_polynomials(PgCellIndex, const PgGraphLimits&, core::ResourcePlan&) const noexcept;
    /// Evaluate a mixed derivative at an exact one-sided location. Association
    /// uses actual owned graph contents, not only the supplied point's digest.
    /// Incompatible content returns StaleGeneration. Negative orders/unknown tags
    /// are InvalidInput; orders above the degree return exact zero without a scale power.
    [[nodiscard]] core::Result<PgGraphValue> evaluate(const PgCellPoint&, PgDerivativeOrders,
                                                      PgDerivativeFrame, const PgGraphLimits&,
                                                      core::ResourcePlan&) const noexcept;
    /// Evaluate the chart anchor s=0 from first s-cell/u=-1. This is not a
    /// physical-space transform; the requested derivative is in chart t only.
    [[nodiscard]] core::Result<PgGraphValue> anchor_at(const ExactRational& t, PgSide side_t,
                                                       core::Index t_order, const PgGraphLimits&,
                                                       core::ResourcePlan&) const noexcept;
    /// Inspect owned diagnostics even after a move; no content proof is inferred.
    [[nodiscard]] const PgGraphCallInfo& preparation() const& noexcept;
    [[nodiscard]] const PgGraphCallInfo& preparation() const&& = delete;

private:
    friend struct detail::PgGraphFactory;
    explicit PgGraph(core::PlannedBuffer<detail::PgGraphData>, PgGraphCallInfo) noexcept;
    core::PlannedBuffer<detail::PgGraphData> data_;
    PgGraphCallInfo info_;
};

/// Owned immutable result; its view borrows only this object's storage.
class PgCellPoint {
public:
    /// Transfer ownership; view() on the moved-from object returns InvalidInput.
    PgCellPoint(PgCellPoint&&) noexcept;
    PgCellPoint& operator=(PgCellPoint&&) noexcept;
    PgCellPoint(const PgCellPoint&) = delete;
    PgCellPoint& operator=(const PgCellPoint&) = delete;
    ~PgCellPoint();
    /// Borrow validated values until this owner moves, is assigned or is destroyed.
    [[nodiscard]] core::Result<PgCellPointView> view() const& noexcept;
    [[nodiscard]] core::Result<PgCellPointView> view() const&& = delete;
    /// Borrow producing-call diagnostics, independent of the ResourcePlan lifetime.
    [[nodiscard]] const PgGraphCallInfo& preparation() const& noexcept;
    [[nodiscard]] const PgGraphCallInfo& preparation() const&& = delete;

private:
    friend struct detail::PgGraphFactory;
    explicit PgCellPoint(core::PlannedBuffer<detail::PgCellPointData>, PgGraphCallInfo) noexcept;
    core::PlannedBuffer<detail::PgCellPointData> data_;
    PgGraphCallInfo info_;
};

/// Owned immutable result; its view borrows only this object's storage.
class PgCellPolynomials {
public:
    /// Transfer ownership; view() on the moved-from object returns InvalidInput.
    PgCellPolynomials(PgCellPolynomials&&) noexcept;
    PgCellPolynomials& operator=(PgCellPolynomials&&) noexcept;
    PgCellPolynomials(const PgCellPolynomials&) = delete;
    PgCellPolynomials& operator=(const PgCellPolynomials&) = delete;
    ~PgCellPolynomials();
    /// Borrow validated values until this owner moves, is assigned or is destroyed.
    [[nodiscard]] core::Result<PgCellPolynomialsView> view() const& noexcept;
    [[nodiscard]] core::Result<PgCellPolynomialsView> view() const&& = delete;
    /// Borrow producing-call diagnostics, independent of the ResourcePlan lifetime.
    [[nodiscard]] const PgGraphCallInfo& preparation() const& noexcept;
    [[nodiscard]] const PgGraphCallInfo& preparation() const&& = delete;

private:
    friend struct detail::PgGraphFactory;
    explicit PgCellPolynomials(core::PlannedBuffer<detail::PgCellPolynomialsData>,
                               PgGraphCallInfo) noexcept;
    core::PlannedBuffer<detail::PgCellPolynomialsData> data_;
    PgGraphCallInfo info_;
};

/// Owned immutable result; its view borrows only this object's storage.
class PgGraphValue {
public:
    /// Transfer ownership; view() on the moved-from object returns InvalidInput.
    PgGraphValue(PgGraphValue&&) noexcept;
    PgGraphValue& operator=(PgGraphValue&&) noexcept;
    PgGraphValue(const PgGraphValue&) = delete;
    PgGraphValue& operator=(const PgGraphValue&) = delete;
    ~PgGraphValue();
    /// Borrow validated values until this owner moves, is assigned or is destroyed.
    [[nodiscard]] core::Result<PgGraphValueView> view() const& noexcept;
    [[nodiscard]] core::Result<PgGraphValueView> view() const&& = delete;
    /// Borrow producing-call diagnostics, independent of the ResourcePlan lifetime.
    [[nodiscard]] const PgGraphCallInfo& preparation() const& noexcept;
    [[nodiscard]] const PgGraphCallInfo& preparation() const&& = delete;

private:
    friend struct detail::PgGraphFactory;
    explicit PgGraphValue(core::PlannedBuffer<detail::PgGraphValueData>, PgGraphCallInfo) noexcept;
    core::PlannedBuffer<detail::PgGraphValueData> data_;
    PgGraphCallInfo info_;
};

/// Owned full-unit-square activity map. It cannot be made from flags or indices,
/// and its coverage proof is not a rank, quadrature or physical-reference proof.
class QualifiedRectangularActivityMap {
public:
    /// Transfer ownership; moved-from maps cannot validate a graph.
    QualifiedRectangularActivityMap(QualifiedRectangularActivityMap&&) noexcept;
    QualifiedRectangularActivityMap& operator=(QualifiedRectangularActivityMap&&) noexcept;
    QualifiedRectangularActivityMap(const QualifiedRectangularActivityMap&) = delete;
    QualifiedRectangularActivityMap& operator=(const QualifiedRectangularActivityMap&) = delete;
    ~QualifiedRectangularActivityMap();
    /// Generate every background cell and exact full coverage from owned contents.
    /// Moved-from inputs are InvalidInput; limits/cancellation have the same failures
    /// as graph preparation. Polygon input is intentionally absent from this API.
    [[nodiscard]] static core::Result<QualifiedRectangularActivityMap>
    prepare(const PgGraph&, const PgGraphLimits&, core::ResourcePlan&) noexcept;
    /// Compare actual immutable graph contents as well as the canonical identity.
    /// Different contents return StaleGeneration; moved-from owners are InvalidInput.
    [[nodiscard]] core::Result<void> validate_for(const PgGraph&) const noexcept;
    /// Borrow the map's independent graph snapshot and all standard-order cells.
    [[nodiscard]] core::Result<PgGraphView> graph() const& noexcept;
    [[nodiscard]] core::Result<PgGraphView> graph() const&& = delete;
    /// Return an owned exact copy of the chart ds*dt area of each listed cell.
    /// This is not a physical surface measure.
    [[nodiscard]] core::Result<ExactRational> cell_area(const PgGraphLimits&,
                                                        core::ResourcePlan&) const noexcept;
    /// Return an owned exact copy of the full chart ds*dt area, equal to one.
    [[nodiscard]] core::Result<ExactRational> total_area(const PgGraphLimits&,
                                                         core::ResourcePlan&) const noexcept;
    /// Borrow producing-call diagnostics.
    [[nodiscard]] const PgGraphCallInfo& preparation() const& noexcept;
    [[nodiscard]] const PgGraphCallInfo& preparation() const&& = delete;

private:
    friend struct detail::PgGraphFactory;
    explicit QualifiedRectangularActivityMap(core::PlannedBuffer<detail::PgRectangleActivityData>,
                                             PgGraphCallInfo) noexcept;
    core::PlannedBuffer<detail::PgRectangleActivityData> data_;
    PgGraphCallInfo info_;
};
} // namespace uvlas::structure
