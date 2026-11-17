#pragma once
#include "exact_polynomial_detail.hpp"

#include <array>
#include <optional>
#include <uvlas/structure/pg_graph.hpp>
namespace uvlas::structure::detail {
struct PgGraphData {
    ExactRational x0, width, skew;
    ExactChartPolynomial height;
    core::Index n = 0;
    PgGraphConvention convention = PgGraphConvention::Version1;
    core::PlannedBuffer<PgCellIndex> cells;
    std::optional<PgGraphContentIdentity> identity;
};
struct PgCellPointData {
    std::optional<PgGraph> source;
    ExactRational s, t, u, v;
    PgCellIndex cell;
    PgSide side_s = PgSide::Interior, side_t = PgSide::Interior;
};
struct PgCellPolynomialsData {
    std::array<ExactCellPolynomial, 3> components;
    PgCellIndex cell;
    std::optional<PgGraphContentIdentity> identity;
};
struct PgGraphValueData {
    std::array<ExactRational, 3> components;
    PgCellIndex cell;
    PgDerivativeOrders orders;
    PgDerivativeFrame frame = PgDerivativeFrame::ChartST;
    std::optional<PgGraphContentIdentity> identity;
};
struct PgRectangleActivityData {
    std::optional<PgGraph> source;
    ExactRational area, total;
};
// Planner/source audit: held rational slots and active exact calls are separate.
// prepare: 18 held; three rational clones + polynomial clone(m).
// locate: 27 held; graph clone + at most 9 scalar-producing calls.
// cell: 35 held; 11 scalar calls, 3*Q(2)+Q(3), substitution U=sum(2^(a+b)).
// evaluate: 42 held; at most 8 scalar calls + two derivatives Q(m).
// anchor: 64 held; locate cumulative payload + evaluate cumulative payload.
// map: 22 held; graph clone + 4 scalar-producing calls.
// Q(m)=m*(sizeof(RawTerm)+sizeof(ExactTerm)+sizeof(RationalStorage))
//      +sizeof(PolynomialStorage); zero outputs can allocate less.
struct PgGraphFactory {
    static const PgGraphData* data(const PgGraph&) noexcept;
    static const PgCellPointData* data(const PgCellPoint&) noexcept;
    static const PgCellPolynomialsData* data(const PgCellPolynomials&) noexcept;
    static const PgGraphValueData* data(const PgGraphValue&) noexcept;
    static const PgRectangleActivityData* data(const QualifiedRectangularActivityMap&) noexcept;
    static PgGraph graph(core::PlannedBuffer<PgGraphData>, PgGraphCallInfo) noexcept;
    static PgCellPoint point(core::PlannedBuffer<PgCellPointData>, PgGraphCallInfo) noexcept;
    static PgCellPolynomials polynomials(core::PlannedBuffer<PgCellPolynomialsData>,
                                         PgGraphCallInfo) noexcept;
    static PgGraphValue value(core::PlannedBuffer<PgGraphValueData>, PgGraphCallInfo) noexcept;
    static QualifiedRectangularActivityMap map(core::PlannedBuffer<PgRectangleActivityData>,
                                               PgGraphCallInfo) noexcept;
    static void info(PgGraph&, PgGraphCallInfo) noexcept;
    static void info(PgCellPoint&, PgGraphCallInfo) noexcept;
    static void info(PgCellPolynomials&, PgGraphCallInfo) noexcept;
    static void info(PgGraphValue&, PgGraphCallInfo) noexcept;
    static void info(QualifiedRectangularActivityMap&, PgGraphCallInfo) noexcept;
};
} // namespace uvlas::structure::detail
