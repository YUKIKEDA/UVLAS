#include "case_wing_sharing_private.hpp"

#include <algorithm>

namespace uvlas::coupling::detail::wing_geometry {
namespace {
namespace pg = structure::detail;
namespace p = aero::detail::pchip;
using namespace wing_chart;
using enum core::FailureCode;
using pg::partition_add;
using pg::partition_fail;
using pg::partition_index;
using pg::partition_take;
using XYZ = std::array<Rational, 3>;

template <class T>
core::PlannedBuffer<T> allocate(Math& m, std::uint64_t count) {
    return partition_take(
        m.plan.allocate<T>(partition_index(count), core::Operation::PrepareModel));
}

bool same_xyz(const XYZ& a, const XYZ& b, Math& m) {
    for (std::size_t axis = 0; axis < a.size(); ++axis)
        if (m.compare(a[axis], b[axis]) != 0)
            return false;
    return true;
}

struct Boundary {
    const Chart& chart;
    std::span<const Point> points;

    Point project(const XYZ& world, Math& m) const {
        auto candidate = chart.graph.chart_projection(world, m);
        if (!same_xyz(world, chart.graph.evaluate_chart(candidate, m, false).position, m))
            partition_fail(GeometryInvalid);
        return candidate;
    }

    XYZ world(const Point& q, Math& m) const {
        return chart.graph.evaluate_chart(q, m, false).position;
    }
};

Rational parameter(const Boundary& boundary, const XYZ& world, Math& m) {
    auto p = boundary_parameter(boundary.points, boundary.project(world, m), m);
    if (!p)
        partition_fail(GeometryInvalid);
    return std::move(*p);
}

bool equal_point(const Point& a, const Point& b, Math& m) {
    return m.compare(a.s, b.s) == 0 && m.compare(a.t, b.t) == 0;
}

void prove_shared_arc(const Boundary& own, const Boundary& other, bool reverse, Math& m,
                      CaseWingDiagnostic& diagnostic) {
    if (own.points.size() < 2 || other.points.size() < 2)
        partition_fail(GeometryInvalid);
    const auto count = partition_add(own.points.size(), other.points.size());
    auto vertices = allocate<Rational>(m, count);
    auto values = vertices.values();
    std::size_t size = 0;
    for (std::size_t i = 0; i < own.points.size(); ++i) {
        m.visit(1);
        values[size++] = integer(i, m);
    }
    const auto other_vertex = [&](std::size_t i) -> const Point& {
        return other.points[reverse ? other.points.size() - 1 - i : i];
    };
    auto last = parameter(own, other.world(other_vertex(0), m), m);
    if (last.sign() != 0)
        partition_fail(GeometryInvalid);
    for (std::size_t i = 1; i < other.points.size(); ++i) {
        m.visit(1);
        auto next = parameter(own, other.world(other_vertex(i), m), m);
        if (m.compare(last, next) >= 0)
            partition_fail(GeometryInvalid);
        values[size++] = m.copy(next);
        last = std::move(next);
    }
    if (m.compare(last, integer(own.points.size() - 1, m)) != 0)
        partition_fail(GeometryInvalid);
    for (std::size_t i = 0; i < own.points.size(); ++i)
        (void)parameter(other, own.world(own.points[i], m), m);
    values = values.first(size);
    std::ranges::sort(values, [&](const auto& a, const auto& b) { return m.compare(a, b) < 0; });
    const auto duplicates = std::ranges::unique(
        values, [&](const auto& a, const auto& b) { return m.compare(a, b) == 0; });
    const auto unique = values.first(static_cast<std::size_t>(duplicates.begin() - values.begin()));
    constexpr std::size_t degree = 16;
    const auto denominator = integer(degree, m);
    std::size_t segment = 0;
    for (std::size_t i = 1; i < unique.size(); ++i) {
        m.visit(1);
        while (m.compare(unique[i - 1], integer(segment + 1, m)) >= 0)
            ++segment;
        diagnostic.strip = segment;
        const auto start = integer(segment, m);
        const auto from = m.sub(unique[i - 1], start);
        const auto to = m.sub(unique[i], start);
        const auto& a = own.points[segment];
        const auto& b = own.points[segment + 1];
        const auto c = other.project(own.world(blend(a, b, from, m), m), m);
        const auto d = other.project(own.world(blend(a, b, to, m), m), m);
        if (equal_point(c, d, m))
            partition_fail(GeometryInvalid);
        const auto segment_present = [&] {
            for (std::size_t j = 1; j < other.points.size(); ++j) {
                m.visit(1);
                if (on_segment(other.points[j - 1], other.points[j], c, m) &&
                    on_segment(other.points[j - 1], other.points[j], d, m)) {
                    diagnostic.other_strip = j - 1;
                    return true;
                }
            }
            return false;
        };
        if (!segment_present())
            partition_fail(GeometryInvalid);

        // Evaluation at 17 distinct rational arguments is an invertible Vandermonde basis
        // for degree <= 16, including the composed height. Every residual coefficient
        // vanishes exactly iff these values vanish. Endpoint membership and injectivity
        // then give the orientation and the portion of the partner segment covered.
        for (std::size_t k = 1; k < degree; ++k) {
            m.visit(1);
            const auto u = m.div(integer(k, m), denominator);
            const auto q = blend(a, b, blend(from, to, u, m), m);
            if (!on_segment(c, d, other.project(own.world(q, m), m), m))
                partition_fail(GeometryInvalid);
        }
    }
}

bool same_words(p::RationalView lhs, p::RationalView rhs, Math& m) {
    m.charge(partition_add(partition_add(lhs.numerator.size(), rhs.numerator.size()),
                           partition_add(lhs.denominator.size(), rhs.denominator.size())));
    return lhs.sign == rhs.sign && std::ranges::equal(lhs.numerator, rhs.numerator) &&
           std::ranges::equal(lhs.denominator, rhs.denominator);
}

bool same_series(std::span<const p::RationalView> lhs, std::span<const p::RationalView> rhs,
                 Math& m) {
    return std::ranges::equal(lhs, rhs, [&](auto a, auto b) { return same_words(a, b, m); });
}

bool zero(const CaseMotionSeriesSnapshotView& series, Math& m) {
    const auto values = std::visit([](const auto& value) { return value.coefficients; }, series);
    return std::ranges::all_of(values, [&](const auto& value) {
        m.visit(1);
        return value.sign == 0;
    });
}

bool same_trajectory(const CaseMotionSeriesSnapshotView& lhs,
                     const CaseMotionSeriesSnapshotView& rhs, Math& m) {
    if (zero(lhs, m) && zero(rhs, m))
        return true;
    if (const auto* a = std::get_if<CaseMotionFormulaSnapshotView>(&lhs)) {
        const auto* b = std::get_if<CaseMotionFormulaSnapshotView>(&rhs);
        if (b == nullptr || a->domain.has_value() != b->domain.has_value())
            return false;
        // The qualified formula domain includes origin and scale. These two
        // fields specify the coefficient parameterization and its derivatives.
        return (!a->domain || (same_words((*a->domain)[0], (*b->domain)[0], m) &&
                               same_words((*a->domain)[1], (*b->domain)[1], m))) &&
               same_series(a->coefficients, b->coefficients, m);
    }
    const auto& a = std::get<p::PchipSnapshotView>(lhs);
    const auto* b = std::get_if<p::PchipSnapshotView>(&rhs);
    if (b == nullptr || a.formula_version != b->formula_version ||
        !same_series(a.input.samples.times, b->input.samples.times, m) ||
        !same_series(a.coefficients, b->coefficients, m) ||
        a.input.boundary.policy.index() != b->input.boundary.policy.index())
        return false;
    const auto* periodic = std::get_if<p::PeriodicBoundary>(&a.input.boundary.policy);
    return periodic == nullptr ||
           same_words(periodic->period,
                      std::get<p::PeriodicBoundary>(b->input.boundary.policy).period, m);
}

void qualify_motion(const Patch& a, const Patch& b, Math& m) {
    const auto snapshot = [](const auto& motion) -> std::optional<CaseMotionSnapshotView> {
        if (!motion)
            return {};
        return partition_take(motion->view());
    };
    const auto lhs = snapshot(a.motion), rhs = snapshot(b.motion);
    const auto stationary = [&](const auto& motion) {
        return !motion || (zero(motion->time.translation, m) && zero(motion->time.rotation, m));
    };
    if (stationary(lhs) && stationary(rhs))
        return;
    if (!lhs || !rhs || lhs->formula_version != rhs->formula_version ||
        !same_series(lhs->frame.components, rhs->frame.components, m) ||
        (lhs->time.pi_power != rhs->time.pi_power &&
         (!zero(lhs->time.rotation, m) || !zero(rhs->time.rotation, m))) ||
        !same_trajectory(lhs->time.translation, rhs->time.translation, m) ||
        !same_trajectory(lhs->time.rotation, rhs->time.rotation, m))
        partition_fail(GeometryInvalid);
}

struct OrderedPoint {
    Rational at;
    std::size_t point{};
};

std::span<OrderedPoint> on_boundary(const Patch& patch, const Boundary& boundary,
                                    std::span<const CaseWingPoint> vertices,
                                    std::span<OrderedPoint> points, Math& m) {
    std::size_t size = 0;
    for (std::size_t i = 0; i <= patch.m(); ++i)
        for (std::size_t j = 0; j <= patch.n(); ++j) {
            m.visit(1);
            const auto index = patch.point(i, j);
            if (auto at = boundary_parameter(boundary.points, vertices[index].chart, m))
                points[size++] = {std::move(*at), index};
        }
    points = points.first(size);
    std::ranges::sort(points,
                      [&](const auto& a, const auto& b) { return m.compare(a.at, b.at) < 0; });
    return points;
}

std::size_t boundary_edge(const Patch& patch, const CaseWingPoint& a, const CaseWingPoint& b) {
    const auto [i, j] = std::pair{std::min(a.chord, b.chord), std::min(a.span, b.span)};
    if (a.chord == b.chord && (i == 0 || i == patch.m()) && std::max(a.span, b.span) == j + 1)
        return patch.span_edge(i, j);
    if (a.span == b.span && (j == 0 || j == patch.n()) && std::max(a.chord, b.chord) == i + 1)
        return patch.chord_edge(i, j);
    partition_fail(GeometryInvalid);
}

} // namespace

std::size_t root(std::span<std::size_t> owner, std::size_t point, Math& m) {
    while (owner[point] != point) {
        m.visit(1);
        owner[point] = owner[owner[point]];
        point = owner[point];
    }
    return point;
}

void unite(std::span<std::size_t> owner, std::size_t a, std::size_t b, Math& m) {
    a = root(owner, a, m);
    b = root(owner, b, m);
    owner[std::max(a, b)] = std::min(a, b);
}

void share(const Patch& lhs, std::uint64_t lhs_id, const Patch& rhs, std::uint64_t rhs_id,
           bool reverse, std::span<const CaseWingPoint> points, std::span<std::size_t> point_owners,
           std::span<std::size_t> edge_owners, std::span<bool> internal, Math& m,
           CaseWingDiagnostic& diagnostic) {
    const Boundary a{*lhs.chart,
                     lhs.chart->boundary(static_cast<std::size_t>(partition_index(lhs_id)))};
    const Boundary b{*rhs.chart,
                     rhs.chart->boundary(static_cast<std::size_t>(partition_index(rhs_id)))};
    diagnostic.strip.reset();
    diagnostic.other_strip.reset();
    prove_shared_arc(a, b, reverse, m, diagnostic);
    diagnostic.strip.reset();
    diagnostic.other_strip.reset();
    qualify_motion(lhs, rhs, m);

    auto space_a = allocate<OrderedPoint>(m, (lhs.m() + 1) * (lhs.n() + 1));
    auto space_b = allocate<OrderedPoint>(m, (rhs.m() + 1) * (rhs.n() + 1));
    const auto pa = on_boundary(lhs, a, points, space_a.values(), m);
    const auto pb = on_boundary(rhs, b, points, space_b.values(), m);
    if (pa.size() < 2 || pa.size() != pb.size())
        partition_fail(GeometryInvalid);
    const auto index_b = [&](std::size_t i) { return pb[reverse ? pb.size() - 1 - i : i].point; };
    for (std::size_t i = 0; i < pa.size(); ++i) {
        m.visit(1);
        const auto ia = pa[i].point, ib = index_b(i);
        const auto& aa = points[ia];
        const auto& bb = points[ib];
        if (!same_xyz(a.world(aa.chart, m), b.world(bb.chart, m), m))
            partition_fail(GeometryInvalid);
        unite(point_owners, ia, ib, m);
        if (i == 0)
            continue;
        const auto ea = boundary_edge(lhs, points[pa[i - 1].point], aa);
        const auto eb = boundary_edge(rhs, points[index_b(i - 1)], bb);
        unite(edge_owners, ea, eb, m);
        internal[ea] = internal[eb] = true;
    }
}

} // namespace uvlas::coupling::detail::wing_geometry
