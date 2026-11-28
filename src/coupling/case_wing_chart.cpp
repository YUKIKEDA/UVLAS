#include "case_wing_chart.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <span>
#include <utility>

namespace uvlas::coupling::detail::wing_chart {
namespace pg = structure::detail;
namespace {
using enum core::FailureCode;
using pg::partition_add;
using pg::partition_fail;
using pg::partition_index;
using pg::partition_mul;
using pg::partition_take;
constexpr auto operation = core::Operation::PrepareModel;

template <typename T>
core::PlannedBuffer<T> allocate(std::uint64_t count, Math& math) {
    return partition_take(math.plan.allocate<T>(partition_index(count), operation));
}

template <typename Field>
Rational number(const io::SiCase& source, Field field, Math& math) {
    const auto quantity = partition_take(source.quantity(field));
    if (quantity.recipe.dimensions.exponents != io::SchemaDimension{}.exponents ||
        quantity.recipe.pi_power != 0)
        partition_fail(InvalidInput);
    const auto digits = std::visit(
        [](const auto& value) -> std::uint64_t {
            if constexpr (std::is_same_v<std::decay_t<decltype(value)>, io::SiDecimalView>)
                return value.digits.size();
            else
                return partition_add(value.numerator.digits.size(),
                                     value.denominator.digits.size());
        },
        quantity.recipe.value);
    math.charge(digits);
    return partition_take(io::materialize_rational(quantity.recipe, math.leaf_limits, math.plan));
}

template <typename Field>
Point point(const io::SiCase& source, Field fields, Math& math) {
    return {number(source, partition_take(fields.at(0)), math),
            number(source, partition_take(fields.at(1)), math)};
}

Point subtract(const Point& a, const Point& b, Math& math) {
    return {math.sub(a.s, b.s), math.sub(a.t, b.t)};
}

int cross_sign(const Point& a, const Point& b, const Point& c, Math& math) {
    const auto ab = subtract(b, a, math);
    const auto ac = subtract(c, a, math);
    return math.sub(math.mul(ab.s, ac.t), math.mul(ab.t, ac.s)).sign();
}

bool equal(const Point& a, const Point& b, Math& math) {
    return math.compare(a.s, b.s) == 0 && math.compare(a.t, b.t) == 0;
}

bool between(const Rational& x, const Rational& a, const Rational& b, Math& math) {
    const auto left = math.compare(a, x);
    const auto right = math.compare(x, b);
    return left == 0 || right == 0 || left == right;
}

Rational ratio(const Point& a, const Point& b, const Point& p, Math& math) {
    if (math.compare(a.t, b.t) != 0)
        return math.div(math.sub(p.t, a.t), math.sub(b.t, a.t));
    return math.div(math.sub(p.s, a.s), math.sub(b.s, a.s));
}

template <typename Field>
core::PlannedBuffer<Rational> grid_axis(const io::SiCase& source, Field values, Math& math) {
    if (values.size() < 2)
        partition_fail(InvalidInput);
    auto result = allocate<Rational>(values.size(), math);
    auto axis = result.values();
    for (std::size_t i = 0; i < axis.size(); ++i) {
        axis[i] = number(source, partition_take(values.at(partition_index(i))), math);
        if (i && math.compare(axis[i - 1], axis[i]) >= 0)
            partition_fail(InvalidInput);
    }
    if (axis.front().sign() != 0 || math.compare(axis.back(), math.constant("1")) != 0)
        partition_fail(InvalidInput);
    return result;
}

template <typename T, typename Project>
std::size_t lower_segment(std::span<const T> values, const Rational& value, Project project,
                          Math& math) {
    const auto upper = std::ranges::upper_bound(
        values, value, [&](const Rational& a, const Rational& b) { return math.compare(a, b) < 0; },
        project);
    return std::min(values.size() - 2, static_cast<std::size_t>(upper - values.begin()) - 1);
}

void validate_span_branch(std::span<const Point> boundary, bool reversed, Math& math) {
    const auto branch = [&](std::size_t i) -> const Point& {
        return boundary[reversed ? boundary.size() - i - 1 : i];
    };
    for (std::size_t i = 1; i < boundary.size(); ++i)
        if (math.compare(branch(i - 1).t, branch(i).t) >= 0)
            partition_fail(InvalidInput);
}

Point at_chart_span(std::span<const Point> boundary, bool reversed, const Rational& t, Math& math) {
    const auto branch = [&](std::size_t i) -> const Point& {
        return boundary[reversed ? boundary.size() - i - 1 : i];
    };
    std::size_t lo = 0, hi = boundary.size() - 1;
    while (hi - lo > 1) {
        const auto mid = lo + (hi - lo) / 2;
        if (math.compare(branch(mid).t, t) <= 0)
            lo = mid;
        else
            hi = mid;
    }
    const auto& a = branch(lo);
    const auto& b = branch(hi);
    const auto u = math.div(math.sub(t, a.t), math.sub(b.t, a.t));
    return {blend(a.s, b.s, u, math), math.copy(t)};
}

struct Parameters {
    Rational eta;
    std::array<Rational, 2> edge;
};

struct KnotPlan {
    core::PlannedBuffer<Knot> knots;
    std::size_t count{};
};

template <typename T>
std::span<T> active_knots(core::PlannedBuffer<T>& values, std::size_t count, Math& math) {
    auto buffer = values.values().first(count);
    std::ranges::sort(
        buffer, [&](const auto& left, const auto& right) { return math.compare(left, right) < 0; });
    const auto unique = std::ranges::unique(buffer, [&](const auto& left, const auto& right) {
        return math.compare(left, right) == 0;
    });
    return buffer.first(static_cast<std::size_t>(unique.begin() - buffer.begin()));
}

KnotPlan chart_span(std::span<const Point> leading, std::span<const Point> trailing, Math& math) {
    const bool left_reversed = math.compare(leading.front().t, leading.back().t) > 0;
    const bool right_reversed = math.compare(trailing.front().t, trailing.back().t) > 0;
    validate_span_branch(leading, left_reversed, math);
    validate_span_branch(trailing, right_reversed, math);
    const auto& low = (left_reversed ? leading.back() : leading.front()).t;
    const auto& high = (left_reversed ? leading.front() : leading.back()).t;
    if (math.compare(low, (right_reversed ? trailing.back() : trailing.front()).t) != 0 ||
        math.compare(high, (right_reversed ? trailing.front() : trailing.back()).t) != 0)
        partition_fail(InvalidInput);
    const auto width = math.sub(high, low);
    auto stations = allocate<Rational>(partition_add(leading.size(), trailing.size()), math);
    std::size_t count = 0;
    for (const auto edge : {leading, trailing})
        for (const auto& station : edge)
            stations.values()[count++] = math.copy(station.t);
    const auto ordered = active_knots(stations, count, math);
    auto knots = allocate<Knot>(ordered.size(), math);
    for (std::size_t i = 0; i < ordered.size(); ++i) {
        const auto& t = ordered[i];
        knots.values()[i] = {math.div(math.sub(t, low), width),
                             at_chart_span(leading, left_reversed, t, math),
                             at_chart_span(trailing, right_reversed, t, math)};
    }
    return {std::move(knots), ordered.size()};
}

Point at_edge_ordinal(std::span<const Point> boundary, const Rational& p, Math& math) {
    std::size_t lo = 0, hi = boundary.size() - 1;
    while (hi - lo > 1) {
        const auto mid = lo + (hi - lo) / 2;
        if (math.compare(integer(mid, math), p) <= 0)
            lo = mid;
        else
            hi = mid;
    }
    return blend(boundary[lo], boundary[hi], math.sub(p, integer(lo, math)), math);
}

KnotPlan explicit_span(const io::SiCase& source, const io::schema::EdgePairing& declaration,
                       std::array<std::span<const Point>, 2> edges, Math& math) {
    const auto size = declaration.knots.size();
    if (size < 2)
        partition_fail(InvalidInput);
    auto records = allocate<Parameters>(size, math);
    const auto entries = records.values();
    auto stations = allocate<Rational>(
        partition_add(size, partition_add(edges[0].size() - 2, edges[1].size() - 2)), math);
    std::size_t count = 0;
    for (std::size_t i = 0; i < entries.size(); ++i) {
        const auto knot = partition_take(declaration.knots.at(partition_index(i))).record();
        entries[i] = {number(source, knot.eta, math),
                      {number(source, knot.leading, math), number(source, knot.trailing, math)}};
        if (i && math.compare(entries[i - 1].eta, entries[i].eta) >= 0)
            partition_fail(InvalidInput);
        stations.values()[count++] = math.copy(entries[i].eta);
    }
    if (entries.front().eta.sign() != 0 ||
        math.compare(entries.back().eta, math.constant("1")) != 0)
        partition_fail(InvalidInput);

    for (std::size_t e = 0; e < edges.size(); ++e) {
        const bool reversed = entries.back().edge[e].sign() == 0;
        const auto& first = (reversed ? entries.back() : entries.front()).edge[e];
        const auto& last = (reversed ? entries.front() : entries.back()).edge[e];
        if (first.sign() != 0 || math.compare(last, integer(edges[e].size() - 1, math)) != 0)
            partition_fail(InvalidInput);
        for (std::size_t i = 1; i < entries.size(); ++i) {
            const auto order = math.compare(entries[i - 1].edge[e], entries[i].edge[e]);
            if (!(reversed ? order > 0 : order < 0))
                partition_fail(InvalidInput);
        }
        std::size_t segment = 0;
        for (std::size_t i = 1; i < edges[e].size() - 1; ++i) {
            const auto vertex = integer(reversed ? edges[e].size() - i - 1 : i, math);
            while (reversed ? math.compare(entries[segment + 1].edge[e], vertex) > 0
                            : math.compare(entries[segment + 1].edge[e], vertex) < 0)
                ++segment;
            const auto& before = entries[segment];
            const auto& after = entries[segment + 1];
            if (math.compare(before.edge[e], vertex) == 0 ||
                math.compare(after.edge[e], vertex) == 0)
                continue;
            const auto fraction =
                math.div(math.sub(vertex, before.edge[e]), math.sub(after.edge[e], before.edge[e]));
            stations.values()[count++] = blend(before.eta, after.eta, fraction, math);
        }
    }
    const auto ordered = active_knots(stations, count, math);
    auto knots = allocate<Knot>(ordered.size(), math);
    const std::span<const Parameters> declared{entries};
    for (std::size_t i = 0; i < ordered.size(); ++i) {
        const auto& eta = ordered[i];
        const auto segment = lower_segment(declared, eta, &Parameters::eta, math);
        const auto& a = declared[segment];
        const auto& b = declared[segment + 1];
        const auto fraction = math.div(math.sub(eta, a.eta), math.sub(b.eta, a.eta));
        knots.values()[i] = {
            math.copy(eta),
            at_edge_ordinal(edges[0], blend(a.edge[0], b.edge[0], fraction, math), math),
            at_edge_ordinal(edges[1], blend(a.edge[1], b.edge[1], fraction, math), math)};
    }
    return {std::move(knots), ordered.size()};
}

using Quad = std::array<const Point*, 4>;
Quad strip(std::span<const Knot> knots, std::size_t i) {
    return {&knots[i].leading, &knots[i].trailing, &knots[i + 1].trailing, &knots[i + 1].leading};
}

bool inside(const Quad& quad, const Point& p, int sign, Math& math) {
    for (std::size_t i = 0; i < quad.size(); ++i)
        if (cross_sign(*quad[i], *quad[(i + 1) % quad.size()], p, math) != sign)
            return false;
    return true;
}

void check_strips(const Quad& first, const Quad& second, bool adjacent, int sign, Math& math) {
    const auto allowed = [&](const Point& contact) {
        return adjacent && on_segment(*first[2], *first[3], contact, math);
    };
    for (std::size_t i = 0; i < first.size(); ++i) {
        const auto& a = *first[i];
        const auto& b = *first[(i + 1) % first.size()];
        if (inside(second, a, sign, math) || inside(first, *second[i], sign, math))
            partition_fail(GeometryInvalid);
        for (std::size_t j = 0; j < second.size(); ++j) {
            const auto& c = *second[j];
            const auto& d = *second[(j + 1) % second.size()];
            const auto ab_c = cross_sign(a, b, c, math), ab_d = cross_sign(a, b, d, math);
            const auto cd_a = cross_sign(c, d, a, math), cd_b = cross_sign(c, d, b, math);
            if ((ab_c * ab_d < 0 && cd_a * cd_b < 0) ||
                (ab_c == 0 && on_segment(a, b, c, math) && !allowed(c)) ||
                (ab_d == 0 && on_segment(a, b, d, math) && !allowed(d)) ||
                (cd_a == 0 && on_segment(c, d, a, math) && !allowed(a)) ||
                (cd_b == 0 && on_segment(c, d, b, math) && !allowed(b)))
                partition_fail(GeometryInvalid);
        }
    }
}

template <typename At>
std::size_t reduce_boundary(std::span<const Point*> output, std::size_t count, At at, Math& math) {
    std::size_t retained = 0;
    for (std::size_t i = 0; i < count; ++i) {
        math.visit(1);
        if (!on_segment(at((i + count - 1) % count), at((i + 1) % count), at(i), math))
            output[retained++] = &at(i);
    }
    return retained;
}

void match_domain(std::span<const Knot> knots, std::span<const Point> domain, Math& math) {
    const auto count = partition_mul(2, knots.size());
    auto outline = allocate<const Point*>(count, math);
    auto target = allocate<const Point*>(domain.size(), math);
    const auto loop_size = reduce_boundary(
        outline.values(), outline.size(),
        [&](std::size_t i) -> const Point& {
            return i < knots.size() ? knots[i].leading : knots[outline.size() - 1 - i].trailing;
        },
        math);
    const auto domain_size = reduce_boundary(
        target.values(), target.size(), [&](std::size_t i) -> const Point& { return domain[i]; },
        math);
    if (loop_size < 3 || loop_size != domain_size)
        partition_fail(GeometryInvalid);
    const auto actual = outline.values().first(loop_size);
    const auto expected = target.values().first(domain_size);
    const auto start = std::ranges::find_if(expected, [&](const Point* coordinate) {
        return equal(*coordinate, *actual.front(), math);
    });
    if (start == expected.end())
        partition_fail(GeometryInvalid);
    const auto offset = static_cast<std::size_t>(start - expected.begin());
    const bool forward = equal(*actual[1], *expected[(offset + 1) % domain_size], math);
    for (std::size_t i = 1; i < loop_size; ++i) {
        const auto at = (offset + (forward ? i : domain_size - i)) % domain_size;
        if (!equal(*actual[i], *expected[at], math))
            partition_fail(GeometryInvalid);
    }
}

void validate_cover(std::span<const Knot> knots, std::span<const Point> domain, Math& math,
                    CaseWingDiagnostic& diagnostic) {
    int sign = 0;
    for (std::size_t k = 0; k + 1 < knots.size(); ++k) {
        diagnostic.strip = k;
        const auto quad = strip(knots, k);
        for (std::size_t i = 0; i < quad.size(); ++i) {
            const auto turn = cross_sign(*quad[i], *quad[(i + 1) % 4], *quad[(i + 2) % 4], math);
            if (!turn || (sign != 0 && sign != turn))
                partition_fail(GeometryInvalid);
            sign = turn;
        }
        // The four turns are the affine Jacobian's four corner signs.
        for (std::size_t i = 0; i < k; ++i) {
            math.visit(1);
            diagnostic.other_strip = i;
            check_strips(strip(knots, i), quad, i + 1 == k, sign, math);
        }
        diagnostic.other_strip.reset();
    }
    diagnostic.strip.reset();
    match_domain(knots, domain, math);
}
} // namespace

Rational integer(std::size_t value, Math& math) {
    std::array<char, 3 * sizeof(value)> text{};
    const auto end = std::to_chars(text.data(), text.data() + text.size(), value).ptr;
    return math.constant({text.data(), static_cast<std::size_t>(end - text.data())});
}

Rational blend(const Rational& a, const Rational& b, const Rational& u, Math& math) {
    return math.add(a, math.mul(u, math.sub(b, a)));
}

Point blend(const Point& a, const Point& b, const Rational& u, Math& math) {
    return {blend(a.s, b.s, u, math), blend(a.t, b.t, u, math)};
}

bool on_segment(const Point& a, const Point& b, const Point& p, Math& math) {
    return cross_sign(a, b, p, math) == 0 && between(p.s, a.s, b.s, math) &&
           between(p.t, a.t, b.t, math);
}

std::optional<Rational> boundary_parameter(std::span<const Point> boundary, const Point& p,
                                           Math& math) {
    for (std::size_t i = 0; i + 1 < boundary.size(); ++i) {
        math.visit(1);
        if (equal(boundary[i], boundary[i + 1], math))
            partition_fail(GeometryInvalid);
        if (on_segment(boundary[i], boundary[i + 1], p, math))
            return math.add(integer(i, math), ratio(boundary[i], boundary[i + 1], p, math));
    }
    return std::nullopt;
}

std::span<const Point> Chart::boundary(std::size_t index) const {
    if (index >= boundaries.size())
        partition_fail(InvalidInput);
    const auto range = boundaries.values()[index];
    return boundary_points.values().subspan(range.offset, range.count);
}

Site Chart::site(const Rational& xi, const Rational& eta, bool derivatives, Math& math,
                 CaseWingDiagnostic& diagnostic) const {
    const auto stations = knots.values().first(knot_count);
    if (math.compare(eta, stations.front().eta) < 0 || math.compare(eta, stations.back().eta) > 0)
        partition_fail(InvalidInput);
    const auto i = lower_segment(stations, eta, &Knot::eta, math);
    diagnostic.strip = i;
    const auto tangent = [&](std::size_t k) {
        const auto& a = stations[k];
        const auto& b = stations[k + 1];
        const auto step = math.sub(b.eta, a.eta);
        const auto selected_a = blend(a.leading, a.trailing, xi, math);
        const auto selected_b = blend(b.leading, b.trailing, xi, math);
        return Point{math.div(math.sub(selected_b.s, selected_a.s), step),
                     math.div(math.sub(selected_b.t, selected_a.t), step)};
    };
    const auto& a = stations[i];
    const auto& b = stations[i + 1];
    const auto fraction = math.div(math.sub(eta, a.eta), math.sub(b.eta, a.eta));
    auto left = blend(a.leading, b.leading, fraction, math);
    auto right = blend(a.trailing, b.trailing, fraction, math);
    auto location = blend(left, right, xi, math);
    if (!derivatives)
        return {std::move(location), {}, {}, i};
    auto span_direction = tangent(i);
    if (i > 0 && math.compare(a.eta, eta) == 0 && !equal(span_direction, tangent(i - 1), math)) {
        diagnostic.other_strip = i - 1;
        partition_fail(GeometryInvalid);
    }
    return {std::move(location), subtract(right, left, math), std::move(span_direction), i};
}

Chart prepare_chart(const io::SiCase& source, const io::schema::AeroSurface& input,
                    const io::schema::GeometrySurfaceGraph& graph, Math& math,
                    CaseWingDiagnostic& diagnostic) {
    diagnostic.phase = CaseWingPhase::Boundaries;
    const auto geometry_id = graph.id.unsigned_value();
    const auto leading = input.leading.unsigned_value(), trailing = input.trailing.unsigned_value();
    diagnostic.boundary = leading;
    diagnostic.other_boundary = trailing;
    if (geometry_id != input.geometry.unsigned_value() || leading == trailing ||
        leading >= graph.boundaries.size() || trailing >= graph.boundaries.size())
        partition_fail(InvalidInput);
    diagnostic.phase = CaseWingPhase::Reference;
    auto graph_limits = math.limits;
    graph_limits.algebra = math.leaf_limits;
    graph_limits.maximum_scalar_operations -= math.scalar_operations;
    graph_limits.maximum_partition_visits -= math.partition_visits;
    auto reference =
        partition_take(CaseReferenceSurface::prepare(source, geometry_id, graph_limits, math.plan));
    const auto graph_work = reference.preparation_work();
    math.charge(graph_work.scalar_operations);
    math.visit(graph_work.partition_visits);
    auto boundary_ranges = allocate<BoundaryRange>(graph.boundaries.size(), math);
    std::uint64_t count = 0;
    for (std::size_t j = 0; j < boundary_ranges.size(); ++j) {
        math.visit(1);
        const auto boundary = partition_take(graph.boundaries.at(partition_index(j))).record();
        const auto offset = count;
        count = partition_add(count, boundary.vertices.size());
        (void)partition_index(count);
        boundary_ranges.values()[j] = {static_cast<std::size_t>(offset),
                                       static_cast<std::size_t>(boundary.vertices.size())};
    }
    auto boundary_points = allocate<Point>(count, math);
    auto domain = allocate<Point>(graph.domain.size(), math);
    for (std::size_t j = 0; j < domain.size(); ++j)
        domain.values()[j] =
            point(source, partition_take(graph.domain.at(partition_index(j))), math);
    for (std::size_t j = 0; j < boundary_ranges.size(); ++j) {
        const auto vertices =
            partition_take(graph.boundaries.at(partition_index(j))).record().vertices;
        const auto range = boundary_ranges.values()[j];
        for (std::size_t i = 0; i < range.count; ++i)
            boundary_points.values()[range.offset + i] =
                point(source, partition_take(vertices.at(partition_index(i))), math);
    }
    diagnostic.phase = CaseWingPhase::Chord;
    auto chord = grid_axis(source, input.chord, math);
    diagnostic.phase = CaseWingPhase::Span;
    auto span = grid_axis(source, input.span, math);
    Chart chart{std::move(reference),
                input.id.unsigned_value(),
                leading,
                trailing,
                input.orientation.text() == "reversed",
                std::move(chord),
                std::move(span),
                {},
                0,
                std::move(boundary_ranges),
                std::move(boundary_points)};
    diagnostic.phase = CaseWingPhase::Boundaries;
    const auto lead = chart.boundary(static_cast<std::size_t>(leading));
    const auto tail = chart.boundary(static_cast<std::size_t>(trailing));
    if (lead.size() < 2 || tail.size() < 2)
        partition_fail(InvalidInput);
    const auto reject_zero = [&](std::span<const Point> edge) {
        for (std::size_t i = 1; i < edge.size(); ++i)
            if (equal(edge[i - 1], edge[i], math))
                partition_fail(GeometryInvalid);
    };
    reject_zero(lead);
    reject_zero(tail);
    diagnostic.phase = CaseWingPhase::Pairing;
    auto prepared = input.edge_pairing
                        ? explicit_span(source, input.edge_pairing->record(), {lead, tail}, math)
                        : chart_span(lead, tail, math);
    diagnostic.phase = CaseWingPhase::Coverage;
    validate_cover(prepared.knots.values().first(prepared.count), domain.values(), math,
                   diagnostic);
    chart.knots = std::move(prepared.knots);
    chart.knot_count = prepared.count;
    return chart;
}
} // namespace uvlas::coupling::detail::wing_chart
