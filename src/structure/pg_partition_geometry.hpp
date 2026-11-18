#pragma once
#include "pg_partition_algebra.hpp"

#include <algorithm>

namespace uvlas::structure::detail {
struct PgPartitionGeometry {
    PgPartitionData& data;
    PgPartitionScratch& scratch;
    PgPartitionAlgebra& op;
    const PgPartitionRequest& input;
    ExactRational zero, one, two, three;
    ExactRational domain_area;
    std::uint64_t raw_line_count = 0, raw_vertex_count = 0;
    std::size_t active_clip = 0;
    std::uint64_t clip_piece_count = 0, clip_edge_count = 0;
    PgPartitionGeometry(PgPartitionData& d, PgPartitionScratch& s, PgPartitionAlgebra& a,
                        const PgPartitionRequest& r)
        : data(d),
          scratch(s),
          op(a),
          input(r),
          one(a.number(1)),
          two(a.number(2)),
          three(a.number(3)) {}
    template <class T>
    T& append(core::PlannedBuffer<T>& array, std::uint64_t& count) {
        op.visit();
        if (count >= array.size())
            partition_fail(core::FailureCode::ResourceExceeded);
        return array.values()[static_cast<std::size_t>(count++)];
    }
    [[nodiscard]] const PgPartitionPoint& original(core::Index i) const {
        return data.source_points.values()[static_cast<std::size_t>(i)];
    }
    [[nodiscard]] const PgPartitionPoint& vertex(core::Index i) const {
        return data.vertices.values()[static_cast<std::size_t>(i)];
    }
    static bool valid_utf8(std::string_view text) noexcept {
        std::uint32_t code = 0, minimum = 0;
        unsigned remaining = 0;
        for (const auto c : text) {
            const auto byte = static_cast<unsigned char>(c);
            if (remaining != 0) {
                if ((byte & 0xc0U) != 0x80U)
                    return false;
                code = (code << 6U) | (byte & 0x3fU);
                if (--remaining == 0 &&
                    (code < minimum || code > 0x10ffffU || (code >= 0xd800U && code <= 0xdfffU)))
                    return false;
            } else if (byte < 0x80U) {
                if (byte == 0)
                    return false;
            } else if (byte >= 0xc2U && byte <= 0xdfU) {
                remaining = 1;
                code = byte & 0x1fU;
                minimum = 0x80U;
            } else if (byte >= 0xe0U && byte <= 0xefU) {
                remaining = 2;
                code = byte & 0x0fU;
                minimum = 0x800U;
            } else if (byte >= 0xf0U && byte <= 0xf4U) {
                remaining = 3;
                code = byte & 0x07U;
                minimum = 0x10000U;
            } else
                return false;
        }
        return remaining == 0;
    }
    void snapshot() {
        const PgGraphLimits gl{op.exact, op.budget.cells, op.caps.cancellation};
        data.graph = partition_take(input.graph.clone(gl, op.payload));
        data.domain_count = input.domain.size();
        data.thickness_offset = input.domain.size();
        std::uint64_t points = 0, sources = 0, text = 0;
        const auto add_point = [&](const PgPartitionPoint& p) {
            append(data.source_points, points) = op.point_clone(p);
        };
        for (const auto& p : input.domain)
            add_point(p);
        for (const auto& p : input.thickness.vertices)
            add_point(p);
        for (std::size_t i = 0; i < input.thickness.triangles.size(); ++i) {
            op.visit();
            data.thickness_triangles.values()[i] = input.thickness.triangles[i];
        }
        for (std::size_t i = 0; i < input.thickness.values.size(); ++i)
            data.thickness_values.values()[i] = op.clone(input.thickness.values[i]);
        for (std::size_t i = 0; i < input.layers.size(); ++i) {
            const auto& x = input.layers[i];
            auto& y = data.layers.values()[i];
            y.id = x.id;
            y.lower = op.clone(x.lower);
            y.upper = op.clone(x.upper);
        }
        for (std::size_t i = 0; i < input.domain.size(); ++i)
            append(data.sources,
                   sources) = {PgPartitionDomainSource{static_cast<core::Index>(i)},
                               {static_cast<core::Index>(i),
                                static_cast<core::Index>((i + 1) % input.domain.size())}};
        for (std::size_t i = 0; i < input.thickness.triangles.size(); ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                const auto& tri = input.thickness.triangles[i];
                append(data.sources,
                       sources) = {PgPartitionThicknessSource{static_cast<core::Index>(i),
                                                              static_cast<core::Index>(j)},
                                   {partition_index(data.thickness_offset) + tri[j],
                                    partition_index(data.thickness_offset) + tri[(j + 1) % 3]}};
            }
        for (std::size_t i = 0; i < input.boundaries.size(); ++i) {
            const auto& b = input.boundaries[i];
            auto& record = data.boundaries.values()[i];
            if (!valid_utf8(b.name))
                partition_fail(core::FailureCode::InvalidInput);
            for (std::size_t j = 0; j < i; ++j) {
                op.visit();
                if (input.boundaries[j].name == b.name)
                    partition_fail(core::FailureCode::InvalidInput);
            }
            record = {partition_index(text), static_cast<core::Index>(b.name.size()),
                      partition_index(points), static_cast<core::Index>(b.vertices.size())};
            for (const auto c : b.name) {
                append(data.text, text) = c;
            }
            for (const auto& p : b.vertices)
                add_point(p);
            for (std::size_t j = 0; j + 1 < b.vertices.size(); ++j)
                append(data.sources,
                       sources) = {PgPartitionNamedSource{static_cast<core::Index>(i),
                                                          static_cast<core::Index>(j)},
                                   {record.first_point + static_cast<core::Index>(j),
                                    record.first_point + static_cast<core::Index>(j + 1)}};
        }
        for (std::size_t i = 0; i < input.refinements.size(); ++i) {
            const auto& r = input.refinements[i];
            const auto first = partition_index(points);
            add_point(r.first);
            add_point(r.second);
            append(data.sources, sources) = {
                PgPartitionRefinementSource{static_cast<core::Index>(i), r.id, r.version},
                {first, first + 1}};
        }
        if (points != data.source_points.size() || sources != data.sources.size() ||
            text != data.text.size())
            partition_fail(core::FailureCode::InvalidInput);
    }
    bool intersects(const PgPartitionSource& x, const PgPartitionSource& y) {
        const auto& a = original(x.endpoints[0]);
        const auto& b = original(x.endpoints[1]);
        const auto& c = original(y.endpoints[0]);
        const auto& d = original(y.endpoints[1]);
        const auto ac = op.orientation(a, b, c).sign(), ad = op.orientation(a, b, d).sign();
        const auto ca = op.orientation(c, d, a).sign(), cb = op.orientation(c, d, b).sign();
        return (ac * ad < 0 && ca * cb < 0) || (ac == 0 && op.on_segment(a, b, c)) ||
               (ad == 0 && op.on_segment(a, b, d)) || (ca == 0 && op.on_segment(c, d, a)) ||
               (cb == 0 && op.on_segment(c, d, b));
    }
    bool positive_overlap(const PgPartitionSource& x, const PgPartitionSource& y) {
        const auto& a = original(x.endpoints[0]);
        const auto& b = original(x.endpoints[1]);
        const auto& c = original(y.endpoints[0]);
        const auto& d = original(y.endpoints[1]);
        if (op.orientation(a, b, c).sign() != 0 || op.orientation(a, b, d).sign() != 0)
            return false;
        const auto& av = a.s == b.s ? a.t : a.s;
        const auto& bv = a.s == b.s ? b.t : b.s;
        const auto& cv = a.s == b.s ? c.t : c.s;
        const auto& dv = a.s == b.s ? d.t : d.s;
        const auto& low1 = op.compare(av, bv) < 0 ? av : bv;
        const auto& hi1 = op.compare(av, bv) < 0 ? bv : av;
        const auto& low2 = op.compare(cv, dv) < 0 ? cv : dv;
        const auto& hi2 = op.compare(cv, dv) < 0 ? dv : cv;
        return op.compare(low1, hi2) < 0 && op.compare(low2, hi1) < 0;
    }
    void validate_geometry() {
        for (std::size_t i = 0; i < input.domain.size(); ++i) {
            const auto& a = original(static_cast<core::Index>(i));
            if (op.compare(a.s, zero) < 0 || op.compare(a.s, one) > 0 ||
                op.compare(a.t, zero) < 0 || op.compare(a.t, one) > 0)
                partition_fail(core::FailureCode::InvalidInput);
        }
        for (const auto& source : data.sources.values()) {
            op.visit();
            if (partition_equal_point(original(source.endpoints[0]), original(source.endpoints[1])))
                partition_fail(
                    std::holds_alternative<PgPartitionRefinementSource>(source.descriptor)
                        ? core::FailureCode::InvalidInput
                        : core::FailureCode::GeometryInvalid);
        }
        const auto ds = data.sources.values().first(static_cast<std::size_t>(data.domain_count));
        for (std::size_t i = 0; i < ds.size(); ++i)
            for (std::size_t j = i + 1; j < ds.size(); ++j) {
                op.visit();
                const auto adjacent = j == i + 1 || (i == 0 && j + 1 == ds.size());
                if (adjacent ? positive_overlap(ds[i], ds[j]) : intersects(ds[i], ds[j]))
                    partition_fail(core::FailureCode::GeometryInvalid);
            }
        const auto area2 = polygon_area2_original();
        if (area2.sign() == 0)
            partition_fail(core::FailureCode::GeometryInvalid);
        domain_area = area2.sign() > 0 ? op.div(area2, two) : op.div(op.sub(zero, area2), two);
        for (const auto& tri : data.thickness_triangles.values()) {
            op.visit();
            if (op.orientation(original(partition_index(data.thickness_offset) + tri[0]),
                               original(partition_index(data.thickness_offset) + tri[1]),
                               original(partition_index(data.thickness_offset) + tri[2]))
                    .sign() == 0)
                partition_fail(core::FailureCode::GeometryInvalid);
        }
        const auto minus = op.number(-1);
        if (!(data.layers.values().front().lower == minus) ||
            !(data.layers.values().back().upper == one))
            partition_fail(core::FailureCode::InvalidInput);
        for (std::size_t i = 0; i < data.layers.size(); ++i) {
            const auto& x = data.layers.values()[i];
            if (op.compare(x.lower, x.upper) >= 0 ||
                (i != 0 && !(data.layers.values()[i - 1].upper == x.lower)))
                partition_fail(core::FailureCode::InvalidInput);
            for (std::size_t j = 0; j < i; ++j) {
                op.visit();
                if (data.layers.values()[j].id == x.id)
                    partition_fail(core::FailureCode::InvalidInput);
            }
        }
        for (std::size_t i = 0; i < data.sources.size(); ++i)
            if (const auto* x =
                    std::get_if<PgPartitionNamedSource>(&data.sources.values()[i].descriptor))
                for (std::size_t j = 0; j < i; ++j)
                    if (const auto* y = std::get_if<PgPartitionNamedSource>(
                            &data.sources.values()[j].descriptor)) {
                        op.visit();
                        if (x->boundary == y->boundary &&
                            positive_overlap(data.sources.values()[i], data.sources.values()[j]))
                            partition_fail(core::FailureCode::GeometryInvalid);
                    }
    }
    ExactRational polygon_area2_original() {
        ExactRational total;
        for (std::size_t i = 0; i < data.domain_count; ++i) {
            const auto& a = original(static_cast<core::Index>(i));
            const auto& b = original(static_cast<core::Index>((i + 1) % data.domain_count));
            const auto det = op.determinant(a.s, a.t, b.s, b.t);
            total = op.add(total, det);
        }
        return total;
    }
    void add_axis(const ExactRational& value, bool first) {
        auto& line = append(scratch.raw_lines, raw_line_count);
        const auto& raw = ExactFactory::value(value).backend();
        line.coefficients[first ? 0 : 1] = op.publish(Rational{Integer{raw.denom()}});
        const Integer negative = -Integer{raw.num()};
        line.coefficients[2] = op.publish(Rational{negative});
    }
    int line_order(const PgPartitionLine& a, const PgPartitionLine& b) {
        for (std::size_t i = 0; i < 3; ++i) {
            const auto cmp = op.compare(a.coefficients[i], b.coefficients[i]);
            if (cmp != 0)
                return cmp;
        }
        return 0;
    }
    void make_lines() {
        for (std::size_t i = 0; i < data.sources.size(); ++i) {
            data.source_lines.values()[i] = partition_index(raw_line_count);
            const auto& source = data.sources.values()[i];
            append(scratch.raw_lines, raw_line_count) =
                op.line(original(source.endpoints[0]), original(source.endpoints[1]));
        }
        for (const auto& p : data.source_points.values()) {
            op.visit();
            add_axis(p.s, true);
            add_axis(p.t, false);
        }
        const auto graph = partition_graph(data);
        for (core::Index j = 0; j <= 2 * graph.n; ++j) {
            const auto x = op.number(j, 2 * graph.n);
            add_axis(x, true);
            add_axis(x, false);
        }
        auto order = scratch.order.values().first(static_cast<std::size_t>(raw_line_count));
        for (std::size_t i = 0; i < order.size(); ++i)
            order[i] = static_cast<core::Index>(i);
        partition_sort(op, order, scratch.merge.values(), [&](core::Index a, core::Index b) {
            return line_order(scratch.raw_lines.values()[static_cast<std::size_t>(a)],
                              scratch.raw_lines.values()[static_cast<std::size_t>(b)]) < 0;
        });
        for (const auto raw : order) {
            op.visit();
            auto& source = scratch.raw_lines.values()[static_cast<std::size_t>(raw)];
            if (data.line_count == 0 ||
                line_order(data.lines.values()[static_cast<std::size_t>(data.line_count - 1)],
                           source) != 0)
                append(data.lines, data.line_count) = std::move(source);
            scratch.remap.values()[static_cast<std::size_t>(raw)] =
                partition_index(data.line_count - 1);
        }
        for (auto& index : data.source_lines.values())
            index = scratch.remap.values()[static_cast<std::size_t>(index)];
    }
    [[nodiscard]] std::size_t pair_index(core::Index a, core::Index b) const {
        if (a > b)
            std::swap(a, b);
        const auto x = static_cast<std::uint64_t>(a), y = static_cast<std::uint64_t>(b),
                   m = data.line_count;
        return static_cast<std::size_t>(x * (2 * m - x - 1) / 2 + y - x - 1);
    }
    void make_vertices() {
        for (auto& x : data.pairs.values())
            x = -1;
        for (core::Index i = 0; i < partition_index(data.line_count); ++i)
            for (core::Index j = i + 1; j < partition_index(data.line_count); ++j) {
                op.visit();
                const auto& a = data.lines.values()[static_cast<std::size_t>(i)];
                const auto& b = data.lines.values()[static_cast<std::size_t>(j)];
                if (op.determinant(a.coefficients[0], a.coefficients[1], b.coefficients[0],
                                   b.coefficients[1])
                        .sign() == 0)
                    continue;
                auto p = op.intersection(a, b);
                if (op.compare(p.s, zero) < 0 || op.compare(p.s, one) > 0 ||
                    op.compare(p.t, zero) < 0 || op.compare(p.t, one) > 0)
                    continue;
                data.pairs.values()[pair_index(i, j)] = partition_index(raw_vertex_count);
                append(scratch.raw_vertices, raw_vertex_count) = std::move(p);
            }
        auto order = scratch.order.values().first(static_cast<std::size_t>(raw_vertex_count));
        for (std::size_t i = 0; i < order.size(); ++i)
            order[i] = static_cast<core::Index>(i);
        partition_sort(op, order, scratch.merge.values(), [&](core::Index a, core::Index b) {
            return op.point_order(scratch.raw_vertices.values()[static_cast<std::size_t>(a)],
                                  scratch.raw_vertices.values()[static_cast<std::size_t>(b)]) < 0;
        });
        for (const auto raw : order) {
            op.visit();
            auto& source = scratch.raw_vertices.values()[static_cast<std::size_t>(raw)];
            if (data.vertex_count == 0 ||
                !partition_equal_point(
                    data.vertices.values()[static_cast<std::size_t>(data.vertex_count - 1)],
                    source))
                append(data.vertices, data.vertex_count) = std::move(source);
            scratch.remap.values()[static_cast<std::size_t>(raw)] =
                partition_index(data.vertex_count - 1);
        }
        for (auto& index : data.pairs.values())
            if (index >= 0)
                index = scratch.remap.values()[static_cast<std::size_t>(index)];
    }
    core::Index find_corner(bool s, bool t) {
        for (std::size_t i = 0; i < data.vertex_count; ++i) {
            op.visit();
            const auto& p = data.vertices.values()[i];
            if (p.s == (s ? one : zero) && p.t == (t ? one : zero))
                return static_cast<core::Index>(i);
        }
        partition_fail(core::FailureCode::GeometryInvalid);
    }
    core::Index find_line(core::Index a, core::Index b) {
        for (std::size_t i = 0; i < data.line_count; ++i) {
            op.visit();
            if (op.line_value(data.lines.values()[i], vertex(a)).sign() == 0 &&
                op.line_value(data.lines.values()[i], vertex(b)).sign() == 0)
                return static_cast<core::Index>(i);
        }
        partition_fail(core::FailureCode::GeometryInvalid);
    }
    void clip_half(std::span<const PgPartitionClipEdge> polygon, core::Index cut, int sign,
                   std::size_t output_side) {
        const auto first = clip_edge_count;
        for (std::size_t i = 0; i < polygon.size(); ++i) {
            const auto& a = polygon[i];
            const auto b = polygon[(i + 1) % polygon.size()].vertex;
            const auto sa = sign * op.line_value(data.lines.values()[static_cast<std::size_t>(cut)],
                                                 vertex(a.vertex))
                                       .sign();
            const auto sb =
                sign *
                op.line_value(data.lines.values()[static_cast<std::size_t>(cut)], vertex(b)).sign();
            if (sa >= 0)
                append(scratch.clip[output_side],
                       clip_edge_count) = {a.vertex, sa == 0 && sb < 0 ? cut : a.line};
            if (sa * sb < 0) {
                const auto p = data.pairs.values()[pair_index(a.line, cut)];
                if (p < 0)
                    partition_fail(core::FailureCode::GeometryInvalid);
                append(scratch.clip[output_side], clip_edge_count) = {p, sa < 0 ? a.line : cut};
            }
        }
        const auto count = clip_edge_count - first;
        if (count < 3)
            partition_fail(core::FailureCode::GeometryInvalid);
        append(scratch.pieces[output_side], clip_piece_count) = {first, count};
    }
    void arrange() {
        const std::array corners{find_corner(false, false), find_corner(true, false),
                                 find_corner(true, true), find_corner(false, true)};
        for (std::size_t i = 0; i < 4; ++i)
            scratch.clip[0].values()[i] = {corners[i], find_line(corners[i], corners[(i + 1) % 4])};
        scratch.pieces[0].values()[0] = {0, 4};
        clip_piece_count = 1;
        clip_edge_count = 4;
        for (core::Index cut = 0; cut < partition_index(data.line_count); ++cut) {
            const auto next = 1 - active_clip;
            const auto old_count = clip_piece_count;
            clip_piece_count = 0;
            clip_edge_count = 0;
            for (std::size_t i = 0; i < old_count; ++i) {
                op.visit();
                const auto old = scratch.pieces[active_clip].values()[i];
                const auto polygon = scratch.clip[active_clip].values().subspan(
                    static_cast<std::size_t>(old.first), static_cast<std::size_t>(old.count));
                bool positive = false, negative = false;
                for (const auto& edge : polygon) {
                    op.visit();
                    const auto value =
                        op.line_value(data.lines.values()[static_cast<std::size_t>(cut)],
                                      vertex(edge.vertex))
                            .sign();
                    positive = positive || value > 0;
                    negative = negative || value < 0;
                }
                if (positive && negative) {
                    clip_half(polygon, cut, 1, next);
                    clip_half(polygon, cut, -1, next);
                } else {
                    const auto first = clip_edge_count;
                    for (const auto edge : polygon)
                        append(scratch.clip[next], clip_edge_count) = edge;
                    append(scratch.pieces[next], clip_piece_count) = {first, old.count};
                }
            }
            active_clip = next;
        }
    }
    bool in_domain(const PgPartitionPoint& p) {
        core::Index winding = 0;
        for (std::size_t i = 0; i < data.domain_count; ++i) {
            op.visit();
            const auto& a = original(static_cast<core::Index>(i));
            const auto& b = original(static_cast<core::Index>((i + 1) % data.domain_count));
            const auto ay = op.compare(a.t, p.t), by = op.compare(b.t, p.t),
                       side = op.orientation(a, b, p).sign();
            if (side == 0)
                partition_fail(core::FailureCode::GeometryInvalid);
            if (ay <= 0 && by > 0 && side > 0)
                ++winding;
            else if (by <= 0 && ay > 0 && side < 0)
                --winding;
        }
        return winding != 0;
    }
    PgPartitionPoint interior(std::span<const PgPartitionClipEdge> polygon) {
        for (std::size_t i = 0; i < polygon.size(); ++i) {
            op.visit();
            const auto& a = vertex(polygon[i].vertex);
            const auto& b = vertex(polygon[(i + 1) % polygon.size()].vertex);
            const auto& c = vertex(polygon[(i + 2) % polygon.size()].vertex);
            if (op.orientation(a, b, c).sign() != 0) {
                const auto ab_s = op.add(a.s, b.s), ab_t = op.add(a.t, b.t);
                return {op.div(op.add(ab_s, c.s), three), op.div(op.add(ab_t, c.t), three)};
            }
        }
        partition_fail(core::FailureCode::GeometryInvalid);
    }
    core::Index containing_thickness(const PgPartitionPoint& p) {
        core::Index selected = -1;
        for (std::size_t i = 0; i < data.thickness_triangles.size(); ++i) {
            op.visit();
            const auto& tri = data.thickness_triangles.values()[i];
            const auto at = [&](std::size_t j) -> const PgPartitionPoint& {
                return original(partition_index(data.thickness_offset) + tri[j]);
            };
            const auto a = op.orientation(at(0), at(1), p).sign(),
                       b = op.orientation(at(1), at(2), p).sign(),
                       c = op.orientation(at(2), at(0), p).sign();
            if ((a > 0 && b > 0 && c > 0) || (a < 0 && b < 0 && c < 0)) {
                if (selected >= 0)
                    partition_fail(core::FailureCode::GeometryInvalid);
                selected = static_cast<core::Index>(i);
            }
        }
        if (selected < 0)
            partition_fail(core::FailureCode::GeometryInvalid);
        return selected;
    }
    PgCellIndex containing_cell(const PgPartitionPoint& p, core::Index axis) {
        PgCellIndex found{-1, -1};
        for (core::Index i = 0; i < axis; ++i) {
            const auto low = op.number(i, axis), high = op.number(i + 1, axis);
            if (op.compare(p.s, low) > 0 && op.compare(p.s, high) < 0)
                found.first = i;
            if (op.compare(p.t, low) > 0 && op.compare(p.t, high) < 0)
                found.second = i;
        }
        if (found.first < 0 || found.second < 0)
            partition_fail(core::FailureCode::GeometryInvalid);
        return found;
    }
    void retain_pieces() {
        const auto graph = partition_graph(data);
        const auto axis = 2 * graph.n;
        for (core::Index i = 0; i < axis; ++i)
            for (core::Index j = 0; j < axis; ++j)
                data.cells.values()[static_cast<std::size_t>(i * axis + j)].cell = {i, j};
        for (std::size_t i = 0; i < clip_piece_count; ++i) {
            op.visit();
            const auto clipped = scratch.pieces[active_clip].values()[i];
            const auto polygon = scratch.clip[active_clip].values().subspan(
                static_cast<std::size_t>(clipped.first), static_cast<std::size_t>(clipped.count));
            const auto p = interior(polygon);
            if (!in_domain(p))
                continue;
            PgPartitionPiece result;
            result.cell = containing_cell(p, axis);
            result.thickness_triangle = containing_thickness(p);
            result.first_vertex = partition_index(data.piece_vertex_count);
            for (std::size_t j = 0; j < polygon.size(); ++j) {
                const auto a = polygon[j].vertex, b = polygon[(j + 1) % polygon.size()].vertex;
                auto order = scratch.order.values();
                std::size_t count = 0;
                for (std::size_t k = 0; k < data.vertex_count; ++k) {
                    op.visit();
                    if (op.on_segment(vertex(a), vertex(b), data.vertices.values()[k]))
                        order[count++] = static_cast<core::Index>(k);
                }
                auto selected = order.first(count);
                partition_sort(op, selected, scratch.merge.values(),
                               [&](core::Index x, core::Index y) { return a < b ? x < y : y < x; });
                if (count < 2 || selected.front() != a || selected.back() != b)
                    partition_fail(core::FailureCode::GeometryInvalid);
                for (std::size_t k = 0; k + 1 < count; ++k)
                    append(data.piece_vertices, data.piece_vertex_count) = selected[k];
            }
            result.vertex_count = partition_index(data.piece_vertex_count) - result.first_vertex;
            auto indices =
                data.piece_vertices.values().subspan(static_cast<std::size_t>(result.first_vertex),
                                                     static_cast<std::size_t>(result.vertex_count));
            std::rotate(indices.begin(), std::min_element(indices.begin(), indices.end()),
                        indices.end());
            ExactRational total;
            for (std::size_t j = 0; j < indices.size(); ++j) {
                const auto& a = vertex(indices[j]);
                const auto& b = vertex(indices[(j + 1) % indices.size()]);
                const auto d = op.determinant(a.s, a.t, b.s, b.t);
                total = op.add(total, d);
            }
            if (total.sign() <= 0)
                partition_fail(core::FailureCode::GeometryInvalid);
            result.area = op.div(total, two);
            auto& cell = data.cells.values()[static_cast<std::size_t>(result.cell.first * axis +
                                                                      result.cell.second)];
            cell.area = op.add(cell.area, result.area);
            cell.active = true;
            append(data.pieces, data.piece_count) = std::move(result);
        }
        ExactRational total;
        for (const auto& cell : data.cells.values())
            total = op.add(total, cell.area);
        if (!(total == domain_area) || data.piece_count == 0)
            partition_fail(core::FailureCode::GeometryInvalid);
    }
};
} // namespace uvlas::structure::detail
