#pragma once
#include "pg_partition_geometry.hpp"

#include <tuple>

namespace uvlas::structure::detail {
struct PgPartitionTopology {
    PgPartitionGeometry& g;
    PgPartitionData& d;
    PgPartitionScratch& s;
    PgPartitionAlgebra& op;
    explicit PgPartitionTopology(PgPartitionGeometry& geometry)
        : g(geometry), d(g.data), s(g.scratch), op(g.op) {}
    ExactRational height(core::Index triangle, const PgPartitionPoint& p) {
        const auto& ids = d.thickness_triangles.values()[static_cast<std::size_t>(triangle)];
        const auto at = [&](std::size_t i) -> const PgPartitionPoint& {
            return g.original(partition_index(d.thickness_offset) + ids[i]);
        };
        const auto denominator = op.orientation(at(0), at(1), at(2));
        ExactRational result;
        for (std::size_t i = 0; i < 3; ++i) {
            const auto area = op.orientation(at((i + 1) % 3), at((i + 2) % 3), p);
            const auto weight = op.div(area, denominator);
            const auto term =
                op.mul(weight, d.thickness_values.values()[static_cast<std::size_t>(ids[i])]);
            result = op.add(result, term);
        }
        return result;
    }
    void make_edges() {
        std::uint64_t count = 0;
        for (std::size_t i = 0; i < d.piece_count; ++i) {
            const auto& piece = d.pieces.values()[i];
            const auto vertices =
                d.piece_vertices.values().subspan(static_cast<std::size_t>(piece.first_vertex),
                                                  static_cast<std::size_t>(piece.vertex_count));
            for (std::size_t j = 0; j < vertices.size(); ++j) {
                const auto a = vertices[j], b = vertices[(j + 1) % vertices.size()];
                if (a == b)
                    partition_fail(core::FailureCode::GeometryInvalid);
                g.append(s.occurrences, count) = {std::min(a, b), std::max(a, b),
                                                  static_cast<core::Index>(i), a < b};
            }
        }
        auto order = s.order.values().first(static_cast<std::size_t>(count));
        for (std::size_t i = 0; i < order.size(); ++i)
            order[i] = static_cast<core::Index>(i);
        const auto key = [&](core::Index i) {
            const auto& x = s.occurrences.values()[static_cast<std::size_t>(i)];
            return std::tuple{x.first, x.second, x.piece};
        };
        partition_sort(op, order, s.merge.values(),
                       [&](core::Index a, core::Index b) { return key(a) < key(b); });
        for (std::size_t i = 0; i < order.size();) {
            const auto first = s.occurrences.values()[static_cast<std::size_t>(order[i])];
            PgPartitionEdge edge;
            edge.vertices = {first.first, first.second};
            edge.first_source = partition_index(d.edge_source_count);
            std::size_t end = i;
            while (end < order.size()) {
                op.visit();
                const auto& occurrence =
                    s.occurrences.values()[static_cast<std::size_t>(order[end])];
                if (occurrence.first != first.first || occurrence.second != first.second)
                    break;
                auto& side = occurrence.forward ? edge.left : edge.right;
                if (side)
                    partition_fail(core::FailureCode::GeometryInvalid);
                side = occurrence.piece;
                ++end;
            }
            if (end - i > 2)
                partition_fail(core::FailureCode::GeometryInvalid);
            bool domain = false;
            for (std::size_t j = 0; j < d.sources.size(); ++j) {
                op.visit();
                const auto& source = d.sources.values()[j];
                if (op.on_segment(g.original(source.endpoints[0]), g.original(source.endpoints[1]),
                                  g.vertex(first.first)) &&
                    op.on_segment(g.original(source.endpoints[0]), g.original(source.endpoints[1]),
                                  g.vertex(first.second))) {
                    g.append(d.edge_sources, d.edge_source_count) = static_cast<core::Index>(j);
                    domain = domain ||
                             std::holds_alternative<PgPartitionDomainSource>(source.descriptor);
                }
            }
            edge.source_count = partition_index(d.edge_source_count) - edge.first_source;
            const auto one_side = edge.left.has_value() != edge.right.has_value();
            if (one_side != domain)
                partition_fail(core::FailureCode::GeometryInvalid);
            edge.external = one_side && domain;
            if (edge.left && edge.right) {
                const auto left =
                    d.pieces.values()[static_cast<std::size_t>(*edge.left)].thickness_triangle;
                const auto right =
                    d.pieces.values()[static_cast<std::size_t>(*edge.right)].thickness_triangle;
                if (left != right)
                    for (const auto vertex : edge.vertices) {
                        const auto a = height(left, g.vertex(vertex)),
                                   b = height(right, g.vertex(vertex));
                        if (!(a == b))
                            partition_fail(core::FailureCode::GeometryInvalid);
                    }
            }
            g.append(d.edges, d.edge_count) = edge;
            i = end;
        }
    }
    core::Index root(core::Index i) {
        while (s.queue.values()[static_cast<std::size_t>(i)] != i) {
            op.visit();
            const auto parent = s.queue.values()[static_cast<std::size_t>(i)];
            s.queue.values()[static_cast<std::size_t>(i)] =
                s.queue.values()[static_cast<std::size_t>(parent)];
            i = parent;
        }
        return i;
    }
    void components() {
        for (std::size_t i = 0; i < d.piece_count; ++i) {
            op.visit();
            s.queue.values()[i] = static_cast<core::Index>(i);
            s.remap.values()[i] = -1;
        }
        for (std::size_t i = 0; i < d.edge_count; ++i) {
            op.visit();
            const auto& edge = d.edges.values()[i];
            if (!edge.left || !edge.right)
                continue;
            const auto& left = d.pieces.values()[static_cast<std::size_t>(*edge.left)];
            const auto& right = d.pieces.values()[static_cast<std::size_t>(*edge.right)];
            if (left.cell != right.cell)
                continue;
            const auto a = root(*edge.left), b = root(*edge.right);
            s.queue.values()[static_cast<std::size_t>(std::max(a, b))] = std::min(a, b);
        }
        const auto graph = partition_graph(d);
        const auto axis = 2 * graph.n;
        for (std::size_t i = 0; i < d.piece_count; ++i) {
            op.visit();
            auto& piece = d.pieces.values()[i];
            const auto parent = root(static_cast<core::Index>(i));
            auto& component = s.remap.values()[static_cast<std::size_t>(parent)];
            auto& cell = d.cells.values()[static_cast<std::size_t>(piece.cell.first * axis +
                                                                   piece.cell.second)];
            if (component < 0)
                component = partition_index(cell.components++);
            piece.component = component;
        }
    }
    void make_domains() {
        for (std::size_t i = 0; i < d.piece_count; ++i) {
            const auto& p = d.pieces.values()[i];
            const auto vertices = d.piece_vertices.values().subspan(
                static_cast<std::size_t>(p.first_vertex), static_cast<std::size_t>(p.vertex_count));
            for (std::size_t j = 1; j + 1 < vertices.size(); ++j) {
                op.visit();
                auto area = op.orientation(g.vertex(vertices[0]), g.vertex(vertices[j]),
                                           g.vertex(vertices[j + 1]));
                if (area.sign() < 0)
                    partition_fail(core::FailureCode::GeometryInvalid);
                if (area.sign() == 0)
                    continue;
                const auto triangle = partition_index(d.triangle_count);
                g.append(d.triangles,
                         d.triangle_count) = {static_cast<core::Index>(i),
                                              {vertices[0], vertices[j], vertices[j + 1]},
                                              std::move(area)};
                for (std::size_t layer = 0; layer < d.layers.size(); ++layer) {
                    g.append(d.volumes, d.volume_count) = {triangle,
                                                           static_cast<core::Index>(layer)};
                    g.append(d.horizontal, d.horizontal_count) = {
                        triangle, static_cast<core::Index>(layer), PgPartitionFaceSide::Lower};
                    g.append(d.horizontal, d.horizontal_count) = {
                        triangle, static_cast<core::Index>(layer), PgPartitionFaceSide::Upper};
                }
            }
        }
        for (std::size_t i = 0; i < d.edge_count; ++i) {
            const auto& edge = d.edges.values()[i];
            if (!edge.external)
                continue;
            // Resolve ownership here; the producer's invariant is not a local presence check.
            const auto adjacent_piece = [&]() {
                if (edge.left)
                    return *edge.left;
                if (edge.right)
                    return *edge.right;
                partition_fail(core::FailureCode::GeometryInvalid);
            }();
            for (std::size_t layer = 0; layer < d.layers.size(); ++layer)
                g.append(d.lateral, d.lateral_count) = {static_cast<core::Index>(i),
                                                        static_cast<core::Index>(layer),
                                                        edge.left.has_value()};
            const auto provenance =
                d.edge_sources.values().subspan(static_cast<std::size_t>(edge.first_source),
                                                static_cast<std::size_t>(edge.source_count));
            for (const auto source : provenance)
                if (std::holds_alternative<PgPartitionNamedSource>(
                        d.sources.values()[static_cast<std::size_t>(source)].descriptor)) {
                    const auto& original = d.sources.values()[static_cast<std::size_t>(source)];
                    const auto forward = op.point_order(g.original(original.endpoints[0]),
                                                        g.original(original.endpoints[1])) < 0;
                    g.append(d.traces, d.trace_count) = {static_cast<core::Index>(i), source,
                                                         adjacent_piece, forward};
                }
        }
        for (std::size_t i = 0; i < d.sources.size(); ++i) {
            const auto& source = d.sources.values()[i];
            if (!std::holds_alternative<PgPartitionNamedSource>(source.descriptor))
                continue;
            const auto& a = g.original(source.endpoints[0]);
            const auto& b = g.original(source.endpoints[1]);
            const auto vertical = a.s == b.s;
            const auto& x = vertical ? a.t : a.s;
            const auto& y = vertical ? b.t : b.s;
            auto expected = op.sub(x, y);
            if (expected.sign() < 0)
                expected = op.sub(g.zero, expected);
            ExactRational total;
            for (std::size_t j = 0; j < d.trace_count; ++j)
                if (d.traces.values()[j].source == static_cast<core::Index>(i)) {
                    const auto& edge =
                        d.edges.values()[static_cast<std::size_t>(d.traces.values()[j].edge)];
                    const auto& e = g.vertex(edge.vertices[0]);
                    const auto& f = g.vertex(edge.vertices[1]);
                    auto length = op.sub(vertical ? e.t : e.s, vertical ? f.t : f.s);
                    if (length.sign() < 0)
                        length = op.sub(g.zero, length);
                    total = op.add(total, length);
                }
            if (!(total == expected))
                partition_fail(core::FailureCode::GeometryInvalid);
        }
    }
    void contacts() {
        const auto graph = partition_graph(d);
        const auto axis = 2 * graph.n;
        const auto lookup_count = partition_mul(d.cells.size(), d.domain_count);
        if (lookup_count > s.queue.size())
            partition_fail(core::FailureCode::ResourceExceeded);
        std::fill_n(s.queue.values().begin(), static_cast<std::size_t>(lookup_count),
                    core::Index{-1});
        for (std::size_t i = 0; i < d.cells.size(); ++i) {
            const auto& cell = d.cells.values()[i];
            if (cell.active)
                continue;
            const auto s0 = op.number(cell.cell.first, axis),
                       s1 = op.number(cell.cell.first + 1, axis),
                       t0 = op.number(cell.cell.second, axis),
                       t1 = op.number(cell.cell.second + 1, axis);
            const auto contains = [&](const PgPartitionPoint& p) {
                return op.compare(p.s, s0) >= 0 && op.compare(p.s, s1) <= 0 &&
                       op.compare(p.t, t0) >= 0 && op.compare(p.t, t1) <= 0;
            };
            for (std::size_t j = 0; j < d.edge_count; ++j) {
                op.visit();
                const auto& edge = d.edges.values()[j];
                if (!edge.external)
                    continue;
                const auto a = contains(g.vertex(edge.vertices[0])),
                           b = contains(g.vertex(edge.vertices[1]));
                if (!a && !b)
                    continue;
                const std::array vertices{a ? edge.vertices[0] : edge.vertices[1],
                                          b ? edge.vertices[1] : edge.vertices[0]};
                const auto sources =
                    d.edge_sources.values().subspan(static_cast<std::size_t>(edge.first_source),
                                                    static_cast<std::size_t>(edge.source_count));
                for (const auto source : sources) {
                    op.visit();
                    if (!std::holds_alternative<PgPartitionDomainSource>(
                            d.sources.values()[static_cast<std::size_t>(source)].descriptor))
                        continue;
                    auto& slot = s.queue.values()[i * static_cast<std::size_t>(d.domain_count) +
                                                  static_cast<std::size_t>(source)];
                    if (slot < 0) {
                        slot = partition_index(d.contact_count);
                        g.append(d.contacts,
                                 d.contact_count) = {cell.cell,
                                                     source,
                                                     {std::min(vertices[0], vertices[1]),
                                                      std::max(vertices[0], vertices[1])}};
                    } else {
                        auto& contact = d.contacts.values()[static_cast<std::size_t>(slot)];
                        contact.vertices[0] =
                            std::min({contact.vertices[0], vertices[0], vertices[1]});
                        contact.vertices[1] =
                            std::max({contact.vertices[1], vertices[0], vertices[1]});
                    }
                }
            }
        }
    }
};
} // namespace uvlas::structure::detail
