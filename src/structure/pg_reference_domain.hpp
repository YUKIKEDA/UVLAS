#pragma once

#include "pg_reference_polynomials.hpp"

#include <array>
#include <span>
#include <type_traits>

namespace uvlas::structure::detail {

template <class T>
[[nodiscard]] const T& reference_at(std::span<const T> values, core::Index index) {
    if (index < 0 || static_cast<std::uint64_t>(index) >= values.size())
        partition_fail(core::FailureCode::InvalidInput);
    return values[static_cast<std::size_t>(index)];
}

[[nodiscard]] inline PgGraphView validate_reference_inputs(const PgReferenceInputs& source) {
    partition_take(source.dictionary.validate_for(source.graph));
    const auto graph = partition_take(source.graph.view());
    const auto partition = partition_take(source.partition.graph());
    if (!partition_same_graph(graph, partition))
        partition_fail(core::FailureCode::StaleGeneration);
    return graph;
}

struct ThicknessPlane {
    const PgPartitionPoint& base;
    const R& base_height;
    std::array<R, 2> gradient;

    [[nodiscard]] R evaluate(ReferenceArithmetic& m, const R& s, const R& t) const {
        return m.add(base_height, m.add(m.mul(gradient[0], m.sub(s, base.s)),
                                        m.mul(gradient[1], m.sub(t, base.t))));
    }
};

[[nodiscard]] inline ThicknessPlane thickness_plane(ReferenceArithmetic& m,
                                                    const PgPartitionThicknessView& thickness,
                                                    core::Index triangle) {
    const auto& ids = reference_at(thickness.triangles, triangle);
    const auto& a = reference_at(thickness.vertices, ids[0]);
    const auto& b = reference_at(thickness.vertices, ids[1]);
    const auto& c = reference_at(thickness.vertices, ids[2]);
    const auto& base = reference_at(thickness.values, ids[0]);
    const auto bs = m.sub(b.s, a.s), bt = m.sub(b.t, a.t);
    const auto cs = m.sub(c.s, a.s), ct = m.sub(c.t, a.t);
    const auto bh = m.sub(reference_at(thickness.values, ids[1]), base);
    const auto ch = m.sub(reference_at(thickness.values, ids[2]), base);
    const auto determinant = m.sub(m.mul(bs, ct), m.mul(bt, cs));
    if (determinant.sign() == 0)
        partition_fail(core::FailureCode::GeometryInvalid);
    return {a,
            base,
            {m.div(m.sub(m.mul(bh, ct), m.mul(bt, ch)), determinant),
             m.div(m.sub(m.mul(bs, ch), m.mul(bh, cs)), determinant)}};
}

struct ReferenceMeasureMap {
    R scale;
    std::array<R, 2> tangent;
    int side = 0; ///< -1 or +1 for horizontal faces, zero for a lateral face.
};

struct ReferenceChartPoint {
    core::Index record;
    core::Index piece;
    const PgPartitionOwnedLayer& layer;
    R s, t, zeta;
    ReferenceMeasureMap map;
};

[[nodiscard]] inline R interpolate(ReferenceArithmetic& m, const R& a, const R& b,
                                   const R& fraction) {
    return m.add(a, m.mul(fraction, m.sub(b, a)));
}

inline void check_fraction(ReferenceArithmetic& m, const R& value, const R& one) {
    if (value.sign() < 0 || m.compare(value, one) > 0)
        partition_fail(core::FailureCode::InvalidInput);
}

template <class Point>
[[nodiscard]] ReferenceChartPoint triangle_point(ReferenceArithmetic& m,
                                                 const PgIntegrationPartition& partition,
                                                 const Point& point) {
    const auto one = m.constant("1");
    if (point.u.sign() < 0 || point.v.sign() < 0 || m.compare(m.add(point.u, point.v), one) > 0)
        partition_fail(core::FailureCode::InvalidInput);
    const auto layers = partition_take(partition.layers());
    const auto triangles = partition_take(partition.triangles());
    const auto vertices = partition_take(partition.vertices());

    const auto locate = [&](core::Index record, const auto& source, auto&& make_depth) {
        const auto& layer = reference_at(layers, source.layer);
        const auto& triangle = reference_at(triangles, source.triangle);
        const auto& a = reference_at(vertices, triangle.vertices[0]);
        const auto& b = reference_at(vertices, triangle.vertices[1]);
        const auto& c = reference_at(vertices, triangle.vertices[2]);
        const auto blend = [&](const R& x, const R& y, const R& z) {
            return m.add(x, m.add(m.mul(point.u, m.sub(y, x)), m.mul(point.v, m.sub(z, x))));
        };
        auto [zeta, map] = make_depth(layer, triangle);
        return ReferenceChartPoint{
            record,          triangle.piece, layer, blend(a.s, b.s, c.s), blend(a.t, b.t, c.t),
            std::move(zeta), std::move(map)};
    };

    if constexpr (std::is_same_v<Point, PgReferenceVolumePoint>) {
        const auto& cell = reference_at(partition_take(partition.volumes()), point.volume);
        check_fraction(m, point.depth, one);
        return locate(point.volume, cell, [&](const auto& layer, const auto& tri) {
            return std::pair{
                interpolate(m, layer.lower, layer.upper, point.depth),
                ReferenceMeasureMap{m.mul(tri.twice_area, m.sub(layer.upper, layer.lower)), {}, 0}};
        });
    } else {
        const auto& face = reference_at(partition_take(partition.horizontal_faces()), point.face);
        return locate(point.face, face, [&](const auto& layer, const auto& tri) {
            const bool lower = face.side == PgPartitionFaceSide::Lower;
            return std::pair{m.copy(lower ? layer.lower : layer.upper),
                             ReferenceMeasureMap{m.copy(tri.twice_area), {}, lower ? -1 : 1}};
        });
    }
}

[[nodiscard]] inline ReferenceChartPoint lateral_point(ReferenceArithmetic& m,
                                                       const PgIntegrationPartition& partition,
                                                       const PgReferenceLateralPoint& point) {
    const auto one = m.constant("1");
    check_fraction(m, point.along, one);
    check_fraction(m, point.depth, one);
    const auto& face = reference_at(partition_take(partition.lateral_faces()), point.face);
    const auto& edge = reference_at(partition_take(partition.edges()), face.edge);
    const auto& layer = reference_at(partition_take(partition.layers()), face.layer);
    const auto vertices = partition_take(partition.vertices());
    const auto owner = face.canonical_direction ? edge.left : edge.right;
    if (!edge.external || !owner)
        partition_fail(core::FailureCode::InvalidInput);
    const auto& a = reference_at(vertices, edge.vertices[face.canonical_direction ? 0 : 1]);
    const auto& b = reference_at(vertices, edge.vertices[face.canonical_direction ? 1 : 0]);
    return {point.face,
            *owner,
            layer,
            interpolate(m, a.s, b.s, point.along),
            interpolate(m, a.t, b.t, point.along),
            interpolate(m, layer.lower, layer.upper, point.depth),
            {m.sub(layer.upper, layer.lower), {m.sub(b.s, a.s), m.sub(b.t, a.t)}, 0}};
}

[[nodiscard]] inline R normal_tube_bound(ReferenceArithmetic& m, const PgGraphView& graph,
                                         const PgIntegrationPartition& partition, const R& length) {
    const auto pieces = partition_take(partition.pieces());
    const auto indices = partition_take(partition.piece_vertex_indices());
    const auto positions = partition_take(partition.vertices());
    const auto thickness = partition_take(partition.thickness());
    m.visit(partition_add(pieces.size(), indices.size()));
    const auto one = m.constant("1");
    const auto two_length = m.mul(m.constant("2"), length);
    const auto determinant_scale = m.div(m.mul(m.square(length), graph.width), m.constant("8"));
    if (pieces.empty())
        partition_fail(core::FailureCode::GeometryInvalid);

    R bound;
    for (const auto& piece : pieces) {
        m.cancel();
        const auto plane = thickness_plane(m, thickness, piece.thickness_triangle);
        auto h_min = R{}, h_max = R{};
        const PgPartitionPoint* max_s = nullptr;
        const PgPartitionPoint* max_t = nullptr;
        for (core::Index i = 0; i < piece.vertex_count; ++i) {
            m.cancel();
            const auto& p = reference_at(positions, reference_at(indices, piece.first_vertex + i));
            auto height = plane.evaluate(m, p.s, p.t);
            if (height.sign() <= 0)
                partition_fail(core::FailureCode::GeometryInvalid);
            if (h_min.sign() == 0 || m.compare(height, h_min) < 0)
                h_min = m.copy(height);
            if (m.compare(height, h_max) > 0)
                h_max = std::move(height);
            if (max_s == nullptr || m.compare(p.s, max_s->s) > 0)
                max_s = &p;
            if (max_t == nullptr || m.compare(p.t, max_t->t) > 0)
                max_t = &p;
        }
        if (max_s == nullptr || max_t == nullptr)
            partition_fail(core::FailureCode::InvalidInput);
        // Each chart coordinate is nonnegative and no larger than this piece's vertex maximum.
        const auto b_squared = physical_hessian_squared_bound(m, graph, max_s->s, max_t->t);
        const auto q = m.mul(m.square(m.div(h_max, two_length)), b_squared);
        if (m.compare(q, one) >= 0)
            partition_fail(core::FailureCode::PrecisionUnresolved);
        auto lower = m.mul(determinant_scale, m.mul(h_min, m.square(m.sub(one, q))));
        if (bound.sign() == 0 || m.compare(lower, bound) < 0)
            bound = std::move(lower);
    }
    return bound;
}

} // namespace uvlas::structure::detail
