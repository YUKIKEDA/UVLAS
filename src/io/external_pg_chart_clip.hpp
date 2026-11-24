#pragma once

#include "../structure/pg_reference_arithmetic.hpp"

#include <array>
#include <span>
#include <utility>

namespace uvlas::io::detail {

// The plane clips a positively oriented triangle against a positively oriented
// convex triangle or quadrilateral. Its intersection is still one convex ring.
namespace external_pg_chart {

using Arithmetic = structure::detail::ReferenceArithmetic;
using Point = structure::PgPartitionPoint;
using Rational = structure::ExactRational;
using Triangle = std::array<const Point*, 3>;

[[nodiscard]] inline Rational cross(Arithmetic& math, const Point& a, const Point& b) {
    return math.sub(math.mul(a.s, b.t), math.mul(a.t, b.s));
}

[[nodiscard]] inline Rational orient(Arithmetic& math, const Point& a, const Point& b,
                                     const Point& c) {
    const Point first{math.sub(b.s, a.s), math.sub(b.t, a.t)};
    const Point second{math.sub(c.s, a.s), math.sub(c.t, a.t)};
    return cross(math, first, second);
}

[[nodiscard]] inline Point crossing(Arithmetic& math, const Point& p, const Point& q,
                                    const Rational& from, const Rational& to) {
    const auto fraction = math.div(from, math.sub(from, to));
    return {math.add(p.s, math.mul(fraction, math.sub(q.s, p.s))),
            math.add(p.t, math.mul(fraction, math.sub(q.t, p.t)))};
}

struct Ring {
    // A convex triangle clipped by four half-planes has at most seven vertices.
    std::array<Point, 7> storage;
    std::size_t count{};

    [[nodiscard]] std::span<const Point> vertices() const& noexcept {
        return std::span{storage}.first(count);
    }
    [[nodiscard]] std::span<const Point> vertices() const&& = delete;

    void append(Point point) {
        if (count == storage.size()) {
            structure::detail::partition_fail(core::FailureCode::GeometryInvalid);
        }
        storage[count++] = std::move(point);
    }
};

[[nodiscard]] inline Ring clip(Arithmetic& math, const Triangle& triangle,
                               std::span<const Point* const> exterior) {
    Ring ring;
    for (const auto* point : triangle) {
        ring.append({math.copy(point->s), math.copy(point->t)});
    }
    for (std::size_t edge = 0; edge < exterior.size(); ++edge) {
        const auto vertices = ring.vertices();
        if (vertices.size() < 3) {
            return {};
        }
        const auto& a = *exterior[edge];
        const auto& b = *exterior[(edge + 1) % exterior.size()];
        Ring clipped;
        auto previous = orient(math, a, b, vertices.back());
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            auto current = orient(math, a, b, vertices[i]);
            if (previous.sign() * current.sign() < 0) {
                const auto prev = (i + vertices.size() - 1) % vertices.size();
                clipped.append(crossing(math, vertices[prev], vertices[i], previous, current));
            }
            if (current.sign() >= 0) {
                clipped.append({math.copy(vertices[i].s), math.copy(vertices[i].t)});
            }
            previous = std::move(current);
        }
        ring = std::move(clipped);
    }
    return ring;
}

[[nodiscard]] inline Rational area2(Arithmetic& math, const Ring& ring) {
    Rational sum;
    const auto vertices = ring.vertices();
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        sum = math.add(sum, cross(math, vertices[i], vertices[(i + 1) % vertices.size()]));
    }
    return sum;
}

} // namespace external_pg_chart
} // namespace uvlas::io::detail
