#include "global_geometry_private.hpp"

namespace uvlas::aero::detail::global_detail {
ExactInteger multiply(const ExactInteger& left, const ExactInteger& right, Run& r) {
    r.step(multiply_cost);
    if (left == 0 || right == 0)
        return 0;
    const bool negative = (left < 0) != (right < 0);
    ExactInteger a = left < 0 ? -left : left, b = right < 0 ? -right : right;
    const auto total_bits = boost::multiprecision::msb(a) + boost::multiprecision::msb(b) + 2U;
    if (total_bits > 32769U)
        reject(core::FailureCode::PrecisionUnsupported);
    const auto ap = boost::multiprecision::lsb(a), bp = boost::multiprecision::lsb(b);
    a >>= ap;
    b >>= bp;
    const auto na = static_cast<std::size_t>(boost::multiprecision::msb(a) / 32 + 1);
    const auto nb = static_cast<std::size_t>(boost::multiprecision::msb(b) / 32 + 1);
    auto& s = r.scratch;
    s.product_left.fill(0);
    s.product_right.fill(0);
    s.product_output.fill(0);
    boost::multiprecision::export_bits(a, s.product_left.begin(), 32, false);
    boost::multiprecision::export_bits(b, s.product_right.begin(), 32, false);
    for (std::size_t i = 0; i < na; ++i) {
        std::uint64_t carry = 0;
        for (std::size_t j = 0; j < nb; ++j) {
            if (i + j >= s.product_output.size())
                reject(core::FailureCode::PrecisionUnsupported);
            // (2^32-1)^2 + digit + carry <= 2^64-1, including the dense-limb case.
            const auto value = static_cast<std::uint64_t>(s.product_left[i]) * s.product_right[j] +
                               s.product_output[i + j] + carry;
            s.product_output[i + j] = static_cast<std::uint32_t>(value);
            carry = value >> 32U;
        }
        std::size_t k = i + nb;
        while (carry != 0) {
            if (k >= s.product_output.size())
                reject(core::FailureCode::PrecisionUnsupported);
            const auto value = static_cast<std::uint64_t>(s.product_output[k]) + carry;
            s.product_output[k++] = static_cast<std::uint32_t>(value);
            carry = value >> 32U;
        }
    }
    ExactInteger value = 0;
    for (auto k = std::min(na + nb, s.product_output.size()); k > 0; --k) {
        value <<= 32U;
        value += s.product_output[k - 1];
    }
    value <<= ap + bp;
    return negative ? -value : value;
}
Vector subtract(const Vector& a, const Vector& b, Run& r) {
    return {sub(a[0], b[0], r), sub(a[1], b[1], r), sub(a[2], b[2], r)};
}
Vector cross(const Vector& a, const Vector& b, Run& r) {
    return {sub(multiply(a[1], b[2], r), multiply(a[2], b[1], r), r),
            sub(multiply(a[2], b[0], r), multiply(a[0], b[2], r), r),
            sub(multiply(a[0], b[1], r), multiply(a[1], b[0], r), r)};
}
ExactInteger dot(const Vector& a, const Vector& b, Run& r) {
    auto value = multiply(a[0], b[0], r);
    value = add(value, multiply(a[1], b[1], r), r);
    return add(value, multiply(a[2], b[2], r), r);
}
namespace {
bool zero(const Vector& a) {
    return a[0] == 0 && a[1] == 0 && a[2] == 0;
}
std::size_t axis(const Vector& v, Run& r) {
    std::size_t best = 0;
    for (std::size_t k = 1; k < 3; ++k) {
        const ExactInteger a = v[k] < 0 ? -v[k] : v[k];
        const ExactInteger b = v[best] < 0 ? -v[best] : v[best];
        if (compare(a, b, r) > 0)
            best = k;
    }
    return best;
}
void set_vertex(RationalPoint& p, const Vector& x) {
    p.numerator = x;
    p.denominator = 1;
}
bool same_point(const RationalPoint& p, const RationalPoint& q, Run& r) {
    for (std::size_t k = 0; k < 3; ++k)
        if (compare(multiply(p.numerator[k], q.denominator, r),
                    multiply(q.numerator[k], p.denominator, r), r) != 0)
            return false;
    return true;
}
int point_order(const RationalPoint& p, const RationalPoint& q, std::size_t k, Run& r) {
    return compare(multiply(p.numerator[k], q.denominator, r),
                   multiply(q.numerator[k], p.denominator, r), r);
}
bool at_vertex(const RationalPoint& p, const Vector& x, Run& r) {
    for (std::size_t k = 0; k < 3; ++k)
        if (compare(p.numerator[k], multiply(x[k], p.denominator, r), r) != 0)
            return false;
    return true;
}
bool on_segment(const RationalPoint& p, const Vector& a, const Vector& b, Run& r) {
    Vector offset;
    for (std::size_t k = 0; k < 3; ++k)
        offset[k] = sub(p.numerator[k], multiply(a[k], p.denominator, r), r);
    const auto direction = subtract(b, a, r);
    if (!zero(cross(offset, direction, r)))
        return false;
    const auto k = axis(direction, r);
    const auto x = multiply(a[k], p.denominator, r), y = multiply(b[k], p.denominator, r);
    const auto reversed = compare(x, y, r) > 0;
    return compare(p.numerator[k], reversed ? y : x, r) >= 0 &&
           compare(p.numerator[k], reversed ? x : y, r) <= 0;
}
std::size_t cut(const std::array<Vector, 3>& vertices, const std::array<ExactInteger, 3>& h,
                std::span<RationalPoint> output, Run& r) {
    std::size_t count = 0;
    for (std::size_t k = 0; k < 3; ++k)
        if (h[k] == 0) {
            if (count == output.size())
                reject(core::FailureCode::PrecisionUnsupported);
            set_vertex(output[count++], vertices[k]);
        }
    for (std::size_t k = 0; k < 3; ++k) {
        const auto j = (k + 1) % 3;
        if ((h[k] < 0 && h[j] > 0) || (h[k] > 0 && h[j] < 0)) {
            if (count == output.size())
                reject(core::FailureCode::PrecisionUnsupported);
            auto& p = output[count++];
            p.denominator = sub(h[k], h[j], r);
            for (std::size_t d = 0; d < 3; ++d)
                p.numerator[d] =
                    sub(multiply(h[k], vertices[j][d], r), multiply(h[j], vertices[k][d], r), r);
            if (p.denominator < 0) {
                p.denominator = -p.denominator;
                for (auto& x : p.numerator)
                    x = -x;
            }
        }
    }
    if (count > 2)
        reject(core::FailureCode::PrecisionUnsupported);
    return count;
}
ExactInteger cross2(const Vector& a, const Vector& b, std::size_t drop, Run& r) {
    const auto x = (drop + 1) % 3, y = (drop + 2) % 3;
    return sub(multiply(a[x], b[y], r), multiply(a[y], b[x], r), r);
}
bool inside(const Vector& p, const std::array<Vector, 3>& t, std::size_t drop, Run& r) {
    bool positive = false, negative = false;
    for (std::size_t k = 0; k < 3; ++k) {
        const auto h = cross2(subtract(t[(k + 1) % 3], t[k], r), subtract(p, t[k], r), drop, r);
        positive = positive || h > 0;
        negative = negative || h < 0;
    }
    return !(positive && negative);
}
bool within_fraction(const ExactInteger& n, const ExactInteger& d, Run& r) {
    return d > 0 ? compare(n, ExactInteger{0}, r) >= 0 && compare(n, d, r) <= 0
                 : compare(n, ExactInteger{0}, r) <= 0 && compare(n, d, r) >= 0;
}
} // namespace

bool intersects_outside_shared(const GeometryData& g, const GlobalTopologyData& topology,
                               std::uint64_t ia, std::uint64_t ib, Run& r) {
    const auto& fa = g.facets.values()[static_cast<std::size_t>(ia)];
    const auto& fb = g.facets.values()[static_cast<std::size_t>(ib)];
    std::array<Vector, 3> a, b;
    for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t d = 0; d < 3; ++d) {
            a[k][d] = coordinate(g, fa.points[k], d, Configuration::Current, r);
            b[k][d] = coordinate(g, fb.points[k], d, Configuration::Current, r);
        }
    const auto ca = cross(subtract(a[1], a[0], r), subtract(a[2], a[0], r), r);
    const auto cb = cross(subtract(b[1], b[0], r), subtract(b[2], b[0], r), r);
    std::array<ExactInteger, 3> ha, hb;
    for (std::size_t k = 0; k < 3; ++k) {
        hb[k] = dot(ca, subtract(b[k], a[0], r), r);
        ha[k] = dot(cb, subtract(a[k], b[0], r), r);
    }
    const auto outside = [](const auto& h) {
        return (h[0] > 0 && h[1] > 0 && h[2] > 0) || (h[0] < 0 && h[1] < 0 && h[2] < 0);
    };
    if (outside(ha) || outside(hb))
        return false;
    auto& candidates = r.scratch.candidates;
    std::size_t count = 0;
    const auto direction = cross(ca, cb, r);
    if (!zero(direction)) {
        const auto k = axis(direction, r);
        auto& cuts = r.scratch.cuts;
        const auto na = cut(a, ha, {cuts.data(), 2}, r), nb = cut(b, hb, {cuts.data() + 2, 2}, r);
        if (na == 0 || nb == 0)
            return false;
        const auto al = na == 2 && point_order(cuts[0], cuts[1], k, r) > 0 ? 1U : 0U;
        const auto ah = na == 2 ? 1U - al : al;
        const auto bl = nb == 2 && point_order(cuts[2], cuts[3], k, r) > 0 ? 3U : 2U;
        const auto bh = nb == 2 ? 5U - bl : bl;
        const auto lo = point_order(cuts[al], cuts[bl], k, r) >= 0 ? al : bl;
        const auto hi = point_order(cuts[ah], cuts[bh], k, r) <= 0 ? ah : bh;
        if (point_order(cuts[lo], cuts[hi], k, r) > 0)
            return false;
        candidates[count++] = cuts[lo];
        if (!same_point(cuts[lo], cuts[hi], r))
            candidates[count++] = cuts[hi];
    } else {
        if (hb[0] != 0 || hb[1] != 0 || hb[2] != 0)
            return false;
        const auto drop = axis(ca, r);
        for (std::size_t k = 0; k < 3; ++k) {
            if (inside(a[k], b, drop, r))
                set_vertex(candidates[count++], a[k]);
            if (inside(b[k], a, drop, r))
                set_vertex(candidates[count++], b[k]);
        }
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                const auto u = subtract(a[(i + 1) % 3], a[i], r),
                           v = subtract(b[(j + 1) % 3], b[j], r);
                const auto delta = subtract(b[j], a[i], r);
                auto denominator = cross2(u, v, drop, r);
                if (denominator == 0)
                    continue; // Collinear extrema are original contained vertices.
                auto numerator = cross2(delta, v, drop, r), other = cross2(delta, u, drop, r);
                if (!within_fraction(numerator, denominator, r) ||
                    !within_fraction(other, denominator, r))
                    continue;
                auto& p = candidates[count++];
                for (std::size_t k = 0; k < 3; ++k)
                    p.numerator[k] =
                        add(multiply(a[i][k], denominator, r), multiply(u[k], numerator, r), r);
                p.denominator = denominator;
                if (p.denominator < 0) {
                    p.denominator = -p.denominator;
                    for (auto& x : p.numerator)
                        x = -x;
                }
            }
    }
    if (count == 0)
        return false;
    // The coplanar construction emits at most 6 contained vertices + 9 edge
    // crossings. Remove exact duplicates before qualifying its convex hull.
    if (count > 15)
        reject(core::FailureCode::PrecisionUnsupported);
    std::size_t unique = 0;
    for (std::size_t i = 0; i < count; ++i) {
        bool duplicate = false;
        for (std::size_t j = 0; j < unique; ++j)
            if (same_point(candidates[i], candidates[j], r)) {
                duplicate = true;
                break;
            }
        if (!duplicate) {
            if (i != unique) {
                r.step(4 * simple_cost);
                candidates[unique] = candidates[i];
            }
            ++unique;
        }
    }
    count = unique;
    // A shared declared edge is convex; endpoint containment certifies its complete hull.
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            const auto ea = topology.edges.values()[static_cast<std::size_t>(fa.edges[i])].owner;
            const auto eb = topology.edges.values()[static_cast<std::size_t>(fb.edges[j])].owner;
            if (ea != eb)
                continue;
            const auto& edge = g.edges.values()[static_cast<std::size_t>(ea)];
            Vector p, q;
            for (std::size_t d = 0; d < 3; ++d) {
                p[d] = coordinate(g, edge.points[0], d, Configuration::Current, r);
                q[d] = coordinate(g, edge.points[1], d, Configuration::Current, r);
            }
            bool all = true;
            for (std::size_t k = 0; k < count; ++k)
                all = all && on_segment(candidates[k], p, q, r);
            if (all)
                return false;
        }
    // Isolated allowed vertices form a nonconvex set: all extrema must be the SAME point.
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            if (fa.points[i] != fb.points[j])
                continue;
            bool all = true;
            for (std::size_t k = 0; k < count; ++k)
                all = all && at_vertex(candidates[k], a[i], r);
            if (all)
                return false;
        }
    return true;
}
} // namespace uvlas::aero::detail::global_detail
