#include "global_geometry_private.hpp"

namespace uvlas::aero::detail::global_detail {
namespace {
std::uint64_t log_height(std::uint64_t n) {
    std::uint64_t h = 0;
    for (auto v = n > 0 ? n - 1 : 0; v != 0; v >>= 1U)
        ++h;
    return h;
}
void leaf_box(const GeometryData& g, Node& node, std::uint64_t ordinal, Configuration config,
              Run& r) {
    const auto& face = g.facets.values()[static_cast<std::size_t>(ordinal)];
    node.count = 1;
    node.facet = ordinal;
    for (std::size_t k = 0; k < 3; ++k) {
        node.extrema[2 * k] = face.points[0];
        node.extrema[2 * k + 1] = face.points[0];
        for (std::size_t j = 1; j < 3; ++j) {
            auto x = coordinate(g, face.points[j], k, config, r);
            auto lo = coordinate(g, node.extrema[2 * k], k, config, r);
            auto hi = coordinate(g, node.extrema[2 * k + 1], k, config, r);
            if (compare(x, lo, r) < 0)
                node.extrema[2 * k] = face.points[j];
            if (compare(x, hi, r) > 0)
                node.extrema[2 * k + 1] = face.points[j];
        }
    }
}
void union_box(const GeometryData& g, Node& output, const Node& a, const Node& b,
               Configuration config, Run& r) {
    for (std::size_t k = 0; k < 3; ++k) {
        for (std::size_t side = 0; side < 2; ++side) {
            const auto index = 2 * k + side;
            auto x = coordinate(g, a.extrema[index], k, config, r),
                 y = coordinate(g, b.extrema[index], k, config, r);
            const auto order = compare(x, y, r);
            output.extrema[index] =
                (side == 0 ? order <= 0 : order >= 0) ? a.extrema[index] : b.extrema[index];
        }
    }
}
template <class Less>
void heap_sort(std::span<std::uint64_t> values, Less less, Run& r) {
    const auto sift = [&](std::size_t root, std::size_t end) {
        while (root < end / 2) {
            auto child = 2 * root + 1;
            if (child + 1 < end && less(values[child], values[child + 1]))
                ++child;
            if (!less(values[root], values[child]))
                break;
            r.step();
            std::swap(values[root], values[child]);
            root = child;
        }
    };
    for (auto n = values.size() / 2; n > 0; --n)
        sift(n - 1, values.size());
    for (auto end = values.size(); end > 1; --end) {
        r.step();
        std::swap(values[0], values[end - 1]);
        sift(0, end - 1);
    }
}
struct Builder {
    GlobalTopologyData& data;
    const GeometryData& geometry;
    Run& run;
    std::uint64_t next;
    std::uint64_t build(std::span<std::uint64_t> range) {
        if (range.size() == 1)
            return range[0];
        Node bounds = data.nodes.values()[static_cast<std::size_t>(range[0])];
        for (std::size_t i = 1; i < range.size(); ++i)
            union_box(geometry, bounds, bounds,
                      data.nodes.values()[static_cast<std::size_t>(range[i])],
                      Configuration::Reference, run);
        std::size_t axis = 0;
        ExactInteger longest = -1;
        for (std::size_t k = 0; k < 3; ++k) {
            const auto lo =
                coordinate(geometry, bounds.extrema[2 * k], k, Configuration::Reference, run);
            const auto hi =
                coordinate(geometry, bounds.extrema[2 * k + 1], k, Configuration::Reference, run);
            auto length = sub(hi, lo, run);
            if (compare(length, longest, run) > 0) {
                longest = std::move(length);
                axis = k;
            }
        }
        heap_sort(
            range,
            [&](auto ia, auto ib) {
                const auto& a = data.nodes.values()[static_cast<std::size_t>(ia)];
                const auto& b = data.nodes.values()[static_cast<std::size_t>(ib)];
                const auto al =
                    coordinate(geometry, a.extrema[2 * axis], axis, Configuration::Reference, run);
                const auto ah = coordinate(geometry, a.extrema[2 * axis + 1], axis,
                                           Configuration::Reference, run);
                const auto bl =
                    coordinate(geometry, b.extrema[2 * axis], axis, Configuration::Reference, run);
                const auto bh = coordinate(geometry, b.extrema[2 * axis + 1], axis,
                                           Configuration::Reference, run);
                const auto order = compare(add(al, ah, run), add(bl, bh, run), run);
                if (order != 0)
                    return order < 0;
                run.step();
                return geometry.facets.values()[static_cast<std::size_t>(ia)].id <
                       geometry.facets.values()[static_cast<std::size_t>(ib)].id;
            },
            run);
        const auto mid = range.size() / 2;
        const auto left = build(range.first(mid)), right = build(range.subspan(mid));
        const auto index = next++;
        auto& node = data.nodes.values()[static_cast<std::size_t>(index)];
        node.left = left;
        node.right = right;
        node.count = range.size();
        node.extrema = bounds.extrema;
        return index;
    }
};
bool box_overlap(const GeometryData& g, const Node& a, const Node& b, Run& r) {
    for (std::size_t k = 0; k < 3; ++k) {
        const auto al = coordinate(g, a.extrema[2 * k], k, Configuration::Current, r);
        const auto bh = coordinate(g, b.extrema[2 * k + 1], k, Configuration::Current, r);
        if (compare(al, bh, r) > 0)
            return false;
        const auto bl = coordinate(g, b.extrema[2 * k], k, Configuration::Current, r);
        const auto ah = coordinate(g, a.extrema[2 * k + 1], k, Configuration::Current, r);
        if (compare(bl, ah, r) > 0)
            return false;
    }
    return true;
}
} // namespace
std::uint64_t tree_shape(std::uint64_t f) {
    return f == 0 ? 0 : plus(f, f - 1);
}
std::uint64_t pair_stack(std::uint64_t f) {
    return plus(1, times(4, log_height(f)));
}
void build_tree(GlobalTopologyData& data, Run& r) {
    const auto& g = *GlobalEngine::source(*data.geometry);
    const auto f = g.facets.size();
    for (std::size_t i = 0; i < f; ++i) {
        data.order.values()[i] = i;
        leaf_box(g, data.nodes.values()[i], i, Configuration::Reference, r);
    }
    if (f > 0) {
        Builder builder{data, g, r, f};
        const auto root = builder.build(data.order.values());
        if (root != data.nodes.size() - 1 || builder.next != data.nodes.size())
            reject(core::FailureCode::PrecisionUnsupported);
    }
}
void refit_tree(const GeometryData& g, std::span<Node> nodes, Run& r) {
    for (std::size_t i = 0; i < g.facets.size(); ++i)
        leaf_box(g, nodes[i], i, Configuration::Current, r);
    for (std::size_t i = g.facets.size(); i < nodes.size(); ++i)
        union_box(g, nodes[i], nodes[static_cast<std::size_t>(nodes[i].left)],
                  nodes[static_cast<std::size_t>(nodes[i].right)], Configuration::Current, r);
}
TraversalResult traverse(const GeometryData& g, const GlobalTopologyData& topology,
                         std::span<const Node> nodes, std::span<NodePair> stack, bool narrow,
                         Run& r) {
    if (nodes.empty())
        return {};
    const auto begin = r.used;
    std::uint64_t narrow_steps = 0;
    std::size_t used = 0;
    std::uint64_t count = 0;
    const auto push = [&](std::uint64_t a, std::uint64_t b) {
        r.step();
        if (used == stack.size())
            reject(core::FailureCode::ResourceExceeded);
        stack[used++] = {a, b};
    };
    push(nodes.size() - 1, nodes.size() - 1);
    std::optional<std::array<std::uint64_t, 2>> invalid;
    while (used != 0) {
        r.step(2);
        const auto pair = stack[--used];
        const auto& a = nodes[static_cast<std::size_t>(pair.first)];
        const auto& b = nodes[static_cast<std::size_t>(pair.second)];
        if (pair.first == pair.second && a.count == 1)
            continue;
        if (!box_overlap(g, a, b, r))
            continue;
        if (pair.first == pair.second) {
            push(a.right, a.right);
            push(a.left, a.right);
            push(a.left, a.left);
        } else if (a.count == 1 && b.count == 1) {
            count = plus(count, 1);
            const auto before_narrow = r.used;
            const auto outside =
                narrow && intersects_outside_shared(g, topology, a.facet, b.facet, r);
            narrow_steps = plus(narrow_steps, r.used - before_narrow);
            if (outside) {
                auto first = g.facets.values()[static_cast<std::size_t>(a.facet)].id;
                auto second = g.facets.values()[static_cast<std::size_t>(b.facet)].id;
                if (first > second)
                    std::swap(first, second);
                const std::array candidate{first, second};
                if (!invalid || candidate < *invalid)
                    invalid = candidate;
            }
        } else if (b.count == 1 ||
                   (a.count != 1 &&
                    (a.count > b.count || (a.count == b.count && pair.first < pair.second)))) {
            push(a.right, pair.second);
            push(a.left, pair.second);
        } else {
            push(pair.first, b.right);
            push(pair.first, b.left);
        }
    }
    if (invalid)
        reject(core::FailureCode::GeometryInvalid, FacetId{(*invalid)[0]}, FacetId{(*invalid)[1]});
    return {count, r.used - begin - narrow_steps};
}
} // namespace uvlas::aero::detail::global_detail
