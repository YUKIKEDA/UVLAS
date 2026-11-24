#include "gmsh_surface_topology.hpp"

#include <algorithm>
#include <limits>
#include <tuple>
#include <utility>

namespace uvlas::io::detail {
namespace {

using Code = core::FailureCode;
using Reason = GmshTopologyFailure;
constexpr auto operation = core::Operation::PrepareModel;
constexpr auto absent = std::numeric_limits<std::size_t>::max();

struct Abort {
    core::Failure failure;
    GmshTopologyDiagnostic diagnostic;
};

[[noreturn]] void reject(GmshTopologyDiagnostic diagnostic, Code code = Code::InvalidInput) {
    throw Abort{{code, operation}, diagnostic};
}

struct Work {
    std::uint64_t limit;
    const std::stop_token& stop;
    std::uint64_t used = 0;

    void charge(std::uint64_t units = 1) {
        if (stop.stop_requested())
            reject({Reason::Limit}, Code::Cancelled);
        if (units > limit - used)
            reject({Reason::Limit}, Code::ResourceExceeded);
        used += units;
    }
};

void add_one(std::size_t& count, std::uint64_t limit, GmshTopologyDiagnostic diagnostic) {
    if (count >= limit || count == std::numeric_limits<std::size_t>::max())
        reject(diagnostic, Code::ResourceExceeded);
    ++count;
}

template <class T>
void admit(std::size_t count, std::uint64_t& remaining, Work& work) {
    if (!std::in_range<core::Index>(count) || count > remaining / sizeof(T))
        reject({Reason::Limit}, Code::ResourceExceeded);
    remaining -= count * sizeof(T);
    work.charge(count);
}

template <class T>
auto allocate(core::ResourcePlan& plan, std::size_t count) {
    auto buffer = plan.allocate<T>(static_cast<core::Index>(count), operation);
    if (!buffer)
        throw Abort{buffer.error(), {Reason::Limit}};
    return std::move(*buffer);
}

const GmshElement& select_face(std::span<const GmshElement> elements, std::uint64_t tag,
                               Work& work) {
    const auto found = std::ranges::lower_bound(
        elements, tag,
        [&work](auto a, auto b) {
            work.charge();
            return a < b;
        },
        &GmshElement::tag);
    if (found == elements.end() || found->tag != tag)
        reject({Reason::UndefinedFace, tag});
    if (found->kind != GmshElementKind::Triangle && found->kind != GmshElementKind::Quadrangle)
        reject({Reason::NonSurfaceElement, tag});
    return *found;
}

std::size_t count_edges(std::span<const GmshElement> elements,
                        std::span<const std::uint64_t> selection, GmshTopologyLimits limits,
                        Work& work) {
    if (selection.size() > limits.faces)
        reject({Reason::Limit, selection[static_cast<std::size_t>(limits.faces)]},
               Code::ResourceExceeded);
    std::size_t count = 0;
    for (const auto tag : selection) {
        work.charge();
        const auto& element = select_face(elements, tag, work);
        const auto ring = element.connectivity();
        for (std::size_t i = 0; i < ring.size(); ++i) {
            work.charge();
            const auto prior = ring.first(i);
            if (std::ranges::find(prior, ring[i]) != prior.end())
                reject({Reason::RepeatedVertex, tag, 0, ring[i]});
            add_one(count, limits.directed_edges, {Reason::Limit, tag, 0, ring[i]});
        }
    }
    return count;
}

struct RingKey {
    std::array<std::uint64_t, 4> vertices{};
    std::size_t face = 0;
};

RingKey canonical_ring(std::span<const std::uint64_t> ring, std::size_t ordinal) noexcept {
    // Each validated ring has distinct nodes, making its smallest tag a unique start.
    const auto start = static_cast<std::size_t>(std::ranges::min_element(ring) - ring.begin());
    auto forward = RingKey{{}, ordinal};
    auto reverse = RingKey{{}, ordinal};
    for (std::size_t i = 0; i < ring.size(); ++i) {
        forward.vertices[i] = ring[(start + i) % ring.size()];
        reverse.vertices[i] = ring[(start + ring.size() - i) % ring.size()];
    }
    return forward.vertices < reverse.vertices ? forward : reverse;
}

void build_faces(std::span<const GmshElement> elements, std::span<const std::uint64_t> selection,
                 std::span<GmshTopologyFace> faces, std::span<GmshTopologyEdge> edges,
                 std::span<RingKey> keys, Work& work) {
    std::size_t start = 0;
    for (std::size_t i = 0; i < faces.size(); ++i) {
        work.charge();
        const auto& element = select_face(elements, selection[i], work);
        const auto ring = element.connectivity();
        faces[i] = {element, start, 0};
        keys[i] = canonical_ring(ring, i);
        for (std::size_t j = 0; j < ring.size(); ++j) {
            work.charge();
            const auto next = (j + 1) % ring.size();
            edges[start + j] = {ring[j], ring[next], i, start + next, {}};
        }
        start += ring.size();
    }
    std::ranges::sort(keys, [&work](const RingKey& a, const RingKey& b) {
        work.charge();
        return std::tie(a.vertices, a.face) < std::tie(b.vertices, b.face);
    });
    for (std::size_t i = 1; i < keys.size(); ++i) {
        work.charge();
        if (keys[i - 1].vertices != keys[i].vertices)
            continue;
        const auto a = faces[keys[i - 1].face].element.tag;
        const auto b = faces[keys[i].face].element.tag;
        reject({a == b ? Reason::RepeatedSelection : Reason::DuplicateFace, b, a});
    }
}

struct Parent {
    std::size_t root = 0, rank = 0;
};

class Sets {
public:
    Sets(std::span<Parent> parents, Work& work) : parents_(parents), work_(work) {
        for (std::size_t i = 0; i < parents_.size(); ++i) {
            work_.charge();
            parents_[i].root = i;
        }
    }

    std::size_t find(std::size_t i) {
        work_.charge();
        while (parents_[i].root != i) {
            work_.charge();
            parents_[i].root = parents_[parents_[i].root].root;
            i = parents_[i].root;
        }
        return i;
    }

    void join(std::size_t a, std::size_t b) {
        a = find(a);
        b = find(b);
        if (a == b)
            return;
        if (parents_[a].rank < parents_[b].rank)
            std::swap(a, b);
        parents_[b].root = a;
        if (parents_[a].rank == parents_[b].rank)
            ++parents_[a].rank;
    }

private:
    std::span<Parent> parents_;
    Work& work_;
};

auto edge_key(const GmshTopologyEdge& edge) noexcept {
    return std::pair{std::min(edge.from, edge.to), std::max(edge.from, edge.to)};
}

std::size_t match_edges(std::span<GmshTopologyEdge> edges, std::span<const GmshTopologyFace> faces,
                        std::span<std::size_t> order, Sets& components, Sets& fans,
                        GmshTopologyLimits limits, Work& work) {
    for (std::size_t i = 0; i < order.size(); ++i) {
        work.charge();
        order[i] = i;
    }
    std::ranges::sort(order, [&](std::size_t a, std::size_t b) {
        work.charge();
        return std::pair{edge_key(edges[a]), a} < std::pair{edge_key(edges[b]), b};
    });

    std::size_t boundary = 0;
    for (std::size_t i = 0; i < order.size();) {
        work.charge();
        const auto first = order[i++];
        auto& a = edges[first];
        const auto key = edge_key(a);
        if (i == order.size() || edge_key(edges[order[i]]) != key) {
            add_one(boundary, limits.boundary_edges,
                    {Reason::Limit, faces[a.face].element.tag, 0, a.from, a.to});
            continue;
        }

        const auto second = order[i++];
        auto& b = edges[second];
        if (i != order.size() && edge_key(edges[order[i]]) == key)
            reject({Reason::NonManifoldEdge, faces[a.face].element.tag,
                    faces[edges[order[i]].face].element.tag, key.first, key.second});
        if (a.from == b.from)
            reject({Reason::InconsistentOrientation, faces[a.face].element.tag,
                    faces[b.face].element.tag, key.first, key.second});

        a.opposite = second;
        b.opposite = first;
        components.join(a.face, b.face);
        // A corner is represented by its outgoing edge. Match both ends of this edge.
        fans.join(first, b.next);
        fans.join(a.next, second);
    }
    return boundary;
}

void check_vertex_fans(std::span<const GmshTopologyEdge> edges,
                       std::span<const GmshTopologyFace> faces, std::span<std::size_t> order,
                       Sets& fans, Work& work) {
    std::ranges::sort(order, [&](std::size_t a, std::size_t b) {
        work.charge();
        return std::pair{edges[a].from, a} < std::pair{edges[b].from, b};
    });
    for (std::size_t i = 1; i < order.size(); ++i) {
        work.charge();
        const auto a = order[i - 1];
        const auto b = order[i];
        if (edges[a].from != edges[b].from)
            continue;
        const auto first = fans.find(a);
        const auto second = fans.find(b);
        if (first != second)
            reject({Reason::DisconnectedVertexFan, faces[edges[a].face].element.tag,
                    faces[edges[b].face].element.tag, edges[a].from});
    }
}

std::size_t label_components(std::span<GmshTopologyFace> faces, std::span<std::size_t> labels,
                             Sets& sets, GmshTopologyLimits limits, Work& work) {
    for (auto& label : labels) {
        work.charge();
        label = absent;
    }
    std::size_t count = 0;
    for (std::size_t i = 0; i < faces.size(); ++i) {
        const auto root = sets.find(i);
        if (labels[root] == absent) {
            labels[root] = count;
            add_one(count, limits.components, {Reason::Limit, faces[i].element.tag});
        }
        faces[i].component = labels[root];
    }
    return count;
}

void record_components(std::span<const GmshTopologyFace> faces,
                       std::span<const GmshTopologyEdge> edges, std::span<std::size_t> boundary,
                       std::span<GmshTopologyComponent> components, Work& work) {
    for (std::size_t i = 0; i < faces.size(); ++i) {
        work.charge();
        auto& row = components[faces[i].component];
        if (row.face_count++ == 0)
            row.first_face = i;
    }
    std::size_t start = 0;
    for (std::size_t i = 0; i < edges.size(); ++i) {
        work.charge();
        const auto& edge = edges[i];
        if (!edge.opposite) {
            boundary[start++] = i;
            ++components[faces[edge.face].component].boundary_count;
        }
    }
}

} // namespace

GmshSurfaceTopology::GmshSurfaceTopology(GmshTopologySource source,
                                         core::PlannedBuffer<GmshTopologyFace>&& faces,
                                         core::PlannedBuffer<GmshTopologyEdge>&& edges,
                                         core::PlannedBuffer<std::size_t>&& boundary,
                                         core::PlannedBuffer<GmshTopologyComponent>&& components,
                                         std::uint64_t work) noexcept
    : source_(source),
      faces_(std::move(faces)),
      edges_(std::move(edges)),
      boundary_(std::move(boundary)),
      components_(std::move(components)),
      work_(work) {}

GmshSurfaceTopology::GmshSurfaceTopology(GmshSurfaceTopology&& other) noexcept
    : source_(std::exchange(other.source_, std::nullopt)),
      faces_(std::move(other.faces_)),
      edges_(std::move(other.edges_)),
      boundary_(std::move(other.boundary_)),
      components_(std::move(other.components_)),
      work_(other.work_) {}

std::optional<GmshTopologyView> GmshSurfaceTopology::view() const& noexcept {
    if (!source_)
        return std::nullopt;
    return GmshTopologyView{*source_,           faces_.values(),      edges_.values(),
                            boundary_.values(), components_.values(), work_};
}

core::Result<GmshSurfaceTopology>
GmshSurfaceTopology::prepare(const GmshMesh& mesh, std::span<const std::uint64_t> selected_tags,
                             GmshTopologyLimits limits, core::ResourcePlan& plan,
                             const std::stop_token& stop,
                             GmshTopologyDiagnostic* diagnostic) noexcept {
    if (diagnostic)
        *diagnostic = {};
    return core::guard_entry(
        operation, Code::InvalidInput, {}, [&]() -> core::Result<GmshSurfaceTopology> {
            try {
                Work work{limits.work, stop};
                work.charge(0);
                const auto input = mesh.view();
                if (!input)
                    reject({Reason::Input});
                const auto face_count = selected_tags.size();
                const auto edge_count = count_edges(input->elements, selected_tags, limits, work);
                auto remaining = plan.remaining_bytes();
                admit<GmshTopologyFace>(face_count, remaining, work);
                admit<GmshTopologyEdge>(edge_count, remaining, work);
                admit<RingKey>(face_count, remaining, work);
                admit<std::size_t>(edge_count, remaining, work);
                admit<Parent>(face_count, remaining, work);
                admit<Parent>(edge_count, remaining, work);

                auto faces = allocate<GmshTopologyFace>(plan, face_count);
                auto edges = allocate<GmshTopologyEdge>(plan, edge_count);
                auto keys = allocate<RingKey>(plan, face_count);
                auto index = allocate<std::size_t>(plan, edge_count);
                auto face_parents = allocate<Parent>(plan, face_count);
                auto fan_parents = allocate<Parent>(plan, edge_count);
                build_faces(input->elements, selected_tags, faces.values(), edges.values(),
                            keys.values(), work);
                Sets components{face_parents.values(), work};
                Sets fans{fan_parents.values(), work};
                const auto boundary_count = match_edges(
                    edges.values(), faces.values(), index.values(), components, fans, limits, work);
                check_vertex_fans(edges.values(), faces.values(), index.values(), fans, work);
                const auto component_count = label_components(
                    faces.values(), index.values().first(face_count), components, limits, work);
                remaining = plan.remaining_bytes();
                admit<std::size_t>(boundary_count, remaining, work);
                admit<GmshTopologyComponent>(component_count, remaining, work);
                auto boundary = allocate<std::size_t>(plan, boundary_count);
                auto records = allocate<GmshTopologyComponent>(plan, component_count);
                record_components(faces.values(), edges.values(), boundary.values(),
                                  records.values(), work);
                work.charge(0);
                return GmshSurfaceTopology{{input->location, input->digest},
                                           std::move(faces),
                                           std::move(edges),
                                           std::move(boundary),
                                           std::move(records),
                                           work.used};
            } catch (const Abort& error) {
                if (diagnostic && error.failure.code != Code::Cancelled)
                    *diagnostic = error.diagnostic;
                return std::unexpected(error.failure);
            }
        });
}

} // namespace uvlas::io::detail
