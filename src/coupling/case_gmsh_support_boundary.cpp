#include "case_gmsh_support_boundary.hpp"

#include "io/external_pg_chart_clip.hpp"
#include "io/gmsh_entity_physical.hpp"
#include "structure/pg_reference_domain.hpp"
#include "structure/pg_trace_polynomial.hpp"

#include <algorithm>

namespace uvlas::coupling::detail {
namespace {
namespace in = io::detail;
namespace ref = structure::detail;
namespace s = io::schema;
using ref::partition_fail;
using ref::partition_index;
using ref::partition_take;
using Code = core::FailureCode;
using Reason = CaseGmshSupportBoundaryFailure;
using Point = structure::PgPartitionPoint;
using Rational = structure::ExactRational;

struct Query {
    ref::ReferenceArithmetic math;
    std::optional<CaseGmshSupportBoundaryDiagnostic>* diagnostic{};
    CaseGmshSupportBoundaryDiagnostic detail{};

    [[noreturn]] void reject(Reason reason, Code code = Code::InvalidInput) {
        if (diagnostic) {
            detail.reason = reason;
            *diagnostic = detail;
        }
        partition_fail(code);
    }
};

[[nodiscard]] bool same_parent(const io::SchemaIdentity& a, const io::SchemaIdentity& b) {
    return a.source == b.source && a.contract_revision == b.contract_revision &&
           a.representation_revision == b.representation_revision;
}

[[nodiscard]] s::Region support_region(Query& q, const s::Model& model, in::MeshSupportId support) {
    for (std::size_t i = 0; i < model.supports.size(); ++i) {
        q.math.visit(1);
        const auto declaration = partition_take(model.supports.at(partition_index(i))).record();
        const auto id =
            std::visit([](const auto& value) { return value.id.unsigned_value(); }, declaration);
        if (id != support.value()) {
            continue;
        }
        const auto* trace = std::get_if<s::SupportTraceConstraint>(&declaration);
        if (trace == nullptr || trace->components.size() != 3) {
            q.reject(Reason::Support, Code::Unsupported);
        }
        const auto region = trace->region.record();
        const auto target = region.target.record();
        if (target.kind.selected() != decltype(target.kind.selected())::ItemShell ||
            region.kind.selected() != decltype(region.kind.selected())::ItemEdge ||
            region.side.selected() != decltype(region.side.selected())::ItemInterior) {
            q.reject(Reason::Support, Code::Unsupported);
        }
        return region;
    }
    q.reject(Reason::Support);
}

[[nodiscard]] const CaseGmshPgSet::Shell&
prepared_shell(Query& q, const CaseGmshPgSet& set, const in::ExternalMeshSidecarBinding& bound,
               in::MeshShellId shell, in::MeshGeometryId pg) {
    for (std::size_t i = 0; i < set.size(); ++i) {
        q.math.visit(1);
        const auto& chosen = partition_take(set.at(partition_index(i))).get();
        const auto& match = chosen.match();
        if (match.shell == shell) {
            if (match.pg != pg || match.sidecar.source != bound.sidecar().source().source ||
                !same_parent(match.sidecar.parent, bound.sidecar().source().parent) ||
                match.sidecar_digest != bound.sidecar().digest()) {
                q.reject(Reason::Binding, Code::StaleGeneration);
            }
            return chosen;
        }
    }
    q.reject(Reason::Support);
}

struct Segment {
    const Point* a{};
    const Point* b{};
};

struct ChartBoundary {
    std::span<const structure::PgPartitionEdge> edges;
    std::span<const structure::PgPartitionSource> sources;
    std::span<const core::Index> source_indices;
    std::span<const Point> vertices;
    core::Index name{};
};

[[nodiscard]] bool belongs(Query& q, const ChartBoundary& chart,
                           const structure::PgPartitionEdge& edge) {
    if (!edge.external) {
        return false;
    }
    for (const auto index :
         chart.source_indices.subspan(static_cast<std::size_t>(edge.first_source),
                                      static_cast<std::size_t>(edge.source_count))) {
        q.math.visit(1);
        const auto& descriptor = ref::reference_at(chart.sources, index).descriptor;
        if (const auto* named = std::get_if<structure::PgPartitionNamedSource>(&descriptor);
            named && named->boundary == chart.name) {
            return true;
        }
    }
    return false;
}

[[nodiscard]] Segment segment(const ChartBoundary& chart, const structure::PgPartitionEdge& edge) {
    return {&ref::reference_at(chart.vertices, edge.vertices[0]),
            &ref::reference_at(chart.vertices, edge.vertices[1])};
}

[[nodiscard]] Rational length(Query& q, Segment line) {
    const auto absolute_difference = [&](const Rational& a, const Rational& b) {
        return q.math.compare(a, b) < 0 ? q.math.sub(b, a) : q.math.sub(a, b);
    };
    return q.math.add(absolute_difference(line.a->s, line.b->s),
                      absolute_difference(line.a->t, line.b->t));
}

[[nodiscard]] Rational overlap(Query& q, Segment a, Segment b) {
    if (in::external_pg_chart::orient(q.math, *a.a, *a.b, *b.a).sign() != 0 ||
        in::external_pg_chart::orient(q.math, *a.a, *a.b, *b.b).sign() != 0) {
        return {};
    }
    const auto less = [&](const Rational& x, const Rational& y) {
        return q.math.compare(x, y) < 0;
    };
    const auto along = [&](auto component) -> Rational {
        const auto [a0, a1] = std::minmax(a.a->*component, a.b->*component, less);
        const auto [b0, b1] = std::minmax(b.a->*component, b.b->*component, less);
        const auto& lo = std::max(a0, b0, less);
        const auto& hi = std::min(a1, b1, less);
        return less(lo, hi) ? q.math.sub(hi, lo) : Rational{};
    };
    return q.math.add(along(&Point::s), along(&Point::t));
}

using Position = std::array<Rational, 3>;

struct Incident {
    Segment chart;
    const Position* first{};
    const Position* second{};
    std::uint64_t face{};
    bool exact{};
};

[[nodiscard]] bool exact_face(const in::ExternalPgOwnedFace&) noexcept {
    return true;
}

[[nodiscard]] bool exact_face(const in::GmshPgFaceBound& face) noexcept {
    return face.distance_upper_metres.sign() == 0;
}

[[nodiscard]] Incident incident(Query& q, const auto faces, std::span<const std::uint64_t> nodes) {
    std::optional<Incident> found;
    for (const auto& face : faces) {
        const auto ring = face.ring();
        for (std::size_t i = 0; i < ring.size(); ++i) {
            q.math.visit(1);
            const auto& a = ring[i];
            const auto& b = ring[(i + 1) % ring.size()];
            if ((a.node == nodes[0] && b.node == nodes[1]) ||
                (a.node == nodes[1] && b.node == nodes[0])) {
                if (found) {
                    q.reject(Reason::Incidence);
                }
                const auto& first = a.node == nodes[0] ? a : b;
                const auto& second = a.node == nodes[0] ? b : a;
                found = Incident{{&first.chart, &second.chart},
                                 &first.position_metres,
                                 &second.position_metres,
                                 face.element,
                                 exact_face(face)};
            }
        }
    }
    if (!found) {
        q.reject(Reason::Incidence);
    }
    return *found;
}

void compare_edge_position(Query& q, const Incident& edge, const structure::PgGraphView& graph,
                           const structure::PgReferenceFrame& frame) {
    auto& math = q.math;
    const auto s = ref::trace_affine(math, edge.chart.a->s, edge.chart.b->s);
    const auto t = ref::trace_affine(math, edge.chart.a->t, edge.chart.b->t);
    // PG attests to a graph of total degree at most four. The two power
    // tables retain at most 30 populated coefficients together, and the
    // remaining degree-four and affine temporaries fit the reference envelope.
    const auto s_powers = ref::trace_powers(math, s, 4);
    const auto t_powers = ref::trace_powers(math, t, 4);
    const auto height = ref::trace_restrict(math, graph.height, s_powers, t_powers);
    auto x = ref::trace_affine(math, graph.x0, graph.x0);
    ref::trace_add_scaled(math, x, s, graph.width);
    ref::trace_add_scaled(math, x, t, graph.skew);

    for (std::size_t axis = 0; axis < 3; ++axis) {
        math.visit(1);
        const auto original = ref::trace_affine(math, (*edge.first)[axis], (*edge.second)[axis]);
        auto world = ref::trace_affine(math, frame.origin_metres[axis], frame.origin_metres[axis]);
        ref::trace_add_scaled(math, world, x,
                              math.mul(frame.length_metres, frame.rotation[3 * axis]));
        ref::trace_add_scaled(math, world, t,
                              math.mul(frame.length_metres, frame.rotation[3 * axis + 1]));
        ref::trace_add_scaled(math, world, height,
                              math.mul(frame.length_metres, frame.rotation[3 * axis + 2]));
        if (!ref::trace_same(math, world, original)) {
            q.detail.world_component = static_cast<std::uint8_t>(axis);
            q.reject(Reason::Geometry);
        }
    }
}

[[nodiscard]] core::PlannedBuffer<CaseGmshSupportLine>
cover(Query& q, const in::GmshMeshView& source, in::GmshKey group, const auto faces,
      const ChartBoundary& chart, const structure::PgGraphView& graph,
      const structure::PgReferenceFrame& frame, std::uint64_t maximum_lines) {
    const auto member = [&](const in::GmshElement& element) {
        q.math.visit(1);
        return in::gmsh_element_in_physical_group(source, element, group, [&] { q.math.visit(1); });
    };
    std::uint64_t count = 0;
    for (const auto& element : source.elements) {
        if (member(element)) {
            if (element.kind != in::GmshElementKind::Line) {
                q.detail.line = element.tag;
                q.reject(Reason::Support, Code::Unsupported);
            }
            if (count == maximum_lines) {
                q.detail.line = element.tag;
                q.reject(Reason::Support, Code::ResourceExceeded);
            }
            ++count;
        }
    }
    auto result = partition_take(q.math.plan.allocate<CaseGmshSupportLine>(
        partition_index(count), ref::partition_operation));
    // Distinct certified free edges have disjoint open chart intervals. The
    // L1 variation is additive on each collinear part of the named boundary.
    Rational total;
    std::size_t written = 0;
    for (const auto& element : source.elements) {
        if (!member(element)) {
            continue;
        }
        q.detail.line = element.tag;
        q.detail.face.reset();
        const auto nodes = element.connectivity();
        for (const auto& previous : result.values().first(written)) {
            q.math.visit(1);
            const auto& pair = previous.nodes;
            if ((pair[0] == nodes[0] && pair[1] == nodes[1]) ||
                (pair[0] == nodes[1] && pair[1] == nodes[0])) {
                q.reject(Reason::DuplicateLine);
            }
        }
        const auto edge = incident(q, faces, nodes);
        q.detail.face = edge.face;
        if (!edge.exact) {
            compare_edge_position(q, edge, graph, frame);
        }
        Rational covered;
        for (const auto& atomic : chart.edges) {
            q.math.visit(1);
            if (belongs(q, chart, atomic)) {
                covered = q.math.add(covered, overlap(q, edge.chart, segment(chart, atomic)));
            }
        }
        if (covered.sign() == 0 || q.math.compare(covered, length(q, edge.chart)) != 0) {
            q.reject(Reason::OutsideBoundary);
        }
        total = q.math.add(total, covered);
        const auto evidence = edge.exact ? CaseGmshSupportPositionEvidence::FaceEquality
                                         : CaseGmshSupportPositionEvidence::EdgePolynomial;
        result.values()[written++] = {element.tag, {nodes[0], nodes[1]}, edge.face, evidence};
    }

    q.detail.line.reset();
    q.detail.face.reset();
    Rational expected;
    for (const auto& edge : chart.edges) {
        q.math.visit(1);
        if (belongs(q, chart, edge)) {
            expected = q.math.add(expected, length(q, segment(chart, edge)));
        }
    }
    if (expected.sign() == 0 || q.math.compare(total, expected) != 0) {
        q.reject(Reason::IncompleteBoundary);
    }
    return result;
}
} // namespace

core::Result<CaseGmshSupportBoundary> certify_case_gmsh_support_boundary(
    const io::SiCase& si, const in::ExternalMeshSidecarBinding& bound, const CaseGmshPgSet& set,
    in::MeshSupportId requested, in::GmshKey physical, const CaseGmshSupportBoundaryLimits& limits,
    core::ResourcePlan& plan,
    std::optional<CaseGmshSupportBoundaryDiagnostic>* diagnostic) noexcept {
    if (diagnostic) {
        diagnostic->reset();
    }
    return ref::partition_entry([&]() -> core::Result<CaseGmshSupportBoundary> {
        Query q{{limits.arithmetic, plan}, diagnostic, {}};
        q.detail.support = requested;
        q.detail.physical = physical;
        q.math.cancel();
        const auto& schema = partition_take(si.source_case()).get();
        if (!same_parent(partition_take(schema.identity()), bound.sidecar().source().parent)) {
            q.reject(Reason::Binding, Code::StaleGeneration);
        }
        const auto* gmsh = std::get_if<in::GmshSidecar>(&bound.sidecar().payload());
        const auto* mesh = std::get_if<in::GmshMesh>(&bound.source());
        if (gmsh == nullptr || mesh == nullptr) {
            q.reject(Reason::Binding, Code::Unsupported);
        }
        const auto mappings = gmsh->supports.values();
        const auto row = std::ranges::find_if(mappings, [&](const auto& candidate) {
            q.math.visit(1);
            return candidate.support == requested && candidate.physical == physical;
        });
        if (row == mappings.end()) {
            q.reject(Reason::Support);
        }
        q.detail.declaration = row->node;
        if (row->physical.dimension != 1) {
            q.reject(Reason::Support, Code::Unsupported);
        }
        const auto region =
            support_region(q, partition_take(schema.root()).model.record(), requested);
        const auto shell = in::MeshShellId{region.target.record().id.unsigned_value()};
        const auto pg = in::MeshGeometryId{region.geometry.unsigned_value()};
        const auto name = partition_index(region.boundary.unsigned_value());
        q.detail.boundary = name;

        const auto& original = prepared_shell(q, set, bound, shell, pg);
        const auto patch = partition_take(original.reference().identity());
        const auto reference = partition_take(original.reference().source());
        const auto graph = partition_take(reference.graph.view());
        const auto frame = partition_take(original.reference().frame());
        const auto& partition = reference.partition;
        const ChartBoundary chart{partition_take(partition.edges()),
                                  partition_take(partition.sources()),
                                  partition_take(partition.edge_source_indices()),
                                  partition_take(partition.vertices()), name};
        const auto source = mesh->view();
        if (!source) {
            q.reject(Reason::Binding);
        }
        q.math.admit();
        auto lines = std::visit(
            [&](const auto& geometry) {
                const auto view = [&] {
                    if constexpr (std::same_as<std::decay_t<decltype(geometry)>,
                                               in::GmshPgConversionEvidence>) {
                        return geometry.geometry.view();
                    } else {
                        return geometry.view();
                    }
                }();
                if (!view || view->patch != patch ||
                    view->source.source != source->location.source ||
                    !same_parent(view->source.parent, source->location.parent) ||
                    view->digest != source->digest) {
                    q.reject(Reason::Binding, Code::StaleGeneration);
                }
                return cover(q, *source, row->physical, view->faces, chart, graph, frame,
                             limits.maximum_lines);
            },
            original.match().geometry);
        q.math.cancel();
        return CaseGmshSupportBoundary{source->location,
                                       source->digest,
                                       bound.sidecar().source(),
                                       bound.sidecar().digest(),
                                       requested,
                                       shell,
                                       row->physical,
                                       row->node,
                                       name,
                                       patch,
                                       std::move(lines)};
    });
}

} // namespace uvlas::coupling::detail
