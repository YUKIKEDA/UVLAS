#include "external_planar_pg.hpp"

#include "../structure/pg_reference_arithmetic.hpp"

#include <algorithm>
#include <compare>
#include <utility>

namespace uvlas::io::detail {
namespace {

namespace geom = structure::detail;
using geom::partition_fail;
using geom::partition_take;
using geom::ReferenceArithmetic;
using Point = structure::PgPartitionPoint;
using Rational = structure::ExactRational;
using Code = core::FailureCode;
using Reason = ExternalPgFailure;
using Ring = std::span<const ExternalPgOwnedCorner>;
using Triangle = std::array<const Point*, 3>;
using Position = std::array<Rational, 3>;
constexpr auto operation = core::Operation::PrepareModel;

struct Corner {
    std::uint64_t tag{};
    std::array<std::string_view, 3> coordinates;
    std::array<std::optional<std::uint64_t>, 3> offsets;
};

struct Face {
    std::uint64_t tag{};
    std::size_t count{};
    std::array<Corner, 4> corners;
};

struct Diagnostic {
    std::optional<ExternalPgDiagnostic>* slot;

    [[noreturn]] void reject(Code code, Reason reason, std::optional<std::uint64_t> face = {},
                             std::optional<std::uint64_t> other = {},
                             std::optional<std::uint64_t> node = {},
                             std::optional<std::uint64_t> offset = {}) const {
        if (slot) {
            *slot = ExternalPgDiagnostic{reason, face, other, node, offset};
        }
        partition_fail(code);
    }
};

[[nodiscard]] bool equal_source(const RootSnapshotLocation& a,
                                const RootSnapshotLocation& b) noexcept {
    return a.source == b.source && a.role == b.role && a.declaration_node == b.declaration_node &&
           a.parent.source == b.parent.source &&
           a.parent.contract_revision == b.parent.contract_revision &&
           a.parent.representation_revision == b.parent.representation_revision;
}

[[nodiscard]] Rational decimal(ReferenceArithmetic& math, std::string_view token) {
    math.charge();
    return partition_take(Rational::from_decimal(token, math.leaf_limits, math.plan));
}

[[nodiscard]] Rational cross(ReferenceArithmetic& math, const Point& a, const Point& b) {
    const auto left = math.mul(a.s, b.t);
    const auto right = math.mul(a.t, b.s);
    return math.sub(left, right);
}

[[nodiscard]] Rational orient(ReferenceArithmetic& math, const Point& a, const Point& b,
                              const Point& c) {
    const Point first{math.sub(b.s, a.s), math.sub(b.t, a.t)};
    const Point second{math.sub(c.s, a.s), math.sub(c.t, a.t)};
    return cross(math, first, second);
}

[[nodiscard]] Point copy(ReferenceArithmetic& math, const Point& p) {
    return {math.copy(p.s), math.copy(p.t)};
}

[[nodiscard]] const Point& positive_corner(const ExternalPgOwnedFace& face,
                                           std::size_t i) noexcept {
    const auto n = static_cast<std::size_t>(face.count);
    return face.corners[face.winding == ExternalPgWinding::Forward ? i : (n - i) % n].chart;
}

[[nodiscard]] Rational area2(ReferenceArithmetic& math, const ExternalPgOwnedFace& face) {
    Rational sum;
    const auto count = static_cast<std::size_t>(face.count);
    for (std::size_t i = 0; i < count; ++i) {
        sum = math.add(
            sum, cross(math, positive_corner(face, i), positive_corner(face, (i + 1) % count)));
    }
    return sum;
}

// This tests the interiors of two convex polygons. A shared edge is not an overlap.
[[nodiscard]] bool separated(ReferenceArithmetic& math, const ExternalPgOwnedFace& first,
                             const ExternalPgOwnedFace& second) {
    for (std::size_t i = 0; i < first.count; ++i) {
        const auto& a = positive_corner(first, i);
        const auto& b = positive_corner(first, (i + 1) % first.count);
        if (std::ranges::none_of(second.ring(), [&](const ExternalPgOwnedCorner& p) {
                return orient(math, a, b, p.chart).sign() > 0;
            })) {
            return true;
        }
    }
    return false;
}

[[nodiscard]] Point crossing(ReferenceArithmetic& math, const Point& p, const Point& q,
                             const Rational& from, const Rational& to) {
    const auto fraction = math.div(from, math.sub(from, to));
    return {math.add(p.s, math.mul(fraction, math.sub(q.s, p.s))),
            math.add(p.t, math.mul(fraction, math.sub(q.t, p.t)))};
}

// A convex triangle clipped by at most four half-planes has at most seven vertices.
[[nodiscard]] Rational intersect_area2(ReferenceArithmetic& math, const Triangle& triangle,
                                       const ExternalPgOwnedFace& face) {
    std::array<Point, 7> ring;
    std::size_t size = triangle.size();
    for (std::size_t i = 0; i < size; ++i) {
        ring[i] = copy(math, *triangle[i]);
    }
    for (std::size_t edge = 0; edge < face.count; ++edge) {
        if (size < 3) {
            return {};
        }
        const auto& a = positive_corner(face, edge);
        const auto& b = positive_corner(face, (edge + 1) % face.count);
        std::array<Point, 7> clipped;
        std::size_t count = 0;
        auto previous = orient(math, a, b, ring[size - 1]);
        for (std::size_t i = 0; i < size; ++i) {
            auto current = orient(math, a, b, ring[i]);
            if (previous.sign() * current.sign() < 0) {
                if (count == clipped.size()) {
                    partition_fail(Code::GeometryInvalid);
                }
                clipped[count++] =
                    crossing(math, ring[(i + size - 1) % size], ring[i], previous, current);
            }
            if (current.sign() >= 0) {
                if (count == clipped.size()) {
                    partition_fail(Code::GeometryInvalid);
                }
                clipped[count++] = copy(math, ring[i]);
            }
            previous = std::move(current);
        }
        ring = std::move(clipped);
        size = count;
    }
    Rational sum;
    for (std::size_t i = 0; i < size; ++i) {
        sum = math.add(sum, cross(math, ring[i], ring[(i + 1) % size]));
    }
    return sum;
}

struct AffineMidSurface {
    Position origin, tangent_s, tangent_t;

    [[nodiscard]] Rational coordinate(ReferenceArithmetic& math, const Point& p,
                                      std::size_t i) const {
        return math.add(origin[i],
                        math.add(math.mul(tangent_s[i], p.s), math.mul(tangent_t[i], p.t)));
    }
};

[[nodiscard]] AffineMidSurface mid_surface(ReferenceArithmetic& math,
                                           const structure::PgGraphView& graph,
                                           const structure::PgReferenceFrame& frame) {
    const Rational zero;
    const auto one = math.constant("1");
    const Rational *f = &zero, *fs = &zero, *ft = &zero;
    for (const auto& term : graph.height.terms()) {
        math.visit(1);
        if (term.powers.first != 0) {
            fs = &term.coefficient;
        } else if (term.powers.second != 0) {
            ft = &term.coefficient;
        } else {
            f = &term.coefficient;
        }
    }
    const auto push = [&](std::size_t i, const Rational& a, const Rational& b, const Rational& c) {
        const auto rotated = math.add(
            math.add(math.mul(frame.rotation[3 * i], a), math.mul(frame.rotation[3 * i + 1], b)),
            math.mul(frame.rotation[3 * i + 2], c));
        return math.mul(frame.length_metres, rotated);
    };
    AffineMidSurface result;
    for (std::size_t i = 0; i < result.origin.size(); ++i) {
        result.origin[i] = math.add(frame.origin_metres[i], push(i, graph.x0, zero, *f));
        result.tangent_s[i] = push(i, graph.width, zero, *fs);
        result.tangent_t[i] = push(i, graph.skew, one, *ft);
    }
    return result;
}

void fill_corner(ReferenceArithmetic& math, const Diagnostic& diagnostic,
                 const AffineMidSurface& physical, const Rational& unit, const Face& source,
                 std::size_t j, const ExternalPgCorner& declared, ExternalPgOwnedCorner& row) {
    const auto& input = source.corners[j];
    if (declared.node != input.tag) {
        diagnostic.reject(Code::InvalidInput, Reason::SourceOrder, source.tag, {}, input.tag);
    }
    row.node = input.tag;
    row.chart = copy(math, declared.point);
    row.source_bytes = input.offsets;
    for (std::size_t i = 0; i < input.coordinates.size(); ++i) {
        auto metres = math.mul(decimal(math, input.coordinates[i]), unit);
        const auto expected = physical.coordinate(math, row.chart, i);
        if (math.compare(metres, expected) != 0) {
            diagnostic.reject(Code::GeometryInvalid, Reason::Position, source.tag, {}, input.tag,
                              input.offsets[i]);
        }
        row.position_metres[i] = std::move(metres);
    }
}

struct Prepared {
    structure::PgReferenceIdentity patch;
    core::PlannedBuffer<ExternalPgOwnedFace> faces;
    structure::PgReferenceCallInfo work;
};

template <class ReadFace>
[[nodiscard]] Prepared match(const structure::PgReferenceGeometry& target, SiLengthUnit length,
                             std::span<const ExternalPgFace> points, const ExternalPgLimits& limits,
                             core::ResourcePlan& plan, const Diagnostic& diagnostic,
                             ReadFace&& read_face) {
    ReferenceArithmetic math{limits.arithmetic, plan};
    math.cancel();
    if (points.empty()) {
        diagnostic.reject(Code::InvalidInput, Reason::Selection);
    }
    if (points.size() > limits.maximum_faces) {
        partition_fail(Code::ResourceExceeded);
    }
    const auto owned_count = geom::partition_index(points.size());
    const auto owner = partition_take(target.source());
    const auto graph = partition_take(owner.graph.view());
    // Check the degree before allocating even when all supplied vertices lie on the graph.
    for (const auto& term : graph.height.terms()) {
        math.visit(1);
        if (term.powers.first + term.powers.second > 1) {
            diagnostic.reject(Code::Unsupported, Reason::CurvedPatch);
        }
    }
    const auto frame = partition_take(target.frame());
    const auto patch = partition_take(target.identity());
    const auto triangles = partition_take(owner.partition.triangles());
    const auto vertices = partition_take(owner.partition.vertices());
    const auto slots =
        geom::partition_add(geom::reference_live_rationals, geom::partition_mul(20, points.size()));
    math.admit(slots);
    const auto physical = mid_surface(math, graph, frame);
    const SiUnits units{length, SiMassUnit::Kilogram, SiTimeUnit::Second, SiAngleUnit::Radian};
    const SiRecipeView scale{SiDecimalView{1, "1", 0}, {}, units, {{1, 0, 0, 0}}, 0};
    math.charge();
    const auto unit = partition_take(materialize_rational(scale, math.leaf_limits, plan));
    auto faces = partition_take(plan.allocate<ExternalPgOwnedFace>(owned_count, operation));
    Rational total;
    for (std::size_t i = 0; i < points.size(); ++i) {
        math.visit(1);
        const auto& given = points[i];
        const auto source = read_face(math, i, given);
        if (source.tag != given.element || source.count != given.corners.size()) {
            diagnostic.reject(Code::InvalidInput, Reason::Selection, given.element);
        }
        auto& row = faces.values()[i];
        row.element = source.tag;
        row.count = static_cast<std::uint8_t>(source.count);
        row.winding = given.winding;
        for (std::size_t j = 0; j < source.count; ++j) {
            fill_corner(math, diagnostic, physical, unit, source, j, given.corners[j],
                        row.corners[j]);
        }
        if (row.winding != ExternalPgWinding::Forward &&
            row.winding != ExternalPgWinding::Reversed) {
            diagnostic.reject(Code::InvalidInput, Reason::Winding, row.element);
        }
        for (std::size_t j = 0; j < source.count; ++j) {
            const auto& a = positive_corner(row, j);
            const auto& b = positive_corner(row, (j + 1) % source.count);
            const auto& c = positive_corner(row, (j + 2) % source.count);
            // These are the corner Jacobian signs of the bilinear map for a quad.
            if (orient(math, a, b, c).sign() <= 0) {
                diagnostic.reject(Code::GeometryInvalid, Reason::Winding, row.element);
            }
        }
        for (std::size_t j = 0; j < i; ++j) {
            math.visit(1);
            const auto& other = faces.values()[j];
            if (!separated(math, row, other) && !separated(math, other, row)) {
                diagnostic.reject(Code::GeometryInvalid, Reason::Overlap, row.element,
                                  other.element);
            }
        }
        const auto area = area2(math, row);
        Rational inside;
        for (const auto& triangle : triangles) {
            math.visit(1);
            const Triangle t{&vertices[static_cast<std::size_t>(triangle.vertices[0])],
                             &vertices[static_cast<std::size_t>(triangle.vertices[1])],
                             &vertices[static_cast<std::size_t>(triangle.vertices[2])]};
            inside = math.add(inside, intersect_area2(math, t, row));
        }
        if (math.compare(inside, area) != 0) {
            diagnostic.reject(Code::GeometryInvalid, Reason::Outside, row.element);
        }
        total = math.add(total, area);
    }

    Rational reference_area;
    for (const auto& triangle : triangles) {
        math.visit(1);
        reference_area = math.add(reference_area, triangle.twice_area);
    }
    if (math.compare(total, reference_area) != 0) {
        diagnostic.reject(Code::GeometryInvalid, Reason::Coverage);
    }
    math.cancel();
    return {patch, std::move(faces), math.info()};
}

[[nodiscard]] Face gmsh_face(ReferenceArithmetic& math, const Diagnostic& diagnostic,
                             const GmshMeshView& mesh, const GmshElement& element) {
    const auto ordered = element.connectivity();
    Face result;
    result.tag = element.tag;
    result.count = ordered.size();
    for (std::size_t j = 0; j < result.count; ++j) {
        const auto node = std::ranges::lower_bound(
            mesh.nodes, ordered[j],
            [&](std::uint64_t a, std::uint64_t b) {
                math.visit(1);
                return a < b;
            },
            &GmshNode::tag);
        if (node == mesh.nodes.end() || node->tag != ordered[j]) {
            diagnostic.reject(Code::InvalidInput, Reason::Source, element.tag, {}, ordered[j]);
        }
        auto& out = result.corners[j];
        out.tag = node->tag;
        for (std::size_t k = 0; k < out.coordinates.size(); ++k) {
            out.coordinates[k] = mesh.text(node->coordinates[k]);
            out.offsets[k] = node->coordinates[k].offset;
        }
    }
    return result;
}

[[nodiscard]] Face nastran_face(ReferenceArithmetic& math, const Diagnostic& diagnostic,
                                const NastranShellView& deck, std::uint64_t id) {
    const auto found = std::ranges::find_if(deck.elements, [&](const NastranShellElement& row) {
        math.visit(1);
        return row.id.value() == id;
    });
    if (found == deck.elements.end()) {
        diagnostic.reject(Code::InvalidInput, Reason::Selection, id);
    }
    if (found->offset && decimal(math, found->offset->value.decimal()).sign() != 0) {
        diagnostic.reject(Code::Unsupported, Reason::Offset, id, {}, {},
                          found->offset->field.position.byte);
    }
    return std::visit(
        [&](const auto& corners) {
            Face result;
            result.tag = id;
            result.count = corners.size();
            for (std::size_t j = 0; j < corners.size(); ++j) {
                const auto& grid = deck.input.grids[corners[j].grid_ordinal];
                auto& out = result.corners[j];
                out.tag = grid.id.value();
                if (grid.position_system.id.value() != 0) {
                    const auto& cp = grid.position_system.field;
                    diagnostic.reject(Code::Unsupported, Reason::PositionSystem, id, {}, out.tag,
                                      cp ? cp->position.byte : grid.id_field.position.byte);
                }
                for (std::size_t k = 0; k < out.coordinates.size(); ++k) {
                    out.coordinates[k] = grid.position[k].decimal();
                    if (const auto& field = grid.position[k].field) {
                        out.offsets[k] = field->position.byte;
                    }
                }
            }
            return result;
        },
        found->corners);
}

} // namespace

ExternalPlanarPg::ExternalPlanarPg(RootSnapshotLocation source, SourceDigest digest,
                                   SiLengthUnit units, structure::PgReferenceIdentity patch,
                                   core::PlannedBuffer<ExternalPgOwnedFace>&& faces,
                                   structure::PgReferenceCallInfo work) noexcept
    : source_(source),
      digest_(digest),
      units_(units),
      patch_(patch),
      faces_(std::move(faces)),
      work_(work) {}

std::optional<ExternalPlanarPgView> ExternalPlanarPg::view() const& noexcept {
    if (faces_.size() == 0) {
        return std::nullopt;
    }
    return ExternalPlanarPgView{source_, digest_, units_, patch_, faces_.values(), work_};
}

core::Result<ExternalPlanarPg>
ExternalPlanarPg::prepare(const structure::PgReferenceGeometry& pg, const GmshMesh& mesh,
                          const GmshSurfaceTopology& topology, SiLengthUnit units,
                          std::span<const ExternalPgFace> charts, const ExternalPgLimits& limits,
                          core::ResourcePlan& plan,
                          std::optional<ExternalPgDiagnostic>* slot) noexcept {
    if (slot) {
        slot->reset();
    }
    const Diagnostic diagnostic{slot};
    return geom::partition_entry([&]() -> core::Result<ExternalPlanarPg> {
        if (limits.arithmetic.cancellation.stop_requested()) {
            partition_fail(Code::Cancelled);
        }
        const auto input = mesh.view();
        const auto selected = topology.view();
        if (!input || !selected) {
            diagnostic.reject(Code::InvalidInput, Reason::Source);
        }
        if (!equal_source(input->location, selected->source.location) ||
            input->digest != selected->source.digest) {
            diagnostic.reject(Code::StaleGeneration, Reason::Source);
        }
        if (selected->faces.size() != charts.size()) {
            diagnostic.reject(Code::InvalidInput, Reason::Selection);
        }
        auto result =
            match(pg, units, charts, limits, plan, diagnostic,
                  [&](ReferenceArithmetic& math, std::size_t i, const ExternalPgFace&) {
                      return gmsh_face(math, diagnostic, *input, selected->faces[i].element);
                  });
        return ExternalPlanarPg{input->location, input->digest,           units,
                                result.patch,    std::move(result.faces), result.work};
    });
}

core::Result<ExternalPlanarPg>
ExternalPlanarPg::prepare(const structure::PgReferenceGeometry& pg, const NastranShellCards& deck,
                          SiLengthUnit units, std::span<const ExternalPgFace> charts,
                          const ExternalPgLimits& limits, core::ResourcePlan& plan,
                          std::optional<ExternalPgDiagnostic>* slot) noexcept {
    if (slot) {
        slot->reset();
    }
    const Diagnostic diagnostic{slot};
    return geom::partition_entry([&]() -> core::Result<ExternalPlanarPg> {
        if (limits.arithmetic.cancellation.stop_requested()) {
            partition_fail(Code::Cancelled);
        }
        const auto input = deck.view();
        if (!input) {
            diagnostic.reject(Code::InvalidInput, Reason::Source);
        }
        auto result =
            match(pg, units, charts, limits, plan, diagnostic,
                  [&](ReferenceArithmetic& math, std::size_t, const ExternalPgFace& chart) {
                      return nastran_face(math, diagnostic, *input, chart.element);
                  });
        return ExternalPlanarPg{input->input.bulk.location,
                                input->input.bulk.digest,
                                units,
                                result.patch,
                                std::move(result.faces),
                                result.work};
    });
}

} // namespace uvlas::io::detail
