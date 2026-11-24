#include "external_planar_pg.hpp"

#include "../structure/pg_conditioning_factor.hpp"
#include "../structure/pg_reference_arithmetic.hpp"
#include "nastran_position_arithmetic.hpp"
#include "nastran_position_enclosures.hpp"

#include <algorithm>
#include <compare>
#include <stdexcept>
#include <tuple>
#include <type_traits>
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
using Triangle = std::array<const Point*, 3>;
using Position = std::array<Rational, 3>;
constexpr auto operation = core::Operation::PrepareModel;
// Each ring retains (s,t) and SI Cartesian coordinates. With ZOFFS the
// deck may also retain the four unshifted GRID positions and the scalar.
constexpr std::uint64_t gmsh_rationals_per_face = std::uint64_t{4} * (2 + 3);
constexpr std::uint64_t nastran_rationals_per_face =
    gmsh_rationals_per_face + std::uint64_t{4} * 3 + 1;

struct Corner {
    std::uint64_t tag{};
    Position grid_metres;
    std::optional<NastranGridPlacement> nastran;
    std::array<std::optional<std::uint64_t>, 3> offsets;
};

struct Face {
    std::uint64_t tag{};
    std::size_t count{};
    std::array<Corner, 4> corners;
    std::optional<NastranReferenceOffset> offset;
    Position shift_metres;
};

struct Diagnostic {
    std::optional<ExternalPgDiagnostic>* slot;

    [[noreturn]] void reject(Code code, Reason reason, std::optional<std::uint64_t> face = {},
                             std::optional<std::uint64_t> other = {},
                             std::optional<std::uint64_t> node = {},
                             std::optional<std::uint64_t> offset = {},
                             std::optional<NastranCoordinateId> coordinate = {}) const {
        if (slot) {
            *slot = ExternalPgDiagnostic{reason, face, other, node, offset, coordinate};
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

[[nodiscard]] Position vector_sub(ReferenceArithmetic& math, const Position& a, const Position& b) {
    return {math.sub(a[0], b[0]), math.sub(a[1], b[1]), math.sub(a[2], b[2])};
}

[[nodiscard]] Position vector_cross(ReferenceArithmetic& math, const Position& a,
                                    const Position& b) {
    Position result;
    for (std::size_t i = 0; i < 3; ++i) {
        const auto j = (i + 1) % 3, k = (i + 2) % 3;
        result[i] = math.sub(math.mul(a[j], b[k]), math.mul(a[k], b[j]));
    }
    return result;
}

[[nodiscard]] Rational dot(ReferenceArithmetic& math, const Position& a, const Position& b) {
    return math.add(math.add(math.mul(a[0], b[0]), math.mul(a[1], b[1])), math.mul(a[2], b[2]));
}

struct DirectionSite {
    const Diagnostic& diagnostic;
    Reason reason = Reason::CoordinateFrame;
    std::uint64_t element{};
    std::optional<std::uint64_t> node, byte;
    std::optional<NastranCoordinateId> coordinate;
};

[[nodiscard]] Position normalize(ReferenceArithmetic& math, const Position& direction,
                                 const DirectionSite& site) {
    const auto squared = dot(math, direction, direction);
    const auto reject = [&](Code code) {
        site.diagnostic.reject(code, site.reason, site.element, {}, site.node, site.byte,
                               site.coordinate);
    };
    if (squared.sign() <= 0) {
        reject(Code::GeometryInvalid);
    }
    // The nearest integer can round above the root. Check its exact square before use.
    const auto n = squared.numerator_bits(), d = squared.denominator_bits();
    if (n > 8191 || d > 8191) {
        reject(Code::ResourceExceeded);
    }
    constexpr auto limb_bits = std::numeric_limits<geom::RankWord>::digits;
    math.charge(2 + n / limb_bits + d / limb_bits);
    const auto fraction = geom::rank_fraction<geom::RankWide>(squared);
    const auto a = geom::conditioning_root(math, fraction.numerator);
    const auto b = geom::conditioning_root(math, fraction.denominator);
    if (geom::rank_product(math, a, a) != fraction.numerator ||
        geom::rank_product(math, b, b) != fraction.denominator) {
        reject(Code::PrecisionUnsupported);
    }
    math.charge();
    const auto radius =
        geom::rank_exact(geom::RankFraction<geom::RankWide>{a, b}, math.leaf_limits, math.plan);
    return {math.div(direction[0], radius), math.div(direction[1], radius),
            math.div(direction[2], radius)};
}

[[nodiscard]] Position tuple(ReferenceArithmetic& math,
                             const std::array<NastranCoordinateValue, 3>& recorded) {
    return {decimal(math, recorded[0].decimal()), decimal(math, recorded[1].decimal()),
            decimal(math, recorded[2].decimal())};
}

[[nodiscard]] std::uint64_t first_byte(const std::array<NastranCoordinateValue, 3>& tuple,
                                       std::uint64_t fallback) noexcept {
    for (const auto& item : tuple) {
        if (const auto& field = item.field) {
            return field->position.byte;
        }
    }
    return fallback;
}

[[nodiscard]] const NastranCoordinateDefinition*
rectangular(ReferenceArithmetic& math, const Diagnostic& diag, const NastranCoordinateView& deck,
            const NastranCoordinateReference& reference, std::uint64_t element, std::uint64_t grid,
            std::uint64_t fallback) {
    if (reference.id.value() == 0) {
        return nullptr;
    }
    const auto found =
        std::ranges::find_if(deck.coordinates, [&](const NastranCoordinateDefinition& row) {
            math.visit(1);
            return row.id == reference.id;
        });
    const auto& token = reference.field;
    const auto offset = token ? token->position.byte : fallback;
    if (found == deck.coordinates.end()) {
        diag.reject(Code::InvalidInput, Reason::PositionSystem, element, {}, grid, offset,
                    reference.id);
    }
    if (found->kind != NastranCoordinateKind::Rectangular) {
        diag.reject(Code::Unsupported, Reason::PositionSystem, element, {}, grid, offset,
                    reference.id);
    }
    return &*found;
}

[[nodiscard]] Position basic_grid(ReferenceArithmetic& math, const Diagnostic& diagnostic,
                                  const NastranCoordinateView& deck, const NastranGridPoint& grid,
                                  std::uint64_t element, const Rational& unit) {
    auto result = tuple(math, grid.position);
    const auto id = grid.id.value();
    auto* frame = rectangular(math, diagnostic, deck, grid.position_system, element, id,
                              grid.id_field.position.byte);
    // Every definition below is interpreted in its rectangular RID dimensions.
    // The owner has already rejected RID cycles. Walking toward basic needs no stack.
    while (frame) {
        math.visit(1);
        const auto fallback = frame->name.position.byte;
        const auto* parent =
            rectangular(math, diagnostic, deck, frame->reference, element, id, fallback);
        const auto a = tuple(math, frame->a);
        const auto b = tuple(math, frame->b);
        const auto c = tuple(math, frame->c);
        const auto z = normalize(math, vector_sub(math, b, a),
                                 {diagnostic, Reason::CoordinateFrame, element, id,
                                  first_byte(frame->b, fallback), frame->id});
        const auto y = normalize(math, vector_cross(math, z, vector_sub(math, c, a)),
                                 {diagnostic, Reason::CoordinateFrame, element, id,
                                  first_byte(frame->c, fallback), frame->id});
        const auto x = vector_cross(math, y, z);
        Position parent_coordinates;
        for (std::size_t i = 0; i < 3; ++i) {
            parent_coordinates[i] = math.add(
                a[i], math.add(math.add(math.mul(x[i], result[0]), math.mul(y[i], result[1])),
                               math.mul(z[i], result[2])));
        }
        result = std::move(parent_coordinates);
        frame = parent;
    }
    for (auto& axis : result) {
        axis = math.mul(axis, unit);
    }
    return result;
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

[[nodiscard]] const Point& positive_corner(const auto& face, std::size_t i) noexcept {
    const auto n = static_cast<std::size_t>(face.count);
    return face.corners[face.winding == ExternalPgWinding::Forward ? i : (n - i) % n].chart;
}

[[nodiscard]] Rational area2(ReferenceArithmetic& math, const auto& face) {
    Rational sum;
    const auto count = static_cast<std::size_t>(face.count);
    for (std::size_t i = 0; i < count; ++i) {
        sum = math.add(
            sum, cross(math, positive_corner(face, i), positive_corner(face, (i + 1) % count)));
    }
    return sum;
}

// This tests the interiors of two convex polygons. A shared edge is not an overlap.
[[nodiscard]] bool separated(ReferenceArithmetic& math, const auto& first, const auto& second) {
    for (std::size_t i = 0; i < first.count; ++i) {
        const auto& a = positive_corner(first, i);
        const auto& b = positive_corner(first, (i + 1) % first.count);
        if (std::ranges::none_of(second.ring(), [&](const auto& p) {
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
                                       const auto& face) {
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
        if (term.powers.first + term.powers.second > 1) {
            continue;
        }
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
                 const AffineMidSurface& physical, Face& source, std::size_t j,
                 const ExternalPgCorner& declared, ExternalPgOwnedCorner& row) {
    auto& input = source.corners[j];
    if (declared.node != input.tag) {
        diagnostic.reject(Code::InvalidInput, Reason::SourceOrder, source.tag, {}, input.tag);
    }
    row.node = input.tag;
    row.chart = copy(math, declared.point);
    row.source_bytes = input.offsets;
    const bool shifted = source.offset && source.offset->signed_metres.sign() != 0;
    for (std::size_t i = 0; i < input.grid_metres.size(); ++i) {
        auto metres = shifted ? math.add(input.grid_metres[i], source.shift_metres[i])
                              : std::move(input.grid_metres[i]);
        const auto expected = physical.coordinate(math, row.chart, i);
        if (math.compare(metres, expected) != 0) {
            auto position = input.offsets[i];
            std::optional<NastranCoordinateId> system;
            if (const auto& defined = input.nastran) {
                system = defined->position_system;
                if (system->value() != 0) {
                    // Any source coordinate may affect this Cartesian component.
                    position = defined->cp_byte;
                }
            }
            if (const auto& offset = source.offset; shifted && offset) {
                position = offset->byte;
            }
            diagnostic.reject(Code::GeometryInvalid, Reason::Position, source.tag, {}, input.tag,
                              position, system);
        }
        row.position_metres[i] = std::move(metres);
    }
    if (input.nastran) {
        if (shifted) {
            input.nastran->unshifted_metres = std::move(input.grid_metres);
        }
        row.nastran = std::move(input.nastran);
    }
}

void validate_chart_row(ReferenceArithmetic& math, const Diagnostic& diagnostic, const auto& row,
                        auto previous, std::span<const structure::PgPartitionTriangle> triangles,
                        std::span<const Point> vertices, Rational& total) {
    if (row.winding != ExternalPgWinding::Forward && row.winding != ExternalPgWinding::Reversed) {
        diagnostic.reject(Code::InvalidInput, Reason::Winding, row.element);
    }
    for (std::size_t j = 0; j < row.count; ++j) {
        const auto& a = positive_corner(row, j);
        const auto& b = positive_corner(row, (j + 1) % row.count);
        const auto& c = positive_corner(row, (j + 2) % row.count);
        // These are the corner Jacobian signs of the bilinear map for a quad.
        if (orient(math, a, b, c).sign() <= 0) {
            diagnostic.reject(Code::GeometryInvalid, Reason::Winding, row.element);
        }
    }
    for (const auto& other : previous) {
        math.visit(1);
        if (!separated(math, row, other) && !separated(math, other, row)) {
            diagnostic.reject(Code::GeometryInvalid, Reason::Overlap, row.element, other.element);
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

void validate_chart_coverage(ReferenceArithmetic& math, const Diagnostic& diagnostic,
                             std::span<const structure::PgPartitionTriangle> triangles,
                             const Rational& total) {
    Rational reference_area;
    for (const auto& triangle : triangles) {
        math.visit(1);
        reference_area = math.add(reference_area, triangle.twice_area);
    }
    if (math.compare(total, reference_area) != 0) {
        diagnostic.reject(Code::GeometryInvalid, Reason::Coverage);
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
                             std::uint64_t rationals_per_face, ReadFace&& read_face) {
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
    const auto slots = geom::partition_add(geom::reference_live_rationals,
                                           geom::partition_mul(rationals_per_face, points.size()));
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
        auto source = read_face(math, i, given, unit);
        if (source.tag != given.element || source.count != given.corners.size()) {
            diagnostic.reject(Code::InvalidInput, Reason::Selection, given.element);
        }
        auto& row = faces.values()[i];
        row.element = source.tag;
        row.count = static_cast<std::uint8_t>(source.count);
        row.winding = given.winding;
        for (std::size_t j = 0; j < source.count; ++j) {
            fill_corner(math, diagnostic, physical, source, j, given.corners[j], row.corners[j]);
        }
        row.nastran_offset = std::move(source.offset);
        validate_chart_row(math, diagnostic, row, faces.values().first(i), triangles, vertices,
                           total);
    }
    validate_chart_coverage(math, diagnostic, triangles, total);
    math.cancel();
    return {patch, std::move(faces), math.info()};
}

[[nodiscard]] const GmshNode& find_gmsh_node(ReferenceArithmetic& math,
                                             const Diagnostic& diagnostic, const GmshMeshView& mesh,
                                             std::uint64_t element, std::uint64_t id) {
    const auto node = std::ranges::lower_bound(
        mesh.nodes, id,
        [&](std::uint64_t a, std::uint64_t b) {
            math.visit(1);
            return a < b;
        },
        &GmshNode::tag);
    if (node == mesh.nodes.end() || node->tag != id) {
        diagnostic.reject(Code::InvalidInput, Reason::Source, element, {}, id);
    }
    return *node;
}

[[nodiscard]] Position gmsh_coordinates(ReferenceArithmetic& math, const GmshMeshView& mesh,
                                        const GmshNode& node, const Rational& unit) {
    Position position;
    for (std::size_t k = 0; k < position.size(); ++k) {
        position[k] = math.mul(decimal(math, mesh.text(node.coordinates[k])), unit);
    }
    return position;
}

[[nodiscard]] Face gmsh_face(ReferenceArithmetic& math, const Diagnostic& diagnostic,
                             const GmshMeshView& mesh, const GmshElement& element,
                             const Rational& unit) {
    const auto ordered = element.connectivity();
    Face result;
    result.tag = element.tag;
    result.count = ordered.size();
    for (std::size_t j = 0; j < result.count; ++j) {
        const auto& node = find_gmsh_node(math, diagnostic, mesh, element.tag, ordered[j]);
        auto& out = result.corners[j];
        out.tag = node.tag;
        out.grid_metres = gmsh_coordinates(math, mesh, node, unit);
        for (std::size_t k = 0; k < out.offsets.size(); ++k) {
            out.offsets[k] = node.coordinates[k].offset;
        }
    }
    return result;
}

void apply_offset(ReferenceArithmetic& math, const Diagnostic& diagnostic, Face& face,
                  const NastranShellReal& offset, const Rational& unit) {
    const auto byte = offset.field.position.byte;
    auto length = math.mul(decimal(math, offset.value.decimal()), unit);
    if (length.sign() != 0) {
        const auto& origin = face.corners[0].grid_metres;
        const auto cross = vector_cross(math, vector_sub(math, face.corners[1].grid_metres, origin),
                                        vector_sub(math, face.corners[2].grid_metres, origin));
        if (face.count == 4 &&
            dot(math, cross, vector_sub(math, face.corners[3].grid_metres, origin)).sign() != 0) {
            diagnostic.reject(Code::GeometryInvalid, Reason::Offset, face.tag, {},
                              face.corners[3].tag, byte);
        }
        const auto direction =
            normalize(math, cross, {diagnostic, Reason::Offset, face.tag, {}, byte, {}});
        for (std::size_t i = 0; i < direction.size(); ++i) {
            face.shift_metres[i] = math.mul(length, direction[i]);
        }
    }
    face.offset = NastranReferenceOffset{std::move(length), byte};
}

[[nodiscard]] const NastranShellElement& find_shell(ReferenceArithmetic& math,
                                                    const Diagnostic& diagnostic,
                                                    const NastranShellView& deck,
                                                    std::uint64_t id) {
    const auto found = std::ranges::find_if(deck.elements, [&](const NastranShellElement& row) {
        math.visit(1);
        return row.id.value() == id;
    });
    if (found == deck.elements.end()) {
        diagnostic.reject(Code::InvalidInput, Reason::Selection, id);
    }
    return *found;
}

[[nodiscard]] Face nastran_face(ReferenceArithmetic& math, const Diagnostic& diagnostic,
                                const NastranShellView& deck, std::uint64_t id,
                                const Rational& unit) {
    const auto& shell = find_shell(math, diagnostic, deck, id);
    auto face = std::visit(
        [&](const auto& corners) {
            Face result;
            result.tag = id;
            result.count = corners.size();
            for (std::size_t j = 0; j < corners.size(); ++j) {
                const auto& grid = deck.input.grids[corners[j].grid_ordinal];
                auto& out = result.corners[j];
                out.tag = grid.id.value();
                out.grid_metres = basic_grid(math, diagnostic, deck.input, grid, id, unit);
                const auto& cp = grid.position_system.field;
                out.nastran = NastranGridPlacement{
                    grid.position_system.id,
                    cp ? std::optional<std::uint64_t>{cp->position.byte} : std::nullopt,
                    {}};
                for (std::size_t k = 0; k < grid.position.size(); ++k) {
                    if (const auto& field = grid.position[k].field) {
                        out.offsets[k] = field->position.byte;
                    }
                }
            }
            return result;
        },
        shell.corners);
    if (const auto& offset = shell.offset) {
        apply_offset(math, diagnostic, face, *offset, unit);
    }
    return face;
}

namespace directed = nastran_directed;
namespace m = directed::m;
using IntervalPosition = directed::Position;

[[nodiscard]] m::Interval dot(const IntervalPosition& a, const IntervalPosition& b, m::Work& work) {
    return m::add(m::add(m::multiply(a[0], b[0], work), m::multiply(a[1], b[1], work), work),
                  m::multiply(a[2], b[2], work), work);
}

void require_positive(const m::Interval& value, const Diagnostic& diagnostic,
                      const NastranPgFaceBound& row) {
    if (value.lower.mantissa <= 0) {
        diagnostic.reject(value.upper.mantissa <= 0 ? Code::GeometryInvalid
                                                    : Code::PrecisionUnresolved,
                          Reason::SourceWinding, row.element, {}, {}, row.element_field.byte);
    }
}

[[nodiscard]] IntervalPosition to_interval(const Position& value, m::Work& work) {
    return {directed::number(value[0], work), directed::number(value[1], work),
            directed::number(value[2], work)};
}

[[nodiscard]] IntervalPosition to_interval(const NastranGridEnclosure& grid, m::Work& work) {
    IntervalPosition result;
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = {directed::number(grid.metres[i].lower, work).lower,
                     directed::number(grid.metres[i].upper, work).upper};
    }
    return result;
}

void prove_planarity(ReferenceArithmetic& math, const Diagnostic& diagnostic,
                     const NastranShellView& deck, const std::array<NastranShellCorner, 4>& corners,
                     const NastranShellElement& shell) {
    const auto& first = deck.input.grids[corners[0].grid_ordinal];
    const auto cp = first.position_system.id;
    if (cp.value() != 0) {
        static_cast<void>(rectangular(math, diagnostic, deck.input, first.position_system,
                                      shell.id.value(), first.id.value(),
                                      first.id_field.position.byte));
    }

    std::array<Position, 4> points;
    for (std::size_t i = 0; i < corners.size(); ++i) {
        const auto& grid = deck.input.grids[corners[i].grid_ordinal];
        if (grid.position_system.id != cp) {
            diagnostic.reject(Code::Unsupported, Reason::Planarity, shell.id.value(), {},
                              grid.id.value(), grid.id_field.position.byte,
                              grid.position_system.id);
        }
        points[i] = tuple(math, grid.position);
    }
    const auto normal = vector_cross(math, vector_sub(math, points[1], points[0]),
                                     vector_sub(math, points[2], points[0]));
    if (dot(math, normal, vector_sub(math, points[3], points[0])).sign() != 0) {
        diagnostic.reject(Code::GeometryInvalid, Reason::Planarity, shell.id.value(), {},
                          corners[3].id.value(), corners[3].field.position.byte);
    }
}

[[nodiscard]] m::Bound upper_norm(const IntervalPosition& value, m::Work& work) {
    const auto squared = directed::norm_squared(value, work);
    return m::square_root({{}, squared.upper}, work).upper;
}

struct PgSurfaceBound {
    const AffineMidSurface& affine;
    const structure::PgGraphView& graph;
    const structure::PgReferenceFrame& frame;
    const IntervalPosition& affine_normal;
    bool curved;
};

[[nodiscard]] m::Interval coordinate_range(const auto ring, const Rational Point::* coordinate,
                                           m::Work& work) {
    auto range = directed::number(ring.front().chart.*coordinate, work);
    for (const auto& corner : ring.subspan(1)) {
        const auto value = directed::number(corner.chart.*coordinate, work);
        if (m::compare(value.lower, range.lower) < 0)
            range.lower = value.lower;
        if (m::compare(value.upper, range.upper) > 0)
            range.upper = value.upper;
    }
    return range;
}

[[nodiscard]] m::Interval height_partial(const structure::PgGraphView& graph, const m::Interval& s,
                                         const m::Interval& t, unsigned ds, unsigned dt,
                                         bool nonlinear_only, m::Work& work) {
    m::Interval sum;
    for (const auto& term : graph.height.terms()) {
        work.step();
        const auto [i, j] = term.powers;
        if (i < ds || j < dt || (nonlinear_only && i + j <= 1))
            continue;
        auto value = directed::number(term.coefficient, work);
        for (unsigned k = 0; k < ds; ++k)
            value = m::multiply(value, m::integer(i - k), work);
        for (unsigned k = 0; k < dt; ++k)
            value = m::multiply(value, m::integer(j - k), work);
        for (auto k = ds; k < i; ++k)
            value = m::multiply(value, s, work);
        for (auto k = dt; k < j; ++k)
            value = m::multiply(value, t, work);
        sum = m::add(sum, value, work);
    }
    return sum;
}

[[nodiscard]] IntervalPosition curved_pg_normal(const PgSurfaceBound& surface, const m::Interval& s,
                                                const m::Interval& t, m::Work& work) {
    const auto fs = height_partial(surface.graph, s, t, 1, 0, false, work);
    const auto ft = height_partial(surface.graph, s, t, 0, 1, false, work);
    const auto width = directed::number(surface.graph.width, work);
    const auto skew = directed::number(surface.graph.skew, work);
    // Cross product of (w,0,f_s) and (a,1,f_t) in the local graph frame.
    const IntervalPosition raw{
        m::subtract({}, fs, work),
        m::subtract(m::multiply(skew, fs, work), m::multiply(width, ft, work), work), width};
    m::Interval length_squared;
    for (const auto& component : raw) {
        const auto absolute = m::absolute(component);
        length_squared = m::add(length_squared, m::multiply(absolute, absolute, work), work);
    }
    const auto length = m::square_root(length_squared, work);
    IntervalPosition world;
    for (std::size_t i = 0; i < world.size(); ++i) {
        for (std::size_t k = 0; k < raw.size(); ++k) {
            const auto rotation = directed::number(surface.frame.rotation[3 * i + k], work);
            world[i] = m::add(world[i], m::multiply(rotation, raw[k], work), work);
        }
        world[i] = m::divide(world[i], length, work);
    }
    return world;
}

[[nodiscard]] m::Interval interpolation_error(const PgSurfaceBound& surface, const m::Interval& s,
                                              const m::Interval& t, m::Work& work) {
    const auto ds = m::subtract({s.upper, s.upper}, {s.lower, s.lower}, work);
    const auto dt = m::subtract({t.upper, t.upper}, {t.lower, t.lower}, work);
    const auto hessian_term = [&](unsigned i, unsigned j, const auto& a, const auto& b) {
        const auto value = m::absolute(height_partial(surface.graph, s, t, i, j, false, work));
        return m::multiply(value, m::multiply(a, b, work), work);
    };
    const auto ss = hessian_term(2, 0, ds, ds);
    const auto st = hessian_term(1, 1, ds, dt);
    const auto tt = hessian_term(0, 2, dt, dt);
    // The interpolation weights have mean (s,t). Their two coordinate
    // variances are at most width^2/4; Taylor contributes the other factor 1/2.
    const auto error = m::scale(m::add(m::add(ss, tt, work), m::scale(st, 1), work), -3);
    return m::multiply(directed::number(surface.frame.length_metres, work), error, work);
}

[[nodiscard]] Rational distance_bound(const PgSurfaceBound& surface, const auto ring,
                                      std::span<const IntervalPosition> corners,
                                      const IntervalPosition& shift, const m::Interval& s,
                                      const m::Interval& t, ReferenceArithmetic& math,
                                      m::Work& work) {
    // Both corner interpolants use the same nonnegative weights on the whole face.
    m::Bound distance;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        IntervalPosition residual;
        const auto& chart = ring[i].chart;
        auto height = m::Interval{};
        if (surface.curved) {
            height = height_partial(surface.graph, directed::number(chart.s, work),
                                    directed::number(chart.t, work), 0, 0, true, work);
            height = m::multiply(directed::number(surface.frame.length_metres, work), height, work);
        }
        for (std::size_t k = 0; k < residual.size(); ++k) {
            const auto actual = m::add(corners[i][k], shift[k], work);
            auto expected = directed::number(surface.affine.coordinate(math, chart, k), work);
            if (surface.curved) {
                const auto direction = directed::number(surface.frame.rotation[3 * k + 2], work);
                expected = m::add(expected, m::multiply(direction, height, work), work);
            }
            residual[k] = m::subtract(actual, expected, work);
        }
        auto corner = upper_norm(residual, work);
        if (m::compare(corner, distance) > 0) {
            distance = std::move(corner);
        }
    }
    if (surface.curved)
        distance = m::add({{}, distance}, interpolation_error(surface, s, t, work), work).upper;
    return directed::retain(distance, math.leaf_limits, math.plan, work);
}

void fill_bounds(const NastranShellCards& deck_owner, const NastranShellView& deck,
                 const NastranShellElement& shell, SiLengthUnit units, const Rational& unit,
                 const ExternalPgFace& declared, const PgSurfaceBound& surface,
                 ReferenceArithmetic& math, m::Work& work, const Diagnostic& diagnostic,
                 NastranPgFaceBound& row) {
    std::array ids{NastranGridId{0}, NastranGridId{0}, NastranGridId{0}, NastranGridId{0}};
    std::visit(
        [&](const auto& corners) {
            if (declared.corners.size() != corners.size()) {
                diagnostic.reject(Code::InvalidInput, Reason::Selection, declared.element);
            }
            row.count = static_cast<std::uint8_t>(corners.size());
            for (std::size_t i = 0; i < corners.size(); ++i) {
                const auto& given = declared.corners[i];
                ids[i] = corners[i].id;
                if (given.node != ids[i].value()) {
                    diagnostic.reject(Code::InvalidInput, Reason::SourceOrder, declared.element, {},
                                      ids[i].value(), corners[i].field.position.byte);
                }
                row.corners[i].node = given.node;
                row.corners[i].chart = copy(math, given.point);
            }
            if constexpr (std::tuple_size_v<std::remove_cvref_t<decltype(corners)>> == 4) {
                // A common Cartesian frame is rigid. Only its raw determinant can
                // certify exact coplanarity in this path, even if CP is irrational.
                prove_planarity(math, diagnostic, deck, corners, shell);
            }
        },
        shell.corners);
    std::optional<NastranPositionDiagnostic> position_diagnostic;
    const auto positions = NastranPositionEnclosures::prepare(
        deck_owner, units, std::span{ids}.first(row.count),
        {math.leaf_limits, row.count, work.limit - work.used, work.stop}, math.plan,
        &position_diagnostic);
    if (!positions) {
        if (position_diagnostic) {
            const auto& error = *position_diagnostic;
            diagnostic.reject(positions.error().code, Reason::Position, row.element, {},
                              error.grid ? std::optional{error.grid->value()} : std::nullopt,
                              error.position ? std::optional{error.position->byte} : std::nullopt,
                              error.coordinate);
        }
        throw geom::PartitionAbort{positions.error()};
    }
    const auto position_view = positions->view();
    if (!position_view) {
        diagnostic.reject(Code::InvalidInput, Reason::Source, row.element);
    }
    const auto source = position_view->grids;
    work.used += position_view->work;
    std::array<IntervalPosition, 4> corners;
    for (std::size_t i = 0; i < row.count; ++i) {
        corners[i] = to_interval(source[i], work);
        row.corners[i].cp = source[i].cp;
        row.corners[i].cp_field = source[i].cp_field;
        row.corners[i].tuple_fields = source[i].tuple_fields;
    }

    const auto raw_normal = directed::cross(directed::subtract(corners[1], corners[0], work),
                                            directed::subtract(corners[2], corners[0], work), work);
    require_positive(directed::norm_squared(raw_normal, work), diagnostic, row);
    const auto normal = directed::normalized(raw_normal, work);
    if (row.count == 4) {
        // The four signed corner Jacobians suffice for this planar bilinear ring.
        for (std::size_t i = 1; i < corners.size(); ++i) {
            const auto turn = directed::cross(
                directed::subtract(corners[(i + 1) % corners.size()], corners[i], work),
                directed::subtract(corners[(i + 2) % corners.size()], corners[i], work), work);
            require_positive(dot(turn, normal, work), diagnostic, row);
        }
    }
    IntervalPosition oriented;
    for (std::size_t i = 0; i < oriented.size(); ++i) {
        oriented[i] = row.winding == ExternalPgWinding::Forward ? normal[i]
                                                                : m::subtract({}, normal[i], work);
    }
    const auto s = surface.curved ? coordinate_range(row.ring(), &Point::s, work) : m::Interval{};
    const auto t = surface.curved ? coordinate_range(row.ring(), &Point::t, work) : m::Interval{};
    const auto pg_normal =
        surface.curved ? curved_pg_normal(surface, s, t, work) : surface.affine_normal;
    require_positive(dot(oriented, pg_normal, work), diagnostic, row);
    const auto normal_upper = upper_norm(directed::subtract(oriented, pg_normal, work), work);
    row.normal_chord_upper = directed::retain(normal_upper, math.leaf_limits, math.plan, work);

    m::Interval offset;
    if (const auto& field = shell.offset) {
        auto metres = math.mul(decimal(math, field->value.decimal()), unit);
        offset = directed::number(metres, work);
        row.nastran_offset = NastranReferenceOffset{std::move(metres), field->field.position.byte};
    }
    IntervalPosition shift;
    for (std::size_t k = 0; k < shift.size(); ++k) {
        shift[k] = m::multiply(offset, normal[k], work);
    }
    row.distance_upper_metres =
        distance_bound(surface, row.ring(), corners, shift, s, t, math, work);
}

void validate_bound_correspondence(ReferenceArithmetic& math, const Diagnostic& diagnostic,
                                   const auto& row, const auto previous,
                                   std::optional<std::uint64_t> element_byte = {}) {
    for (const auto& other : previous) {
        math.visit(1);
        if (row.element == other.element) {
            diagnostic.reject(Code::InvalidInput, Reason::Selection, row.element, {}, {},
                              element_byte);
        }
        for (const auto& corner : row.ring()) {
            for (const auto& old : other.ring()) {
                const bool same_chart = math.compare(corner.chart.s, old.chart.s) == 0 &&
                                        math.compare(corner.chart.t, old.chart.t) == 0;
                if ((corner.node == old.node) != same_chart) {
                    diagnostic.reject(Code::GeometryInvalid, Reason::SourceOrder, row.element,
                                      other.element, corner.node, element_byte);
                }
            }
        }
    }
}

[[nodiscard]] bool bounded_pg_is_curved(ReferenceArithmetic& math, const Diagnostic& diagnostic,
                                        const structure::PgGraphView& graph) {
    bool curved = false;
    for (const auto& term : graph.height.terms()) {
        math.visit(1);
        const auto degree = term.powers.first + term.powers.second;
        if (degree > 4) {
            diagnostic.reject(Code::Unsupported, Reason::CurvedPatch);
        }
        curved = curved || degree > 1;
    }
    return curved;
}

[[nodiscard]] Rational source_si_scale(ReferenceArithmetic& math, SiLengthUnit units) {
    const SiUnits system{units, SiMassUnit::Kilogram, SiTimeUnit::Second, SiAngleUnit::Radian};
    const SiRecipeView scale{SiDecimalView{1, "1", 0}, {}, system, {{1, 0, 0, 0}}, 0};
    math.charge();
    return partition_take(materialize_rational(scale, math.leaf_limits, math.plan));
}

[[nodiscard]] IntervalPosition source_normal(const GmshPgFaceBound& face,
                                             const IntervalPosition& pg_normal,
                                             ReferenceArithmetic& math, m::Work& work,
                                             const Diagnostic& diagnostic) {
    IntervalPosition hull;
    m::Bound projected_lower;
    // For a bilinear map the cross product of its derivatives is affine.
    // At each parameter corner, next-current crossed with previous-current
    // gives that cross product with the original ring orientation.
    const std::size_t count = face.count;
    const auto checks = count == 3 ? std::size_t{1} : count;
    for (std::size_t i = 0; i < checks; ++i) {
        const auto& here = face.corners[i];
        const auto next =
            vector_sub(math, face.corners[(i + 1) % count].position_metres, here.position_metres);
        const auto previous = vector_sub(
            math, face.corners[(i + count - 1) % count].position_metres, here.position_metres);
        auto normal = to_interval(vector_cross(math, next, previous), work);
        if (face.winding == ExternalPgWinding::Reversed) {
            for (auto& component : normal) {
                component = m::subtract({}, component, work);
            }
        }
        const auto projected = dot(normal, pg_normal, work);
        if (projected.lower.mantissa <= 0) {
            diagnostic.reject(
                projected.upper.mantissa <= 0 ? Code::GeometryInvalid : Code::PrecisionUnresolved,
                Reason::SourceWinding, face.element, {}, here.node, here.source_bytes.front());
        }
        if (i == 0) {
            hull = std::move(normal);
            projected_lower = projected.lower;
        } else {
            for (std::size_t k = 0; k < hull.size(); ++k) {
                if (m::compare(normal[k].lower, hull[k].lower) < 0) {
                    hull[k].lower = normal[k].lower;
                }
                if (m::compare(normal[k].upper, hull[k].upper) > 0) {
                    hull[k].upper = normal[k].upper;
                }
            }
            if (m::compare(projected.lower, projected_lower) < 0) {
                projected_lower = projected.lower;
            }
        }
    }
    // Every PG normal being compared has unit length. Its positive projection
    // bounds the source cross-product length below on the entire parameter cell.
    const auto squared = directed::norm_squared(hull, work);
    auto length = m::square_root(
        {squared.lower.mantissa > 0 ? squared.lower : m::Bound{}, squared.upper}, work);
    if (m::compare(projected_lower, length.lower) > 0) {
        length.lower = std::move(projected_lower);
    }
    return {m::divide(hull[0], length, work), m::divide(hull[1], length, work),
            m::divide(hull[2], length, work)};
}

void store_gmsh_corners(const GmshMeshView& mesh, const GmshElement& element, const Rational& unit,
                        const ExternalPgFace& declared, ReferenceArithmetic& math,
                        const Diagnostic& diagnostic, GmshPgFaceBound& row) {
    const auto ordered = element.connectivity();
    if (element.tag != declared.element || ordered.size() != declared.corners.size()) {
        diagnostic.reject(Code::InvalidInput, Reason::Selection, declared.element);
    }
    row.element = element.tag;
    row.entity = element.entity;
    row.winding = declared.winding;
    row.count = static_cast<std::uint8_t>(ordered.size());
    for (std::size_t i = 0; i < ordered.size(); ++i) {
        const auto& node = find_gmsh_node(math, diagnostic, mesh, element.tag, ordered[i]);
        const auto& corner = declared.corners[i];
        if (corner.node != node.tag) {
            diagnostic.reject(Code::InvalidInput, Reason::SourceOrder, element.tag, {}, node.tag,
                              node.coordinates.front().offset);
        }
        auto& out = row.corners[i];
        out.node = node.tag;
        out.chart = copy(math, corner.point);
        out.position_metres = gmsh_coordinates(math, mesh, node, unit);
        for (std::size_t k = 0; k < out.source_bytes.size(); ++k) {
            out.source_bytes[k] = node.coordinates[k].offset;
        }
    }
}

void measure_gmsh_face(const PgSurfaceBound& surface, ReferenceArithmetic& math, m::Work& work,
                       const Diagnostic& diagnostic, GmshPgFaceBound& row) {
    const auto s = surface.curved ? coordinate_range(row.ring(), &Point::s, work) : m::Interval{};
    const auto t = surface.curved ? coordinate_range(row.ring(), &Point::t, work) : m::Interval{};
    const auto pg_normal =
        surface.curved ? curved_pg_normal(surface, s, t, work) : surface.affine_normal;
    const auto normal = source_normal(row, pg_normal, math, work, diagnostic);
    const auto normal_upper = upper_norm(directed::subtract(normal, pg_normal, work), work);
    row.normal_chord_upper = directed::retain(normal_upper, math.leaf_limits, math.plan, work);

    std::array<IntervalPosition, 4> positions;
    for (std::size_t i = 0; i < row.count; ++i) {
        positions[i] = to_interval(row.corners[i].position_metres, work);
    }
    row.distance_upper_metres =
        distance_bound(surface, row.ring(), positions, {}, s, t, math, work);
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
            match(pg, units, charts, limits, plan, diagnostic, gmsh_rationals_per_face,
                  [&](ReferenceArithmetic& math, std::size_t i, const ExternalPgFace&,
                      const Rational& unit) {
                      return gmsh_face(math, diagnostic, *input, selected->faces[i].element, unit);
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
        auto result = match(pg, units, charts, limits, plan, diagnostic, nastran_rationals_per_face,
                            [&](ReferenceArithmetic& math, std::size_t, const ExternalPgFace& chart,
                                const Rational& unit) {
                                return nastran_face(math, diagnostic, *input, chart.element, unit);
                            });
        return ExternalPlanarPg{input->input.bulk.location,
                                input->input.bulk.digest,
                                units,
                                result.patch,
                                std::move(result.faces),
                                result.work};
    });
}

NastranPlanarPgBounds::NastranPlanarPgBounds(RootSnapshotLocation source, SourceDigest digest,
                                             SiLengthUnit units,
                                             structure::PgReferenceIdentity patch,
                                             core::PlannedBuffer<NastranPgFaceBound>&& faces,
                                             structure::PgReferenceCallInfo chart_work,
                                             std::uint64_t directed_work) noexcept
    : source_(source),
      digest_(digest),
      units_(units),
      patch_(patch),
      faces_(std::move(faces)),
      chart_work_(chart_work),
      directed_work_(directed_work) {}

std::optional<NastranPlanarPgBoundView> NastranPlanarPgBounds::view() const& noexcept {
    if (faces_.size() == 0) {
        return std::nullopt;
    }
    return NastranPlanarPgBoundView{source_,         digest_,     units_,        patch_,
                                    faces_.values(), chart_work_, directed_work_};
}

core::Result<NastranPlanarPgBounds> NastranPlanarPgBounds::prepare(
    const structure::PgReferenceGeometry& pg, const NastranShellCards& deck, SiLengthUnit units,
    std::span<const ExternalPgFace> charts, const ExternalPgBoundLimits& limits,
    core::ResourcePlan& plan, std::optional<ExternalPgDiagnostic>* slot) noexcept {
    if (slot) {
        slot->reset();
    }
    const Diagnostic diagnostic{slot};
    return geom::partition_entry([&]() -> core::Result<NastranPlanarPgBounds> {
        try {
            ReferenceArithmetic math{limits.region.arithmetic, plan};
            math.cancel();
            const auto input = deck.view();
            if (!input) {
                diagnostic.reject(Code::InvalidInput, Reason::Source);
            }
            if (charts.empty()) {
                diagnostic.reject(Code::InvalidInput, Reason::Selection);
            }
            if (charts.size() > limits.region.maximum_faces) {
                partition_fail(Code::ResourceExceeded);
            }
            const auto owner = partition_take(pg.source());
            const auto graph = partition_take(owner.graph.view());
            const bool curved = bounded_pg_is_curved(math, diagnostic, graph);
            const auto patch = partition_take(pg.identity());
            const auto triangles = partition_take(owner.partition.triangles());
            const auto vertices = partition_take(owner.partition.vertices());
            // Four chart points, two bounds and one optional offset per face.
            const auto slots = geom::partition_add(geom::reference_live_rationals,
                                                   geom::partition_mul(11, charts.size()));
            math.admit(slots);
            const auto frame = partition_take(pg.frame());
            const auto affine = mid_surface(math, graph, frame);
            m::Work work{0, limits.maximum_directed_work, limits.region.arithmetic.cancellation};
            const auto pq = vector_cross(math, affine.tangent_s, affine.tangent_t);
            const auto affine_normal = directed::normalized(to_interval(pq, work), work);
            const PgSurfaceBound surface{affine, graph, frame, affine_normal, curved};
            const auto unit = source_si_scale(math, units);
            auto faces = partition_take(
                plan.allocate<NastranPgFaceBound>(geom::partition_index(charts.size()), operation));
            Rational total;
            for (std::size_t i = 0; i < charts.size(); ++i) {
                math.visit(1);
                const auto& chart = charts[i];
                const auto& shell = find_shell(math, diagnostic, *input, chart.element);
                auto& row = faces.values()[i];
                row.element = chart.element;
                row.element_field = shell.id_field.position;
                row.winding = chart.winding;
                if (row.winding != ExternalPgWinding::Forward &&
                    row.winding != ExternalPgWinding::Reversed) {
                    diagnostic.reject(Code::InvalidInput, Reason::Winding, row.element);
                }
                fill_bounds(deck, *input, shell, units, unit, chart, surface, math, work,
                            diagnostic, row);
                const auto previous = faces.values().first(i);
                validate_bound_correspondence(math, diagnostic, row, previous,
                                              row.element_field.byte);
                validate_chart_row(math, diagnostic, row, previous, triangles, vertices, total);
            }
            validate_chart_coverage(math, diagnostic, triangles, total);
            math.cancel();
            return NastranPlanarPgBounds{input->input.bulk.location,
                                         input->input.bulk.digest,
                                         units,
                                         patch,
                                         std::move(faces),
                                         math.info(),
                                         work.used};
        } catch (const m::ArithmeticFailure& failure) {
            partition_fail(failure.code);
        } catch (const std::overflow_error&) {
            partition_fail(Code::PrecisionUnresolved);
        } catch (const std::range_error&) {
            partition_fail(Code::PrecisionUnresolved);
        }
    });
}
GmshPgBounds::GmshPgBounds(RootSnapshotLocation source, SourceDigest digest, SiLengthUnit units,
                           structure::PgReferenceIdentity patch,
                           core::PlannedBuffer<GmshPgFaceBound>&& faces,
                           structure::PgReferenceCallInfo chart_work,
                           std::uint64_t directed_work) noexcept
    : source_(source),
      digest_(digest),
      units_(units),
      patch_(patch),
      faces_(std::move(faces)),
      chart_work_(chart_work),
      directed_work_(directed_work) {}

std::optional<GmshPgBoundView> GmshPgBounds::view() const& noexcept {
    if (faces_.size() == 0) {
        return std::nullopt;
    }
    return GmshPgBoundView{source_,         digest_,     units_,        patch_,
                           faces_.values(), chart_work_, directed_work_};
}

core::Result<GmshPgBounds>
GmshPgBounds::prepare(const structure::PgReferenceGeometry& pg, const GmshMesh& mesh,
                      const GmshSurfaceTopology& topology, SiLengthUnit units,
                      std::span<const ExternalPgFace> charts, const ExternalPgBoundLimits& limits,
                      core::ResourcePlan& plan,
                      std::optional<ExternalPgDiagnostic>* slot) noexcept {
    if (slot) {
        slot->reset();
    }
    const Diagnostic diagnostic{slot};
    return geom::partition_entry([&]() -> core::Result<GmshPgBounds> {
        try {
            ReferenceArithmetic math{limits.region.arithmetic, plan};
            math.cancel();
            const auto input = mesh.view();
            const auto selected = topology.view();
            if (!input || !selected) {
                diagnostic.reject(Code::InvalidInput, Reason::Source);
            }
            if (!equal_source(input->location, selected->source.location) ||
                input->digest != selected->source.digest) {
                diagnostic.reject(Code::StaleGeneration, Reason::Source);
            }
            if (charts.empty() || selected->faces.size() != charts.size()) {
                diagnostic.reject(Code::InvalidInput, Reason::Selection);
            }
            if (charts.size() > limits.region.maximum_faces) {
                partition_fail(Code::ResourceExceeded);
            }
            const auto owner = partition_take(pg.source());
            const auto graph = partition_take(owner.graph.view());
            const bool curved = bounded_pg_is_curved(math, diagnostic, graph);
            const auto patch = partition_take(pg.identity());
            const auto triangles = partition_take(owner.partition.triangles());
            const auto vertices = partition_take(owner.partition.vertices());
            // Own all source positions, chart points and two retained bounds per face.
            const auto slots = geom::partition_add(
                geom::reference_live_rationals,
                geom::partition_mul(gmsh_rationals_per_face + 2, charts.size()));
            math.admit(slots);
            const auto frame = partition_take(pg.frame());
            const auto affine = mid_surface(math, graph, frame);
            m::Work work{0, limits.maximum_directed_work, limits.region.arithmetic.cancellation};
            const auto pq = vector_cross(math, affine.tangent_s, affine.tangent_t);
            const auto affine_normal = directed::normalized(to_interval(pq, work), work);
            const PgSurfaceBound surface{affine, graph, frame, affine_normal, curved};
            const auto unit = source_si_scale(math, units);
            auto faces = partition_take(
                plan.allocate<GmshPgFaceBound>(geom::partition_index(charts.size()), operation));
            Rational total;
            for (std::size_t i = 0; i < charts.size(); ++i) {
                math.visit(1);
                const auto& element = selected->faces[i].element;
                auto& row = faces.values()[i];
                store_gmsh_corners(*input, element, unit, charts[i], math, diagnostic, row);
                const auto previous = faces.values().first(i);
                validate_bound_correspondence(math, diagnostic, row, previous);
                validate_chart_row(math, diagnostic, row, previous, triangles, vertices, total);
                measure_gmsh_face(surface, math, work, diagnostic, row);
            }
            validate_chart_coverage(math, diagnostic, triangles, total);
            math.cancel();
            return GmshPgBounds{input->location,  input->digest, units,    patch,
                                std::move(faces), math.info(),   work.used};
        } catch (const m::ArithmeticFailure& failure) {
            partition_fail(failure.code);
        } catch (const std::overflow_error&) {
            partition_fail(Code::PrecisionUnresolved);
        } catch (const std::range_error&) {
            partition_fail(Code::PrecisionUnresolved);
        }
    });
}

} // namespace uvlas::io::detail
