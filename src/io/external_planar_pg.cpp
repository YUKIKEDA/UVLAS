#include "external_planar_pg.hpp"

#include "../structure/pg_conditioning_factor.hpp"
#include "../structure/pg_reference_arithmetic.hpp"
#include "nastran_position_arithmetic.hpp"
#include "nastran_position_enclosures.hpp"

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

[[nodiscard]] Face gmsh_face(ReferenceArithmetic& math, const Diagnostic& diagnostic,
                             const GmshMeshView& mesh, const GmshElement& element,
                             const Rational& unit) {
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
        for (std::size_t k = 0; k < out.grid_metres.size(); ++k) {
            out.grid_metres[k] = math.mul(decimal(math, mesh.text(node->coordinates[k])), unit);
            out.offsets[k] = node->coordinates[k].offset;
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

void fill_bounds(const NastranShellCards& deck_owner, const NastranShellView& deck,
                 const NastranShellElement& shell, SiLengthUnit units, const Rational& unit,
                 const ExternalPgFace& declared, const AffineMidSurface& surface,
                 const IntervalPosition& pg_normal, ReferenceArithmetic& math, m::Work& work,
                 const Diagnostic& diagnostic, NastranPgFaceBound& row) {
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
    require_positive(dot(oriented, pg_normal, work), diagnostic, row);
    const auto normal_upper = upper_norm(directed::subtract(oriented, pg_normal, work), work);
    row.normal_chord_upper = directed::retain(normal_upper, math.leaf_limits, math.plan, work);

    m::Interval offset;
    if (const auto& field = shell.offset) {
        auto metres = math.mul(decimal(math, field->value.decimal()), unit);
        offset = directed::number(metres, work);
        row.nastran_offset = NastranReferenceOffset{std::move(metres), field->field.position.byte};
    }
    // PG is affine in its chart, so also its corner interpolant. Both
    // residual maps use nonnegative weights summing to one over the whole face.
    m::Bound distance;
    for (std::size_t i = 0; i < row.count; ++i) {
        IntervalPosition residual;
        for (std::size_t k = 0; k < residual.size(); ++k) {
            const auto actual = m::add(corners[i][k], m::multiply(offset, normal[k], work), work);
            const auto expected = surface.coordinate(math, row.corners[i].chart, k);
            residual[k] = m::subtract(actual, directed::number(expected, work), work);
        }
        auto corner = upper_norm(residual, work);
        if (m::compare(corner, distance) > 0) {
            distance = std::move(corner);
        }
    }
    row.distance_upper_metres = directed::retain(distance, math.leaf_limits, math.plan, work);
}

void validate_bound_correspondence(ReferenceArithmetic& math, const Diagnostic& diagnostic,
                                   const NastranPgFaceBound& row,
                                   std::span<const NastranPgFaceBound> previous) {
    for (const auto& other : previous) {
        math.visit(1);
        if (row.element == other.element) {
            diagnostic.reject(Code::InvalidInput, Reason::Selection, row.element, {}, {},
                              row.element_field.byte);
        }
        for (const auto& corner : row.ring()) {
            for (const auto& old : other.ring()) {
                const bool same_chart = math.compare(corner.chart.s, old.chart.s) == 0 &&
                                        math.compare(corner.chart.t, old.chart.t) == 0;
                if ((corner.node == old.node) != same_chart) {
                    diagnostic.reject(Code::GeometryInvalid, Reason::SourceOrder, row.element,
                                      other.element, corner.node, row.element_field.byte);
                }
            }
        }
    }
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
    std::span<const ExternalPgFace> charts, const NastranPgBoundLimits& limits,
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
            for (const auto& term : graph.height.terms()) {
                math.visit(1);
                if (term.powers.first + term.powers.second > 1) {
                    diagnostic.reject(Code::Unsupported, Reason::CurvedPatch);
                }
            }
            const auto patch = partition_take(pg.identity());
            const auto triangles = partition_take(owner.partition.triangles());
            const auto vertices = partition_take(owner.partition.vertices());
            // Four chart points, two bounds and one optional offset per face.
            const auto slots = geom::partition_add(geom::reference_live_rationals,
                                                   geom::partition_mul(11, charts.size()));
            math.admit(slots);
            const auto surface = mid_surface(math, graph, partition_take(pg.frame()));
            m::Work work{0, limits.maximum_directed_work, limits.region.arithmetic.cancellation};
            const auto pq = vector_cross(math, surface.tangent_s, surface.tangent_t);
            const auto pg_normal = directed::normalized(to_interval(pq, work), work);
            const SiUnits system{units, SiMassUnit::Kilogram, SiTimeUnit::Second,
                                 SiAngleUnit::Radian};
            const SiRecipeView scale{SiDecimalView{1, "1", 0}, {}, system, {{1, 0, 0, 0}}, 0};
            math.charge();
            const auto unit = partition_take(materialize_rational(scale, math.leaf_limits, plan));
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
                fill_bounds(deck, *input, shell, units, unit, chart, surface, pg_normal, math, work,
                            diagnostic, row);
                const auto previous = faces.values().first(i);
                validate_bound_correspondence(math, diagnostic, row, previous);
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
} // namespace uvlas::io::detail
