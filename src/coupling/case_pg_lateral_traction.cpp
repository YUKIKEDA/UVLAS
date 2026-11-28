#include "case_pg_lateral_traction.hpp"

#include "io/si_case_detail.hpp"
#include "structure/pg_material_grid.hpp"
#include "structure/pg_reference_domain.hpp"
#include "structure/pg_reference_projection.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace uvlas::coupling::detail {
namespace {
namespace pg = structure;
namespace s = io::schema;
using core::FailureCode;
using namespace pg::detail;
using R = pg::ExactRational;
using Bound = pg::PgMaterialEnclosure;
constexpr auto operation = core::Operation::PrepareModel;

void limit(std::uint64_t count, std::uint64_t maximum) {
    if (count > maximum)
        partition_fail(FailureCode::ResourceExceeded);
}

template <class T>
auto allocate(ReferenceArithmetic& math, std::uint64_t count) {
    math.visit(count);
    return partition_take(math.plan.allocate<T>(partition_index(count), operation));
}

struct Selection {
    std::uint64_t id{};
    std::uint64_t boundary{};
    core::Index first{}, count{};
};

struct Atom {
    core::Index segment{}, edge{};
    bool canonical_direction{};
    R lower, upper;
};

struct BoundaryView {
    std::span<const pg::PgPartitionPoint> vertices, original;
    std::span<const pg::PgPartitionEdge> edges;
    std::span<const pg::PgPartitionEdgeTrace> traces;
    std::span<const pg::PgPartitionSource> sources;
    std::span<const pg::PgPartitionBoundaryRecord> names;
    std::span<const pg::PgPartitionLateralFace> faces;
    std::span<const pg::PgPartitionOwnedLayer> layers;
};

BoundaryView boundary_view(const pg::PgPreparedPatch& patch) {
    const auto basis = partition_take(patch.conditioning().reference_basis());
    const auto& partition = partition_take(basis->source()).partition;
    return {partition_take(partition.vertices()),      partition_take(partition.source_points()),
            partition_take(partition.edges()),         partition_take(partition.edge_traces()),
            partition_take(partition.sources()),       partition_take(partition.boundaries()),
            partition_take(partition.lateral_faces()), partition_take(partition.layers())};
}

template <class Visit>
void visit_loads(const s::Shell& shell, const s::AnalysisStatic& analysis,
                 ReferenceArithmetic& math, Visit visit) {
    const auto visit_list = [&](const auto& list, bool member) {
        for (core::Index i = 0; std::cmp_less(i, list.size()); ++i) {
            math.visit(1);
            const auto entry = partition_take(list.at(i)).record();
            const auto* traction = std::get_if<s::LoadTraction>(&entry);
            if (traction == nullptr)
                continue;
            const auto target = traction->region.record().target.record();
            const bool selected =
                target.kind.selected() == decltype(target.kind.selected())::ItemShell &&
                target.id.unsigned_value() == shell.id.unsigned_value();
            if (member || selected)
                visit(*traction, selected);
        }
    };
    visit_list(shell.loads, true);
    visit_list(analysis.loads, false);
}

template <class Field>
R read_quantity(const io::SiCase& source, Field value, io::SchemaDimension units,
                ReferenceArithmetic& math) {
    math.charge();
    const auto recipe = partition_take(source.quantity(value)).recipe;
    if (recipe.dimensions.exponents != units.exponents || recipe.pi_power != 0)
        partition_fail(FailureCode::InvalidInput);
    return partition_take(io::materialize_rational(recipe, math.leaf_limits, math.plan));
}

std::array<R, 3> world_traction(const io::SiCase& source, const s::LoadTraction& traction,
                                ReferenceArithmetic& math) {
    const auto spatial = traction.spatial.record();
    const auto time = traction.multiplier.record();
    const auto* vector = std::get_if<s::ReferenceFieldConstant>(&spatial);
    const auto* multiplier = std::get_if<s::TimeSeriesConstant>(&time);
    if (traction.area.selected() != decltype(traction.area.selected())::ItemReference ||
        traction.frame.selected() != decltype(traction.frame.selected())::ItemWorld ||
        vector == nullptr || multiplier == nullptr)
        partition_fail(FailureCode::Unsupported);
    if (vector->value.size() != 3 || multiplier->value.size() != 1)
        partition_fail(FailureCode::InvalidInput);
    constexpr io::SchemaDimension pressure{{-1, 1, -2, 0}}, dimensionless{{0, 0, 0, 0}};
    const auto scale =
        read_quantity(source, partition_take(multiplier->value.at(0)), dimensionless, math);
    std::array<R, 3> result;
    for (core::Index i = 0; i < 3; ++i)
        result[static_cast<std::size_t>(i)] = math.mul(
            scale, read_quantity(source, partition_take(vector->value.at(i)), pressure, math));
    return result;
}

void validate_rule(std::span<const pg::PgLineQuadratureNode> nodes, const R& one,
                   ReferenceArithmetic& math) {
    if (nodes.empty())
        partition_fail(FailureCode::InvalidInput);
    R sum;
    const R* previous = nullptr;
    for (const auto& node : nodes) {
        math.visit(1);
        if (node.coordinate.sign() < 0 || math.compare(node.coordinate, one) > 0 ||
            node.weight.sign() <= 0 ||
            (previous != nullptr && math.compare(*previous, node.coordinate) >= 0))
            partition_fail(FailureCode::InvalidInput);
        sum = math.add(sum, node.weight);
        previous = &node.coordinate;
    }
    if (math.compare(sum, one) != 0)
        partition_fail(FailureCode::InvalidInput);
}

const pg::PgPartitionNamedSource* named(const BoundaryView& view,
                                        const pg::PgPartitionEdgeTrace& trace, core::Index id) {
    const auto& source = reference_at(view.sources, trace.source);
    const auto* name = std::get_if<pg::PgPartitionNamedSource>(&source.descriptor);
    return name != nullptr && name->boundary == id ? name : nullptr;
}

std::uint64_t count_atoms(const BoundaryView& view, core::Index id, ReferenceArithmetic& math) {
    std::uint64_t count = 0;
    for (const auto& trace : view.traces) {
        math.visit(1);
        if (named(view, trace, id) != nullptr)
            ++count;
    }
    if (count == 0)
        partition_fail(FailureCode::GeometryInvalid);
    return count;
}

Atom locate_atom(const BoundaryView& view, const pg::PgPartitionBoundaryRecord& record,
                 const pg::PgPartitionNamedSource& name, const pg::PgPartitionEdgeTrace& trace,
                 ReferenceArithmetic& math) {
    if (name.segment < 0 || name.segment >= record.point_count - 1)
        partition_fail(FailureCode::GeometryInvalid);
    const auto& edge = reference_at(view.edges, trace.edge);
    if (!edge.external || edge.left.has_value() == edge.right.has_value())
        partition_fail(FailureCode::GeometryInvalid);
    const auto& start = reference_at(view.original, record.first_point + name.segment);
    const auto& stop = reference_at(view.original, record.first_point + name.segment + 1);
    const auto& first = reference_at(view.vertices, edge.vertices[0]);
    const auto& last = reference_at(view.vertices, edge.vertices[1]);
    const bool use_s = math.compare(start.s, stop.s) != 0;
    const auto coordinate = [&](const auto& point) -> const R& {
        return use_s ? point.s : point.t;
    };
    const auto span = math.sub(coordinate(stop), coordinate(start));
    if (span.sign() == 0)
        partition_fail(FailureCode::GeometryInvalid);
    auto x = math.div(math.sub(coordinate(first), coordinate(start)), span);
    auto y = math.div(math.sub(coordinate(last), coordinate(start)), span);
    const auto orientation = math.compare(x, y);
    if (orientation == 0 || (orientation < 0) != trace.canonical_direction)
        partition_fail(FailureCode::GeometryInvalid);
    if (orientation > 0)
        std::swap(x, y);
    return {name.segment, trace.edge, trace.canonical_direction, std::move(x), std::move(y)};
}

void prepare_atoms(const BoundaryView& view, core::Index id, std::span<Atom> atoms,
                   std::span<std::uint8_t> seen, const R& one, ReferenceArithmetic& math,
                   CasePgLateralTractionDiagnostic& where) {
    const auto& record = reference_at(view.names, id);
    std::size_t i = 0;
    for (const auto& trace : view.traces) {
        math.visit(1);
        if (const auto* segment = named(view, trace, id)) {
            where.segment = segment->segment;
            where.atom = trace.edge;
            // This also rejects positive overlap between different original segments in one Load.
            auto& used = seen[static_cast<std::size_t>(trace.edge)];
            if (used != 0)
                partition_fail(FailureCode::GeometryInvalid);
            used = 1;
            atoms[i++] = locate_atom(view, record, *segment, trace, math);
        }
    }
    std::ranges::sort(atoms, [&](const auto& a, const auto& b) {
        math.visit(1);
        if (a.segment != b.segment)
            return a.segment < b.segment;
        const auto order = math.compare(a.lower, b.lower);
        return order == 0 ? a.edge < b.edge : order < 0;
    });
    std::size_t cursor = 0;
    const R zero;
    where.atom.reset();
    for (core::Index segment = 0; segment < record.point_count - 1; ++segment) {
        math.visit(1);
        where.segment = segment;
        const R* endpoint = &zero;
        while (cursor < atoms.size() && atoms[cursor].segment == segment) {
            const auto& atom = atoms[cursor++];
            where.atom = atom.edge;
            if (math.compare(*endpoint, atom.lower) != 0)
                partition_fail(FailureCode::GeometryInvalid);
            endpoint = &atom.upper;
            seen[static_cast<std::size_t>(atom.edge)] = 0;
        }
        if (math.compare(*endpoint, one) != 0)
            partition_fail(FailureCode::GeometryInvalid);
        where.atom.reset();
    }
    if (cursor != atoms.size())
        partition_fail(FailureCode::GeometryInvalid);
    where.segment.reset();
}

std::span<const pg::PgPartitionLateralFace> atom_faces(const BoundaryView& view, const Atom& atom,
                                                       ReferenceArithmetic& math) {
    // The partition stores physical faces in edge order, with contiguous layers.
    const auto first = std::ranges::lower_bound(
        view.faces, atom.edge,
        [&](core::Index a, core::Index b) {
            math.visit(1);
            return a < b;
        },
        &pg::PgPartitionLateralFace::edge);
    const auto offset = static_cast<std::size_t>(first - view.faces.begin());
    if (view.faces.size() - offset < view.layers.size())
        partition_fail(FailureCode::GeometryInvalid);
    const auto end = offset + view.layers.size();
    if (end < view.faces.size() && view.faces[end].edge == atom.edge)
        partition_fail(FailureCode::GeometryInvalid);
    return view.faces.subspan(offset, view.layers.size());
}

pg::PgReferenceLimits remaining(const ReferenceArithmetic& math, std::uint64_t held_rationals) {
    auto limits = math.limits;
    limits.maximum_scalar_operations -= math.scalar_operations;
    limits.maximum_partition_visits -= math.partition_visits;
    math.reserve_nested_rationals(limits.algebra, held_rationals);
    return limits;
}

template <class Comparator>
Bound two_word_bound(const Comparator& comparator, const R& base) {
    auto& math = comparator.math;
    const auto high = reference_binary64(math, pg_round_word(comparator, base));
    const auto residual = math.sub(base, high);
    if (comparator.compare(residual, {}) == 0)
        return {math.copy(high), math.copy(high)};
    const auto low = pg_round_word(comparator, residual);
    const auto low_exact = reference_binary64(math, low);
    auto middle = math.add(high, low_exact);
    const auto sign = comparator.compare(residual, low_exact);
    if (sign == 0)
        return {math.copy(middle), std::move(middle)};
    auto endpoint = math.add(
        high, reference_binary64(
                  math, std::nextafter(low, sign * std::numeric_limits<double>::infinity())));
    return sign < 0 ? Bound{std::move(endpoint), std::move(middle)}
                    : Bound{std::move(middle), std::move(endpoint)};
}

Bound positive_root(const pg::PgReferenceScalar& square, const R& g, FailureCode nonpositive,
                    ReferenceArithmetic& math) {
    const auto sign = square.radical.sign();
    const auto squared_term = sign == 0 ? R{} : math.mul(g, math.square(square.radical));
    const PgProjectionComparator squared{math, squared_term, sign};
    if (squared.compare(square.rational, {}) <= 0)
        partition_fail(nonpositive);
    const PgPositiveRootComparator root{math, squared, square.rational};
    auto result = two_word_bound(root, {});
    if (result.lower.sign() <= 0)
        partition_fail(FailureCode::PrecisionUnresolved);
    return result;
}

Bound scaled_shape(const pg::PgReferenceScalar& shape, const Bound& root_g, const R& weight,
                   ReferenceArithmetic& math) {
    auto rational = math.mul(shape.rational, weight);
    if (shape.radical.sign() == 0)
        return {math.copy(rational), std::move(rational)};
    const auto radical = math.mul(shape.radical, weight);
    const bool negative = radical.sign() < 0;
    return {math.add(rational, math.mul(radical, negative ? root_g.upper : root_g.lower)),
            math.add(rational, math.mul(radical, negative ? root_g.lower : root_g.upper))};
}

void append_force(Bound& sum, const Bound& shape, const Bound& area, const R& pascals,
                  std::uint32_t grid, ReferenceArithmetic& math) {
    if (pascals.sign() == 0)
        return;
    // Area is strictly positive; each extremum chooses its magnitude by the sign of the shape.
    const auto lower = math.mul(shape.lower, shape.lower.sign() < 0 ? area.upper : area.lower);
    const auto upper = math.mul(shape.upper, shape.upper.sign() < 0 ? area.lower : area.upper);
    const bool negative = pascals.sign() < 0;
    const auto lo = math.mul(negative ? upper : lower, pascals);
    const auto hi = math.mul(negative ? lower : upper, pascals);
    sum.lower = math.add(sum.lower, outward_material_grid(lo, grid, GridEnd::Lower, math));
    sum.upper = math.add(sum.upper, outward_material_grid(hi, grid, GridEnd::Upper, math));
}

} // namespace

core::Result<CasePgLateralTraction>
CasePgLateralTraction::prepare(const io::SiCase& source, const CaseStaticPgPointPatch& bound,
                               const CasePgLateralTractionRule& rule,
                               const CasePgLateralTractionLimits& limits, core::ResourcePlan& plan,
                               CasePgLateralTractionDiagnostic* diagnostic) noexcept {
    CasePgLateralTractionDiagnostic where;
    auto result = partition_entry([&]() -> core::Result<CasePgLateralTraction> {
        ReferenceArithmetic math{limits.reference, plan};
        math.cancel();
        const auto& typed = partition_take(source.source_case()).get();
        const auto identity = partition_take(typed.identity());
        if (!io::detail::same_identity(identity, bound.source()))
            partition_fail(FailureCode::StaleGeneration);
        const auto root = partition_take(typed.root());
        const auto analysis = root.analysis.record();
        const auto* statics = std::get_if<s::AnalysisStatic>(&analysis);
        if (statics == nullptr)
            partition_fail(FailureCode::Unsupported);
        const auto model = root.model.record();
        std::optional<s::Shell> shell;
        const auto shell_id = bound.shell().value();
        for (core::Index i = 0; std::cmp_less(i, model.shells.size()); ++i) {
            math.visit(1);
            const auto member = partition_take(model.shells.at(i)).record();
            if (member.id.unsigned_value() == shell_id) {
                shell.emplace(member);
                break;
            }
        }
        if (!shell)
            partition_fail(FailureCode::InvalidInput);
        const auto info = partition_take(bound.patch().info());
        const auto view = boundary_view(bound.patch());
        limit(partition_add(rule.edge.size(), rule.depth.size()), limits.maximum_rule_nodes);
        limit(view.layers.size(), limits.maximum_layers);
        math.admit();
        const auto one = math.constant("1");
        validate_rule(rule.edge, one, math);
        validate_rule(rule.depth, one, math);

        std::uint64_t load_count = 0;
        visit_loads(*shell, *statics, math, [&](const auto& load, bool) {
            where.load = load.id.unsigned_value();
            limit(++load_count, limits.maximum_loads);
        });
        auto loads = allocate<Load>(math, load_count);
        auto selection = allocate<Selection>(math, load_count);
        math.leaf_limits = limits.reference.algebra;
        math.admit(partition_add(reference_live_rationals, partition_mul(load_count, 3)));
        std::size_t i = 0;
        std::uint64_t atom_count = 0;
        const auto points_per_atom =
            partition_mul(view.layers.size(), partition_mul(rule.edge.size(), rule.depth.size()));
        visit_loads(*shell, *statics, math, [&](const auto& traction, bool selected) {
            const auto region = traction.region.record();
            const auto id = traction.id.unsigned_value(),
                       boundary = region.boundary.unsigned_value();
            where = {.load = id, .boundary = boundary};
            if (!selected || region.geometry.unsigned_value() != shell->geometry.unsigned_value())
                partition_fail(FailureCode::InvalidInput);
            if (region.kind.selected() != decltype(region.kind.selected())::ItemEdge ||
                region.side.selected() != decltype(region.side.selected())::ItemInterior)
                partition_fail(FailureCode::Unsupported);
            if (boundary >= view.names.size())
                partition_fail(FailureCode::InvalidInput);
            const auto count = count_atoms(view, partition_index(boundary), math);
            selection.values()[i] = {id, boundary, partition_index(atom_count),
                                     partition_index(count)};
            atom_count = partition_add(atom_count, count);
            limit(atom_count, limits.maximum_atoms);
            limit(partition_mul(atom_count, points_per_atom), limits.maximum_points);
            loads.values()[i++] = {world_traction(source, traction, math),
                                   partition_index(partition_mul(count, points_per_atom))};
        });
        const auto point_count = partition_mul(atom_count, points_per_atom);
        auto points = allocate<CasePgLateralTractionPoint>(math, point_count);
        if (point_count == 0)
            return CasePgLateralTraction{identity, bound.shell(), info.identity, std::move(loads),
                                         std::move(points)};
        auto atoms = allocate<Atom>(math, atom_count);
        auto seen = allocate<std::uint8_t>(math, view.edges.size());
        math.leaf_limits = limits.reference.algebra;
        const auto held = partition_add(
            partition_mul(3, load_count),
            partition_add(partition_mul(2, atom_count), partition_mul(4, point_count)));
        math.admit(partition_add(reference_live_rationals, held));
        std::size_t cursor = 0;
        for (const auto& load : selection.values()) {
            where = {.load = load.id, .boundary = load.boundary};
            auto group = atoms.values().subspan(static_cast<std::size_t>(load.first),
                                                static_cast<std::size_t>(load.count));
            prepare_atoms(view, partition_index(load.boundary), group, seen.values(), one, math,
                          where);
            for (const auto& atom : group) {
                where.segment = atom.segment;
                where.atom = atom.edge;
                where.face.reset();
                where.layer.reset();
                const auto faces = atom_faces(view, atom, math);
                for (std::size_t j = 0; j < faces.size(); ++j) {
                    math.visit(1);
                    const auto& face = faces[j];
                    const auto face_index =
                        partition_index(static_cast<std::size_t>(&face - view.faces.data()));
                    where.face = face_index;
                    where.layer = view.layers[j].id;
                    if (face.edge != atom.edge || face.layer != partition_index(j))
                        partition_fail(FailureCode::GeometryInvalid);
                    const bool forward = atom.canonical_direction == face.canonical_direction;
                    for (const auto& u : rule.edge)
                        for (const auto& v : rule.depth) {
                            math.visit(1);
                            points.values()[cursor++] = {load.id,
                                                         load.boundary,
                                                         atom.segment,
                                                         atom.edge,
                                                         face_index,
                                                         view.layers[j].id,
                                                         math.copy(u.coordinate),
                                                         forward ? math.copy(u.coordinate)
                                                                 : math.sub(one, u.coordinate),
                                                         math.copy(v.coordinate),
                                                         math.mul(u.weight, v.weight)};
                        }
                }
            }
        }
        math.cancel();
        return CasePgLateralTraction{identity, bound.shell(), info.identity, std::move(loads),
                                     std::move(points)};
    });
    if (diagnostic != nullptr)
        *diagnostic = result ? CasePgLateralTractionDiagnostic{} : where;
    return result;
}

core::Result<CasePgLateralTractionForce>
CasePgLateralTraction::evaluate(const CaseStaticPgPointPatch& bound, double time_seconds,
                                const CasePgLateralTractionAccuracy& accuracy,
                                const CasePgLateralTractionLimits& limits, core::ResourcePlan& plan,
                                CasePgLateralTractionDiagnostic* diagnostic) const noexcept {
    CasePgLateralTractionDiagnostic where;
    auto result = partition_entry([&]() -> core::Result<CasePgLateralTractionForce> {
        ReferenceArithmetic math{limits.reference, plan};
        math.cancel();
        const auto info = partition_take(bound.patch().info());
        if (!source_ || !io::detail::same_identity(*source_, bound.source()) ||
            shell_ != bound.shell() || patch_ != info.identity)
            partition_fail(FailureCode::StaleGeneration);
        if (!std::isfinite(time_seconds))
            partition_fail(FailureCode::Nonfinite);
        if (accuracy.maximum_force_width_newtons.sign() < 0)
            partition_fail(FailureCode::InvalidInput);
        if (accuracy.sum_fraction_bits >= 32768)
            partition_fail(FailureCode::PrecisionUnsupported);
        limit(loads_.size(), limits.maximum_loads);
        limit(points_.size(), limits.maximum_points);
        const auto retained =
            partition_add(partition_mul(loads_.size(), 3), partition_mul(points_.size(), 4));
        const auto assembly = partition_add(retained, partition_mul(info.basis.rank, 6));
        const auto nested_payload = partition_add(assembly, reference_live_rationals);
        const auto batch =
            partition_add(reference_live_rationals, partition_mul(info.basis.rank, 8));
        math.admit(partition_add(assembly, batch));
        auto force = allocate<pg::PgShellCartesianEnclosure>(math, info.basis.rank);
        std::size_t cursor = 0;
        for (const auto& load : loads_.values()) {
            for (core::Index i = 0; i < load.point_count; ++i) {
                const auto& point = points_.values()[cursor++];
                where = {point.load, point.boundary, point.segment,
                         point.atom, point.face,     point.layer};
                math.visit(1);
                const pg::PgReferencePoint query =
                    pg::PgReferenceLateralPoint{point.face, point.face_along, point.depth};
                auto sample = partition_take(bound.patch().evaluate_reference(
                    std::span{&query, 1}, remaining(math, nested_payload), plan));
                math.visit(sample.call_info.partition_visits);
                math.charge(sample.call_info.scalar_operations);
                const auto& reference = sample.points.values().front();
                const auto* lateral =
                    std::get_if<pg::PgReferenceLateralMeasure>(&reference.measure);
                if (lateral == nullptr || reference.layer != point.layer ||
                    sample.column_count != info.basis.rank)
                    partition_fail(FailureCode::StaleGeneration);
                const auto area = positive_root(lateral->area.squared_norm, reference.g,
                                                FailureCode::GeometryInvalid, math);
                const auto columns = sample.columns.values();
                const auto root_g = std::ranges::any_of(columns,
                                                        [&](const auto& entry) {
                                                            math.visit(1);
                                                            return entry.value.radical.sign() != 0;
                                                        })
                                        ? positive_root({math.copy(reference.g), {}}, reference.g,
                                                        FailureCode::GeometryInvalid, math)
                                        : Bound{};
                for (std::size_t column = 0; column < columns.size(); ++column) {
                    math.visit(1);
                    const auto shape =
                        scaled_shape(columns[column].value, root_g, point.weight, math);
                    for (std::size_t axis = 0; axis < 3; ++axis)
                        append_force(force.values()[column][axis], shape, area, load.pascals[axis],
                                     accuracy.sum_fraction_bits, math);
                }
            }
        }
        where = {};
        for (const auto& row : force.values()) {
            math.visit(1);
            for (const auto& axis : row)
                if (math.compare(math.sub(axis.upper, axis.lower),
                                 accuracy.maximum_force_width_newtons) > 0)
                    partition_fail(FailureCode::PrecisionUnresolved);
        }
        math.cancel();
        return CasePgLateralTractionForce{patch_, shell_, std::move(force), math.info()};
    });
    if (diagnostic != nullptr)
        *diagnostic = result ? CasePgLateralTractionDiagnostic{} : where;
    return result;
}

} // namespace uvlas::coupling::detail
