#include "case_nastran_static_point_support.hpp"

#include "case_material_reference.hpp"
#include "io/nastran_support_selection.hpp"
#include "io/si_case_detail.hpp"
#include "structure/pg_reference_arithmetic.hpp"

#include <algorithm>
#include <bit>

namespace uvlas::coupling::detail {
namespace {
namespace in = io::detail;
namespace ref = structure::detail;
namespace s = io::schema;
using ref::partition_index;
using ref::partition_take;
using Code = core::FailureCode;
using Reason = CaseNastranPointSupportFailure;
using Rational = structure::ExactRational;

struct Query {
    ref::ReferenceArithmetic math;
    std::optional<CaseNastranPointSupportDiagnostic>* diagnostic{};
    CaseNastranPointSupportDiagnostic detail{};

    [[noreturn]] void reject(Reason reason, Code code = Code::InvalidInput) {
        detail.reason = reason;
        if (diagnostic)
            *diagnostic = detail;
        ref::partition_fail(code);
    }

    [[nodiscard]] Rational exact(const io::SiRecipeView& recipe, Reason reason) {
        if (recipe.pi_power != 0)
            reject(reason, Code::Unsupported);
        math.charge();
        return partition_take(io::materialize_rational(recipe, math.leaf_limits, math.plan));
    }

    [[nodiscard]] bool equal(const Rational& a, const Rational& b) {
        const auto words = [](std::uint64_t bits) { return bits / 32 + (bits % 32 != 0 ? 1 : 0); };
        const auto size = [&](const Rational& value) {
            return words(value.numerator_bits()) + words(value.denominator_bits());
        };
        math.charge(1 + size(a) + size(b));
        return a == b;
    }
};

[[nodiscard]] std::optional<in::NastranBulkPosition>
position(const std::optional<in::NastranBulkLexeme>& field) {
    return field ? std::optional{field->position} : std::nullopt;
}

[[nodiscard]] s::SupportConstraint select_support(Query& q, const s::Model& model,
                                                  in::MeshSupportId requested) {
    for (std::size_t i = 0; i < model.supports.size(); ++i) {
        q.math.visit(1);
        const auto declaration = partition_take(model.supports.at(partition_index(i))).record();
        if (const auto* support = std::get_if<s::SupportConstraint>(&declaration);
            support && support->id.unsigned_value() == requested.value())
            return *support;
    }
    q.reject(Reason::Support, Code::Unsupported);
}

[[nodiscard]] const CaseNastranPlanarPgSet::Shell&
select_shell(Query& q, const CaseNastranPlanarPgSet& set, const in::ParsedMeshSidecar& sidecar,
             in::MeshShellId shell) {
    for (std::size_t i = 0; i < set.size(); ++i) {
        q.math.visit(1);
        const auto& prepared = partition_take(set.at(partition_index(i))).get();
        const auto& match = prepared.match().planar;
        if (match.shell == shell) {
            if (!in::same_identity(match.sidecar.parent, sidecar.source().parent) ||
                match.sidecar.source != sidecar.source().source ||
                match.sidecar_digest != sidecar.digest())
                q.reject(Reason::Binding, Code::StaleGeneration);
            return prepared;
        }
    }
    q.reject(Reason::Position, Code::Unsupported);
}

struct PointMatch {
    std::uint64_t element;
    structure::PgPartitionPoint chart;
    structure::PgCartesianExactPosition world_metres;
};

[[nodiscard]] PointMatch compare_point(Query& q, const io::SiCase& si,
                                       const s::MaterialPoint& target,
                                       const in::ExternalPlanarPgView& geometry) {
    const auto coordinate = [&](core::Index i) {
        return q.exact(partition_take(si.quantity(partition_take(target.coordinates.at(i)))).recipe,
                       Reason::Position);
    };
    if (coordinate(2).sign() != 0)
        q.reject(Reason::Position);
    structure::PgPartitionPoint chart{coordinate(0), coordinate(1)};
    for (const auto& face : geometry.faces) {
        for (const auto& corner : face.ring()) {
            q.math.visit(1);
            if (corner.node != q.detail.grid.value() || !q.equal(corner.chart.s, chart.s) ||
                !q.equal(corner.chart.t, chart.t))
                continue;
            q.detail.element = face.element;
            if (corner.nastran && corner.nastran->unshifted_metres &&
                !std::ranges::equal(*corner.nastran->unshifted_metres, corner.position_metres,
                                    [&](const auto& original, const auto& current) {
                                        return q.equal(original, current);
                                    })) {
                q.detail.offset_byte =
                    face.nastran_offset.transform([](const auto& offset) { return offset.byte; });
                q.reject(Reason::Position);
            }
            const auto& xyz = corner.position_metres;
            return {face.element,
                    std::move(chart),
                    {q.math.copy(xyz[0]), q.math.copy(xyz[1]), q.math.copy(xyz[2])}};
        }
    }
    q.reject(Reason::Position);
}

[[nodiscard]] int sign(const io::SiRecipeView& recipe) {
    return std::visit(
        [](const auto& value) {
            if constexpr (requires { value.sign; })
                return value.sign;
            else
                return value.numerator.sign;
        },
        recipe.value);
}

void require_static(Query& q, const io::SiCase& si, const s::Motion& motion) {
    const auto zero = [&](const auto& series) {
        const auto data = series.record();
        const auto* constant = std::get_if<s::TimeSeriesConstant>(&data);
        if (!constant)
            q.reject(Reason::Motion, Code::Unsupported);
        for (core::Index i = 0; i < 3; ++i) {
            q.math.charge();
            const auto recipe =
                partition_take(si.quantity(partition_take(constant->value.at(i)))).recipe;
            if (sign(recipe) != 0)
                q.reject(Reason::Motion, Code::Unsupported);
        }
    };
    zero(motion.translation);
    zero(motion.rotation_vector);
    const auto matrix = motion.origin.record().rotation;
    const auto one = q.math.constant("1");
    for (core::Index i = 0; i < 3; ++i) {
        const auto row = partition_take(matrix.at(i));
        for (core::Index j = 0; j < 3; ++j) {
            const auto field = partition_take(row.at(j));
            const auto entry = q.exact(partition_take(si.quantity(field)).recipe, Reason::Motion);
            if (i == j ? entry != one : entry.sign() != 0)
                q.reject(Reason::Motion, Code::Unsupported);
        }
    }
}

[[nodiscard]] unsigned world_mask(Query& q, const s::SupportConstraint& support) {
    unsigned result = 0;
    for (std::size_t i = 0; i < support.components.size(); ++i) {
        q.math.visit(1);
        const auto component =
            partition_take(support.components.at(partition_index(i))).unsigned_value();
        if (component >= 3)
            q.reject(Reason::Components, Code::Unsupported);
        result |= 1U << component;
    }
    return result;
}

using Fields = std::array<std::optional<in::NastranBulkPosition>, 3>;

void compare_components(Query& q, const in::NastranShellCards& shells, unsigned local,
                        unsigned world, const Fields& fields) {
    if (std::popcount(local) != std::popcount(world))
        q.reject(Reason::Components);
    const auto& limits = q.math.limits;
    // Borrow no support or chart rationals in this nested calculation.
    const in::NastranSpcDirectionLimits directions{
        limits.algebra, limits.maximum_partition_visits - q.math.partition_visits,
        limits.cancellation};
    auto result = in::prepare_nastran_displacement_axes(
        shells, q.detail.grid, static_cast<std::uint8_t>(local), directions, q.math.plan,
        &q.detail.direction);
    if (!result) {
        q.detail.reason = Reason::Components;
        if (q.diagnostic)
            *q.diagnostic = q.detail;
        throw ref::PartitionAbort{result.error()};
    }
    q.math.visit(result->work);
    // There are equally many independent directions on each side. Containment gives equality.
    bool unresolved = false;
    for (std::size_t i = 0; i < 3; ++i) {
        if ((local & (1U << i)) == 0)
            continue;
        for (std::size_t j = 0; j < 3; ++j) {
            q.math.charge();
            if ((world & (1U << j)) != 0)
                continue;
            const auto& interval = result->directions[i][j];
            if (interval.lower.sign() > 0 || interval.upper.sign() < 0) {
                q.detail.field = fields[i];
                q.reject(Reason::Components);
            }
            if (interval.lower.sign() != 0 || interval.upper.sign() != 0) {
                q.detail.field = fields[i];
                unresolved = true;
            }
        }
    }
    if (unresolved)
        q.reject(Reason::Components, Code::PrecisionUnresolved);
}

} // namespace

core::Result<CaseNastranPointSupport> certify_case_nastran_static_point_support(
    const io::SiCase& si, const in::ExternalMeshSidecarBinding& bound,
    const CaseNastranPlanarPgSet& prepared, in::NastranSupportSetId set, in::NastranGridId grid,
    const CaseNastranPointSupportLimits& limits, core::ResourcePlan& plan,
    std::optional<CaseNastranPointSupportDiagnostic>* diagnostic) noexcept {
    if (diagnostic)
        diagnostic->reset();
    return ref::partition_entry([&]() -> core::Result<CaseNastranPointSupport> {
        Query q{{limits.arithmetic, plan}, diagnostic, {}};
        q.detail.set = set;
        q.detail.grid = grid;
        q.math.cancel();
        const auto& schema = partition_take(si.source_case()).get();
        const auto& sidecar = bound.sidecar();
        if (!in::same_identity(partition_take(schema.identity()), sidecar.source().parent))
            q.reject(Reason::Binding, Code::StaleGeneration);
        const auto* input = std::get_if<in::NastranSidecar>(&sidecar.payload());
        const auto* shells = std::get_if<in::NastranShellCards>(&bound.source());
        if (input == nullptr || shells == nullptr)
            q.reject(Reason::Binding, Code::Unsupported);
        const auto mappings = input->supports.values();
        const auto step = [&] { q.math.visit(1); };
        const auto row = std::ranges::find_if(mappings, [&](const auto& candidate) {
            step();
            return candidate.set == set && candidate.grid == grid;
        });
        if (row == mappings.end())
            q.reject(Reason::Support);
        q.detail.declaration = row->node;
        const auto* source = bound.nastran_supports();
        if (source == nullptr)
            q.reject(Reason::Binding, Code::StaleGeneration);
        const auto terms = source->view();
        const auto deck = shells->view();
        if (!terms || !deck ||
            !deck->input.bulk.same_snapshot(terms->source.location, terms->source.digest))
            q.reject(Reason::Binding, Code::StaleGeneration);
        const auto grids = deck->input.grids;
        const auto node = std::ranges::find_if(grids, [&](const auto& candidate) {
            step();
            return candidate.id == grid;
        });
        if (node == grids.end())
            q.reject(Reason::Binding, Code::StaleGeneration);
        q.detail.cd_field = position(node->displacement_system.field);
        q.detail.field = node->id_field.position;

        const auto support =
            select_support(q, partition_take(schema.root()).model.record(), row->support);
        const auto target = support.target.record();
        const auto member = target.target.record();
        if (member.kind.selected() != decltype(member.kind.selected())::ItemShell ||
            target.side.selected() != decltype(target.side.selected())::ItemInterior)
            q.reject(Reason::Position, Code::Unsupported);
        const auto shell = in::MeshShellId{member.id.unsigned_value()};
        const auto& selected = select_shell(q, prepared, sidecar, shell);
        const auto patch = partition_take(selected.reference().identity());
        const auto geometry = selected.match().planar.geometry.view();
        if (!geometry || geometry->patch != patch ||
            !deck->input.bulk.same_snapshot(geometry->source, geometry->digest))
            q.reject(Reason::Binding, Code::StaleGeneration);

        q.math.admit(ref::reference_live_rationals + 5);
        auto point = compare_point(q, si, target, *geometry);
        require_static(q, si, support.motion.record());
        const auto world = world_mask(q, support);
        std::size_t count = 0;
        Fields fields{};
        in::visit_nastran_support_terms(*terms, set, grid, step, [&](const auto& term) {
            if (count == limits.maximum_terms)
                ref::partition_fail(Code::ResourceExceeded);
            q.detail.field = term.component_field.position;
            if ((term.components & ~7U) != 0)
                q.reject(Reason::Components, Code::Unsupported);
            if (term.prescribed_motion) {
                q.math.charge(term.prescribed_motion->decimal().size());
                const auto value = partition_take(Rational::from_decimal(
                    term.prescribed_motion->decimal(), q.math.leaf_limits, plan));
                if (value.sign() != 0) {
                    q.detail.field = position(term.motion_field);
                    q.reject(Reason::Motion);
                }
            }
            for (std::size_t i = 0; i < fields.size(); ++i) {
                if (!fields[i] && (term.components & (1U << i)) != 0)
                    fields[i] = term.component_field.position;
            }
            ++count;
        });
        compare_components(q, *shells, row->component_bits, world, fields);
        auto owned = partition_take(plan.allocate<CaseNastranPointSupportTerm>(
            partition_index(count), ref::partition_operation));
        std::size_t next = 0;
        in::visit_nastran_support_terms(*terms, set, grid, step, [&](const auto& term) {
            owned.values()[next++] = {term.set, term.components, term.card.position,
                                      term.component_field.position, position(term.motion_field)};
        });
        return CaseNastranPointSupport{terms->source,
                                       sidecar.source(),
                                       sidecar.digest(),
                                       *row,
                                       shell,
                                       patch,
                                       point.element,
                                       std::move(point.chart),
                                       std::move(point.world_metres),
                                       node->id_field.position,
                                       node->displacement_system.id,
                                       q.detail.cd_field,
                                       static_cast<std::uint8_t>(world),
                                       std::move(owned)};
    });
}

core::Result<CaseNastranPointSupportRow> assemble_case_nastran_static_point_row(
    const CaseNastranPointSupport& point, const structure::PgPreparedPatch& patch,
    const structure::PgReferenceLimits& limits, core::ResourcePlan& plan) noexcept {
    return ref::partition_entry([&]() -> core::Result<CaseNastranPointSupportRow> {
        namespace pg = structure;
        ref::ReferenceArithmetic math{limits, plan};
        math.cancel();
        const auto info = partition_take(patch.info());
        if (info.geometry != point.patch)
            ref::partition_fail(Code::StaleGeneration);
        if (point.world_components == 0 || (point.world_components & ~7U) != 0)
            ref::partition_fail(Code::InvalidInput);
        const auto* basis = partition_take(patch.conditioning().reference_basis());
        const auto source = partition_take(basis->source());
        const auto vertices = partition_take(source.partition.vertices());
        const auto triangles = partition_take(source.partition.triangles());
        const auto layers = partition_take(source.partition.layers());
        const auto volumes = partition_take(source.partition.volumes());
        const auto coordinates = patch.reference_coordinates();
        const auto columns = static_cast<std::uint64_t>(coordinates.size());
        // One live batch holds N and three physical-gradient scalars per column.
        // Each scalar has two rationals. Also admit the retained row and RHS.
        const auto local_slots = ref::partition_add(ref::basis_live_rationals, 12);
        math.admit(ref::partition_add(local_slots, ref::partition_mul(9, columns)));
        std::optional<CaseNastranPointSupportRow> result;
        for (std::size_t i = 0; i < volumes.size(); ++i) {
            math.visit(1);
            const auto& volume = volumes[i];
            const auto& layer = layers[static_cast<std::size_t>(volume.layer)];
            if (layer.lower.sign() > 0 || layer.upper.sign() < 0)
                continue;
            const auto& triangle = triangles[static_cast<std::size_t>(volume.triangle)];
            const auto uv = case_material_triangle_coordinates(math, triangle, vertices,
                                                               point.chart.s, point.chart.t);
            if (!uv)
                continue;
            const auto u = math.div((*uv)[0], triangle.twice_area);
            const auto v = math.div((*uv)[1], triangle.twice_area);
            const auto depth = math.div(math.neg(layer.lower), math.sub(layer.upper, layer.lower));
            const std::array<pg::PgReferencePoint, 1> query{
                pg::PgReferenceVolumePoint{partition_index(i), u, v, depth}};
            auto nested = limits;
            nested.maximum_scalar_operations -= math.scalar_operations;
            nested.maximum_partition_visits -= math.partition_visits;
            const auto saved = result ? result->terms.size() + point.world_metres.size() : 0;
            math.reserve_nested_rationals(nested.algebra, ref::partition_add(saved, 7));
            auto evaluated = partition_take(patch.evaluate_reference(query, nested, plan));
            math.charge(evaluated.call_info.scalar_operations);
            math.visit(evaluated.call_info.partition_visits);
            auto entries = evaluated.columns.values();
            if (entries.size() != columns)
                ref::partition_fail(Code::InvalidInput);
            if (result) {
                const auto terms = result->terms.values();
                std::size_t next = 0;
                for (std::size_t j = 0; j < entries.size(); ++j) {
                    math.visit(1);
                    const auto& value = entries[j].value;
                    if (value.radical.sign() != 0)
                        ref::partition_fail(Code::PrecisionUnresolved);
                    if (value.rational.sign() == 0)
                        continue;
                    if (next >= terms.size() || terms[next].column != partition_index(j) ||
                        math.compare(terms[next].coefficient, value.rational) != 0)
                        ref::partition_fail(Code::GeometryInvalid);
                    ++next;
                }
                if (next != terms.size())
                    ref::partition_fail(Code::GeometryInvalid);
                continue;
            }
            std::uint64_t nonzeros = 0;
            for (const auto& entry : entries) {
                math.visit(1);
                if (entry.value.radical.sign() != 0)
                    ref::partition_fail(Code::PrecisionUnresolved);
                nonzeros += entry.value.rational.sign() != 0 ? 1U : 0U;
            }
            if (nonzeros == 0)
                ref::partition_fail(Code::GeometryInvalid);
            auto terms = partition_take(plan.allocate<pg::PgTraceSparseEntry>(
                partition_index(nonzeros), ref::partition_operation));
            pg::PgCartesianExactPosition actual;
            std::size_t next = 0;
            for (std::size_t j = 0; j < entries.size(); ++j) {
                math.visit(1);
                auto& value = entries[j].value.rational;
                if (value.sign() == 0)
                    continue;
                for (std::size_t axis = 0; axis < actual.size(); ++axis)
                    actual[axis] = math.add(actual[axis], math.mul(value, coordinates[j][axis]));
                terms.values()[next++] = {partition_index(j), std::move(value)};
            }
            for (std::size_t axis = 0; axis < actual.size(); ++axis)
                if (math.compare(actual[axis], point.world_metres[axis]) != 0)
                    ref::partition_fail(Code::GeometryInvalid);
            result.emplace(CaseNastranPointSupportRow{info.identity, point.mapping.support,
                                                      point.shell, point.world_components,
                                                      std::move(terms), std::move(actual)});
        }
        if (!result)
            ref::partition_fail(Code::GeometryInvalid);
        math.cancel();
        return std::move(*result);
    });
}
} // namespace uvlas::coupling::detail
