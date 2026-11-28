#include "case_material_reference.hpp"

#include "structure/pg_reference_arithmetic.hpp"

#include <optional>

namespace uvlas::coupling::detail {
namespace {
namespace pg = structure::detail;
using core::FailureCode;
using pg::partition_fail;
using pg::partition_index;
using pg::partition_take;
using R = structure::ExactRational;
using Point = structure::PgPartitionPoint;

constexpr auto operation = core::Operation::EvaluateStep;

R orientation(pg::ReferenceArithmetic& m, const R& ax, const R& ay, const R& bx, const R& by,
              const R& px, const R& py) {
    return m.sub(m.mul(m.sub(bx, ax), m.sub(py, ay)), m.mul(m.sub(by, ay), m.sub(px, ax)));
}

} // namespace

std::optional<std::array<R, 2>>
case_material_triangle_coordinates(pg::ReferenceArithmetic& m,
                                   const structure::PgPartitionTriangle& cell,
                                   std::span<const Point> vertices, const R& s, const R& t) {
    m.visit(1);
    const auto& a = vertices[static_cast<std::size_t>(cell.vertices[0])];
    const auto& b = vertices[static_cast<std::size_t>(cell.vertices[1])];
    const auto& c = vertices[static_cast<std::size_t>(cell.vertices[2])];
    auto u = orientation(m, a.s, a.t, s, t, c.s, c.t);
    auto v = orientation(m, a.s, a.t, b.s, b.t, s, t);
    if (u.sign() < 0 || v.sign() < 0 || m.compare(m.add(u, v), cell.twice_area) > 0)
        return std::nullopt;
    return std::array{std::move(u), std::move(v)};
}

namespace {
std::array<structure::PgReferenceVector, 2> chain(pg::ReferenceArithmetic& m,
                                                  const structure::PgReferenceValue& value,
                                                  const CaseMaterialChartPoint& chart) {
    pg::ReferenceRadicals algebra{m, value.g};
    std::array<structure::PgReferenceVector, 2> result;
    const std::array directions{&chart.chord_tangent, &chart.span_tangent};
    for (std::size_t direction = 0; direction < result.size(); ++direction)
        for (std::size_t axis = 0; axis < result[direction].size(); ++axis)
            for (std::size_t component = 0; component < value.jacobian_columns_metres.size();
                 ++component)
                result[direction][axis] =
                    algebra.add(result[direction][axis],
                                algebra.scale(value.jacobian_columns_metres[component][axis],
                                              (*directions[direction])[component]));
    return result;
}

bool equal(pg::ReferenceRadicals& algebra, const structure::PgReferenceVector& a,
           const structure::PgReferenceVector& b) {
    for (std::size_t axis = 0; axis < a.size(); ++axis)
        if (algebra.sign(algebra.sub(a[axis], b[axis])) != 0)
            return false;
    return true;
}

structure::PgReferencePoint
locate_query(pg::ReferenceArithmetic& m, CaseMaterialSide side, core::Index triangle,
             core::Index layer, const std::array<R, 3>& local,
             std::span<const structure::PgPartitionVolume> volumes,
             std::span<const structure::PgPartitionHorizontalFace> horizontal) {
    if (side == CaseMaterialSide::Interior) {
        for (std::size_t i = 0; i < volumes.size(); ++i) {
            m.visit(1);
            if (volumes[i].triangle == triangle && volumes[i].layer == layer)
                return structure::PgReferenceVolumePoint{partition_index(i), local[0], local[1],
                                                         local[2]};
        }
    } else {
        const auto expected = side == CaseMaterialSide::Lower
                                  ? structure::PgPartitionFaceSide::Lower
                                  : structure::PgPartitionFaceSide::Upper;
        for (std::size_t i = 0; i < horizontal.size(); ++i) {
            m.visit(1);
            const auto& face = horizontal[i];
            if (face.triangle == triangle && face.layer == layer && face.side == expected)
                return structure::PgReferenceHorizontalPoint{partition_index(i), local[0],
                                                             local[1]};
        }
    }
    partition_fail(FailureCode::InvalidInput);
}

template <class Build>
core::Result<CaseMaterialReferenceSample>
evaluate(const CaseMaterialChart& chart, const structure::PgReferenceLimits& limits,
         core::ResourcePlan& plan, bool requires_tangents, Build build) noexcept {
    return core::guard_entry(
        operation, FailureCode::PrecisionUnresolved, {},
        [&]() -> core::Result<CaseMaterialReferenceSample> {
            try {
                const auto initial_bytes = plan.allocated_bytes();
                auto point = partition_take(build());
                pg::ReferenceArithmetic m{limits, plan};
                m.charge(point.call_info.scalar_operations);
                m.visit(point.call_info.partition_visits);
                m.admit(pg::reference_live_rationals + 128);
                const auto* support = partition_take(chart.support());
                const auto source = partition_take(support->source());
                const auto vertices = partition_take(source.partition.vertices());
                const auto triangles = partition_take(source.partition.triangles());
                const auto layers = partition_take(source.partition.layers());
                const auto volumes = partition_take(source.partition.volumes());
                const auto horizontal = partition_take(source.partition.horizontal_faces());
                std::optional<std::size_t> layer;
                for (std::size_t i = 0; i < layers.size(); ++i) {
                    m.visit(1);
                    if (layers[i].id == point.layer) {
                        layer = i;
                        break;
                    }
                }
                if (!layer)
                    partition_fail(FailureCode::InvalidInput);
                const auto& layer_range = layers[*layer];
                const auto depth = m.div(m.sub(point.material[2], layer_range.lower),
                                         m.sub(layer_range.upper, layer_range.lower));
                // Count exact containments before retaining any alternate records.
                std::size_t candidates = 0;
                for (const auto& cell : triangles)
                    if (case_material_triangle_coordinates(m, cell, vertices, point.material[0],
                                                           point.material[1]))
                        ++candidates;
                if (candidates == 0)
                    partition_fail(FailureCode::GeometryInvalid);
                std::optional<CaseMaterialReferenceSample> result;
                for (std::size_t i = 0; i < triangles.size(); ++i) {
                    if (result && result->alternate_count == candidates - 1)
                        break;
                    const auto& cell = triangles[i];
                    const auto uv = case_material_triangle_coordinates(
                        m, cell, vertices, point.material[0], point.material[1]);
                    if (!uv)
                        continue;
                    std::array<R, 3> local{m.div((*uv)[0], cell.twice_area),
                                           m.div((*uv)[1], cell.twice_area), m.copy(depth)};
                    const auto query =
                        locate_query(m, point.side, partition_index(i), partition_index(*layer),
                                     local, volumes, horizontal);
                    auto remaining = limits;
                    remaining.maximum_scalar_operations -= m.scalar_operations;
                    remaining.maximum_partition_visits -= m.partition_visits;
                    const auto retained = result ? result->alternate_count : 0;
                    m.reserve_nested_rationals(
                        remaining.algebra, pg::partition_add(128, pg::partition_mul(3, retained)));
                    auto reference =
                        partition_take(support->evaluate(query, std::nullopt, remaining, plan));
                    m.charge(reference.call_info.scalar_operations);
                    m.visit(reference.call_info.partition_visits);
                    std::optional<std::array<structure::PgReferenceVector, 2>> tangents;
                    if (requires_tangents)
                        tangents.emplace(chain(m, reference, point));
                    if (!result) {
                        result.emplace(CaseMaterialReferenceSample{{},
                                                                   std::move(reference),
                                                                   std::move(tangents),
                                                                   partition_index(i),
                                                                   std::move(local),
                                                                   {},
                                                                   0,
                                                                   {}});
                        continue;
                    }
                    // An interior aerodynamic point has a two-dimensional neighbourhood.
                    // At an actual PG cut both limits must produce the same physical tangent.
                    const auto& first = result->reference;
                    pg::ReferenceRadicals algebra{m, first.g};
                    if (reference.g != first.g ||
                        !equal(algebra, first.position_metres, reference.position_metres) ||
                        result->tangents.has_value() != tangents.has_value())
                        partition_fail(FailureCode::GeometryInvalid);
                    if (tangents && result->tangents &&
                        (!equal(algebra, (*result->tangents)[0], (*tangents)[0]) ||
                         !equal(algebra, (*result->tangents)[1], (*tangents)[1])))
                        partition_fail(FailureCode::GeometryInvalid);
                    if (result->alternate.size() == 0)
                        result->alternate =
                            partition_take(plan.allocate<CaseMaterialReferenceCandidate>(
                                partition_index(candidates - 1), operation));
                    result->alternate.values()[result->alternate_count++] = {
                        reference.partition_record, std::move(local)};
                    m.reserve_nested_rationals(m.leaf_limits, 3);
                }
                if (!result)
                    partition_fail(FailureCode::GeometryInvalid);
                m.cancel();
                result->chart = std::move(point);
                result->call_info = m.info();
                result->call_info.allocated_owned_payload_bytes =
                    plan.allocated_bytes() - initial_bytes;
                return std::move(*result);
            } catch (const pg::PartitionAbort& rejected) {
                return std::unexpected(
                    core::Failure{rejected.failure.code, operation, rejected.failure.details});
            }
        });
}
} // namespace

core::Result<CaseMaterialReferenceSample> evaluate_case_material_collocation_reference(
    const CaseMaterialChart& chart, core::Index chord, core::Index span,
    const structure::PgReferenceLimits& limits, core::ResourcePlan& plan) noexcept {
    return evaluate(chart, limits, plan, true,
                    [&] { return chart.collocation(chord, span, limits, plan); });
}

core::Result<CaseMaterialReferenceSample> evaluate_case_material_vertex_reference(
    const CaseMaterialChart& chart, core::Index chord, core::Index span, core::Index owning_span,
    const structure::PgReferenceLimits& limits, core::ResourcePlan& plan) noexcept {
    return evaluate(chart, limits, plan, false,
                    [&] { return chart.vertex(chord, span, owning_span, limits, plan); });
}

namespace {
bool same_schema(const io::SchemaIdentity& a, const io::SchemaIdentity& b) noexcept {
    return a.source == b.source && a.contract_revision == b.contract_revision &&
           a.representation_revision == b.representation_revision;
}

void verify_cartesian_reference(pg::ReferenceArithmetic& m,
                                const CaseMaterialReferenceSample& actual,
                                const CaseReferenceSurfaceSample& expected, const R& tolerance) {
    pg::ReferenceRadicals algebra{m, actual.reference.g};
    const structure::PgReferenceScalar squared{m.square(tolerance), {}};
    const auto difference = [&](const structure::PgReferenceVector& left,
                                const std::array<R, 3>& right) {
        structure::PgReferenceVector result;
        for (std::size_t axis = 0; axis < result.size(); ++axis)
            result[axis] = {m.sub(left[axis].rational, right[axis]), m.copy(left[axis].radical)};
        return result;
    };
    const auto check_error = [&](const structure::PgReferenceVector& value, const R& scale) {
        const auto distance = algebra.scale(algebra.dot(value, value), m.square(scale));
        if (algebra.sign(algebra.sub(distance, squared)) > 0)
            partition_fail(FailureCode::GeometryInvalid);
    };
    check_error(difference(actual.reference.position_metres, expected.position), m.constant("1"));
    const auto& tangents = actual.tangents;
    if (!tangents)
        return;
    check_error(difference((*tangents)[0], expected.chord_tangent), actual.chart.panel_extent[0]);
    check_error(difference((*tangents)[1], expected.span_tangent), actual.chart.panel_extent[1]);
    structure::PgReferenceVector x, y;
    for (std::size_t axis = 0; axis < x.size(); ++axis) {
        x[axis].rational = m.copy(expected.chord_tangent[axis]);
        y[axis].rational = m.copy(expected.span_tangent[axis]);
    }
    const auto actual_normal = algebra.cross((*tangents)[0], (*tangents)[1]);
    if (algebra.sign(algebra.dot(actual_normal, algebra.cross(x, y))) <= 0)
        partition_fail(FailureCode::GeometryInvalid);
}

template <class Build>
core::Result<CaseMaterialReferenceSample> qualify(const CaseMaterialChart& chart,
                                                  const CaseReferenceSurface& thin,
                                                  const structure::PgReferenceLimits& limits,
                                                  core::ResourcePlan& plan, Build build) noexcept {
    return core::guard_entry(
        operation, FailureCode::PrecisionUnresolved, {},
        [&]() -> core::Result<CaseMaterialReferenceSample> {
            try {
                const auto initial_bytes = plan.allocated_bytes();
                const auto material_identity = partition_take(chart.identity());
                const auto reference_identity = partition_take(thin.identity());
                if (material_identity.geometry_id != reference_identity.geometry_id ||
                    !same_schema(material_identity.source, reference_identity.source))
                    partition_fail(FailureCode::StaleGeneration);
                auto result = partition_take(build());
                pg::ReferenceArithmetic m{limits, plan};
                m.charge(result.call_info.scalar_operations);
                m.visit(result.call_info.partition_visits);
                m.admit(pg::reference_live_rationals + 128);
                auto remaining = limits;
                remaining.maximum_scalar_operations -= m.scalar_operations;
                remaining.maximum_partition_visits -= m.partition_visits;
                m.reserve_nested_rationals(remaining.algebra, 128);
                const auto independent = partition_take(thin.evaluate(
                    result.chart.aerodynamic[0], result.chart.aerodynamic[1], remaining, plan));
                m.charge(independent.call_info.scalar_operations);
                m.visit(independent.call_info.partition_visits);
                verify_cartesian_reference(m, result, independent,
                                           *partition_take(chart.geometry_error_metres()));
                m.cancel();
                result.call_info = m.info();
                result.call_info.allocated_owned_payload_bytes =
                    plan.allocated_bytes() - initial_bytes;
                return result;
            } catch (const pg::PartitionAbort& rejected) {
                return std::unexpected(
                    core::Failure{rejected.failure.code, operation, rejected.failure.details});
            }
        });
}
} // namespace

core::Result<CaseMaterialReferenceSample> qualify_case_material_collocation_reference(
    const CaseMaterialChart& chart, const CaseReferenceSurface& thin, core::Index chord,
    core::Index span, const structure::PgReferenceLimits& limits,
    core::ResourcePlan& plan) noexcept {
    return qualify(chart, thin, limits, plan, [&] {
        return evaluate_case_material_collocation_reference(chart, chord, span, limits, plan);
    });
}

core::Result<CaseMaterialReferenceSample> qualify_case_material_vertex_reference(
    const CaseMaterialChart& chart, const CaseReferenceSurface& thin, core::Index chord,
    core::Index span, core::Index owning_span, const structure::PgReferenceLimits& limits,
    core::ResourcePlan& plan) noexcept {
    return qualify(chart, thin, limits, plan, [&] {
        return evaluate_case_material_vertex_reference(chart, chord, span, owning_span, limits,
                                                       plan);
    });
}

} // namespace uvlas::coupling::detail
