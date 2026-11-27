#include "case_material_rows.hpp"

#include "aero/local_geometry_bounds_private.hpp"
#include "aero/p1_interval_private.hpp"
#include "aero/p1_shared_kernel.hpp"
#include "case_material_alias.hpp"
#include "structure/pg_reference_arithmetic.hpp"

#include <cmath>
#include <memory>
#include <memory_resource>
#include <utility>

namespace uvlas::coupling::detail {
namespace {
namespace pg = structure::detail;
namespace s = structure;
using core::FailureCode;
using pg::partition_add;
using pg::partition_fail;
using pg::partition_index;
using pg::partition_mul;
using pg::partition_take;

template <class F>
auto run(core::Operation operation, F action) noexcept -> std::invoke_result_t<F> {
    return core::guard_entry(
        operation, FailureCode::PrecisionUnresolved, {}, [&]() -> std::invoke_result_t<F> {
            try {
                return action();
            } catch (const pg::PartitionAbort& rejected) {
                return std::unexpected(
                    core::Failure{rejected.failure.code, operation, rejected.failure.details});
            }
        });
}

template <class T>
std::uint64_t storage_bytes(std::size_t entries) {
    return partition_add(partition_mul(entries, sizeof(T)), alignof(T) - 1);
}

void match_reference(const CaseMaterialChart& chart, const s::PgPreparedPatch& patch,
                     pg::ReferenceArithmetic& m) {
    const auto* support = partition_take(chart.support());
    const auto info = partition_take(patch.info());
    if (info.geometry != partition_take(support->identity()))
        partition_fail(FailureCode::StaleGeneration);
    const auto* basis = partition_take(patch.conditioning().reference_basis());
    const auto a_frame = partition_take(support->frame());
    const auto b_frame = partition_take(basis->frame());
    m.charge(13);
    if (a_frame.rotation != b_frame.rotation || a_frame.origin_metres != b_frame.origin_metres ||
        a_frame.length_metres != b_frame.length_metres)
        partition_fail(FailureCode::StaleGeneration);
    const auto a = partition_take(support->source());
    const auto b = partition_take(basis->source());
    const auto edges = partition_take(a.partition.sources());
    const auto points = partition_take(a.partition.source_points());
    const auto thickness = partition_take(a.partition.thickness());
    const auto layers = partition_take(a.partition.layers());
    const auto names = partition_take(a.partition.boundaries());
    const auto text = partition_take(a.partition.boundary_text());
    std::size_t domain_count = 0;
    for (const auto& edge : edges) {
        m.visit(1);
        if (std::holds_alternative<s::PgPartitionDomainSource>(edge.descriptor))
            ++domain_count;
    }
    auto scratch = partition_take(m.plan.allocate<std::byte>(
        partition_index(storage_bytes<s::PgPartitionLayer>(layers.size())),
        core::Operation::PrepareModel));
    std::pmr::monotonic_buffer_resource resource{scratch.values().data(), scratch.size(),
                                                 std::pmr::null_memory_resource()};
    auto* layer_view =
        std::pmr::polymorphic_allocator<s::PgPartitionLayer>{&resource}.allocate(layers.size());
    auto boundary_view = partition_take(m.plan.allocate<s::PgPartitionBoundary>(
        partition_index(names.size()), core::Operation::PrepareModel));
    for (std::size_t i = 0; i < layers.size(); ++i) {
        m.visit(1);
        const auto& layer = layers[i];
        std::construct_at(layer_view + i, layer.id, layer.lower, layer.upper);
    }
    for (std::size_t i = 0; i < names.size(); ++i) {
        const auto& name = names[i];
        m.visit(partition_add(1, static_cast<std::uint64_t>(name.point_count)));
        const auto label = text.subspan(static_cast<std::size_t>(name.text_offset),
                                        static_cast<std::size_t>(name.text_count));
        boundary_view.values()[i] = {{label.data(), label.size()},
                                     points.subspan(static_cast<std::size_t>(name.first_point),
                                                    static_cast<std::size_t>(name.point_count))};
    }
    m.visit(partition_add(domain_count,
                          partition_add(thickness.vertices.size(), thickness.triangles.size())));
    m.visit(partition_add(text.size(), layers.size()));
    // A Case-created PG reference has no separate computational cuts. Compare
    // the actual chart inputs; the reference-content stamp is only an early check.
    const s::PgPartitionRequest expected{
        a.graph,
        points.first(domain_count),
        {thickness.vertices, thickness.triangles, thickness.values},
        {layer_view, layers.size()},
        boundary_view.values(),
        {}};
    partition_take(b.partition.validate_for(expected));
    m.cancel();
}

s::PgReferenceLimits remaining(const pg::ReferenceArithmetic& m) {
    auto active = m.limits;
    active.maximum_scalar_operations -= m.scalar_operations;
    active.maximum_partition_visits -= m.partition_visits;
    return active;
}

void use(pg::ReferenceArithmetic& m, const s::PgReferenceCallInfo& counts) {
    m.charge(counts.scalar_operations);
    m.visit(counts.partition_visits);
}

s::PgReferencePoint query_for(CaseMaterialSide side, core::Index record,
                              const std::array<s::ExactRational, 3>& local) {
    if (side == CaseMaterialSide::Interior)
        return s::PgReferenceVolumePoint{record, local[0], local[1], local[2]};
    return s::PgReferenceHorizontalPoint{record, local[0], local[1]};
}

void qualify_shared_rows(const s::PgPreparedPatch& patch, const CaseMaterialReferenceSample& point,
                         const std::array<s::PgReferencePoint, 1>& canonical,
                         pg::ReferenceArithmetic& m) {
    if (point.alternate_count == 0)
        return;
    const auto info = partition_take(patch.info());
    // Each exact conditioned column retains N plus three radical gradient components.
    const auto held = partition_add(512, partition_mul(8, info.basis.rank));
    m.admit(held);
    auto caps = remaining(m);
    m.reserve_nested_rationals(caps.algebra, 256);
    const auto original = partition_take(patch.evaluate_reference(canonical, caps, m.plan));
    use(m, original.call_info);
    const auto original_columns = original.columns.values();
    pg::ReferenceRadicals algebra{m, point.reference.g};
    for (const auto& alternate : point.alternate.values().first(point.alternate_count)) {
        const std::array query{
            query_for(point.chart.side, alternate.partition_record, alternate.local)};
        caps = remaining(m);
        m.reserve_nested_rationals(caps.algebra, held);
        const auto comparison = partition_take(patch.evaluate_reference(query, caps, m.plan));
        use(m, comparison.call_info);
        const auto columns = comparison.columns.values();
        if (columns.size() != original_columns.size())
            partition_fail(FailureCode::StaleGeneration);
        for (std::size_t i = 0; i < columns.size(); ++i) {
            m.visit(1);
            const auto& a = original_columns[i];
            const auto& b = columns[i];
            if (algebra.sign(algebra.sub(a.value, b.value)) != 0)
                partition_fail(FailureCode::PrecisionUnresolved);
            if (const auto& tangents = point.tangents)
                for (const auto& tangent : *tangents)
                    if (algebra.sign(
                            algebra.sub(algebra.dot(a.physical_gradient_per_metre, tangent),
                                        algebra.dot(b.physical_gradient_per_metre, tangent))) != 0)
                        partition_fail(FailureCode::PrecisionUnresolved);
        }
    }
}
} // namespace

CaseMaterialRows::CaseMaterialRows(CaseMaterialChartIdentity surface, CaseMaterialSite site,
                                   s::PgProjectedReferenceBatch&& position,
                                   std::optional<s::PgMaterialTangentReference>&& tangents,
                                   std::optional<CaseMaterialLoadKey> original_load) noexcept
    : surface_(surface),
      site_(site),
      position_(std::move(position)),
      tangents_(std::move(tangents)),
      original_load_(original_load) {}

core::Result<CaseMaterialRows>
CaseMaterialRows::prepare(const CaseMaterialChart& chart, const CaseReferenceSurface& reference,
                          const s::PgPreparedPatch& patch, CaseMaterialSite site,
                          const s::PgPreparedReferenceProjection& projection,
                          double maximum_tangent_row_error, const s::PgReferenceLimits& limits,
                          core::ResourcePlan& plan,
                          std::optional<CaseMaterialLoadKey> original_load) noexcept {
    return run(core::Operation::PrepareModel, [&]() -> core::Result<CaseMaterialRows> {
        pg::ReferenceArithmetic m{limits, plan};
        m.cancel();
        if (original_load && (!std::holds_alternative<CaseMaterialVertexSite>(site) ||
                              original_load->local_vertex >= 3))
            partition_fail(FailureCode::InvalidInput);
        match_reference(chart, patch, m);
        auto selected = partition_take(std::visit(
            [&](const auto& key) -> core::Result<CaseMaterialReferenceSample> {
                if constexpr (std::is_same_v<std::remove_cvref_t<decltype(key)>,
                                             CaseMaterialCollocationSite>)
                    return qualify_case_material_collocation_reference(
                        chart, reference, key.chord_cell, key.span_cell, remaining(m), plan);
                else
                    return qualify_case_material_vertex_reference(
                        chart, reference, key.chord_column, key.span_knot, key.original_span_cell,
                        remaining(m), plan);
            },
            site));
        use(m, selected.call_info);
        // Coincident physical coordinates alone do not prove identical basis rows.
        const std::array points{query_for(selected.chart.side, selected.reference.partition_record,
                                          selected.triangle_coordinates)};
        qualify_shared_rows(patch, selected, points, m);
        auto position =
            partition_take(patch.project_reference(points, projection, remaining(m), plan));
        use(m, position.call_info());
        std::optional<s::PgMaterialTangentReference> tangents;
        if (std::holds_alternative<CaseMaterialCollocationSite>(site)) {
            const std::array<s::PgMaterialDirections, 1> directions{{
                {std::move(selected.chart.chord_tangent), std::move(selected.chart.span_tangent)},
            }};
            tangents.emplace(partition_take(s::PgMaterialTangentReference::prepare(
                patch, points, directions, {projection.representation, maximum_tangent_row_error},
                remaining(m), plan)));
            use(m, tangents->call_info());
        }
        m.cancel();
        return CaseMaterialRows{partition_take(chart.identity()), site, std::move(position),
                                std::move(tangents), original_load};
    });
}

core::Result<CaseMaterialPointKinematics>
CaseMaterialRows::evaluate(s::PgPreparedPatchIdentity patch, s::PgCartesianWords q,
                           s::PgCartesianWords v, const CaseMaterialActionLimits& limits,
                           core::Workspace& workspace, core::ResourcePlan& plan) const noexcept {
    return run(core::Operation::EvaluateStep, [&]() -> core::Result<CaseMaterialPointKinematics> {
        if (position_.points().empty())
            partition_fail(FailureCode::InvalidInput);
        if (patch != position_.source())
            partition_fail(FailureCode::StaleGeneration);
        if (limits.cancellation.stop_requested())
            partition_fail(FailureCode::Cancelled);
        const auto precision = position_.projection().representation;
        auto position = partition_take(s::PgVolumePositionResult::allocate(3, precision, plan));
        auto velocity = partition_take(s::PgVolumePositionResult::allocate(3, precision, plan));
        std::optional<s::PgVolumePositionResult> tangents;
        if (tangents_)
            tangents.emplace(
                partition_take(s::PgVolumePositionResult::allocate(6, precision, plan)));
        const auto action = s::PgVolumePositionDirection::ToPoints;
        auto steps = limits.maximum_steps;
        partition_take(s::apply_pg_prepared_reference_position(
            position_, patch, action, q, {steps, limits.position_error, limits.cancellation},
            workspace, position));
        steps -= position.used_steps();
        partition_take(s::apply_pg_prepared_reference_position(
            position_, patch, action, v, {steps, limits.velocity_error, limits.cancellation},
            workspace, velocity));
        steps -= velocity.used_steps();
        if (tangents && tangents_)
            partition_take(s::apply_pg_material_tangents(
                *tangents_, patch, action, q, {steps, limits.tangent_error, limits.cancellation},
                workspace, *tangents));
        return CaseMaterialPointKinematics{std::move(position), std::move(velocity),
                                           std::move(tangents)};
    });
}

core::Result<CaseMaterialOriginalForce>
CaseMaterialRows::pullback(s::PgPreparedPatchIdentity patch, s::PgCartesianWords force,
                           const s::PgVolumePositionLimits& limits, core::Workspace& workspace,
                           core::ResourcePlan& plan) const noexcept {
    return run(core::Operation::EvaluateStep, [&]() -> core::Result<CaseMaterialOriginalForce> {
        if (position_.points().empty())
            partition_fail(FailureCode::InvalidInput);
        if (patch != position_.source())
            partition_fail(FailureCode::StaleGeneration);
        if (limits.cancellation.stop_requested())
            partition_fail(FailureCode::Cancelled);
        auto result = partition_take(s::PgVolumePositionResult::allocate(
            partition_index(partition_mul(3, position_.column_count())),
            position_.projection().representation, plan));
        partition_take(s::apply_pg_prepared_reference_position(
            position_, patch, s::PgVolumePositionDirection::TransposeToColumns, force, limits,
            workspace, result));
        return CaseMaterialOriginalForce{surface_, site_, patch, std::move(result)};
    });
}

core::Result<CaseMaterialOriginalEnclosedForce> CaseMaterialRows::pullback_enclosed(
    s::PgPreparedPatchIdentity patch, const DiscreteOriginalLoadContribution& load,
    const s::PgVolumePositionLimits& limits, core::ResourcePlan& plan) const noexcept {
    namespace aero = uvlas::aero::detail;
    namespace interval = aero::p1;
    constexpr auto operation = core::Operation::EvaluateStep;
    return core::guard_entry(
        operation, FailureCode::PrecisionUnresolved, {},
        [&]() -> core::Result<CaseMaterialOriginalEnclosedForce> {
            try {
                interval::Work work{0, limits.maximum_steps, limits.cancellation};
                work.step();
                const auto* vertex = std::get_if<CaseMaterialVertexSite>(&site_);
                if (!vertex || position_.points().size() != 1)
                    partition_fail(FailureCode::InvalidInput);
                if (patch != position_.source())
                    partition_fail(FailureCode::StaleGeneration);
                if (!original_load_)
                    partition_fail(FailureCode::InvalidInput);
                if (*original_load_ !=
                    CaseMaterialLoadKey{load.face, load.facet, load.local_vertex, load.canonical})
                    partition_fail(FailureCode::StaleGeneration);
                const auto& force = load.force;
                if (!std::isfinite(limits.maximum_absolute_error))
                    partition_fail(FailureCode::Nonfinite);
                if (limits.maximum_absolute_error < 0)
                    partition_fail(FailureCode::InvalidInput);
                const auto unpack = [&](const auto& words) {
                    return std::visit(
                        [&](const auto& value) {
                            if constexpr (std::is_same_v<std::remove_cvref_t<decltype(value)>,
                                                         aero::SingleWord>)
                                return aero::decode_p1_word(value.value);
                            else
                                return interval::add(aero::decode_p1_word(value.high),
                                                     aero::decode_p1_word(value.low), work);
                        },
                        words);
                };
                const auto widen = [&](const interval::Interval& center, double radius) {
                    const auto width = aero::decode_p1_word(radius);
                    return interval::Interval{interval::subtract(center, width, work).lower,
                                              interval::add(center, width, work).upper};
                };
                std::array<interval::Interval, 3> loads;
                for (std::size_t i = 0; i < loads.size(); ++i)
                    loads[i] = interval::scale(widen(unpack(force[i].center()), force[i].radius()),
                                               force[i].exponent());
                const auto high = position_.high_words(), low = position_.low_words();
                const auto projection = position_.projection();
                const auto allowance = aero::decode_p1_word(limits.maximum_absolute_error);
                auto output = partition_take(plan.allocate<std::optional<DiscreteLoadVector>>(
                    partition_index(position_.column_count()), operation));
                aero::LocalRootScratch scratch;
                for (std::size_t column = 0; column < output.size(); ++column) {
                    work.step();
                    // The stored row starts with X0 and measure, then holds each N and gradient.
                    const auto index = 4 * (column + 1);
                    const auto center =
                        low.empty()
                            ? aero::decode_p1_word(high[index])
                            : unpack(aero::ScalarWords{aero::DoubleWord{high[index], low[index]}});
                    const auto coefficient = widen(center, projection.maximum_shape_error);
                    const auto axis = [&](std::size_t i) {
                        const auto product = interval::multiply(coefficient, loads[i], work);
                        const auto word =
                            aero::enclose_p1_value(product, projection.representation, scratch);
                        const auto error =
                            interval::scale(aero::decode_p1_word(word.radius()), word.exponent());
                        if (interval::compare(error.upper, allowance.lower) > 0)
                            partition_fail(FailureCode::PrecisionUnresolved);
                        return word;
                    };
                    output.values()[column].emplace(DiscreteLoadVector{axis(0), axis(1), axis(2)});
                }
                work.step();
                return CaseMaterialOriginalEnclosedForce{surface_, *vertex, patch,
                                                         std::move(output), work.used};
            } catch (const interval::ArithmeticFailure& rejected) {
                return std::unexpected(core::Failure{rejected.code, operation});
            } catch (const pg::PartitionAbort& rejected) {
                return std::unexpected(
                    core::Failure{rejected.failure.code, operation, rejected.failure.details});
            }
        });
}

CaseMaterialAliasRow::CaseMaterialAliasRow(
    CaseMaterialChartIdentity surface, s::PgPreparedPatchIdentity patch,
    s::ExactRational&& radicand, core::PlannedBuffer<s::PgConditionedReferenceColumn>&& row,
    s::PgReferenceCallInfo info) noexcept
    : info_(info),
      surface_(surface),
      patch_(patch),
      radicand_(std::move(radicand)),
      row_(std::move(row)) {}

core::Result<CaseMaterialAliasRow>
CaseMaterialAliasRow::prepare(const CaseMaterialChart& chart, const CaseReferenceSurface& reference,
                              const s::PgPreparedPatch& patch, CaseMaterialVertexSite site,
                              const s::PgReferenceLimits& limits,
                              core::ResourcePlan& plan) noexcept {
    return run(core::Operation::PrepareModel, [&]() -> core::Result<CaseMaterialAliasRow> {
        pg::ReferenceArithmetic m{limits, plan};
        m.cancel();
        match_reference(chart, patch, m);
        auto selected = partition_take(qualify_case_material_vertex_reference(
            chart, reference, site.chord_column, site.span_knot, site.original_span_cell,
            remaining(m), plan));
        use(m, selected.call_info);
        const std::array points{query_for(selected.chart.side, selected.reference.partition_record,
                                          selected.triangle_coordinates)};
        qualify_shared_rows(patch, selected, points, m);
        auto exact = partition_take(patch.evaluate_reference(points, remaining(m), plan));
        use(m, exact.call_info);
        m.cancel();
        return CaseMaterialAliasRow{
            partition_take(chart.identity()), partition_take(patch.info()).identity,
            std::move(selected.reference.g), std::move(exact.columns), m.info()};
    });
}

} // namespace uvlas::coupling::detail
