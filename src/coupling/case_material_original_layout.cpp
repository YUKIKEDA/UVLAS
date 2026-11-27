#include "case_material_original_layout.hpp"

#include "case_material_alias.hpp"
#include "case_pg_support_motion.hpp"
#include "structure/pg_reference_arithmetic.hpp"

#include <algorithm>

namespace uvlas::coupling::detail {
namespace {
namespace a = aero::detail;
namespace pg = structure::detail;
using core::FailureCode;
using pg::partition_add;
using pg::partition_fail;
using pg::partition_index;
using pg::partition_mul;
using pg::partition_take;
using R = structure::ExactRational;

template <class T, class E>
T aero_result(std::expected<T, E> result) {
    if (!result) {
        if constexpr (requires { result.error().cause; })
            throw pg::PartitionAbort{result.error().cause};
        else
            throw pg::PartitionAbort{result.error().origin.cause};
    }
    return std::move(*result);
}

template <class F>
auto entry(core::Operation operation, F action) noexcept -> std::invoke_result_t<F> {
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

R exact_word(pg::ReferenceArithmetic& m, const a::ScalarWords& word) {
    const auto from = [&](double x) {
        m.charge();
        return partition_take(R::from_binary64(x, m.leaf_limits, m.plan));
    };
    if (const auto* single = std::get_if<a::SingleWord>(&word))
        return from(single->value);
    const auto pair = std::get<a::DoubleWord>(word);
    return m.add(from(pair.high), from(pair.low));
}

void match_corner(pg::ReferenceArithmetic& m, const structure::PgReferenceValue& physical,
                  const std::array<a::ScalarWords, 3>& retained, const R& allowance) {
    pg::ReferenceRadicals algebra{m, physical.g};
    structure::PgReferenceVector difference;
    for (std::size_t axis = 0; axis < difference.size(); ++axis)
        difference[axis] = {
            m.sub(physical.position_metres[axis].rational, exact_word(m, retained[axis])),
            m.copy(physical.position_metres[axis].radical)};
    const auto squared = algebra.dot(difference, difference);
    const structure::PgReferenceScalar limit{m.square(allowance), {}};
    if (algebra.sign(algebra.sub(squared, limit)) > 0)
        partition_fail(FailureCode::GeometryInvalid);
}

bool same_surface(const CaseMaterialChartIdentity& a, const CaseMaterialChartIdentity& b) noexcept {
    return a.surface_id == b.surface_id && a.shell_id == b.shell_id &&
           a.geometry_id == b.geometry_id && a.reversed == b.reversed &&
           a.source.source == b.source.source &&
           a.source.contract_revision == b.source.contract_revision &&
           a.source.representation_revision == b.source.representation_revision;
}

auto key(const CaseMaterialOriginalCorner& row) noexcept {
    return std::pair{row.facet, row.local_vertex};
}

} // namespace

CaseMaterialOriginalLayout::CaseMaterialOriginalLayout(
    CaseMaterialChartIdentity surface, a::GeometryContext context,
    std::array<std::byte, 32> geometry,
    core::PlannedBuffer<CaseMaterialOriginalCorner> corners) noexcept
    : surface_(surface), context_(context), geometry_(geometry), corners_(std::move(corners)) {}

core::Result<CaseMaterialOriginalLayout> CaseMaterialOriginalLayout::prepare(
    const CaseMaterialChart& chart, const CaseReferenceSurface& thin,
    const a::GlobalContactState& contact, std::span<const a::PhysicalFaceId> faces,
    const structure::PgReferenceLimits& limits, core::ResourcePlan& plan) noexcept {
    constexpr auto operation = core::Operation::PrepareModel;
    return entry(operation, [&]() -> core::Result<CaseMaterialOriginalLayout> {
        pg::ReferenceArithmetic math{limits, plan};
        math.cancel();
        const auto dims = partition_take(chart.panel_shape());
        const auto nc = static_cast<std::uint64_t>(dims[0]),
                   ns = static_cast<std::uint64_t>(dims[1]);
        if (faces.size() != partition_mul(nc, ns))
            partition_fail(FailureCode::InvalidInput);
        const auto surface = partition_take(chart.identity());
        const auto snapshot = aero_result(contact.geometry());
        const auto& geometry = snapshot.get();
        const auto& tolerance = *partition_take(chart.geometry_error_metres());
        const auto columns = partition_add(nc, 1), width = partition_add(ns, 1);
        auto grid = partition_take(plan.allocate<std::optional<a::PointId>>(
            partition_index(partition_mul(columns, width)), operation));
        auto corners = partition_take(plan.allocate<CaseMaterialOriginalCorner>(
            partition_index(partition_mul(6, faces.size())), operation));
        math.admit();
        struct Panel {
            a::PhysicalFaceId face{0};
            std::size_t ordinal{};
            std::optional<a::UvlmOriginalQuad> quad;
        };
        auto panels =
            partition_take(plan.allocate<Panel>(partition_index(faces.size()), operation));
        for (std::size_t i = 0; i < faces.size(); ++i) {
            math.visit(1);
            panels.values()[i] = {faces[i], i, {}};
        }
        std::ranges::sort(panels.values(), [&](const Panel& a, const Panel& b) {
            math.visit(1);
            return a.face < b.face;
        });
        for (std::size_t i = 1; i < panels.size(); ++i) {
            math.visit(1);
            if (panels.values()[i - 1].face == panels.values()[i].face)
                partition_fail(FailureCode::InvalidInput);
        }
        const auto count = aero_result(geometry.counts()).quads;
        for (std::size_t i = 0; i < count; ++i) {
            const auto quad = aero_result(contact.uvlm_original_quad(
                i, limits.maximum_partition_visits - math.partition_visits, limits.cancellation));
            math.visit(quad.used_steps);
            const auto records = panels.values();
            const auto panel = std::ranges::lower_bound(
                records, quad.face,
                [&](auto a, auto b) {
                    math.visit(1);
                    return a < b;
                },
                &Panel::face);
            if (panel != records.end() && panel->face == quad.face) {
                if (quad.role != a::AeroFaceRole::Bound || panel->quad)
                    partition_fail(FailureCode::InvalidInput);
                panel->quad = quad;
            }
        }
        constexpr std::array<std::array<std::size_t, 3>, 2> children{{{0, 2, 1}, {0, 3, 2}}};
        constexpr std::array<std::size_t, 4> chord_offsets{0, 0, 1, 1};
        const std::array<std::size_t, 4> span_offsets =
            surface.reversed ? std::array<std::size_t, 4>{1, 0, 0, 1}
                             : std::array<std::size_t, 4>{0, 1, 1, 0};
        for (const auto& panel : panels.values()) {
            math.visit(1);
            if (!panel.quad)
                partition_fail(FailureCode::InvalidInput);
            const auto f = panel.ordinal;
            const auto chord = f / ns, span = f % ns;
            const auto& quad = *panel.quad;
            const auto current_limits = [&] {
                auto caps = limits;
                caps.maximum_scalar_operations -= math.scalar_operations;
                caps.maximum_partition_visits -= math.partition_visits;
                return caps;
            };
            const auto collocation_work =
                partition_take(qualify_case_material_collocation_reference(
                                   chart, thin, partition_index(chord), partition_index(span),
                                   current_limits(), plan))
                    .call_info;
            math.charge(collocation_work.scalar_operations);
            math.visit(collocation_work.partition_visits);

            std::array<a::PointId, 4> owners{a::PointId{0}, a::PointId{0}, a::PointId{0},
                                             a::PointId{0}};
            std::array<CaseMaterialVertexSite, 4> sites;
            for (std::size_t i = 0; i < sites.size(); ++i) {
                const auto c = chord + chord_offsets[i], s = span + span_offsets[i];
                sites[i] = {partition_index(c), partition_index(s), partition_index(span)};
                owners[i] = aero_result(geometry.owner(quad.corners[i]));
                auto& shared = grid.values()[static_cast<std::size_t>(c * width + s)];
                if (shared && *shared != owners[i])
                    partition_fail(FailureCode::GeometryInvalid);
                shared = owners[i];
                const auto point = partition_take(qualify_case_material_vertex_reference(
                    chart, thin, sites[i].chord_column, sites[i].span_knot,
                    sites[i].original_span_cell, current_limits(), plan));
                math.charge(point.call_info.scalar_operations);
                math.visit(point.call_info.partition_visits);
                match_corner(
                    math, point.reference,
                    aero_result(geometry.point_words(quad.corners[i], a::Configuration::Reference)),
                    tolerance);
            }
            for (std::size_t child = 0; child < children.size(); ++child) {
                const auto facet = quad.facets[child];
                const auto actual = aero_result(geometry.facet_points(facet));
                for (std::size_t local = 0; local < actual.size(); ++local) {
                    math.visit(1);
                    const auto i = children[child][local];
                    if (aero_result(geometry.owner(actual[local])) != owners[i])
                        partition_fail(FailureCode::GeometryInvalid);
                    corners.values()[6 * f + 3 * child + local] = {
                        faces[f], facet, static_cast<std::uint8_t>(local), owners[i], sites[i]};
                }
            }
        }
        std::ranges::sort(grid.values(), [&](const auto& first, const auto& second) {
            math.visit(1);
            return first < second;
        });
        for (std::size_t i = 0; i < grid.size(); ++i) {
            math.visit(1);
            if (!grid.values()[i] || (i != 0 && grid.values()[i - 1] == grid.values()[i]))
                partition_fail(FailureCode::GeometryInvalid);
        }
        std::ranges::sort(corners.values(), [&](const auto& first, const auto& second) {
            math.visit(1);
            return key(first) < key(second);
        });
        for (std::size_t i = 1; i < corners.size(); ++i) {
            math.visit(1);
            if (key(corners.values()[i - 1]) == key(corners.values()[i]))
                partition_fail(FailureCode::InvalidInput);
        }
        math.cancel();
        return CaseMaterialOriginalLayout{surface, geometry.context(),
                                          aero_result(geometry.identity()), std::move(corners)};
    });
}

core::Result<CaseMaterialOriginalCorner>
CaseMaterialOriginalLayout::corner(a::FacetId facet, std::uint8_t local) const noexcept {
    const auto rows = corners_.values();
    const std::pair expected{facet, local};
    const auto at = std::lower_bound(
        rows.begin(), rows.end(), expected,
        [](const auto& row, const auto& candidate) { return key(row) < candidate; });
    if (at == rows.end() || key(*at) != expected)
        return std::unexpected(
            core::Failure{FailureCode::InvalidInput, core::Operation::EvaluateStep});
    return *at;
}

core::Result<CaseMaterialOriginalEnclosedForce> CaseMaterialOriginalLayout::pullback(
    const DiscreteCirculationLoads& loads, std::size_t index, const CaseMaterialChart& chart,
    const CaseReferenceSurface& thin, const structure::PgPreparedPatch& patch,
    const structure::PgPreparedReferenceProjection& projection,
    const structure::PgReferenceLimits& preparation, const structure::PgVolumePositionLimits& apply,
    core::ResourcePlan& plan) const noexcept {
    return entry(
        core::Operation::EvaluateStep, [&]() -> core::Result<CaseMaterialOriginalEnclosedForce> {
            if (loads.source().context != context_ || loads.source().geometry != geometry_ ||
                !same_surface(surface_, partition_take(chart.identity())))
                partition_fail(FailureCode::StaleGeneration);
            const auto load = partition_take(loads.original(index));
            const auto row = partition_take(corner(load.facet, load.local_vertex));
            if (row.face != load.face || row.canonical != load.canonical)
                partition_fail(FailureCode::StaleGeneration);
            const auto material = partition_take(CaseMaterialRows::prepare(
                chart, thin, patch, row.site, projection, 0., preparation, plan,
                CaseMaterialLoadKey{row.face, row.facet, row.local_vertex, row.canonical}));
            return material.pullback_enclosed(partition_take(patch.info()).identity, load, apply,
                                              plan);
        });
}

core::Result<void> qualify_case_material_original_alias(const CaseMaterialOriginalAliasPoint& left,
                                                        const CaseMaterialOriginalAliasPoint& right,
                                                        const CaseStructuralSupport& support,
                                                        const structure::PgReferenceLimits& limits,
                                                        core::Workspace& workspace,
                                                        core::ResourcePlan& plan) noexcept {
    return entry(core::Operation::PrepareModel, [&]() -> core::Result<void> {
        pg::ReferenceArithmetic math{limits, plan};
        math.cancel();
        const auto source = support.source();
        const auto same_case = [&](const io::SchemaIdentity& candidate) {
            return candidate.source == source.source &&
                   candidate.contract_revision == source.contract_revision &&
                   candidate.representation_revision == source.representation_revision;
        };
        if (left.layout.context_ != right.layout.context_ ||
            left.layout.geometry_ != right.layout.geometry_ ||
            !same_case(left.layout.surface_.source) || !same_case(right.layout.surface_.source))
            partition_fail(FailureCode::StaleGeneration);
        const auto a = partition_take(left.layout.corner(left.facet, left.local_vertex));
        const auto b = partition_take(right.layout.corner(right.facet, right.local_vertex));
        if (a.canonical != b.canonical)
            partition_fail(FailureCode::GeometryInvalid);
        const auto remaining = [&](std::uint64_t retained) {
            auto caps = limits;
            caps.maximum_scalar_operations -= math.scalar_operations;
            caps.maximum_partition_visits -= math.partition_visits;
            math.reserve_nested_rationals(caps.algebra, retained);
            return caps;
        };
        const auto prepare = [&](const CaseMaterialOriginalAliasPoint& point, const auto& corner,
                                 std::uint64_t retained) {
            if (!same_surface(point.layout.surface_, partition_take(point.chart.identity())))
                partition_fail(FailureCode::StaleGeneration);
            auto row = partition_take(CaseMaterialAliasRow::prepare(
                point.chart, point.thin, point.patch, corner.site, remaining(retained), plan));
            math.charge(row.call_info().scalar_operations);
            math.visit(row.call_info().partition_visits);
            return row;
        };
        const auto payload = [](const CaseMaterialAliasRow& row) {
            // The conditioned column retains value and its exact gradient.
            return partition_add(1, partition_mul(8, row.coefficients().size()));
        };
        const auto a_row = prepare(left, a, 0);
        const auto b_row = prepare(right, b, payload(a_row));
        return qualify_case_material_constraint_alias(
            a_row, b_row, support.support(),
            remaining(partition_add(payload(a_row), payload(b_row))), workspace, plan);
    });
}

} // namespace uvlas::coupling::detail
