#include "case_material_alias.hpp"

#include "structure/pg_reference_arithmetic.hpp"

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
using R = s::ExactRational;

struct Duals {
    s::PgSparseChartDual rational, radical;
};

struct View {
    const R& constant;
    const R& radical;
};

bool equal(pg::ReferenceArithmetic& m, View a, const R& g, View b, const R& h) {
    auto x = m.sub(a.constant, b.constant);
    pg::ReferenceRadicals algebra{m, g};
    if (g == h)
        return algebra.sign({std::move(x), m.sub(a.radical, b.radical)}) == 0;
    // x + a.radical sqrt(g) = b.radical sqrt(h), with h positive.
    // Check the sign before squaring, which would otherwise admit the opposite root.
    if (algebra.sign({m.copy(x), m.copy(a.radical)}) != b.radical.sign())
        return false;
    auto difference =
        m.sub(m.add(m.square(x), m.mul(m.square(a.radical), g)), m.mul(m.square(b.radical), h));
    return algebra.sign({std::move(difference), m.mul(m.constant("2"), m.mul(x, a.radical))}) == 0;
}

s::PgConstraintPatchLayout layout_of(pg::ReferenceArithmetic& m, const CaseMaterialAliasRow& row,
                                     const s::PgGlobalConstraints& constraints) {
    for (const auto& patch : constraints.layout()) {
        m.visit(1);
        if (patch.id.value() == row.surface().shell_id) {
            if (patch.basis != row.patch() ||
                static_cast<std::uint64_t>(patch.count) != row.coefficients().size())
                partition_fail(FailureCode::StaleGeneration);
            return patch;
        }
    }
    partition_fail(FailureCode::StaleGeneration);
}

Duals transpose(pg::ReferenceArithmetic& m, const CaseMaterialAliasRow& point,
                const s::PgSparseChart& chart, const s::PgGlobalConstraintInfo& info,
                s::PgSparseChartIdentity identity, std::uint64_t retained,
                core::Workspace& workspace) {
    const auto layout = layout_of(m, point, chart.original());
    auto force = partition_take(
        m.plan.allocate<s::PgCartesianExactPosition>(info.columns, core::Operation::PrepareModel));
    const auto apply = [&](bool radical, std::uint64_t held) {
        const auto row = point.coefficients();
        for (std::size_t i = 0; i < row.size(); ++i) {
            m.visit(1);
            auto& entry = force.values()[static_cast<std::size_t>(layout.first) + i][0];
            entry = m.copy(radical ? row[i].value.radical : row[i].value.rational);
        }
        auto caps = m.limits;
        caps.maximum_scalar_operations -= m.scalar_operations;
        caps.maximum_partition_visits -= m.partition_visits;
        m.reserve_nested_rationals(caps.algebra,
                                   partition_add(held, partition_mul(3, force.size())));
        auto result =
            partition_take(chart.transpose({identity, force.values()}, caps, workspace, m.plan));
        m.charge(result.call_info.scalar_operations);
        m.visit(result.call_info.partition_visits);
        return result;
    };
    auto rational = apply(false, retained);
    const auto count =
        partition_add(rational.independent_force.size(), rational.original_multipliers.size());
    auto radical = apply(true, partition_add(retained, partition_mul(3, count)));
    return {std::move(rational), std::move(radical)};
}

struct PrescribedAction {
    s::PgReferenceVector constant;
    struct Generators {
        s::PgReferenceVector rotation;
        s::PgReferenceScalar translation;
    };
    core::PlannedBuffer<Generators> motion;
};

PrescribedAction prescribed(pg::ReferenceArithmetic& m, const s::PgRigidSupportRows& support,
                            const Duals& dual) {
    auto groups = partition_take(m.plan.allocate<PrescribedAction::Generators>(
        partition_index(support.group_count()), core::Operation::PrepareModel));
    PrescribedAction result{{}, std::move(groups)};
    const auto plain = support.chart().original().original_rows();
    const auto active = support.rows();
    const auto real = dual.rational.original_multipliers.values();
    const auto radial = dual.radical.original_multipliers.values();
    if (plain.size() != active.size() || plain.size() != real.size() ||
        plain.size() != radial.size())
        partition_fail(FailureCode::StaleGeneration);
    for (std::size_t i = 0; i < active.size(); ++i) {
        m.visit(1);
        const auto append = [&](s::PgReferenceScalar& target, const R& source) {
            target.rational = m.add(target.rational, m.mul(real[i][0], source));
            target.radical = m.add(target.radical, m.mul(radial[i][0], source));
        };
        const auto& item = active[i];
        if (item.group) {
            auto& group = result.motion.values()[*item.group];
            append(group.translation, item.translation);
            for (std::size_t axis = 0; axis < group.rotation.size(); ++axis)
                append(group.rotation[axis], item.reference_position[axis]);
        } else
            for (std::size_t axis = 0; axis < result.constant.size(); ++axis)
                append(result.constant[axis], plain[i].prescribed[axis]);
    }
    return result;
}
} // namespace

core::Result<void> qualify_case_material_constraint_alias(const CaseMaterialAliasRow& left,
                                                          const CaseMaterialAliasRow& right,
                                                          const s::PgRigidSupportRows& support,
                                                          const s::PgReferenceLimits& limits,
                                                          core::Workspace& workspace,
                                                          core::ResourcePlan& plan) noexcept {
    constexpr auto operation = core::Operation::PrepareModel;
    return core::guard_entry(
        operation, FailureCode::PrecisionUnresolved, {}, [&]() -> core::Result<void> {
            try {
                pg::ReferenceArithmetic m{limits, plan};
                m.cancel();
                const auto a = left.surface().source, b = right.surface().source;
                if (a.source != b.source || a.contract_revision != b.contract_revision ||
                    a.representation_revision != b.representation_revision ||
                    left.coefficients().empty() || right.coefficients().empty())
                    partition_fail(FailureCode::StaleGeneration);
                const auto& chart = support.chart();
                const auto info = partition_take(chart.original().info());
                const auto identity = partition_take(chart.identity());
                const auto one_dual =
                    partition_mul(3, partition_add(chart.free_columns().size(),
                                                   static_cast<std::uint64_t>(info.rows)));
                const auto retained = partition_add(partition_mul(4, one_dual),
                                                    partition_mul(16, support.group_count()));
                m.admit(partition_add(pg::reference_live_rationals, retained));
                const auto left_dual = transpose(m, left, chart, info, identity, 0, workspace);
                const auto right_dual = transpose(m, right, chart, info, identity,
                                                  partition_mul(2, one_dual), workspace);
                const auto& g = left.radicand();
                const auto& h = right.radicand();
                const auto compare = [&](View x, View y) {
                    if (!equal(m, x, g, y, h))
                        partition_fail(FailureCode::GeometryInvalid);
                };
                const auto x_free = left_dual.rational.independent_force.values();
                const auto y_free = right_dual.rational.independent_force.values();
                const auto x_rad = left_dual.radical.independent_force.values();
                const auto y_rad = right_dual.radical.independent_force.values();
                for (std::size_t i = 0; i < x_free.size(); ++i) {
                    m.visit(1);
                    compare({x_free[i][0], x_rad[i][0]}, {y_free[i][0], y_rad[i][0]});
                }
                const auto x = prescribed(m, support, left_dual);
                const auto y = prescribed(m, support, right_dual);
                const auto compare_row = [&](const s::PgReferenceScalar& u,
                                             const s::PgReferenceScalar& v) {
                    compare({u.rational, u.radical}, {v.rational, v.radical});
                };
                for (std::size_t axis = 0; axis < x.constant.size(); ++axis)
                    compare_row(x.constant[axis], y.constant[axis]);
                for (std::size_t group = 0; group < x.motion.size(); ++group) {
                    m.visit(1);
                    const auto& l = x.motion.values()[group];
                    const auto& r = y.motion.values()[group];
                    compare_row(l.translation, r.translation);
                    for (std::size_t axis = 0; axis < l.rotation.size(); ++axis)
                        compare_row(l.rotation[axis], r.rotation[axis]);
                }
                m.cancel();
                return {};
            } catch (const pg::PartitionAbort& rejected) {
                return std::unexpected(
                    core::Failure{rejected.failure.code, operation, rejected.failure.details});
            }
        });
}

} // namespace uvlas::coupling::detail
