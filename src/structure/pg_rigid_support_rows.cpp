#include "pg_reference_arithmetic.hpp"

#include <algorithm>
#include <type_traits>
#include <uvlas/structure/pg_rigid_support_rows.hpp>

namespace uvlas::structure {
namespace {
using namespace detail;
using Position = PgCartesianExactPosition;

template <class T>
auto allocate(ReferenceArithmetic& math, std::uint64_t n) {
    return partition_take(math.plan.allocate<T>(partition_index(n), partition_operation));
}

void bind_patches(ReferenceArithmetic& math, std::span<const PgConstraintPatchLayout> layout,
                  std::span<const PgConstraintPatchBinding> bindings,
                  std::span<const PgPreparedPatch*> output) {
    if (bindings.size() != layout.size())
        partition_fail(core::FailureCode::InvalidInput);
    for (const auto& binding : bindings) {
        const auto match = std::ranges::find_if(layout, [&](const auto& x) {
            math.visit(1);
            return x.id == binding.id;
        });
        if (match == layout.end())
            partition_fail(core::FailureCode::InvalidInput);
        auto& slot = output[static_cast<std::size_t>(match - layout.begin())];
        if (slot != nullptr)
            partition_fail(core::FailureCode::InvalidInput);
        const auto info = partition_take(binding.patch.info());
        const auto count = static_cast<std::size_t>(match->count);
        if (info.identity != match->basis ||
            binding.patch.reference_coordinates().size() != count ||
            binding.patch.translation().size() != count)
            partition_fail(core::FailureCode::StaleGeneration);
        slot = &binding.patch;
    }
}

template <class Binding>
void bind_supports(ReferenceArithmetic& math, std::span<const PgConstraintSource> sources,
                   std::span<const Binding> moving, std::span<bool> used,
                   std::span<std::optional<std::size_t>> output) {
    if (moving.size() > sources.size())
        partition_fail(core::FailureCode::InvalidInput);
    for (const auto& binding : moving) {
        const auto entry = [&] {
            if constexpr (std::is_same_v<Binding, PgConstraintId>)
                return PgRigidSupportBinding{binding, 0};
            else
                return binding;
        }();
        if (entry.group >= used.size())
            partition_fail(core::FailureCode::InvalidInput);
        const auto match = std::ranges::find_if(sources, [&](const auto& x) {
            math.visit(1);
            return x.kind == PgTraceKind::Support && x.id == entry.support;
        });
        if (match == sources.end())
            partition_fail(core::FailureCode::InvalidInput);
        auto& selected = output[static_cast<std::size_t>(match - sources.begin())];
        if (selected)
            partition_fail(core::FailureCode::InvalidInput);
        selected = entry.group;
        used[entry.group] = true;
    }
    for (const auto seen : used) {
        math.visit(1);
        if (!seen)
            partition_fail(core::FailureCode::InvalidInput);
    }
}

void fill_rows(ReferenceArithmetic& math, const PgGlobalConstraints& original,
               std::span<const PgPreparedPatch* const> patches,
               std::span<const std::optional<std::size_t>> groups,
               std::span<PgRigidSupportRow> result) {
    const auto sources = original.sources();
    const auto layout = original.layout();
    for (std::size_t i = 0; i < sources.size(); ++i) {
        math.visit(1);
        if (!groups[i])
            continue;
        const auto& source = sources[i];
        const auto it = std::ranges::find_if(layout, [&](const auto& x) {
            math.visit(1);
            return x.id == source.left;
        });
        if (it == layout.end())
            partition_fail(core::FailureCode::InvalidInput);
        const auto& patch = *patches[static_cast<std::size_t>(it - layout.begin())];
        const auto xyz = patch.reference_coordinates();
        const auto b = patch.translation();
        const auto first_row = static_cast<std::size_t>(source.first_row);
        const auto count = static_cast<std::size_t>(source.rows);
        for (std::size_t local = 0; local < count; ++local) {
            math.visit(1);
            const auto row = first_row + local;
            auto& output = result[row];
            output.group = groups[i];
            for (const auto& term : partition_take(original.row(partition_index(row)))) {
                math.visit(1);
                const auto column = term.column - it->first;
                if (column < 0 || column >= it->count)
                    partition_fail(core::FailureCode::InvalidInput);
                const auto j = static_cast<std::size_t>(column);
                for (std::size_t axis = 0; axis < 3; ++axis)
                    output.reference_position[axis] = math.add(
                        output.reference_position[axis], math.mul(term.coefficient, xyz[j][axis]));
                output.translation = math.add(output.translation, math.mul(term.coefficient, b[j]));
            }
        }
    }
}

enum class Generator : std::uint8_t { Rotation, Translation, Fixed };

void verify_generators(ReferenceArithmetic& math, const PgSparseChart& chart,
                       std::span<const PgRigidSupportRow> map, std::size_t groups,
                       PgReferenceLimits checking, core::Workspace& workspace) {
    checking.cancellation = math.limits.cancellation;
    math.reserve_nested_rationals(checking.algebra, partition_mul(7, map.size()));
    const auto identity = partition_take(chart.identity());
    const auto independent = allocate<Position>(math, chart.free_columns().size());
    const auto original = chart.original().original_rows();
    auto rhs = allocate<Position>(math, map.size());
    // C acts identically in each Cartesian lane. Test the three columns of X together.
    const auto check = [&](Generator generator, std::optional<std::size_t> group) {
        for (std::size_t i = 0; i < map.size(); ++i) {
            math.visit(1);
            auto& value = rhs.values()[i];
            value = {};
            if (map[i].group != group)
                continue;
            if (generator == Generator::Translation) {
                value[0] = math.copy(map[i].translation);
            } else {
                const auto& source = generator == Generator::Rotation ? map[i].reference_position
                                                                      : original[i].prescribed;
                for (std::size_t axis = 0; axis < 3; ++axis)
                    value[axis] = math.copy(source[axis]);
            }
        }
        static_cast<void>(partition_take(chart.lift_rhs({identity, independent.values()},
                                                        {identity, rhs.values()}, checking,
                                                        workspace, math.plan)));
    };
    check(Generator::Fixed, std::nullopt);
    for (std::size_t group = 0; group < groups; ++group) {
        math.visit(1);
        check(Generator::Rotation, group);
        check(Generator::Translation, group);
    }
}

template <class Binding>
core::PlannedBuffer<PgRigidSupportRow>
prepare_rows(const PgSparseChart& chart, std::span<const PgConstraintPatchBinding> bindings,
             std::span<const Binding> moving, std::size_t groups,
             const PgRigidSupportLimits& limits, core::Workspace& workspace,
             core::ResourcePlan& plan) {
    ReferenceArithmetic math{limits.arithmetic, plan};
    math.cancel();
    if (groups > moving.size())
        partition_fail(core::FailureCode::InvalidInput);
    static_cast<void>(partition_take(chart.identity()));
    const auto& original = chart.original();
    const auto rows = original.original_rows().size();
    if (rows > limits.maximum_original_rows)
        partition_fail(core::FailureCode::ResourceExceeded);
    math.admit(partition_add(reference_live_rationals, partition_mul(7, rows)));
    core::PlannedBuffer<PgRigidSupportRow> result;
    {
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        const auto layout = original.layout();
        auto patches = allocate<const PgPreparedPatch*>(math, layout.size());
        bind_patches(math, layout, bindings, patches.values());
        const auto sources = original.sources();
        auto selected = allocate<std::optional<std::size_t>>(math, sources.size());
        auto used = allocate<bool>(math, groups);
        bind_supports(math, sources, moving, used.values(), selected.values());
        result = allocate<PgRigidSupportRow>(math, rows);
        fill_rows(math, original, patches.values(), selected.values(), result.values());
    }
    verify_generators(math, chart, result.values(), groups, limits.chart_verification, workspace);
    math.cancel();
    return result;
}
} // namespace

PgRigidSupportRows::PgRigidSupportRows(PgSparseChart&& chart,
                                       core::PlannedBuffer<PgRigidSupportRow> rows,
                                       std::size_t group_count) noexcept
    : chart_(std::move(chart)), rows_(std::move(rows)), group_count_(group_count) {}

core::Result<PgRigidSupportRows> PgRigidSupportRows::prepare(
    PgSparseChart&& chart, std::span<const PgConstraintPatchBinding> bindings,
    std::span<const PgConstraintId> moving, const PgRigidSupportLimits& limits,
    core::Workspace& workspace, core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgRigidSupportRows> {
        const std::size_t groups = moving.empty() ? 0 : 1;
        auto rows = prepare_rows(chart, bindings, moving, groups, limits, workspace, plan);
        return PgRigidSupportRows{std::move(chart), std::move(rows), groups};
    });
}

core::Result<PgRigidSupportRows>
PgRigidSupportRows::prepare_groups(PgSparseChart&& chart,
                                   std::span<const PgConstraintPatchBinding> bindings,
                                   std::span<const PgRigidSupportBinding> moving,
                                   std::size_t groups, const PgRigidSupportLimits& limits,
                                   core::Workspace& workspace, core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgRigidSupportRows> {
        auto rows = prepare_rows(chart, bindings, moving, groups, limits, workspace, plan);
        return PgRigidSupportRows{std::move(chart), std::move(rows), groups};
    });
}

const PgSparseChart& PgRigidSupportRows::chart() const& noexcept {
    return chart_;
}
std::span<const PgRigidSupportRow> PgRigidSupportRows::rows() const& noexcept {
    return rows_.values();
}

} // namespace uvlas::structure
