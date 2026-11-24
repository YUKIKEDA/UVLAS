#include "pg_reference_arithmetic.hpp"

#include <algorithm>
#include <utility>
#include <uvlas/structure/pg_global_constraints.hpp>

namespace uvlas::structure {
namespace {
using namespace detail;

void require_capacity(std::uint64_t n, std::uint64_t limit) {
    if (n > limit)
        partition_fail(core::FailureCode::ResourceExceeded);
}

std::size_t checked_index(core::Index i, std::size_t size) {
    if (i < 0 || static_cast<std::uint64_t>(i) >= size)
        partition_fail(core::FailureCode::InvalidInput);
    return static_cast<std::size_t>(i);
}

template <class Binding, class Key>
void unique_ids(std::span<const Binding> bindings, ReferenceArithmetic& math, Key key) {
    if (bindings.size() < 2)
        return;
    math.visit(bindings.size());
    using Identity = decltype(key(bindings.front()));
    auto scratch = partition_take(
        math.plan.allocate<Identity>(partition_index(bindings.size()), partition_operation));
    auto ids = scratch.values();
    for (std::size_t i = 0; i < bindings.size(); ++i) {
        math.cancel();
        ids[i] = key(bindings[i]);
    }
    const auto less = [&](const Identity& a, const Identity& b) {
        math.visit(1);
        return a < b;
    };
    std::ranges::sort(ids, less);
    for (std::size_t i = 1; i < ids.size(); ++i) {
        math.visit(1);
        if (ids[i - 1] == ids[i])
            partition_fail(core::FailureCode::InvalidInput);
    }
}

PgTraceDepth copy_depth(ReferenceArithmetic& math, const PgTraceDepth& x) {
    return {math.copy(x.lower_first), math.copy(x.upper_first), math.copy(x.lower_last),
            math.copy(x.upper_last)};
}

PgTraceOverlap copy_overlap(ReferenceArithmetic& math, const PgTraceOverlap& x) {
    return {x.left_edge,
            x.right_edge,
            x.left_layer,
            x.right_layer,
            math.copy(x.first),
            math.copy(x.last),
            math.copy(x.active_first),
            math.copy(x.active_last),
            copy_depth(math, x.left_depth),
            copy_depth(math, x.right_depth)};
}

PgConstraintSource bind_source(const PgConstraintTraceBinding& connection, const PgTraceInfo& info,
                               std::span<const PgConstraintPatchLayout> layout,
                               core::Index first_row, core::Index first_overlap) {
    const auto& left = layout[checked_index(connection.left_patch, layout.size())];
    if (left.basis != info.left || left.count != info.left_columns)
        partition_fail(core::FailureCode::StaleGeneration);
    std::optional<numerics::PatchId> right;
    if (info.kind == PgTraceKind::Connection) {
        if (!connection.right_patch || *connection.right_patch == connection.left_patch)
            partition_fail(core::FailureCode::InvalidInput);
        const auto& patch = layout[checked_index(*connection.right_patch, layout.size())];
        if (patch.basis != info.right || patch.count != info.right_columns)
            partition_fail(core::FailureCode::StaleGeneration);
        right = patch.id;
    } else if (connection.right_patch) {
        partition_fail(core::FailureCode::InvalidInput);
    }
    return {connection.id, info.kind,
            left.id,       right,
            first_row,     partition_index(info.original_rows),
            first_overlap, partition_index(info.overlaps)};
}

/// Local columns are ordered left then right. Emit the two contiguous spans in global order.
void scatter_row(ReferenceArithmetic& math, std::span<const PgTraceSparseEntry> local,
                 const ExactRational& scale, const PgConstraintPatchLayout& left,
                 const PgConstraintPatchLayout* right, std::span<PgTraceSparseEntry> output) {
    std::size_t out = 0;
    const auto copy = [&](std::span<const PgTraceSparseEntry> block, core::Index shift) {
        for (const auto& term : block) {
            math.visit(1);
            output[out++] = {term.column + shift, math.mul(scale, term.coefficient)};
        }
    };
    if (right == nullptr) {
        copy(local, left.first);
        return;
    }
    const auto boundary = std::ranges::partition_point(local, [&](const PgTraceSparseEntry& x) {
        math.visit(1);
        return x.column < left.count;
    });
    const auto split = static_cast<std::size_t>(boundary - local.begin());
    if (left.first < right->first) {
        copy(local.first(split), left.first);
        copy(local.subspan(split), right->first - left.count);
    } else {
        copy(local.subspan(split), right->first - left.count);
        copy(local.first(split), left.first);
    }
}

void accumulate(ReferenceArithmetic& math, PgCartesianExactPosition& value,
                const ExactRational& coefficient, const PgCartesianExactPosition& input) {
    for (std::size_t axis = 0; axis < value.size(); ++axis)
        value[axis] = math.add(value[axis], math.mul(coefficient, input[axis]));
}
} // namespace

struct PgGlobalConstraints::Data {
    PgGlobalConstraintInfo info{};
    core::PlannedBuffer<PgConstraintPatchLayout> layout;
    core::PlannedBuffer<PgConstraintSource> sources;
    core::PlannedBuffer<PgTraceOverlap> overlaps;
    core::PlannedBuffer<PgConstraintOriginalRow> rows;
    core::PlannedBuffer<std::size_t> offsets;
    core::PlannedBuffer<PgTraceSparseEntry> terms;

    [[nodiscard]] std::span<const PgTraceSparseEntry> row(std::size_t i) const noexcept {
        const auto ranges = offsets.values();
        return terms.values().subspan(ranges[i], ranges[i + 1] - ranges[i]);
    }
};

PgGlobalConstraints::PgGlobalConstraints(core::PlannedBuffer<Data> data) noexcept
    : data_(std::move(data)) {}

PgGlobalConstraints::PgGlobalConstraints(PgGlobalConstraints&&) noexcept = default;
PgGlobalConstraints& PgGlobalConstraints::operator=(PgGlobalConstraints&&) noexcept = default;
PgGlobalConstraints::~PgGlobalConstraints() = default;

core::Result<PgGlobalConstraints>
PgGlobalConstraints::prepare(std::span<const PgConstraintPatchBinding> patches,
                             std::span<const PgConstraintTraceBinding> bindings,
                             const PgGlobalConstraintLimits& limits, core::Workspace& workspace,
                             core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgGlobalConstraints> {
        ReferenceArithmetic math{limits.arithmetic, plan};
        math.cancel();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        require_capacity(patches.size(), limits.maximum_patches);
        require_capacity(bindings.size(), limits.maximum_sources);
        unique_ids(patches, math, [](const auto& binding) { return binding.id.value(); });
        unique_ids(bindings, math, [](const auto& binding) {
            return std::pair{partition_take(binding.trace.info()).kind, binding.id.value()};
        });

        auto owner = partition_take(plan.allocate<Data>(1, partition_operation));
        auto& data = owner.values().front();
        data.layout = partition_take(plan.allocate<PgConstraintPatchLayout>(
            partition_index(patches.size()), partition_operation));
        data.sources = partition_take(plan.allocate<PgConstraintSource>(
            partition_index(bindings.size()), partition_operation));
        std::uint64_t columns = 0, rows = 0, overlaps = 0, nonzeros = 0;
        for (std::size_t i = 0; i < patches.size(); ++i) {
            math.visit(1);
            const auto patch = partition_take(patches[i].patch.info());
            const auto count = patch.basis.rank;
            data.layout.values()[i] = {patches[i].id, patch.identity, partition_index(columns),
                                       partition_index(count)};
            columns = partition_add(columns, count);
        }
        static_cast<void>(partition_index(partition_mul(3, columns)));
        for (std::size_t i = 0; i < bindings.size(); ++i) {
            math.visit(1);
            const auto& source = bindings[i];
            const auto info = partition_take(source.trace.info());
            data.sources.values()[i] =
                bind_source(source, info, data.layout.values(), partition_index(rows),
                            partition_index(overlaps));
            rows = partition_add(rows, info.original_rows);
            overlaps = partition_add(overlaps, info.overlaps);
            require_capacity(rows, limits.maximum_rows);
            for (const auto& original : source.trace.original_rows()) {
                math.visit(1);
                if (original.unique_row < 0)
                    continue;
                const auto local = partition_take(source.trace.unique_row(original.unique_row));
                nonzeros = partition_add(nonzeros, local.size());
                require_capacity(nonzeros, limits.maximum_nonzeros);
            }
        }
        static_cast<void>(partition_index(partition_mul(3, rows)));
        const auto offset_count = partition_add(rows, 1);
        const auto retained = partition_add(
            nonzeros, partition_add(partition_mul(4, rows), partition_mul(12, overlaps)));
        math.admit(partition_add(reference_live_rationals, retained));
        data.overlaps = partition_take(
            plan.allocate<PgTraceOverlap>(partition_index(overlaps), partition_operation));
        data.rows = partition_take(
            plan.allocate<PgConstraintOriginalRow>(partition_index(rows), partition_operation));
        data.offsets = partition_take(
            plan.allocate<std::size_t>(partition_index(offset_count), partition_operation));
        data.terms = partition_take(
            plan.allocate<PgTraceSparseEntry>(partition_index(nonzeros), partition_operation));

        std::size_t overlap_at = 0, row_at = 0, term_at = 0;
        for (std::size_t i = 0; i < bindings.size(); ++i) {
            math.visit(1);
            const auto& source = bindings[i];
            const auto& left = data.layout.values()[static_cast<std::size_t>(source.left_patch)];
            const auto* right =
                source.right_patch
                    ? &data.layout.values()[static_cast<std::size_t>(*source.right_patch)]
                    : nullptr;
            for (const auto& overlap : source.trace.overlaps()) {
                math.visit(1);
                data.overlaps.values()[overlap_at++] = copy_overlap(math, overlap);
            }
            const auto originals = source.trace.original_rows();
            for (std::size_t j = 0; j < originals.size(); ++j) {
                math.visit(1);
                const auto& original = originals[j];
                auto& row = data.rows.values()[row_at];
                row = {static_cast<core::Index>(i),
                       static_cast<core::Index>(j),
                       original.overlap,
                       original.field,
                       original.power,
                       original.unique_row,
                       math.copy(original.scale),
                       {}};
                for (std::size_t axis = 0; axis < row.prescribed.size(); ++axis)
                    row.prescribed[axis] = math.copy(original.prescribed[axis]);
                if (original.unique_row >= 0) {
                    const auto local = partition_take(source.trace.unique_row(original.unique_row));
                    scatter_row(math, local, original.scale, left, right,
                                data.terms.values().subspan(term_at, local.size()));
                    term_at += local.size();
                }
                data.offsets.values()[++row_at] = term_at;
            }
        }
        math.cancel();
        data.info = {partition_index(columns), partition_index(rows), partition_index(nonzeros),
                     math.info()};
        return PgGlobalConstraints{std::move(owner)};
    });
}

core::Result<PgGlobalConstraintInfo> PgGlobalConstraints::info() const noexcept {
    if (data_.size() == 0)
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, partition_operation});
    return data_.values().front().info;
}

std::span<const PgConstraintPatchLayout> PgGlobalConstraints::layout() const& noexcept {
    return data_.size() == 0 ? std::span<const PgConstraintPatchLayout>{}
                             : data_.values().front().layout.values();
}

std::span<const PgConstraintSource> PgGlobalConstraints::sources() const& noexcept {
    return data_.size() == 0 ? std::span<const PgConstraintSource>{}
                             : data_.values().front().sources.values();
}

std::span<const PgTraceOverlap> PgGlobalConstraints::overlaps() const& noexcept {
    return data_.size() == 0 ? std::span<const PgTraceOverlap>{}
                             : data_.values().front().overlaps.values();
}

std::span<const PgConstraintOriginalRow> PgGlobalConstraints::original_rows() const& noexcept {
    return data_.size() == 0 ? std::span<const PgConstraintOriginalRow>{}
                             : data_.values().front().rows.values();
}

core::Result<std::span<const PgTraceSparseEntry>>
PgGlobalConstraints::row(core::Index i) const& noexcept {
    return partition_entry([&]() -> core::Result<std::span<const PgTraceSparseEntry>> {
        const auto index = checked_index(i, original_rows().size());
        return data_.values().front().row(index);
    });
}

core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
PgGlobalConstraints::apply_impl(PgGlobalConstraintCoordinates input, bool subtract_prescribed,
                                const PgReferenceLimits& limits, core::Workspace& workspace,
                                core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<core::PlannedBuffer<PgCartesianExactPosition>> {
        ReferenceArithmetic math{limits, plan};
        math.cancel();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        const auto current = partition_take(info());
        if (input.values.size() != static_cast<std::size_t>(current.columns) ||
            input.layout.size() != layout().size())
            partition_fail(core::FailureCode::InvalidInput);
        for (std::size_t i = 0; i < input.layout.size(); ++i) {
            math.visit(1);
            if (input.layout[i] != layout()[i])
                partition_fail(core::FailureCode::StaleGeneration);
        }
        math.admit(partition_add(reference_live_rationals,
                                 partition_mul(3, static_cast<std::uint64_t>(current.rows))));
        math.visit(static_cast<std::uint64_t>(current.rows));
        auto output = partition_take(
            plan.allocate<PgCartesianExactPosition>(current.rows, partition_operation));
        const auto& data = data_.values().front();
        for (std::size_t i = 0; i < output.size(); ++i) {
            math.visit(1);
            auto& value = output.values()[i];
            for (const auto& term : data.row(i)) {
                math.visit(1);
                accumulate(math, value, term.coefficient,
                           input.values[static_cast<std::size_t>(term.column)]);
            }
            if (subtract_prescribed)
                for (std::size_t axis = 0; axis < value.size(); ++axis)
                    value[axis] = math.sub(value[axis], data.rows.values()[i].prescribed[axis]);
        }
        math.cancel();
        return output;
    });
}

core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
PgGlobalConstraints::apply(PgGlobalConstraintCoordinates input, const PgReferenceLimits& limits,
                           core::Workspace& workspace, core::ResourcePlan& plan) const noexcept {
    return apply_impl(input, false, limits, workspace, plan);
}

core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
PgGlobalConstraints::residual(PgGlobalConstraintCoordinates input, const PgReferenceLimits& limits,
                              core::Workspace& workspace, core::ResourcePlan& plan) const noexcept {
    return apply_impl(input, true, limits, workspace, plan);
}

core::Result<core::PlannedBuffer<PgCartesianExactPosition>>
PgGlobalConstraints::transpose(std::span<const PgCartesianExactPosition> lambda,
                               const PgReferenceLimits& limits, core::Workspace& workspace,
                               core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<core::PlannedBuffer<PgCartesianExactPosition>> {
        ReferenceArithmetic math{limits, plan};
        math.cancel();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        const auto current = partition_take(info());
        if (lambda.size() != static_cast<std::size_t>(current.rows))
            partition_fail(core::FailureCode::InvalidInput);
        math.admit(partition_add(reference_live_rationals,
                                 partition_mul(3, static_cast<std::uint64_t>(current.columns))));
        math.visit(static_cast<std::uint64_t>(current.columns));
        auto output = partition_take(
            plan.allocate<PgCartesianExactPosition>(current.columns, partition_operation));
        const auto& data = data_.values().front();
        for (std::size_t i = 0; i < lambda.size(); ++i) {
            math.visit(1);
            for (const auto& term : data.row(i)) {
                math.visit(1);
                accumulate(math, output.values()[static_cast<std::size_t>(term.column)],
                           term.coefficient, lambda[i]);
            }
        }
        math.cancel();
        return output;
    });
}

} // namespace uvlas::structure
