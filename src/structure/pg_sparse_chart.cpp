#include "pg_reference_hash.hpp"

#include <algorithm>
#include <chrono>
#include <uvlas/structure/pg_sparse_chart.hpp>

namespace uvlas::structure {
namespace {
using namespace detail;
using Position = PgCartesianExactPosition;
using Term = PgTraceSparseEntry;

struct FactorRow {
    std::size_t start{}, upper_begin{}, end{};
};

void subtract_scaled(ReferenceArithmetic& math, Position& from, const ExactRational& scale,
                     const Position& value) {
    for (std::size_t k = 0; k < from.size(); ++k)
        from[k] = math.sub(from[k], math.mul(scale, value[k]));
}

void divide_by(ReferenceArithmetic& math, Position& value, const ExactRational& divisor) {
    for (auto& component : value)
        component = math.div(component, divisor);
}

Position copy_position(ReferenceArithmetic& math, const Position& value) {
    return {math.copy(value[0]), math.copy(value[1]), math.copy(value[2])};
}

std::size_t merge_eliminated(ReferenceArithmetic& math, std::span<const Term> active,
                             std::span<const Term> upper, const ExactRational& multiplier,
                             std::span<Term> output) {
    std::size_t a = 0, b = 0, size = 0;
    while (a < active.size() || b < upper.size()) {
        math.visit(1);
        const bool use_a =
            b == upper.size() || (a < active.size() && active[a].column < upper[b].column);
        const bool use_b =
            a == active.size() || (b < upper.size() && upper[b].column < active[a].column);
        const auto column = use_b ? upper[b].column : active[a].column;
        ExactRational value;
        if (use_a) {
            value = math.copy(active[a++].coefficient);
        } else if (use_b) {
            value = math.neg(math.mul(multiplier, upper[b++].coefficient));
        } else {
            value = math.sub(active[a++].coefficient, math.mul(multiplier, upper[b++].coefficient));
        }
        if (value.sign() == 0)
            continue;
        if (size == output.size())
            partition_fail(core::FailureCode::ResourceExceeded);
        output[size++] = {column, std::move(value)};
    }
    return size;
}

template <class Less>
void ordered_indices(ReferenceArithmetic& math, std::span<core::Index> result, Less less) {
    for (std::size_t i = 0; i < result.size(); ++i) {
        math.visit(1);
        result[i] = static_cast<core::Index>(i);
    }
    std::ranges::sort(result, [&](core::Index a, core::Index b) {
        math.visit(1);
        return less(static_cast<std::size_t>(a), static_cast<std::size_t>(b));
    });
}

PgSparseChartIdentity chart_identity(ReferenceArithmetic& math, const PgGlobalConstraints& original,
                                     std::span<const PgSparseChartPivot> pivots) {
    ReferenceHash hash{math};
    constexpr std::string_view version = "uvlas.pg.sparse-chart.support-first.v1";
    hash.bytes(version.data(), version.size());
    const auto integer = [&](auto x) { hash.word(static_cast<std::uint64_t>(x)); };
    const auto& blocks = original.layout();
    integer(blocks.size());
    for (const auto& patch : blocks) {
        math.visit(1);
        hash.word(patch.id.value());
        hash.bytes(patch.basis.bytes.data(), patch.basis.bytes.size());
        integer(patch.first);
        integer(patch.count);
    }
    const auto sources = original.sources();
    integer(sources.size());
    for (const auto& source : sources) {
        math.visit(1);
        hash.word(source.id.value());
        integer(source.kind);
        hash.word(source.left.value());
        integer(source.right.has_value());
        if (source.right)
            hash.word(source.right->value());
        integer(source.first_row);
        integer(source.rows);
        integer(source.first_overlap);
        integer(source.overlaps);
    }
    const auto overlaps = original.overlaps();
    integer(overlaps.size());
    for (const auto& span : overlaps) {
        math.visit(1);
        integer(span.left_edge);
        integer(span.right_edge);
        integer(span.left_layer.has_value());
        if (span.left_layer)
            integer(*span.left_layer);
        integer(span.right_layer.has_value());
        if (span.right_layer)
            integer(*span.right_layer);
        for (const auto* scalar : {&span.first, &span.last, &span.active_first, &span.active_last,
                                   &span.left_depth.lower_first, &span.left_depth.upper_first,
                                   &span.left_depth.lower_last, &span.left_depth.upper_last,
                                   &span.right_depth.lower_first, &span.right_depth.upper_first,
                                   &span.right_depth.lower_last, &span.right_depth.upper_last})
            hash.rational(*scalar);
    }
    const auto rows = original.original_rows();
    integer(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
        math.visit(1);
        const auto& row = rows[i];
        integer(row.source);
        integer(row.source_row);
        integer(row.source_overlap);
        integer(row.field);
        integer(row.power);
        for (const auto& value : row.prescribed)
            hash.rational(value);
        const auto terms = partition_take(original.row(static_cast<core::Index>(i)));
        integer(terms.size());
        for (const auto& term : terms) {
            math.visit(1);
            integer(term.column);
            hash.rational(term.coefficient);
        }
    }
    integer(pivots.size());
    for (const auto& pivot : pivots) {
        math.visit(1);
        integer(pivot.original_row);
        integer(pivot.dependent_column);
    }
    return {hash.result().bytes};
}

template <class T>
auto allocate(ReferenceArithmetic& math, std::uint64_t n) {
    math.visit(n);
    return partition_take(math.plan.allocate<T>(partition_index(n), partition_operation));
}
} // namespace

struct PgSparseChart::Data {
    PgSparseChartIdentity identity;
    PgSparseChartPreparation report{};
    core::PlannedBuffer<PgSparseChartPivot> pivots;
    core::PlannedBuffer<FactorRow> ranges;
    core::PlannedBuffer<Term> terms;
    core::PlannedBuffer<core::Index> stable_to_global;
    core::PlannedBuffer<core::Index> free;
    std::size_t rank = 0, free_size = 0;

    [[nodiscard]] std::span<const Term> lower(std::size_t i) const noexcept {
        const auto at = ranges.values()[i];
        return terms.values().subspan(at.start, at.upper_begin - at.start);
    }

    [[nodiscard]] std::span<const Term> upper(std::size_t i) const noexcept {
        const auto at = ranges.values()[i];
        return terms.values().subspan(at.upper_begin, at.end - at.upper_begin);
    }

    [[nodiscard]] std::size_t global(core::Index stable) const noexcept {
        return static_cast<std::size_t>(
            stable_to_global.values()[static_cast<std::size_t>(stable)]);
    }

    template <class Rhs>
    [[nodiscard]] core::PlannedBuffer<Position>
    solve(ReferenceArithmetic& math, const PgGlobalConstraints& original,
          std::span<const Position> independent, Rhs rhs) const {
        auto q = allocate<Position>(math, stable_to_global.size());
        auto transformed = allocate<Position>(math, rank);
        for (std::size_t j = 0; j < independent.size(); ++j) {
            math.visit(1);
            q.values()[static_cast<std::size_t>(free.values()[j])] =
                copy_position(math, independent[j]);
        }
        // C_selected = L*U. First solve the unit triangular equations for the RHS.
        for (std::size_t i = 0; i < rank; ++i) {
            math.visit(1);
            transformed.values()[i] =
                copy_position(math, rhs(static_cast<std::size_t>(pivots.values()[i].original_row)));
            for (const auto& term : lower(i)) {
                math.visit(1);
                subtract_scaled(math, transformed.values()[i], term.coefficient,
                                transformed.values()[static_cast<std::size_t>(term.column)]);
            }
        }
        // U may couple to the free coordinates, which have already been put in q.
        for (std::size_t i = rank; i > 0;) {
            --i;
            math.visit(1);
            auto& value = transformed.values()[i];
            const auto row = upper(i);
            for (const auto& term : row.subspan(1)) {
                math.visit(1);
                subtract_scaled(math, value, term.coefficient, q.values()[global(term.column)]);
            }
            divide_by(math, value, row.front().coefficient);
            q.values()[global(row.front().column)] = std::move(value);
        }
        // Check all original rows, including every skipped duplicate or dependent equation.
        const auto rows = original.original_rows();
        for (std::size_t i = 0; i < rows.size(); ++i) {
            math.visit(1);
            auto residual = copy_position(math, rhs(i));
            const auto row = partition_take(original.row(static_cast<core::Index>(i)));
            for (const auto& term : row) {
                math.visit(1);
                subtract_scaled(math, residual, term.coefficient,
                                q.values()[static_cast<std::size_t>(term.column)]);
            }
            if (std::ranges::any_of(residual, [](const auto& x) { return x.sign() != 0; }))
                partition_fail(core::FailureCode::InvalidInput);
        }
        math.cancel();
        return q;
    }
};

PgSparseChart::PgSparseChart(PgGlobalConstraints&& original,
                             core::PlannedBuffer<Data> data) noexcept
    : original_(std::move(original)), data_(std::move(data)) {}

PgSparseChart::PgSparseChart(PgSparseChart&&) noexcept = default;
PgSparseChart& PgSparseChart::operator=(PgSparseChart&&) noexcept = default;
PgSparseChart::~PgSparseChart() = default;

core::Result<PgSparseChart> PgSparseChart::prepare(PgGlobalConstraints&& original,
                                                   const PgSparseChartLimits& limits,
                                                   core::Workspace& workspace,
                                                   core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgSparseChart> {
        ReferenceArithmetic math{limits.arithmetic, plan};
        math.cancel();
        std::optional<std::chrono::steady_clock::time_point> started;
        if (limits.measure_preparation)
            started = std::chrono::steady_clock::now();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        const auto info = partition_take(original.info());
        const auto columns = static_cast<std::uint64_t>(info.columns);
        const auto rows = static_cast<std::uint64_t>(info.rows);
        const auto pc = std::min({limits.maximum_pivots, columns, rows});
        const auto wc = std::min(limits.maximum_row_terms, columns);
        const auto fc = limits.maximum_factor_terms;
        const auto live =
            partition_add(partition_add(fc, partition_mul(2, wc)),
                          partition_add(pc, partition_mul(3, partition_add(columns, pc))));
        math.admit(partition_add(reference_live_rationals, live));

        auto owner = allocate<Data>(math, 1);
        auto& data = owner.values().front();
        data.pivots = allocate<PgSparseChartPivot>(math, pc);
        data.ranges = allocate<FactorRow>(math, pc);
        data.terms = allocate<Term>(math, fc);
        data.stable_to_global = allocate<core::Index>(math, columns);
        auto global_to_stable = allocate<core::Index>(math, columns);
        auto pivot_lookup = allocate<core::Index>(math, columns);
        for (auto& slot : pivot_lookup.values()) {
            math.visit(1);
            slot = -1;
        }
        const auto patches = original.layout();
        auto patch_order = allocate<core::Index>(math, patches.size());
        ordered_indices(math, patch_order.values(), [&](std::size_t a, std::size_t b) {
            return patches[a].id.value() < patches[b].id.value();
        });
        std::size_t stable = 0;
        for (const auto block : patch_order.values()) {
            const auto& patch = patches[static_cast<std::size_t>(block)];
            for (core::Index j = 0; j < patch.count; ++j) {
                math.visit(1);
                const auto global = patch.first + j;
                data.stable_to_global.values()[stable] = global;
                global_to_stable.values()[static_cast<std::size_t>(global)] =
                    static_cast<core::Index>(stable++);
            }
        }
        const auto sources = original.sources();
        auto order = allocate<core::Index>(math, sources.size());
        ordered_indices(math, order.values(), [&](std::size_t a, std::size_t b) {
            const auto& x = sources[a];
            const auto& y = sources[b];
            return std::pair{x.kind != PgTraceKind::Support, x.id.value()} <
                   std::pair{y.kind != PgTraceKind::Support, y.id.value()};
        });

        auto first = allocate<Term>(math, wc);
        auto second = allocate<Term>(math, wc);
        auto lower_row = allocate<Term>(math, pc);
        std::size_t retained = 0, retained_lower = 0;
        for (const auto source : order.values()) {
            math.visit(1);
            const auto& block = sources[static_cast<std::size_t>(source)];
            for (core::Index row = block.first_row; row < block.first_row + block.rows; ++row) {
                math.visit(1);
                const auto raw = partition_take(original.row(row));
                if (raw.size() > wc)
                    partition_fail(core::FailureCode::ResourceExceeded);
                auto active = first.values();
                auto scratch = second.values();
                auto size = raw.size();
                for (std::size_t i = 0; i < size; ++i) {
                    math.visit(1);
                    const auto index = static_cast<std::size_t>(raw[i].column);
                    active[i] = {global_to_stable.values()[index], math.copy(raw[i].coefficient)};
                }
                std::ranges::sort(active.first(size), [&](const Term& x, const Term& y) {
                    math.visit(1);
                    return x.column < y.column;
                });

                std::size_t lower_size = 0;
                while (size != 0) {
                    // Earlier U rows never contain an earlier pivot. Eliminate the
                    // least existing pivot that is actually present, from the sparse row.
                    core::Index earliest = static_cast<core::Index>(data.rank);
                    const Term* selected = nullptr;
                    for (const auto& term : active.first(size)) {
                        math.visit(1);
                        const auto owner_index =
                            pivot_lookup.values()[static_cast<std::size_t>(term.column)];
                        if (owner_index >= 0 && owner_index < earliest) {
                            earliest = owner_index;
                            selected = &term;
                        }
                    }
                    if (selected == nullptr)
                        break;
                    const auto upper = data.upper(static_cast<std::size_t>(earliest));
                    auto multiplier = math.div(selected->coefficient, upper.front().coefficient);
                    size = merge_eliminated(math, active.first(size), upper, multiplier, scratch);
                    lower_row.values()[lower_size++] = {earliest, std::move(multiplier)};
                    std::swap(active, scratch);
                }
                if (size == 0)
                    continue;
                if (data.rank == data.pivots.size() ||
                    size + lower_size > data.terms.size() - retained)
                    partition_fail(core::FailureCode::ResourceExceeded);
                const auto pivot = active.front().column;
                const auto global = data.global(pivot);
                pivot_lookup.values()[static_cast<std::size_t>(pivot)] =
                    static_cast<core::Index>(data.rank);
                data.pivots.values()[data.rank] = {row, static_cast<core::Index>(global)};
                const auto start = retained;
                for (auto& term : lower_row.values().first(lower_size)) {
                    math.visit(1);
                    data.terms.values()[retained++] = std::move(term);
                }
                const auto upper_begin = retained;
                retained_lower += lower_size;
                for (auto& term : active.first(size)) {
                    math.visit(1);
                    data.terms.values()[retained++] = std::move(term);
                }
                data.ranges.values()[data.rank++] = {start, upper_begin, retained};
            }
        }
        data.free = allocate<core::Index>(math, columns - data.rank);
        for (std::size_t i = 0; i < pivot_lookup.size(); ++i) {
            math.visit(1);
            if (pivot_lookup.values()[i] < 0)
                data.free.values()[data.free_size++] = data.stable_to_global.values()[i];
        }
        const auto rows_view = original.original_rows();
        static_cast<void>(data.solve(math, original, {}, [&](std::size_t i) -> const Position& {
            return rows_view[i].prescribed;
        }));
        data.identity = chart_identity(math, original, data.pivots.values().first(data.rank));
        std::optional<core::StageTiming> timing;
        if (started) {
            const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                     std::chrono::steady_clock::now() - *started)
                                     .count();
            timing = core::StageTiming{1, static_cast<std::uint64_t>(elapsed)};
        }
        data.report = {data.rank,
                       static_cast<std::uint64_t>(info.nonzeros),
                       retained_lower,
                       retained - retained_lower,
                       pc,
                       wc,
                       fc,
                       math.info(),
                       timing};
        math.cancel();
        return PgSparseChart{std::move(original), std::move(owner)};
    });
}

core::Result<PgSparseChartIdentity> PgSparseChart::identity() const noexcept {
    if (data_.size() == 0)
        return std::unexpected(
            core::Failure{core::FailureCode::StaleGeneration, partition_operation});
    return data_.values().front().identity;
}

core::Result<PgSparseChartPreparation> PgSparseChart::preparation() const noexcept {
    if (data_.size() == 0)
        return std::unexpected(
            core::Failure{core::FailureCode::StaleGeneration, partition_operation});
    return data_.values().front().report;
}

std::span<const PgSparseChartPivot> PgSparseChart::pivots() const& noexcept {
    if (data_.size() == 0)
        return {};
    const auto& data = data_.values().front();
    return data.pivots.values().first(data.rank);
}

std::span<const core::Index> PgSparseChart::free_columns() const& noexcept {
    return data_.size() == 0 ? std::span<const core::Index>{}
                             : data_.values().front().free.values();
}

const PgGlobalConstraints& PgSparseChart::original() const& noexcept {
    return original_;
}

core::Result<core::PlannedBuffer<Position>>
PgSparseChart::lift_impl(PgSparseChartInput r, PgSparseChartInput rhs, RhsKind kind,
                         const PgReferenceLimits& limits, core::Workspace& workspace,
                         core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<core::PlannedBuffer<Position>> {
        ReferenceArithmetic math{limits, plan};
        math.cancel();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        if (data_.size() == 0)
            partition_fail(core::FailureCode::StaleGeneration);
        const auto info = partition_take(original_.info());
        const auto& data = data_.values().front();
        if (r.identity != data.identity ||
            (kind == RhsKind::Explicit && rhs.identity != data.identity))
            partition_fail(core::FailureCode::StaleGeneration);
        if (r.values.size() != data.free_size ||
            (kind == RhsKind::Explicit && rhs.values.size() != original_.original_rows().size()))
            partition_fail(core::FailureCode::InvalidInput);
        const auto work =
            partition_mul(3, partition_add(static_cast<std::uint64_t>(info.columns), data.rank));
        math.admit(partition_add(reference_live_rationals, work));
        const Position zero{};
        const auto original = original_.original_rows();
        return data.solve(math, original_, r.values, [&](std::size_t i) -> const Position& {
            if (kind == RhsKind::Fixed)
                return original[i].prescribed;
            if (kind == RhsKind::Explicit)
                return rhs.values[i];
            return zero;
        });
    });
}

core::Result<core::PlannedBuffer<Position>>
PgSparseChart::lift(PgSparseChartInput r, const PgReferenceLimits& limits,
                    core::Workspace& workspace, core::ResourcePlan& plan) const noexcept {
    return lift_impl(r, {}, RhsKind::Fixed, limits, workspace, plan);
}

core::Result<core::PlannedBuffer<Position>>
PgSparseChart::tangent(PgSparseChartInput r, const PgReferenceLimits& limits,
                       core::Workspace& workspace, core::ResourcePlan& plan) const noexcept {
    return lift_impl(r, {}, RhsKind::Zero, limits, workspace, plan);
}

core::Result<core::PlannedBuffer<Position>>
PgSparseChart::lift_rhs(PgSparseChartInput r, PgSparseChartInput rhs,
                        const PgReferenceLimits& limits, core::Workspace& workspace,
                        core::ResourcePlan& plan) const noexcept {
    return lift_impl(r, rhs, RhsKind::Explicit, limits, workspace, plan);
}

core::Result<PgSparseChartDual> PgSparseChart::transpose(PgSparseChartInput input,
                                                         const PgReferenceLimits& limits,
                                                         core::Workspace& workspace,
                                                         core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<PgSparseChartDual> {
        ReferenceArithmetic math{limits, plan};
        math.cancel();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        if (data_.size() == 0)
            partition_fail(core::FailureCode::StaleGeneration);
        const auto info = partition_take(original_.info());
        const auto& data = data_.values().front();
        if (input.identity != data.identity)
            partition_fail(core::FailureCode::StaleGeneration);
        const auto force = input.values;
        if (force.size() != data.stable_to_global.size())
            partition_fail(core::FailureCode::InvalidInput);
        const auto live =
            partition_add(partition_mul(2, force.size()), static_cast<std::uint64_t>(info.rows));
        math.admit(partition_add(reference_live_rationals, partition_mul(3, live)));
        auto work = allocate<Position>(math, force.size());
        for (std::size_t i = 0; i < force.size(); ++i) {
            math.visit(1);
            work.values()[i] = copy_position(math, force[i]);
        }
        auto reduced = allocate<Position>(math, data.free.size());
        auto multipliers = allocate<Position>(math, static_cast<std::uint64_t>(info.rows));
        // U^T, scattering to the still unsolved pivots and to the free coordinates.
        for (std::size_t i = 0; i < data.rank; ++i) {
            math.visit(1);
            const auto upper = data.upper(i);
            auto& pivot = work.values()[data.global(upper.front().column)];
            divide_by(math, pivot, upper.front().coefficient);
            for (const auto& term : upper.subspan(1)) {
                math.visit(1);
                subtract_scaled(math, work.values()[data.global(term.column)], term.coefficient,
                                pivot);
            }
        }
        // L^T followed by injection into the original row numbering.
        for (std::size_t i = data.rank; i > 0;) {
            --i;
            math.visit(1);
            const auto pivot = data.pivots.values()[i];
            auto& value = work.values()[static_cast<std::size_t>(pivot.dependent_column)];
            for (const auto& term : data.lower(i)) {
                math.visit(1);
                const auto earlier = data.pivots.values()[static_cast<std::size_t>(term.column)];
                subtract_scaled(math,
                                work.values()[static_cast<std::size_t>(earlier.dependent_column)],
                                term.coefficient, value);
            }
            multipliers.values()[static_cast<std::size_t>(pivot.original_row)] = std::move(value);
        }
        for (std::size_t i = 0; i < reduced.size(); ++i) {
            math.visit(1);
            reduced.values()[i] =
                std::move(work.values()[static_cast<std::size_t>(data.free.values()[i])]);
        }
        math.cancel();
        return PgSparseChartDual{std::move(reduced), std::move(multipliers)};
    });
}

} // namespace uvlas::structure
