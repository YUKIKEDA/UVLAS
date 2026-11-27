#include "pg_reference_hash.hpp"

#include <algorithm>
#include <bit>
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
                     const Position& value, std::uint8_t components = 7) {
    for (std::size_t k = 0; k < from.size(); ++k)
        if ((components & (1U << k)) != 0)
            from[k] = math.sub(from[k], math.mul(scale, value[k]));
}

void divide_by(ReferenceArithmetic& math, Position& value, const ExactRational& divisor,
               std::uint8_t components) {
    for (std::size_t k = 0; k < value.size(); ++k)
        if ((components & (1U << k)) != 0)
            value[k] = math.div(value[k], divisor);
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

    [[nodiscard]] std::array<std::size_t, 2>
    factor(ReferenceArithmetic&, const PgGlobalConstraints&, std::uint64_t columns,
           std::uint64_t pivots_capacity, std::uint64_t work_capacity, std::uint64_t terms_capacity,
           std::uint8_t components);

    template <class Rhs>
    void solve_pivots(ReferenceArithmetic& math, std::span<Position> q,
                      std::span<Position> transformed, const Rhs& rhs,
                      std::uint8_t components) const {
        // C_selected = L*U. First solve the unit triangular equations for the RHS.
        for (std::size_t i = 0; i < rank; ++i) {
            math.visit(1);
            const auto& value = rhs(static_cast<std::size_t>(pivots.values()[i].original_row));
            for (std::size_t axis = 0; axis < value.size(); ++axis)
                if ((components & (1U << axis)) != 0)
                    transformed[i][axis] = math.copy(value[axis]);
            for (const auto& term : lower(i)) {
                math.visit(1);
                subtract_scaled(math, transformed[i], term.coefficient,
                                transformed[static_cast<std::size_t>(term.column)], components);
            }
        }
        // U may couple to the free coordinates, which have already been put in q.
        for (std::size_t i = rank; i > 0;) {
            --i;
            math.visit(1);
            auto& value = transformed[i];
            const auto row = upper(i);
            for (const auto& term : row.subspan(1)) {
                math.visit(1);
                subtract_scaled(math, value, term.coefficient, q[global(term.column)], components);
            }
            divide_by(math, value, row.front().coefficient, components);
            auto& q_pivot = q[global(row.front().column)];
            for (std::size_t axis = 0; axis < value.size(); ++axis)
                if ((components & (1U << axis)) != 0)
                    q_pivot[axis] = std::move(value[axis]);
        }
    }

    void transpose_pivots(ReferenceArithmetic& math, std::span<Position> work,
                          std::span<Position> multipliers, std::uint8_t components) const {
        // U^T, scattering to the still unsolved pivots and to the free coordinates.
        for (std::size_t i = 0; i < rank; ++i) {
            math.visit(1);
            const auto active = upper(i);
            auto& pivot = work[global(active.front().column)];
            divide_by(math, pivot, active.front().coefficient, components);
            for (const auto& term : active.subspan(1)) {
                math.visit(1);
                subtract_scaled(math, work[global(term.column)], term.coefficient, pivot,
                                components);
            }
        }
        // L^T followed by injection into the original row numbering.
        for (std::size_t i = rank; i > 0;) {
            --i;
            math.visit(1);
            const auto pivot = pivots.values()[i];
            auto& value = work[static_cast<std::size_t>(pivot.dependent_column)];
            for (const auto& term : lower(i)) {
                math.visit(1);
                const auto earlier = pivots.values()[static_cast<std::size_t>(term.column)];
                subtract_scaled(math, work[static_cast<std::size_t>(earlier.dependent_column)],
                                term.coefficient, value, components);
            }
            auto& original = multipliers[static_cast<std::size_t>(pivot.original_row)];
            for (std::size_t axis = 0; axis < value.size(); ++axis)
                if ((components & (1U << axis)) != 0)
                    original[axis] = std::move(value[axis]);
        }
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
        solve_pivots(math, q.values(), transformed.values(), rhs, 7);
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

std::array<std::size_t, 2> PgSparseChart::Data::factor(ReferenceArithmetic& math,
                                                       const PgGlobalConstraints& original,
                                                       std::uint64_t columns, std::uint64_t pc,
                                                       std::uint64_t wc, std::uint64_t fc,
                                                       std::uint8_t components) {
    auto& data = *this;
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
        if ((block.world_components & components) == 0)
            continue;
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
            if (data.rank == data.pivots.size() || size + lower_size > data.terms.size() - retained)
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
    return {retained_lower, retained - retained_lower};
}

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
        if (!info.common_cartesian)
            partition_fail(core::FailureCode::Unsupported);
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
        const auto [retained_lower, retained_upper] =
            data.factor(math, original, columns, pc, wc, fc, 7);
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
                       retained_upper,
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
        data.transpose_pivots(math, work.values(), multipliers.values(), 7);
        for (std::size_t i = 0; i < reduced.size(); ++i) {
            math.visit(1);
            reduced.values()[i] =
                std::move(work.values()[static_cast<std::size_t>(data.free.values()[i])]);
        }
        math.cancel();
        return PgSparseChartDual{std::move(reduced), std::move(multipliers), math.info()};
    });
}

struct PgComponentSparseChart::Data {
    PgSparseChartIdentity identity;
    PgComponentSparsePreparation report;
    std::array<PgSparseChart::Data, 3> factors;
    std::array<std::uint8_t, 3> components{};
    std::size_t groups{}, maximum_rank{};
    core::PlannedBuffer<core::Index> free;

    [[nodiscard]] core::PlannedBuffer<Position>
    solve(ReferenceArithmetic& math, const PgGlobalConstraints& original, std::uint64_t columns,
          std::span<const ExactRational> independent, bool homogeneous) const {
        auto q = allocate<Position>(math, columns);
        auto transformed = allocate<Position>(math, maximum_rank);
        for (std::size_t i = 0; i < independent.size(); ++i) {
            math.visit(1);
            const auto flat = static_cast<std::size_t>(free.values()[i]);
            q.values()[flat / 3][flat % 3] = math.copy(independent[i]);
        }
        const auto rows = original.original_rows();
        const auto sources = original.sources();
        const Position zero{};
        const auto rhs = [&](std::size_t i) -> const Position& {
            return homogeneous ? zero : rows[i].prescribed;
        };
        for (std::size_t i = 0; i < groups; ++i)
            factors[i].solve_pivots(math, q.values(), transformed.values(), rhs, components[i]);
        for (std::size_t i = 0; i < rows.size(); ++i) {
            math.visit(1);
            const auto mask = sources[static_cast<std::size_t>(rows[i].source)].world_components;
            auto residual = copy_position(math, rhs(i));
            for (const auto& term : partition_take(original.row(partition_index(i)))) {
                math.visit(1);
                subtract_scaled(math, residual, term.coefficient,
                                q.values()[static_cast<std::size_t>(term.column)], mask);
            }
            if (std::ranges::any_of(residual, [](const auto& value) { return value.sign() != 0; }))
                partition_fail(core::FailureCode::InvalidInput);
        }
        math.cancel();
        return q;
    }
};

PgComponentSparseChart::PgComponentSparseChart(PgGlobalConstraints&& original,
                                               core::PlannedBuffer<Data> data) noexcept
    : original_(std::move(original)), data_(std::move(data)) {}

PgComponentSparseChart::PgComponentSparseChart(PgComponentSparseChart&&) noexcept = default;
PgComponentSparseChart&
PgComponentSparseChart::operator=(PgComponentSparseChart&&) noexcept = default;
PgComponentSparseChart::~PgComponentSparseChart() = default;

core::Result<PgComponentSparseChart>
PgComponentSparseChart::prepare(PgGlobalConstraints&& original, const PgSparseChartLimits& limits,
                                core::Workspace& workspace, core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgComponentSparseChart> {
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
        std::array<std::uint8_t, 3> components{};
        std::size_t groups = 0;
        const auto sources = original.sources();
        for (unsigned axis = 0; axis < 3; ++axis) {
            const auto bit = static_cast<std::uint8_t>(1U << axis);
            auto group = std::span{components}.first(groups);
            const auto matching = std::ranges::find_if(group, [&](std::uint8_t previous) {
                return std::ranges::all_of(sources, [&](const auto& source) {
                    math.visit(1);
                    const auto mask = source.world_components;
                    return ((mask & bit) != 0) == ((mask & previous) != 0);
                });
            });
            if (matching == group.end())
                components[groups++] = bit;
            else
                *matching |= bit;
        }
        std::array<std::uint64_t, 3> active_rows{};
        std::uint64_t active_factors = 0;
        for (std::size_t i = 0; i < groups; ++i) {
            for (const auto& source : sources) {
                math.visit(1);
                if ((source.world_components & components[i]) != 0)
                    active_rows[i] =
                        partition_add(active_rows[i], static_cast<std::uint64_t>(source.rows));
            }
            if (active_rows[i] != 0)
                ++active_factors;
        }
        const auto factor_live = partition_mul(active_factors, fc);
        const auto row_scratch = partition_add(partition_mul(2, wc), pc);
        const auto lift_scratch = partition_mul(3, partition_add(columns, pc));
        const auto live = partition_add(factor_live, std::max(row_scratch, lift_scratch));
        math.admit(partition_add(reference_live_rationals, live));
        auto owner = allocate<Data>(math, 1);
        auto& data = owner.values().front();
        data.components = components;
        data.groups = groups;
        auto free_mask = allocate<std::uint8_t>(math, columns);
        std::uint64_t total_free = 0, lower_terms = 0, upper_terms = 0;
        for (std::size_t i = 0; i < groups; ++i) {
            auto& factor = data.factors[i];
            const auto [lower, upper] = factor.factor(
                math, original, columns, std::min(pc, active_rows[i]), active_rows[i] == 0 ? 0 : wc,
                active_rows[i] == 0 ? 0 : fc, components[i]);
            lower_terms = partition_add(lower_terms, lower);
            upper_terms = partition_add(upper_terms, upper);
            data.maximum_rank = std::max(data.maximum_rank, factor.rank);
            const auto axes = static_cast<std::uint64_t>(std::popcount(components[i]));
            total_free = partition_add(total_free, partition_mul(axes, factor.free.size()));
            for (const auto column : factor.free.values()) {
                math.visit(1);
                free_mask.values()[static_cast<std::size_t>(column)] |= components[i];
            }
        }
        data.free = allocate<core::Index>(math, total_free);
        std::size_t next = 0;
        for (std::size_t i = 0; i < free_mask.size(); ++i) {
            math.visit(1);
            for (unsigned axis = 0; axis < 3; ++axis)
                if ((free_mask.values()[i] & (1U << axis)) != 0)
                    data.free.values()[next++] = partition_index(3 * i + axis);
        }
        static_cast<void>(data.solve(math, original, columns, {}, false));
        ReferenceHash hash{math};
        constexpr std::string_view version = "uvlas.pg.component-sparse-chart.v1";
        hash.bytes(version.data(), version.size());
        const auto original_id = chart_identity(math, original, {});
        hash.bytes(original_id.bytes.data(), original_id.bytes.size());
        for (const auto& source : sources) {
            math.visit(1);
            hash.word(source.world_components);
        }
        hash.word(groups);
        for (std::size_t i = 0; i < groups; ++i) {
            const auto& factor = data.factors[i];
            hash.word(data.components[i]);
            hash.word(factor.rank);
            for (const auto& pivot : factor.pivots.values().first(factor.rank)) {
                math.visit(1);
                hash.word(static_cast<std::uint64_t>(pivot.original_row));
                hash.word(static_cast<std::uint64_t>(pivot.dependent_column));
            }
        }
        data.identity = {hash.result().bytes};
        std::optional<core::StageTiming> timing;
        if (started) {
            const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                     std::chrono::steady_clock::now() - *started)
                                     .count();
            timing = core::StageTiming{1, static_cast<std::uint64_t>(elapsed)};
        }
        const auto rank = partition_mul(3, columns) - total_free;
        data.report = {active_factors, rank, lower_terms, upper_terms, math.info(), timing};
        math.cancel();
        return PgComponentSparseChart{std::move(original), std::move(owner)};
    });
}

core::Result<PgSparseChartIdentity> PgComponentSparseChart::identity() const noexcept {
    if (data_.size() == 0)
        return std::unexpected(
            core::Failure{core::FailureCode::StaleGeneration, partition_operation});
    return data_.values().front().identity;
}

core::Result<PgComponentSparsePreparation> PgComponentSparseChart::preparation() const noexcept {
    if (data_.size() == 0)
        return std::unexpected(
            core::Failure{core::FailureCode::StaleGeneration, partition_operation});
    return data_.values().front().report;
}

const PgGlobalConstraints& PgComponentSparseChart::original() const& noexcept {
    return original_;
}

std::span<const core::Index> PgComponentSparseChart::free_columns() const& noexcept {
    return data_.size() == 0 ? std::span<const core::Index>{}
                             : data_.values().front().free.values();
}

core::Result<core::PlannedBuffer<Position>>
PgComponentSparseChart::lift_impl(PgComponentSparseInput input, bool homogeneous,
                                  const PgReferenceLimits& limits, core::Workspace& workspace,
                                  core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<core::PlannedBuffer<Position>> {
        ReferenceArithmetic math{limits, plan};
        math.cancel();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        if (input.identity != partition_take(identity()))
            partition_fail(core::FailureCode::StaleGeneration);
        const auto& data = data_.values().front();
        if (input.values.size() != data.free.size())
            partition_fail(core::FailureCode::InvalidInput);
        const auto columns = static_cast<std::uint64_t>(partition_take(original_.info()).columns);
        const auto live = partition_mul(3, partition_add(columns, data.maximum_rank));
        math.admit(partition_add(reference_live_rationals, live));
        return data.solve(math, original_, columns, input.values, homogeneous);
    });
}

core::Result<core::PlannedBuffer<Position>>
PgComponentSparseChart::lift(PgComponentSparseInput input, const PgReferenceLimits& limits,
                             core::Workspace& workspace, core::ResourcePlan& plan) const noexcept {
    return lift_impl(input, false, limits, workspace, plan);
}

core::Result<core::PlannedBuffer<Position>>
PgComponentSparseChart::tangent(PgComponentSparseInput input, const PgReferenceLimits& limits,
                                core::Workspace& workspace,
                                core::ResourcePlan& plan) const noexcept {
    return lift_impl(input, true, limits, workspace, plan);
}

core::Result<PgComponentSparseDual>
PgComponentSparseChart::transpose(PgSparseChartInput input, const PgReferenceLimits& limits,
                                  core::Workspace& workspace,
                                  core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<PgComponentSparseDual> {
        ReferenceArithmetic math{limits, plan};
        math.cancel();
        const auto lease = partition_take(workspace.acquire());
        static_cast<void>(lease);
        if (input.identity != partition_take(identity()))
            partition_fail(core::FailureCode::StaleGeneration);
        const auto info = partition_take(original_.info());
        const auto columns = static_cast<std::uint64_t>(info.columns);
        const auto rows = static_cast<std::uint64_t>(info.rows);
        if (input.values.size() != columns)
            partition_fail(core::FailureCode::InvalidInput);
        const auto& data = data_.values().front();
        const auto live =
            partition_add(partition_mul(3, partition_add(columns, rows)), data.free.size());
        math.admit(partition_add(reference_live_rationals, live));
        auto work = allocate<Position>(math, columns);
        auto multipliers = allocate<Position>(math, rows);
        auto reduced = allocate<ExactRational>(math, data.free.size());
        for (std::size_t i = 0; i < input.values.size(); ++i) {
            math.visit(1);
            work.values()[i] = copy_position(math, input.values[i]);
        }
        for (std::size_t i = 0; i < data.groups; ++i)
            data.factors[i].transpose_pivots(math, work.values(), multipliers.values(),
                                             data.components[i]);
        for (std::size_t i = 0; i < data.free.size(); ++i) {
            math.visit(1);
            const auto flat = static_cast<std::size_t>(data.free.values()[i]);
            reduced.values()[i] = std::move(work.values()[flat / 3][flat % 3]);
        }
        math.cancel();
        return PgComponentSparseDual{std::move(reduced), std::move(multipliers), math.info()};
    });
}

} // namespace uvlas::structure
