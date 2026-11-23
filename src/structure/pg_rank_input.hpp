#pragma once

#include "pg_rank_storage.hpp"

namespace uvlas::structure::detail {

struct RankShape {
    std::uint32_t cells_axis = 0;
    std::uint32_t cells = 0;
    std::uint32_t columns = 0;
    std::uint32_t originals = 0;
    std::uint32_t candidates = 0;
    std::uint32_t rows = 0;
    std::uint32_t trace_rows = 0;
    std::uint32_t rank_ceiling = 0;
    std::uint64_t activity_terms = 0;
    std::uint64_t candidate_terms = 0;
    std::uint64_t map_terms = 0;
    std::array<bool, rank_max_cells> active{};
    std::array<std::uint32_t, rank_key_count> row_ids{};
    std::array<std::uint32_t, rank_max_columns> raw_counts{};
    std::array<std::uint32_t, rank_max_candidates> map_counts{};
};

template <class F>
void rank_map_terms(const PgCoefficientDictionary& dictionary, const RankShape& shape, F&& apply) {
    for (const auto& entry : partition_take(dictionary.required_map()))
        apply(entry.target, RankInputMap{static_cast<std::uint32_t>(entry.source.value), false,
                                         &entry.coefficient});
    for (std::uint32_t raw = 0; raw < shape.columns; ++raw)
        apply(rank_required + raw, RankInputMap{raw});
    const auto anchor_start = rank_required + shape.columns;
    for (std::uint32_t raw = 0; raw < shape.originals; ++raw)
        apply(anchor_start + raw, RankInputMap{raw});
    const auto columns = partition_take(dictionary.columns());
    for (std::uint32_t raw = shape.originals; raw < shape.columns; ++raw) {
        const auto& source = columns[raw];
        if (source.kind != PgCoefficientKind::AnchorR && source.kind != PgCoefficientKind::AnchorP)
            partition_fail(core::FailureCode::RankUnresolved);
        const auto original = static_cast<std::uint32_t>(source.original.value);
        if (original >= shape.originals)
            partition_fail(core::FailureCode::RankUnresolved);
        apply(anchor_start + original, RankInputMap{raw, true, nullptr});
    }
    for (const auto& entry : partition_take(dictionary.anchor_c_map())) {
        const auto source = static_cast<std::uint32_t>(entry.source.value);
        const auto target = static_cast<std::uint32_t>(entry.target.value);
        if (source >= shape.originals || target >= shape.columns)
            partition_fail(core::FailureCode::RankUnresolved);
        apply(anchor_start + source, RankInputMap{target, true, &entry.coefficient});
    }
}

struct RankKey {
    std::uint32_t activity;
    std::uint32_t trace;
    std::uint32_t cell;
    bool anchor;
};

[[nodiscard]] inline RankKey rank_key(const RankShape& shape, const PgCoefficientRecord& record,
                                      const ExactPowers& powers) {
    const auto field = static_cast<std::size_t>(record.field);
    if (field >= rank_degree.size() || record.cell.first < 0 || record.cell.second < 0 ||
        record.cell.first >= shape.cells_axis || record.cell.second >= shape.cells_axis ||
        powers.first < 0 || powers.second < 0 || powers.first > rank_degree[field] ||
        powers.second > rank_degree[field] - powers.first)
        partition_fail(core::FailureCode::RankUnresolved);
    const auto degree = rank_degree[field];
    const auto a = static_cast<std::uint32_t>(powers.first);
    const auto b = static_cast<std::uint32_t>(powers.second);
    const auto cell =
        static_cast<std::uint32_t>(record.cell.first * shape.cells_axis + record.cell.second);
    const auto slot = a * (degree + 1) - a * (a == 0 ? 0 : a - 1) / 2 + b;
    const auto active = rank_field_slots[field] * shape.cells +
                        (rank_field_slots[field + 1] - rank_field_slots[field]) * cell + slot;
    const auto trace = rank_trace_slots[field] * shape.cells_axis +
                       (degree + 1) * static_cast<std::uint32_t>(record.cell.second) + b;
    return {active, trace, cell, record.cell.first == 0};
}

[[nodiscard]] inline RankShape count_rank_input(const PgCoefficientDictionary& dictionary,
                                                const PgIntegrationPartition& partition,
                                                RankBudget& budget) {
    const auto graph = partition_take(dictionary.graph());
    if (!partition_same_graph(graph, partition_take(partition.graph()),
                              [&](std::uint64_t units) { budget.charge(units); }))
        partition_fail(core::FailureCode::StaleGeneration);
    RankShape result;
    result.cells_axis = static_cast<std::uint32_t>(2 * graph.n);
    result.cells = result.cells_axis * result.cells_axis;
    result.trace_rows = rank_trace_slots.back() * result.cells_axis;
    const auto columns = partition_take(dictionary.columns());
    const auto originals = partition_take(dictionary.original_columns());
    const auto activities = partition_take(partition.activity_map());
    if ((graph.n != 1 && graph.n != 2 && graph.n != 4) || activities.size() != result.cells ||
        columns.empty() || columns.size() > rank_max_columns || originals <= 0 ||
        originals > rank_max_original || static_cast<std::uint64_t>(originals) > columns.size())
        partition_fail(core::FailureCode::InvalidInput);
    result.columns = static_cast<std::uint32_t>(columns.size());
    result.originals = static_cast<std::uint32_t>(originals);
    result.candidates = rank_required + result.columns + result.originals;
    for (std::uint32_t cell = 0; cell < result.cells; ++cell) {
        budget.charge();
        const auto& activity = activities[cell];
        const auto expected = PgCellIndex{cell / result.cells_axis, cell % result.cells_axis};
        if (activity.cell != expected || activity.area.sign() < 0 ||
            activity.active != (activity.area.sign() > 0))
            partition_fail(core::FailureCode::RankUnresolved);
        result.active[cell] = activity.active;
    }
    result.row_ids.fill(rank_missing);
    for (const auto& record : partition_take(dictionary.records())) {
        const auto raw = static_cast<std::uint32_t>(record.column.value);
        if (raw >= result.columns)
            partition_fail(core::FailureCode::RankUnresolved);
        for (const auto& term : record.polynomial.terms()) {
            budget.charge();
            const auto key = rank_key(result, record, term.powers);
            if (result.active[key.cell]) {
                result.row_ids[key.activity] = 0;
                ++result.raw_counts[raw];
                ++result.activity_terms;
            }
        }
    }
    for (auto& id : result.row_ids)
        if (id != rank_missing)
            id = result.rows++;
    if (result.rows == 0)
        partition_fail(core::FailureCode::RankUnresolved);
    result.rank_ceiling = std::min(result.rows, result.columns);
    std::array<std::uint64_t, rank_max_candidates> upper{};
    rank_map_terms(dictionary, result, [&](std::uint32_t j, const RankInputMap& term) {
        budget.charge();
        if (j >= result.candidates || term.source >= result.columns)
            partition_fail(core::FailureCode::RankUnresolved);
        upper[j] += result.raw_counts[term.source];
        ++result.map_counts[j];
        ++result.map_terms;
    });
    for (std::uint32_t j = 0; j < result.candidates; ++j)
        result.candidate_terms += std::min<std::uint64_t>(result.rows, upper[j]);
    return result;
}

inline void fill_rank_rows(const RankShape& shape, std::span<PgRankRow> rows) {
    for (std::size_t field = 0; field < rank_degree.size(); ++field) {
        const auto degree = rank_degree[field];
        for (std::uint32_t cell = 0; cell < shape.cells; ++cell) {
            PgCoefficientRecord record;
            record.cell = {cell / shape.cells_axis, cell % shape.cells_axis};
            record.field = static_cast<PgCoefficientField>(field);
            for (std::uint32_t a = 0; a <= degree; ++a)
                for (std::uint32_t b = 0; b <= degree - a; ++b) {
                    const auto id = shape.row_ids[rank_key(shape, record, {a, b}).activity];
                    if (id != rank_missing)
                        rows[id] = {record.field, record.cell, {a, b}};
                }
        }
    }
}

struct RankColumns {
    RankCsc values;
    core::PlannedBuffer<RankCell> row_scale, column_scale;
    core::PlannedBuffer<std::uint64_t> map_offsets;
    core::PlannedBuffer<RankInputMap> map_entries;

    [[nodiscard]] std::span<const RankInputMap> map(std::uint32_t j) const noexcept {
        const auto offsets = map_offsets.values();
        return map_entries.values().subspan(static_cast<std::size_t>(offsets[j]),
                                            static_cast<std::size_t>(offsets[j + 1] - offsets[j]));
    }
};

inline void count_rank_columns(RankSizer& s, const RankShape& shape) {
    s.add<std::uint64_t>(shape.columns + 1);
    s.add<std::uint64_t>((shape.candidates + std::uint64_t{1}) * 2);
    s.add<RankCellEntry>(shape.activity_terms);
    s.add<RankCellEntry>(shape.candidate_terms);
    s.add<RankCell>(shape.rows * 2 + shape.candidates + shape.trace_rows);
    s.add<RankCell>(partition_mul(shape.columns, shape.trace_rows));
    s.add<RankInputMap>(shape.map_terms);
}

[[nodiscard]] inline RankColumns assemble_rank_columns(const PgCoefficientDictionary& dictionary,
                                                       const RankShape& shape, RankBudget& budget,
                                                       core::ResourcePlan& plan) {
    auto columns = RankColumns{rank_csc(plan, shape.candidates, shape.candidate_terms),
                               rank_allocate<RankCell>(plan, shape.rows),
                               rank_allocate<RankCell>(plan, shape.candidates),
                               rank_allocate<std::uint64_t>(plan, shape.candidates + 1),
                               rank_allocate<RankInputMap>(plan, shape.map_terms)};
    const auto records = partition_take(dictionary.records());
    auto raw = rank_csc(plan, shape.columns, shape.activity_terms);
    auto anchor_scale = rank_allocate<RankCell>(plan, shape.trace_rows);
    auto trace = rank_allocate<RankCell>(plan, partition_mul(shape.columns, shape.trace_rows));
    auto work = rank_allocate<RankCell>(plan, shape.rows);
    std::fill(columns.row_scale.values().begin(), columns.row_scale.values().end(), RankCell{1});
    std::fill(anchor_scale.values().begin(), anchor_scale.values().end(), RankCell{1});

    for (const auto& record : records)
        for (const auto& term : record.polynomial.terms()) {
            budget.charge();
            const auto key = rank_key(shape, record, term.powers);
            const auto id = shape.row_ids[key.activity];
            if (id == rank_missing && !key.anchor)
                continue;
            const auto coefficient = rank_fraction<RankCell>(term.coefficient);
            if (id != rank_missing) {
                auto& scale = columns.row_scale.values()[id];
                scale = rank_lcm(budget, scale, coefficient.denominator);
            }
            if (key.anchor) {
                auto& scale = anchor_scale.values()[key.trace];
                scale = rank_lcm(budget, scale, coefficient.denominator);
            }
        }

    std::array<std::uint64_t, rank_max_columns> raw_next{};
    for (std::uint32_t j = 0; j < shape.columns; ++j) {
        raw_next[j] = raw.offsets.values()[j];
        raw.offsets.values()[j + 1] = raw_next[j] + shape.raw_counts[j];
    }
    for (const auto& record : records)
        for (const auto& term : record.polynomial.terms()) {
            budget.charge();
            const auto key = rank_key(shape, record, term.powers);
            const auto id = shape.row_ids[key.activity];
            if (id == rank_missing && !key.anchor)
                continue;
            const auto j = static_cast<std::uint32_t>(record.column.value);
            const auto value = rank_fraction<RankCell>(term.coefficient);
            if (id != rank_missing) {
                const RankCell multiple = columns.row_scale.values()[id] / value.denominator;
                raw.entries.values()[static_cast<std::size_t>(raw_next[j]++)] = {
                    id, rank_product(budget, multiple, value.numerator)};
            }
            if (key.anchor) {
                const RankCell multiple = anchor_scale.values()[key.trace] / value.denominator;
                auto contribution = rank_product(budget, multiple, value.numerator);
                if (term.powers.first % 2 != 0)
                    contribution = -contribution;
                trace.values()[static_cast<std::size_t>(j) * shape.trace_rows + key.trace] +=
                    contribution;
            }
        }

    std::array<std::uint64_t, rank_max_candidates> map_next{};
    for (std::uint32_t j = 0; j < shape.candidates; ++j) {
        map_next[j] = columns.map_offsets.values()[j];
        columns.map_offsets.values()[j + 1] = map_next[j] + shape.map_counts[j];
    }
    rank_map_terms(dictionary, shape, [&](std::uint32_t j, const RankInputMap& term) {
        budget.charge();
        columns.map_entries.values()[static_cast<std::size_t>(map_next[j]++)] = term;
    });

    std::array<RankCell, rank_trace_count> anchor_residual{};
    std::uint64_t position = 0;
    for (std::uint32_t j = 0; j < shape.candidates; ++j) {
        budget.charge(shape.rows);
        std::fill(work.values().begin(), work.values().end(), RankCell{0});
        const bool anchor = j >= rank_required + shape.columns;
        auto scale = RankCell{1};
        for (const auto& entry : columns.map(j))
            if (entry.coefficient != nullptr)
                scale = rank_lcm(budget, scale,
                                 rank_fraction<RankCell>(*entry.coefficient).denominator);
        columns.column_scale.values()[j] = scale;
        for (const auto& entry : columns.map(j)) {
            budget.charge();
            auto coefficient = entry.coefficient != nullptr
                                   ? rank_fraction<RankCell>(*entry.coefficient)
                                   : RankFraction<RankCell>{RankCell{1}, RankCell{1}};
            const RankCell quotient = scale / coefficient.denominator;
            auto weight = rank_product(budget, quotient, coefficient.numerator);
            if (entry.negative)
                weight = -weight;
            for (const auto& term : raw.column(entry.source))
                work.values()[term.row] += rank_product(budget, term.value, weight);
            if (anchor) {
                const auto terms = trace.values().subspan(
                    static_cast<std::size_t>(entry.source) * shape.trace_rows, shape.trace_rows);
                for (std::uint32_t i = 0; i < shape.trace_rows; ++i)
                    anchor_residual[i] += rank_product(budget, terms[i], weight);
            }
        }
        if (anchor) {
            budget.charge(shape.trace_rows);
            for (std::uint32_t i = 0; i < shape.trace_rows; ++i)
                if (anchor_residual[i] != 0)
                    partition_fail(core::FailureCode::RankUnresolved);
        }
        budget.charge(shape.rows);
        for (std::uint32_t i = 0; i < shape.rows; ++i)
            if (work.values()[i] != 0) {
                if (position == shape.candidate_terms)
                    partition_fail(core::FailureCode::RankUnresolved);
                columns.values.entries.values()[static_cast<std::size_t>(position++)] = {
                    i, std::move(work.values()[i])};
            }
        columns.values.offsets.values()[j + 1] = position;
    }
    return columns;
}

} // namespace uvlas::structure::detail
