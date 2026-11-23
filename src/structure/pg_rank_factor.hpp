#pragma once

#include "pg_rank_input.hpp"

namespace uvlas::structure::detail {

constexpr std::array<std::uint32_t, 4> rank_primes{2147483647, 2147483629, 2147483587, 2147483579};

struct RankModTerm {
    std::uint32_t row = 0;
    std::uint32_t value = 0;
};

[[nodiscard]] inline std::uint32_t rank_mod_product(std::uint32_t a, std::uint32_t b,
                                                    std::uint32_t prime) noexcept {
    return static_cast<std::uint32_t>((static_cast<std::uint64_t>(a) * b) % prime);
}

[[nodiscard]] inline std::uint32_t rank_mod_subtract(std::uint32_t a, std::uint32_t b,
                                                     std::uint32_t prime) noexcept {
    return a >= b ? a - b : prime - (b - a);
}

[[nodiscard]] inline std::uint32_t rank_inverse(RankBudget& budget, std::uint32_t x,
                                                std::uint32_t prime) {
    if (x == 0)
        partition_fail(core::FailureCode::RankUnresolved);
    auto result = std::uint32_t{1};
    for (auto exponent = prime - 2; exponent > 0; exponent >>= 1U) {
        budget.charge(5);
        if (exponent % 2 != 0)
            result = rank_mod_product(result, x, prime);
        x = rank_mod_product(x, x, prime);
    }
    return result;
}

[[nodiscard]] inline std::uint32_t rank_mod(const RankCell& value, std::uint32_t prime,
                                            RankBudget& budget) {
    budget.charge(value.backend().size());
    const RankCell remainder = value % prime;
    return static_cast<std::uint32_t>(remainder < 0 ? RankCell{remainder + prime} : remainder);
}

struct RankFactor {
    std::size_t stride;
    core::PlannedBuffer<std::uint32_t> work, degree, row_order, buckets;
    core::PlannedBuffer<std::uint32_t> prefix, selected, pivots, inverse, row_to_pivot;
    core::PlannedBuffer<std::uint32_t> upper, lower;
    core::PlannedBuffer<std::uint64_t> offsets;
    core::PlannedBuffer<RankModTerm> reduced;
    std::uint32_t rank = 0;

    [[nodiscard]] std::span<const RankModTerm> reduced_column(std::uint32_t k) const noexcept {
        const auto start = offsets.values()[k];
        return reduced.values().subspan(static_cast<std::size_t>(start),
                                        static_cast<std::size_t>(offsets.values()[k + 1] - start));
    }
    [[nodiscard]] std::size_t square_index(std::uint32_t i, std::uint32_t j) const noexcept {
        return static_cast<std::size_t>(i) * stride + j;
    }
};

inline void count_rank_factor(RankSizer& s, const RankShape& shape) {
    const auto r = shape.rank_ceiling;
    s.add<std::uint32_t>(shape.rows * 4 + 2 * shape.candidates + 1 + 3 * r);
    s.add<std::uint32_t>(2 * partition_mul(r, r));
    s.add<std::uint64_t>(r + 1);
    s.add<RankModTerm>(partition_mul(shape.rows, r));
}

[[nodiscard]] inline RankFactor rank_factor_storage(core::ResourcePlan& plan,
                                                    const RankShape& shape) {
    const auto r = shape.rank_ceiling;
    const auto square = partition_mul(r, r);
    const auto row = [&] { return rank_allocate<std::uint32_t>(plan, shape.rows); };
    const auto basis = [&] { return rank_allocate<std::uint32_t>(plan, r); };
    return {r,
            row(),
            row(),
            row(),
            rank_allocate<std::uint32_t>(plan, shape.candidates + 1),
            rank_allocate<std::uint32_t>(plan, shape.candidates),
            basis(),
            basis(),
            basis(),
            row(),
            rank_allocate<std::uint32_t>(plan, square),
            rank_allocate<std::uint32_t>(plan, square),
            rank_allocate<std::uint64_t>(plan, r + 1),
            rank_allocate<RankModTerm>(plan, partition_mul(shape.rows, r))};
}

inline void rank_order_rows(const RankShape& shape, const RankCsc& csc, RankFactor& factor,
                            RankBudget& budget) {
    const auto used = csc.offsets.values()[shape.candidates];
    budget.charge(
        partition_add(std::uint64_t{2} * shape.rows, partition_add(used, shape.candidates + 1)));
    for (const auto& entry : csc.entries.values().first(static_cast<std::size_t>(used)))
        ++factor.degree.values()[entry.row];
    auto buckets = factor.buckets.values();
    for (const auto degree : factor.degree.values())
        ++buckets[degree];
    std::uint32_t pos = 0;
    for (auto& bucket : buckets) {
        const auto count = bucket;
        bucket = pos;
        pos += count;
    }
    for (std::uint32_t row = 0; row < shape.rows; ++row)
        factor.row_order.values()[buckets[factor.degree.values()[row]]++] = row;
}

inline void rank_discover(const RankShape& shape, const RankColumns& columns, RankFactor& factor,
                          std::uint32_t prime, RankBudget& budget) {
    for (const auto& d : columns.row_scale.values())
        if (rank_mod(d, prime, budget) == 0)
            partition_fail(core::FailureCode::RankUnresolved);
    for (const auto& e : columns.column_scale.values())
        if (rank_mod(e, prime, budget) == 0)
            partition_fail(core::FailureCode::RankUnresolved);
    auto work = factor.work.values();
    auto reduced = factor.reduced.values();
    auto offsets = factor.offsets.values();
    auto selected = factor.selected.values();
    auto pivots = factor.pivots.values();
    auto upper = factor.upper.values();
    auto lower = factor.lower.values();
    auto inverse = factor.inverse.values();
    auto row_to_pivot = factor.row_to_pivot.values();
    budget.charge(partition_add(row_to_pivot.size(), lower.size()));
    std::fill(row_to_pivot.begin(), row_to_pivot.end(), rank_missing);
    std::fill(lower.begin(), lower.end(), std::uint32_t{0});
    factor.rank = 0;
    offsets[0] = 0;
    std::uint64_t pos = 0;
    for (std::uint32_t j = 0; j < shape.candidates; ++j) {
        const auto r = factor.rank;
        budget.charge(shape.rows);
        std::fill(work.begin(), work.end(), std::uint32_t{0});
        for (const auto& term : columns.values.column(j))
            work[term.row] = rank_mod(term.value, prime, budget);
        for (std::uint32_t k = 0; k < r; ++k) {
            budget.charge();
            const auto rhs = work[pivots[k]];
            if (r < shape.rank_ceiling)
                upper[factor.square_index(k, r)] = rhs;
            if (rhs == 0)
                continue;
            const auto terms = factor.reduced_column(k);
            budget.charge(3 * terms.size());
            for (const auto& term : terms)
                work[term.row] = rank_mod_subtract(work[term.row],
                                                   rank_mod_product(rhs, term.value, prime), prime);
        }
        factor.prefix.values()[j] = r;
        budget.charge(shape.rows);
        const auto order = factor.row_order.values();
        const auto pivot = std::find_if(order.begin(), order.end(),
                                        [&](std::uint32_t row) { return work[row] != 0; });
        if (pivot == order.end())
            continue;
        if (j >= rank_required + shape.columns || r >= shape.rank_ceiling ||
            row_to_pivot[*pivot] != rank_missing)
            partition_fail(core::FailureCode::RankUnresolved);

        const auto pivot_value = work[*pivot];
        const auto inv = rank_inverse(budget, pivot_value, prime);
        if (rank_mod_product(pivot_value, inv, prime) != 1)
            partition_fail(core::FailureCode::RankUnresolved);
        inverse[r] = inv;
        upper[factor.square_index(r, r)] = pivot_value;
        selected[r] = j;
        pivots[r] = *pivot;
        row_to_pivot[*pivot] = r;
        budget.charge(shape.rows);
        for (std::uint32_t row = 0; row < shape.rows; ++row)
            if (work[row] != 0) {
                budget.charge(2);
                reduced[static_cast<std::size_t>(pos++)] = {
                    row, rank_mod_product(work[row], inv, prime)};
            }
        offsets[r + 1] = pos;
        ++factor.rank;
    }

    if (factor.rank == 0)
        partition_fail(core::FailureCode::RankUnresolved);
    for (std::uint32_t k = 0; k < factor.rank; ++k)
        for (const auto& term : factor.reduced_column(k)) {
            budget.charge();
            const auto pivot = row_to_pivot[term.row];
            if (pivot != rank_missing && pivot > k)
                lower[factor.square_index(pivot, k)] = term.value;
        }
}

inline void rank_solve(const RankFactor& factor, std::span<std::uint32_t> v, std::uint32_t prime,
                       RankBudget& budget) {
    const auto lower = factor.lower.values();
    const auto upper = factor.upper.values();
    const auto inverse = factor.inverse.values();
    for (std::uint32_t i = 0; i < factor.rank; ++i) {
        auto accum = v[i];
        budget.charge(2 * i + 1);
        for (std::uint32_t j = 0; j < i; ++j)
            accum = rank_mod_subtract(
                accum, rank_mod_product(lower[factor.square_index(i, j)], v[j], prime), prime);
        v[i] = accum;
    }
    for (auto k = factor.rank; k > 0; --k) {
        const auto i = k - 1;
        auto accum = v[i];
        budget.charge(2 * (factor.rank - k) + 2);
        for (auto j = k; j < factor.rank; ++j)
            accum = rank_mod_subtract(
                accum, rank_mod_product(upper[factor.square_index(i, j)], v[j], prime), prime);
        v[i] = rank_mod_product(accum, inverse[i], prime);
    }
}

} // namespace uvlas::structure::detail
