#pragma once

#include "pg_rank_factor.hpp"

namespace uvlas::structure::detail {

struct RankLift {
    RankCsc pivot_matrix;
    core::PlannedBuffer<RankCell> current, row_sum, bound;
    core::PlannedBuffer<std::uint32_t> digit;
    core::PlannedBuffer<RankWide> residue, scaled, residual;
    core::PlannedBuffer<RankFraction<RankWide>> recovered;
};

inline void count_rank_lift(RankSizer& s, const RankShape& shape) {
    const auto r = shape.rank_ceiling;
    s.add<std::uint64_t>(r + 1);
    s.add<RankCellEntry>(std::min(shape.candidate_terms, partition_mul(r, r)));
    s.add<RankCell>(std::uint64_t{3} * r);
    s.add<std::uint32_t>(r);
    s.add<RankWide>(2 * r + shape.rows);
    s.add<RankFraction<RankWide>>(r);
}

[[nodiscard]] inline RankLift rank_lift_storage(core::ResourcePlan& plan, const RankShape& shape) {
    const auto r = shape.rank_ceiling;
    return {rank_csc(plan, r, std::min(shape.candidate_terms, partition_mul(r, r))),
            rank_allocate<RankCell>(plan, r),
            rank_allocate<RankCell>(plan, r),
            rank_allocate<RankCell>(plan, r),
            rank_allocate<std::uint32_t>(plan, r),
            rank_allocate<RankWide>(plan, r),
            rank_allocate<RankWide>(plan, r),
            rank_allocate<RankWide>(plan, shape.rows),
            rank_allocate<RankFraction<RankWide>>(plan, r)};
}

inline void rank_pivot_matrix(const RankColumns& columns, const RankFactor& factor, RankLift& lift,
                              RankBudget& budget) {
    auto sums = lift.row_sum.values();
    budget.charge(sums.size());
    std::fill(sums.begin(), sums.end(), RankCell{0});
    std::uint64_t position = 0;
    for (std::uint32_t k = 0; k < factor.rank; ++k) {
        for (const auto& term : columns.values.column(factor.selected.values()[k])) {
            budget.charge();
            const auto pivot = factor.row_to_pivot.values()[term.row];
            if (pivot == rank_missing)
                continue;
            sums[pivot] += rank_abs(term.value);
            lift.pivot_matrix.entries.values()[static_cast<std::size_t>(position++)] = {pivot,
                                                                                        term.value};
        }
        lift.pivot_matrix.offsets.values()[k + 1] = position;
    }
}

[[nodiscard]] inline bool rank_reconstruct(const RankWide& residue, const RankWide& modulus,
                                           const RankWide& bound, RankFraction<RankWide>& out,
                                           RankBudget& budget) {
    RankWide r0 = modulus, r1 = residue, t0 = 0, t1 = 1;
    while (r1 > bound) {
        budget.charge(2 * partition_mul(r0.backend().size(), r1.backend().size()));
        const RankWide quotient = r0 / r1;
        auto next_r = RankWide{r0 % r1};
        auto next_t = RankWide{t0 - rank_product(budget, quotient, t1)};
        r0 = std::move(r1);
        r1 = std::move(next_r);
        t0 = std::move(t1);
        t1 = std::move(next_t);
    }
    if (t1 == 0 || rank_abs(t1) > bound)
        return false;
    const RankWide check = r1 - rank_product(budget, t1, residue);
    budget.charge(partition_mul(check.backend().size(), modulus.backend().size()));
    if (check % modulus != 0)
        return false;
    out = normalized_rank_fraction(budget, std::move(r1), std::move(t1));
    return true;
}

[[nodiscard]] inline RankWide rank_proof_scale(const RankFactor& factor, RankLift& lift,
                                               RankBudget& budget) {
    auto common = RankWide{1};
    const auto recovered = lift.recovered.values().first(factor.rank);
    for (const auto& pair : recovered) {
        budget.charge();
        common = rank_lcm(budget, common, pair.denominator);
    }
    for (std::uint32_t k = 0; k < factor.rank; ++k) {
        budget.charge(
            partition_mul(common.backend().size(), recovered[k].denominator.backend().size()));
        const RankWide scale = common / recovered[k].denominator;
        lift.scaled.values()[k] = rank_product(budget, recovered[k].numerator, scale);
    }
    return common;
}

[[nodiscard]] inline bool rank_prove(const RankColumns& columns, const RankFactor& factor,
                                     RankLift& lift, std::uint32_t j, const RankWide& common,
                                     std::span<RankWide> residual, bool pivot_only,
                                     RankBudget& budget) {
    budget.charge(residual.size());
    std::fill(residual.begin(), residual.end(), RankWide{0});
    for (std::uint32_t k = 0; k < factor.rank; ++k) {
        budget.charge();
        const auto& coordinate = lift.scaled.values()[k];
        if (coordinate == 0)
            continue;
        const auto terms = pivot_only ? lift.pivot_matrix.column(k)
                                      : columns.values.column(factor.selected.values()[k]);
        for (const auto& term : terms)
            residual[term.row] += rank_product(budget, coordinate, RankWide{term.value});
    }
    for (const auto& term : columns.values.column(j)) {
        budget.charge();
        const auto row = pivot_only ? factor.row_to_pivot.values()[term.row] : term.row;
        if (row != rank_missing)
            residual[row] -= rank_product(budget, common, RankWide{term.value});
    }
    budget.charge(residual.size());
    return std::all_of(residual.begin(), residual.end(), [](const RankWide& x) { return x == 0; });
}

inline void rank_lift_column(const RankColumns& columns, const RankFactor& factor, RankLift& lift,
                             std::uint32_t j, std::uint32_t prime, RankBudget& budget) {
    const auto r = factor.rank;
    auto current = lift.current.values().first(r);
    auto sizes = lift.bound.values().first(r);
    auto residue = lift.residue.values().first(r);
    auto digit = lift.digit.values().first(r);
    auto recovered = lift.recovered.values().first(r);
    budget.charge(std::uint64_t{3} * r);
    std::fill(current.begin(), current.end(), RankCell{0});
    std::fill(residue.begin(), residue.end(), RankWide{0});
    for (const auto& term : columns.values.column(j)) {
        budget.charge();
        const auto row = factor.row_to_pivot.values()[term.row];
        if (row != rank_missing)
            current[row] = term.value;
    }
    for (std::uint32_t row = 0; row < r; ++row)
        sizes[row] = std::max(rank_abs(current[row]), lift.row_sum.values()[row]);

    auto modulus = RankWide{1};
    for (std::uint32_t step = 1; step <= budget.limits.maximum_digits; ++step) {
        for (std::uint32_t row = 0; row < r; ++row)
            digit[row] = rank_mod(current[row], prime, budget);
        rank_solve(factor, digit, prime, budget);
        for (std::uint32_t k = 0; k < r; ++k) {
            budget.charge();
            if (digit[k] == 0)
                continue;
            for (const auto& term : lift.pivot_matrix.column(k))
                current[term.row] -= rank_product(budget, term.value, RankCell{digit[k]});
            residue[k] += rank_product(budget, modulus, RankWide{digit[k]});
        }
        modulus = rank_product(budget, modulus, RankWide{prime});
        for (std::uint32_t row = 0; row < r; ++row) {
            budget.charge(2 * current[row].backend().size());
            if (current[row] % prime != 0)
                partition_fail(core::FailureCode::RankUnresolved);
            current[row] /= prime;
            if (rank_abs(current[row]) > sizes[row])
                partition_fail(core::FailureCode::RankUnresolved);
        }
        if (!std::has_single_bit(step))
            continue;

        budget.charge(partition_mul(modulus.backend().size(), modulus.backend().size()));
        const RankWide bound = boost::multiprecision::sqrt(RankWide{modulus / 2});
        bool ready = true;
        for (std::uint32_t k = 0; k < r && ready; ++k) {
            ready = rank_reconstruct(residue[k], modulus, bound, recovered[k], budget);
            if (ready && k >= factor.prefix.values()[j])
                ready = recovered[k].numerator == 0;
        }
        if (!ready)
            continue;
        const auto common = rank_proof_scale(factor, lift, budget);
        auto all_rows = lift.residual.values();
        if (rank_prove(columns, factor, lift, j, common, all_rows.first(r), true, budget) &&
            rank_prove(columns, factor, lift, j, common, all_rows, false, budget))
            return;
    }
    partition_fail(core::FailureCode::RankUnresolved);
}

} // namespace uvlas::structure::detail
