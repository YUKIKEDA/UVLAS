#include "pg_rank_lift.hpp"
#include "pg_rank_witness.hpp"

#include <limits>
#include <utility>

namespace uvlas::structure::detail {
struct PgRankState {
    PgRankPreparationInfo info;
    RankWitness witness;
    std::uint32_t columns = 0;
    std::uint32_t originals = 0;
    core::PlannedBuffer<PgRankRow> rows;
    core::PlannedBuffer<std::uint32_t> selected;
    core::PlannedBuffer<std::uint64_t> basis_offsets, source_offsets;
    core::PlannedBuffer<std::uint32_t> basis_rows;
    core::PlannedBuffer<RankFraction<RankCell>> basis_values;
    core::PlannedBuffer<RankSourceMap> source_values;
    core::PlannedBuffer<RankFraction<RankWide>> coordinates;
};
} // namespace uvlas::structure::detail

namespace uvlas::structure {
namespace {
using namespace detail;

inline void count_rank_state(RankSizer& size, const RankShape& shape) {
    size.add<PgRankState>(1);
    size.add<PgRankRow>(shape.rows);
    size.add<std::uint32_t>(partition_add(shape.rank_ceiling, shape.candidate_terms));
    size.add<std::uint64_t>(std::uint64_t{2} * (shape.rank_ceiling + 1));
    size.add<RankFraction<RankCell>>(shape.candidate_terms);
    size.add<RankSourceMap>(shape.map_terms);
    size.add<RankFraction<RankWide>>(partition_mul(shape.rank_ceiling, shape.candidates));
}

[[nodiscard]] auto rank_state_storage(core::ResourcePlan& plan, const RankShape& shape) {
    auto state = rank_allocate<PgRankState>(plan, 1);
    auto& out = state.values().front();
    out.rows = rank_allocate<PgRankRow>(plan, shape.rows);
    out.selected = rank_allocate<std::uint32_t>(plan, shape.rank_ceiling);
    out.coordinates = rank_allocate<RankFraction<RankWide>>(
        plan, partition_mul(shape.rank_ceiling, shape.candidates));
    out.columns = shape.columns;
    out.originals = shape.originals;
    fill_rank_rows(shape, out.rows.values());
    return state;
}

void rank_proof_and_coordinates(const RankShape& shape, const RankColumns& columns,
                                const RankFactor& factor, RankLift& lift, PgRankState& out,
                                std::uint32_t prime, RankBudget& budget) {
    const auto r = factor.rank;
    const auto selected = factor.selected.values().first(r);
    std::copy(selected.begin(), selected.end(), out.selected.values().begin());
    rank_pivot_matrix(columns, factor, lift, budget);
    for (std::uint32_t j = 0; j < shape.candidates; ++j) {
        const auto begin = static_cast<std::size_t>(j) * r;
        auto to = out.coordinates.values().subspan(begin, r);
        const auto index = factor.prefix.values()[j];
        budget.charge(r);
        if (index < r && selected[index] == j) {
            // This column is the actual basis column. Check the same full-row
            // identity through the general source proof, including zeros.
            auto scaled = lift.scaled.values().first(r);
            std::fill(scaled.begin(), scaled.end(), RankWide{0});
            scaled[index] = 1;
            if (!rank_prove(columns, factor, lift, j, RankWide{1}, lift.residual.values(), false,
                            budget))
                partition_fail(core::FailureCode::RankUnresolved);
            std::fill(to.begin(), to.end(), RankFraction<RankWide>{});
            to[index].numerator = 1;
            continue;
        }
        rank_lift_column(columns, factor, lift, j, prime, budget);
        const auto from = lift.recovered.values();
        for (std::uint32_t k = 0; k < r; ++k) {
            budget.charge();
            if (from[k].numerator == 0) {
                to[k] = {};
                continue;
            }
            // M_b y = M_j, therefore C_b (e_b/e_j)y = C_j.
            auto n = rank_product(budget, from[k].numerator,
                                  RankWide{columns.column_scale.values()[selected[k]]});
            auto d = rank_product(budget, from[k].denominator,
                                  RankWide{columns.column_scale.values()[j]});
            to[k] = normalized_rank_fraction(budget, std::move(n), std::move(d));
        }
    }
}

void rank_own_basis(const RankShape& shape, const RankColumns& columns, const RankFactor& factor,
                    PgRankState& out, RankBudget& budget, core::ResourcePlan& plan) {
    const auto selected = factor.selected.values().first(factor.rank);
    std::uint64_t terms = 0, mapped = 0;
    for (const auto j : selected) {
        budget.charge();
        terms = partition_add(terms, columns.values.column(j).size());
        mapped = partition_add(mapped, columns.map(j).size());
    }
    out.basis_offsets = rank_allocate<std::uint64_t>(plan, factor.rank + 1);
    out.source_offsets = rank_allocate<std::uint64_t>(plan, factor.rank + 1);
    out.basis_rows = rank_allocate<std::uint32_t>(plan, terms);
    out.basis_values = rank_allocate<RankFraction<RankCell>>(plan, terms);
    out.source_values = rank_allocate<RankSourceMap>(plan, mapped);

    std::uint64_t position = 0, map_position = 0;
    for (std::uint32_t k = 0; k < factor.rank; ++k) {
        const auto j = selected[k];
        budget.charge();
        if (j < rank_required)
            ++out.info.independent_required;
        bool has_c = false, has_rp = false;
        for (const auto& entry : columns.values.column(j)) {
            budget.charge();
            auto denominator = rank_product(budget, RankWide{columns.row_scale.values()[entry.row]},
                                            RankWide{columns.column_scale.values()[j]});
            const auto value =
                normalized_rank_fraction(budget, RankWide{entry.value}, std::move(denominator));
            out.basis_rows.values()[static_cast<std::size_t>(position)] = entry.row;
            out.basis_values.values()[static_cast<std::size_t>(position++)] = {
                RankCell{value.numerator}, RankCell{value.denominator}};
            if (out.rows.values()[entry.row].field == PgCoefficientField::C)
                has_c = true;
            else
                has_rp = true;
        }
        if (has_c && has_rp)
            partition_fail(core::FailureCode::RankUnresolved);
        if (has_c)
            ++out.info.rank_c;
        out.basis_offsets.values()[k + 1] = position;
        for (const auto& entry : columns.map(j)) {
            budget.charge();
            auto value = entry.coefficient != nullptr
                             ? rank_fraction<RankCell>(*entry.coefficient)
                             : RankFraction<RankCell>{RankCell{1}, RankCell{1}};
            if (entry.negative)
                value.numerator = -value.numerator;
            out.source_values.values()[static_cast<std::size_t>(map_position++)] = {
                entry.source, std::move(value)};
        }
        out.source_offsets.values()[k + 1] = map_position;
    }
    out.info.active_rows = shape.rows;
    out.info.anchor_rows = shape.trace_rows;
    out.info.candidates = shape.candidates;
    out.info.rank = factor.rank;
    out.info.normalized_terms = columns.values.offsets.values()[shape.candidates];
    out.info.proof_rows = partition_mul(shape.rows, shape.candidates);
    out.info.verified_anchor_rows = partition_mul(shape.originals, shape.trace_rows);
}

template <class I>
[[nodiscard]] ExactRational rank_exact(const RankFraction<I>& fraction, const ExactLimits& limits,
                                       core::ResourcePlan& plan) {
    if (fraction.numerator == 0)
        return {};
    const auto width = std::max(rank_bits(fraction.numerator), rank_bits(fraction.denominator));
    if (width > limits.maximum_value_bits)
        partition_fail(core::FailureCode::ResourceExceeded);
    if (sizeof(RationalStorage) > plan.remaining_bytes())
        partition_fail(core::FailureCode::ResourceExceeded);
    const auto admission = partition_take(exact_rational_admission(limits, width, 0));
    try {
        Rational raw;
        raw.backend().num() = fraction.numerator.backend();
        raw.backend().denom() = fraction.denominator.backend();
        return partition_take(ExactFactory::rational(std::move(raw), limits, plan, admission));
    } catch (const std::overflow_error&) {
        partition_fail(core::FailureCode::ResourceExceeded);
    } catch (const std::range_error&) {
        partition_fail(core::FailureCode::ResourceExceeded);
    }
}

[[nodiscard]] const PgRankState& rank_state(const core::PlannedBuffer<PgRankState>& buffer) {
    if (buffer.size() != 1)
        partition_fail(core::FailureCode::InvalidInput);
    return buffer.values().front();
}

[[nodiscard]] std::size_t rank_index(core::Index index, std::uint64_t count) {
    if (index < 0 || static_cast<std::uint64_t>(index) >= count)
        partition_fail(core::FailureCode::InvalidInput);
    return static_cast<std::size_t>(index);
}

[[nodiscard]] ExactRational rank_read_coordinate(const PgRankState& state, PgRankBasisColumn basis,
                                                 std::size_t candidate, const ExactLimits& limits,
                                                 core::ResourcePlan& plan) {
    const auto k = rank_index(basis.value, state.info.rank);
    return rank_exact(state.coordinates.values()[candidate * state.info.rank + k], limits, plan);
}

} // namespace

PgRankPreparation::PgRankPreparation(core::PlannedBuffer<PgRankState> state) noexcept
    : state_(std::move(state)) {}
PgRankPreparation::PgRankPreparation(PgRankPreparation&&) noexcept = default;
PgRankPreparation& PgRankPreparation::operator=(PgRankPreparation&&) noexcept = default;
PgRankPreparation::~PgRankPreparation() = default;

core::Result<PgRankPreparation>
PgRankPreparation::prepare(const PgCoefficientDictionary& dictionary,
                           const PgIntegrationPartition& partition,
                           const PgRankPreparationLimits& limits, core::Workspace& workspace,
                           core::ResourcePlan& plan) noexcept {
    return partition_entry([&]() -> core::Result<PgRankPreparation> {
        try {
            const auto lease = partition_take(workspace.acquire());
            static_cast<void>(lease);
            if (limits.maximum_digits == 0 || limits.maximum_digits > 128)
                partition_fail(core::FailureCode::ResourceExceeded);

            auto budget = RankBudget{limits};
            const auto shape = count_rank_input(dictionary, partition, budget);
            const auto domain_count = count_rank_domain(partition, budget);
            RankSizer forecast;
            count_rank_columns(forecast, shape);
            count_rank_factor(forecast, shape);
            count_rank_lift(forecast, shape);
            count_rank_state(forecast, shape);
            count_rank_witness(forecast, domain_count);
            if (forecast.bytes > limits.maximum_payload_bytes ||
                forecast.bytes > plan.remaining_bytes())
                partition_fail(core::FailureCode::ResourceExceeded);

            const auto starting = plan.allocated_bytes();
            auto state = rank_state_storage(plan, shape);
            auto& out = state.values().front();
            fill_rank_witness(out.witness, dictionary, partition, domain_count, budget, plan);
            const auto columns = assemble_rank_columns(dictionary, shape, budget, plan);
            auto factor = rank_factor_storage(plan, shape);
            auto lift = rank_lift_storage(plan, shape);
            rank_order_rows(shape, columns.values, factor, budget);
            for (const auto prime : rank_primes) {
                try {
                    rank_discover(shape, columns, factor, prime, budget);
                    rank_proof_and_coordinates(shape, columns, factor, lift, out, prime, budget);
                } catch (const PartitionAbort& error) {
                    if (error.failure.code == core::FailureCode::RankUnresolved)
                        continue;
                    throw;
                }
                rank_own_basis(shape, columns, factor, out, budget, plan);
                budget.charge();
                out.info.prime = prime;
                out.info.steps = budget.steps;
                out.info.preflight_bytes = forecast.bytes;
                out.info.allocated_payload_bytes = plan.allocated_bytes() - starting;
                if (out.info.allocated_payload_bytes > forecast.bytes)
                    partition_fail(core::FailureCode::ResourceExceeded);
                return PgRankPreparation{std::move(state)};
            }
            partition_fail(core::FailureCode::RankUnresolved);
        } catch (const std::overflow_error&) {
            partition_fail(core::FailureCode::ResourceExceeded);
        } catch (const std::range_error&) {
            partition_fail(core::FailureCode::ResourceExceeded);
        }
    });
}

core::Result<std::span<const PgRankRow>> PgRankPreparation::rows() const& noexcept {
    return partition_entry([&]() -> core::Result<std::span<const PgRankRow>> {
        return rank_state(state_).rows.values();
    });
}
core::Result<std::span<const std::uint32_t>> PgRankPreparation::selected() const& noexcept {
    return partition_entry([&]() -> core::Result<std::span<const std::uint32_t>> {
        const auto& state = rank_state(state_);
        return state.selected.values().first(state.info.rank);
    });
}
core::Result<PgRankPreparationInfo> PgRankPreparation::preparation() const noexcept {
    return partition_entry(
        [&]() -> core::Result<PgRankPreparationInfo> { return rank_state(state_).info; });
}
core::Result<PgRankWitnessInfo> PgRankPreparation::witness() const noexcept {
    return partition_entry([&]() -> core::Result<PgRankWitnessInfo> {
        const auto& source = rank_state(state_).witness;
        return PgRankWitnessInfo{source.n,
                                 source.convention,
                                 source.recipe_version,
                                 source.reduction_version,
                                 source.height_count,
                                 source.domain.size(),
                                 source.activity_count,
                                 source.graph_identity};
    });
}

core::Result<PgRankChartParameters>
PgRankPreparation::source_chart(const ExactLimits& limits,
                                core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<PgRankChartParameters> {
        const auto& source = rank_state(state_).witness;
        return PgRankChartParameters{rank_exact(source.x0, limits, plan),
                                     rank_exact(source.width, limits, plan),
                                     rank_exact(source.skew, limits, plan)};
    });
}

core::Result<PgRankHeightTerm>
PgRankPreparation::source_height(core::Index index, const ExactLimits& limits,
                                 core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<PgRankHeightTerm> {
        const auto& source = rank_state(state_).witness;
        const auto& value = source.height[rank_index(index, source.height_count)];
        return PgRankHeightTerm{value.powers, rank_exact(value.coefficient, limits, plan)};
    });
}

core::Result<PgRankDomainPoint>
PgRankPreparation::source_domain(core::Index index, const ExactLimits& limits,
                                 core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<PgRankDomainPoint> {
        const auto& source = rank_state(state_).witness;
        const auto& value = source.domain.values()[rank_index(index, source.domain.size())];
        return PgRankDomainPoint{rank_exact(value.s, limits, plan),
                                 rank_exact(value.t, limits, plan)};
    });
}

core::Result<PgRankActivityCell>
PgRankPreparation::source_activity(core::Index index, const ExactLimits& limits,
                                   core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<PgRankActivityCell> {
        const auto& source = rank_state(state_).witness;
        const auto& value = source.activity[rank_index(index, source.activity_count)];
        return PgRankActivityCell{value.cell, rank_exact(value.area, limits, plan), value.active};
    });
}

core::Result<ExactRational> PgRankPreparation::basis(PgRankBasisColumn column, core::Index row,
                                                     const ExactLimits& limits,
                                                     core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<ExactRational> {
        const auto& state = rank_state(state_);
        const auto k = rank_index(column.value, state.info.rank);
        const auto i = rank_index(row, state.info.active_rows);
        const auto begin = static_cast<std::size_t>(state.basis_offsets.values()[k]);
        const auto end = static_cast<std::size_t>(state.basis_offsets.values()[k + 1]);
        const auto rows = state.basis_rows.values().subspan(begin, end - begin);
        const auto found = std::lower_bound(rows.begin(), rows.end(), i);
        if (found == rows.end() || *found != i)
            return ExactRational{};
        return rank_exact(
            state.basis_values.values()[begin + static_cast<std::size_t>(found - rows.begin())],
            limits, plan);
    });
}
core::Result<ExactRational> PgRankPreparation::source_map(PgRankBasisColumn column,
                                                          PgRankOriginalColumn source,
                                                          const ExactLimits& limits,
                                                          core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<ExactRational> {
        const auto& state = rank_state(state_);
        const auto k = rank_index(column.value, state.info.rank);
        const auto i = rank_index(source.value, state.columns);
        const auto begin = static_cast<std::size_t>(state.source_offsets.values()[k]);
        const auto end = static_cast<std::size_t>(state.source_offsets.values()[k + 1]);
        const auto values = state.source_values.values().subspan(begin, end - begin);
        ExactRational result;
        for (const auto& entry : values)
            if (entry.source == i) {
                auto value = rank_exact(entry.coefficient, limits, plan);
                result = partition_take(result.add(value, limits, plan));
            }
        return result;
    });
}
core::Result<ExactRational> PgRankPreparation::original(PgRankBasisColumn basis,
                                                        PgRankOriginalColumn original,
                                                        const ExactLimits& limits,
                                                        core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<ExactRational> {
        const auto& state = rank_state(state_);
        const auto j = rank_index(original.value, state.columns);
        return rank_read_coordinate(state, basis, rank_required + j, limits, plan);
    });
}
core::Result<ExactRational> PgRankPreparation::required(PgRankBasisColumn basis,
                                                        PgRankRequiredColumn required,
                                                        const ExactLimits& limits,
                                                        core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<ExactRational> {
        const auto& state = rank_state(state_);
        return rank_read_coordinate(state, basis, rank_index(required.value, rank_required), limits,
                                    plan);
    });
}
core::Result<ExactRational> PgRankPreparation::anchor(PgRankBasisColumn basis,
                                                      PgRankAnchorColumn anchor,
                                                      const ExactLimits& limits,
                                                      core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<ExactRational> {
        const auto& state = rank_state(state_);
        const auto j = rank_index(anchor.value, state.originals);
        return rank_read_coordinate(state, basis, rank_required + state.columns + j, limits, plan);
    });
}

} // namespace uvlas::structure
