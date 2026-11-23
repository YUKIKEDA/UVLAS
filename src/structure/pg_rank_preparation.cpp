#include "pg_rank_lift.hpp"
#include "pg_rank_witness.hpp"
#include "pg_reference_polynomials.hpp"

#include <limits>
#include <tuple>
#include <utility>

namespace uvlas::structure::detail {
struct PgRankState {
    PgRankPreparationInfo info;
    RankWitness witness;
    std::array<PgCoefficientSource, rank_max_columns> provenance{};
    std::uint32_t columns = 0;
    std::uint32_t originals = 0;
    core::PlannedBuffer<PgRankRow> rows;
    core::PlannedBuffer<std::uint32_t> selected;
    core::PlannedBuffer<std::uint64_t> basis_offsets, recipe_offsets;
    core::PlannedBuffer<std::uint32_t> basis_rows;
    core::PlannedBuffer<RankFraction<RankCell>> basis_values;
    core::PlannedBuffer<RankSourceMap> recipe_values;
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
    size.add<std::uint64_t>(std::uint64_t{shape.rank_ceiling} + 1);
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

[[nodiscard]] RankFraction<RankCell> rank_recipe_value(const RankInputMap& entry,
                                                       RankBudget& budget) {
    budget.charge();
    if (entry.coefficient == nullptr)
        return {RankCell{entry.negative ? -1 : 1}, RankCell{1}};
    const auto& native = ExactFactory::value(*entry.coefficient).backend();
    budget.charge(std::uint64_t{native.num().size()} + native.denom().size());
    auto value = rank_fraction<RankCell>(*entry.coefficient);
    if (entry.negative)
        value.numerator = -value.numerator;
    return value;
}

void own_rank_recipe(PgRankState& out, RankColumns& columns,
                     const PgCoefficientDictionary& dictionary, const RankShape& shape,
                     RankBudget& budget, core::ResourcePlan& plan) {
    const auto sources = partition_take(dictionary.columns());
    budget.charge(sources.size());
    std::copy(sources.begin(), sources.end(), out.provenance.begin());
    out.recipe_values = rank_allocate<RankSourceMap>(plan, shape.map_terms);
    const auto borrowed = columns.map_entries.values();
    auto owned = out.recipe_values.values();
    for (std::size_t i = 0; i < borrowed.size(); ++i) {
        const auto& entry = borrowed[i];
        owned[i] = {entry.source, rank_recipe_value(entry, budget)};
    }
    out.recipe_offsets = std::move(columns.map_offsets);
}

[[nodiscard]] std::span<const RankSourceMap> rank_recipe(const PgRankState& out,
                                                         std::size_t candidate) {
    const auto offsets = out.recipe_offsets.values();
    const auto start = static_cast<std::size_t>(offsets[candidate]);
    const auto stop = static_cast<std::size_t>(offsets[candidate + 1]);
    return out.recipe_values.values().subspan(start, stop - start);
}

void rank_match_recipe(const PgRankState& previous, const PgCoefficientDictionary& dictionary,
                       const RankShape& shape, RankBudget& budget) {
    const auto sources = partition_take(dictionary.columns());
    budget.charge();
    if (shape.columns != previous.columns || shape.originals != previous.originals)
        partition_fail(core::FailureCode::StaleGeneration);
    for (std::size_t i = 0; i < sources.size(); ++i) {
        budget.charge();
        const auto& a = previous.provenance[i];
        const auto& b = sources[i];
        if (a.kind != b.kind || a.first != b.first || a.second != b.second ||
            a.component != b.component || a.axis != b.axis || a.original != b.original)
            partition_fail(core::FailureCode::StaleGeneration);
    }

    std::array<std::size_t, rank_max_candidates> counts{};
    rank_map_terms(dictionary, shape, [&](std::uint32_t candidate, const RankInputMap& entry) {
        budget.charge();
        if (candidate >= shape.candidates)
            partition_fail(core::FailureCode::RankUnresolved);
        const auto saved = rank_recipe(previous, candidate);
        const auto i = counts[candidate]++;
        if (i >= saved.size() || saved[i].source != entry.source ||
            !rank_equal_fraction(saved[i].coefficient, rank_recipe_value(entry, budget), budget))
            partition_fail(core::FailureCode::StaleGeneration);
    });
    for (std::size_t candidate = 0; candidate < shape.candidates; ++candidate) {
        budget.charge();
        if (counts[candidate] != rank_recipe(previous, candidate).size())
            partition_fail(core::FailureCode::StaleGeneration);
    }
}

void rank_own_basis(const RankShape& shape, const RankColumns& columns, const RankFactor& factor,
                    PgRankState& out, RankBudget& budget, core::ResourcePlan& plan) {
    const auto selected = factor.selected.values().first(factor.rank);
    std::uint64_t terms = 0;
    for (const auto j : selected) {
        budget.charge();
        terms = partition_add(terms, columns.values.column(j).size());
    }
    out.basis_offsets = rank_allocate<std::uint64_t>(plan, factor.rank + 1);
    out.basis_rows = rank_allocate<std::uint32_t>(plan, terms);
    out.basis_values = rank_allocate<RankFraction<RankCell>>(plan, terms);

    std::uint64_t position = 0;
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

[[nodiscard]] core::PlannedBuffer<PgRankState>
prepare_rank_state(const PgCoefficientDictionary& dictionary,
                   const PgIntegrationPartition& partition, RankBudget& budget,
                   core::ResourcePlan& plan, const PgRankState* previous) {
    const auto shape = count_rank_input(dictionary, partition, budget);
    const auto domain_count = count_rank_domain(partition, budget);
    if (previous != nullptr) {
        rank_match_witness(previous->witness, partition_take(dictionary.graph()), partition,
                           domain_count, budget);
        rank_match_recipe(*previous, dictionary, shape, budget);
    }
    RankSizer forecast;
    count_rank_columns(forecast, shape);
    count_rank_factor(forecast, shape);
    count_rank_lift(forecast, shape);
    count_rank_state(forecast, shape);
    count_rank_witness(forecast, domain_count);
    if (forecast.bytes > budget.limits.maximum_payload_bytes ||
        forecast.bytes > plan.remaining_bytes())
        partition_fail(core::FailureCode::ResourceExceeded);

    const auto starting = plan.allocated_bytes();
    auto state = rank_state_storage(plan, shape);
    auto& out = state.values().front();
    fill_rank_witness(out.witness, dictionary, partition, domain_count, budget, plan);
    auto columns = assemble_rank_columns(dictionary, shape, budget, plan);
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
        own_rank_recipe(out, columns, dictionary, shape, budget, plan);
        budget.charge();
        out.info.prime = prime;
        out.info.steps = budget.steps;
        out.info.preflight_bytes = forecast.bytes;
        out.info.allocated_payload_bytes = plan.allocated_bytes() - starting;
        if (out.info.allocated_payload_bytes > forecast.bytes)
            partition_fail(core::FailureCode::ResourceExceeded);
        return state;
    }
    partition_fail(core::FailureCode::RankUnresolved);
}

template <class T>
void rank_compare(std::span<const T> a, std::span<const T> b, RankBudget& budget) {
    budget.charge(a.size());
    if (!std::ranges::equal(a, b))
        partition_fail(core::FailureCode::RankUnresolved);
}

void rank_compare_result(const PgRankState& a, const PgRankState& b, RankBudget& budget) {
    if (a.info.rank != b.info.rank || a.info.rank_c != b.info.rank_c ||
        a.info.independent_required != b.info.independent_required ||
        a.info.active_rows != b.info.active_rows || a.info.candidates != b.info.candidates)
        partition_fail(core::FailureCode::RankUnresolved);
    const auto r = static_cast<std::size_t>(a.info.rank);
    rank_compare(a.selected.values().first(r), b.selected.values().first(r), budget);
    rank_compare(a.basis_offsets.values(), b.basis_offsets.values(), budget);
    rank_compare(a.basis_rows.values(), b.basis_rows.values(), budget);
    const auto rows = a.rows.values();
    for (std::size_t i = 0; i < rows.size(); ++i) {
        budget.charge();
        const auto& x = rows[i];
        const auto& y = b.rows.values()[i];
        if (x.field != y.field || x.cell != y.cell || x.powers.first != y.powers.first ||
            x.powers.second != y.powers.second)
            partition_fail(core::FailureCode::RankUnresolved);
    }
    for (std::size_t i = 0; i < a.basis_values.size(); ++i)
        if (!rank_equal_fraction(a.basis_values.values()[i], b.basis_values.values()[i], budget))
            partition_fail(core::FailureCode::RankUnresolved);
    const auto coordinates = r * static_cast<std::size_t>(a.info.candidates);
    for (std::size_t i = 0; i < coordinates; ++i)
        if (!rank_equal_fraction(a.coordinates.values()[i], b.coordinates.values()[i], budget))
            partition_fail(core::FailureCode::RankUnresolved);
    // Every source recipe was matched against the current input before replaying.
    // The identical selected ordinals therefore also certify the same T.
}

template <class T, class F>
core::Result<T> rank_entry(F&& fn, const PgRankPreparationLimits& limits,
                           core::Workspace& workspace) noexcept {
    return partition_entry([&]() -> core::Result<T> {
        try {
            const auto lease = partition_take(workspace.acquire());
            static_cast<void>(lease);
            auto budget = RankBudget{limits};
            budget.charge(0);
            if (limits.maximum_digits == 0 || limits.maximum_digits > 128)
                partition_fail(core::FailureCode::ResourceExceeded);
            return fn(budget);
        } catch (const std::overflow_error&) {
            partition_fail(core::FailureCode::ResourceExceeded);
        } catch (const std::range_error&) {
            partition_fail(core::FailureCode::ResourceExceeded);
        }
    });
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
    return rank_entry<PgRankPreparation>(
        [&](RankBudget& budget) -> core::Result<PgRankPreparation> {
            return PgRankPreparation{
                prepare_rank_state(dictionary, partition, budget, plan, nullptr)};
        },
        limits, workspace);
}

core::Result<PgRankValidationInfo>
PgRankPreparation::validate_for(const PgCoefficientDictionary& dictionary,
                                const PgIntegrationPartition& partition,
                                const PgRankPreparationLimits& limits, core::Workspace& workspace,
                                core::ResourcePlan& plan) const noexcept {
    return rank_entry<PgRankValidationInfo>(
        [&](RankBudget& budget) -> core::Result<PgRankValidationInfo> {
            const auto& previous = rank_state(state_);
            const auto replay = prepare_rank_state(dictionary, partition, budget, plan, &previous);
            const auto& current = rank_state(replay);
            rank_compare_result(previous, current, budget);
            budget.charge(0);
            return PgRankValidationInfo{budget.steps, current.info.preflight_bytes,
                                        current.info.allocated_payload_bytes};
        },
        limits, workspace);
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

core::Result<std::span<const PgCoefficientSource>> PgRankPreparation::provenance() const& noexcept {
    return partition_entry([&]() -> core::Result<std::span<const PgCoefficientSource>> {
        const auto& state = rank_state(state_);
        return std::span{state.provenance}.first(state.columns);
    });
}

core::Result<std::uint64_t>
PgRankPreparation::recipe_terms(PgRankRequiredColumn target) const noexcept {
    return partition_entry([&]() -> core::Result<std::uint64_t> {
        return rank_recipe(rank_state(state_), rank_index(target.value, rank_required)).size();
    });
}

core::Result<std::uint64_t>
PgRankPreparation::recipe_terms(PgRankAnchorColumn target) const noexcept {
    return partition_entry([&]() -> core::Result<std::uint64_t> {
        const auto& state = rank_state(state_);
        return rank_recipe(state, rank_required + state.columns +
                                      rank_index(target.value, state.originals))
            .size();
    });
}

core::Result<PgRankRecipeTerm> PgRankPreparation::recipe(PgRankRequiredColumn target,
                                                         core::Index term,
                                                         const ExactLimits& limits,
                                                         core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<PgRankRecipeTerm> {
        const auto values =
            rank_recipe(rank_state(state_), rank_index(target.value, rank_required));
        const auto& entry = values[rank_index(term, values.size())];
        return PgRankRecipeTerm{{entry.source}, rank_exact(entry.coefficient, limits, plan)};
    });
}

core::Result<PgRankRecipeTerm> PgRankPreparation::recipe(PgRankAnchorColumn target,
                                                         core::Index term,
                                                         const ExactLimits& limits,
                                                         core::ResourcePlan& plan) const noexcept {
    return partition_entry([&]() -> core::Result<PgRankRecipeTerm> {
        const auto& state = rank_state(state_);
        const auto values = rank_recipe(state, rank_required + state.columns +
                                                   rank_index(target.value, state.originals));
        const auto& entry = values[rank_index(term, values.size())];
        return PgRankRecipeTerm{{entry.source}, rank_exact(entry.coefficient, limits, plan)};
    });
}

core::Result<std::optional<PgRankTranslationAxis>>
PgRankPreparation::translation_axis(PgRankRequiredColumn target) const noexcept {
    return partition_entry([&]() -> core::Result<std::optional<PgRankTranslationAxis>> {
        static_cast<void>(rank_state(state_));
        const auto index = rank_index(target.value, rank_required);
        constexpr std::array axes{PgRankTranslationAxis::X, PgRankTranslationAxis::Y,
                                  PgRankTranslationAxis::Z};
        if (index < rank_required - axes.size())
            return std::nullopt;
        return axes[index - (rank_required - axes.size())];
    });
}

std::array<std::array<ExactRational, 3>, 3>
PgRankPreparation::basis_cell_jet(ReferenceArithmetic& math, PgRankBasisColumn column,
                                  PgCellIndex cell, const ExactRational& u, const ExactRational& v,
                                  const ExactRational& axis_scale) const {
    const auto& state = rank_state(state_);
    math.visit(1);
    const auto k = rank_index(column.value, state.info.rank);
    const auto begin = static_cast<std::size_t>(state.basis_offsets.values()[k]);
    const auto end = static_cast<std::size_t>(state.basis_offsets.values()[k + 1]);
    const auto indices = state.basis_rows.values().subspan(begin, end - begin);
    const auto key = [&](std::uint32_t ordinal) {
        math.visit(1);
        const auto& row = state.rows.values()[ordinal];
        return std::tuple{row.field, row.cell.first, row.cell.second};
    };
    constexpr std::array fields{PgCoefficientField::R, PgCoefficientField::P,
                                PgCoefficientField::C};
    std::array<std::array<ExactRational, 3>, 3> jet{};
    for (std::size_t field = 0; field < fields.size(); ++field) {
        const auto target = std::tuple{fields[field], cell.first, cell.second};
        const auto first = std::partition_point(indices.begin(), indices.end(),
                                                [&](auto i) { return key(i) < target; });
        for (auto term = first; term != indices.end() && key(*term) == target; ++term) {
            const auto& row = state.rows.values()[*term];
            const auto offset = begin + static_cast<std::size_t>(term - indices.begin());
            math.charge();
            const auto coefficient =
                rank_exact(state.basis_values.values()[offset], math.leaf_limits, math.plan);
            for (std::size_t i = 0; i < jet[field].size(); ++i) {
                const auto order = reference_orders[i];
                if (order.first > row.powers.first || order.second > row.powers.second)
                    continue;
                auto value =
                    times_monomial(math, math.copy(coefficient), row.powers.first, order.first, u);
                value = times_monomial(math, std::move(value), row.powers.second, order.second, v);
                jet[field][i] = math.add(jet[field][i], value);
            }
        }
        for (std::size_t i = 1; i < jet[field].size(); ++i)
            jet[field][i] = math.mul(jet[field][i], axis_scale);
    }
    return jet;
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
        const auto values = rank_recipe(state, state.selected.values()[k]);
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
