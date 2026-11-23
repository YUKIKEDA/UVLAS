#pragma once

#include "pg_rank_storage.hpp"

namespace uvlas::structure::detail {

struct RankWitnessTerm {
    ExactPowers powers{};
    RankFraction<RankWide> coefficient{};
};

struct RankWitnessPoint {
    RankFraction<RankWide> s{}, t{};
};

struct RankWitnessCell {
    PgCellIndex cell{};
    RankFraction<RankWide> area{};
    bool active = false;
};

inline constexpr std::uint8_t rank_recipe_version = 1;
inline constexpr std::uint8_t rank_reduction_version = 1;

struct RankWitness {
    RankFraction<RankWide> x0{}, width{}, skew{};
    std::array<RankWitnessTerm, 15> height{};
    std::size_t height_count = 0;
    core::Index n = 0;
    PgGraphConvention convention = PgGraphConvention::Version1;
    PgGraphContentIdentity graph_identity{{}};
    std::uint8_t recipe_version = rank_recipe_version;
    std::uint8_t reduction_version = rank_reduction_version;
    std::array<RankWitnessCell, rank_max_cells> activity{};
    std::size_t activity_count = 0;
    core::PlannedBuffer<RankWitnessPoint> domain;
};

[[nodiscard]] inline std::uint64_t count_rank_domain(const PgIntegrationPartition& partition,
                                                     RankBudget& budget) {
    const auto sources = partition_take(partition.sources());
    const auto points = partition_take(partition.source_points());
    std::uint64_t count = 0;
    for (const auto& source : sources) {
        budget.charge();
        if (std::holds_alternative<PgPartitionDomainSource>(source.descriptor))
            ++count;
    }
    if (count < 3 || count > points.size())
        partition_fail(core::FailureCode::RankUnresolved);

    // The partition owns the original polygon in the initial point slice.
    // Check every directed edge before using that convention to copy it.
    std::uint64_t ordinal = 0;
    for (const auto& source : sources) {
        budget.charge();
        const auto* edge = std::get_if<PgPartitionDomainSource>(&source.descriptor);
        if (edge == nullptr)
            continue;
        const auto start = partition_index(ordinal);
        const auto end = partition_index((ordinal + 1) % count);
        if (edge->edge != start || source.endpoints[0] != start || source.endpoints[1] != end)
            partition_fail(core::FailureCode::RankUnresolved);
        ++ordinal;
    }
    return count;
}

inline void count_rank_witness(RankSizer& size, std::uint64_t domain_count) {
    size.add<RankWitnessPoint>(domain_count);
}

[[nodiscard]] inline RankFraction<RankWide> rank_witness_value(const ExactRational& value,
                                                               RankBudget& budget) {
    constexpr std::uint64_t max_bits = 8192;
    constexpr auto limb_bits = IntegerBackend::limb_bits;
    const auto numerator_bits = value.numerator_bits();
    const auto denominator_bits = value.denominator_bits();
    if (numerator_bits > max_bits || denominator_bits > max_bits)
        partition_fail(core::FailureCode::ResourceExceeded);
    budget.charge((numerator_bits + limb_bits - 1) / limb_bits +
                  (denominator_bits + limb_bits - 1) / limb_bits + 1);
    return rank_fraction<RankWide>(value);
}

inline void fill_rank_witness(RankWitness& witness, const PgCoefficientDictionary& dictionary,
                              const PgIntegrationPartition& partition, std::uint64_t domain_count,
                              RankBudget& budget, core::ResourcePlan& plan) {
    const auto graph = partition_take(dictionary.graph());
    const auto height = graph.height.terms();
    const auto cells = partition_take(partition.activity_map());
    if (height.size() > witness.height.size() || cells.size() > witness.activity.size())
        partition_fail(core::FailureCode::RankUnresolved);
    witness.x0 = rank_witness_value(graph.x0, budget);
    witness.width = rank_witness_value(graph.width, budget);
    witness.skew = rank_witness_value(graph.skew, budget);
    witness.n = graph.n;
    witness.convention = graph.convention;
    witness.graph_identity = graph.identity;
    witness.height_count = height.size();
    for (std::size_t i = 0; i < height.size(); ++i) {
        witness.height[i].powers = height[i].powers;
        witness.height[i].coefficient = rank_witness_value(height[i].coefficient, budget);
    }
    witness.activity_count = cells.size();
    for (std::size_t i = 0; i < cells.size(); ++i) {
        witness.activity[i].cell = cells[i].cell;
        witness.activity[i].area = rank_witness_value(cells[i].area, budget);
        witness.activity[i].active = cells[i].active;
    }

    const auto points =
        partition_take(partition.source_points()).first(static_cast<std::size_t>(domain_count));
    witness.domain = rank_allocate<RankWitnessPoint>(plan, domain_count);
    for (std::size_t i = 0; i < points.size(); ++i) {
        auto& point = witness.domain.values()[i];
        point.s = rank_witness_value(points[i].s, budget);
        point.t = rank_witness_value(points[i].t, budget);
    }
}

inline void rank_match_witness(const RankWitness& witness, const PgGraphView& graph,
                               const PgIntegrationPartition& partition, std::uint64_t domain_count,
                               RankBudget& budget) {
    const auto height = graph.height.terms();
    const auto cells = partition_take(partition.activity_map());
    const auto match = [&](const RankFraction<RankWide>& old, const ExactRational& current) {
        return rank_equal_fraction(old, rank_witness_value(current, budget), budget);
    };
    budget.charge();
    if (witness.n != graph.n || witness.convention != graph.convention ||
        witness.recipe_version != rank_recipe_version ||
        witness.reduction_version != rank_reduction_version ||
        height.size() != witness.height_count || cells.size() != witness.activity_count ||
        domain_count != witness.domain.size() || !match(witness.x0, graph.x0) ||
        !match(witness.width, graph.width) || !match(witness.skew, graph.skew))
        partition_fail(core::FailureCode::StaleGeneration);
    for (std::size_t i = 0; i < height.size(); ++i) {
        budget.charge();
        const auto& saved = witness.height[i];
        if (saved.powers.first != height[i].powers.first ||
            saved.powers.second != height[i].powers.second ||
            !match(saved.coefficient, height[i].coefficient))
            partition_fail(core::FailureCode::StaleGeneration);
    }
    for (std::size_t i = 0; i < cells.size(); ++i) {
        budget.charge();
        const auto& saved = witness.activity[i];
        if (saved.cell != cells[i].cell || saved.active != cells[i].active ||
            !match(saved.area, cells[i].area))
            partition_fail(core::FailureCode::StaleGeneration);
    }
    const auto points = partition_take(partition.source_points());
    for (std::size_t i = 0; i < witness.domain.size(); ++i) {
        budget.charge();
        const auto& saved = witness.domain.values()[i];
        if (!match(saved.s, points[i].s) || !match(saved.t, points[i].t))
            partition_fail(core::FailureCode::StaleGeneration);
    }
}

} // namespace uvlas::structure::detail
