#pragma once

#include "pg_rank_storage.hpp"

#include <algorithm>
#include <array>
#include <span>

namespace uvlas::structure::detail {

inline constexpr std::size_t cut_degree = 2 * std::size_t{rank_degree[1]};
inline constexpr std::size_t cut_slots = (cut_degree + 1) * (cut_degree + 2) / 2;
inline constexpr std::size_t cut_edge_degree = cut_degree + 1;
inline constexpr std::size_t cut_denominator_power = cut_degree + 2;
static_assert(cut_slots == 946 && cut_denominator_power == 44);

struct GramCutLayout {
    std::array<core::Index, rank_max_cells> slots{};
    std::uint64_t count = 0;
    core::Index axis = 0;

    explicit GramCutLayout(core::Index axis_cells) noexcept : axis(axis_cells) {
        slots.fill(-1);
    }
    [[nodiscard]] std::size_t index(PgCellIndex cell) const noexcept {
        return static_cast<std::size_t>(cell.first * axis + cell.second);
    }
    [[nodiscard]] core::Index slot(PgCellIndex cell) const noexcept {
        return slots[index(cell)];
    }
};

[[nodiscard]] inline std::size_t cut_index(std::size_t u, std::size_t v) noexcept {
    const auto degree = u + v;
    return degree * (degree + 1) / 2 + u;
}

template <class Budget>
[[nodiscard]] RankWide gram_divide(Budget& work, const RankWide& a, const RankWide& b) {
    work.charge(partition_mul(a.backend().size(), b.backend().size()));
    return RankWide{a / b};
}

template <class Budget>
[[nodiscard]] RankFraction<RankWide> cut_fraction(Budget& work, const ExactRational& exact) {
    const auto& source = ExactFactory::value(exact).backend();
    work.charge(std::uint64_t{source.num().size()} + source.denom().size());
    return rank_fraction<RankWide>(exact);
}

template <class Budget, class Visit>
void visit_cut_boundary(Budget& work, const PgIntegrationPartition& partition,
                        const GramCutLayout& layout, const Visit& visit) {
    const auto edges = partition_take(partition.edges());
    const auto pieces = partition_take(partition.pieces());
    const auto vertices = partition_take(partition.vertices());
    for (const auto& edge : edges) {
        work.charge();
        const auto left = edge.left.value_or(-1), right = edge.right.value_or(-1);
        const auto left_slot =
            left < 0 ? -1 : layout.slot(pieces[static_cast<std::size_t>(left)].cell);
        const auto right_slot =
            right < 0 ? -1 : layout.slot(pieces[static_cast<std::size_t>(right)].cell);
        if (left_slot == right_slot)
            continue;
        const auto& a = vertices[static_cast<std::size_t>(edge.vertices.front())];
        const auto& b = vertices[static_cast<std::size_t>(edge.vertices.back())];
        if (left_slot >= 0)
            visit(left_slot, pieces[static_cast<std::size_t>(left)].cell, a, b, true);
        if (right_slot >= 0)
            visit(right_slot, pieces[static_cast<std::size_t>(right)].cell, a, b, false);
    }
}

[[nodiscard]] inline std::span<RankWide> cut_cell(std::span<RankWide> moments, core::Index slot) {
    return moments.subspan(static_cast<std::size_t>(slot) * cut_slots, cut_slots);
}
[[nodiscard]] inline std::span<const RankWide> cut_cell(std::span<const RankWide> moments,
                                                        core::Index slot) {
    return moments.subspan(static_cast<std::size_t>(slot) * cut_slots, cut_slots);
}

template <class Budget>
[[nodiscard]] RankWide cut_coordinate(Budget& work, const ExactRational& input,
                                      const RankWide& scale, core::Index four_n,
                                      core::Index cell_index) {
    const auto term = cut_fraction(work, input);
    const auto scaled =
        rank_product(work, term.numerator, gram_divide(work, scale, term.denominator));
    return rank_product(work, scaled, RankWide{four_n}) -
           rank_product(work, scale, RankWide{2 * cell_index + 1});
}

struct CutPolynomial {
    std::array<RankWide, cut_edge_degree + 1> terms{};

    template <class Budget>
    void extend(Budget& work, std::size_t degree, const RankWide& origin, const RankWide& delta) {
        terms[degree + 1] = rank_product(work, terms[degree], delta);
        for (auto k = degree; k > 0; --k) {
            work.charge();
            terms[k] =
                rank_product(work, terms[k], origin) + rank_product(work, terms[k - 1], delta);
        }
        terms[0] = rank_product(work, terms[0], origin);
    }
};

struct CutScale {
    std::array<RankWide, cut_denominator_power + 1> powers{};
    std::array<RankWide, cut_edge_degree + 1> integration{};
    std::array<RankWide, cut_degree + 1> green{};
    RankWide denominator;
};

template <class Budget>
void accumulate_cut_edge(Budget& work, const CutScale& scale, std::span<RankWide> moments,
                         const RankWide& u0, const RankWide& u_delta, const RankWide& v0,
                         const RankWide& v_delta) {
    if (v_delta == 0)
        return;
    CutPolynomial u{}, uv{};
    u.terms[0] = 1;
    for (std::size_t a = 0; a <= cut_degree; ++a) {
        u.extend(work, a, u0, u_delta);
        work.charge(a + 2);
        std::ranges::copy(std::span{u.terms}.first(a + 2), uv.terms.begin());
        for (std::size_t b = 0; b <= cut_degree - a; ++b) {
            const auto degree = a + b + 1;
            RankWide integral{};
            for (std::size_t k = 0; k <= degree; ++k) {
                work.charge();
                integral += rank_product(work, uv.terms[k], scale.integration[k]);
            }
            auto scaled = rank_product(work, integral, v_delta);
            scaled = rank_product(work, scaled, scale.green[a]);
            scaled = rank_product(work, scaled, scale.powers[cut_denominator_power - degree - 1]);
            work.charge();
            moments[cut_index(a, b)] += scaled;
            if (b < cut_degree - a)
                uv.extend(work, degree, v0, v_delta);
        }
    }
}

template <class Budget>
[[nodiscard]] RankWide generate_cut_moments(Budget& work, const PgIntegrationPartition& partition,
                                            const GramCutLayout& layout,
                                            std::span<RankWide> moments) {
    RankWide common{1};
    std::uint64_t side_count = 0;
    visit_cut_boundary(
        work, partition, layout,
        [&](core::Index, PgCellIndex, const PgPartitionPoint& a, const PgPartitionPoint& b, bool) {
            side_count = partition_add(side_count, 1);
            for (const auto* x : {&a.s, &a.t, &b.s, &b.t}) {
                const auto& source = ExactFactory::value(*x).backend().denom();
                work.charge(source.size());
                common = rank_lcm(work, common, rank_copy_integer<RankWide>(source));
            }
        });
    RankWide green_denominator{1};
    for (std::size_t k = 1; k <= cut_edge_degree; ++k)
        green_denominator = rank_lcm(work, green_denominator, RankWide{k});
    const auto integration_denominator =
        rank_lcm(work, green_denominator, RankWide{cut_denominator_power});
    // Endpoint local numerators are at most S in magnitude and increments at
    // most 2*S. The coefficient L1-norm is bounded by (3*S)^degree.
    // Include the edge derivative, both LCMs, and accumulation across all sides.
    const auto power_bits = partition_mul(cut_denominator_power, rank_bits(common) + 2);
    const auto extra_bits = rank_bits(green_denominator) + rank_bits(integration_denominator) +
                            partition_log(std::max(side_count, std::uint64_t{1}));
    if (partition_add(power_bits, extra_bits) > 8192)
        partition_fail(core::FailureCode::ResourceExceeded);

    CutScale scale{};
    scale.powers.front() = 1;
    for (std::size_t k = 1; k < scale.powers.size(); ++k)
        scale.powers[k] = rank_product(work, scale.powers[k - 1], common);
    for (std::size_t k = 0; k < scale.integration.size(); ++k)
        scale.integration[k] = gram_divide(work, integration_denominator, RankWide{k + 1});
    for (std::size_t k = 0; k < scale.green.size(); ++k)
        scale.green[k] = gram_divide(work, green_denominator, RankWide{k + 1});
    scale.denominator = rank_product(
        work, rank_product(work, green_denominator, integration_denominator), scale.powers.back());

    const auto four_n = 2 * layout.axis;
    visit_cut_boundary(
        work, partition, layout,
        [&](core::Index slot, PgCellIndex cell, const PgPartitionPoint& a,
            const PgPartitionPoint& b, bool ccw) {
            // Reverse endpoint order for the piece to the right of a canonical edge.
            const auto& first = ccw ? a : b;
            const auto& last = ccw ? b : a;
            const auto u0 = cut_coordinate(work, first.s, common, four_n, cell.first);
            const auto v0 = cut_coordinate(work, first.t, common, four_n, cell.second);
            const auto u_delta = cut_coordinate(work, last.s, common, four_n, cell.first) - u0;
            const auto v_delta = cut_coordinate(work, last.t, common, four_n, cell.second) - v0;
            accumulate_cut_edge(work, scale, cut_cell(moments, slot), u0, u_delta, v0, v_delta);
        });

    // The constant moment independently reconciles the edge roster with the
    // certified chart area for every affected background cell.
    const auto activity = partition_take(partition.activity_map());
    for (std::size_t c = 0; c < activity.size(); ++c) {
        work.charge();
        const auto slot = layout.slots[c];
        if (slot < 0)
            continue;
        const auto area = cut_fraction(work, activity[c].area);
        const auto local_num = rank_product(work, area.numerator, RankWide{four_n * four_n});
        const auto factor = gram_divide(work, scale.denominator, area.denominator);
        if (rank_product(work, area.denominator, factor) != scale.denominator ||
            rank_product(work, local_num, factor) != cut_cell(moments, slot).front())
            partition_fail(core::FailureCode::RankUnresolved);
    }
    return std::move(scale.denominator);
}

} // namespace uvlas::structure::detail
