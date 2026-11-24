#pragma once

#include "pg_flat_gram_cut_moments.hpp"
#include "pg_reference_domain.hpp"

namespace uvlas::structure::detail {

inline constexpr std::array<std::size_t, 3> thickness_gram_degree{2 * std::size_t{rank_degree[0]},
                                                                  2 * std::size_t{rank_degree[1]},
                                                                  2 * std::size_t{rank_degree[2]}};
inline constexpr std::size_t thickness_polygon_degree = thickness_gram_degree[1] + 3;
inline constexpr std::array<std::size_t, 4> thickness_gram_offset = [] {
    std::array<std::size_t, 4> offsets{};
    for (std::size_t j = 0; j < thickness_gram_degree.size(); ++j)
        offsets[j + 1] = offsets[j] + cut_index(0, thickness_gram_degree[j] + 1);
    return offsets;
}();
inline constexpr auto thickness_gram_slots = thickness_gram_offset.back();
inline constexpr auto thickness_piece_slots =
    (thickness_polygon_degree + 1) * (thickness_polygon_degree + 2) / 2;
static_assert(thickness_gram_slots == 1339 && thickness_piece_slots == 1081);

struct GramThicknessDomain {
    PgPartitionThicknessView thickness;
    std::span<const PgPartitionPiece> pieces;
    std::span<const PgPartitionPoint> vertices;
    std::span<const core::Index> indices;
};

[[nodiscard]] inline GramThicknessDomain
gram_thickness_domain(const PgIntegrationPartition& partition) {
    return {partition_take(partition.thickness()), partition_take(partition.pieces()),
            partition_take(partition.vertices()), partition_take(partition.piece_vertex_indices())};
}

[[nodiscard]] inline std::array<R, 3> local_half_thickness(ReferenceArithmetic& math,
                                                           PgPartitionThicknessView thickness,
                                                           const PgPartitionPiece& piece,
                                                           const R& four_n, const R& two_length,
                                                           const R& slope_scale) {
    math.visit(1);
    const auto plane = thickness_plane(math, thickness, piece.thickness_triangle);
    const auto s = math.div(index_value(math, 2 * piece.cell.first + 1), four_n);
    const auto t = math.div(index_value(math, 2 * piece.cell.second + 1), four_n);
    return {math.div(plane.evaluate(math, s, t), two_length),
            math.div(plane.gradient[0], slope_scale), math.div(plane.gradient[1], slope_scale)};
}

[[nodiscard]] inline std::span<const core::Index>
gram_piece_vertices(const GramThicknessDomain& domain, const PgPartitionPiece& piece) {
    return domain.indices.subspan(static_cast<std::size_t>(piece.first_vertex),
                                  static_cast<std::size_t>(piece.vertex_count));
}

struct GramThicknessScale {
    RankWide chart{1}, height{1};
    std::uint64_t maximum_coefficient_bits = 0;
};

template <class Budget>
[[nodiscard]] GramThicknessScale
gram_thickness_scale(Budget& work, const GramThicknessDomain& domain, const R& four_n,
                     const R& two_length, const R& slope_scale) {
    GramThicknessScale scale{};
    for (const auto& piece : domain.pieces) {
        work.charge();
        for (const auto& value : local_half_thickness(work.math, domain.thickness, piece, four_n,
                                                      two_length, slope_scale)) {
            const auto fraction = cut_fraction(work, value);
            scale.height = rank_lcm(work, scale.height, fraction.denominator);
            scale.maximum_coefficient_bits =
                std::max(scale.maximum_coefficient_bits, rank_bits(fraction.numerator));
        }
        for (const auto index : gram_piece_vertices(domain, piece)) {
            work.charge();
            const auto& vertex = domain.vertices[static_cast<std::size_t>(index)];
            for (const auto* coordinate : {&vertex.s, &vertex.t}) {
                const auto& denominator = ExactFactory::value(*coordinate).backend().denom();
                work.charge(denominator.size());
                scale.chart = rank_lcm(work, scale.chart, rank_copy_integer<RankWide>(denominator));
            }
        }
    }
    return scale;
}

struct GramThicknessPower {
    std::array<RankWide, 21> coefficients = [] {
        std::array<RankWide, 21> values{};
        values.front() = 1;
        return values;
    }();

    template <class Budget>
    void extend(Budget& work, std::size_t new_degree, const std::array<RankWide, 3>& affine) {
        for (auto d = new_degree + 1; d > 0; --d) {
            const auto degree = d - 1;
            for (std::size_t u = 0; u <= degree; ++u) {
                work.charge();
                const auto v = degree - u;
                auto& value = coefficients[cut_index(u, v)];
                auto next = rank_product(work, value, affine[0]);
                if (u > 0)
                    next += rank_product(work, coefficients[cut_index(u - 1, v)], affine[1]);
                if (v > 0)
                    next += rank_product(work, coefficients[cut_index(u, v - 1)], affine[2]);
                value = std::move(next);
            }
        }
    }
};

template <class Budget>
void combine_thickness_moments(Budget& work, std::span<RankWide> target,
                               std::span<const RankWide> piece,
                               const std::array<RankWide, 3>& affine,
                               const std::array<RankWide, 6>& height_powers) {
    GramThicknessPower h;
    for (std::size_t p = 1; p <= 5; ++p) {
        h.extend(work, p, affine);
        if (p % 2 == 0)
            continue;
        const auto mode = p / 2;
        const auto offset = thickness_gram_offset[mode];
        const auto maximum = thickness_gram_degree[mode];
        auto power = h.coefficients;
        work.charge(power.size());
        for (auto& term : std::span{power}.first(cut_index(0, p + 1)))
            term = rank_product(work, term, height_powers[5 - p]);
        for (std::size_t a = 0; a <= maximum; ++a)
            for (std::size_t b = 0; b <= maximum - a; ++b) {
                RankWide sum{};
                for (std::size_t u = 0; u <= p; ++u)
                    for (std::size_t v = 0; v <= p - u; ++v) {
                        work.charge();
                        const auto& coefficient = power[cut_index(u, v)];
                        if (coefficient != 0)
                            sum += rank_product(work, coefficient, piece[cut_index(a + u, b + v)]);
                    }
                target[offset + cut_index(a, b)] += sum;
            }
    }
}

template <class Budget>
[[nodiscard]] RankWide
generate_thickness_moments(Budget& work, const PgIntegrationPartition& partition,
                           const GramCutLayout& layout, const R& length, const R& four_n,
                           std::span<RankWide> table, std::span<RankWide> scratch) {
    auto& math = work.math;
    const auto two_length = math.mul(math.constant("2"), length);
    const auto slope_scale = math.mul(two_length, four_n);
    const auto domain = gram_thickness_domain(partition);
    const auto common = gram_thickness_scale(work, domain, four_n, two_length, slope_scale);
    // The L1-norm of the scaled affine coefficients is bounded by
    // 3 * max(abs(numerator)) * common.height. Include all five increments
    // in addition to the shared geometric bound and all visited edges.
    const auto height_bits =
        partition_mul(5, common.maximum_coefficient_bits + rank_bits(common.height) + 2);
    const auto scale = setup_cut_scale<thickness_polygon_degree>(
        work, common.chart, domain.indices.size(), height_bits);
    std::array<RankWide, 6> height_powers{};
    height_powers.front() = 1;
    for (std::size_t p = 1; p < height_powers.size(); ++p)
        height_powers[p] = rank_product(work, height_powers[p - 1], common.height);
    const auto axis_scale = 2 * layout.axis;

    for (const auto& piece : domain.pieces) {
        work.charge(scratch.size());
        std::ranges::fill(scratch, RankWide{});
        const auto indices = gram_piece_vertices(domain, piece);
        for (std::size_t j = 0; j < indices.size(); ++j) {
            work.charge();
            const auto& a = domain.vertices[static_cast<std::size_t>(indices[j])];
            const auto& b =
                domain.vertices[static_cast<std::size_t>(indices[(j + 1) % indices.size()])];
            const auto u0 = cut_coordinate(work, a.s, common.chart, axis_scale, piece.cell.first);
            const auto v0 = cut_coordinate(work, a.t, common.chart, axis_scale, piece.cell.second);
            const auto du =
                cut_coordinate(work, b.s, common.chart, axis_scale, piece.cell.first) - u0;
            const auto dv =
                cut_coordinate(work, b.t, common.chart, axis_scale, piece.cell.second) - v0;
            accumulate_cut_edge(work, scale, scratch, u0, du, v0, dv);
        }
        verify_cut_area(work, piece.area, axis_scale, scale.denominator, scratch.front());
        const auto local =
            local_half_thickness(math, domain.thickness, piece, four_n, two_length, slope_scale);
        std::array<RankWide, 3> affine{};
        for (std::size_t j = 0; j < affine.size(); ++j) {
            const auto fraction = cut_fraction(work, local[j]);
            const auto factor = gram_divide(work, common.height, fraction.denominator);
            affine[j] = rank_product(work, fraction.numerator, factor);
        }
        combine_thickness_moments(
            work,
            table.subspan(layout.index(piece.cell) * thickness_gram_slots, thickness_gram_slots),
            scratch, affine, height_powers);
    }
    return rank_product(work, scale.denominator, height_powers.back());
}

} // namespace uvlas::structure::detail
