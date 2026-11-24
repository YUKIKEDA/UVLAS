#pragma once

#include "pg_rank_storage.hpp"

namespace uvlas::structure::detail {

using PgConditioningFactor = RankInteger<3072>;

[[nodiscard]] inline std::size_t packed(std::size_t i, std::size_t j) noexcept {
    return j * (j + 1) / 2 + i;
}

template <class Work>
[[nodiscard]] RankWide conditioning_divide(Work& work, const RankWide& a, const RankWide& b) {
    work.charge(partition_mul(a.backend().size(), b.backend().size()));
    return RankWide{a / b};
}

template <class Work>
[[nodiscard]] RankWide conditioning_nearest(Work& work, const RankWide& a, const RankWide& b) {
    work.charge();
    const auto absolute = rank_abs(a);
    auto quotient = conditioning_divide(work, absolute, b);
    const auto twice_remainder = (absolute - rank_product(work, quotient, b)) << 1U;
    if (twice_remainder > b ||
        (twice_remainder == b && boost::multiprecision::bit_test(quotient, 0)))
        ++quotient;
    return a < 0 ? -quotient : quotient;
}

template <class Work>
[[nodiscard]] RankWide conditioning_rescale(Work& work, const RankWide& units, std::uint32_t from,
                                            std::uint32_t to) {
    work.charge();
    if (from > to)
        return conditioning_nearest(work, units, RankWide{1} << (from - to));
    auto shifted = rank_abs(units) << (to - from);
    return units < 0 ? -shifted : shifted;
}

/// Integer Newton iteration, started above the root. Charge each exact divide,
/// then use the exact integer halfway test to round a square root.
template <class Work>
[[nodiscard]] RankWide conditioning_root(Work& work, const RankWide& positive) {
    auto root = RankWide{1} << ((rank_bits(positive) + 1) / 2);
    for (;;) {
        work.charge();
        auto next = (root + conditioning_divide(work, positive, root)) >> 1U;
        if (next >= root)
            break;
        root = std::move(next);
    }
    if (positive > rank_product(work, root, root) + root)
        ++root;
    return root;
}

template <class Work, class Interval>
void conditioning_factor(Work& work, std::span<const Interval> gram, std::uint32_t gram_bits,
                         std::uint32_t factor_bits, std::size_t count,
                         std::span<PgConditioningFactor> factor) {
    for (std::size_t j = 0; j < count; ++j)
        for (std::size_t i = 0; i <= j; ++i) {
            work.charge();
            const auto index = packed(i, j);
            const auto& interval = gram[index];
            auto residual = conditioning_rescale(work, RankWide{interval.lower} + interval.upper,
                                                 gram_bits + 1, 2 * factor_bits);
            for (std::size_t k = 0; k < i; ++k) {
                work.charge();
                residual -= rank_product(work, RankWide{factor[packed(k, i)]},
                                         RankWide{factor[packed(k, j)]});
            }
            if (i == j) {
                if (residual <= 0)
                    throw PartitionAbort{
                        {core::FailureCode::PrecisionUnresolved, core::Operation::CertifyError}};
                factor[index] = PgConditioningFactor{conditioning_root(work, residual)};
            } else {
                const auto diagonal = RankWide{factor[packed(i, i)]};
                factor[index] =
                    PgConditioningFactor{conditioning_nearest(work, residual, diagonal)};
            }
        }
}

template <class Work>
void conditioning_inverse(Work& work, std::span<const PgConditioningFactor> factor,
                          std::uint32_t factor_bits, std::uint32_t candidate_bits,
                          std::span<PgConditioningFactor> column, std::span<RankCell> candidate) {
    const auto one = RankWide{1} << (2 * factor_bits);
    for (std::size_t j = 0; j < column.size(); ++j)
        for (std::size_t rows = j + 1; rows > 0; --rows) {
            work.charge();
            const auto i = rows - 1;
            auto residual = i == j ? one : RankWide{};
            for (std::size_t k = i + 1; k <= j; ++k) {
                work.charge();
                residual -= rank_product(work, RankWide{factor[packed(i, k)]}, RankWide{column[k]});
            }
            column[i] = PgConditioningFactor{
                conditioning_nearest(work, residual, RankWide{factor[packed(i, i)]})};
            candidate[packed(i, j)] = RankCell{
                conditioning_rescale(work, RankWide{column[i]}, factor_bits, candidate_bits)};
        }
}

} // namespace uvlas::structure::detail
