#pragma once

#include "pg_partition_common.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <uvlas/structure/pg_rank_preparation.hpp>

namespace uvlas::structure::detail {

template <unsigned Bits>
using RankInteger = boost::multiprecision::number<
    boost::multiprecision::cpp_int_backend<Bits, Bits, boost::multiprecision::signed_magnitude,
                                           boost::multiprecision::checked, void>,
    boost::multiprecision::et_off>;

using RankCell = RankInteger<1536>;
using RankWide = RankInteger<8192>;
using RankWord = boost::multiprecision::limb_type;
constexpr std::uint32_t rank_required = 33;
constexpr std::uint32_t rank_max_columns = 933;
constexpr std::uint32_t rank_max_original = 530;
constexpr std::uint32_t rank_max_candidates = rank_required + rank_max_columns + rank_max_original;
constexpr std::uint32_t rank_max_cells = 64;
constexpr std::array<std::uint32_t, 3> rank_degree{13, 21, 2};
constexpr std::array<std::uint32_t, 4> rank_field_slots{0, 105, 358, 364};
constexpr std::array<std::uint32_t, 4> rank_trace_slots{0, 14, 36, 39};
constexpr std::uint32_t rank_key_count = rank_field_slots.back() * rank_max_cells;
constexpr std::uint32_t rank_trace_count = rank_trace_slots.back() * 8;
constexpr std::uint32_t rank_missing = UINT32_MAX;

template <class I>
struct RankFraction {
    I numerator{};
    I denominator{1};
};

struct RankBudget {
    const PgRankPreparationLimits& limits;
    std::uint64_t steps = 0;

    void charge(std::uint64_t units = 1) {
        if (limits.cancellation.stop_requested())
            partition_fail(core::FailureCode::Cancelled);
        if (units > limits.maximum_steps - steps)
            partition_fail(core::FailureCode::ResourceExceeded);
        steps += units;
    }
};

template <class I>
[[nodiscard]] std::uint64_t rank_bits(const I& value) {
    if (value == 0)
        return 0;
    const auto& backend = value.backend();
    const auto high = backend.limbs()[backend.size() - 1];
    constexpr auto limb_bits = static_cast<std::uint64_t>(std::numeric_limits<RankWord>::digits);
    return std::uint64_t{backend.size() - 1} * limb_bits +
           static_cast<std::uint64_t>(std::bit_width(high));
}

template <class I>
[[nodiscard]] I rank_abs(const I& value) {
    return value < 0 ? I{-value} : I{value};
}

/// Horner over native limbs. Only a single native limb is multiplied by Boost
/// at a time, so neither this path nor the one-limb path can use Karatsuba.
template <class I>
[[nodiscard]] I rank_product(RankBudget& budget, const I& left, const I& right) {
    const auto* a = &left;
    const auto* b = &right;
    if (a->backend().size() < b->backend().size())
        std::swap(a, b);
    const auto count = b->backend().size();
    budget.charge(partition_mul(a->backend().size(), count));
    if (count == 1) {
        auto result = *a;
        result *= b->backend().limbs()[0];
        if (*b < 0)
            result = -result;
        return result;
    }
    const I magnitude = rank_abs(*a);
    I result{0};
    for (auto i = count; i > 0; --i) {
        result <<= std::numeric_limits<RankWord>::digits;
        auto next = magnitude;
        next *= b->backend().limbs()[i - 1];
        result += next;
    }
    if ((*a < 0) != (*b < 0))
        result = -result;
    return result;
}

template <class I>
[[nodiscard]] I rank_gcd(RankBudget& budget, I a, I b) {
    a = rank_abs(a);
    b = rank_abs(b);
    while (b != 0) {
        budget.charge(partition_mul(a.backend().size(), b.backend().size()));
        auto r = a % b;
        a = std::move(b);
        b = std::move(r);
    }
    return a;
}

template <class I>
[[nodiscard]] I rank_lcm(RankBudget& budget, const I& left, const I& right) {
    if (left == 1 || right == 1)
        return left == 1 ? right : left;
    const I divisor = rank_gcd(budget, left, right);
    const I reduced = left / divisor;
    return rank_product(budget, reduced, right);
}

template <class I>
[[nodiscard]] RankFraction<I> normalized_rank_fraction(RankBudget& budget, I n, I d) {
    if (d == 0)
        partition_fail(core::FailureCode::RankUnresolved);
    if (n == 0)
        return {};
    const I common = rank_gcd(budget, n, d);
    budget.charge(partition_mul(common.backend().size(),
                                partition_add(n.backend().size(), d.backend().size())));
    n /= common;
    d /= common;
    if (d < 0) {
        n = -n;
        d = -d;
    }
    return {std::move(n), std::move(d)};
}

template <class I>
[[nodiscard]] I rank_copy_integer(const IntegerBackend& source) {
    I result;
    auto& target = result.backend();
    target.resize(source.size(), source.size());
    std::copy_n(source.limbs(), source.size(), target.limbs());
    target.sign(source.sign());
    target.normalize();
    return result;
}

template <class I>
[[nodiscard]] RankFraction<I> rank_fraction(const ExactRational& exact) {
    const auto& source = ExactFactory::value(exact).backend();
    return {rank_copy_integer<I>(source.num()), rank_copy_integer<I>(source.denom())};
}

struct RankSizer {
    std::uint64_t bytes = 0;

    template <class T>
    void add(std::uint64_t count) {
        bytes = partition_add(bytes, partition_mul(count, sizeof(T)));
        static_cast<void>(partition_index(count));
    }
};

template <class T>
[[nodiscard]] core::PlannedBuffer<T> rank_allocate(core::ResourcePlan& plan, std::uint64_t count) {
    return partition_take(plan.allocate<T>(partition_index(count), partition_operation));
}

template <class T>
[[nodiscard]] std::span<T> rank_span(core::PlannedBuffer<T>& values) noexcept {
    return values.values();
}
template <class T>
[[nodiscard]] std::span<const T> rank_span(const core::PlannedBuffer<T>& values) noexcept {
    return values.values();
}

struct RankCellEntry {
    std::uint32_t row = 0;
    RankCell value;
};

struct RankCsc {
    core::PlannedBuffer<std::uint64_t> offsets;
    core::PlannedBuffer<RankCellEntry> entries;

    [[nodiscard]] std::span<const RankCellEntry> column(std::uint32_t j) const noexcept {
        const auto begin = offsets.values()[j];
        return entries.values().subspan(static_cast<std::size_t>(begin),
                                        static_cast<std::size_t>(offsets.values()[j + 1] - begin));
    }
};

[[nodiscard]] inline RankCsc rank_csc(core::ResourcePlan& plan, std::uint64_t columns,
                                      std::uint64_t terms) {
    return {rank_allocate<std::uint64_t>(plan, columns + 1),
            rank_allocate<RankCellEntry>(plan, terms)};
}

struct RankInputMap {
    std::uint32_t source = 0;
    bool negative = false;
    const ExactRational* coefficient = nullptr; ///< Null represents unit magnitude.
};

struct RankSourceMap {
    std::uint32_t source = 0;
    RankFraction<RankCell> coefficient;
};

} // namespace uvlas::structure::detail
