#pragma once

#include "pg_reference_arithmetic.hpp"

#include <bit>
#include <limits>
#include <uvlas/core/error_certificate.hpp>

namespace uvlas::structure::detail {

inline constexpr std::uint64_t pg_binary64_sign = 0x8000000000000000ULL;
inline constexpr std::uint64_t pg_binary64_infinity = 0x7ff0000000000000ULL;
inline constexpr std::uint64_t pg_binary64_max_finite = pg_binary64_infinity - 1;
static_assert(std::numeric_limits<double>::is_iec559 && sizeof(double) == sizeof(std::uint64_t));

struct PgProjectedWord {
    double high{};
    double low{};
};

[[nodiscard]] inline ExactRational reference_binary64(ReferenceArithmetic& math, double value) {
    math.charge();
    return partition_take(ExactRational::from_binary64(value, math.leaf_limits, math.plan));
}

/// Compare a + b sqrt(g) with a rational dyadic threshold. The square of the
/// radical term is computed once and is compared only when signs oppose.
struct PgProjectionComparator {
    ReferenceArithmetic& math;
    const ExactRational& square_radical;
    int radical_sign;

    [[nodiscard]] int compare(const ExactRational& a, const ExactRational& threshold) const {
        if (radical_sign == 0) {
            const auto comparison = math.compare(a, threshold);
            return comparison == 0 ? 0 : (comparison > 0 ? 1 : -1);
        }
        const auto delta = math.sub(a, threshold);
        const auto sign = delta.sign();
        if (sign == 0 || sign == radical_sign)
            return radical_sign;
        const auto comparison = math.compare(math.square(delta), square_radical);
        return comparison == 0 ? 0 : (comparison > 0 ? sign : radical_sign);
    }
};

template <class Comparator>
[[nodiscard]] double pg_round_word(const Comparator& comparator, const ExactRational& a) {
    const auto sign = comparator.compare(a, {});
    if (sign == 0)
        return 0.0;
    auto& math = comparator.math;
    const auto word = [mask = sign < 0 ? pg_binary64_sign : std::uint64_t{0}](std::uint64_t n) {
        return std::bit_cast<double>(n | mask);
    };
    const auto compare_word = [&](std::uint64_t value) {
        return sign * comparator.compare(a, reference_binary64(math, word(value)));
    };
    if (compare_word(pg_binary64_max_finite) >= 0) {
        const auto max = reference_binary64(math, word(pg_binary64_max_finite));
        const auto gap = reference_binary64(math, word(std::bit_cast<std::uint64_t>(0x1p970)));
        if (sign * comparator.compare(a, math.add(max, gap)) >= 0)
            partition_fail(core::FailureCode::PrecisionUnresolved);
        return word(pg_binary64_max_finite);
    }
    std::uint64_t lower = 0;
    auto upper = pg_binary64_max_finite;
    while (upper - lower > 1) {
        const auto middle = lower + (upper - lower) / 2;
        if (compare_word(middle) < 0)
            upper = middle;
        else
            lower = middle;
    }
    const auto x = reference_binary64(math, word(lower));
    const auto y = reference_binary64(math, word(upper));
    const auto midpoint = math.div(math.add(x, y), math.constant("2"));
    const auto relation = sign * comparator.compare(a, midpoint);
    return word(relation < 0 || (relation == 0 && lower % 2 == 0) ? lower : upper);
}

template <class Comparator>
[[nodiscard]] PgProjectedWord pg_project_words(const Comparator& comparator, const ExactRational& a,
                                               const ExactRational& tolerance,
                                               core::Precision representation) {
    PgProjectedWord words;
    auto& math = comparator.math;
    words.high = pg_round_word(comparator, a);
    auto residual = math.sub(a, reference_binary64(math, words.high));
    if (representation == core::Precision::TwoWord) {
        words.low = pg_round_word(comparator, residual);
        residual = math.sub(residual, reference_binary64(math, words.low));
    }
    if (comparator.compare(residual, tolerance) > 0 ||
        comparator.compare(residual, math.neg(tolerance)) < 0)
        partition_fail(core::FailureCode::PrecisionUnresolved);
    return words;
}

[[nodiscard]] inline PgProjectedWord
pg_project_scalar(ReferenceArithmetic& math, const PgReferenceScalar& value, const ExactRational& g,
                  const ExactRational& tolerance, core::Precision representation) {
    const auto radical_sign = value.radical.sign();
    if (value.rational.sign() == 0 && radical_sign == 0)
        return {};
    const auto square =
        radical_sign == 0 ? ExactRational{} : math.mul(g, math.square(value.radical));
    const PgProjectionComparator comparator{math, square, radical_sign};
    return pg_project_words(comparator, value.rational, tolerance, representation);
}

/// Compare offset + sqrt(squared_norm) to a rational without rounding the area.
struct PgPositiveRootComparator {
    ReferenceArithmetic& math;
    const PgProjectionComparator& squared;
    const ExactRational& rational;

    [[nodiscard]] int compare(const ExactRational& offset, const ExactRational& threshold) const {
        const auto delta = math.sub(threshold, offset);
        return delta.sign() <= 0 ? 1 : squared.compare(rational, math.square(delta));
    }
};

[[nodiscard]] inline PgProjectedWord pg_project_positive_root(ReferenceArithmetic& math,
                                                              const PgReferenceScalar& square,
                                                              const ExactRational& g,
                                                              const ExactRational& tolerance,
                                                              core::Precision representation) {
    const auto sign = square.radical.sign();
    const auto squared_radical =
        sign == 0 ? ExactRational{} : math.mul(g, math.square(square.radical));
    const PgProjectionComparator squared{math, squared_radical, sign};
    if (squared.compare(square.rational, {}) <= 0)
        partition_fail(core::FailureCode::PrecisionUnresolved);
    const PgPositiveRootComparator root{math, squared, square.rational};
    return pg_project_words(root, {}, tolerance, representation);
}

} // namespace uvlas::structure::detail
