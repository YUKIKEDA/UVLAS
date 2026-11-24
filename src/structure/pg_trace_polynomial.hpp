#pragma once

#include "pg_reference_polynomials.hpp"

#include <algorithm>
#include <array>

namespace uvlas::structure::detail {

inline constexpr std::size_t pg_trace_max_degree = 30;
inline constexpr std::array<std::size_t, 3> pg_trace_field_slots{14, 31, 3};
inline constexpr std::size_t pg_trace_row_count = 48;
inline constexpr std::size_t pg_trace_basis_max_degree = 21;

/// Ascending powers in the local parameter of one common external segment.
struct PgTracePolynomial {
    std::array<ExactRational, pg_trace_max_degree + 1> coefficients{};
    std::size_t degree = 0;

    [[nodiscard]] bool zero() const noexcept {
        return degree == 0 && coefficients[0].sign() == 0;
    }
    void trim() noexcept {
        while (degree != 0 && coefficients[degree].sign() == 0)
            --degree;
    }
};

[[nodiscard]] inline PgTracePolynomial
trace_affine(ReferenceArithmetic& math, const ExactRational& begin, const ExactRational& end) {
    PgTracePolynomial output;
    output.coefficients[0] = math.copy(begin);
    output.coefficients[1] = math.sub(end, begin);
    output.degree = 1;
    output.trim();
    return output;
}

inline void trace_add_scaled(ReferenceArithmetic& math, PgTracePolynomial& sum,
                             const PgTracePolynomial& value, const ExactRational& scale) {
    if (scale.sign() == 0 || value.zero())
        return;
    sum.degree = std::max(sum.degree, value.degree);
    for (std::size_t k = 0; k <= value.degree; ++k) {
        math.visit(1);
        if (value.coefficients[k].sign() != 0)
            sum.coefficients[k] =
                math.add(sum.coefficients[k], math.mul(value.coefficients[k], scale));
    }
    sum.trim();
}

[[nodiscard]] inline PgTracePolynomial trace_product(ReferenceArithmetic& math,
                                                     const PgTracePolynomial& left,
                                                     const PgTracePolynomial& right) {
    PgTracePolynomial result;
    if (left.zero() || right.zero())
        return result;
    result.degree = left.degree + right.degree;
    if (result.degree > pg_trace_max_degree)
        partition_fail(core::FailureCode::ResourceExceeded);
    for (std::size_t i = 0; i <= left.degree; ++i) {
        math.visit(1);
        if (left.coefficients[i].sign() == 0)
            continue;
        for (std::size_t j = 0; j <= right.degree; ++j) {
            math.visit(1);
            if (right.coefficients[j].sign() != 0)
                result.coefficients[i + j] =
                    math.add(result.coefficients[i + j],
                             math.mul(left.coefficients[i], right.coefficients[j]));
        }
    }
    result.trim();
    return result;
}

[[nodiscard]] inline bool trace_same(ReferenceArithmetic& math, const PgTracePolynomial& left,
                                     const PgTracePolynomial& right) {
    if (left.degree != right.degree)
        return false;
    for (std::size_t i = 0; i <= left.degree; ++i) {
        math.visit(1);
        if (math.compare(left.coefficients[i], right.coefficients[i]) != 0)
            return false;
    }
    return true;
}

using PgTracePowers = std::array<PgTracePolynomial, pg_trace_basis_max_degree + 1>;

[[nodiscard]] inline PgTracePowers
trace_powers(ReferenceArithmetic& math, const PgTracePolynomial& affine, std::size_t through) {
    PgTracePowers result;
    result[0].coefficients[0] = math.constant("1");
    for (std::size_t i = 1; i <= through; ++i)
        result[i] = trace_product(math, result[i - 1], affine);
    return result;
}

[[nodiscard]] inline PgTracePolynomial
trace_restrict(ReferenceArithmetic& math, const ExactChartPolynomial& polynomial,
               const PgTracePowers& s, const PgTracePowers& t, ExactPowers order = {}) {
    PgTracePolynomial result;
    for (const auto& term : polynomial.terms()) {
        math.visit(1);
        if (term.powers.first < order.first || term.powers.second < order.second)
            continue;
        auto multiplier = math.copy(term.coefficient);
        if (order.first == 1)
            multiplier = math.mul(multiplier, index_value(math, term.powers.first));
        if (order.second == 1)
            multiplier = math.mul(multiplier, index_value(math, term.powers.second));
        const auto a = static_cast<std::size_t>(term.powers.first - order.first);
        const auto b = static_cast<std::size_t>(term.powers.second - order.second);
        trace_add_scaled(math, result, trace_product(math, s[a], t[b]), multiplier);
    }
    return result;
}

} // namespace uvlas::structure::detail
