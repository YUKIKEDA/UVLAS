#pragma once

#include "pg_reference_arithmetic.hpp"

#include <array>
#include <charconv>
#include <limits>

namespace uvlas::structure::detail {

using R = ExactRational;
using Jet2 = std::array<R, 6>;
inline constexpr std::array<ExactPowers, 15> reference_orders{{
    {0, 0},
    {1, 0},
    {0, 1},
    {2, 0},
    {1, 1},
    {0, 2},
    {3, 0},
    {2, 1},
    {1, 2},
    {0, 3},
    {4, 0},
    {3, 1},
    {2, 2},
    {1, 3},
    {0, 4},
}};

[[nodiscard]] inline R index_value(ReferenceArithmetic& m, core::Index value) {
    std::array<char, std::numeric_limits<core::Index>::digits10 + 3> chars{};
    const auto converted = std::to_chars(chars.data(), chars.data() + chars.size(), value);
    if (converted.ec != std::errc{})
        partition_fail(core::FailureCode::ResourceExceeded);
    return m.constant(
        std::string_view(chars.data(), static_cast<std::size_t>(converted.ptr - chars.data())));
}

[[nodiscard]] inline R times_monomial(ReferenceArithmetic& m, R value, core::Index exponent,
                                      core::Index order, const R& coordinate) {
    for (core::Index i = 0; i < order; ++i)
        value = m.mul(value, index_value(m, exponent - i));
    for (core::Index i = order; i < exponent; ++i)
        value = m.mul(value, coordinate);
    return value;
}

template <class Polynomial>
[[nodiscard]] R polynomial_derivative(ReferenceArithmetic& m, const Polynomial& polynomial,
                                      ExactPowers order, const R& s, const R& t,
                                      bool coefficient_bound = false) {
    R result;
    for (const auto& term : polynomial.terms()) {
        m.cancel();
        if (term.powers.first < order.first || term.powers.second < order.second)
            continue;
        auto value = coefficient_bound ? m.abs(term.coefficient) : m.copy(term.coefficient);
        value = times_monomial(m, std::move(value), term.powers.first, order.first, s);
        value = times_monomial(m, std::move(value), term.powers.second, order.second, t);
        result = m.add(result, value);
    }
    return result;
}

[[nodiscard]] inline std::array<R, 3>
physical_rotate(ReferenceArithmetic& m, const std::array<R, 9>& q, const std::array<R, 3>& x) {
    std::array<R, 3> result;
    for (std::size_t i = 0; i < result.size(); ++i)
        for (std::size_t j = 0; j < x.size(); ++j)
            result[i] = m.add(result[i], m.mul(q[i * x.size() + j], x[j]));
    return result;
}

[[nodiscard]] inline R determinant(ReferenceArithmetic& m, const std::array<R, 9>& q) {
    R result;
    for (std::size_t i = 0; i < 3; ++i) {
        const auto j = (i + 1) % 3, k = (i + 2) % 3;
        result =
            m.add(result, m.mul(q[i], m.sub(m.mul(q[3 + j], q[6 + k]), m.mul(q[3 + k], q[6 + j]))));
    }
    return result;
}

inline void validate_frame(ReferenceArithmetic& m, const PgReferenceFrame& frame) {
    if (frame.length_metres.sign() <= 0)
        partition_fail(core::FailureCode::GeometryInvalid);
    const auto one = m.constant("1");
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = i; j < 3; ++j) {
            R value;
            for (std::size_t k = 0; k < 3; ++k)
                value = m.add(value, m.mul(frame.rotation[k * 3 + i], frame.rotation[k * 3 + j]));
            if (i == j ? value != one : value.sign() != 0)
                partition_fail(core::FailureCode::GeometryInvalid);
        }
    if (determinant(m, frame.rotation).sign() <= 0)
        partition_fail(core::FailureCode::GeometryInvalid);
}

[[nodiscard]] inline R physical_hessian_squared_bound(ReferenceArithmetic& m,
                                                      const PgGraphView& graph, const R& max_s,
                                                      const R& max_t) {
    const auto ss = polynomial_derivative(m, graph.height, {2, 0}, max_s, max_t, true);
    const auto st = polynomial_derivative(m, graph.height, {1, 1}, max_s, max_t, true);
    const auto tt = polynomial_derivative(m, graph.height, {0, 2}, max_s, max_t, true);
    const auto inv_w = m.div(m.constant("1"), graph.width);
    const auto skew = m.mul(m.abs(graph.skew), inv_w);
    const auto xx = m.mul(ss, m.square(inv_w));
    const auto xy = m.mul(inv_w, m.add(st, m.mul(skew, ss)));
    const auto yy =
        m.add(m.add(tt, m.mul(m.constant("2"), m.mul(skew, st))), m.mul(m.square(skew), ss));
    const auto xy_squared = m.square(xy);
    return m.add(m.add(m.square(xx), m.add(xy_squared, xy_squared)), m.square(yy));
}

} // namespace uvlas::structure::detail
