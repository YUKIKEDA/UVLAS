#pragma once

#include "pg_coefficients_detail.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cstdint>
#include <span>
#include <string_view>
#include <utility>

namespace uvlas::structure::detail {
using PgPoly = ExactCellPolynomial;
struct PgCoefficientAbort {
    core::Failure failure;
};
template <class T>
T pg_take(core::Result<T> value) {
    if (!value)
        throw PgCoefficientAbort{value.error()};
    return std::move(*value);
}
inline void pg_take(core::Result<void> value) {
    if (!value)
        throw PgCoefficientAbort{value.error()};
}
/// All native and exact leaf calls pass through this preparation-local adapter.
/// Nested timing is disabled; no caller-owned input or mutable global cache is used.
struct PgCoefficientAlgebra {
    PgCoefficientAlgebra(const PgCoefficientLimits& caps, core::ResourcePlan& payload)
        : limits(caps.algebra), plan(payload), cancellation(caps.cancellation) {
        limits.timing = ExactTimingMode::Disabled;
    }
    ExactLimits limits;
    core::ResourcePlan& plan;
    std::stop_token cancellation;
    ExactRational zero;
    ExactPreparationInfo leaf_maximum;
    std::array<std::uint64_t, 10> calls{};
    void check() const {
        if (cancellation.stop_requested())
            throw PgCoefficientAbort{{core::FailureCode::Cancelled, core::Operation::PrepareModel}};
    }
    template <class T>
    T finish(core::Result<T> value, std::size_t operation) {
        auto result = pg_take(std::move(value));
        const auto& info = result.preparation();
        leaf_maximum.algebraic_work_bits =
            std::max(leaf_maximum.algebraic_work_bits, info.algebraic_work_bits);
        leaf_maximum.planned_work_bits =
            std::max(leaf_maximum.planned_work_bits, info.planned_work_bits);
        leaf_maximum.planned_raw_terms =
            std::max(leaf_maximum.planned_raw_terms, info.planned_raw_terms);
        leaf_maximum.estimated_external_payload_bytes = std::max(
            leaf_maximum.estimated_external_payload_bytes, info.estimated_external_payload_bytes);
        ++calls[operation];
        check();
        return result;
    }
    ExactRational number(std::int64_t numerator, std::int64_t denominator = 1) {
        check();
        std::array<char, 32> n{}, d{};
        const auto nr = std::to_chars(n.data(), n.data() + n.size(), numerator);
        const auto dr = std::to_chars(d.data(), d.data() + d.size(), denominator);
        if (nr.ec != std::errc{} || dr.ec != std::errc{})
            throw PgCoefficientAbort{
                {core::FailureCode::ResourceExceeded, core::Operation::PrepareModel}};
        return finish(ExactRational::from_integer_pair(
                          {n.data(), static_cast<std::size_t>(nr.ptr - n.data())},
                          {d.data(), static_cast<std::size_t>(dr.ptr - d.data())}, limits, plan),
                      0);
    }
    ExactRational rm(const ExactRational& a, const ExactRational& b) {
        return finish(a.multiply(b, limits, plan), 0);
    }
    ExactRational rd(const ExactRational& a, const ExactRational& b) {
        return finish(a.divide(b, limits, plan), 0);
    }
    ExactRational ra(const ExactRational& a, const ExactRational& b) {
        return finish(a.add(b, limits, plan), 0);
    }
    ExactRational rs(const ExactRational& a, const ExactRational& b) {
        return finish(a.subtract(b, limits, plan), 0);
    }
    ExactRational rk(const ExactRational& a) {
        return finish(a.clone(limits, plan), 0);
    }
    PgPoly make(std::span<const ExactInputTerm> terms) {
        return finish(PgPoly::make(terms, limits, plan), 1);
    }
    PgPoly a(const PgPoly& x, const PgPoly& y) {
        return finish(x.add(y, limits, plan), 2);
    }
    PgPoly sub(const PgPoly& x, const PgPoly& y) {
        return finish(x.subtract(y, limits, plan), 2);
    }
    PgPoly m(const PgPoly& x, const PgPoly& y) {
        return finish(x.multiply(y, limits, plan), 3);
    }
    PgPoly s(const PgPoly& x, const ExactRational& y) {
        return finish(x.scale(y, limits, plan), 4);
    }
    PgPoly d(const PgPoly& x, ExactAxis axis, core::Index order = 1) {
        return finish(x.differentiate(axis, order, limits, plan), 5);
    }
    PgPoly i(const PgPoly& x, ExactAxis axis, const ExactRational& lower) {
        return finish(x.integrate(axis, lower, limits, plan), 6);
    }
    ExactRational e(const PgPoly& x, const ExactRational& u, const ExactRational& v) {
        return finish(x.evaluate(u, v, limits, plan), 7);
    }
    PgPoly k(const PgPoly& x) {
        return finish(x.clone(limits, plan), 8);
    }
    PgPoly constant(const ExactRational& value) {
        if (value.sign() == 0)
            return {};
        const std::array terms{ExactInputTerm{{0, 0}, value}};
        return make(terms);
    }
    PgPoly constant(std::int64_t value) {
        const auto r = number(value);
        return constant(r);
    }
    PgPoly affine(const ExactRational& origin, const ExactRational& slope, ExactAxis axis) {
        const ExactPowers powers = axis == ExactAxis::First ? ExactPowers{1, 0} : ExactPowers{0, 1};
        const std::array terms{ExactInputTerm{{0, 0}, origin}, ExactInputTerm{powers, slope}};
        return make(terms);
    }
};
template <std::size_t N, class F, std::size_t... I>
std::array<ExactInputTerm, N> pg_inputs(F&& get, std::index_sequence<I...>) {
    return {get(I)...};
}
template <std::size_t N, class F>
std::array<ExactInputTerm, N> pg_inputs(F&& get) {
    return pg_inputs<N>(std::forward<F>(get), std::make_index_sequence<N>{});
}
inline bool pg_equal(const PgPoly& a, const PgPoly& b) noexcept {
    const auto x = a.terms(), y = b.terms();
    if (x.size() != y.size())
        return false;
    for (std::size_t i = 0; i < x.size(); ++i)
        if (x[i].powers.first != y[i].powers.first || x[i].powers.second != y[i].powers.second ||
            !(x[i].coefficient == y[i].coefficient))
            return false;
    return true;
}
/// Exact numerator restriction at the original chart anchor. It leaves v symbolic.
/// Coefficients remain borrowed until make() returns; all negated values are owned.
inline PgPoly pg_anchor(PgCoefficientAlgebra& op, const PgPoly& value) {
    const auto terms = value.terms();
    if (terms.empty())
        return {};
    if (terms.size() > 253)
        throw PgCoefficientAbort{
            {core::FailureCode::ResourceExceeded, core::Operation::PrepareModel}};
    std::array<ExactRational, 253> negatives;
    const auto minus = op.number(-1);
    for (std::size_t j = 0; j < terms.size(); ++j)
        if ((terms[j].powers.first & 1) != 0)
            negatives[j] = op.rm(terms[j].coefficient, minus);
    const auto inputs = pg_inputs<253>([&](std::size_t j) -> ExactInputTerm {
        if (j >= terms.size())
            return {{0, 0}, op.zero};
        return {{0, terms[j].powers.second},
                (terms[j].powers.first & 1) != 0 ? negatives[j] : terms[j].coefficient};
    });
    return op.make(std::span{inputs}.first(terms.size()));
}
inline std::array<PgPoly, 5> pg_powers(PgCoefficientAlgebra& op, const PgPoly& x) {
    std::array<PgPoly, 5> result;
    result[0] = op.constant(1);
    result[1] = op.k(x);
    for (std::size_t j = 2; j < result.size(); ++j)
        result[j] = op.m(result[j - 1], x);
    return result;
}
inline std::array<PgPoly, 8> pg_legendre(PgCoefficientAlgebra& op, ExactAxis axis) {
    std::array<PgPoly, 8> result;
    const auto one = op.number(1);
    const auto x = op.affine(op.zero, one, axis);
    result[0] = op.constant(one);
    result[1] = op.k(x);
    for (std::int64_t degree = 1; degree < 7; ++degree) {
        const auto left =
            op.s(op.m(x, result[static_cast<std::size_t>(degree)]), op.number(2 * degree + 1));
        const auto right = op.s(result[static_cast<std::size_t>(degree - 1)], op.number(degree));
        result[static_cast<std::size_t>(degree + 1)] =
            op.s(op.sub(left, right), op.number(1, degree + 1));
    }
    return result;
}
struct PgMomentRow {
    std::array<PgPoly, 5> active, completed;
};
inline PgMomentRow pg_moments(PgCoefficientAlgebra& op, core::Index n, core::Index origin,
                              const PgPoly& curvature, ExactAxis axis) {
    const auto h = op.number(1, 4 * n), center = op.number(2 * origin + 1, 4 * n);
    const auto coordinate = op.affine(center, h, axis);
    const auto powers = pg_powers(op, coordinate);
    const auto lower = op.number(-1), upper = op.number(1);
    PgMomentRow result;
    for (std::size_t j = 0; j < 5; ++j) {
        result.active[j] = op.i(op.s(op.m(powers[j], curvature), h), axis, lower);
        const auto endpoint = axis == ExactAxis::First ? op.e(result.active[j], upper, op.zero)
                                                       : op.e(result.active[j], op.zero, upper);
        result.completed[j] = op.constant(endpoint);
    }
    return result;
}
/// Vector lifting from the original global f coefficients and local coordinate powers.
/// Moment carry is already selected by origin/cell; no per-cell constant reset occurs.
inline std::array<PgPoly, 3>
pg_lifting(PgCoefficientAlgebra& op, const ExactChartPolynomial& global_f, const PgPoly& local_f,
           const PgPoly& along, const std::array<PgPoly, 5>& transverse,
           const std::array<PgPoly, 5>& moments, PgPairedAxis axis,
           const ExactRational& inverse_width, const ExactRational& skew) {
    PgPoly J, Jt, Jst;
    for (const auto& term : global_f.terms()) {
        const auto i = static_cast<std::size_t>(axis == PgPairedAxis::S ? term.powers.first
                                                                        : term.powers.second);
        const auto j = static_cast<std::size_t>(axis == PgPairedAxis::S ? term.powers.second
                                                                        : term.powers.first);
        J = op.a(J, op.s(op.m(transverse[j], moments[i]), term.coefficient));
        if (j != 0) {
            const auto coefficient =
                op.rm(term.coefficient, op.number(static_cast<std::int64_t>(j)));
            Jt = op.a(Jt, op.s(op.m(transverse[j - 1], moments[i]), coefficient));
            Jst = op.a(Jst, op.s(op.m(transverse[j - 1], moments[i + 1]), coefficient));
        }
    }
    auto first = op.sub(J, op.m(local_f, moments[0]));
    auto second = op.sub(Jst, op.m(along, Jt));
    auto third = op.sub(op.m(along, moments[0]), moments[1]);
    if (axis == PgPairedAxis::T)
        std::swap(first, second);
    first = op.s(first, inverse_width);
    second = op.sub(second, op.s(first, skew));
    return {std::move(first), std::move(second), std::move(third)};
}
/// Numerators use the same physical-plane derivatives as the adopted v field.
/// g is not square-rooted, and no normalized director is fitted at sample points.
inline std::array<PgPoly, 3> pg_director(PgCoefficientAlgebra& op, const std::array<PgPoly, 3>& v,
                                         const PgPoly& fx, const PgPoly& fy, const PgPoly& U,
                                         const PgPoly& V, const PgPoly& xy,
                                         const ExactRational& chart_scale,
                                         const ExactRational& inverse_width,
                                         const ExactRational& skew) {
    const auto scale_x = op.rm(chart_scale, inverse_width);
    std::array<PgPoly, 3> dx, dy;
    for (std::size_t j = 0; j < 3; ++j) {
        dx[j] = op.s(op.d(v[j], ExactAxis::First), scale_x);
        dy[j] = op.sub(op.s(op.d(v[j], ExactAxis::Second), chart_scale), op.s(dx[j], skew));
    }
    const auto tx = op.sub(op.sub(dx[2], op.m(fx, dx[0])), op.m(fy, dx[1]));
    const auto ty = op.sub(op.sub(dy[2], op.m(fx, dy[0])), op.m(fy, dy[1]));
    auto px = op.sub(op.m(xy, ty), op.m(V, tx));
    auto py = op.sub(op.m(xy, tx), op.m(U, ty));
    const PgPoly zero;
    auto pz = op.sub(zero, op.a(op.m(fx, tx), op.m(fy, ty)));
    return {std::move(px), std::move(py), std::move(pz)};
}
} // namespace uvlas::structure::detail
