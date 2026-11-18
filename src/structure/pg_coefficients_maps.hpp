#pragma once

#include "pg_coefficients_anchor.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace uvlas::structure::detail {
inline const ExactRational* pg_power_coefficient(const PgPoly& value, core::Index degree) noexcept {
    for (const auto& term : value.terms())
        if (term.powers.first == degree && term.powers.second == 0)
            return &term.coefficient;
    return nullptr;
}
inline PgPoly pg_global_curvature(PgCoefficientAlgebra& op, const PgPoly& legendre, core::Index n,
                                  core::Index origin) {
    const auto coordinate =
        op.affine(op.number(2 * origin + 1 - 2 * n, 2 * n), op.number(1, 2 * n), ExactAxis::First);
    std::array<PgPoly, 8> powers;
    powers[0] = op.constant(1);
    powers[1] = op.k(coordinate);
    for (std::size_t j = 2; j < 8; ++j)
        powers[j] = op.m(powers[j - 1], coordinate);
    PgPoly result;
    for (const auto& term : legendre.terms())
        result = op.a(result,
                      op.s(powers[static_cast<std::size_t>(term.powers.first)], term.coefficient));
    return result;
}
struct PgCoefficientMapBuilder {
    PgCoefficientGeneration& generated;
    void add(core::Index source, std::uint8_t target, ExactRational coefficient) {
        auto& g = generated;
        auto& data = g.data;
        if (coefficient.sign() == 0)
            return;
        if (data.required_count >= data.required.size())
            g.exceeded();
        data.required.values()[data.required_count++] = {{source}, target, std::move(coefficient)};
    }
    void required(const PgGraphView& graph) {
        auto& g = generated;
        auto& data = g.data;
        auto& op = g.op;
        const auto n = data.n;
        const auto nr = (n + 2) * (n + 2) + 1, nd = (2 * n + 1) * (2 * n + 1) + 3;
        std::array<core::Index, 9> knots{};
        std::size_t size = 3;
        for (core::Index j = 1; j < n; ++j)
            knots[size++] = j;
        for (std::size_t j = 0; j < 3; ++j)
            knots[size++] = n;
        std::array<ExactRational, 6> greville;
        for (core::Index j = 0; j < n + 2; ++j)
            greville[static_cast<std::size_t>(j)] = op.number(
                knots[static_cast<std::size_t>(j + 1)] + knots[static_cast<std::size_t>(j + 2)],
                2 * n);
        const auto one = op.number(1);
        for (core::Index i = 0; i < n + 2; ++i)
            for (core::Index j = 0; j < n + 2; ++j) {
                const auto source = i * (n + 2) + j;
                const auto x = op.ra(
                    op.ra(graph.x0, op.rm(graph.width, greville[static_cast<std::size_t>(i)])),
                    op.rm(graph.skew, greville[static_cast<std::size_t>(j)]));
                add(source, 0, op.rk(x));
                add(source, 1, op.rk(greville[static_cast<std::size_t>(j)]));
                for (std::uint8_t target = 30; target < 33; ++target)
                    add(source, target, op.rk(one));
            }
        add(nr - 1, 2, op.rk(one));
        for (std::uint8_t j = 0; j < 3; ++j)
            add(nr + nd - 3 + j, static_cast<std::uint8_t>(27 + j), op.rk(one));
        const auto legendre = pg_legendre(op, ExactAxis::First);
        for (core::Index origin = 0; origin < 2 * n; ++origin)
            for (core::Index degree = 0; degree < 8; ++degree) {
                auto remainder =
                    pg_global_curvature(op, legendre[static_cast<std::size_t>(degree)], n, origin);
                for (core::Index local = degree; local >= 0; --local) {
                    const auto* numerator = pg_power_coefficient(remainder, local);
                    if (numerator == nullptr || numerator->sign() == 0)
                        continue;
                    const auto& basis = legendre[static_cast<std::size_t>(local)];
                    const auto* denominator = pg_power_coefficient(basis, local);
                    pg_require_identity(denominator != nullptr);
                    const auto coefficient = op.rd(*numerator, *denominator);
                    remainder = op.sub(remainder, op.s(basis, coefficient));
                    for (std::uint8_t j = 0; j < 3; ++j) {
                        const auto source =
                            nr + nd + static_cast<core::Index>(j) * (32 * n) + origin * 8 + local;
                        const auto target = static_cast<std::uint8_t>(3 + j * 8 + degree);
                        add(source, target, op.rk(coefficient));
                    }
                }
                pg_require_identity(remainder.terms().empty());
            }
        auto entries = data.required.values().first(data.required_count);
        std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
            return a.target < b.target || (a.target == b.target && a.source.value < b.source.value);
        });
    }
    PgPoly reconstruct(std::uint8_t target, PgCellIndex cell, PgCoefficientField field) {
        auto& g = generated;
        auto& op = g.op;
        const auto entries = g.data.required.values().first(g.data.required_count);
        const auto start =
            std::lower_bound(entries.begin(), entries.end(), target,
                             [](const auto& entry, std::uint8_t id) { return entry.target < id; });
        PgPoly value;
        for (auto it = start; it != entries.end() && it->target == target; ++it) {
            const auto& coefficient = g.get(it->source.value, cell, field);
            if (!coefficient.terms().empty())
                value = op.a(value, op.s(coefficient, it->coefficient));
        }
        return value;
    }
    void check_field(std::uint8_t target, PgCellIndex cell, PgCoefficientField field,
                     const PgPoly& expected) {
        const auto actual = reconstruct(target, cell, field);
        pg_require_identity(pg_equal(actual, expected));
    }
    void verify_basic(const PgGraphView& graph) {
        auto& g = generated;
        auto& data = g.data;
        auto& op = g.op;
        const auto axis = 2 * data.n;
        const auto minus = op.number(-1);
        const auto unit = op.constant(1);
        const auto x0 = op.constant(graph.x0);
        const PgPoly zero;
        for (core::Index cs = 0; cs < axis; ++cs)
            for (core::Index ct = 0; ct < axis; ++ct) {
                const PgCellIndex cell{cs, ct};
                const auto index = static_cast<std::size_t>(cs * axis + ct);
                const auto& c = g.cells.values()[index];
                const auto x = op.a(op.a(x0, op.s(c.s, graph.width)), op.s(c.t, graph.skew));
                check_field(0, cell, PgCoefficientField::R, x);
                check_field(1, cell, PgCoefficientField::R, c.t);
                check_field(2, cell, PgCoefficientField::R, c.f);
                for (std::uint8_t target = 0; target < 3; ++target) {
                    check_field(target, cell, PgCoefficientField::P, zero);
                    check_field(target, cell, PgCoefficientField::C, zero);
                }
                const auto& den = data.denominator.values()[index];
                const std::array normal{op.s(op.m(den, c.fx), minus), op.s(op.m(den, c.fy), minus),
                                        op.k(den)};
                for (std::uint8_t j = 0; j < 3; ++j) {
                    check_field(static_cast<std::uint8_t>(27 + j), cell, PgCoefficientField::P,
                                normal[j]);
                    check_field(static_cast<std::uint8_t>(27 + j), cell, PgCoefficientField::R,
                                zero);
                    check_field(static_cast<std::uint8_t>(27 + j), cell, PgCoefficientField::C,
                                zero);
                    check_field(static_cast<std::uint8_t>(30 + j), cell, PgCoefficientField::R,
                                unit);
                    check_field(static_cast<std::uint8_t>(30 + j), cell, PgCoefficientField::P,
                                zero);
                    check_field(static_cast<std::uint8_t>(30 + j), cell, PgCoefficientField::C,
                                zero);
                }
            }
    }
    /// Construct a separate all-domain curvature and carry its moments across
    /// every interval. Reconstruction coordinates are not used in this path.
    void verify_curvatures(const PgGraphView& graph, const ExactRational& inverse_width,
                           const ExactRational& chart_scale) {
        auto& g = generated;
        auto& data = g.data;
        auto& op = g.op;
        const auto n = data.n, axis = 2 * n;
        const auto legs = pg_legendre(op, ExactAxis::First);
        const auto h = op.number(1, 4 * n), lower = op.number(-1), upper = op.number(1);
        const PgPoly zero;
        for (core::Index degree = 0; degree < 8; ++degree) {
            std::array<ExactRational, 5> carry;
            for (core::Index cs = 0; cs < axis; ++cs) {
                const auto curvature =
                    pg_global_curvature(op, legs[static_cast<std::size_t>(degree)], n, cs);
                const auto& along = g.cells.values()[static_cast<std::size_t>(cs * axis)];
                std::array<PgPoly, 5> moments;
                for (std::size_t m = 0; m < 5; ++m) {
                    const auto active =
                        op.i(op.s(op.m(along.s_power[m], curvature), h), ExactAxis::First, lower);
                    moments[m] = op.a(active, op.constant(carry[m]));
                    carry[m] = op.e(moments[m], upper, op.zero);
                }
                for (core::Index ct = 0; ct < axis; ++ct) {
                    const PgCellIndex cell{cs, ct};
                    const auto& c = g.cells.values()[static_cast<std::size_t>(cs * axis + ct)];
                    const auto v = pg_lifting(op, graph.height, c.f, c.s, c.t_power, moments,
                                              PgPairedAxis::S, inverse_width, graph.skew);
                    const auto p = pg_director(op, v, c.fx, c.fy, c.U, c.V, c.xy, chart_scale,
                                               inverse_width, graph.skew);
                    for (std::uint8_t j = 0; j < 3; ++j) {
                        const auto target = static_cast<std::uint8_t>(3 + j * 8 + degree);
                        check_field(target, cell, PgCoefficientField::R, v[j]);
                        check_field(target, cell, PgCoefficientField::P, p[j]);
                        check_field(target, cell, PgCoefficientField::C, zero);
                    }
                }
            }
        }
    }
};
} // namespace uvlas::structure::detail
