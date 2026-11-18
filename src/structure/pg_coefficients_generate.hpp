#pragma once

#include "pg_coefficients_basis.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>

namespace uvlas::structure::detail {
struct PgCoefficientCellCache {
    PgPoly s, t, f, fx, fy, U, V, xy;
    std::array<PgPoly, 5> s_power, t_power;
};
struct PgCoefficientGeneration {
    PgCoefficientGeneration(PgCoefficientData& target, PgCoefficientAlgebra& algebra,
                            const PgCoefficientLimits& caps)
        : data(target), op(algebra), limits(caps) {}
    PgCoefficientData& data;
    PgCoefficientAlgebra& op;
    const PgCoefficientLimits& limits;
    core::PlannedBuffer<PgCoefficientCellCache> cells;
    std::uint64_t retained_terms = 0;
    std::size_t sorted_record_count = 0;
    [[noreturn]] static void exceeded() {
        throw PgCoefficientAbort{
            {core::FailureCode::ResourceExceeded, core::Operation::PrepareModel}};
    }
    static core::Index source_count(core::Index n) {
        return (n + 2) * (n + 2) + 1 + (2 * n + 1) * (2 * n + 1) + 3 + 96 * n + (n + 1) * (n + 1);
    }
    void source(PgCoefficientKind kind, core::Index first = 0, core::Index second = 0,
                std::uint8_t component = 0, PgPairedAxis axis = PgPairedAxis::S) {
        if (data.column_count >= data.columns.size())
            exceeded();
        const auto id = static_cast<core::Index>(data.column_count);
        data.columns.values()[data.column_count++] = {kind, first, second, component, axis, {id}};
    }
    void original_sources() {
        const auto n = data.n;
        for (core::Index i = 0; i < n + 2; ++i)
            for (core::Index j = 0; j < n + 2; ++j)
                source(PgCoefficientKind::BasicR, i, j);
        source(PgCoefficientKind::ReferenceR);
        for (core::Index i = 0; i < 2 * n + 1; ++i)
            for (core::Index j = 0; j < 2 * n + 1; ++j)
                source(PgCoefficientKind::BasicD, i, j);
        for (std::uint8_t j = 0; j < 3; ++j)
            source(PgCoefficientKind::NormalD, 0, 0, j);
        for (std::uint8_t j = 0; j < 3; ++j)
            for (const auto axis : {PgPairedAxis::S, PgPairedAxis::T})
                for (core::Index origin = 0; origin < 2 * n; ++origin)
                    for (core::Index degree = 0; degree < 8; ++degree)
                        source(PgCoefficientKind::Paired, origin, degree, j, axis);
        for (core::Index i = 0; i < n + 1; ++i)
            for (core::Index j = 0; j < n + 1; ++j)
                source(PgCoefficientKind::BasicC, i, j);
        data.original = static_cast<core::Index>(data.column_count);
    }
    void record(core::Index column, PgCellIndex cell, PgCoefficientField field, PgPoly value) {
        op.check();
        if (value.terms().empty())
            return;
        if (data.record_count >= data.records.size())
            exceeded();
        const auto size = static_cast<std::uint64_t>(value.terms().size());
        if (retained_terms > limits.maximum_stored_terms ||
            size > limits.maximum_stored_terms - retained_terms)
            exceeded();
        retained_terms += size;
        data.records.values()[data.record_count++] = {{column}, cell, field, std::move(value)};
    }
    [[nodiscard]] std::uint64_t key(const PgCoefficientRecord& entry) const noexcept {
        return (static_cast<std::uint64_t>(entry.column.value) *
                    static_cast<std::uint64_t>(4 * data.n * data.n) +
                static_cast<std::uint64_t>(entry.cell.first * (2 * data.n) + entry.cell.second)) *
                   3 +
               static_cast<std::uint64_t>(entry.field);
    }
    void sort_records() {
        auto all = data.records.values().first(data.record_count);
        std::sort(all.begin(), all.end(),
                  [&](const auto& a, const auto& b) { return key(a) < key(b); });
        sorted_record_count = data.record_count;
    }
    [[nodiscard]] const PgPoly& get(core::Index column, PgCellIndex cell,
                                    PgCoefficientField field) const noexcept {
        const auto target =
            (static_cast<std::uint64_t>(column) * static_cast<std::uint64_t>(4 * data.n * data.n) +
             static_cast<std::uint64_t>(cell.first * (2 * data.n) + cell.second)) *
                3 +
            static_cast<std::uint64_t>(field);
        const auto all = data.records.values().first(sorted_record_count);
        const auto found = std::lower_bound(
            all.begin(), all.end(), target,
            [&](const auto& entry, std::uint64_t value) { return key(entry) < value; });
        if (found != all.end() && key(*found) == target)
            return found->polynomial;
        return data.zero;
    }
    void geometry(const PgGraphView& graph, const ExactRational& inverse_width,
                  const ExactRational& chart_scale) {
        const auto axis = 2 * data.n;
        cells = pg_take(
            op.plan.allocate<PgCoefficientCellCache>(axis * axis, core::Operation::PrepareModel));
        const auto h = op.number(1, 4 * data.n);
        const auto unit = op.constant(1);
        const auto sx = op.rm(chart_scale, inverse_width);
        for (core::Index cs = 0; cs < axis; ++cs)
            for (core::Index ct = 0; ct < axis; ++ct) {
                const auto index = static_cast<std::size_t>(cs * axis + ct);
                auto& c = cells.values()[index];
                c.s = op.affine(op.number(2 * cs + 1, 4 * data.n), h, ExactAxis::First);
                c.t = op.affine(op.number(2 * ct + 1, 4 * data.n), h, ExactAxis::Second);
                c.f = op.finish(graph.height.substitute(c.s, c.t, op.limits, op.plan), 9);
                c.fx = op.s(op.d(c.f, ExactAxis::First), sx);
                c.fy =
                    op.sub(op.s(op.d(c.f, ExactAxis::Second), chart_scale), op.s(c.fx, graph.skew));
                const auto xx = op.m(c.fx, c.fx), yy = op.m(c.fy, c.fy);
                c.U = op.a(unit, xx);
                c.V = op.a(unit, yy);
                c.xy = op.m(c.fx, c.fy);
                data.denominator.values()[index] = op.a(c.U, yy);
                c.s_power = pg_powers(op, c.s);
                c.t_power = pg_powers(op, c.t);
            }
    }
    void basic() {
        const auto n = data.n, axis = 2 * n;
        const auto nr = (n + 2) * (n + 2) + 1, nd = (2 * n + 1) * (2 * n + 1) + 3;
        const auto c_start = nr + nd + 96 * n;
        const auto minus = op.number(-1);
        struct BasisSet {
            PgUnivariateRow r, d, c;
        };
        std::array<BasisSet, 8> first, second;
        for (core::Index cell = 0; cell < axis; ++cell) {
            auto& a = first[static_cast<std::size_t>(cell)];
            auto& b = second[static_cast<std::size_t>(cell)];
            a.r = pg_basis(op, n, cell, ExactAxis::First, PgUnivariateKind::QuadraticC1);
            a.d = pg_basis(op, n, cell, ExactAxis::First, PgUnivariateKind::QuadraticC0);
            a.c = pg_basis(op, n, cell, ExactAxis::First, PgUnivariateKind::LinearC0);
            b.r = pg_basis(op, n, cell, ExactAxis::Second, PgUnivariateKind::QuadraticC1);
            b.d = pg_basis(op, n, cell, ExactAxis::Second, PgUnivariateKind::QuadraticC0);
            b.c = pg_basis(op, n, cell, ExactAxis::Second, PgUnivariateKind::LinearC0);
        }
        for (core::Index cs = 0; cs < axis; ++cs) {
            const auto& row = first[static_cast<std::size_t>(cs)];
            const auto& rs = row.r;
            const auto& ds = row.d;
            const auto& hs = row.c;
            for (core::Index ct = 0; ct < axis; ++ct) {
                const PgCellIndex cell{cs, ct};
                const auto index = static_cast<std::size_t>(cs * axis + ct);
                const auto& cache = cells.values()[index];
                const auto& transverse = second[static_cast<std::size_t>(ct)];
                const auto& rt = transverse.r;
                const auto& dt = transverse.d;
                const auto& ht = transverse.c;
                const auto& g = data.denominator.values()[index];
                for (core::Index i = 0; i < rs.count; ++i)
                    for (core::Index j = 0; j < rt.count; ++j)
                        if (!rs.values[static_cast<std::size_t>(i)].terms().empty() &&
                            !rt.values[static_cast<std::size_t>(j)].terms().empty())
                            record(i * (n + 2) + j, cell, PgCoefficientField::R,
                                   op.m(rs.values[static_cast<std::size_t>(i)],
                                        rt.values[static_cast<std::size_t>(j)]));
                record(nr - 1, cell, PgCoefficientField::R, op.k(cache.f));
                for (core::Index i = 0; i < ds.count; ++i)
                    for (core::Index j = 0; j < dt.count; ++j)
                        if (!ds.values[static_cast<std::size_t>(i)].terms().empty() &&
                            !dt.values[static_cast<std::size_t>(j)].terms().empty())
                            record(nr + i * (2 * n + 1) + j, cell, PgCoefficientField::P,
                                   op.m(g, op.m(ds.values[static_cast<std::size_t>(i)],
                                                dt.values[static_cast<std::size_t>(j)])));
                record(nr + nd - 3, cell, PgCoefficientField::P, op.s(op.m(g, cache.fx), minus));
                record(nr + nd - 2, cell, PgCoefficientField::P, op.s(op.m(g, cache.fy), minus));
                record(nr + nd - 1, cell, PgCoefficientField::P, op.k(g));
                for (core::Index i = 0; i < hs.count; ++i)
                    for (core::Index j = 0; j < ht.count; ++j)
                        if (!hs.values[static_cast<std::size_t>(i)].terms().empty() &&
                            !ht.values[static_cast<std::size_t>(j)].terms().empty())
                            record(c_start + i * (n + 1) + j, cell, PgCoefficientField::C,
                                   op.m(hs.values[static_cast<std::size_t>(i)],
                                        ht.values[static_cast<std::size_t>(j)]));
            }
        }
    }
    void paired(const PgGraphView& graph, const ExactRational& inverse_width,
                const ExactRational& chart_scale) {
        const auto n = data.n, axis = 2 * n;
        const auto nr = (n + 2) * (n + 2) + 1, nd = (2 * n + 1) * (2 * n + 1) + 3;
        for (const auto paired_axis : {PgPairedAxis::S, PgPairedAxis::T}) {
            const auto derivative_axis =
                paired_axis == PgPairedAxis::S ? ExactAxis::First : ExactAxis::Second;
            const auto legs = pg_legendre(op, derivative_axis);
            for (core::Index origin = 0; origin < axis; ++origin)
                for (core::Index degree = 0; degree < 8; ++degree) {
                    const auto moments = pg_moments(
                        op, n, origin, legs[static_cast<std::size_t>(degree)], derivative_axis);
                    for (core::Index along = origin; along < axis; ++along) {
                        const auto& selected = along == origin ? moments.active : moments.completed;
                        for (core::Index transverse = 0; transverse < axis; ++transverse) {
                            const PgCellIndex cell = paired_axis == PgPairedAxis::S
                                                         ? PgCellIndex{along, transverse}
                                                         : PgCellIndex{transverse, along};
                            const auto index =
                                static_cast<std::size_t>(cell.first * axis + cell.second);
                            const auto& c = cells.values()[index];
                            const auto& coordinate = paired_axis == PgPairedAxis::S ? c.s : c.t;
                            const auto& powers =
                                paired_axis == PgPairedAxis::S ? c.t_power : c.s_power;
                            auto v = pg_lifting(op, graph.height, c.f, coordinate, powers, selected,
                                                paired_axis, inverse_width, graph.skew);
                            auto p = pg_director(op, v, c.fx, c.fy, c.U, c.V, c.xy, chart_scale,
                                                 inverse_width, graph.skew);
                            for (core::Index j = 0; j < 3; ++j) {
                                const auto column = nr + nd + j * (32 * n) +
                                                    (paired_axis == PgPairedAxis::S ? 0 : 16 * n) +
                                                    origin * 8 + degree;
                                record(column, cell, PgCoefficientField::R,
                                       std::move(v[static_cast<std::size_t>(j)]));
                                record(column, cell, PgCoefficientField::P,
                                       std::move(p[static_cast<std::size_t>(j)]));
                            }
                        }
                    }
                }
        }
        sort_records();
    }
};
} // namespace uvlas::structure::detail
