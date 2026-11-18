#pragma once

#include "pg_coefficients_generate.hpp"

#include <array>
#include <cstddef>

namespace uvlas::structure::detail {
inline void pg_require_identity(bool valid) {
    if (!valid)
        throw PgCoefficientAbort{
            {core::FailureCode::PrecisionUnresolved, core::Operation::PrepareModel}};
}
struct PgCoefficientCompletion {
    explicit PgCoefficientCompletion(PgCoefficientGeneration& value) : generated(value) {}
    PgCoefficientGeneration& generated;
    std::array<PgPoly, 8> ell;
    void anchor() {
        auto& g = generated;
        auto& data = g.data;
        auto& op = g.op;
        const auto axis = 2 * data.n;
        const auto one = op.constant(1);
        for (core::Index cs = 0; cs < axis; ++cs)
            ell[static_cast<std::size_t>(cs)] =
                op.sub(one, g.cells.values()[static_cast<std::size_t>(cs * axis)].s);
        for (const auto field : {PgCoefficientField::R, PgCoefficientField::P})
            for (core::Index column = 0; column < data.original; ++column) {
                const auto& source = data.columns.values()[static_cast<std::size_t>(column)];
                const bool allowed = field == PgCoefficientField::R
                                         ? (source.kind == PgCoefficientKind::BasicR ||
                                            source.kind == PgCoefficientKind::ReferenceR ||
                                            source.kind == PgCoefficientKind::Paired)
                                         : (source.kind == PgCoefficientKind::BasicD ||
                                            source.kind == PgCoefficientKind::NormalD ||
                                            source.kind == PgCoefficientKind::Paired);
                if (!allowed)
                    continue;
                std::array<PgPoly, 8> traces;
                bool nonzero = false;
                for (core::Index ct = 0; ct < axis; ++ct) {
                    traces[static_cast<std::size_t>(ct)] =
                        pg_anchor(op, g.get(column, {0, ct}, field));
                    nonzero = nonzero || !traces[static_cast<std::size_t>(ct)].terms().empty();
                }
                const bool structural_zero =
                    (source.kind == PgCoefficientKind::Paired && source.axis == PgPairedAxis::S) ||
                    ((source.kind == PgCoefficientKind::BasicR ||
                      source.kind == PgCoefficientKind::BasicD) &&
                     source.first != 0);
                if (structural_zero) {
                    pg_require_identity(!nonzero);
                    continue;
                }
                if (!nonzero)
                    continue;
                if (data.column_count >= data.columns.size())
                    g.exceeded();
                auto origin = source;
                origin.kind = field == PgCoefficientField::R ? PgCoefficientKind::AnchorR
                                                             : PgCoefficientKind::AnchorP;
                origin.original = {column};
                const auto target = static_cast<core::Index>(data.column_count);
                data.columns.values()[data.column_count++] = origin;
                for (core::Index cs = 0; cs < axis; ++cs)
                    for (core::Index ct = 0; ct < axis; ++ct)
                        if (!traces[static_cast<std::size_t>(ct)].terms().empty())
                            g.record(target, {cs, ct}, field,
                                     op.m(ell[static_cast<std::size_t>(cs)],
                                          traces[static_cast<std::size_t>(ct)]));
            }
        g.sort_records();
    }
    void c_map() {
        auto& g = generated;
        auto& data = g.data;
        auto& op = g.op;
        const auto n = data.n;
        const auto start = data.original - (n + 1) * (n + 1);
        for (core::Index it = 0; it < n + 1; ++it)
            for (core::Index is = 0; is < n; ++is) {
                if (data.c_map_count >= data.c_map.size())
                    g.exceeded();
                data.c_map.values()[data.c_map_count++] = {
                    {start + it}, {start + is * (n + 1) + it}, op.number(n - is, n)};
            }
    }
    void verify_anchor() {
        auto& g = generated;
        auto& data = g.data;
        auto& op = g.op;
        const auto axis = 2 * data.n;
        for (core::Index column = data.original;
             column < static_cast<core::Index>(data.column_count); ++column) {
            const auto& source = data.columns.values()[static_cast<std::size_t>(column)];
            const auto field = source.kind == PgCoefficientKind::AnchorR ? PgCoefficientField::R
                                                                         : PgCoefficientField::P;
            for (core::Index ct = 0; ct < axis; ++ct) {
                const auto actual = pg_anchor(op, g.get(column, {0, ct}, field));
                const auto expected = pg_anchor(op, g.get(source.original.value, {0, ct}, field));
                pg_require_identity(pg_equal(actual, expected));
            }
        }
        const auto c_count = (data.n + 1) * (data.n + 1), c_start = data.original - c_count;
        for (core::Index original = c_start; original < data.original; ++original)
            for (core::Index ct = 0; ct < axis; ++ct) {
                const auto trace = pg_anchor(op, g.get(original, {0, ct}, PgCoefficientField::C));
                for (core::Index cs = 0; cs < axis; ++cs) {
                    const auto expected = op.m(ell[static_cast<std::size_t>(cs)], trace);
                    PgPoly actual;
                    for (const auto& entry : data.c_map.values().first(data.c_map_count))
                        if (entry.source.value == original) {
                            const auto& basis =
                                g.get(entry.target.value, {cs, ct}, PgCoefficientField::C);
                            if (!basis.terms().empty())
                                actual = op.a(actual, op.s(basis, entry.coefficient));
                        }
                    pg_require_identity(pg_equal(actual, expected));
                }
            }
    }
};
} // namespace uvlas::structure::detail
