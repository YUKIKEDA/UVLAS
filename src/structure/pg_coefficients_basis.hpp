#pragma once

#include "pg_coefficients_algebra.hpp"

#include <array>
#include <cstddef>
#include <utility>

namespace uvlas::structure::detail {
enum class PgUnivariateKind : std::uint8_t { QuadraticC1, QuadraticC0, LinearC0 };
struct PgUnivariateRow {
    std::array<PgPoly, 9> values;
    core::Index count = 0;
};
/// Cell-interior knot membership chooses the polynomial piece, never an averaged
/// endpoint value. Endpoints are then evaluated from this same one-sided piece.
/// Knots use integer multiples of 1/n, avoiding rounded comparisons entirely.
inline PgUnivariateRow pg_basis(PgCoefficientAlgebra& op, core::Index n, core::Index cell,
                                ExactAxis axis, PgUnivariateKind kind) {
    PgUnivariateRow result;
    if (kind == PgUnivariateKind::LinearC0) {
        result.count = n + 1;
        const auto coarse = cell / 2;
        const auto left_origin = op.number(4 * (coarse + 1) - (2 * cell + 1), 4);
        const auto right_origin = op.number((2 * cell + 1) - 4 * coarse, 4);
        const auto quarter = op.number(1, 4), minus_quarter = op.number(-1, 4);
        result.values[static_cast<std::size_t>(coarse)] =
            op.affine(left_origin, minus_quarter, axis);
        result.values[static_cast<std::size_t>(coarse + 1)] =
            op.affine(right_origin, quarter, axis);
        return result;
    }
    std::array<core::Index, 12> knots{};
    std::size_t size = 3;
    for (core::Index j = 1; j < n; ++j) {
        knots[size++] = j;
        if (kind == PgUnivariateKind::QuadraticC0)
            knots[size++] = j;
    }
    for (std::size_t j = 0; j < 3; ++j)
        knots[size++] = n;
    std::array<PgPoly, 11> current;
    for (std::size_t j = 0; j + 1 < size; ++j)
        if (4 * knots[j] <= 2 * cell + 1 && 2 * cell + 1 < 4 * knots[j + 1])
            current[j] = op.constant(1);
    for (std::size_t degree = 1; degree <= 2; ++degree) {
        std::array<PgPoly, 11> next;
        for (std::size_t j = 0; j + degree + 1 < size; ++j) {
            const auto left_denominator = knots[j + degree] - knots[j];
            const auto right_denominator = knots[j + degree + 1] - knots[j + 1];
            PgPoly left, right;
            if (left_denominator != 0 && !current[j].terms().empty()) {
                const auto origin = op.number(2 * cell + 1 - 4 * knots[j], 4 * left_denominator);
                const auto slope = op.number(1, 4 * left_denominator);
                const auto factor = op.affine(origin, slope, axis);
                left = op.m(factor, current[j]);
            }
            if (right_denominator != 0 && !current[j + 1].terms().empty()) {
                const auto origin =
                    op.number(4 * knots[j + degree + 1] - (2 * cell + 1), 4 * right_denominator);
                const auto slope = op.number(-1, 4 * right_denominator);
                const auto factor = op.affine(origin, slope, axis);
                right = op.m(factor, current[j + 1]);
            }
            if (left.terms().empty())
                next[j] = std::move(right);
            else if (right.terms().empty())
                next[j] = std::move(left);
            else
                next[j] = op.a(left, right);
        }
        current = std::move(next);
    }
    result.count = static_cast<core::Index>(size - 3);
    for (std::size_t j = 0; j < static_cast<std::size_t>(result.count); ++j)
        result.values[j] = std::move(current[j]);
    return result;
}
} // namespace uvlas::structure::detail
