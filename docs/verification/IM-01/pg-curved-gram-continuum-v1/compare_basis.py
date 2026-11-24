"""Measure the curved Gram convergence on the exported product basis."""

import argparse
from collections import defaultdict
from collections.abc import Sequence
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from pathlib import Path

from reference import FIELD_TERMS, analytic, decimal, unit_gauss

RANK = 130
FIELD_PAIRS = (("RR", "RP", "RC"), ("RP", "PP", "PC"), ("RC", "PC", "CC"))
FIELD_FACTORS = {name: (depth, exponent, scale) for name, depth, exponent, scale in FIELD_TERMS}
Matrix = list[list[D]]


def read_basis(path: Path):
    columns = [defaultdict(list) for _ in range(RANK)]
    volumes = {}
    term_keys = set()
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        words = line.split(",")
        if words[0] == "volume":
            _, ordinal, s_cell, t_cell, *coordinates = words
            index = int(ordinal)
            if index in volumes or len(coordinates) != 8:
                raise ValueError("Repeated or incomplete volume record")
            cell = int(s_cell), int(t_cell)
            x0, y0, x1, y1, x2, y2, lower, upper = map(F, coordinates)
            volumes[index] = cell, (x0, y0, x1, y1, x2, y2), (lower, upper)
        elif len(words) == 7 and words[0].isdigit():
            column, s_cell, t_cell, field, s_power, t_power = map(int, words[:6])
            key = tuple(words[:6])
            if key in term_keys or not (0 <= column < RANK and 0 <= field < len(FIELD_PAIRS)):
                raise ValueError("Unexpected basis term")
            term_keys.add(key)
            columns[column][(s_cell, t_cell)].append((field, s_power, t_power, F(words[6])))

    if set(volumes) != set(range(40)) or len(term_keys) != 4008 or not all(columns):
        raise ValueError("Basis export lacks expected columns or partition records")
    layers_by_triangle = defaultdict(set)
    for cell, vertices, layer in volumes.values():
        layers_by_triangle[(cell, vertices)].add(layer)
    expected_layers = {(F(-1), F(-1, 3)), (F(-1, 3), F(1))}
    if len(layers_by_triangle) != 20 or any(
        layers != expected_layers for layers in layers_by_triangle.values()
    ):
        raise ValueError("Layer coverage in the export is not the tested full thickness")

    triangles = defaultdict(list)
    for cell, vertices in sorted(layers_by_triangle):
        triangles[cell].append(vertices)
    if set(triangles) != {(0, 0), (0, 1), (1, 0), (1, 1)}:
        raise ValueError("Inconsistent native cells")
    return columns, triangles


def field_pairs(columns, cell):
    rows = [
        (i, [(field, a, b, decimal(coefficient)) for field, a, b, coefficient in parts[cell]])
        for i, parts in enumerate(columns)
        if cell in parts
    ]
    pairs, keys = {}, set()
    for row, (i, left) in enumerate(rows):
        for j, right in rows[row:]:
            combined = defaultdict(D)
            for field_i, a_i, b_i, value_i in left:
                for field_j, a_j, b_j, value_j in right:
                    field = FIELD_PAIRS[field_i][field_j]
                    combined[(field, a_i + a_j, b_i + b_j)] += value_i * value_j
            active = {key: coefficient for key, coefficient in combined.items() if coefficient}
            pairs[(i, j)] = active
            keys.update(active)
    return pairs, sorted(keys)


def continuous_moments(cell, keys):
    return {
        (field, a, b): decimal(analytic(*cell, a, b, *FIELD_FACTORS[field])[0])
        for field, a, b in keys
    }


def powers(base: D, maximum: int) -> list[D]:
    values = [D(1)]
    for _ in range(maximum):
        values.append(values[-1] * base)
    return values


def sampled_moments(cell, triangles, keys, order):
    # Three exact Gauss points integrate a degree-five depth polynomial.
    # Apply its exact integral across the two verified layers so that the
    # measured error isolates the in-plane quadrature and not depth rounding.
    moments = dict.fromkeys(keys, D(0))
    s_cell, t_cell = cell
    max_s = max(a for _, a, _ in keys)
    max_t = max(b for _, _, b in keys)
    for triangle in triangles:
        s0, t0, s1, t1, s2, t2 = map(decimal, triangle)
        ds_u, dt_u = s1 - s0, t1 - t0
        ds_v, dt_v = s2 - s0, t2 - t0
        determinant = ds_u * dt_v - dt_u * ds_v
        if determinant <= 0:
            raise ValueError("Nonpositive triangle orientation")
        for u, u_weight in unit_gauss(order):
            for raw_v, v_weight in unit_gauss(order):
                v = (1 - u) * raw_v
                s = s0 + u * ds_u + v * ds_v
                t = t0 + u * dt_u + v * dt_v
                s_values = powers(4 * s - (2 * s_cell + 1), max_s)
                t_values = powers(4 * t - (2 * t_cell + 1), max_t)
                g = 1 + s * s / 64
                root = g.sqrt()
                z_half = D(1) / 4 + s / 8
                weight = determinant * (1 - u) * u_weight * v_weight
                fields = {}
                for name, depth_power, exponent, scale in FIELD_TERMS:
                    metric = root if exponent == F(1, 2) else 1 / (root * g * g)
                    thickness = 2 * z_half ** (depth_power + 1) / (depth_power + 1)
                    fields[name] = weight * metric * thickness * decimal(scale)
                for field, a, b in keys:
                    moments[(field, a, b)] += fields[field] * s_values[a] * t_values[b]
    return moments


def matrix(cell_terms, triangles, order: int | None = None) -> Matrix:
    result = [[D(0)] * RANK for _ in range(RANK)]
    for cell, (pairs, keys) in cell_terms.items():
        if order is None:
            moments = continuous_moments(cell, keys)
        else:
            moments = sampled_moments(cell, triangles[cell], keys, order)
        for (i, j), terms in pairs.items():
            result[i][j] += sum(coefficient * moments[key] for key, coefficient in terms.items())
    for j in range(RANK):
        for i in range(j):
            result[j][i] = result[i][j]
    return result


def cholesky(gram: Matrix) -> Matrix:
    lower = [[D(0)] * RANK for _ in range(RANK)]
    for i in range(RANK):
        for j in range(i + 1):
            remainder = gram[i][j] - sum(lower[i][k] * lower[j][k] for k in range(j))
            lower[i][j] = remainder.sqrt() if i == j else remainder / lower[j][j]
    return lower


def forward(lower: Matrix, rhs: Sequence[Sequence[D]]) -> Matrix:
    result = [[D(0)] * RANK for _ in range(RANK)]
    for i in range(RANK):
        inverse_diagonal = 1 / lower[i][i]
        for j in range(RANK):
            value = rhs[i][j] - sum(lower[i][k] * result[k][j] for k in range(i))
            result[i][j] = value * inverse_diagonal
    return result


def row_norm(values: Sequence[Sequence[D]]) -> D:
    return max(sum(map(abs, row), start=D(0)) for row in values)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("basis_export", type=Path)
    args = parser.parse_args()

    with localcontext() as context:
        context.prec = 85
        columns, triangles = read_basis(args.basis_export)
        cell_terms = {cell: field_pairs(columns, cell) for cell in triangles}
        continuous = matrix(cell_terms, triangles)
        factor = cholesky(continuous)
        smallest = min(factor[i][i] for i in range(RANK))
        print(f"Continuous Gram minimum Cholesky diagonal {smallest:.8E}", flush=True)
        for order in (9, 11, 13, 15, 17):
            approximate = matrix(cell_terms, triangles, order)
            delta = [
                [a - b for a, b in zip(x, y, strict=True)]
                for x, y in zip(approximate, continuous, strict=True)
            ]
            left = forward(factor, delta)
            normalized_transpose = forward(factor, list(zip(*left, strict=True)))
            normalized = list(zip(*normalized_transpose, strict=True))
            print(
                f"order={order} raw_row_norm={row_norm(delta):.8E} "
                f"continuum_normalized_row_norm={row_norm(normalized):.8E}",
                flush=True,
            )


if __name__ == "__main__":
    main()
