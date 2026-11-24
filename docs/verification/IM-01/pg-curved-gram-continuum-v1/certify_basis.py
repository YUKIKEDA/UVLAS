"""Check the cut-parabola continuum and finite Gram through the returned exact transform."""

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

from compare_basis import FIELD_FACTORS, RANK, cholesky, field_pairs, forward, read_basis
from finite_moments import MOMENT_SCALE as GRAM_SCALE
from finite_moments import sampled_bounds
from reference import analytic

TRANSFORM_SCALE = 2**216


def floor_ratio(numerator: int, denominator: int) -> int:
    return numerator // denominator


def ceil_ratio(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def on_grid(value: str, scale: int) -> int:
    scaled = Fraction(value) * scale
    if scaled.denominator != 1:
        raise ValueError("Value is off the declared dyadic grid")
    return scaled.numerator


def read_preparation(path: Path):
    transform = [[0] * RANK for _ in range(RANK)]
    finite = []
    recovered = []
    bound = None
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        label, *words = line.split(",")
        if label == "gram":
            lower, upper = (on_grid(word, GRAM_SCALE) for word in words)
            if lower > upper:
                raise ValueError("Reversed finite Gram interval")
            finite.append((lower, upper))
        if label == "transform":
            i, j = map(int, words[:2])
            if not 0 <= i <= j < RANK:
                raise ValueError("Out-of-range transform entry")
            recovered.append((i, j))
            transform[i][j] = on_grid(words[2], TRANSFORM_SCALE)
        if label == "congruence":
            if bound is not None:
                raise ValueError("Repeated finite certificate")
            bound = Fraction(words[0])

    ordered = [(i, j) for j in range(RANK) for i in range(j + 1)]
    if len(finite) != len(ordered) or recovered != ordered or bound is None:
        raise ValueError("Incomplete finite Gram or transform")
    if not 0 <= bound <= Fraction(1, 10**20) or any(transform[i][i] <= 0 for i in range(RANK)):
        raise ValueError("Finite transform is outside the qualification bounds")
    return finite, transform, bound


def continuum_bounds(cell, keys):
    moments = {}
    for field, a, b in keys:
        centre, radius = analytic(cell, (a, b), *FIELD_FACTORS[field])
        low = (centre - radius) * GRAM_SCALE
        high = (centre + radius) * GRAM_SCALE
        moments[(field, a, b)] = (
            floor_ratio(low.numerator, low.denominator),
            ceil_ratio(high.numerator, high.denominator),
        )
    return moments


def gram_enclosure(cell_terms, moment_provider):
    lower = [[0] * RANK for _ in range(RANK)]
    upper = [[0] * RANK for _ in range(RANK)]
    for cell, (pairs, keys) in cell_terms.items():
        moments = moment_provider(cell, keys)
        for (i, j), terms in pairs.items():
            for key, coefficient in terms.items():
                lo, hi = moments[key]
                numerator = coefficient.numerator
                denominator = coefficient.denominator
                if numerator < 0:
                    lo, hi = hi, lo
                lower[i][j] += floor_ratio(lo * numerator, denominator)
                upper[i][j] += ceil_ratio(hi * numerator, denominator)

    return [(lower[i][j], upper[i][j]) for j in range(RANK) for i in range(j + 1)]


def congruence_bounds(
    entries, transform: list[list[int]], subtract_identity: bool
) -> tuple[Fraction, Fraction]:
    size = len(transform)
    midpoints = [[0] * size for _ in range(size)]
    radii = [[0] * size for _ in range(size)]
    positions = ((i, j) for j in range(size) for i in range(j + 1))
    for (i, j), (lo, hi) in zip(positions, entries, strict=True):
        # The doubled midpoint and radius share an additional factor of two.
        midpoints[i][j] = midpoints[j][i] = hi + lo
        radii[i][j] = radii[j][i] = hi - lo

    transposed = list(zip(*transform, strict=True))
    magnitudes = [tuple(map(abs, row)) for row in transposed]
    scale = 2 * GRAM_SCALE * TRANSFORM_SCALE**2

    def dot(left, right) -> int:
        return sum(a * b for a, b in zip(left, right, strict=True))

    row_lower = []
    row_upper = []
    for i in range(size):
        left_midpoints = [dot(transposed[i], row) for row in midpoints]
        left_radii = [dot(magnitudes[i], row) for row in radii]
        lo = 0
        hi = 0
        for j in range(size):
            midpoint = dot(left_midpoints, transposed[j])
            if subtract_identity and i == j:
                midpoint -= scale
            radius = dot(left_radii, magnitudes[j])
            lo += max(abs(midpoint) - radius, 0)
            hi += abs(midpoint) + radius
        row_lower.append(lo)
        row_upper.append(hi)
    return Fraction(max(row_lower), scale), Fraction(max(row_upper), scale)


def propose(finite):
    with localcontext() as context:
        context.prec = 120
        matrix = [[Decimal(0)] * RANK for _ in range(RANK)]
        positions = ((i, j) for j in range(RANK) for i in range(j + 1))
        for (i, j), (lo, hi) in zip(positions, finite, strict=True):
            matrix[i][j] = matrix[j][i] = Decimal(lo + hi) / (2 * GRAM_SCALE)
        identity = [[Decimal(i == j) for j in range(RANK)] for i in range(RANK)]
        inverse = forward(cholesky(matrix), identity)
        return [
            [int((v * TRANSFORM_SCALE).to_integral_value()) for v in row]
            for row in zip(*inverse, strict=True)
        ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("basis_export", type=Path)
    parser.add_argument("--finite-preparation", type=Path)
    parser.add_argument("--order", type=int, choices=(13, 15, 17), default=17)
    args = parser.parse_args()

    columns, triangles = read_basis(args.basis_export)
    cell_terms = {cell: field_pairs(columns, cell) for cell in triangles}
    continuum = gram_enclosure(cell_terms, continuum_bounds)
    finite = gram_enclosure(
        cell_terms, lambda cell, keys: sampled_bounds(cell, triangles[cell], keys, args.order)
    )
    if args.finite_preparation is not None:
        product, transform, certificate = read_preparation(args.finite_preparation)
        if any(
            a_hi < b_lo or b_hi < a_lo
            for (a_lo, a_hi), (b_lo, b_hi) in zip(finite, product, strict=True)
        ):
            raise ValueError("Independent and production finite Gram intervals disagree")
        decimal_certificate = Decimal(certificate.numerator) / Decimal(certificate.denominator)
        print(f"product_finite_certificate={decimal_certificate:.12E}", flush=True)
    else:
        transform = propose(finite)
    difference = [
        (lo - finite_hi, hi - finite_lo)
        for (lo, hi), (finite_lo, finite_hi) in zip(continuum, finite, strict=True)
    ]
    limit = Fraction(1, 10**20)
    for label, (lo, hi) in (
        ("independent_finite_congruence", congruence_bounds(finite, transform, True)),
        ("continuous_difference", congruence_bounds(difference, transform, False)),
        ("continuous_congruence", congruence_bounds(continuum, transform, True)),
    ):
        lower = Decimal(lo.numerator) / Decimal(lo.denominator)
        upper = Decimal(hi.numerator) / Decimal(hi.denominator)
        print(f"{label}=[{lower:.12E}, {upper:.12E}]", flush=True)
        if lo > limit:
            raise ValueError(f"{label} exceeds {limit} on the cut-parabola fixture")
        if hi > limit:
            raise ValueError(f"{label} cannot be qualified at {limit} on the cut-parabola fixture")


if __name__ == "__main__":
    main()
