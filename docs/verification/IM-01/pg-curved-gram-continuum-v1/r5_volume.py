"""Bound the full curved swept-taper R5 volume using exact polynomial moments."""

from decimal import Decimal as D
from decimal import localcontext
from fractions import Fraction as Q
from functools import cache
from math import comb

from reference import decimal, unit_gauss

Poly = dict[tuple[int, int], int]
# Scale factors: slopes by 100, second derivatives by 100, thickness shape by 20.
SLOPE_S = {(1, 0): 40, (2, 0): 10, (1, 2): 5}
SLOPE_T = {(0, 1): 20, (0, 2): -6, (2, 1): 5}
SECOND_SS = {(0, 0): 40, (1, 0): 20, (0, 2): 5}
SECOND_TT = {(0, 0): 20, (0, 1): -12, (2, 0): 5}
SECOND_ST = {(1, 1): 10}
THICKNESS = {(0, 0): 20, (1, 0): 2, (0, 1): 1}
AREA = Q(4, 5)
SLOPE_SCALE_SQUARED = 10_000
SERIES_TERMS = 36


def add(left: Poly, right: Poly, sign: int = 1) -> Poly:
    result = left.copy()
    for key, value in right.items():
        result[key] = result.get(key, 0) + sign * value
    return {key: value for key, value in result.items() if value}


def multiply(left: Poly, right: Poly) -> Poly:
    result: Poly = {}
    for (i, j), a in left.items():
        for (k, m), b in right.items():
            key = (i + k, j + m)
            result[key] = result.get(key, 0) + a * b
    return {key: value for key, value in result.items() if value}


@cache
def moment(a: int, b: int) -> Q:
    # Integrate s^a t^b over 0.3s <= t <= 1 - 0.1s in the unit chart.
    degree = b + 1
    top = sum(
        (Q((-1) ** k * comb(degree, k), 10**k * (a + k + 1)) for k in range(degree + 1)),
        start=Q(0),
    )
    bottom = Q(3, 10) ** degree / (a + degree + 1)
    return (top - bottom) / degree


def series_volume_term(polynomial: Poly, perturbation: Poly, exponent: Q) -> tuple[Q, Q]:
    base = polynomial
    value, coefficient = Q(0), Q(1)
    for order in range(SERIES_TERMS):
        value += coefficient * sum((c * moment(*key) for key, c in base.items()), start=Q(0))
        coefficient *= (exponent - order) / ((order + 1) * SLOPE_SCALE_SQUARED)
        if order + 1 < SERIES_TERMS:
            base = multiply(base, perturbation)

    slope_bound = sum(map(abs, SLOPE_S.values())) ** 2 + sum(map(abs, SLOPE_T.values())) ** 2
    ratio = (1 + abs(exponent) / SERIES_TERMS) * Q(slope_bound, SLOPE_SCALE_SQUARED)
    remainder = abs(coefficient) * slope_bound**SERIES_TERMS / (1 - ratio)
    return value, AREA * sum(map(abs, polynomial.values())) * remainder


def pointwise_volume(h: D, order: int) -> D:
    # Independent integration of the exact depth antiderivative along vertical chart strips.
    result = D(0)
    for s, ws in unit_gauss(order):
        lower = 3 * s / 10
        span = 1 - 4 * s / 10
        for v, wv in unit_gauss(order):
            t = lower + span * v
            fs = (40 * s + 10 * s * s + 5 * s * t * t) / 100
            ft = (20 * t - 6 * t * t + 5 * s * s * t) / 100
            fss = (40 + 20 * s + 5 * t * t) / 100
            ftt = (20 - 12 * t + 5 * s * s) / 100
            fst = s * t / 10
            gd = 1 + fs * fs + ft * ft
            root = gd.sqrt()
            thickness = h * (1 + s / 10 + t / 20)
            value = thickness * root + thickness**3 * (fss * ftt - fst * fst) / (12 * gd * root)
            result += ws * wv * span * value
    return result


def main() -> None:
    assert moment(0, 0) == AREA
    assert moment(1, 0) == Q(11, 30)
    assert moment(0, 1) == Q(131, 300)
    metric = add(multiply(SLOPE_S, SLOPE_S), multiply(SLOPE_T, SLOPE_T))
    gaussian = add(multiply(SECOND_SS, SECOND_TT), multiply(SECOND_ST, SECOND_ST), -1)
    thickness_cubed = multiply(multiply(THICKNESS, THICKNESS), THICKNESS)
    primary, primary_error = series_volume_term(THICKNESS, metric, Q(1, 2))
    secondary, secondary_error = series_volume_term(
        multiply(thickness_cubed, gaussian), metric, Q(-3, 2)
    )

    with localcontext() as context:
        context.prec = 85
        for h in (Q(1, 100), Q(1, 10000)):
            scale = h**3 / (12 * 20**3 * SLOPE_SCALE_SQUARED)
            centre = primary * h / 20 + secondary * scale
            error = primary_error * h / 20 + secondary_error * scale
            direct = pointwise_volume(decimal(h), 32)
            inverse_degree = pointwise_volume(decimal(h), 28)
            if not decimal(centre - error) <= direct <= decimal(centre + error):
                raise ValueError("Independent chart quadrature fell outside the volume enclosure")
            if abs(direct - inverse_degree) > decimal(error):
                raise ValueError("The two chart quadratures differ beyond the enclosure radius")
            print(
                f"thickness={decimal(h)} volume={direct:.34E} "
                f"exact_series_radius={decimal(error):.8E}",
                flush=True,
            )


if __name__ == "__main__":
    main()
