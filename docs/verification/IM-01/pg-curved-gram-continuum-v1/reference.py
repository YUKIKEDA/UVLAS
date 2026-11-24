"""Continuum moments for the cut parabola in the product-volume regression."""

from collections.abc import Sequence
from decimal import Decimal as D
from decimal import localcontext
from fractions import Fraction as F
from functools import cache
from math import cos, pi

# Each tuple is the surviving term after integrating over the symmetric
# normalized depth. The last two entries are the exponent of g and the multiplier.
FIELD_TERMS = (
    ("RR", 0, F(1, 2), F(1)),
    ("RP", 2, F(-5, 2), F(-1, 8)),
    ("RC", 2, F(1, 2), F(1, 2)),
    ("PP", 2, F(-5, 2), F(1)),
    ("PC", 4, F(-5, 2), F(-1, 16)),
    ("CC", 4, F(1, 2), F(1, 4)),
)
SHAPE_POWERS = ((0, 0), (3, 5), (7, 6), (14, 4), (9, 9), (4, 14))
STRESS_POWERS = ((0, 0), (3, 5), (7, 6), (14, 14))
DEPTH_HALF_WIDTH = (F(1, 4), F(1, 8))
SERIES_TERMS = 56


def decimal(value: F) -> D:
    return D(value.numerator) / D(value.denominator)


def multiply(left: Sequence[F], right: Sequence[F]) -> list[F]:
    product = [F(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            product[i + j] += a * b
    return product


def polynomial_power(base: Sequence[F], exponent: int) -> list[F]:
    product = [F(1)]
    for _ in range(exponent):
        product = multiply(product, base)
    return product


def strip_integral(t_cell: int, t_power: int) -> list[F]:
    """Integrate the local t monomial and return a polynomial in global s."""
    if t_cell == 0:
        lower, upper = (F(0), F(1, 4)), (F(1, 2), F(0))
    else:
        lower, upper = (F(1, 2), F(0)), (F(1), F(-1, 8))
    offset = 2 * t_cell + 1
    local_lower = (4 * lower[0] - offset, 4 * lower[1])
    local_upper = (4 * upper[0] - offset, 4 * upper[1])
    exponent = t_power + 1
    values_at_lower = polynomial_power(local_lower, exponent)
    values_at_upper = polynomial_power(local_upper, exponent)
    return [
        (upper_value - lower_value) / (4 * exponent)
        for lower_value, upper_value in zip(values_at_lower, values_at_upper, strict=True)
    ]


def integrate_polynomial(polynomial: Sequence[F], lower: F, upper: F, shift: int) -> F:
    return sum(
        (
            value
            * (upper ** (degree + shift + 1) - lower ** (degree + shift + 1))
            / (degree + shift + 1)
            for degree, value in enumerate(polynomial)
        ),
        start=F(0),
    )


def analytic(
    cell: tuple[int, int], powers: tuple[int, int], depth_power: int, exponent: F, scale: F
) -> tuple[F, F]:
    s_cell, t_cell = cell
    s_power, t_power = powers
    strip = strip_integral(t_cell, t_power)
    s_monomial = polynomial_power((F(-(2 * s_cell + 1)), F(4)), s_power)
    depth = [
        value * scale * 2 / (depth_power + 1)
        for value in polynomial_power(DEPTH_HALF_WIDTH, depth_power + 1)
    ]
    polynomial = multiply(multiply(strip, s_monomial), depth)
    lower, upper = F(s_cell, 2), F(s_cell + 1, 2)
    coefficient, central = F(1), F(0)
    for term in range(SERIES_TERMS):
        central += coefficient * integrate_polynomial(polynomial, lower, upper, 2 * term)
        coefficient *= (exponent - term) / (64 * (term + 1))

    # On [0,1] the omitted monomials are bounded by the coefficient one-norm.
    ratio = (1 + abs(exponent) / SERIES_TERMS) / 64
    radius = abs(coefficient) * sum(map(abs, polynomial)) / (1 - ratio)
    return central, radius


@cache
def unit_gauss(order: int) -> tuple[tuple[D, D], ...]:
    def evaluate(z: D) -> tuple[D, D]:
        previous, current = D(1), z
        for degree in range(2, order + 1):
            previous, current = (
                current,
                ((2 * degree - 1) * z * current - (degree - 1) * previous) / degree,
            )
        derivative = order * (z * current - previous) / (z * z - 1)
        return current, derivative

    nodes = []
    for root in range(1, order + 1):
        z = D(cos(pi * (root - 0.25) / (order + 0.5)))
        for _ in range(11):
            value, derivative = evaluate(z)
            z -= value / derivative
        _, derivative = evaluate(z)
        nodes.append(((z + 1) / 2, 1 / ((1 - z * z) * derivative * derivative)))
    return tuple(nodes)


def chart_strip(s: D, t_cell: int) -> tuple[D, D]:
    if t_cell == 0:
        return s / 4, D("0.5")
    return D("0.5"), 1 - s / 8


def decimal_power(value: D, exponent: int) -> D:
    # Decimal treats 0**0 as invalid, but the constant monomial is one.
    return value**exponent if exponent else D(1)


def sweep(
    cell: tuple[int, int], powers: tuple[int, int], depth_power: int, exponent: F, scale: F
) -> D:
    """Cross-check the series by integrating the analytic vertical sections."""
    s_cell, t_cell = cell
    s_power, t_power = powers
    total = D(0)
    for node, weight in unit_gauss(48):
        s = (D(s_cell) + node) / 2
        lower, upper = chart_strip(s, t_cell)
        t_offset = 2 * t_cell + 1
        t_exponent = t_power + 1
        strip = ((4 * upper - t_offset) ** t_exponent - (4 * lower - t_offset) ** t_exponent) / (
            4 * t_exponent
        )
        g = 1 + s * s / 64
        root = g.sqrt()
        metric = root if exponent == F(1, 2) else 1 / (root * g * g)
        depth = (D(1) / 4 + s / 8) ** (depth_power + 1)
        s_monomial = decimal_power(4 * s - (2 * s_cell + 1), s_power)
        total += weight * strip * s_monomial * metric * depth * decimal(scale) / (depth_power + 1)
    return total


def product_moments(
    s_cell: int, t_cell: int, order: int, powers: Sequence[tuple[int, int]]
) -> dict[tuple[str, int, int], D]:
    start, end = D(s_cell) / 2, D(s_cell + 1) / 2
    lower_start, upper_start = chart_strip(start, t_cell)
    lower_end, upper_end = chart_strip(end, t_cell)
    a, b = (start, lower_start), (end, lower_end)
    c, d = (end, upper_end), (start, upper_start)
    integrals = {(field, i, j): D(0) for field, *_ in FIELD_TERMS for i, j in powers}

    for origin, first, second in ((a, b, c), (a, c, d)):
        us, ut = first[0] - origin[0], first[1] - origin[1]
        vs, vt = second[0] - origin[0], second[1] - origin[1]
        determinant = us * vt - ut * vs
        for u, u_weight in unit_gauss(order):
            for raw_v, v_weight in unit_gauss(order):
                v = (1 - u) * raw_v
                s = origin[0] + u * us + v * vs
                t = origin[1] + u * ut + v * vt
                local_s, local_t = 4 * s - (2 * s_cell + 1), 4 * t - (2 * t_cell + 1)
                g = 1 + s * s / 64
                root = g.sqrt()
                h = D(1) / 4 + s / 8
                surface_weight = determinant * (1 - u) * u_weight * v_weight
                for depth, depth_weight in unit_gauss(3):
                    z = (2 * depth - 1) * h
                    jacobian = root - z / (8 * g)
                    weight = surface_weight * (2 * h) * depth_weight * jacobian
                    director, quadratic = z / (root * g), z * z / 2
                    terms = {
                        "RR": D(1),
                        "RP": director,
                        "RC": quadratic,
                        "PP": director * director,
                        "PC": director * quadratic,
                        "CC": quadratic * quadratic,
                    }
                    for i, j in powers:
                        monomial = decimal_power(local_s, i) * decimal_power(local_t, j)
                        for field, term in terms.items():
                            integrals[(field, i, j)] += weight * term * monomial
    return integrals


def compare(powers: Sequence[tuple[int, int]], orders: Sequence[int]) -> None:
    reference = {}
    largest_tail, largest_sweep_error = D(0), D(0)
    for field, depth_power, exponent, scale in FIELD_TERMS:
        for s_cell in range(2):
            for t_cell in range(2):
                for i, j in powers:
                    exact, tail = analytic((s_cell, t_cell), (i, j), depth_power, exponent, scale)
                    cross_check = sweep((s_cell, t_cell), (i, j), depth_power, exponent, scale)
                    largest_tail = max(largest_tail, decimal(tail))
                    largest_sweep_error = max(
                        largest_sweep_error, abs(cross_check - decimal(exact))
                    )
                    reference[(s_cell, t_cell, field, i, j)] = exact, tail
    if largest_tail >= D("1e-50") or largest_sweep_error >= D("1e-40"):
        raise ValueError("Independent continuum reference exceeded its precision budget")

    print(
        f"moments={len(reference)} exact_tail={largest_tail:.8E} "
        f"independent_sweep_gap={largest_sweep_error:.8E}"
    )
    for order in orders:
        worst = D(0)
        controlling_term = None
        for s_cell in range(2):
            for t_cell in range(2):
                samples = product_moments(s_cell, t_cell, order, powers)
                for (field, i, j), value in samples.items():
                    key = (s_cell, t_cell, field, i, j)
                    exact, tail = reference[key]
                    gap = abs(value - decimal(exact)) + decimal(tail)
                    if gap > worst:
                        worst = gap
                        controlling_term = key
        print(
            f"order={order} points={8 * order**2 * 3} largest_gap={worst:.8E} "
            f"controlling_term={controlling_term} threshold_1e-12={worst < D('1e-12')}"
        )


def main() -> None:
    with localcontext() as context:
        context.prec = 85
        print("Native half-cells, products up to summed powers 18")
        compare(SHAPE_POWERS, (5, 7, 9, 11, 13, 15))
        print("Higher-degree stress comparison, powers (14, 14)")
        compare(STRESS_POWERS, (9, 12, 14, 16, 18, 20))


if __name__ == "__main__":
    main()
