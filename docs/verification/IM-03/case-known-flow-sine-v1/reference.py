"""Independent rational bounds for the analytic 45 degree known-flow fixture."""

from fractions import Fraction
from math import isqrt


def arctan_bounds(denominator: int) -> tuple[Fraction, Fraction]:
    """The first even number of terms and the next term bound arctan(1/d)."""
    terms = 100
    lower = sum(
        (Fraction((-1) ** i, (2 * i + 1) * denominator ** (2 * i + 1)) for i in range(terms)),
        Fraction(0),
    )
    return lower, lower + Fraction(1, (2 * terms + 1) * denominator ** (2 * terms + 1))


def decimal_test_data() -> list[list[str]]:
    atan3 = arctan_bounds(3)
    atan7 = arctan_bounds(7)
    # tan(2 atan(1/3) + atan(1/7)) = 1 in (0, pi/2).
    pi = tuple(8 * u + 4 * v for u, v in zip(atan3, atan7, strict=True))
    resolution = 10**85
    root = isqrt(2 * resolution**2)
    sine = (Fraction(root, 2 * resolution), Fraction(root + 1, 2 * resolution))
    result = []
    output_resolution = 10**64
    uncertainty = Fraction(1, 2**205)
    for n in (3, 4):
        e = Fraction(n, 5)
        scale = Fraction(7, 2) * e
        velocity = tuple(n - 1 + scale * s for s in sine)
        spatial = tuple(scale * p * s / 2 for p, s in zip(pi, sine, strict=True))
        time = tuple(e * 3 * s - scale * p * s / 2 for p in pi for s in sine)
        row = []
        for values in (velocity, time, spatial):
            lo, hi = min(values), max(values)
            digits = round((lo + hi) * output_resolution / 2)
            reference = Fraction(digits, output_resolution)
            assert max(abs(reference - lo), abs(reference - hi)) <= uncertainty
            row.append(str(digits))
        result.append(row)
    return result


if __name__ == "__main__":
    for row in decimal_test_data():
        for value in row:
            print(value)