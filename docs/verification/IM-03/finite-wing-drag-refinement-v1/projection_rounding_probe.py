"""Certify the closest Binary64 projection of four product-geometry diagonals."""

from fractions import Fraction
import math


class Interval:
    def __init__(self, lower, upper=None):
        self.lower = Fraction(lower)
        self.upper = self.lower if upper is None else Fraction(upper)

    def __neg__(self):
        return Interval(-self.upper, -self.lower)

    def __add__(self, other):
        return Interval(self.lower + other.lower, self.upper + other.upper)

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        products = (self.lower * other.lower, self.lower * other.upper,
                    self.upper * other.lower, self.upper * other.upper)
        return Interval(min(products), max(products))

    def __truediv__(self, other):
        if other.lower * other.upper <= 0:
            raise ValueError("Division by an interval containing zero")
        return self * Interval(1 / other.upper, 1 / other.lower)


def root(value):
    scale = 1 << 256
    floor = math.isqrt(value.numerator * scale * scale // value.denominator)
    lower = Fraction(floor, scale)
    return Interval(lower, lower if lower * lower == value else Fraction(floor + 1, scale))


def atan(value):
    partial = sum((-1) ** k * value ** (2 * k + 1) / (2 * k + 1) for k in range(96))
    return Interval(partial, partial + value ** 193 / 193)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b, strict=True))


def difference(a, b):
    return tuple(x - y for x, y in zip(a, b, strict=True))


def cross(a, b):
    return tuple(a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3]
                 for i in range(3))


def coefficient(query, tangents, vertices, radius):
    normal = cross(tangents[:3], tangents[3:])
    four_pi = Interval(64) * atan(Fraction(1, 5)) - Interval(16) * atan(Fraction(1, 239))
    factor = four_pi * root(dot(normal, normal))
    radius2 = radius * radius
    value = Interval(0)
    for a, b in zip(vertices, vertices[1:] + vertices[:1], strict=True):
        e, ra, rb = difference(b, a), difference(query, a), difference(query, b)
        n = cross(ra, rb)
        trace = (Interval(dot(e, ra)) / root(dot(ra, ra) + radius2) -
                 Interval(dot(e, rb)) / root(dot(rb, rb) + radius2))
        value += (trace * Interval(dot(n, normal)) /
                  (factor * Interval(dot(n, n) + radius2 * dot(e, e))))
    return value


def nearest_error(value):
    center = float((value.lower + value.upper) / 2)
    words = tuple(map(Fraction, (math.nextafter(center, -math.inf), center,
                                 math.nextafter(center, math.inf))))
    if not words[0] <= value.lower <= value.upper <= words[-1]:
        raise ValueError("The neighbouring Binary64 words do not bracket the coefficient")
    return Interval(
        min(max(value.lower - word, word - value.upper, Fraction(0)) for word in words),
        min(max(abs(value.lower - word), abs(value.upper - word)) for word in words),
    )


# First chord strip, zero-based indices. Product coordinates from the Case lattice
# used by finite_wing_steady_test.cpp, in query, tangents, and directed Quad order.
CASES = (
    (0, 33, "1.2850e-13", "1.2851e-13", """
        0x1.0527e40d2a436p+1 -0x1.74f3a4e1feba7p+1 0x0p+0
        0x1p+2 0x0p+0 0x0p+0 -0x1.66819a3a0bf7bp+1 0x1p+2 0x0p+0
        0x1.07c1aa7f8a6c8p+1 -0x1.78ad74e01bd8fp+1 0x0p+0
        0x1.0289f2e4bdb38p+1 -0x1.7139d4e3e19bfp+1 0x0p+0
        0x1.02903f355a7ccp+1 -0x1.7139d4e3e19bfp+1 0x0p+0
        0x1.07c7ddfed3dd5p+1 -0x1.78ad74e01bd8fp+1 0x0p+0"""),
    (0, 162, "2.5005e-13", "2.5006e-13", """
        0x1.1943913070f7bp+1 0x1.91ab8df812cdbp+1 0x0p+0
        0x1p+2 0x0p+0 0x0p+0 0x1.66819a3a0bf7bp+1 0x1p+2 0x0p+0
        0x1.16d83d783cd13p+1 0x1.8e39d9cd73464p+1 0x0p+0
        0x1.1baaf9f7e0ba5p+1 0x1.951d4222b2552p+1 0x0p+0
        0x1.1bb0cec1dbc6ap+1 0x1.951d4222b2552p+1 0x0p+0
        0x1.16de29333aaeap+1 0x1.8e39d9cd73464p+1 0x0p+0"""),
    (5, 33, "1.6972e-13", "1.6973e-13", """
        0x1.0527e40d2a436p+1 -0x1.74f3a4e1feba7p+1 0x1.05084ac01faedp-2
        0x1p+2 0x0p+0 0x0p+0 -0x1.66819a3a0bf7bp+1 0x1p+2 -0x1.665a8349d55e1p-2
        0x1.07c1aa7f8a6c8p+1 -0x1.78ad74e01bd8fp+1 0x1.07a3dfccf1f97p-2
        0x1.0289f2e4bdb38p+1 -0x1.7139d4e3e19bfp+1 0x1.026cb5b34d644p-2
        0x1.02903f355a7ccp+1 -0x1.7139d4e3e19bfp+1 0x1.026cb5b34d644p-2
        0x1.07c7ddfed3dd5p+1 -0x1.78ad74e01bd8fp+1 0x1.07a3dfccf1f97p-2"""),
    (5, 162, "1.9910e-13", "1.9911e-13", """
        0x1.1943913070f7bp+1 0x1.91ab8df812cdbp+1 0x1.1921f66bbaed8p-2
        0x1p+2 0x0p+0 0x0p+0 0x1.66819a3a0bf7bp+1 0x1p+2 0x1.665a8349d55e1p-2
        0x1.16d83d783cd13p+1 0x1.8e39d9cd73464p+1 0x1.16b8d9932b254p-2
        0x1.1baaf9f7e0ba5p+1 0x1.951d4222b2552p+1 0x1.1b8b13444ab5cp-2
        0x1.1bb0cec1dbc6ap+1 0x1.951d4222b2552p+1 0x1.1b8b13444ab5cp-2
        0x1.16de29333aaeap+1 0x1.8e39d9cd73464p+1 0x1.16b8d9932b254p-2"""),
)


def main():
    radius = Fraction(float.fromhex("0x1.75e469fd8d403p-17"))
    budget = Fraction(1e-13)
    for angle, row, lower, upper, data in CASES:
        query, tangents, *vertices = (tuple(Fraction(float.fromhex(word)) for word in line.split())
                                      for line in data.strip().splitlines())
        error = nearest_error(coefficient(query, tangents, vertices, radius))
        if not budget < Fraction(lower) <= error.lower <= error.upper <= Fraction(upper):
            raise ValueError(f"Unexpected nearest error for {angle} degrees, row {row}")
        print(f"dihedral={angle} row={row} minimum_error={float(error.lower):.17g}")


if __name__ == "__main__":
    main()
