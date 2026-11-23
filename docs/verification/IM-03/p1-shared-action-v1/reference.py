# SPDX-License-Identifier: Apache-2.0
"""Generate the shared P1 velocity and directional derivative fixture on stdout.

The reference evaluates the regularized surface integral directly. It shares
neither the product boundary/sheet decomposition nor its derivative code.
"""

import json
import math
from decimal import Decimal, localcontext
from fractions import Fraction

COORDINATES = (
    ("0", "0", "0"),
    ("0", "1", "0"),
    ("1", "1", "1/4"),
    ("1", "0", "0"),
    ("-3/4", "1/2", "-1/8"),
)
COORDINATE_DIRECTION = (
    ("1/16", "0", "-1/32"),
    ("0", "1/16", "1/32"),
    ("0", "-1/32", "1/16"),
    ("1/32", "-1/16", "0"),
    ("-1/32", "0", "1/16"),
)
QUERY = ("3/8", "7/16", "3/4")
QUERY_DIRECTION = ("-1/16", "1/32", "1/16")
RADIUS = "1/8"
RADIUS_DIRECTION = "1/64"
STRENGTH = (("3/4", "1/4"), ("-3/8", "1/8"), ("1/2", "-1/4"))
STRENGTH_DIRECTION = (("1/8", "-1/8"), ("1/4", "1/16"), ("-1/8", "1/8"))
FACETS = ((0, 2, 1), (0, 3, 2), (0, 1, 4))
TRACES = {
    "equal_trace": ((0, 1, 0), (0, 1, 1), (0, 0, 2)),
    "unequal_trace": ((0, 1, 0), (0, 1, 1), (1, 0, 2)),
}
ORDERS = (24, 36, 48)
DECIMAL_DIGITS = 60
PI_TIMES_FOUR = "12.566370614359172953850573533118011536788677597500423283899778"


def rational(text):
    number = Fraction(text)
    return Decimal(number.numerator) / Decimal(number.denominator)


def vector(values):
    return tuple(rational(value) for value in values)


def add(left, right):
    return tuple(x + y for x, y in zip(left, right, strict=True))


def subtract(left, right):
    return tuple(x - y for x, y in zip(left, right, strict=True))


def scale(values, factor):
    return tuple(value * factor for value in values)


def cross(left, right):
    x, y, z = left
    u, v, w = right
    return y * w - z * v, z * u - x * w, x * v - y * u


def dot(left, right):
    return sum(x * y for x, y in zip(left, right, strict=True))


def interpolate(weights, rows):
    return tuple(dot(weights, component) for component in zip(*rows, strict=True))


def legendre_pair(order, x):
    previous, current = Decimal(1), x
    for degree in range(2, order + 1):
        following = ((2 * degree - 1) * x * current - (degree - 1) * previous) / degree
        previous, current = current, following
    return previous, current


def gauss_unit_interval(order):
    """Find the Gauss nodes in Decimal. Float cosines only seed Newton iterations."""
    result = []
    for index in range(1, order // 2 + 1):
        angle = math.pi * (index - 0.25) / (order + 0.5)
        point = Decimal(str(math.cos(angle)))
        for _ in range(40):
            previous, polynomial = legendre_pair(order, point)
            slope = order * (previous - point * polynomial) / (1 - point * point)
            step = polynomial / slope
            point -= step
            if abs(step) < Decimal("1e-55"):
                break
        else:
            raise RuntimeError(f"Gauss node {index} of order {order} failed to converge")
        previous, polynomial = legendre_pair(order, point)
        slope = order * (previous - point * polynomial) / (1 - point * point)
        weight = 1 / ((1 - point * point) * slope * slope)
        result.extend((((1 - point) / 2, weight), ((1 + point) / 2, weight)))
    return result


def integrate(order, traces):
    query = vector(QUERY)
    query_direction = vector(QUERY_DIRECTION)
    radius = rational(RADIUS)
    radius_direction = rational(RADIUS_DIRECTION)
    points = tuple(map(vector, COORDINATES))
    points_direction = tuple(map(vector, COORDINATE_DIRECTION))
    strengths = tuple(map(vector, STRENGTH))
    strengths_direction = tuple(map(vector, STRENGTH_DIRECTION))

    velocity = [[Decimal(0), Decimal(0)] for _ in range(3)]
    derivative = [[Decimal(0), Decimal(0)] for _ in range(3)]
    rule = gauss_unit_interval(order)
    for vertices, indices in zip(FACETS, traces, strict=True):
        facet = [points[index] for index in vertices]
        facet_direction = [points_direction[index] for index in vertices]
        coefficients = [strengths[index] for index in indices]
        coefficients_direction = [strengths_direction[index] for index in indices]
        edge_a = subtract(facet[1], facet[0])
        edge_b = subtract(facet[2], facet[0])
        edge_a_direction = subtract(facet_direction[1], facet_direction[0])
        edge_b_direction = subtract(facet_direction[2], facet_direction[0])
        normal = cross(edge_a, edge_b)
        normal_direction = add(cross(edge_a_direction, edge_b), cross(edge_a, edge_b_direction))
        for u, weight_u in rule:
            for v, weight_v in rule:
                barycentric = ((1 - u) * (1 - v), u, (1 - u) * v)
                offset = subtract(query, interpolate(barycentric, facet))
                offset_direction = subtract(
                    query_direction, interpolate(barycentric, facet_direction)
                )
                squared_distance = dot(offset, offset)
                squared_direction = 2 * dot(offset, offset_direction)
                core_squared = radius * radius
                core_squared_direction = 2 * radius * radius_direction
                distance = squared_distance + core_squared
                distance_direction = squared_direction + core_squared_direction
                normal_dot_offset = dot(offset, normal)
                normal_dot_direction = dot(offset_direction, normal) + dot(offset, normal_direction)
                radial_factor = squared_distance - 2 * core_squared
                radial_direction = squared_direction - 2 * core_squared_direction

                numerator = subtract(
                    scale(normal, radial_factor), scale(offset, 3 * normal_dot_offset)
                )
                numerator_direction = subtract(
                    add(scale(normal, radial_direction), scale(normal_direction, radial_factor)),
                    scale(
                        add(
                            scale(offset, normal_dot_direction),
                            scale(offset_direction, normal_dot_offset),
                        ),
                        3,
                    ),
                )
                weight = weight_u * weight_v * (1 - u)
                weight /= Decimal(PI_TIMES_FOUR) * distance * distance * distance.sqrt()
                numerator_direction = subtract(
                    numerator_direction,
                    scale(numerator, Decimal("2.5") * distance_direction / distance),
                )
                value = interpolate(barycentric, coefficients)
                direction = interpolate(barycentric, coefficients_direction)
                for axis in range(3):
                    for plane in range(2):
                        velocity[axis][plane] += weight * value[plane] * numerator[axis]
                        derivative[axis][plane] += weight * (
                            direction[plane] * numerator[axis]
                            + value[plane] * numerator_direction[axis]
                        )
    return {"velocity": velocity, "jvp": derivative}


def maximum_delta(first, later):
    return max(
        abs(x - y)
        for key, vectors in first.items()
        for a, b in zip(vectors, later[key], strict=True)
        for x, y in zip(a, b, strict=True)
    )


def maximum_precision_records():
    cases = {}
    convergence = {}
    for name, traces in TRACES.items():
        ordered = [integrate(order, traces) for order in ORDERS]
        convergence[name] = [
            f"{maximum_delta(first, later):.8E}"
            for first, later in zip(ordered[:-1], ordered[1:], strict=True)
        ]
        cases[name] = {
            key: [[f"{value:.24g}" for value in pair] for pair in coordinates]
            for key, coordinates in ordered[-1].items()
        }
    return {
        "schema": "uvlas.p1-shared-action-oracle/1",
        "quadrature": {"decimal_digits": DECIMAL_DIGITS, "orders": ORDERS},
        "input": {
            "vertices": COORDINATES,
            "vertex_jvp": COORDINATE_DIRECTION,
            "query": QUERY,
            "query_jvp": QUERY_DIRECTION,
            "radius": RADIUS,
            "radius_jvp": RADIUS_DIRECTION,
            "strength": STRENGTH,
            "strength_jvp": STRENGTH_DIRECTION,
            "facets": FACETS,
            "traces": TRACES,
        },
        "convergence_max_component": convergence,
        "cases": cases,
    }


def main():
    with localcontext() as context:
        context.prec = DECIMAL_DIGITS
        print(json.dumps(maximum_precision_records(), indent=2))


if __name__ == "__main__":
    main()
