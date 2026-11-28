"""Check a rectangular full-thickness surface load using rational coordinates."""

import json
from fractions import Fraction
from itertools import pairwise, product

type Vector = tuple[Fraction, Fraction, Fraction]
type Interval = tuple[Fraction, Fraction]
type Rectangle = tuple[Interval, Interval]

THICKNESS = Fraction(127, 50_000)
CENTRE_X = Fraction(1)
LOADED_END = Fraction(4)
FIXED_END = Fraction(2)
ORIGIN = (Fraction(2), Fraction(1), Fraction(3))
SHIFT = (Fraction(1, 2), Fraction(-3, 4), Fraction(1))
ROTATION = (Fraction(1, 4), Fraction(-1, 2), Fraction(3, 4))
RULE = ((Fraction(), Fraction(1, 6)), (Fraction(1, 2), Fraction(2, 3)),
        (Fraction(1), Fraction(1, 6)))


def cross(a: Vector, b: Vector) -> Vector:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def dot(a: Vector, b: Vector) -> Fraction:
    return sum((x * y for x, y in zip(a, b, strict=True)), Fraction())


def require(condition: bool, observation: str) -> None:
    if not condition:
        raise ValueError(observation)


def integrate(rectangles: list[Rectangle], traction: Fraction, end_y: Fraction,
              reverse_edge: bool = False) -> tuple[Fraction, ...]:
    """Integrate force, moments, rigid work, and two placement-sensitive probes."""
    area = force_y = moment_x = moment_z = work = axial = depth = Fraction()
    for (x0, x1), (z0, z1) in rectangles:
        for (u, w_u), (v, w_v) in product(RULE, repeat=2):
            x = x0 + (x1 - x0) * u
            z = z0 + (z1 - z0) * (1 - v if reverse_edge else v)
            weight = (x1 - x0) * (z1 - z0) * w_u * w_v
            radius = (x - ORIGIN[0], end_y - ORIGIN[1], z - ORIGIN[2])
            force_at_point = weight * traction
            area += weight
            force_y += force_at_point
            moment_x -= force_at_point * radius[2]
            moment_z += force_at_point * radius[0]
            work += force_at_point * (SHIFT[1] + cross(ROTATION, radius)[1])
            axial += force_at_point * (end_y - FIXED_END) / (LOADED_END - FIXED_END)
            depth += force_at_point * (2 * (x - CENTRE_X) / THICKNESS) ** 2
    return area, force_y, moment_x, moment_z, work, axial, depth


def edge_rule() -> list[tuple[Fraction, Fraction]]:
    numerators = (0, 5856, 23293, 51921, 91101, 139958, 197400, 262144, 332744,
                  407623, 485108, 563468, 640953, 715832, 786432, 851176, 908618,
                  957475, 996655, 1025283, 1042720, 1048576)
    nodes = [Fraction(numerator, 2**20) for numerator in numerators]
    weights = []
    for node in nodes:
        coefficients = [Fraction(1)]
        divisor = Fraction(1)
        for other in nodes:
            if node == other:
                continue
            product_coefficients = [Fraction()] * (len(coefficients) + 1)
            for power, coefficient in enumerate(coefficients):
                product_coefficients[power] -= coefficient * other
                product_coefficients[power + 1] += coefficient
            coefficients = product_coefficients
            divisor *= node - other
        # Integrate the Lagrange cardinal polynomial on [0, 1].
        weights.append(sum(coefficient / (power + 1)
                           for power, coefficient in enumerate(coefficients)) / divisor)
    rule = list(zip(nodes, weights, strict=True))
    require(all(weight > 0 for _, weight in rule), 'positive edge weights')
    for power in range(len(nodes)):
        moment = sum(weight * node**power for node, weight in rule)
        require(moment == Fraction(1, power + 1), f'edge moment of degree {power}')
    return rule


def main() -> None:
    young = Fraction(200_000_000_000)
    stretch = Fraction(1025, 1024)
    traction = young * stretch * (stretch**2 - 1) / 2
    require(traction == Fraction(102550048828125, 524288), 'traction')

    area = 2 * THICKNESS
    force = (Fraction(), area * traction, Fraction())
    require(force[1] == Fraction(4167633984375, 4194304), 'force')
    centre = (CENTRE_X, LOADED_END, Fraction(4))
    radius = (centre[0] - ORIGIN[0], centre[1] - ORIGIN[1], centre[2] - ORIGIN[2])
    moment = cross(radius, force)
    require(moment == (-force[1], Fraction(), -force[1]), 'moment')
    work = dot(SHIFT, force) + dot(ROTATION, moment)
    require(work == -Fraction(7, 4) * force[1], 'rigid work')

    x0, x1 = CENTRE_X - THICKNESS / 2, CENTRE_X + THICKNESS / 2
    z0, z1 = Fraction(3), Fraction(5)
    cells = list(product(pairwise((x0, x0 + THICKNESS / 3, x1)),
                         pairwise((z0, Fraction(17, 5), z1))))
    require(len(cells) == 4, 'number of side pieces')
    expected = (area, force[1], moment[0], moment[2], work, force[1], force[1] / 3)
    whole = integrate([((x0, x1), (z0, z1))], traction, LOADED_END)
    require(whole == expected, 'whole side')
    require(integrate(cells, traction, LOADED_END) == expected, 'partitioned side')
    require(integrate(cells, traction, LOADED_END, reverse_edge=True) == expected,
            'reversed edge')

    require(integrate(cells[1:], traction, LOADED_END)[1] < force[1], 'missing piece')
    require(integrate(cells + cells[:1], traction, LOADED_END)[1] > force[1], 'repeated piece')
    require(integrate(cells, -traction, LOADED_END)[1] == -force[1], 'flipped force')
    require((z1 - z0) * traction != force[1], 'line measure')
    opposite_end = integrate(cells, traction, FIXED_END)
    require(opposite_end[:5] == expected[:5] and opposite_end[5] == 0 < expected[5],
            'opposite end')
    centreline_work = force[1] * (2 * (centre[0] - CENTRE_X) / THICKNESS) ** 2
    require(centreline_work == 0 < expected[6], 'depth centreline')
    require(dot(ROTATION, cross(centre, force)) != dot(ROTATION, moment), 'wrong moment origin')

    authored = Fraction(1, 3)
    authored_reversed_z = z1 + (z0 - z1) * authored
    owner_z = z0 + (z1 - z0) * (1 - authored)
    require(authored_reversed_z == owner_z == Fraction(13, 3), 'reversed point')
    require(z0 + (z1 - z0) * authored == Fraction(11, 3) != owner_z, 'omitted reversal')

    extension = (LOADED_END - FIXED_END) * (stretch - 1)
    rule = edge_rule()
    result = {
        'area_m2': str(area),
        'traction_y_pa': str(traction),
        'force_n': list(map(str, force)),
        'moment_nm': list(map(str, moment)),
        'rigid_work_j': str(work),
        'axially_affine_work_j': str(expected[5]),
        'quadratic_depth_work_j': str(expected[6]),
        'end_extension_m': str(extension),
        'uniform_end_extension_work_j': str(force[1] * extension),
        'reversed_point_z_m': str(owner_z),
        'exact_edge_rule': [[str(node), str(weight)] for node, weight in rule],
        'negative_controls': ['missing_piece', 'repeated_piece', 'flipped_force',
                              'line_measure', 'opposite_end', 'depth_centreline',
                              'wrong_moment_origin', 'omitted_edge_reversal'],
    }
    print(json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
