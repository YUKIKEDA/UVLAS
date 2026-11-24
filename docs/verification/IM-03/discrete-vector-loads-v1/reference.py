# SPDX-License-Identifier: Apache-2.0
"""Independent Decimal line integrals, wake-area integrals and midpoint/triangle loads."""

from decimal import Decimal as D, localcontext
from pathlib import Path
from runpy import run_path

helper = run_path(str(Path(__file__).resolve().parents[1] / "p1-shared-action-v1/reference.py"))
cross, dot, subtract = (helper[k] for k in ("cross", "dot", "subtract"))
gauss = helper["gauss_unit_interval"]


def vector(*entries):
    return tuple(D(str(value)) for value in entries)


def plus(a, b):
    return tuple(x + y for x, y in zip(a, b, strict=True))


def scale(value, factor):
    return tuple(x * factor for x in value)


class Fixture:
    def __init__(self, order):
        self.rule = gauss(order)
        self.eps2 = D(1) / 64
        self.pi4 = D(helper["PI_TIMES_FOUR"])
        self.points = [vector(x, y, 0) for x, y in
                       [(0, 0), (0, 1), (1, 1), (1, 0), (2, 1), (2, 0)]]
        self.rings = [[0, 1, 2, 3], [3, 2, 4, 5]]
        self.origin = vector(.25, -.125, .5)
        self.rho = D("1.25")
        self.freestream = vector(0, 0, .508)

    def line_velocity(self, query, start, end):
        tangent = subtract(end, start)
        r_a, r_b = subtract(query, start), subtract(query, end)
        numerator = cross(tangent, r_a)
        scalar = dot(tangent, r_a) / (dot(r_a, r_a) + self.eps2).sqrt()
        scalar -= dot(tangent, r_b) / (dot(r_b, r_b) + self.eps2).sqrt()
        denominator = self.pi4 * (dot(numerator, numerator) + self.eps2 * dot(tangent, tangent))
        return scale(numerator, scalar / denominator)

    def ring_velocity(self, point, ring, skip=None):
        total = vector(0, 0, 0)
        for left, right in zip(ring, ring[1:] + ring[:1], strict=True):
            if {left, right} != skip:
                total = plus(total, self.line_velocity(point, self.points[left], self.points[right]))
        return total

    def wake_velocity(self, point):
        # The reference fixture and all its load queries lie in the disk's plane.
        front, history = D(0), D(0)
        for x, wx in self.rule:
            for y, wy in self.rule:
                r = subtract(point, (1 + x, y, D(0)))
                r2 = dot(r, r)
                den = r2 + self.eps2
                k = (r2 - 2 * self.eps2) * wx * wy / (self.pi4 * den ** 2 * den.sqrt())
                front += (1 - x) * k
                history += x * k
        return front, history

    def circulation(self, wake):
        p0 = vector(.75, .5, 0)
        a = self.ring_velocity(p0, self.rings[0])[2]
        if wake:
            front, history = self.wake_velocity(p0)
            return ((D("-.008") + history / 4) / (a + front),)
        p1 = vector(1.625, .375, 0)
        b = self.ring_velocity(p0, self.rings[1])[2]
        c = self.ring_velocity(p1, self.rings[0])[2]
        d = self.ring_velocity(p1, self.rings[1])[2]
        f, g = D("-.008"), D(".242")
        det = a * d - b * c
        return (f * d - g * b) / det, (a * g - c * f) / det

    def forces(self, wake):
        gamma = self.circulation(wake)
        rings = self.rings[:len(gamma)]
        supplied_rates = [D(".375") + D(2) ** -57, -D(".125")]
        motions = {10: vector(.125, -.25, .375), 11: vector(.25, .125, .125),
                   12: vector(.25, 0, .125), 13: vector(-.125, .25, -.25),
                   15: vector(.125, -.25, .125), 16: vector(0, .5, .25),
                   17: vector(.125, .25, -.125)}
        edges = {10: (0, 1), 11: (1, 2), 12: (2, 3), 13: (0, 3),
                 15: (2, 4), 16: (4, 5), 17: (3, 5)}
        force = moment = vector(0, 0, 0)
        nodes = [vector(0, 0, 0) for _ in range(6)]

        def distribute(position, value, vertices):
            nonlocal force, moment
            force = plus(force, value)
            moment = plus(moment, cross(subtract(position, self.origin), value))
            for vertex in vertices:
                nodes[vertex] = plus(nodes[vertex], scale(value, D(1) / len(vertices)))

        edge_rows = {}
        for edge, endpoints in edges.items():
            if wake and edge not in [10, 11, 13]:
                continue
            left, right = endpoints
            position = scale(plus(self.points[left], self.points[right]), D(".5"))
            strength = D(0)
            velocity = self.freestream
            for ring, explained in zip(rings, gamma, strict=True):
                segments = list(zip(ring, ring[1:] + ring[:1], strict=True))
                strength += explained * (int((left, right) in segments) - int((right, left) in segments))
                velocity = plus(velocity, scale(self.ring_velocity(position, ring, {left, right}), explained))
            if wake:
                front, history = self.wake_velocity(position)
                velocity = plus(velocity, (D(0), D(0), front * gamma[0] - history / 4))
            relative = subtract(velocity, motions[edge])
            dl = subtract(self.points[right], self.points[left])
            value = scale(cross(relative, dl), self.rho * strength)
            edge_rows[edge] = value
            distribute(position, value, endpoints)
        for panel, corners in enumerate(rings):
            for indices in [(0, 2, 1), (0, 3, 2)]:
                vertices = tuple(corners[j] for j in indices)
                a, b, c = (self.points[j] for j in vertices)
                position = scale(plus(plus(a, b), c), D(1) / 3)
                area = scale(cross(subtract(b, a), subtract(c, a)), D(".5"))
                distribute(position, scale(area, self.rho * supplied_rates[panel]), vertices)
        return edge_rows, force, moment, nodes


with localcontext() as context:
    context.prec = 60
    for order in (48, 72, 96):
        model = Fixture(order)
        for wake in (True, False):
            if not wake and order != 96:
                continue
            edge_rows, force, moment, nodes = model.forces(wake)
            print("wake" if wake else "two bound", order)
            for i, row in edge_rows.items():
                print(i, *(f"{x:.26g}" for x in row))
            print("total", *(f"{x:.26g}" for x in force))
            print("moment", *(f"{x:.26g}" for x in moment))
            print("node 1", *(f"{x:.26g}" for x in nodes[0]))
