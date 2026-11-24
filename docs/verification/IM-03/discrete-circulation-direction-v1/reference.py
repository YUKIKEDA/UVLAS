# SPDX-License-Identifier: Apache-2.0
"""Independent direct finite-core area and line formulas for the circulation direction."""

from decimal import Decimal as D, localcontext
from pathlib import Path
from runpy import run_path

policy = run_path(str(Path(__file__).resolve().parents[1] / "p1-shared-action-v1/reference.py"))
cross, dot, subtract = (policy[name] for name in ("cross", "dot", "subtract"))
gauss = policy["gauss_unit_interval"]


def vector(*values):
    return tuple(D(str(value)) for value in values)


def perturb(base, direction, h):
    return tuple(x + h * dx for x, dx in zip(base, direction, strict=True))


def affine_direction(point):
    x, y, _ = point
    return (x / 8 + y / 16, -x / 32 + y / 8, x / 16 + y / 32)


class Fixture:
    def __init__(self, order):
        self.rule = gauss(order)
        self.eps2 = D(1) / 64
        self.pi4 = D(policy["PI_TIMES_FOUR"])
        self.point = vector(0.75, 0.5, 0)

    def ring(self, point, h, normal):
        loop = [vector(x, y, 0) for x, y in [(0, 0), (0, 1), (1, 1), (1, 0)]]
        loop = [perturb(p, affine_direction(p), h) for p in loop]
        total = D(0)
        for a, b in zip(loop, loop[1:] + loop[:1], strict=True):
            dl = subtract(b, a)
            start, end = subtract(point, a), subtract(point, b)
            num = cross(dl, start)
            factor = dot(dl, start) / (dot(start, start) + self.eps2).sqrt()
            factor -= dot(dl, end) / (dot(end, end) + self.eps2).sqrt()
            total += dot(normal, num) * factor / (
                (dot(num, num) + self.eps2 * dot(dl, dl)) * self.pi4
            )
        return total

    def equation(self, h, material):
        query_direction = affine_direction(self.point) if material else vector(-.03125, .0625, .125)
        point = perturb(self.point, query_direction, h)
        tx = perturb(vector(1, 0, 0), vector(.0625, -.0625, .125), h)
        ty = perturb(vector(0, 2, 0), vector(-.125, .125, -.0625), h)
        area = cross(tx, ty)
        length = dot(area, area).sqrt()
        normal = tuple(x / length for x in area)
        tx = perturb(vector(1, 0, 0), affine_direction(vector(1, 0, 0)), h)
        ty = perturb(vector(0, 1, 0), affine_direction(vector(0, 1, 0)), h)
        area = cross(tx, ty)
        front, history = D(0), D(0)
        for y, wy in self.rule:
            for x, wx in self.rule:
                vertex = vector(1 + x, y, 0)
                r = subtract(point, perturb(vertex, affine_direction(vertex), h))
                d2 = dot(r, r)
                den = d2 + self.eps2
                numerator = ((d2 - 2 * self.eps2) * dot(normal, area)
                             - 3 * dot(r, normal) * dot(r, area))
                kernel = numerator * wx * wy / (self.pi4 * den ** 2 * den.sqrt())
                front += (1 - x) * kernel
                history += x * kernel
        relative = perturb(vector(-.25, -.25, .008), vector(.125, -.0625, .0625), h)
        return self.ring(point, h, normal) + front, dot(normal, relative) - history / 4

    def direction(self, h, material):
        a, b = self.equation(D(0), material)
        gamma = -b / a
        ap, bp = self.equation(h, material)
        am, bm = self.equation(-h, material)
        forcing = ((ap - am) * gamma + bp - bm) / (2 * h)
        tangent = (-bp / ap + bm / am) / (2 * h)
        return forcing, tangent


with localcontext() as context:
    context.prec = 65
    for quadrature, step in ((48, "1e-12"), (72, "1e-15"), (96, "1e-18")):
        fixture = Fixture(quadrature)
        for material in (False, True):
            print(quadrature, step, "material" if material else "free",
                  *fixture.direction(D(step), material))
