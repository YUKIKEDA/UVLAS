# SPDX-License-Identifier: Apache-2.0
"""Independent load differences with circulation and circulation-rate directions supplied."""

from contextlib import redirect_stdout
from decimal import Decimal as D, localcontext
from io import StringIO
from pathlib import Path
from runpy import run_path

with redirect_stdout(StringIO()):
    predecessor = run_path(
        str(Path(__file__).resolve().parents[1] / "discrete-vector-loads-v1/reference.py")
    )
Fixture, vector, plus, scale, cross, dot, subtract = (
    predecessor[name] for name in
    ("Fixture", "vector", "plus", "scale", "cross", "dot", "subtract")
)


def point_direction(point):
    x, y, _ = point
    return x / 8 + y / 16, -x / 32 + y / 8, x / 16 + y / 32


def at(point, direction, h):
    return plus(point, scale(direction, h))


def wake_induction(model, h, point, gamma):
    x = at(vector(1, 0, 0), point_direction(vector(1, 0, 0)), h)
    y = at(vector(0, 1, 0), point_direction(vector(0, 1, 0)), h)
    oriented_area = cross(x, y)
    result = vector(0, 0, 0)
    for u, wu in model.rule:
        for v, wv in model.rule:
            w = (D(1) + u, v, D(0))
            r = subtract(point, at(w, point_direction(w), h))
            r2 = dot(r, r)
            den = r2 + model.eps2
            kernel = subtract(scale(oriented_area, r2 - 2 * model.eps2),
                              scale(r, 3 * dot(r, oriented_area)))
            trace = (1 - u) * gamma - u / 4
            result = plus(result, scale(
                kernel, trace * wu * wv / (model.pi4 * den ** 2 * den.sqrt())
            ))
    return result


def sample(model, h, use_wake, base_gamma):
    points = [at(p, point_direction(p), h) for p in model.points]
    rings = model.rings[:len(base_gamma)]
    dgamma = [D(".1875") + D(2) ** -56, -D(".15625")]
    gamma = [g + h * dgamma[i] for i, g in enumerate(base_gamma)]
    rates = [D(".375") + D(2) ** -57 + h * (-D(".25") + D(2) ** -58),
             -D(".125") + h * D(".0625")]
    ambient = at(model.freestream, vector(.125, -.0625, .09375), h)
    edges = {10: (0, 1), 11: (1, 2), 12: (2, 3), 13: (0, 3),
             15: (2, 4), 16: (4, 5), 17: (3, 5)}
    motions = {10: vector(.125, -.25, .375), 11: vector(.25, .125, .125),
               12: vector(.25, 0, .125), 13: vector(-.125, .25, -.25),
               15: vector(.125, -.25, .125), 16: vector(0, .5, .25),
               17: vector(.125, .25, -.125)}
    result = {}
    force = moment = vector(0, 0, 0)
    nodes = [vector(0, 0, 0) for _ in points]

    def distribute(position, value, vertices):
        nonlocal force, moment
        force = plus(force, value)
        moment = plus(moment, cross(subtract(position, model.origin), value))
        for vertex in vertices:
            nodes[vertex] = plus(nodes[vertex], scale(value, D(1) / len(vertices)))

    for edge, (left, right) in edges.items():
        if use_wake and edge not in (10, 11, 13):
            continue
        position = scale(plus(points[left], points[right]), D(".5"))
        circulation = D(0)
        velocity = ambient
        for ring, g in zip(rings, gamma, strict=True):
            for a, b in zip(ring, ring[1:] + ring[:1], strict=True):
                if {a, b} == {left, right}:
                    circulation += g if a == left else -g
                else:
                    induced = model.line_velocity(position, points[a], points[b])
                    velocity = plus(velocity, scale(induced, g))
        if use_wake:
            velocity = plus(velocity, wake_induction(model, h, position, gamma[0]))
        motion = at(motions[edge], vector(.03125, -.0625, D(edge) / 256), h)
        load = scale(cross(subtract(velocity, motion), subtract(points[right], points[left])),
                     model.rho * circulation)
        distribute(position, load, (left, right))
        result[f"edge {edge} position"] = position
        result[f"edge {edge} force"] = load

    for panel, ring in enumerate(rings):
        for facet, local in enumerate(((0, 2, 1), (0, 3, 2))):
            ids = tuple(ring[i] for i in local)
            x, y, z = (points[i] for i in ids)
            center = scale(plus(plus(x, y), z), D(1) / 3)
            area = scale(cross(subtract(y, x), subtract(z, x)), D(".5"))
            load = scale(area, model.rho * rates[panel])
            distribute(center, load, ids)
            result[f"facet {panel}:{facet} position"] = center
            result[f"facet {panel}:{facet} force"] = load
    result["total"] = force
    result["moment"] = moment
    result["node 1"] = nodes[0]
    return result


def difference(model, use_wake, step):
    gamma = model.circulation(use_wake)
    upper, lower = (sample(model, step * sign, use_wake, gamma) for sign in (1, -1))
    return {key: scale(subtract(upper[key], lower[key]), 1 / (2 * step)) for key in upper}


with localcontext() as context:
    context.prec = 65
    for use_wake in (True, False):
        coarse = difference(Fixture(72), use_wake, D("1e-15"))
        reference = difference(Fixture(96), use_wake, D("1e-18"))
        print("wake" if use_wake else "two bound")
        for key, value in reference.items():
            print(key, *(f"{x:.24g}" for x in value))
        for key in ("total", "moment"):
            print(key, "change", max(abs(x) for x in subtract(reference[key], coarse[key])))
