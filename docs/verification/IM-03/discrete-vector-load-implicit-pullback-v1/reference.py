# SPDX-License-Identifier: Apache-2.0
"""Differentiate a load objective after re-solving the perturbed discrete circulation."""

from contextlib import redirect_stdout
from decimal import Decimal as D, localcontext
from io import StringIO
from pathlib import Path
from runpy import run_path

with redirect_stdout(StringIO()):
    previous = run_path(str(
        Path(__file__).resolve().parents[1] / "discrete-vector-load-pullback-v1/reference.py"
    ))
    direction = previous["source"]
Fixture, sample, objective = (previous[key] for key in ("Fixture", "sample", "objective"))
at, cross, dot, point_direction, plus, scale, subtract, vector, wake_induction = (
    direction[key] for key in
    ("at", "cross", "dot", "point_direction", "plus", "scale", "subtract", "vector", "wake_induction")
)


def circulation(model, wake, h):
    rows = (
        (vector(.75, .5, 0), vector(.03125, -.046875, .0234375),
         vector(.03125, .0625, .125), vector(-.09375, .015625, -.0625),
         vector(0, 0, .5), vector(.0625, -.03125, .125)),
        (vector(1.625, .375, 0), vector(-.046875, .0625, -.0390625),
         vector(-.0625, .03125, -.125), vector(.125, -.03125, .0625),
         vector(0, 0, .75), vector(-.03125, .0625, -.09375)),
    )
    points = [at(point, point_direction(point), h) for point in model.points]
    ambient = at(model.freestream, vector(.125, -.0625, .09375), h)
    matrix, right = [], []
    for q, dq, dx, dy, velocity, dv in rows[:1 if wake else 2]:
        query = at(q, dq, h)
        normal = cross(at(vector(1, 0, 0), dx, h), at(vector(0, 1, 0), dy, h))
        relative = subtract(ambient, at(velocity, dv, h))
        coefficients = []
        for ring in model.rings[:1 if wake else 2]:
            induced = vector(0, 0, 0)
            for left, right_point in zip(ring, ring[1:] + ring[:1], strict=True):
                induced = plus(induced, model.line_velocity(query, points[left], points[right_point]))
            coefficients.append(dot(normal, induced))
        if wake:
            inherited = wake_induction(model, h, query, D(0))
            leading = subtract(wake_induction(model, h, query, D(1)), inherited)
            relative = plus(relative, inherited)
            coefficients[0] += dot(normal, leading)
        matrix.append(coefficients)
        right.append(-dot(normal, relative))
    if wake:
        return (right[0] / matrix[0][0],)
    (a, b), (c, d) = matrix
    f, g = right
    det = a * d - b * c
    return ((f * d - g * b) / det, (a * g - c * f) / det)


def value(model, wake, h):
    # Cancel the earlier explicit dGamma increment in sample so the load formula
    # evaluates at the circulation solved independently for this perturbation.
    explicit = (D(".1875") + D(2) ** -56, -D(".15625"))
    solved = tuple(gamma - h * explicit[i] for i, gamma in enumerate(circulation(model, wake, h)))
    return objective(sample(model, h, wake, solved))


def difference(wake, order, step):
    model = Fixture(order)
    return (value(model, wake, step) - value(model, wake, -step)) / (2 * step)


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 65
        for with_wake in (True, False):
            coarse = difference(with_wake, 72, D("1e-15"))
            fine = difference(with_wake, 96, D("1e-18"))
            print("wake" if with_wake else "two bound", fine, "change", abs(fine - coarse))