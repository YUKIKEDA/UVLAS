# SPDX-License-Identifier: Apache-2.0
"""Independent scalar objective differences on the direct vector-load reference."""

from contextlib import redirect_stdout
from decimal import Decimal as D, localcontext
from io import StringIO
from pathlib import Path
from runpy import run_path

with redirect_stdout(StringIO()):
    source = run_path(str(
        Path(__file__).resolve().parents[1] / "discrete-vector-load-direction-v1/reference.py"
    ))
Fixture, sample, vector, dot = (source[name] for name in ("Fixture", "sample", "vector", "dot"))


def objective(values):
    weights = {
        "edge 10 position": vector(.125, .0625, -.125),
        "edge 10 force": vector(.25, -.125, .5),
        "facet 0:0 position": vector(-.25, .125, .0625),
        "facet 0:0 force": vector(.0625, .25, -.125),
        "node 1": vector(-.375, .125, .25),
        "total": vector(.125, -.25, D(".375") + D(2) ** -57),
        "moment": vector(-.125, .25, .0625),
    }
    return sum(dot(values[key], dual) for key, dual in weights.items())


def difference(wake, count, step):
    model = Fixture(count)
    gamma = model.circulation(wake)
    left = objective(sample(model, -step, wake, gamma))
    right = objective(sample(model, step, wake, gamma))
    return (right - left) / (2 * step)


with localcontext() as context:
    context.prec = 65
    for wake in (True, False):
        first = difference(wake, 72, D("1e-15"))
        second = difference(wake, 96, D("1e-18"))
        print("wake" if wake else "two bound", second, "change", abs(second - first))
