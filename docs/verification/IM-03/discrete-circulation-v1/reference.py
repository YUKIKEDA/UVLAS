# SPDX-License-Identifier: Apache-2.0
"""Independent line and wake surface integrals for the circulation solve fixture."""

from decimal import Decimal as D, localcontext
from pathlib import Path
from runpy import run_path

helper = run_path(str(Path(__file__).resolve().parents[1] / "p1-shared-action-v1/reference.py"))
cross = helper["cross"]
dot = helper["dot"]
subtract = helper["subtract"]
gauss = helper["gauss_unit_interval"]


def ring(point, y0):
    """Regularized straight contours using unscaled physical endpoint coordinates."""
    loop = [(D(x), D(y + y0), D(0)) for x, y in [(0, 0), (0, 1), (1, 1), (1, 0)]]
    total = D(0)
    for a, b in zip(loop, loop[1:] + loop[:1], strict=True):
        dl = subtract(b, a)
        start, end = subtract(point, a), subtract(point, b)
        numerator = cross(dl, start)
        factor = dot(dl, start) / (dot(start, start) + EPS2).sqrt()
        factor -= dot(dl, end) / (dot(end, end) + EPS2).sqrt()
        total += numerator[2] * factor / ((dot(numerator, numerator) + EPS2 * dot(dl, dl)) * PI4)
    return total


with localcontext() as ctx:
    ctx.prec = 60
    EPS2 = D(1) / 64
    PI4 = D(helper["PI_TIMES_FOUR"])
    p0 = (D("0.75"), D("0.5"), D(0))
    p1 = (D("0.625"), D("2.375"), D(0))
    a00, a01 = ring(p0, 0), ring(p0, 2)
    a10, a11 = ring(p1, 0), ring(p1, 2)
    b0, b1 = D("-0.008"), D("0.242")
    determinant = a00 * a11 - a01 * a10
    print("two bound panels", (b0 * a11 - b1 * a01) / determinant,
          (a00 * b1 - a10 * b0) / determinant)

    for order in (48, 72, 96):
        front, historical = D(0), D(0)
        rule = gauss(order)
        for x, wx in rule:
            for y, wy in rule:
                r2 = (p0[0] - (1 + x)) ** 2 + (p0[1] - y) ** 2
                den = r2 + EPS2
                k = (r2 - 2 * EPS2) * wx * wy / (PI4 * den ** 2 * den.sqrt())
                front += (1 - x) * k
                historical += x * k
        print("attached wake", order, (b0 + historical / 4) / (a00 + front))
