# SPDX-License-Identifier: Apache-2.0
"""Independent regularized P0 line and W-P1 surface oracle at wake material vertices."""

from decimal import Decimal as D, localcontext
from pathlib import Path
from runpy import run_path

helper = run_path(str(Path(__file__).resolve().parents[1] / "p1-shared-action-v1/reference.py"))
gauss = helper["gauss_unit_interval"]
PI4 = D(helper["PI_TIMES_FOUR"])
POINTS = tuple(tuple(map(D, pair)) for pair in ((0, 0), (0, 1), (1, 1), (1, 0)))


def ring(point, radius, corners=POINTS):
    result = D(0)
    for start, end in zip(corners, corners[1:] + corners[:1], strict=True):
        dx, dy = (b - a for a, b in zip(start, end, strict=True))
        xa, ya = (q - a for q, a in zip(point, start, strict=True))
        xb, yb = (q - b for q, b in zip(point, end, strict=True))
        c = dx * ya - dy * xa
        edge = (dx * xa + dy * ya) / (xa * xa + ya * ya + radius**2).sqrt()
        edge -= (dx * xb + dy * yb) / (xb * xb + yb * yb + radius**2).sqrt()
        result += c * edge / (PI4 * (c * c + radius**2 * (dx * dx + dy * dy)))
    return result


def sheet(q, radius, rule):
    leading, trailing, dleading, dtrailing = (D(0),) * 4
    for x, wx in rule:
        for y, wy in rule:
            r2 = (q[0] - 1 - x) ** 2 + (q[1] - y) ** 2
            den = r2 + radius**2
            k = (r2 - 2 * radius**2) / (PI4 * den**2 * den.sqrt())
            dk = (-4 * radius - 5 * radius * (r2 - 2 * radius**2) / den) / (
                PI4 * den**2 * den.sqrt()
            )
            leading += (1 - x) * k * wx * wy
            trailing += x * k * wx * wy
            dleading += (1 - x) * dk * wx * wy
            dtrailing += x * dk * wx * wy
    return leading, trailing, dleading, dtrailing


def run(order):
    rule = gauss(order)
    eps, gamma, old = D(1) / 8, D(1) / 2, -D(1) / 4
    deps, dgamma, dold = D(1) / 32, D(1) / 8, D(5) / 16
    rows = []
    for i, q in enumerate(((1, 0), (2, 1))):
        q = tuple(map(D, q))
        lead, trail, dlead, dtrail = sheet(q, eps, rule)
        a = ring(q, eps)
        h = D("1e-22")
        da = (ring(q, eps + h) - ring(q, eps - h)) / (2 * h)
        common = a + lead
        p = (da + dlead) * gamma + dtrail * old
        ambient = (D(1), D(2))[i]
        phase = -PI4 / 64 if i == 0 else D(0)
        v = ambient + common * gamma + trail * old
        jvp = D(1) / 8 + phase + common * dgamma + trail * dold + p * deps
        rows.append((v, jvp, common, trail, p))
    return rows


def receiver_alias_dual(order):
    """Move canonical TE point 4 in x; its aliased receiving point moves with it."""
    rule = gauss(order)
    eps, gamma, old = D(1) / 8, D(1) / 2, -D(1) / 4
    values = []
    for qx, qy, dq in ((D(1), D(0), D(1)), (D(2), D(1), D(0))):
        h = D("1e-22")
        plus = POINTS[:3] + ((D(1) + h, D(0)),)
        minus = POINTS[:3] + ((D(1) - h, D(0)),)
        bound = (ring((qx + h * dq, qy), eps, plus) - ring((qx - h * dq, qy), eps, minus)) / (2 * h)
        wake = D(0)
        for u, wu in rule:
            for t, wt in rule:
                for right in (False, True):
                    x, y = (u, u * t) if right else (u * t, u)
                    rx, ry = qx - 1 - x, qy - y
                    r2 = rx**2 + ry**2
                    den = r2 + eps**2
                    k = (r2 - 2 * eps**2) / (PI4 * den**2 * den.sqrt())
                    dk = (1 - D(5) / 2 * (r2 - 2 * eps**2) / den) / (PI4 * den**2 * den.sqrt())
                    # Both affine triangles have N_0=1-u. Only x>=y has d(det)/d(x4)=-1.
                    dr2 = 2 * rx * (dq - (1 - u))
                    density = (1 - x) * gamma + x * old
                    wake += density * (dk * dr2 - (k if right else 0)) * u * wu * wt
        values.append(bound * gamma + wake + dq * PI4 / 8)
    return values[0] / 2 - 3 * values[1] / 8


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 75
        old = run(64)
        result = run(88)
        coarse = receiver_alias_dual(72)
        fine = receiver_alias_dual(96)
        assert abs(coarse - fine) < D("1e-20")
        print("canonical TE point 4 x dual", fine)
        print("value, JVP, bound dual, old-trace dual, core-radius dual:", result)
        print(
            "max refinement difference",
            max(
                abs(a - b)
                for r, s in zip(old, result, strict=True)
                for a, b in zip(r, s, strict=True)
            ),
        )
