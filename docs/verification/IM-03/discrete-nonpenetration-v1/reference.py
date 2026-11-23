# SPDX-License-Identifier: Apache-2.0
"""Independent finite-core line and wake surface reference for a planar UVLM row."""

from decimal import Decimal, localcontext
from pathlib import Path
from runpy import run_path

helper = run_path(str(Path(__file__).resolve().parents[1] / "p1-shared-action-v1" / "reference.py"))
gauss = helper["gauss_unit_interval"]
vector = helper["vector"]
cross = helper["cross"]
dot = helper["dot"]
subtract = helper["subtract"]

with localcontext() as ctx:
    ctx.prec = 60
    radius2 = Decimal(1) / 64
    pi4 = Decimal(helper["PI_TIMES_FOUR"])
    point = vector(("3/4", "1/2", "0"))
    corners = [vector(p) for p in [("0", "0", "0"), ("0", "1", "0"),
                                   ("1", "1", "0"), ("1", "0", "0")]]
    bound = Decimal(0)
    for a, b in zip(corners, corners[1:] + corners[:1], strict=True):
        chord = subtract(b, a)
        start, end = subtract(point, a), subtract(point, b)
        normal = cross(chord, start)
        denominator = dot(normal, normal) + radius2 * dot(chord, chord)
        edge = dot(chord, start) / (dot(start, start) + radius2).sqrt()
        edge -= dot(chord, end) / (dot(end, end) + radius2).sqrt()
        bound += normal[2] * edge / (pi4 * denominator)

    gamma = Decimal(3) / 4 + Decimal(2) ** -60
    history = -Decimal(1) / 4
    flow_minus_body = Decimal(1) / 125
    for order in (48, 72, 96):
        rule = gauss(order)
        wake_front, wake_old = Decimal(0), Decimal(0)
        for x, wx in rule:
            for y, wy in rule:
                r2 = (point[0] - (1 + x)) ** 2 + (point[1] - y) ** 2
                denominator = r2 + radius2
                value = (r2 - 2 * radius2) / (denominator ** 2 * denominator.sqrt() * pi4)
                wake_front += wx * wy * (1 - x) * value
                wake_old += wx * wy * x * value
        residual = (bound + wake_front) * gamma + wake_old * history + flow_minus_body
        action = -(bound + wake_front) / 2
        print(order, f"bound={bound:.40g}", f"front={wake_front:.40g}", f"old={wake_old:.40g}")
        print(f"residual={residual:.40g} action={action:.40g}")
