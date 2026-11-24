# SPDX-License-Identifier: Apache-2.0
"""Direct same-core finite-strip induction at both end receivers."""

import json
from decimal import Decimal as D
from decimal import localcontext
from pathlib import Path
from runpy import run_path

REFERENCE = Path(__file__).resolve().parents[1] / "discrete-free-wake-velocity-v1/reference.py"
previous = run_path(str(REFERENCE))
EPSILON = D(1) / 8


def velocity(count, receiver, rule):
    value = D(1) / 8 + previous["ring"](receiver, EPSILON) / 2
    leading = D(1) / 2
    for row in range(count):
        trailing = -D(1) / 4 + D(row) / 64 + D(2) ** -58
        x, y = receiver
        first, second, *_ = previous["sheet"]((x - row, y), EPSILON, rule)
        value += leading * first + trailing * second
        leading = trailing
    return value


def main():
    with localcontext() as context:
        context.prec = 75
        low, high = (previous["gauss"](order) for order in (64, 88))
        cases = {}
        for count in (2, 4, 8):
            receivers = [(D(2), D(0)), (D(count + 1), D(1))]
            values = [velocity(count, receiver, high) for receiver in receivers]
            nearby = [velocity(count, receiver, low) for receiver in receivers]
            cases[str(count)] = {
                "receivers_m": [list(map(str, q)) for q in receivers],
                "velocity_z_m_s": list(map(str, values)),
                "max_gl64_to_gl88_difference": str(
                    max(abs(a - b) for a, b in zip(values, nearby, strict=True))
                ),
            }
        return {
            "finite_kernel": "regularized P0 ring and W-P1 wake sheet, radius 1/8 m",
            "decimal_digits": 75,
            "quadrature_orders": [64, 88],
            "cases": cases,
        }


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
