# SPDX-License-Identifier: Apache-2.0
"""Independent exact stencil and inherited spatial quadrature for old-particle transport."""

import json
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path


def rate(samples, h):
    current, previous, older = samples
    return (3 * current - 4 * previous + older) / (2 * h)


def solve_position(previous, older, velocity, h):
    return (2 * h * velocity + 4 * previous - older) / 3


def main():
    spatial = run_path(
        str(Path(__file__).resolve().parents[1] / "discrete-free-wake-velocity-v1/reference.py")
    )
    t, h = Q(5, 4), Q(1, 8)
    x = (Q(2), Q(15, 8), Q(7, 4))
    dx = (Q(1, 8), Q(1, 4), Q(3, 8))
    prescribed_v = 1 + 2 * (t - 1)
    prescribed_dv = 2 * Q(1, 4)
    seed = (-Q(1, 8), Q(1, 4), -Q(3, 8))
    with localcontext() as context:
        context.prec = 75
        value, jvp, *_ = spatial["run"](88)[1]

        def decimal(x):
            return Decimal(x.numerator) / Decimal(x.denominator)

        residual = (decimal(rate(x, h)), Decimal(0), -value)
        tangent = (decimal(rate(dx, h) - Q(1, 8)), decimal(Q(1, 16)), -jvp)
        work = sum(decimal(s) * d for s, d in zip(seed, tangent, strict=True))
        return {
            "times_s": [str(t - i * h) for i in range(3)],
            "step_s": str(h),
            "x_m": list(map(str, x)),
            "dx_m": list(map(str, dx)),
            "free": {
                "velocity_source": (
                    "fixed geometry from discrete-free-wake-velocity-v1, second receiver, GL88"
                ),
                "residual_m_s": list(map(str, residual)),
                "jvp_m_s": list(map(str, tangent)),
                "transpose_work": str(work),
            },
            "prescribed": {
                "residual_x_m_s": str(rate(x, h) - prescribed_v),
                "jvp_x_m_s": str(rate(dx, h) - prescribed_dv),
                "position_x_m": str(solve_position(*x[1:], prescribed_v, h)),
                "position_jvp_x_m": str(solve_position(*dx[1:], prescribed_dv, h)),
                "time_dual": str(-2 * 2),
                "transpose_work": str(2 * (rate(dx, h) - prescribed_dv)),
            },
        }


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
