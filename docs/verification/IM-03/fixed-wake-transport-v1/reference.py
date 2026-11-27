# SPDX-License-Identifier: Apache-2.0
"""Independent exact stencils and spatial quadrature for free and prescribed transport."""

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


def previous_trailing_velocity(spatial, order):
    """Independent W-P1 integral on the actual parent quadrilateral after moving one corner."""
    q = (Decimal(1), Decimal(1))
    radius, bound, old = Decimal(1) / 8, Decimal(1) / 2, -Decimal(1) / 4
    pi4 = spatial["PI4"]
    rule = spatial["gauss"](order)
    velocity = bound * spatial["ring"](q, radius)
    kradius = (
        spatial["ring"](q, radius + Decimal("1e-22"))
        - spatial["ring"](q, radius - Decimal("1e-22"))
    ) / Decimal("2e-22") * bound
    bound_coefficient, old_coefficient = spatial["ring"](q, radius), Decimal(0)
    for u, wu in rule:
        for t, wt in rule:
            points = (
                (1 + 7 * u * t / 8, u, u * t, 7 * u / 8),
                (1 + (7 + t) * u / 8, u * (1 - t), u, u),
            )
            for x, y, w, jacobian in points:
                r2 = (q[0] - x) ** 2 + (q[1] - y) ** 2
                den = r2 + radius**2
                k = (r2 - 2 * radius**2) / (pi4 * den**2 * den.sqrt())
                dk = (-4 * radius - 5 * radius * (r2 - 2 * radius**2) / den) / (
                    pi4 * den**2 * den.sqrt()
                )
                weight = jacobian * wu * wt
                rho = (1 - w) * bound + w * old
                velocity += rho * k * weight
                kradius += rho * dk * weight
                bound_coefficient += (1 - w) * k * weight
                old_coefficient += w * k * weight
    phase = pi4 / 64
    sine = sin_term = phase
    cosine = cos_term = Decimal(1)
    for n in range(1, 40):
        sin_term *= -phase**2 / (2 * n * (2 * n + 1))
        cos_term *= -phase**2 / (2 * n * (2 * n - 1))
        sine += sin_term
        cosine += cos_term
    return (
        velocity + 1 + sine,
        Decimal(1) / 8
        - phase * cosine
        + bound_coefficient / 8
        + 5 * old_coefficient / 16
        + kradius / 32,
    )


def main():
    spatial = run_path(
        str(Path(__file__).resolve().parents[1] / "discrete-free-wake-velocity-v1/reference.py")
    )
    t, h = Q(5, 4), Q(1, 8)
    x = (Q(2), Q(15, 8), Q(13, 8))
    dx = (Q(1, 8), Q(1, 4), Q(1, 2))
    backward_velocity = (x[0] - x[1]) / h
    backward_direction = (dx[0] - dx[1]) / h
    assert rate(x, h) != backward_velocity
    assert rate(dx, h) != backward_direction
    prescribed_v = 1 + 2 * (t - 1)
    prescribed_dv = 2 * Q(1, 4)
    seed = (-Q(1, 8), Q(1, 4), -Q(3, 8))
    with localcontext() as context:
        context.prec = 75
        value, jvp, *_ = spatial["run"](88)[1]
        parent = previous_trailing_velocity(spatial, 88)
        coarser = previous_trailing_velocity(spatial, 64)
        assert max(abs(x - y) for x, y in zip(parent, coarser, strict=True)) < Decimal("1e-20")
        assert abs(value - parent[0]) > Decimal("0.1")

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
            "backward_euler_time_derivative": {
                "position_x_m_s": str(backward_velocity),
                "direction_x_m_s": str(backward_direction),
            },
            "free": {
                "velocity_source": (
                    "fixed geometry from discrete-free-wake-velocity-v1, second receiver, GL88"
                ),
                "residual_m_s": list(map(str, residual)),
                "jvp_m_s": list(map(str, tangent)),
                "transpose_work": str(work),
            },
            "free_birth": {
                "parent_velocity_z_m_s": str(parent[0]),
                "parent_velocity_jvp_z_m_s": str(parent[1]),
                "residual_z_m_s": str(-(value + parent[0]) / 2),
                "jvp_z_m_s": str(-(jvp + parent[1]) / 2),
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
