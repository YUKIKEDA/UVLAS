"""Limited exact design witness; not the product implementation or allocator test."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
from typing import Literal

Mode = Literal[
    "control",
    "rounded_times",
    "arithmetic",
    "no_clip",
    "no_time_scale",
    "trunc_period",
    "allow_hold_kink",
    "skip_period_slope",
]


def sign(value: F) -> int:
    return (value > 0) - (value < 0)


def endpoint(h: F, other_h: F, delta: F, other_delta: F, mode: Mode) -> F:
    value = ((2 * h + other_h) * delta - h * other_delta) / (h + other_h)
    if sign(value) != sign(delta):
        return F(0)
    if mode != "no_clip" and sign(delta) != sign(other_delta) and abs(value) > 3 * abs(delta):
        return 3 * delta
    return value


def slopes(times: tuple[F, ...], values: tuple[F, ...], mode: Mode) -> tuple[F, ...]:
    if len(times) < 2 or len(times) != len(values):
        raise ValueError("shape")
    widths = tuple(b - a for a, b in itertools.pairwise(times))
    if any(width <= 0 for width in widths):
        raise ValueError("time_order")
    secants = tuple(
        (b - a) / width for a, b, width in zip(values, values[1:], widths, strict=False)
    )
    if len(times) == 2:
        return (secants[0], secants[0])
    answer = [endpoint(widths[0], widths[1], secants[0], secants[1], mode)]
    for i in range(1, len(times) - 1):
        left, right = secants[i - 1], secants[i]
        if sign(left) == 0 or sign(left) != sign(right):
            answer.append(F(0))
        elif mode == "arithmetic":
            answer.append((left + right) / 2)
        else:
            w1, w2 = 2 * widths[i] + widths[i - 1], widths[i] + 2 * widths[i - 1]
            answer.append((w1 + w2) / (w1 / left + w2 / right))
    answer.append(endpoint(widths[-1], widths[-2], secants[-1], secants[-2], mode))
    return tuple(answer)


@dataclass(frozen=True)
class Curve:
    times: tuple[F, ...]
    values: tuple[F, ...]
    derivatives: tuple[F, ...]
    coefficients: tuple[tuple[F, F, F, F], ...]
    policy: str
    mode: Mode

    @classmethod
    def prepare(
        cls,
        times: tuple[F, ...],
        values: tuple[F, ...],
        policy: str = "reject",
        mode: Mode = "control",
        *,
        period: F | None = None,
    ) -> Curve:
        if mode == "rounded_times":
            times = tuple(F(float(time)) for time in times)
        derivatives = slopes(times, values, mode)
        if policy == "periodic":
            if period is not None and period != times[-1] - times[0]:
                raise ValueError("period_width")
            if values[0] != values[-1]:
                raise ValueError("period_value")
            if mode != "skip_period_slope" and derivatives[0] != derivatives[-1]:
                raise ValueError("period_slope")
        if policy not in {"reject", "hold", "periodic"}:
            raise ValueError("policy")
        coefficients: list[tuple[F, F, F, F]] = []
        for i, width in enumerate(b - a for a, b in itertools.pairwise(times)):
            a, b, left, right = values[i], values[i + 1], derivatives[i], derivatives[i + 1]
            coefficients.append(
                (
                    a,
                    width * left,
                    3 * (b - a) - width * (2 * left + right),
                    2 * (a - b) + width * (left + right),
                )
            )
        return cls(times, values, derivatives, tuple(coefficients), policy, mode)

    def demand_c1(self, a: F, b: F) -> None:
        if a > b:
            raise ValueError("demand_order")
        if self.policy == "reject" and (a < self.times[0] or b > self.times[-1]):
            raise ValueError("demand_range")
        if self.policy == "hold" and self.mode != "allow_hold_kink":
            for edge, slope in (
                (self.times[0], self.derivatives[0]),
                (self.times[-1], self.derivatives[-1]),
            ):
                if a <= edge <= b and slope != 0:
                    raise ValueError("hold_not_c1")

    def evaluate(self, time: F, side: str = "auto", order: int = 2) -> tuple[F, F, F]:
        if side not in {"left", "right", "auto"} or order not in {0, 1, 2}:
            raise ValueError("request")
        if self.policy == "periodic":
            period = self.times[-1] - self.times[0]
            quotient = (time - self.times[0]) / period
            turns = (
                int(quotient)
                if self.mode == "trunc_period"
                else quotient.numerator // quotient.denominator
            )
            time -= turns * period
            if time == self.times[0] and side == "left":
                time = self.times[-1]
        if time < self.times[0] or time > self.times[-1]:
            if self.policy == "hold":
                return (self.values[0] if time < self.times[0] else self.values[-1], F(0), F(0))
            raise ValueError("range")
        at_left, at_right = time == self.times[0], time == self.times[-1]
        outside_side = (at_left and side == "left") or (at_right and side == "right")
        if self.policy == "hold" and outside_side:
            return (self.values[0] if at_left else self.values[-1], F(0), F(0))
        index = len(self.times) - 2
        for i in range(len(self.times) - 1):
            if time < self.times[i + 1] or (time == self.times[i + 1] and side == "left"):
                index = i
                break
        width = self.times[index + 1] - self.times[index]
        u = (time - self.times[index]) / width
        a, b, c, d = self.coefficients[index]
        first = (3 * d * u + 2 * c) * u + b
        second = 6 * d * u + 2 * c
        if self.mode != "no_time_scale":
            first /= width
            second /= width * width
        result = (((d * u + c) * u + b) * u + a, first, second)
        if (
            side == "auto"
            and self.policy == "hold"
            and (at_left or at_right)
            and ((order >= 1 and first != 0) or (order == 2 and second != 0))
        ):
            raise ValueError("side_required")
        if (
            side == "auto"
            and order == 2
            and (time in self.times[1:-1] or (self.policy == "periodic" and at_left))
            and self.evaluate(time, "left")[2] != self.evaluate(time, "right")[2]
        ):
            raise ValueError("side_required")
        return result


def hermite(a: F, b: F, sa: F, sb: F, *, width: F, u: F) -> tuple[F, F, F]:
    # Independent endpoint data/basis, rather than the coefficient producer above.
    value = (
        (2 * u**3 - 3 * u**2 + 1) * a
        + (u**3 - 2 * u**2 + u) * width * sa
        + (-2 * u**3 + 3 * u**2) * b
        + (u**3 - u**2) * width * sb
    )
    first = (
        (6 * u**2 - 6 * u) * a
        + (3 * u**2 - 4 * u + 1) * width * sa
        + (-6 * u**2 + 6 * u) * b
        + (3 * u**2 - 2 * u) * width * sb
    ) / width
    second = (
        (12 * u - 6) * a + (6 * u - 4) * width * sa + (-12 * u + 6) * b + (6 * u - 2) * width * sb
    ) / width**2
    return (value, first, second)


def run(mode: Mode) -> dict[str, object]:
    fixtures = (
        ("two_points", (0, 2), (1, 5), (2, 2)),
        ("increasing", (0, 1, 2), (0, 1, 4), (0, F(3, 2), 4)),
        ("nonuniform", (0, 1, 3), (0, 1, 5), (F(2, 3), F(9, 7), F(8, 3))),
        ("endpoint_clip", (0, 1, 2), (0, 1, -10), (3, 0, -17)),
        ("right_clip", (0, 1, 2), (-10, 1, 0), (17, 0, -3)),
        ("flat_ends", (0, 1, 2, 3), (0, 0, 1, 1), (0, 0, 0, 0)),
        ("clip_equality", (0, 1, 2), (0, 1, -2), (3, 0, -5)),
        ("zero_secant", (0, 1, 2), (1, 1, 2), (0, 0, F(3, 2))),
        ("end_wrong_sign", (0, 1, 2), (0, 1, 10), (0, F(9, 5), 13)),
        ("constant", (0, 1, 3), (4, 4, 4), (0, 0, 0)),
        ("negative", (0, 1, 3), (0, -1, -5), (-F(2, 3), -F(9, 7), -F(8, 3))),
        (
            "low_time_words",
            (F(1), F(1) + F(1, 2**200), F(1) + F(1, 2**199)),
            (0, 1, 2),
            (2**200, 2**200, 2**200),
        ),
    )
    checks = 0
    encoded: list[dict[str, object]] = []
    for name, source_times, source_values, source_slopes in fixtures:
        times, values, expected = (
            tuple(map(F, source_times)),
            tuple(map(F, source_values)),
            tuple(map(F, source_slopes)),
        )
        curve = Curve.prepare(times, values, mode=mode)
        assert curve.derivatives == expected, name
        checks += 1
        for index, coefficients in enumerate(curve.coefficients):
            width = times[index + 1] - times[index]
            for u in (F(0), F(1, 4), F(1, 2), F(3, 4), F(1)):
                query = times[index] + u * width
                side = "left" if u == 1 else "right"
                assert curve.evaluate(query, side) == hermite(
                    values[index],
                    values[index + 1],
                    expected[index],
                    expected[index + 1],
                    width=width,
                    u=u,
                ), (name, index, u)
                checks += 1
            integral = width * sum(
                (coefficient / (power + 1) for power, coefficient in enumerate(coefficients)), F(0)
            )
            assert integral == width * (
                (values[index] + values[index + 1]) / 2
                + width * (expected[index] - expected[index + 1]) / 12
            ), name
            checks += 1
        encoded.append(
            {
                "name": name,
                "times": list(map(str, times)),
                "values": list(map(str, values)),
                "slopes": list(map(str, expected)),
                "coefficients": [[str(value) for value in row] for row in curve.coefficients],
            }
        )
    held = Curve.prepare((F(0), F(2)), (F(1), F(5)), "hold", mode)
    assert held.evaluate(F(-1), order=1) == (F(1), F(0), F(0))
    assert held.evaluate(F(0), "left") == (F(1), F(0), F(0))
    assert held.evaluate(F(0), "right")[1] == 2
    held.demand_c1(F(1, 2), F(3, 2))
    smooth_hold = Curve.prepare(
        tuple(map(F, (0, 1, 2, 3))), tuple(map(F, (0, 0, 1, 1))), "hold", mode
    )
    smooth_hold.demand_c1(F(-1), F(4))
    assert smooth_hold.evaluate(F(0)) == (F(0), F(0), F(0))
    rejected: list[str] = []
    negatives: tuple[tuple[str, Callable[[], object]], ...] = (
        ("short", lambda: Curve.prepare((F(0),), (F(1),), mode=mode)),
        (
            "period_width",
            lambda: Curve.prepare(
                tuple(map(F, (0, 1, 2, 3))),
                tuple(map(F, (0, 1, -1, 0))),
                "periodic",
                mode,
                period=F(2),
            ),
        ),
        ("duplicate", lambda: Curve.prepare((F(0), F(0)), (F(1), F(2)), mode=mode)),
        ("reversed", lambda: Curve.prepare((F(1), F(0)), (F(1), F(2)), mode=mode)),
        (
            "reject_range",
            lambda: Curve.prepare((F(0), F(2)), (F(1), F(5)), mode=mode).evaluate(F(-1)),
        ),
        ("hold_kink", lambda: held.demand_c1(F(-1), F(1))),
        ("hold_side", lambda: held.evaluate(F(0), order=1)),
        (
            "period_value",
            lambda: Curve.prepare((F(0), F(1), F(2)), (F(0), F(1), F(2)), "periodic", mode),
        ),
        (
            "period_slope",
            lambda: Curve.prepare((F(0), F(1), F(2)), (F(0), F(1), F(0)), "periodic", mode),
        ),
    )
    for name, call in negatives:
        try:
            call()
        except ValueError:
            rejected.append(name)
        else:
            raise AssertionError(name)
    periodic = Curve.prepare(
        tuple(map(F, (0, 1, 2, 3))), tuple(map(F, (0, 1, -1, 0))), "periodic", mode
    )
    for query in (F(-1, 2), F(-13, 4), F(25, 4), F(3 * 2**200) + F(1, 4)):
        reduced = query % 3
        assert periodic.evaluate(query, "right") == periodic.evaluate(reduced, "right"), query
        checks += 1
    assert periodic.evaluate(F(0), "left", 1)[:2] == periodic.evaluate(F(0), "right", 1)[:2]
    try:
        periodic.evaluate(F(0), order=2)
    except ValueError:
        rejected.append("period_second_side")
    else:
        raise AssertionError("period_second_side")
    # Positive time scaling and a common origin must preserve the physical chain rule.
    base = Curve.prepare(tuple(map(F, (0, 1, 3))), tuple(map(F, (0, 1, 5))), mode=mode)
    for scale in (F(1, 2**200), F(2**200)):
        origin = F(2**900)
        changed = Curve.prepare(
            tuple(origin + scale * time for time in base.times), base.values, mode=mode
        )
        for query in (F(1, 4), F(3, 2)):
            value = changed.evaluate(origin + scale * query, "right")
            reference = base.evaluate(query, "right")
            assert value == (reference[0], reference[1] / scale, reference[2] / scale**2)
            checks += 1
    assert float(F(1) + F(1, 2**200)) == float(F(1)), (
        "binary64 mutation must collapse the separate knot"
    )
    checks += 1
    gcd_cases = 0
    for a in range(1, 128):
        for b in range(1, 128):
            x, y, count = a, b, 0
            while y:
                x, y = y, x % y
                count += 1
            assert count <= 2 * max(a.bit_length(), b.bit_length()) + 1
            gcd_cases += 1
    return {
        "scalar_fixtures": encoded,
        "exact_equalities": checks,
        "rejected": rejected,
        "gcd_bound_cases": gcd_cases,
        "claim": "limited design witness only",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Refusing to overwrite output")
    result = run("control")
    mutations: list[dict[str, object]] = []
    for mode in (
        "rounded_times",
        "arithmetic",
        "no_clip",
        "no_time_scale",
        "trunc_period",
        "allow_hold_kink",
        "skip_period_slope",
    ):
        try:
            run(mode)
        except (AssertionError, ValueError) as error:
            mutations.append({"mode": mode, "detected": True, "reason": str(error)})
        else:
            raise AssertionError(f"survived mutation: {mode}")
    result["mutations"] = mutations
    result["producer_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    print(
        json.dumps(
            {
                key: result[key]
                for key in ("exact_equalities", "rejected", "gcd_bound_cases", "mutations")
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
