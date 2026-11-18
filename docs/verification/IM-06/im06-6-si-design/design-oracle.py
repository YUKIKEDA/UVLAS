"""Check SI design mathematics without importing UVLAS product code."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from pathlib import Path


def pair(value: F) -> list[str]:
    return [str(value.numerator), str(value.denominator)]


def atan_bounds(inverse: int) -> tuple[F, F]:
    total = sum(
        (F((-1) ** j, (2 * j + 1) * inverse ** (2 * j + 1)) for j in range(128)),
        start=F(0),
    )
    return total, total + F(1, 257 * inverse**257)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    low5, high5 = atan_bounds(5)
    low239, high239 = atan_bounds(239)
    lo, hi = 16 * low5 - 4 * high239, 16 * high5 - 4 * low239
    assert F("3.14159265358979323846264338327950288419716939937510") < lo
    assert hi < F("3.14159265358979323846264338327950288419716939937511")
    assert hi - lo < F(1, 2**580)
    # tan(4 atan(1/5) - atan(1/239))=1, with angle in (0, pi/2).
    twice = 2 * F(1, 5) / (1 - F(1, 5) ** 2)
    four = 2 * twice / (1 - twice**2)
    assert four == F(120, 119)
    assert (four - F(1, 239)) / (1 + four * F(1, 239)) == 1
    units = {
        "m": F(1), "mm": F(1, 1000), "cm": F(1, 100), "km": F(1000),
        "in": F(127, 5000), "ft": F(381, 1250),
        "kg": F(1), "g": F(1, 1000), "t": F(1000),
        "lbm": F(45359237, 100000000), "s": F(1), "ms": F(1, 1000), "min": F(60),
    }
    cases = [
        ("millimetre", "1000", ("mm", "kg", "s"), (1, 0, 0), "1"),
        ("centimetre", "100", ("cm", "kg", "s"), (1, 0, 0), "1"),
        ("kilometre", "1", ("km", "kg", "s"), (1, 0, 0), "1000"),
        ("inch", "100", ("in", "kg", "s"), (1, 0, 0), "2.54"),
        ("foot", "10", ("ft", "kg", "s"), (1, 0, 0), "3.048"),
        ("gram", "1000", ("m", "g", "s"), (0, 1, 0), "1"),
        ("tonne", "1", ("m", "t", "s"), (0, 1, 0), "1000"),
        ("pound", "2", ("m", "lbm", "s"), (0, 1, 0), "0.90718474"),
        ("millisecond", "1000", ("m", "kg", "ms"), (0, 0, 1), "1"),
        ("minute", "2", ("m", "kg", "min"), (0, 0, 1), "120"),
        ("velocity", "2", ("mm", "kg", "ms"), (1, 0, -1), "2"),
        ("density", "1", ("mm", "g", "s"), (-3, 1, 0), "1000000"),
        ("pressure", "1", ("mm", "g", "ms"), (-1, 1, -2), "1000000"),
        ("beam_curvature", "1", ("mm", "kg", "s"), (-1, 0, 0), "1000"),
        ("beam_axial_K", "1", ("mm", "g", "ms"), (1, 1, -2), "1"),
        ("beam_bending_K", "1", ("mm", "g", "ms"), (3, 1, -2), "0.000001"),
        ("beam_bending_D", "1", ("mm", "g", "ms"), (3, 1, -1), "0.000000001"),
        ("exact_cancellation", "0.1", ("km", "kg", "s"), (1, 0, 0), "100"),
        ("exact_half", "500", ("mm", "kg", "s"), (1, 0, 0), "0.5"),
        ("non_dyadic", "1.20", ("in", "kg", "s"), (1, 0, 0), "0.03048"),
    ]
    conversions = []
    for name, raw, roots, dimensions, expected in cases:
        factor = F(1)
        for root, power in zip(roots, dimensions, strict=True):
            factor *= units[root] ** power
        actual = F(raw) * factor
        assert actual == F(expected), name
        dyadic = actual.denominator & (actual.denominator - 1) == 0
        conversions.append({"name": name, "actual": pair(actual), "dyadic": dyadic})
    assert conversions[-3]["dyadic"] and conversions[-2]["dyadic"]
    assert not conversions[-1]["dyadic"]
    angles = []
    for name, raw, power, coefficient, pi_power in [
        ("180deg", F(180), 1, F(1), 1),
        ("360deg", F(360), 1, F(2), 1),
        ("zero", F(0), 1, F(0), 0),
        ("support_rotational_K", F(1), -2, F(32400), -2),
    ]:
        assert raw * F(1, 180) ** power == coefficient
        assert (power if raw else 0) == pi_power
        angles.append({"name": name, "coefficient": pair(coefficient), "pi_power": pi_power})
    assert F(1, 180) ** 2 * F(180) ** 2 == 1  # Rotation work invariance.
    budgets = []
    for name, y, scale, absolute, relative, error, expected in [
        ("equal", F(7), F(2), F(3), F(0), F(6), True),
        ("above", F(7), F(2), F(3), F(0), F(6) + F(1, 2**100), False),
        ("relative_negative", F(-7), F(2), F(1), F(1, 7), F(3), True),
        ("zero", F(0), F(1), F(0), F(0), F(0), True),
        ("zero_not_exact", F(0), F(1), F(0), F(0), F(1, 10**400), False),
        ("tiny", F(1, 10**400), F(1, 10**400), F(1, 10), F(0), F(1, 10**401), True),
        ("huge", F(10**400), F(10**400), F(1, 10), F(0), F(10**399), True),
    ]:
        accepted = error <= absolute * scale + relative * abs(y)
        assert accepted == expected
        budgets.append({"name": name, "accepted": accepted})
    first, second = F(1), F(1) + F(1, 10**30)
    assert first < second and float(first) == float(second)
    assert first / 1000 < second / 1000
    assert F("1e-400") != 0 and float(F("1e-400")) == 0
    assert F(1, 10) <= F(3, 20) < 2 * F(1, 10)
    assert not (F(3, 2) <= 1) and F(3, 2) <= 2
    assert not (lo <= F("3.14") <= hi)
    assert F(1000) * units["mm"] == 1 != F(1000) * units["mm"] ** 2
    assert len({("source-A", 7), ("source-A", 7)}) == 1
    assert len({("source-A", 7), ("source-A", 8)}) == 2
    result = {
        "schema": "uvlas.si-design-oracle/1",
        "scope": "Design mathematics only; no product API, interval kernel or QoI qualification.",
        "unit_controls": conversions, "angle_recipes": angles, "budget_boundaries": budgets,
        "pi": {"terms_each": 128, "width_less_than_2_pow_minus_580": True,
               "lower": pair(lo), "upper": pair(hi)},
        "exact_order_and_ledger_checks": "passed",
        "negative_witnesses": [
            "float64 preparse collapses 1 and 1+1e-30",
            "float64 preparse flushes nonzero 1e-400",
            "relative RHS upper end incorrectly admits E=1.5 against [1,2]",
            "double charging E=0.1 violates otherwise satisfied budget 0.15",
            "rational pi surrogate 3.14 lies outside the proved enclosure",
            "repeated millimetre conversion turns 1 m into 0.001 m",
        ],
        "product_mutations_executed": False,
        "product_tests_executed": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("20 unit controls, 4 angle recipes, 7 budget boundaries; Machin/order/ledger: passed")


if __name__ == "__main__":
    main()
