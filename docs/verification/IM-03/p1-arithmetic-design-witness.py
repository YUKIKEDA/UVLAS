"""Bounded scalar design witnesses; no product P1 implementation is imported."""
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path
import platform
import sys


def decimal_value(value: Fraction) -> Decimal:
    return Decimal(value.numerator) / Decimal(value.denominator)


def log_bounds(t: Fraction, terms: int = 128) -> tuple[Fraction, Fraction]:
    total = sum((2 * t ** (2 * j + 1) / (2 * j + 1) for j in range(terms)), Fraction())
    remainder = 2 * abs(t) ** (2 * terms + 1) / ((2 * terms + 1) * (1 - t * t))
    return total - remainder, total + remainder


def atan_bounds(x: Fraction, terms: int) -> tuple[Fraction, Fraction]:
    total = sum(((-1) ** j * x ** (2 * j + 1) / (2 * j + 1) for j in range(terms)), Fraction())
    remainder = abs(x) ** (2 * terms + 1) / (2 * terms + 1)
    return total - remainder, total + remainder


def sqrt_bracket(value: Fraction) -> tuple[Fraction, Fraction]:
    # Independent direct integer sqrt, rather than the proposed bit trials.
    integer = isqrt((value.numerator << 512) // value.denominator)
    lower = Fraction(integer, 1 << 256)
    upper = lower if lower * lower == value else Fraction(integer + 1, 1 << 256)
    assert lower * lower <= value <= upper * upper
    return lower, upper


rows: list[dict[str, object]] = []
with localcontext() as context:
    context.prec = 240
    for value in (Fraction(1), Fraction(2), Fraction(3), Fraction(1, 3), Fraction(7, 5)):
        lower, upper = sqrt_bracket(value)
        want = decimal_value(value).sqrt()
        assert decimal_value(lower) <= want <= decimal_value(upper)
        rows.append({"case": f"sqrt-{value}", "bounded": True, "exact": lower == upper})

    for t in (Fraction(1, 3), Fraction(1, 7), Fraction(1, 100), Fraction(1, 1 << 60)):
        lower, upper = log_bounds(t)
        y = (1 + t) / (1 - t)
        rounded_y = decimal_value(y)
        # Bound BOTH conversion of the exact argument and correctly rounded ln.
        # Cancellation near y=1 magnifies the former relative to ln's own ulp.
        oracle_lower = rounded_y.next_minus().ln().next_minus()
        oracle_upper = rounded_y.next_plus().ln().next_plus()
        assert decimal_value(lower) <= oracle_upper
        assert decimal_value(upper) >= oracle_lower
        rows.append({"case": f"log-transform-{t}", "overlaps_independent_oracle_interval": True})

    # Exact tan double-angle identity: atan(3/4) = 2*atan(1/3).
    # The independent direct series uses another argument and 1024 terms.
    direct_lower, direct_upper = atan_bounds(Fraction(3, 4), 1024)
    half_lower, half_upper = atan_bounds(Fraction(1, 3), 128)
    assert 2 * half_lower <= direct_lower <= direct_upper <= 2 * half_upper
    # The missing half-angle factor must be detected by this witness.
    assert half_upper < direct_lower
    rows.append({"case": "atan-half-angle-identity", "bounded": True, "missing_factor_detected": True})

    # The same identity also bounds atan(1) from the exact sqrt(2) bracket.
    lower_root, upper_root = sqrt_bracket(Fraction(2))
    z_lower = 1 / (1 + upper_root)
    z_upper = 1 / (1 + lower_root)
    half_lower, _ = atan_bounds(z_lower, 128)
    _, half_upper = atan_bounds(z_upper, 128)
    first_lower, first_upper = atan_bounds(Fraction(1, 5), 128)
    second_lower, second_upper = atan_bounds(Fraction(1, 239), 128)
    pi_lower = 16 * first_lower - 4 * second_upper
    pi_upper = 16 * first_upper - 4 * second_lower
    assert 2 * half_lower <= pi_lower / 4 <= pi_upper / 4 <= 2 * half_upper
    rows.append({"case": "atan-one-versus-machin", "bounded": True})

    # The analytic remainder bounds cover their entire stated domains.
    assert 2 * Fraction(1, 2) ** 257 / 257 < Fraction(1, 1 << 256)
    assert 2 * Fraction(1, 3) ** 257 / (257 * (1 - Fraction(1, 9))) < Fraction(1, 1 << 400)
    rows.append({"case": "uniform-series-remainders", "bounded": True})

print(json.dumps({
    "id": "p1-arithmetic-design-witness-v1",
    "scope": "scalar design identities and remainder bounds, not product P1 qualification",
    "python": sys.version.split()[0],
    "platform": platform.system(),
    "source": {"id": "p1-arithmetic-design-witness-source-v1", "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "hint": Path(__file__).name},
    "cases": rows,
    "passed": len(rows) == 12,
}, indent=2))