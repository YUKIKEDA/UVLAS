"""Evaluate exact interval witnesses with independent cubic Hermite basis functions."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
from typing import cast

SOURCE_PATH = "tests/aero/fixtures/pchip-coefficients-v1.json"
SOURCE_SHA = "c8b8afb2b25acd691eef0060b8a38afe6dedf56eb9e6bbe8e10cd877c903c1e8"
PRODUCER_PATH = "docs/verification/IM-03/pchip-interval-evaluation-v1/reference.py"
POINTS = (F(0), F(1, 3), F(1))


def rationals(value: object) -> list[F]:
    if not isinstance(value, list):
        raise ValueError("Expected an exact rational array")
    result: list[F] = []
    for item in cast("list[object]", value):
        if not isinstance(item, str):
            raise ValueError("Exact rational values must be strings")
        result.append(F(item))
    return result


def hermite(y0: F, y1: F, m0: F, m1: F, *, h: F, u: F) -> tuple[F, F, F]:
    """Use endpoint data and physical-time slopes, without coefficient construction."""
    if h <= 0:
        raise ValueError("Interval length must be positive")
    # The value basis is factored. Derivatives below are analytic basis derivatives,
    # with dt/du=h; no power-basis or product coefficient is read or generated.
    h00 = (1 + 2 * u) * (1 - u) ** 2
    h10 = u * (1 - u) ** 2
    h01 = u**2 * (3 - 2 * u)
    h11 = u**2 * (u - 1)
    value = h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1
    dt = 6 * u * (u - 1) * y0 / h + (1 - u) * (1 - 3 * u) * m0
    dt += 6 * u * (1 - u) * y1 / h + u * (3 * u - 2) * m1
    dtt = (12 * u - 6) * y0 / h**2 + (6 * u - 4) * m0 / h
    dtt += (6 - 12 * u) * y1 / h**2 + (6 * u - 2) * m1 / h
    return value, dt, dtt


def closed_form_checks() -> int:
    """Compare against direct line/cubic/constant formulas, not another Hermite evaluator."""
    count = 0
    for u in POINTS:
        t = 2 * u
        witnesses = (
            (hermite(F(1), F(5), F(2), F(2), h=F(2), u=u), (1 + 2 * t, F(2), F(0))),
            (hermite(F(0), F(8), F(0), F(12), h=F(2), u=u), (t**3, 3 * t**2, 6 * t)),
            (hermite(F(-7), F(-7), F(0), F(0), h=F(2), u=u), (F(-7), F(0), F(0))),
        )
        for actual, expected in witnesses:
            if actual != expected:
                raise ValueError("Independent closed-form identity failed")
            count += 1
    return count


def generate(raw: bytes) -> dict[str, object]:
    if hashlib.sha256(raw).hexdigest() != SOURCE_SHA:
        raise ValueError("Frozen coefficient fixture hash mismatch")
    source = cast("dict[str, object]", json.loads(raw))
    if source.get("schema") != "uvlas.pchip-coefficient-reference/1":
        raise ValueError("Source schema mismatch")
    records = cast("list[dict[str, object]]", source["cases"])
    if len(records) != 60:
        raise ValueError("Expected all sixty frozen sources")
    names: set[str] = set()
    cases: list[dict[str, object]] = []
    for record in records:
        name, components = record["name"], record["components"]
        if not isinstance(name, str) or not name or name in names:
            raise ValueError("Source names must be unique nonempty strings")
        if isinstance(components, bool) or not isinstance(components, int) or components < 1:
            raise ValueError("Source component count must be a positive integer")
        names.add(name)
        times = rationals(record["times"])
        values = rationals(record["values"])
        slopes = rationals(record["slopes"])
        # Flattened source arrays are point-major, then component-major.
        if len(times) < 2 or len(values) != len(times) * components or len(slopes) != len(values):
            raise ValueError("Source endpoint shape mismatch")
        for interval in range(len(times) - 1):
            h = times[interval + 1] - times[interval]
            if h <= 0:
                raise ValueError("Source times must be strictly increasing")
            for u in POINTS:
                samples: list[dict[str, str]] = []
                for component in range(components):
                    left = interval * components + component
                    right = left + components
                    value, dt, dtt = hermite(
                        values[left], values[right], slopes[left], slopes[right], h=h, u=u
                    )
                    # Endpoint identities are checked separately on the original input data.
                    if (u == 0 and (value != values[left] or dt != slopes[left])) or (
                        u == 1 and (value != values[right] or dt != slopes[right])
                    ):
                        raise ValueError("Endpoint Hermite identity failed")
                    samples.append({"value": str(value), "dt": str(dt), "dtt": str(dtt)})
                cases.append(
                    {
                        "source_name": name,
                        "interval": interval,
                        "time": str(times[interval] + h * u),
                        "u": str(u),
                        "values": samples,
                    }
                )
    return {
        "schema": "uvlas.pchip-interval-reference/1",
        "source_path": SOURCE_PATH,
        "source_sha256": SOURCE_SHA,
        "producer_path": PRODUCER_PATH,
        "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "method": "Fraction endpoint data and independent factored cubic Hermite basis",
        "independent_closed_form_checks": closed_form_checks(),
        "source_count": len(names),
        "cases": cases,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source_path = cast("Path", args.source)
    output_path = cast("Path", args.output)
    if output_path.exists():
        raise SystemExit("Refusing to overwrite an existing fixture")
    result = generate(source_path.read_bytes())
    encoded = (json.dumps(result, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    try:
        with output_path.open("xb") as stream:
            stream.write(encoded)
    except FileExistsError as error:
        raise SystemExit("Refusing to overwrite an existing fixture") from error
    print(json.dumps({"output_sha256": hashlib.sha256(encoded).hexdigest(), "bytes": len(encoded)}))


if __name__ == "__main__":
    main()
