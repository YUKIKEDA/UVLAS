"""Import frozen scalar PCHIP witnesses and solve Hermite conditions independently."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
from typing import cast

SOURCE_SHA = "497287324f6084d162a41aa98bfad9993adb91287fb78e6e5467adef9a2e22c5"


def rational_list(value: object) -> list[F]:
    if not isinstance(value, list):
        raise ValueError("expected exact rational list")
    result: list[F] = []
    for item in cast("list[object]", value):
        if not isinstance(item, str):
            raise ValueError("exact rational must be a string")
        result.append(F(item))
    return result


def hermite(y0: F, y1: F, h: F, s0: F, s1: F) -> list[F]:
    # Independent row elimination of endpoint value/derivative conditions in u.
    rows = [[F(x) for x in row] for row in ([1, 0, 0, 0], [1, 1, 1, 1], [0, 1, 0, 0], [0, 1, 2, 3])]
    rhs = [y0, y1, h * s0, h * s1]
    for col in range(4):
        pivot = next(i for i in range(col, 4) if rows[i][col])
        rows[pivot], rows[col] = rows[col], rows[pivot]
        rhs[pivot], rhs[col] = rhs[col], rhs[pivot]
        divisor = rows[col][col]
        rows[col] = [v / divisor for v in rows[col]]
        rhs[col] /= divisor
        for row in range(4):
            if row == col:
                continue
            factor = rows[row][col]
            rows[row] = [a - factor * b for a, b in zip(rows[row], rows[col], strict=True)]
            rhs[row] -= factor * rhs[col]
    return rhs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Refusing to overwrite fixture")
    raw = args.source.read_bytes()
    if hashlib.sha256(raw).hexdigest() != SOURCE_SHA:
        raise SystemExit("Frozen source hash mismatch")
    source = cast("dict[str, object]", json.loads(raw))
    records = cast("list[dict[str, object]]", source["scalar_fixtures"])
    cases: list[dict[str, object]] = []
    checks = 0
    for record in records:
        name = record["name"]
        if not isinstance(name, str):
            raise ValueError("case name")
        times = rational_list(record["times"])
        values = rational_list(record["values"])
        slopes = rational_list(record["slopes"])
        original = cast("list[object]", record["coefficients"])
        if len(times) < 2 or len(values) != len(times) or len(slopes) != len(times):
            raise ValueError("source shape")
        for i in range(len(times) - 1):
            solved = hermite(
                values[i], values[i + 1], times[i + 1] - times[i], slopes[i], slopes[i + 1]
            )
            if solved != rational_list(original[i]):
                raise ValueError("source coefficients violate independent Hermite equations")
            checks += 4
        transforms = [
            ("base", F(1), F(0), F(1), F(0)),
            ("origin", F(1), F(2**300), F(1), F(0)),
            ("time_scale", F(1, 2**100), F(0), F(1), F(0)),
            ("value_scale", F(1), F(0), F(-3), F(5)),
        ]
        for label, ts, to, ys, yo in transforms:
            t = [v * ts + to for v in times]
            y = [v * ys + yo for v in values]
            s = [v * ys / ts for v in slopes]
            coefficients = [
                hermite(y[i], y[i + 1], t[i + 1] - t[i], s[i], s[i + 1]) for i in range(len(t) - 1)
            ]
            cases.append(
                {
                    "name": name + "/" + label,
                    "times": list(map(str, t)),
                    "values": list(map(str, y)),
                    "components": 1,
                    "slopes": list(map(str, s)),
                    "coefficients": [[str(v) for v in row] for row in coefficients],
                }
            )
        # Three mathematically ordered components: original, affine image, exact zero.
        y3 = [v for y in values for v in (y, 2 * y + 1, F(0))]
        s3 = [v for s in slopes for v in (s, 2 * s, F(0))]
        c3: list[list[str]] = []
        for i in range(len(times) - 1):
            for component in range(3):
                c3.append(
                    list(
                        map(
                            str,
                            hermite(
                                y3[3 * i + component],
                                y3[3 * (i + 1) + component],
                                times[i + 1] - times[i],
                                s3[3 * i + component],
                                s3[3 * (i + 1) + component],
                            ),
                        )
                    )
                )
        cases.append(
            {
                "name": name + "/components3",
                "times": list(map(str, times)),
                "values": list(map(str, y3)),
                "components": 3,
                "slopes": list(map(str, s3)),
                "coefficients": c3,
            }
        )
    output = {
        "schema": "uvlas.pchip-coefficient-reference/1",
        "source_sha256": SOURCE_SHA,
        "independent_hermite_scalar_checks": checks,
        "cases": cases,
    }
    encoded = (json.dumps(output, ensure_ascii=False, sort_keys=True, indent=2) + "\n").encode()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("xb") as stream:
        stream.write(encoded)
    print(
        json.dumps(
            {
                "cases": len(cases),
                "sha256": hashlib.sha256(encoded).hexdigest(),
                "hermite_checks": checks,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
