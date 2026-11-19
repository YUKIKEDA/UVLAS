"""Secondary floating comparison; exact acceptance remains the Fraction witness."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from collections.abc import Callable, Sequence
from fractions import Fraction
from pathlib import Path
from typing import Protocol, cast


class PchipView(Protocol):
    c: Sequence[Sequence[float]]

    def derivative(self) -> Callable[[list[float]], Sequence[float]]: ...


def exact_list(value: object) -> list[Fraction]:
    assert isinstance(value, list)
    result: list[Fraction] = []
    for element in cast("list[object]", value):
        assert isinstance(element, str)
        result.append(Fraction(element))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Refusing to overwrite output")
    data = cast("dict[str, object]", json.loads(args.reference.read_text(encoding="utf-8")))
    interpolate = importlib.import_module("scipy.interpolate")
    constructor = cast(
        "Callable[[list[float], list[float]], PchipView]", interpolate.PchipInterpolator
    )
    records = cast("list[dict[str, object]]", data["scalar_fixtures"])
    maximum = 0.0
    compared = 0
    for record in records:
        if record["name"] == "low_time_words":
            continue  # This case intentionally collapses when converted to binary64.
        times = exact_list(record["times"])
        values = exact_list(record["values"])
        expected = exact_list(record["slopes"])
        numeric_times = list(map(float, times))
        candidate = constructor(numeric_times, list(map(float, values)))
        actual = candidate.derivative()(numeric_times)
        for left, right in zip(actual, expected, strict=True):
            error = abs(float(left) - float(right)) / max(1.0, abs(float(right)))
            maximum = max(maximum, error)
            assert error <= 1e-12, record["name"]
            compared += 1
        for i, coefficients in enumerate(cast("list[object]", record["coefficients"])):
            width = float(times[i + 1] - times[i])
            for power, coefficient in enumerate(exact_list(coefficients)):
                expected_value = float(coefficient)
                value = float(candidate.c[3 - power][i]) * width**power
                error = abs(value - expected_value) / max(1.0, abs(expected_value))
                maximum = max(maximum, error)
                assert error <= 1e-12, record["name"]
                compared += 1
    scipy = importlib.import_module("scipy")
    source = importlib.import_module("scipy.interpolate._cubic")
    result = {
        "scipy_version": str(scipy.__version__),
        "scipy_cubic_sha256": hashlib.sha256(Path(str(source.__file__)).read_bytes()).hexdigest(),
        "normal_fixtures": len(records) - 1,
        "compared_values": compared,
        "maximum_scaled_error": maximum,
        "fixed_secondary_threshold": 1e-12,
        "not_product_error_budget": True,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
