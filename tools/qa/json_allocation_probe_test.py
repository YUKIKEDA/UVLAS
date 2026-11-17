"""Check fixed JSON SAX values and isolate every observed C++ allocation failure."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import cast


@dataclass(frozen=True)
class Expected:
    name: str
    accepted: bool
    events: bytes
    count: int
    error_id: int = 0
    error_byte: int = 0


# Hand-written event literals: no product parser generates the expected values.
# Encoding includes signed/unsigned kinds, float bits/raw spelling, byte lengths,
# embedded NUL and decoded UTF-8. Object order is deliberately preserved.
EXPECTED = (
    Expected(
        "all-kinds",
        True,
        b"O;K1:a;A;N;T;F;I-1;U0;U9223372036854775807;U18446744073709551615;"
        b"D3ff8000000000000:3:1.5;S3:a\x00b;Z;K1:b;O;K1:k;S9:" + "日本語".encode() + b";C;C;",
        19,
    ),
    Expected("surrogate", False, b"A;", 1, 101, 9),
    Expected("range", False, b"A;", 1, 406, 6),
    Expected("deep", True, b"A;" * 80 + b"U0;" + b"Z;" * 80, 161),
    Expected("long-token", False, b"A;S600:" + b"x" * 600 + b";Z;", 3, 101, 608),
)


def invoke(
    probe: Path, name: str, timeout: float, fault: int | None = None
) -> subprocess.CompletedProcess[str]:
    command = [str(probe), name]
    if fault is not None:
        command.append(str(fault))
    return subprocess.run(
        command, capture_output=True, text=True, encoding="utf-8", timeout=timeout, check=False
    )


def control(probe: Path, expected: Expected, timeout: float) -> tuple[int, bool]:
    result = invoke(probe, expected.name, timeout)
    if result.returncode != 0:
        raise ValueError(f"{expected.name}: control exit {result.returncode}: {result.stderr}")
    value: object = json.loads(result.stdout)
    if not isinstance(value, dict):
        raise ValueError(f"{expected.name}: control response is not an object")
    payload = cast("dict[str, object]", value)
    required: dict[str, object] = {
        "schema": 1,
        "accepted": expected.accepted,
        "events": expected.count,
        "events_hex": expected.events.hex(),
        "error_id": expected.error_id,
        "error_byte": expected.error_byte,
        "error_prefix": expected.error_id != 0,
        "complete": True,
    }
    for key, wanted in required.items():
        actual = payload.get(key)
        if type(actual) is not type(wanted) or actual != wanted:
            raise ValueError(f"{expected.name}: {key}: expected {wanted!r}, got {actual!r}")
    count = payload.get("new_count")
    debug = payload.get("debug_proxies")
    if type(count) is not int or count < 0 or type(debug) is not bool:
        raise ValueError(f"{expected.name}: invalid allocation count or configuration")
    return count, debug


def require_faults(probe: Path, name: str, count: int, exit_code: int, timeout: float) -> None:
    for fault in range(count):
        result = invoke(probe, name, timeout, fault)
        if result.returncode != exit_code:
            raise ValueError(
                f"{name}: fault {fault}/{count}: expected exit {exit_code}, "
                f"got {result.returncode}: {result.stdout} {result.stderr}"
            )


def run(probe: Path, timeout: float) -> dict[str, object]:
    observations: list[dict[str, object]] = []
    configuration: bool | None = None
    total = 0
    for expected in EXPECTED:
        count, debug = control(probe, expected, timeout)
        if count == 0:
            raise ValueError(f"{expected.name}: zero allocation sites cannot verify propagation")
        if configuration is not None and configuration != debug:
            raise ValueError("Configuration changed between child processes")
        configuration = debug
        require_faults(probe, expected.name, count, 42, timeout)
        total += count
        observations.append({"case": expected.name, "new_count": count, "caught": count})
    constructor_observations: list[dict[str, object]] = []
    for name in ("default-vector", "explicit-bool", "explicit-char", "explicit-string"):
        count, debug = control(probe, Expected(name, True, b"", 0), timeout)
        if debug != configuration:
            raise ValueError("Configuration changed in constructor controls")
        if debug and count == 0:
            raise ValueError(f"{name}: MSVC Debug proxy negative control observed no allocation")
        exit_code = 77 if debug and name == "default-vector" else 42
        require_faults(probe, name, count, exit_code, timeout)
        constructor_observations.append(
            {"case": name, "new_count": count, "fault_exit": exit_code if count else None}
        )
    return {
        "schema": "uvlas.json-allocation-probe/1",
        "debug_proxies": configuration,
        "sax_cases": observations,
        "sax_faults_caught": total,
        "constructor_controls": constructor_observations,
        "scope": "Selected JSON SAX C++new sites; no CRT malloc injection or all-format guarantee",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=float, default=10.0)
    arguments = parser.parse_args()
    probe = cast("Path", arguments.probe).resolve()
    timeout = cast("float", arguments.timeout_seconds)
    if not math.isfinite(timeout) or timeout <= 0 or timeout > 60:
        parser.error("--timeout-seconds must be finite and in (0, 60]")
    try:
        summary = run(probe, timeout)
    except (ValueError, OSError, subprocess.TimeoutExpired) as error:
        print(f"JSON allocation regression failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
