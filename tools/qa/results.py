"""Baseline results and result change control (QA-03, QA-04, QA-07).

    python tools/qa/results.py check  --case numeric-canary --run <executable>
    python tools/qa/results.py accept --case numeric-canary --run <executable> [--record NNNN]
    python tools/qa/results.py platforms

A case program prints "platform <key>" and then one "name value" line per
quantity, the value in hexadecimal floating point. `check` compares that output
bit for bit with qa/baselines/<case>/<platform>.json and fails when the
platform has no baseline or any value differs.

`accept` writes a baseline, and only in one of these situations:

  * the case has no baseline yet on any platform (a new case);
  * the platform has no baseline and its values equal those of a qualified
    platform bit for bit (the new platform is recorded as identical to it);
  * --record names a result change record in docs/quality/result-changes/
    that mentions this case.

There is no other way to change a baseline. Standard library only (Python 3.10+).
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
BASELINES = REPO / "qa" / "baselines"
RECORDS = REPO / "docs" / "quality" / "result-changes"
# std::format("{:a}") writes hexadecimal floating point without the 0x prefix.
HEX = re.compile(r"^-?(?:(?:0x)?[0-9a-f]+(?:\.[0-9a-f]+)?p[+-]\d+|inf|nan)$")


class ResultError(Exception):
    pass


def run_case(executable: Path) -> tuple[str, dict[str, str]]:
    done = subprocess.run(
        [str(executable)], capture_output=True, text=True, timeout=600, check=False
    )
    if done.returncode != 0:
        raise ResultError(f"{executable} exited with {done.returncode}: {done.stderr.strip()}")
    platform = ""
    values: dict[str, str] = {}
    for line in done.stdout.splitlines():
        name, _, value = line.strip().partition(" ")
        if not name:
            continue
        if name == "platform":
            platform = value
        elif not HEX.match(value):
            raise ResultError(f"{name}: {value!r} is not a hexadecimal floating-point value")
        elif name in values:
            raise ResultError(f"{name} is printed twice")
        else:
            values[name] = value
    if not platform or not values:
        raise ResultError("the case printed no platform key or no values")
    return platform, values


def baseline_path(case: str, platform: str) -> Path:
    return BASELINES / case / f"{platform}.json"


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def differences(expected: dict[str, str], actual: dict[str, str]) -> list[str]:
    lines: list[str] = []
    for name in sorted(expected.keys() | actual.keys()):
        old, new = expected.get(name), actual.get(name)
        if old == new:
            continue
        if old is None:
            lines.append(f"  {name}: new quantity {new}")
        elif new is None:
            lines.append(f"  {name}: quantity removed (was {old})")
        else:
            change = ""
            try:
                a, b = float.fromhex(old), float.fromhex(new)
                if a != 0.0:
                    change = f"  relative change {abs(b - a) / abs(a):.3e}"
            except ValueError:
                pass
            lines.append(f"  {name}: {old} -> {new}{change}")
    return lines


def check(case: str, executable: Path) -> int:
    platform, values = run_case(executable)
    path = baseline_path(case, platform)
    if not path.is_file():
        known = (
            sorted(p.stem for p in (BASELINES / case).glob("*.json"))
            if (BASELINES / case).is_dir()
            else []
        )
        print(f"{case}: platform {platform} is not qualified (QA-07).")
        print(f"  qualified platforms: {', '.join(known) or 'none'}")
        print(
            "  qualify it with `results.py accept`; if its values differ from every "
            "qualified platform, that needs a result change record."
        )
        return 1
    found = differences(load(path)["values"], values)
    if found:
        print(
            f"{case} on {platform}: {len(found)} of {len(values)} values differ "
            "from the baseline (QA-04)."
        )
        print("\n".join(found))
        print(
            "  A change of results needs a record in docs/quality/result-changes/ in the same "
            "pull request; then `results.py accept --record NNNN`."
        )
        return 1
    print(f"{case} on {platform}: {len(values)} values equal the baseline bit for bit")
    return 0


def accept(case: str, executable: Path, record: str | None) -> int:
    platform, values = run_case(executable)
    path = baseline_path(case, platform)
    others: dict[str, dict[str, Any]] = (
        {p.stem: load(p) for p in sorted((BASELINES / case).glob("*.json")) if p != path}
        if (BASELINES / case).is_dir()
        else {}
    )
    entry: dict[str, Any] = {"schema": "uvlas.baseline/1", "case": case, "platform": platform}
    if record is not None:
        matches = sorted(RECORDS.glob(f"{record}-*.md")) if RECORDS.is_dir() else []
        if len(matches) != 1:
            raise ResultError(f"no single result change record {record}-*.md in {RECORDS}")
        if case not in matches[0].read_text(encoding="utf-8"):
            raise ResultError(f"{matches[0].name} does not mention the case {case}")
        entry["result_change_record"] = matches[0].name
    elif path.is_file():
        if load(path)["values"] == values:
            print(f"{case} on {platform}: already equal to the baseline; nothing written")
            return 0
        raise ResultError(
            "the values differ from the baseline: a result change record is "
            "required (--record NNNN)"
        )
    elif others:
        same = [name for name, data in others.items() if data["values"] == values]
        if not same:
            raise ResultError(
                f"platform {platform} differs from every qualified platform "
                f"({', '.join(others)}): a result change record is required"
            )
        entry["identical_to"] = same[0]
    else:
        entry["note"] = "first baseline of a new case"
    entry["values"] = dict(sorted(values.items()))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(entry, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {path.relative_to(REPO).as_posix()}: {len(values)} values")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("mode", choices=("check", "accept", "platforms"))
    parser.add_argument("--case")
    parser.add_argument("--run", type=Path, help="The case program built from this tree")
    parser.add_argument("--record", help="Number of the result change record, for example 0001")
    args = parser.parse_args(argv)
    try:
        if args.mode == "platforms":
            for case in (
                sorted(p for p in BASELINES.iterdir() if p.is_dir()) if BASELINES.is_dir() else []
            ):
                print(f"{case.name}: " + ", ".join(sorted(p.stem for p in case.glob("*.json"))))
            return 0
        if not args.case or not args.run:
            parser.error("--case and --run are required")
        if args.mode == "check":
            return check(args.case, args.run)
        return accept(args.case, args.run, args.record)
    except (ResultError, OSError, ValueError, subprocess.TimeoutExpired) as error:
        print(f"results: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
