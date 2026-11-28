# SPDX-License-Identifier: Apache-2.0
"""Regenerate pinned positive 1e-7 m V-MULTI companion fixtures."""

import argparse
import hashlib
import json
import runpy
from pathlib import Path

import numpy as np

ORACLE_SHA = "86a5226879f6a39e93dbc285d1ece4ba80023b37ad453f120ea0061f62f6c332"
COMPANION_SHA = "a0d44fd777be6a42909bdad846b068e0766930fff523dec95adf0eb4e211346f"
RADIUS = 1e-7


def checked_bytes(path, expected):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError(f"Pinned input differs: {path.name}")
    return raw


def build(oracle, companion, key):
    original = next(row for row in companion["records"] if tuple(row["case"]) == key)
    vertices, points, normals, rhs, quads = oracle["received_context"](*key)
    age = oracle["age_partition"](points, vertices)
    panels = oracle["subdivide"](age, 4)
    if original["panel_starts_seconds"] != panels[:-1].tolist():
        raise ValueError(f"Pinned age panels differ for {key}")
    expected_rhs = np.asarray(original["rhs"])
    if not (np.array_equal(rhs.real, expected_rhs[..., 0])
            and np.array_equal(rhs.imag, expected_rhs[..., 1])):
        raise ValueError(f"Pinned driving rows differ for {key}")

    initial = oracle["ring_projection"](points, normals, quads, RADIUS)
    finite, _ = oracle["evaluate"](initial, vertices, points, normals, rhs, RADIUS, panels, 48)
    row = {
        "case": key,
        **{name: original[name] for name in
           ("panel_starts_seconds", "vertices_m", "collocation_m", "normals")},
        **{name: np.stack((finite[name].real, finite[name].imag), axis=-1).tolist()
           for name in ("A", "gamma", "dual", "rhs")},
    }
    return {**{name: value for name, value in companion.items() if name != "records"},
            "core_radius_m": RADIUS, "records": [row]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--companion", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--case", type=int, nargs=3, metavar=("DIHEDRAL", "M", "N"))
    args = parser.parse_args()
    path = Path(__file__).with_name("reference.py")
    checked_bytes(path, ORACLE_SHA)
    oracle = runpy.run_path(str(path))
    companion = json.loads(checked_bytes(args.companion, COMPANION_SHA))
    cases = oracle["CASES"] if args.case is None else (tuple(args.case),)
    if not set(cases) <= set(oracle["CASES"]):
        parser.error("Case is outside the frozen six-condition roster")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    targets = [args.output_dir / f"d{d}-m{m}-n{n}-reference.json" for d, m, n in cases]
    for target in targets:
        if target.exists():
            raise FileExistsError(target)
    for key, target in zip(cases, targets, strict=True):
        body = build(oracle, companion, key)
        with target.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(body, stream, allow_nan=False, separators=(",", ":"))
            stream.write("\n")
        print(f"{hashlib.sha256(target.read_bytes()).hexdigest()}  {target.name}", flush=True)


if __name__ == "__main__":
    main()
