# SPDX-License-Identifier: Apache-2.0
"""Extend the frozen finite-wing fixture with all closed steady columns and loads."""

import argparse
import hashlib
import importlib.util
import json
from decimal import Decimal as D
from decimal import localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
ORIGINAL = ROOT / "tests/aero/fixtures/finite-wing-case-v1.json"
ORIGINAL_HASH = "0fa3597ffecbddba915652af463b5475a98e32a625d0cba2c3ab5b56093c6e01"
GENERATOR = ROOT / "docs/verification/IM-03/finite-wing-case-v1/reference.py"
GENERATOR_HASH = "18d84c7def0ee7a40426a4dcab71e6537b686e9f2addec5c5e972d81a3835950"
DENSITY = D("1.2")
ORIGIN = (D(".375"), D(0), D(0))


def checked(path, expected):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f"Frozen reference changed: {path.name}")
    return data


def original_oracle():
    checked(GENERATOR, GENERATOR_HASH)
    spec = importlib.util.spec_from_file_location("frozen_finite_wing", GENERATOR)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    checked(ROOT / module.ORIGIN, module.SHA256)
    return module


def all_observations(oracle, data, frozen):
    # Use the incoming frozen geometry and ambient values, rather than
    # approximating the inputs with a different wing or flow convention.
    points = [tuple(map(D, map(str, point))) for point in frozen["vertices_m"]]
    collocation = [tuple(map(D, map(str, point))) for point in frozen["collocations_m"]]
    tangents = [tuple(map(D, map(str, point))) for point in frozen["tangents_y"]]
    flow = tuple(map(D, map(str, frozen["flow_m_s"])))
    radius = D(float(data["core_radius_m"]))

    def corners(chord, span):
        first = 5 * chord + span
        return (first, first + 1, first + 6, first + 5)

    bound = [corners(c, s) for c in range(2) for s in range(4)]
    closed = [*points[:10], *points[15:]]
    rings = [tuple(closed[index] for index in face) for face in bound]
    norms = []
    for tangent in tangents:
        normal = oracle.cross((D(1), D(0), D(0)), tangent)
        norm = oracle.dot(normal, normal).sqrt()
        norms.append(tuple(value / norm for value in normal))
    matrix = [
        [oracle.ring(query, normal, ring, radius) for ring in rings]
        for query, normal in zip(collocation, norms, strict=True)
    ]
    rhs = [-oracle.dot(normal, flow) for normal in norms]
    gamma = oracle.solve(matrix, rhs)
    residual = [oracle.dot(row, gamma) - b for row, b in zip(matrix, rhs, strict=True)]
    assert max(map(abs, residual)) < D("1e-65")
    for row, expected in zip(matrix, frozen["first_bound_column"], strict=True):
        assert abs(row[0] - D(expected)) < D("1e-13")
    for value, expected in zip(gamma, frozen["circulation_m2_s"], strict=True):
        assert abs(value - D(expected)) < D("1e-13")

    # Signed circulation on each physical bound edge. The actual attached
    # edge never carries force. A straight edge's midpoint self-action is zero.
    strengths = {}
    for panel, circulation in zip(bound, gamma, strict=True):
        for start, end in zip(panel, panel[1:] + panel[:1], strict=True):
            key = (min(start, end), max(start, end))
            strengths[key] = strengths.get(key, D(0)) + (circulation if start < end else -circulation)

    force, moment = [D(0)] * 3, [D(0)] * 3
    for (start, end), circulation in strengths.items():
        if start >= 10 and end >= 10:
            continue
        midpoint = tuple((a + b) / 2 for a, b in zip(points[start], points[end], strict=True))
        velocity = list(flow)
        for ring, strength in zip(rings, gamma, strict=True):
            for a, b in zip(ring, ring[1:] + ring[:1], strict=True):
                term = oracle.segment(midpoint, a, b, radius)
                for axis, value in enumerate(term):
                    velocity[axis] += strength * value
        edge = oracle.subtract(points[end], points[start])
        edge_force = tuple(DENSITY * circulation * x for x in oracle.cross(velocity, edge))
        edge_moment = oracle.cross(oracle.subtract(midpoint, ORIGIN), edge_force)
        for axis in range(3):
            force[axis] += edge_force[axis]
            moment[axis] += edge_moment[axis]
    return {
        "dihedral_degrees": frozen["dihedral_degrees"],
        "steady_columns_per_m": [list(map(float, row)) for row in matrix],
        "rhs_m_s": list(map(float, rhs)),
        "circulation_m2_s": list(map(float, gamma)),
        "force_N": list(map(float, force)),
        "moment_Nm": list(map(float, moment)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    frozen = json.loads(checked(ORIGINAL, ORIGINAL_HASH))
    oracle = original_oracle()
    with localcontext() as context:
        context.prec = 75
        result = {
            "id": "finite-wing-steady-v1",
            "original_fixture_sha256": ORIGINAL_HASH,
            "density_kg_m3": float(DENSITY),
            "moment_origin_m": list(map(float, ORIGIN)),
            "cases": [all_observations(oracle, frozen, case) for case in frozen["cases"]],
        }
    with args.output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
