# SPDX-License-Identifier: Apache-2.0
"""Independent closed-ring solution and midpoint drag on the frozen FW-V01/V02 geometry."""

import argparse
import hashlib
import importlib.util
import json
from decimal import Decimal as D
from decimal import localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
CASE = "docs/verification/IM-03/finite-wing-case-v1/reference.py"
LINES = "docs/verification/IM-03/finite-wing-span-lines-v1/reference.py"
RESULT = "rd/aero/finite-wing/records/drag-common-limit-result.json"
POLICY = "rd/aero/finite-wing/records/drag-common-limit-assessment-policy.md"
SOURCE_SHA256 = {
    CASE: "18d84c7def0ee7a40426a4dcab71e6537b686e9f2addec5c5e972d81a3835950",
    LINES: "6cdc1161e8e6a011fdf3a642e5bdc5f019730c6c3c5083937c01dadf8869c3a8",
    RESULT: "c875d9a995b4129b61e50a7d0fb1184108a3387f17fda4fd28eeebaaddfc9aa5",
    POLICY: "f573ce4f677f175161aca26c1e26ad678b9e88e58f798cdff9c0dbb517a38ea0",
}
DENSITY = D("1.2")
SPEED = D(20)
RADIUS = D("0.00001114285714285714")
CLOSURE = D(1024)


def checked(name):
    path = ROOT / name
    if hashlib.sha256(path.read_bytes()).hexdigest() != SOURCE_SHA256[name]:
        raise ValueError(f"Frozen source has changed: {path}")
    return path


def import_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def cosine_breakpoints(oracle, count):
    result = [D(0)]
    for index in range(1, count):
        _, cosine = oracle.trig(D(180) * index / count)
        result.append((1 - cosine) / 2)
    return [*result, D(1)]


def grids(oracle, chord, span, inclination, sweep):
    chord_breaks = cosine_breakpoints(oracle, chord)
    half_span = cosine_breakpoints(oracle, span // 2)
    y_breaks = [4 * (x - 1) for x in half_span] + [4 * x for x in half_span[1:]]

    def point(x, y):
        return abs(y) * sweep + x * (D("1.5") - D("0.225") * abs(y)), y, abs(y) * inclination

    bound = [(3 * x + y) / 4 for x, y in zip(chord_breaks[:-1], chord_breaks[1:], strict=True)]
    xyz = [point(x, y) for x in [*bound, D(1)] for y in y_breaks]
    xyz += [(x + CLOSURE, y, z) for x, y, z in xyz[-span-1:]]
    collocation = [
        point((x0 + 3 * x1) / 4, (y0 + y1) / 2)
        for x0, x1 in zip(chord_breaks[:-1], chord_breaks[1:], strict=True)
        for y0, y1 in zip(y_breaks[:-1], y_breaks[1:], strict=True)
    ]
    return chord_breaks, half_span, xyz, collocation


def panels(chord, span):
    pitch = span + 1
    for c in range(chord + 1):
        for s in range(span):
            first = c * pitch + s
            yield first, first + 1, first + pitch + 1, first + pitch


def ring_values(oracle, positions, vertices, normals, chord, span):
    rows = list(panels(chord, span))
    bound = chord * span
    values = []
    for point, normal in zip(positions, normals, strict=True):
        influences = [oracle.ring(point, normal, [vertices[k] for k in row], RADIUS) for row in rows]
        for s in range(span):
            influences[bound - span + s] += influences[bound + s]
        values.append(influences[:bound])
    return values


def induced_edge_drag(oracle, vertices, gamma, chord, span):
    signed_edges = {}
    for row, strength in zip(panels(chord, span), gamma + gamma[-span:], strict=True):
        for a, b in zip(row, row[1:] + row[:1], strict=True):
            pair = min(a, b), max(a, b)
            signed_edges[pair] = signed_edges.get(pair, D(0)) + (strength if a < b else -strength)
    te_start = chord * (span + 1)
    result = D(0)
    for (a, b), strength in sorted(signed_edges.items()):
        if a >= te_start or not strength:
            continue
        tangent = oracle.subtract(vertices[b], vertices[a])
        midpoint = tuple((u + v) / 2 for u, v in zip(vertices[a], vertices[b], strict=True))
        velocity = [D(0)] * 3
        for (i, j), source in signed_edges.items():
            if not source or (a, b) == (i, j):
                continue
            induced = oracle.segment(midpoint, vertices[i], vertices[j], RADIUS)
            for axis, value in enumerate(induced):
                velocity[axis] += source * value
        result += DENSITY * strength * oracle.cross(velocity, tangent)[0]
    return result


def observe(oracle, line, chord, span, angle, limit):
    sin, cos = oracle.trig(35)
    sweep = sin / cos
    sin, cos = oracle.trig(angle)
    inclination = sin / cos
    chord_breaks, half_span, xyz, collocation = grids(oracle, chord, span, inclination, sweep)
    normal_z = 1 / (1 + inclination**2).sqrt()
    normals = [(D(0), normal_z * inclination * (1 if p[1] < 0 else -1), normal_z)
               for p in collocation]
    matrix = ring_values(oracle, collocation, xyz, normals, chord, span)
    rhs = [-SPEED * n[2] for n in normals]
    gamma = oracle.solve(matrix, rhs)
    residual = max(abs(oracle.dot(row, gamma) - b) for row, b in zip(matrix, rhs, strict=True)) / SPEED
    assert residual < D("1e-64")
    near = induced_edge_drag(oracle, xyz, gamma, chord, span)
    te = chord * (span + 1)
    far = line.trefftz(oracle, xyz[te:te+span+1], gamma[-span:], DENSITY)
    limit_value = D(str(limit["reference_limit"]))
    uncertainty = D(str(limit["reference_relative_uncertainty"]))
    result = {
        "dihedral_degrees": angle,
        "chord_panels": chord,
        "span_panels": span,
        "chord_fractions": chord_breaks,
        "half_span_fractions": half_span,
        "vertices_m": xyz,
        "collocations_m": collocation,
        "gamma_derivative_m2_s_per_rad": gamma,
        "near_drag_N_per_rad2": near,
        "trefftz_drag_N_per_rad2": far,
        "relative_near_far_gap": abs(near / far - 1),
        "relative_near_to_limit_plus_uncertainty": abs(near / limit_value - 1) + uncertainty,
        "relative_original_residual": residual,
    }
    print(angle, chord, span, f"{near:.10f}", f"{far:.10f}", f"{abs(near / far - 1):.6%}", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--grid", action="append", choices=["2x4", "4x8", "4x12", "8x24"])
    parser.add_argument("--summary", action="store_true")
    args = parser.parse_args()
    oracle = import_path("frozen_case", checked(CASE))
    line = import_path("frozen_lines", checked(LINES))
    checked(POLICY)
    original = json.loads(checked(RESULT).read_bytes())["records"]
    levels = [tuple(map(int, x.split("x"))) for x in args.grid or ["2x4", "4x8"]]
    with localcontext() as decimal_context:
        decimal_context.prec = 75
        cases = [observe(oracle, line, m, n, record["dihedral"], record)
                 for m, n in levels for record in original]
    if args.summary:
        geometry = ("chord_fractions", "half_span_fractions", "vertices_m",
                    "collocations_m", "gamma_derivative_m2_s_per_rad")
        cases = [{key: value for key, value in case.items() if key not in geometry}
                 for case in cases]
    with args.output.open("x", encoding="utf-8", newline="\n") as output:
        json.dump({
            "id": "finite-wing-drag-refinement-v1", "sources": SOURCE_SHA256,
            "closure_m": float(CLOSURE), "core_m": float(RADIUS), "cases": cases,
        }, output, default=float, indent=2, allow_nan=False)
        output.write("\n")


if __name__ == "__main__":
    main()
