# SPDX-License-Identifier: Apache-2.0
"""Frozen two-wing steady first-order independent comparison, decimal edition."""

import argparse
import hashlib
from decimal import Decimal as D, localcontext
import importlib.util
import json
from pathlib import Path

root = Path(__file__).resolve().parents[4]
FROZEN_REFERENCE = root / "rd/aero/multi-surface/records/multiple-surface-independent-result.json"
FROZEN_SHA256 = "9375e2825da70fd9da4d057e2c36c6129bd7a5a47ebfb0d2e01160d4fadb4048"
STEADY_SHA256 = "7cd47141282ee114900b939405185feacb61fc782a269ca7c8b7d172aa7ead43"
assert hashlib.sha256((root / "docs/verification/IM-03/finite-wing-steady-v1/reference.py").read_bytes()).hexdigest() == STEADY_SHA256
spec = importlib.util.spec_from_file_location("steady", root / "docs/verification/IM-03/finite-wing-steady-v1/reference.py")
assert spec and spec.loader
steady = importlib.util.module_from_spec(spec)
spec.loader.exec_module(steady)
oracle = steady.original_oracle()
shapes = (
    (D(4), D("1.5"), D(".6"), D(35), (D(0), D(0), D(0))),
    (D("2.5"), D(1), D(".4"), D(20), (D(5), D(0), D(".6"))),
)
normal = (D(0), D(0), D(1))

def geometry(m, n, shape):
    semispan, root_chord, tip_chord, sweep, origin = shape
    sin, cos = oracle.trig(sweep)
    tangent = sin / cos
    half = n // 2
    span_fraction = [(1 - oracle.trig(D(180) * i / half)[1]) / 2 for i in range(half + 1)]
    span_fraction[0], span_fraction[-1] = D(0), D(1)
    ys = [semispan * (value - 1) for value in span_fraction] + [semispan * value for value in span_fraction[1:]]
    xs = [(1 - oracle.trig(D(180) * i / m)[1]) / 2 for i in range(m + 1)]
    xs[0], xs[-1] = D(0), D(1)
    v = [(3 * a + b) / 4 for a, b in zip(xs[:-1], xs[1:], strict=True)] + [D(1)]
    c = [(a + 3 * b) / 4 for a, b in zip(xs[:-1], xs[1:], strict=True)]
    def position(x, y):
        return (origin[0] + abs(y) * tangent + x * (root_chord + (tip_chord - root_chord) * abs(y) / semispan),
                origin[1] + y, origin[2])
    vertices = [[position(x, y) for y in ys] for x in v]
    probes = [position(x, (a + b) / 2) for x in c for a, b in zip(ys[:-1], ys[1:], strict=True)]
    closed = vertices[:-1] + [[(p[0] + 1024, p[1], p[2]) for p in vertices[-1]]]
    rings = [(closed[i][j], closed[i][j+1], closed[i+1][j+1], closed[i+1][j])
             for i in range(m) for j in range(n)]
    return vertices, probes, rings

def wrench(gamma, vertices, m, n):
    force, moment = D(0), D(0)
    for i in range(m):
        for j in range(n):
            strength = gamma[i*n+j] - (gamma[(i-1)*n+j] if i else D(0))
            a, b = vertices[i][j:j+2]
            fz = D("1.2") * strength * 20 * (b[1] - a[1])
            force += fz
            moment -= (a[0] + b[0]) / 2 * fz
    return [force, moment]

def run(m, n, radius):
    left, right = (geometry(m, n, shape) for shape in shapes)
    probes = left[1] + right[1]
    rings = left[2] + right[2]
    matrix = [[oracle.ring(point, normal, ring, radius) for ring in rings] for point in probes]
    size = m*n
    gamma = oracle.solve(matrix, [D(-20)] * len(probes))
    returned = []
    for i, wing in enumerate((left, right)):
        begin, end = i * size, (i + 1) * size
        alone = oracle.solve([row[begin:end] for row in matrix[begin:end]], [D(-20)] * size)
        returned.append((wrench(gamma[begin:end], wing[0], m, n),
                         wrench(alone, wing[0], m, n)))
    residual = max(abs(sum(a*b for a,b in zip(row,gamma,strict=True))+20) for row in matrix)
    return {"m":m, "n":n, "radius_m":str(radius),"vertices_m":[left[0],right[0]],"queries_m":probes,
            "matrix_per_m":matrix,"circulation_slope_m2_s":gamma,"wrenches":returned,
            "residual_m_s":str(residual)}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    reference_bytes = FROZEN_REFERENCE.read_bytes()
    assert hashlib.sha256(reference_bytes).hexdigest() == FROZEN_SHA256
    reference = json.loads(reference_bytes)
    with localcontext() as context:
        context.prec = 75
        result = run(2, 12, D(".00001"))
        for i, (coupled, isolated) in enumerate(result["wrenches"]):
            frozen = reference["records"][0]["independent_wrenches"][i]
            defects = reference["records"][0]["coupling_change"][i]
            for both, single, target, defect in zip(coupled, isolated, frozen, defects, strict=True):
                error = abs((both-D(str(target)))/D(str(target)))
                assert error < D(".01"), (i, both, target)
                assert abs(abs((both-single)/single)-D(str(defect))) < D("1e-5")
                assert abs(abs((both-single)/single)-D(str(defect))) <= abs(D(str(defect))) / 100
        assert D(result["residual_m_s"]) < D("1e-65")
    with args.output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(result, handle, default=str, indent=2)
        handle.write("\n")
