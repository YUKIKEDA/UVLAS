# SPDX-License-Identifier: Apache-2.0
"""Independent Decimal closed-ring specialization of the frozen finite-wing proof."""

import argparse
import hashlib
import json
from decimal import Decimal as D
from decimal import localcontext
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
ORIGIN = "rd/aero/finite-wing/records/finite-wing-specialization-result.json"
SHA256 = "56da441111e16b149f6030916cb3decd3c091ba1744fa9073c22077188dca719"
PI = D("3.141592653589793238462643383279502884197169399375105820974944592307816406286")


def trig(degrees):
    x = PI * D(degrees) / 180
    s, c, si, co = x, D(1), x, D(1)
    for k in range(1, 42):
        si = -si * x * x / ((2 * k) * (2 * k + 1))
        co = -co * x * x / ((2 * k - 1) * (2 * k))
        s += si
        c += co
    return s, c


def dot(a, b):
    return sum(x * y for x, y in zip(a, b, strict=True))


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b, strict=True))


def cross(a, b):
    return tuple(
        a[(k + 1) % 3] * b[(k + 2) % 3] - a[(k + 2) % 3] * b[(k + 1) % 3] for k in range(3)
    )


def segment(p, a, b, radius):
    ra, rb, edge = subtract(p, a), subtract(p, b), subtract(b, a)
    n = cross(ra, rb)
    strength = (
        dot(edge, ra) / (dot(ra, ra) + radius * radius).sqrt()
        - dot(edge, rb) / (dot(rb, rb) + radius * radius).sqrt()
    ) / (4 * PI * (dot(n, n) + radius * radius * dot(edge, edge)))
    return tuple(strength * value for value in n)


def ring(p, normal, four, radius):
    return sum(
        dot(segment(p, a, b, radius), normal)
        for a, b in zip(four, four[1:] + four[:1], strict=True)
    )


def solve(matrix, values):
    ns = len(values)
    rows = [[*list(row), value] for row, value in zip(matrix, values, strict=True)]
    for i in range(ns):
        pivot = max(range(i, ns), key=lambda j: abs(rows[j][i]))
        rows[i], rows[pivot] = rows[pivot], rows[i]
        divisor = rows[i][i]
        rows[i] = [item / divisor for item in rows[i]]
        for j, row in enumerate(rows):
            if i != j:
                rows[j] = [value - row[i] * lead for value, lead in zip(row, rows[i], strict=True)]
    return [row[-1] for row in rows]


def rational(value):
    x = Fraction(value)
    return {"numerator": str(x.numerator), "denominator": str(x.denominator)}


def chart(side, sweep, inclination):
    corner = [D(".15"), D(".375")] if side == 0 else [D(".375"), D(".15")]
    shift = sweep if side else -sweep
    slope = inclination if side else -inclination
    lower, upper = [(rational(0), rational(0)), (rational(0), rational(1))]
    lower_te, upper_te = [(rational(chord), rational(t)) for t, chord in enumerate(corner)]
    identification = str(side + 1)
    height = []
    if inclination:
        if side == 0:
            height.append({"powers": [0, 0], "coefficient": rational(inclination)})
        height.append({"powers": [0, 1], "coefficient": rational(slope)})
    geometry = {
        "kind": "surface_graph",
        "id": identification,
        "frame": {
            "origin": [0, -4 if side == 0 else 0, 0],
            "rotation": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        },
        "reference_length": 4,
        "x0": rational(0 if side else sweep),
        "width": rational(1),
        "skew": rational(shift),
        "degree": [0, 1] if inclination else [0, 0],
        "coefficients": height,
        "domain": [lower, lower_te, upper_te, upper],
        "boundaries": [
            {"name": "leading", "vertices": [lower, upper]},
            {"name": "trailing", "vertices": [lower_te, upper_te]},
            {"name": "span-start", "vertices": [lower, lower_te]},
            {"name": "span-end", "vertices": [upper, upper_te]},
        ],
        "continuity": [],
    }
    surface = {
        "id": identification,
        "geometry": identification,
        "attachment": {"kind": "fixed_reference"},
        "chord": [0, 0.5, 1],
        "span": [0, 0.5, 1],
        "leading": 0,
        "trailing": 1,
        "orientation": "positive",
        "sharing": [
            {
                "other": str(2 - side),
                "own_boundary": 2 if side else 3,
                "other_boundary": 3 if side else 2,
                "orientation": "same",
            }
        ],
    }
    return geometry, surface


def check_chart(graph, xi, y, expected):
    def decode(value):
        return D(value["numerator"]) / D(value["denominator"])

    x0, width, skew = (decode(graph[name]) for name in ("x0", "width", "skew"))
    xref, y0, z0 = graph["frame"]["origin"]
    length = D(graph["reference_length"])
    t = (y - y0) / length
    bottom, top = (decode(graph["domain"][i][0]) for i in (1, 2))
    s = xi * (bottom + (top - bottom) * t)
    height = D(0)
    for term in graph["coefficients"]:
        assert term["powers"][0] == 0
        height += decode(term["coefficient"]) * (t if term["powers"][1] else D(1))
    actual = (xref + length * (x0 + width * s + skew * t), y, z0 + length * height)
    assert max(abs(a - b) for a, b in zip(actual, expected, strict=True)) < D("1e-28")


def build(dihedral):
    opening, closing = trig(35)
    sweep = opening / closing
    opening, closing = trig(dihedral)
    inclination = opening / closing

    def point(xi, y):
        y = D(y)
        return (abs(y) * sweep + xi * (D("1.5") - D(".225") * abs(y)), y, abs(y) * inclination)

    ys = (-4, -2, 0, 2, 4)
    grid = [[point(D(xi), y) for y in ys] for xi in (".125", ".625", "1")]
    open_end = [(x + 32, y, z) for x, y, z in grid[-1]]
    closed = [*grid[:2], open_end]
    panel_rings = [
        (closed[i][j], closed[i][j + 1], closed[i + 1][j + 1], closed[i + 1][j])
        for i in range(2)
        for j in range(4)
    ]
    collocation = [point(D(xi), y) for xi in (".375", ".875") for y in (-3, -1, 1, 3)]
    tangents_y = [subtract(grid[i][j + 1], grid[i][j]) for i in range(2) for j in range(4)]
    normals = []
    for tangent in tangents_y:
        n = cross((D(1), D(0), D(0)), tangent)
        length = dot(n, n).sqrt()
        normals.append(tuple(x / length for x in n))

    # Match the exact binary64 radius in the frozen specialization.
    radius = D("0.001")
    matrix = [
        [ring(q, normal, panel, radius) for panel in panel_rings]
        for q, normal in zip(collocation, normals, strict=True)
    ]
    sin, cos = trig(1)
    flow = tuple(D(str(float(x))) for x in (20 * cos, D(0), 20 * sin))
    gamma = solve(matrix, [-dot(n, flow) for n in normals])
    residual = [dot(row, gamma) + dot(n, flow) for row, n in zip(matrix, normals, strict=True)]
    assert max(map(abs, residual)) < D("1e-65")
    graphs, surfaces = zip(
        *(chart(i, round(sweep, 30), round(inclination, 30)) for i in range(2)), strict=True
    )

    for xi in (D(".125"), D(".375"), D(".625"), D(".875"), D(1)):
        for side, graph in enumerate(graphs):
            for y in range(-4 + 4 * side, 1 + 4 * side):
                check_chart(graph, xi, D(y), point(xi, y))

    def numbers(points):
        return [[float(x) for x in p] for p in points]

    return {
        "dihedral_degrees": dihedral,
        "case_geometry": graphs,
        "case_surfaces": surfaces,
        "flow_m_s": list(map(float, flow)),
        "vertices_m": numbers([p for row in grid for p in row] + open_end),
        "collocations_m": numbers(collocation),
        "tangents_y": numbers(tangents_y),
        "circulation_m2_s": list(map(float, gamma)),
        "first_bound_column": [float(row[0]) for row in matrix],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if hashlib.sha256((ROOT / ORIGIN).read_bytes()).hexdigest() != SHA256:
        raise ValueError("Frozen specialization record hash mismatch")
    with localcontext() as context:
        context.prec = 75
        data = {
            "id": "finite-wing-case-v1",
            "source": ORIGIN,
            "source_sha256": SHA256,
            "chord_panels": 2,
            "span_panels": 4,
            "wake_length_m": 32,
            "core_radius_m": 0.001,
            "angle_degrees": 1,
            "cases": [build(0), build(5)],
        }
    with args.output.open("x", encoding="utf-8", newline="\n") as target:
        json.dump(data, target, ensure_ascii=False, indent=2)
        target.write("\n")


if __name__ == "__main__":
    main()
