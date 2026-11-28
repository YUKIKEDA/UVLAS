# SPDX-License-Identifier: Apache-2.0
"""Exact geometry witnesses for the FD-06 wing-edge pairing design."""

import argparse
import hashlib
import json
import math
from bisect import bisect_right
from copy import deepcopy
from decimal import Decimal
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
RECORD = "rd/aero/multi-surface/records/multiple-surface-harmonic-refined-result.json"
SOURCE_HASHES = {
    RECORD: "72427831e87dd0a25f98ebccd2d0e3dac6f3bcbeb3ffa4c7181f6ba5a5b537bf",
    "rd/aero/multi-surface/experiments/run_multiple_surface_harmonic.py":
        "5c2815bc91a897a5ed4b82605aed29e7d06411a1b0d85d1e92f6de33827456a9",
    "rd/aero/finite-wing/experiments/uvlm_steady_reference.py":
        "36c6c6d03b6ef00bfff11697b9036c0130c364b8f9ac4a48534646d3d71e8ad1",
    "rd/aero/unsteady-reference/experiments/ptera_frequency_reference.py":
        "9664faf5e404fb84421e6728fe1ae6a4df8b389fc16b3e6086f3f9a4075f72d3",
}
CASE_HASHES = {
    (0, 2, 8): "3c1b172e594b47f781bdff7cb74b55270e3311bc76c35e8241f6f555e1c669ec",
    (0, 4, 16): "b7dabb8177107ee975d1c5cb5a1f1e9fabb02dbc7b0f2f7d9d8a8a34466df81d",
    (0, 8, 32): "f6e9000507c674ea88c3106cb80d3833b96816c0e63d6ccf01aaeee6f4456af5",
    (5, 2, 8): "a96dddf14e8e3bc54beffe9f05b26df570e773682837b3d4c9990b17186e3819",
    (5, 4, 16): "ca69ce26c62e848dfb0d80e5c23cc4d7cbda2a5b8e834adcafc2455ef1a9d96e",
    (5, 8, 32): "4cc3106d19eead36d6ac8080741f0223549aa91162595ea3f6c3c19fabaae926",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def checked_bytes(path, sha256):
    data = path.read_bytes()
    require(hashlib.sha256(data).hexdigest() == sha256, f"Changed reference: {path}")
    return data


def rational(value):
    if isinstance(value, dict):
        return Q(int(value["numerator"]), int(value["denominator"]))
    return Q(value)


def between(a, b, ratio):
    return tuple(Q(x) + ratio * (y - x) for x, y in zip(a, b, strict=True))


def series_at(knots, points, coordinate):
    i = min(bisect_right(knots, coordinate) - 1, len(knots) - 2)
    u = (coordinate - knots[i]) / (knots[i + 1] - knots[i])
    return between(points[i], points[i + 1], u)


def pair_by_chart(leading, trailing):
    if leading[0][1] > leading[-1][1]:
        leading = tuple(reversed(leading))
    if trailing[0][1] > trailing[-1][1]:
        trailing = tuple(reversed(trailing))
    lt, rt = (tuple(Q(t) for _, t in edge) for edge in (leading, trailing))
    require(lt[0] == rt[0] and lt[-1] == rt[-1] and all(
        a < b for knots in (lt, rt) for a, b in zip(knots, knots[1:])
    ), "Explicit pairing required")
    start, end = lt[0], lt[-1]
    return tuple(
        ((t - start) / (end - start), series_at(lt, leading, t), series_at(rt, trailing, t))
        for t in sorted(set(lt) | set(rt))
    )


def pair_explicitly(leading, trailing, records):
    eta, left, right = tuple(zip(*records, strict=True))
    require(eta[0] == 0 and eta[-1] == 1 and all(a < b for a, b in zip(eta, eta[1:])),
            "Invalid common knot sequence")
    refined = set(eta)
    for edge, parameters in ((leading, left), (trailing, right)):
        end = len(edge) - 1
        sign = 1 if parameters[0] < parameters[-1] else -1
        require({parameters[0], parameters[-1]} == {0, end} and all(
            sign * (b - a) > 0 for a, b in zip(parameters, parameters[1:])
        ), "Invalid edge parameter sequence")
        for i, (a, b) in enumerate(zip(parameters, parameters[1:])):
            for k in range(math.floor(min(a, b)) + 1, math.ceil(max(a, b))):
                refined.add(eta[i] + (eta[i + 1] - eta[i]) * (k - a) / (b - a))
    return tuple(
        (r, series_at(tuple(map(Q, range(len(leading)))), leading,
                      series_at(eta, tuple((p,) for p in left), r)[0]),
         series_at(tuple(map(Q, range(len(trailing)))), trailing,
                   series_at(eta, tuple((p,) for p in right), r)[0]))
        for r in sorted(refined)
    )


def edges(loop):
    return tuple(zip(loop, loop[1:] + loop[:1], strict=True))


def turn(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def on_segment(a, b, point):
    return turn(a, b, point) == 0 and all(
        min(x, y) <= z <= max(x, y) for x, y, z in zip(a, b, point, strict=True)
    )


def contacts(first, second):
    a, b = first
    c, d = second
    require(not (turn(a, b, c) * turn(a, b, d) < 0
                 and turn(c, d, a) * turn(c, d, b) < 0), "Proper edge crossing")
    return {point for point in (a, b, c, d)
            if on_segment(a, b, point) and on_segment(c, d, point)}


def boundary(loop):
    # Remove only collinear subdivision points that do not change direction.
    return tuple(b for a, b, c in zip(loop[-1:] + loop[:-1], loop,
                                     loop[1:] + loop[:1], strict=True)
                 if not on_segment(a, c, b))


def same_boundary(a, b):
    a, b = boundary(a), boundary(b)
    return len(a) == len(b) and any(
        a == cycle[i:] + cycle[:i] for cycle in (b, b[::-1]) for i in range(len(b))
    )


def verify_coverage(knots, domain):
    quads = tuple((left, right, next_right, next_left)
                  for (_, left, right), (_, next_left, next_right) in zip(knots, knots[1:]))
    turns = [turn(previous, current, following) for quad in quads
             for previous, current, following in zip(quad[-1:] + quad[:-1], quad,
                                                    quad[1:] + quad[:1], strict=True)]
    require(all(value > 0 for value in turns) or all(value < 0 for value in turns),
            "Fold, collapse, or changed orientation")
    sign = 1 if turns[0] > 0 else -1
    for i, j in combinations(range(len(quads)), 2):
        first, second = quads[i], quads[j]
        expected = {first[2], first[3]} if j == i + 1 else set()
        actual = set().union(*(contacts(a, b) for a in edges(first) for b in edges(second)))
        require(actual == expected, "Extra or missing strip contact")
        require(not any(all(sign * turn(a, b, p) > 0 for a, b in edges(target))
                        for source, target in ((first, second), (second, first)) for p in source),
                "Overlapping strips")
    left, right = (tuple(row[column] for row in knots) for column in (1, 2))
    require(same_boundary(left + right[::-1], domain), "Different outer boundary")


def point_at(knots, xi, eta):
    locations, left, right = tuple(zip(*knots, strict=True))
    return between(series_at(locations, left, eta), series_at(locations, right, eta), xi)


def rejects(operation, phrase):
    try:
        operation()
    except ValueError as error:
        require(phrase in str(error), f"Wrong rejection: {error}")
    else:
        raise ValueError(f"Negative control accepted: {phrase}")


def counterexamples():
    # Use lengths in metres for this witness. Dividing both coordinates by 5
    # yields the same witness inside the unit Case chart.
    leading = ((0, 0), (0, 5))
    trailing = ((1, 0), (1, 1), (4, 5))
    perimeter = leading + trailing[::-1]
    chart = pair_by_chart(leading[::-1], trailing)
    equal_arc = pair_explicitly(leading, trailing, (
        (Q(0), Q(0), Q(0)), (Q(1, 6), Q(1, 6), Q(1)), (Q(1), Q(1), Q(2)),
    ))
    for knots in (chart, equal_arc):
        verify_coverage(knots, perimeter)
    new, old = (point_at(knots, Q(3, 4), Q(1, 2)) for knots in (chart, equal_arc))
    require(new == (Q(51, 32), Q(5, 2)) and old == (Q(33, 20), Q(103, 40)),
            "Changed folded-edge witness")
    require(((old[0] - new[0]), (old[1] - new[1])) == (Q(9, 160), Q(3, 40)),
            "Changed position difference")
    require((old[1] - new[1]) / 10 == Q(3, 400), "Changed normal-flow difference")
    print(f"Controls: common span {new}, separate arclength {old}, normal difference 3/400 m/s")

    shifted = ((1, 1), (1, 6))
    rejects(lambda: pair_by_chart(leading, shifted), "Explicit pairing required")
    linear = ((Q(0), Q(0), Q(0)), (Q(1), Q(1), Q(1)))
    verify_coverage(pair_explicitly(leading, shifted, linear), leading + shifted[::-1])
    swapped = ((Q(0), Q(0), Q(1)), (Q(1), Q(1), Q(0)))
    rejects(lambda: verify_coverage(pair_explicitly(leading, shifted, swapped),
                                    leading + shifted[::-1]), "Fold")
    skips_vertex = ((Q(0), Q(0), Q(0)), (Q(1), Q(1), Q(2)))
    require(pair_explicitly(leading, trailing, skips_vertex)[1][0] == Q(1, 2),
            "Missing exact preimage of a boundary vertex")
    duplicate = (linear[0], (Q(1, 2), Q(1, 2), Q(0)), linear[1])
    rejects(lambda: pair_explicitly(leading, shifted, duplicate), "Invalid edge parameter")
    notched = ((0, 0), (0, 5), (4, 5), (1, 1), (1, 0), (Q(1, 2), Q(1, 4)))
    rejects(lambda: verify_coverage(chart, notched), "Different outer boundary")
    # Strictly convex strips turn once around the origin, then overlap.
    arc = ((Q(-1, 2), 1), (1, 1), (1, -1), (-1, -1), (-1, 1), (Q(1, 2), 1))
    winding = tuple((Q(i, len(arc) - 1), (x, y), (2 * x, 2 * y))
                    for i, (x, y) in enumerate(arc))
    rejects(lambda: verify_coverage(winding, arc), "Extra or missing strip contact")


def frozen_constants(side, dihedral):
    main = side < 2
    length, root = (Q(4), Q(3, 2)) if main else (Q(5, 2), Q(1))
    taper = Q(0.9 if main else 0.6)
    sweep = Q(float.fromhex("0x1.66819a3a0bf7ap-1" if main else "0x1.74b49cf3902d4p-2"))
    bank = Q(float.fromhex("0x1.665a8349d55e1p-4")) if main and dihedral else Q(0)
    base = (0, Q(0)) if main else (5, Q(0.6))
    return length, root, taper, sweep, bank, base


def frozen_position(side, dihedral, xi, eta):
    length, root, taper, sweep, bank, (x0, z0) = frozen_constants(side, dihedral)
    y = length * (eta - (side % 2 == 0))
    return x0 + sweep * abs(y) + xi * (root - taper * abs(y) / length), y, z0 + bank * abs(y)


def chart_declarations(dihedral, m, n):
    # These minimal declarations reproduce the graph fields used in the Case.
    def spacing(count):
        return [Q((1 - math.cos(math.pi * j / count)) / 2) for j in range(count + 1)]

    geometries, surfaces = [], []
    for side in range(4):
        length, root, taper, sweep, bank, (x0, z0) = frozen_constants(side, dihedral)
        negative = side % 2 == 0
        sign = -1 if negative else 1
        ends = (root, root - taper)[::sign]
        left, right = ((0, 0), (0, 1)), ((ends[0] / length, 0), (ends[1] / length, 1))
        identity = str(side + 1)
        sides = ((left[0], right[0]), (left[1], right[1]))
        geometries.append({
            "id": identity, "kind": "surface_graph",
            "frame": {"origin": (x0, -length if negative else 0, z0),
                      "rotation": ((1, 0, 0), (0, 1, 0), (0, 0, 1))},
            "reference_length": length, "width": 1,
            "x0": sweep if negative else 0, "skew": sign * sweep,
            "coefficients": [{"powers": (0, 0), "coefficient": bank if negative else 0},
                             {"powers": (0, 1), "coefficient": sign * bank}],
            "domain": left + right[::-1],
            "boundaries": [{"vertices": points} for points in (left, right, *sides)],
        })
        surfaces.append({"id": identity, "geometry": identity, "leading": 0, "trailing": 1,
                         "chord": spacing(m), "span": spacing(n // 2), "orientation": "positive",
                         "attachment": {"kind": "fixed_reference"},
                         "sharing": [{"other": str((side ^ 1) + 1),
                                      "own_boundary": 3 if negative else 2,
                                      "other_boundary": 2 if negative else 3,
                                      "orientation": "same"}]})
    return {"units": {"length": "m"}, "model": {"geometry": geometries, "aero_surfaces": surfaces}}


def chart_polygon(vertices):
    return tuple(tuple(map(rational, vertex)) for vertex in vertices)


def graph_position(geometry, chart):
    s, t = chart
    x0, width, skew = (rational(geometry[name]) for name in ("x0", "width", "skew"))
    z = sum(rational(term["coefficient"]) * s ** term["powers"][0] * t ** term["powers"][1]
            for term in geometry["coefficients"])
    local = (x0 + width * s + skew * t, t, z)
    length = rational(geometry["reference_length"])
    frame = geometry["frame"]
    return tuple(rational(translation) + length * sum(
        rational(coefficient) * value for coefficient, value in zip(row, local, strict=True)
    ) for translation, row in zip(frame["origin"], frame["rotation"], strict=True))


def verify_case(case, dihedral, m, n):
    require(case["units"]["length"] == "m", "Expected the frozen SI case")
    registry = {entry["id"]: entry for entry in case["model"]["geometry"]}
    surfaces = case["model"]["aero_surfaces"]
    require([entry["id"] for entry in surfaces] == ["1", "2", "3", "4"], "Wrong owner roster")
    seam = {}
    vertices = queries = bands = 0
    for side, surface in enumerate(surfaces):
        require(surface["orientation"] == "positive"
                and surface["attachment"]["kind"] == "fixed_reference"
                and "edge_pairing" not in surface, "Changed frozen variant")
        geometry = registry[surface["geometry"]]
        boundary_rows = [chart_polygon(row["vertices"]) for row in geometry["boundaries"]]
        knots = pair_by_chart(*(boundary_rows[surface[key]] for key in ("leading", "trailing")))
        verify_coverage(knots, chart_polygon(geometry["domain"]))
        chord, span = (tuple(map(rational, surface[key])) for key in ("chord", "span"))
        require(len(chord) == m + 1 and len(span) == n // 2 + 1
                and all(axis[0] == 0 and axis[-1] == 1
                        and all(a < b for a, b in zip(axis, axis[1:]))
                        for axis in (chord, span)), "Wrong lattice shape")
        quarter = tuple((3 * a + b) / 4 for a, b in zip(chord, chord[1:])) + (Q(1),)
        three_quarter = tuple((a + 3 * b) / 4 for a, b in zip(chord, chord[1:]))
        midpoint = tuple((a + b) / 2 for a, b in zip(span, span[1:]))

        def position(xi, eta):
            return graph_position(geometry, point_at(knots, xi, eta))

        for coordinates in ((quarter, span), (three_quarter, midpoint)):
            for xi in coordinates[0]:
                for eta in coordinates[1]:
                    require(position(xi, eta) == frozen_position(side, dihedral, xi, eta),
                            f"Incorrect source wing geometry for AeroSurface {surface['id']}")
        vertices += len(quarter) * len(span)
        queries += len(three_quarter) * len(midpoint)
        bands += len(span) - 1

        require(len(surface["sharing"]) == 1, "Wrong root-sharing declaration")
        registered, = surface["sharing"]
        root_eta = Q(1 if side % 2 == 0 else 0)
        expected = {"other": str((side ^ 1) + 1), "own_boundary": 2 + int(root_eta),
                    "other_boundary": 3 - int(root_eta), "orientation": "same"}
        require(expected == registered, "Wrong root-sharing provenance")
        shared_boundary = tuple(graph_position(geometry, point)
                                for point in boundary_rows[expected["own_boundary"]])
        require(shared_boundary == tuple(frozen_position(side, dihedral, xi, root_eta)
                                         for xi in (Q(0), Q(1))), "Wrong named root boundary")
        seam[surface["id"]] = tuple(position(xi, root_eta) for xi in quarter)
    require(seam["1"] == seam["2"] and seam["3"] == seam["4"] and seam["1"] != seam["3"],
            "Incorrect root owner or shared order")
    return vertices, queries, bands


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=Path, help="Also verify the six byte-pinned typed Case files")
    args = parser.parse_args()
    pinned = {name: checked_bytes(ROOT / name, digest) for name, digest in SOURCE_HASHES.items()}
    frozen = json.loads(pinned[RECORD])
    keys = tuple((row["dihedral"], row["m"], row["n"]) for row in frozen["records"])
    require(keys == tuple(CASE_HASHES), "Frozen series changed")
    counterexamples()
    for use_files in (False, True) if args.cases else (False,):
        totals = [0, 0, 0]
        for d, m, n in keys:
            label = f"d{d}-m{m}-n{n}"
            case = (json.loads(checked_bytes(args.cases / f"{label}.json", CASE_HASHES[d, m, n]),
                               parse_float=Decimal) if use_files else chart_declarations(d, m, n))
            result = verify_case(case, d, m, n)
            require(result == (4 * (m + 1) * (n // 2 + 1), 2 * m * n, 2 * n), "Lost grid rows")
            totals = [a + b for a, b in zip(totals, result, strict=True)]
            print(f"{'Case files' if use_files else 'Chart declarations'} {label}: "
                  f"{result[0]} vertices, {result[1]} queries, {result[2]} trailing bands")
            if (d, m, n) == keys[0]:
                swapped = deepcopy(case)
                for surface in swapped["model"]["aero_surfaces"]:
                    surface["geometry"] = str(((int(surface["geometry"]) - 1) ^ 2) + 1)
                rejects(lambda: verify_case(swapped, d, m, n), "Incorrect source wing geometry")
        require(totals == [1704, 1344, 224], "Frozen totals changed")
        print(f"Totals: {totals[0]} vertices, {totals[1]} queries, {totals[2]} trailing bands")


if __name__ == "__main__":
    main()
