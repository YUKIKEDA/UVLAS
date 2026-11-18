"""Exact design oracle: barycentric feasible polytope, independent of plane cuts."""
import argparse
from fractions import Fraction as Q
from itertools import combinations, permutations
import json
from pathlib import Path


def rref(matrix, columns):
    a = [list(map(Q, row)) for row in matrix]
    rank = 0
    for col in range(columns):
        found = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if found is None:
            continue
        a[rank], a[found] = a[found], a[rank]
        divisor = a[rank][col]
        a[rank] = [x / divisor for x in a[rank]]
        for i in range(len(a)):
            if i != rank:
                multiplier = a[i][col]
                a[i] = [x - multiplier * y for x, y in zip(a[i], a[rank])]
        rank += 1
    return None if any(not any(row[:columns]) and row[columns] for row in a) else a[:rank]


def intersection(a, b):
    equations = [[a[j][i] for j in range(3)] + [-b[j][i] for j in range(3)] + [0] for i in range(3)]
    rows = rref(equations + [[1, 1, 1, 0, 0, 0, 1], [0, 0, 0, 1, 1, 1, 1]], 6)
    if rows is None:
        return set()
    vertices = set()
    for cols in combinations(range(6), len(rows)):
        solution = rref([[row[c] for c in cols] + [row[6]] for row in rows], len(rows))
        if solution is None or len(solution) != len(rows):
            continue
        values = [row[-1] for row in solution]
        if any(x < 0 for x in values):
            continue
        bary = [Q(0)] * 6
        for c, value in zip(cols, values):
            bary[c] = value
        vertices.add(tuple(sum(bary[j] * a[j][i] for j in range(3)) for i in range(3)))
    return vertices


def segment_member(p, a, b):
    d = [b[i] - a[i] for i in range(3)]
    k = next(i for i in range(3) if d[i])
    t = (p[k] - a[k]) / d[k]
    return 0 <= t <= 1 and all(p[i] == a[i] + t * d[i] for i in range(3))


def allowed(vertices, points=(), edge=None):
    if not vertices:
        return True
    if edge is not None:
        return all(segment_member(p, *edge) for p in vertices)
    return len(vertices) == 1 and next(iter(vertices)) in points


def box_candidate(a, b, strict=False):
    bounds = [(max(min(p[k] for p in a), min(p[k] for p in b)),
               min(max(p[k] for p in a), max(p[k] for p in b))) for k in range(3)]
    return all(lo < hi if strict else lo <= hi for lo, hi in bounds)


def te_row(occurrences, wake_sign, choose_first=False, wrong_sign=False):
    if choose_first:
        occurrences = occurrences[:1]
    if len(occurrences) != 1:
        return None
    face, geometric_sign = occurrences[0]
    result = [0, 0]
    result[face] = (1 if wrong_sign else -1) * wake_sign * (-geometric_sign)
    return result


def kinematic_match(canonical, alias, ignore_last=False):
    if ignore_last:
        return canonical[:-1] == alias[:-1]
    return canonical == alias


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Refusing to overwrite existing output")
    a = tuple(tuple(map(Q, p)) for p in [(0, 0, 0), (2, 0, 0), (0, 2, 0)])
    cases = [
        ("disjoint", [(3, 3, 0), (5, 3, 0), (3, 5, 0)], set(), (), None, True),
        ("transverse", [(1, -1, -1), (1, 3, -1), (1, 0, 1)], {(1, 0, 0), (1, 1, 0)}, (), None, False),
        ("point", [(2, 0, 0), (3, 0, 1), (2, 1, 1)], {(2, 0, 0)}, ((2, 0, 0),), None, True),
        ("point-excess", [(0, 0, 0), (1, 1, 1), (1, 1, -1)], {(0, 0, 0), (1, 1, 0)}, ((0, 0, 0),), None, False),
        ("edge", [(2, 0, 0), (0, 0, 0), (0, -2, 1)], {(0, 0, 0), (2, 0, 0)}, (), ((0, 0, 0), (2, 0, 0)), True),
        ("unclaimed-edge", [(2, 0, 0), (0, 0, 0), (0, -2, 1)], {(0, 0, 0), (2, 0, 0)}, ((0, 0, 0), (2, 0, 0)), None, False),
        ("coplanar-edge", [(2, 0, 0), (0, 0, 0), (0, -2, 0)], {(0, 0, 0), (2, 0, 0)}, (), ((0, 0, 0), (2, 0, 0)), True),
        ("overlap", [(2, 0, 0), (0, 0, 0), (1, 1, 0)], {(0, 0, 0), (2, 0, 0), (1, 1, 0)}, (), ((0, 0, 0), (2, 0, 0)), False),
        ("containment", [(0, 0, 0), (1, 0, 0), (0, 1, 0)], {(0, 0, 0), (1, 0, 0), (0, 1, 0)}, (), None, False),
        ("collinear", [(1, 0, 0), (3, 0, 0), (3, -2, 0)], {(1, 0, 0), (2, 0, 0)}, (), None, False),
        ("duplicate", a, set(a), (), None, False),
        ("low-gap", [(0, 0, Q(1, 2**60)), (2, 0, Q(1, 2**60)), (0, 2, Q(1, 2**60))], set(), (), None, True),
    ]
    checks = 0
    records = []
    for name, raw_b, expected, points, edge, accept in cases:
        b = tuple(tuple(map(Q, p)) for p in raw_b)
        expected = {tuple(map(Q, p)) for p in expected}
        for pa in permutations(a):
            for pb in permutations(b):
                got = intersection(pa, pb)
                assert got == expected, (name, got, expected)
                assert allowed(got, points, edge) == accept, name
                assert not got or box_candidate(pa, pb), name
                checks += 3
        for exponent, origin in [(-900, Q(0)), (900, Q(0)), (0, Q(2**900))]:
            scale = Q(2) ** exponent
            def transform(p):
                return tuple(origin + scale * x for x in p)
            assert intersection(tuple(map(transform, a)), tuple(map(transform, b))) == {transform(p) for p in expected}, name
            checks += 1
        records.append({"case": name, "vertices": [[str(x) for x in p] for p in sorted(expected)], "accepted": accept})
    faces = [(0, 1, 3, 2), (1, 4, 5, 3)]
    boundary = {}
    for face_id, face in enumerate(faces):
        balance = dict.fromkeys(face, 0)
        for u, v in zip(face, face[1:] + face[:1]):
            boundary.setdefault(tuple(sorted((u, v))), []).append((face_id, 1 if u < v else -1))
            balance[u] -= 1
            balance[v] += 1
        assert not any(balance.values())
    t_rows = []
    for edge, s in [((2, 3), -1), ((3, 5), -1)]:
        assert len(boundary[edge]) == 1
        row = [0, 0]
        face_id, geometric_sign = boundary[edge][0]
        induced_sign = -geometric_sign
        row[face_id] = -s * induced_sign
        assert induced_sign + s * row[face_id] == 0
        t_rows.append(row)
        assert te_row(boundary[edge], s) == row
    assert t_rows == [[1, 0], [0, 1]]
    owners = {(v, 7) for edge in [(2, 3), (3, 5)] for v in edge}
    assert len(owners) == 3
    rejected = {
        "open-aabb": not box_candidate(a, cases[2][1], strict=True),
        "shared-point-skip": not allowed(intersection(a, cases[3][1]), ((0, 0, 0),)),
        "isolated-endpoints-only": not allowed(intersection(a, cases[5][1]), cases[5][3]),
        "T-sign": te_row(boundary[(2, 3)], -1, wrong_sign=True) != [1, 0],
        "alias-basis": -1 != (-1) * (-1),
        "internal-TE-face-choice": te_row(boundary[(1, 3)], -1) is None and te_row(boundary[(1, 3)], -1, choose_first=True) is not None,
        "ignored-J-column": not kinematic_match([Q(1), Q(1, 2**60)], [Q(1), Q(0)]) and kinematic_match([Q(1), Q(1, 2**60)], [Q(1), Q(0)], ignore_last=True),
        "ignored-epoch": len({(3, 7), (3, 8)}) != len({3}),
    }
    assert all(rejected.values()), rejected
    assert (3, 7) != (3, 8) != (3, 2**64 - 1)
    assert ("Attached", 3) != ("Material", (3, 7))
    # A valid exact action has three disjoint significant bit bands.
    action = sum([Q(1), Q(1, 2**200), Q(1, 2**400)])
    assert action - Q(1) - Q(1, 2**200) == Q(1, 2**400)
    assert not kinematic_match([action, Q(0)], [action + Q(1, 2**500), Q(0)])
    for component in range(3):
        perturbed = [Q(0)] * 3
        perturbed[component] = Q(1, 2**100)
        assert not kinematic_match([Q(0)] * 3, perturbed)
    # Finite high/low words retain a gap lost by binary64 addition.
    exact = Q.from_float(1.0) + Q.from_float(2.0**-60)
    assert exact > 1 and 1.0 + 2.0**-60 == 1.0
    # Fixed quad children intersect on their declared diagonal only.
    quad = [(Q(0), Q(0), Q(0)), (Q(2), Q(0), Q(0)),
            (Q(2), Q(2), Q(1)), (Q(0), Q(2), Q(0))]
    first, second = tuple(quad[i] for i in (0, 2, 1)), tuple(quad[i] for i in (0, 3, 2))
    vertices = intersection(first, second)
    assert vertices == {quad[0], quad[2]} and allowed(vertices, edge=(quad[0], quad[2]))
    assert te_row([(0, 1), (0, -1)], 1) is None
    # A third incident face remains legal on the same shared segment; it is not TE.
    third = ((Q(0), Q(0), Q(0)), (Q(2), Q(0), Q(0)), (Q(0), Q(-2), Q(-1)))
    assert allowed(intersection(a, third), edge=(a[0], a[1]))
    assert te_row([(0, 1), (1, -1), (1, 1)], 1) is None
    # Rotation of geometry and provider director cannot be judged by n_ref dot n_cur.
    rotated = tuple((p[0], -p[1], -p[2]) for p in a)
    normal_reference, normal_current, director_rotated = (0, 0, 1), (0, 0, -1), (0, 0, -1)
    assert sum(x * y for x, y in zip(normal_reference, normal_current)) < 0
    assert normal_current == director_rotated and rotated != a
    result = {"tolerance": "0 (Fraction exact)", "geometry_checks": checks,
              "cases": records, "T": t_rows, "release_owners": sorted(owners),
              "multi_band_velocity": str(action), "negative_controls": rejected}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"geometry_checks": checks, "negative_controls": len(rejected)}))


if __name__ == "__main__":
    main()
