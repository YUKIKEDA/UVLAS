# SPDX-License-Identifier: Apache-2.0
"""Finite closure and core sensitivity on the byte-identified coarse wings."""

import argparse
import hashlib
import importlib.util
import json
from decimal import Decimal as D
from decimal import localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
STEADY = ROOT / "docs/verification/IM-03/finite-wing-steady-v1/reference.py"
SOURCES = {
    STEADY: "7cd47141282ee114900b939405185feacb61fc782a269ca7c8b7d172aa7ead43",
    ROOT / "tests/aero/fixtures/finite-wing-steady-v1.json":
        "e4c0246d5ef9ddd978294d753e03ddc8831f44d0aac62824f04ca3707da9d667",
    ROOT / "rd/aero/finite-wing/records/finite-wing-span-extension-protocol.md":
        "0cf5e7dc4ff0483af6fe16dec289cbf5f50f0bf708eb1a7711b1126b23d40886",
    ROOT / "rd/aero/finite-wing/records/finite-wing-extended-campaign-result.json":
        "046bff5c6d808dc130c5bc993b3a0e9f56ae26d48a3e4380edd90c01a11e74bc",
}


def segment(oracle, query, start, end, core):
    if not core:
        start_to_query = oracle.subtract(query, start)
        edge = oracle.subtract(end, start)
        if not any(oracle.cross(start_to_query, edge)):
            along = oracle.dot(start_to_query, edge)
            if 0 <= along <= oracle.dot(edge, edge):
                raise ValueError("Ideal kernel queried on a non-excluded segment")
            return (D(0),) * 3
    return oracle.segment(query, start, end, core)


def observe(steady, oracle, wing, length, core):
    def vector(values):
        return tuple(D(str(x)) for x in values)

    physical = list(map(vector, wing["vertices_m"][:15]))
    closure = [(x + length, y, z) for x, y, z in physical[10:15]]
    if length == 32:
        for calculated, frozen in zip(closure, wing["vertices_m"][15:], strict=True):
            assert max(map(abs, oracle.subtract(calculated, vector(frozen)))) < D("1e-14")
    points = physical[:10] + closure
    collocations = list(map(vector, wing["collocations_m"]))
    flow = vector(wing["flow_m_s"])
    normals = []
    for t in map(vector, wing["tangents_y"]):
        normal = oracle.cross((D(1), D(0), D(0)), t)
        norm = oracle.dot(normal, normal).sqrt()
        normals.append(tuple(x / norm for x in normal))
    panels = [(i, i + 1, i + 6, i + 5) for i in (0, 1, 2, 3, 5, 6, 7, 8)]
    rings = [list(zip(panel, panel[1:] + panel[:1], strict=True)) for panel in panels]
    matrix = [
        [
            sum(oracle.dot(normal, segment(oracle, p, points[a], points[b], core)) for a, b in ring)
            for ring in rings
        ]
        for p, normal in zip(collocations, normals, strict=True)
    ]
    rhs = [-oracle.dot(normal, flow) for normal in normals]
    gamma = oracle.solve(matrix, rhs)
    assert max(abs(oracle.dot(row, gamma) - value) for row, value in zip(matrix, rhs, strict=True)) < D(
        "1e-65"
    )

    # Keep the physical trailing vertices. Combining a bound edge and its
    # wake continuation is valid for circulation, but cannot move a load point.
    load_points = physical + closure
    trailing_panels = [(i, i + 1, i + 6, i + 5) for i in (10, 11, 12, 13)]
    strengths = {}
    for panel, strength in zip(panels + trailing_panels, gamma + gamma[4:], strict=True):
        for a, b in zip(panel, panel[1:] + panel[:1], strict=True):
            key = (min(a, b), max(a, b))
            strengths[key] = strengths.get(key, D(0)) + (strength if a < b else -strength)
    force, moment = [D(0)] * 3, [D(0)] * 3
    band_force = [[D(0)] * 3 for _ in range(20)]
    band_moment = [[D(0)] * 3 for _ in range(20)]
    for (a, b), strength in strengths.items():
        if a >= 10 or b >= 15:
            continue  # Forces are observed only on physical bound edges.
        left, right = load_points[a], load_points[b]
        midpoint = tuple((x + y) / 2 for x, y in zip(left, right, strict=True))
        velocity = list(flow)
        for (i, j), source in strengths.items():
            if not source or (i, j) == (a, b):
                continue  # Physical straight self-edge, including every occurrence of it.
            induced = segment(oracle, midpoint, load_points[i], load_points[j], core)
            for axis, component in enumerate(induced):
                velocity[axis] += source * component
        f = tuple(steady.DENSITY * strength * x for x in oracle.cross(velocity, oracle.subtract(right, left)))
        m = oracle.cross(oracle.subtract(midpoint, steady.ORIGIN), f)
        band = min(int((midpoint[1] + 4) / D(".4")), 19)
        assert 0 <= band < 20
        for axis in range(3):
            force[axis] += f[axis]
            moment[axis] += m[axis]
            band_force[band][axis] += f[axis]
            band_moment[band][axis] += m[axis]
    sine, cosine = oracle.trig(1)
    lift = (-sine, D(0), cosine)
    return {
        "wake_length_m": length,
        "core_radius_m": core,
        "steady_columns_per_m": matrix,
        "rhs_m_s": rhs,
        "circulation_m2_s": gamma,
        "force_N": force,
        "moment_Nm": moment,
        "angles": [{
            "alpha_degrees": 1,
            "flow_m_s": flow,
            "circulation_m2_s": gamma,
            "force_N": force,
            "moment_Nm": moment,
            "drag_lift_N": [oracle.dot(force, (cosine, D(0), sine)), oracle.dot(force, lift)],
            "band_forces_N": band_force,
            "band_moments_Nm": band_moment,
            "band_lift_N": [oracle.dot(f, lift) for f in band_force],
        }],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    checked = {}
    for path, digest in SOURCES.items():
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != digest:
            raise ValueError(f"Frozen source changed: {path.name}")
        checked[path] = data
    spec = importlib.util.spec_from_file_location("wing_steady", STEADY)
    assert spec is not None and spec.loader is not None
    steady = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(steady)
    oracle = steady.original_oracle()
    wings = json.loads(steady.checked(steady.ORIGINAL, steady.ORIGINAL_HASH))
    prior = json.loads(checked[ROOT / "tests/aero/fixtures/finite-wing-steady-v1.json"])
    with localcontext() as context:
        context.prec = 75
        core = D(".001")
        cases = []
        for wing, previous in zip(wings["cases"], prior["cases"], strict=True):
            derived = observe(steady, oracle, wing, 32, core)
            for row, frozen in zip(derived["steady_columns_per_m"], previous["steady_columns_per_m"], strict=True):
                assert max(abs(a - D(b)) for a, b in zip(row, frozen, strict=True)) < D("1e-13")
            for key in ("circulation_m2_s", "force_N", "moment_Nm"):
                assert max(abs(a - D(b)) for a, b in zip(derived[key], previous[key], strict=True)) < D("1e-12")
            cases.append({
                "dihedral_degrees": wing["dihedral_degrees"],
                "baseline": derived,
                "closure": [observe(steady, oracle, wing, length, core) for length in (256, 512, 1024)],
                "core": [observe(steady, oracle, wing, 512, r) for r in (D(".0001"), D(".00001"))],
                "ideal_same_finite_closure": observe(steady, oracle, wing, 512, D(0)),
            })
    with args.output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump({"id": "finite-wing-wake-core-v1", "cases": cases}, handle, default=float, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
