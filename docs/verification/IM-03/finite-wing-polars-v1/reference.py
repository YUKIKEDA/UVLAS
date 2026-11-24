# SPDX-License-Identifier: Apache-2.0
"""Independent closed-ring reference for the five frozen finite-wing angles."""

import argparse
import importlib.util
import json
from decimal import Decimal as D
from decimal import localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
STEADY_SOURCE = ROOT / "docs/verification/IM-03/finite-wing-steady-v1/reference.py"
STEADY_HASH = "7cd47141282ee114900b939405185feacb61fc782a269ca7c8b7d172aa7ead43"
STEADY_FIXTURE = ROOT / "tests/aero/fixtures/finite-wing-steady-v1.json"
FIXTURE_HASH = "e4c0246d5ef9ddd978294d753e03ddc8831f44d0aac62824f04ca3707da9d667"
PROTOCOL = ROOT / "rd/aero/finite-wing/records/finite-wing-load-reference-specification.md"
PROTOCOL_HASH = "0b731ec1bcd8294c893665f41a00f109de14573cdf21dacd1ffda41148f59459"
ANGLES = (-2, -1, 0, 1, 2)
BOUNDS = tuple(D("-4") + D(".4") * i for i in range(21))


def prerequisites():
    import hashlib

    assert hashlib.sha256(STEADY_SOURCE.read_bytes()).hexdigest() == STEADY_HASH
    spec = importlib.util.spec_from_file_location("frozen_steady", STEADY_SOURCE)
    assert spec is not None and spec.loader is not None
    steady = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(steady)
    steady.checked(PROTOCOL, PROTOCOL_HASH)
    shape = json.loads(steady.checked(steady.ORIGINAL, steady.ORIGINAL_HASH))
    previous = json.loads(steady.checked(STEADY_FIXTURE, FIXTURE_HASH))
    return steady, shape, previous


def wing_data(steady, oracle, frozen, one_degree, radius):
    def vector(values):
        return tuple(D(str(x)) for x in values)

    points = list(map(vector, frozen["vertices_m"]))
    collocations = list(map(vector, frozen["collocations_m"]))
    normals = []
    for tangent in map(vector, frozen["tangents_y"]):
        normal = oracle.cross((D(1), D(0), D(0)), tangent)
        length = oracle.dot(normal, normal).sqrt()
        normals.append(tuple(component / length for component in normal))
    bound = [(i, i + 1, i + 6, i + 5) for i in (0, 1, 2, 3, 5, 6, 7, 8)]
    closed = points[:10] + points[15:]
    rings = [tuple(closed[i] for i in quad) for quad in bound]
    matrix = [
        [oracle.ring(point, normal, ring, radius) for ring in rings]
        for point, normal in zip(collocations, normals, strict=True)
    ]
    for row, expected_row in zip(matrix, one_degree["steady_columns_per_m"], strict=True):
        assert max(abs(got - D(expected)) for got, expected in zip(row, expected_row, strict=True)) < D(
            "1e-13"
        )

    def observation(degrees):
        sine, cosine = oracle.trig(degrees)
        flow = vector(frozen["flow_m_s"]) if degrees == 1 else vector(map(float, (20 * cosine, 0, 20 * sine)))
        gamma = oracle.solve(matrix, [-oracle.dot(normal, flow) for normal in normals])
        strengths = {}
        for panel, strength in zip(bound, gamma, strict=True):
            for start, end in zip(panel, panel[1:] + panel[:1], strict=True):
                key = (min(start, end), max(start, end))
                strengths[key] = strengths.get(key, D(0)) + (strength if start < end else -strength)

        force, moment = [D(0)] * 3, [D(0)] * 3
        band_force = [[D(0)] * 3 for _ in range(20)]
        band_moment = [[D(0)] * 3 for _ in range(20)]
        for (start, end), circulation in strengths.items():
            if start >= 10 and end >= 10:
                continue
            left, right = points[start], points[end]
            midpoint = tuple((x + y) / 2 for x, y in zip(left, right, strict=True))
            velocity = list(flow)
            for ring, strength in zip(rings, gamma, strict=True):
                for a, b in zip(ring, ring[1:] + ring[:1], strict=True):
                    for k, value in enumerate(oracle.segment(midpoint, a, b, radius)):
                        velocity[k] += strength * value
            f = tuple(steady.DENSITY * circulation * v for v in oracle.cross(velocity, oracle.subtract(right, left)))
            m = oracle.cross(oracle.subtract(midpoint, steady.ORIGIN), f)
            band = next(
                (i for i in range(20) if BOUNDS[i] <= midpoint[1] < BOUNDS[i + 1]),
                19,
            )
            assert BOUNDS[band] <= midpoint[1] <= BOUNDS[band + 1]
            for k in range(3):
                force[k] += f[k]
                moment[k] += m[k]
                band_force[band][k] += f[k]
                band_moment[band][k] += m[k]
        drag_axis = (cosine, D(0), sine)
        lift_axis = (-sine, D(0), cosine)
        return {
            "alpha_degrees": degrees,
            "flow_m_s": flow,
            "circulation_m2_s": gamma,
            "force_N": force,
            "moment_Nm": moment,
            "drag_lift_N": [oracle.dot(force, drag_axis), oracle.dot(force, lift_axis)],
            "band_forces_N": band_force,
            "band_moments_Nm": band_moment,
            "band_lift_N": [oracle.dot(f, lift_axis) for f in band_force],
        }

    samples = [observation(angle) for angle in ANGLES]
    plus_one = samples[3]
    for name in ("circulation_m2_s", "force_N", "moment_Nm"):
        assert max(abs(a - D(b)) for a, b in zip(plus_one[name], one_degree[name], strict=True)) < D("1e-12")
    qs = steady.DENSITY * 20**2 * D("8.4") / 2
    cref = D("1.114285714285714")
    slopes = []
    for delta in (1, 2):
        minus, plus = samples[2 - delta], samples[2 + delta]
        denominator = 2 * delta * oracle.PI / 180
        slopes.append(
            {
                "delta_degrees": delta,
                "cl_per_rad": (plus["drag_lift_N"][1] - minus["drag_lift_N"][1]) / (qs * denominator),
                "cmy_per_rad": (plus["moment_Nm"][1] - minus["moment_Nm"][1]) / (qs * cref * denominator),
                "band_cl_per_rad": [
                    (right - left) / (qs * denominator)
                    for left, right in zip(minus["band_lift_N"], plus["band_lift_N"], strict=True)
                ],
            }
        )
    return {"dihedral_degrees": frozen["dihedral_degrees"], "angles": samples, "slopes": slopes}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    steady, inputs, prior = prerequisites()
    oracle = steady.original_oracle()
    with localcontext() as ctx:
        ctx.prec = 75
        output = {
            "id": "finite-wing-polars-v1",
            "prior_fixture_sha256": FIXTURE_HASH,
            "protocol_sha256": PROTOCOL_HASH,
            "band_bounds_m": BOUNDS,
            "cases": [
                wing_data(steady, oracle, frozen, first, D(inputs["core_radius_m"]))
                for frozen, first in zip(inputs["cases"], prior["cases"], strict=True)
            ],
        }
    with args.output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(output, handle, default=float, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
