# SPDX-License-Identifier: Apache-2.0
"""Independent edge loads, span-coordinate integrals and Trefftz observations."""

import argparse
import hashlib
import importlib.util
import json
from decimal import Decimal as D
from decimal import localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
POLARS = ROOT / "docs/verification/IM-03/finite-wing-polars-v1/reference.py"
PRIOR = ROOT / "tests/aero/fixtures/finite-wing-polars-v1.json"
SOURCES = {
    POLARS: "89554ba4e854ef78838fa36ef1dcca337ed65e1ffcbd065acd922a365ad32835",
    PRIOR: "50bf3906dfc16c81eb1eb8bbf53b0c2424cdb4b6d45b5f7db7c4566962146a13",
    ROOT / "rd/aero/avl-reference/spanwise-load-observation-contract.md":
        "3d1a590132c8f7a38f72a3b36491460e34d9685aab7b57c56ba65798a0e79595",
    ROOT / "rd/aero/finite-wing/records/induced-drag-common-limit-protocol.md":
        "662a9773f3378615a92d07f19fba1c1cf7bb918d0b97fd8650b434c936579c6c",
    ROOT / "rd/aero/finite-wing/records/trefftz-midpoint-operator-result.json":
        "331a99b7020ee04ab46bd4fdb20dff5c2cdf5122e46ba96f39358653fb034e35",
    ROOT / "rd/aero/finite-wing/experiments/trefftz_observations.py":
        "c69b5085a57476652bde1254683d0ccb89d52427275fb2be03cbe1bd4890776e",
}
BOUNDS = tuple(D(-4) + D(".4") * i for i in range(21))


def reference_version():
    for path, expected in SOURCES.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Frozen input has changed: {path}")
    spec = importlib.util.spec_from_file_location("frozen_polars", POLARS)
    assert spec is not None and spec.loader is not None
    polars = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(polars)
    return polars


def span_integrals(oracle, start, end, position, force, origin):
    """Integrate the moment in spatial y, independently of the test's t clips."""
    delta = oracle.subtract(end, start)
    midpoint = tuple((a + b) / 2 for a, b in zip(start, end, strict=True))
    couple = oracle.cross(oracle.subtract(position, midpoint), force)
    if not delta[1]:
        band = next(
            (i for i, (low, high) in enumerate(zip(BOUNDS[:-1], BOUNDS[1:], strict=True))
             if low <= position[1] < high),
            19,
        )
        assert BOUNDS[band] <= position[1] <= BOUNDS[band + 1]
        yield band, force, oracle.cross(oracle.subtract(position, origin), force)
        return
    lower, upper = sorted((start[1], end[1]))
    for i, (edge_low, edge_high) in enumerate(zip(BOUNDS[:-1], BOUNDS[1:], strict=True)):
        low, high = max(lower, edge_low), min(upper, edge_high)
        if high <= low:
            continue
        length = high - low
        # Integral of y - start_y, in the positive physical y direction.
        first_moment = ((high - start[1])**2 - (low - start[1])**2) / 2
        arm_integral = tuple(
            ((a - o) * length + distance * first_moment / delta[1]) / abs(delta[1])
            for a, o, distance in zip(start, origin, delta, strict=True)
        )
        weight = length / abs(delta[1])
        f = tuple(weight * v for v in force)
        moment = tuple(
            x + weight * c for x, c in zip(oracle.cross(arm_integral, force), couple, strict=True)
        )
        yield i, f, moment


def trefftz(oracle, te, traces, rho):
    jumps = [x - y for x, y in zip([D(0)] + traces, traces + [D(0)], strict=True)]
    drag = D(0)
    for trace, start, end in zip(traces, te[:-1], te[1:], strict=True):
        delta = oracle.subtract(end[1:], start[1:])
        center = tuple((a + b) / 2 for a, b in zip(start[1:], end[1:], strict=True))
        assert any(delta)
        projected = D(0)
        for jump, vertex in zip(jumps, te, strict=True):
            radius = oracle.subtract(center, vertex[1:])
            denominator = oracle.dot(radius, radius)
            if not denominator:
                raise ValueError("Ideal Trefftz query touches a streamwise vortex")
            projected += jump * oracle.dot(radius, delta) / denominator
        drag -= trace * projected
    return rho * drag / (4 * oracle.PI)


def edge_loads(steady, oracle, points, gamma, flow, radius):
    panels = [(i, i + 1, i + 6, i + 5) for i in (0, 1, 2, 3, 5, 6, 7, 8)]
    wake = [(i, i + 1, i + 6, i + 5) for i in (10, 11, 12, 13)]
    circulation = {}
    for panel, strength in zip(panels + wake, gamma + gamma[4:], strict=True):
        for first, second in zip(panel, panel[1:] + panel[:1], strict=True):
            edge = (min(first, second), max(first, second))
            circulation[edge] = circulation.get(edge, D(0)) + (strength if first < second else -strength)
    for (i, j), strength in sorted(circulation.items()):
        if i >= 10 or j >= 15:
            continue
        start, end = points[i], points[j]
        midpoint = tuple((a + b) / 2 for a, b in zip(start, end, strict=True))
        air = list(flow)
        for (a, b), other in circulation.items():
            if not other or (a, b) == (i, j):
                continue
            induced = oracle.segment(midpoint, points[a], points[b], radius)
            for axis, value in enumerate(induced):
                air[axis] += other * value
        force = tuple(
            steady.DENSITY * strength * v
            for v in oracle.cross(air, oracle.subtract(end, start))
        )
        yield i, j, midpoint, force


def observe(steady, oracle, wing, polar, radius):
    points = [tuple(D(str(x)) for x in vertex) for vertex in wing["vertices_m"]]
    samples = []
    for sample in polar["angles"]:
        gamma = sample["circulation_m2_s"]
        by_band = [[D(0)] * 3 for _ in range(20)]
        moment_by_band = [[D(0)] * 3 for _ in range(20)]
        edges = []
        for i, j, midpoint, force in edge_loads(
            steady, oracle, points, gamma, sample["flow_m_s"], radius
        ):
            edges.append({
                "endpoint_ids": [i + 1, j + 1], "point_m": midpoint, "force_N": force,
            })
            for band, f, m in span_integrals(
                oracle, points[i], points[j], midpoint, force, steady.ORIGIN
            ):
                for axis in range(3):
                    by_band[band][axis] += f[axis]
                    moment_by_band[band][axis] += m[axis]
        assert len(edges) == 18
        for quantities, expected in (
            (by_band, sample["force_N"]), (moment_by_band, sample["moment_Nm"])
        ):
            for component, target in zip(zip(*quantities, strict=True), expected, strict=True):
                assert abs(sum(component) - target) < D("1e-63")
        sine, cosine = oracle.trig(sample["alpha_degrees"])
        samples.append({
            "alpha_degrees": sample["alpha_degrees"],
            "edges": edges,
            "span_forces_N": by_band,
            "span_moments_Nm": moment_by_band,
            "span_lift_N": [oracle.dot((-sine, D(0), cosine), f) for f in by_band],
            "trefftz_drag_N": trefftz(oracle, points[10:15], gamma[4:], steady.DENSITY),
        })
    slopes = []
    for delta in (1, 2):
        radians = delta * oracle.PI / 180
        low, high, zero = samples[2 - delta], samples[2 + delta], samples[2]
        qs = steady.DENSITY * 20**2 * D("8.4") / 2
        near = [s["drag_lift_N"][0] for s in polar["angles"]]
        slopes.append({
            "delta_degrees": delta,
            "span_lift_per_rad": [
                (right - left) / (2 * radians * qs)
                for left, right in zip(low["span_lift_N"], high["span_lift_N"], strict=True)
            ],
            "near_drag_N_per_rad2": (near[2-delta] + near[2+delta] - 2*near[2]) / (2*radians**2),
            "trefftz_drag_N_per_rad2": (
                low["trefftz_drag_N"] + high["trefftz_drag_N"] - 2*zero["trefftz_drag_N"]
            ) / (2*radians**2),
        })
    plus_one = polar["angles"][3]
    gamma_prime = [20 * value / plus_one["flow_m_s"][2] for value in plus_one["circulation_m2_s"]]
    near_coefficient = sum(
        f[0] for _, _, _, f in edge_loads(steady, oracle, points, gamma_prime, (D(0),) * 3, radius)
    )
    return {
        "dihedral_degrees": wing["dihedral_degrees"], "angles": samples, "slopes": slopes,
        "quadratic_near_drag_N_per_rad2": near_coefficient,
        "quadratic_trefftz_drag_N_per_rad2": trefftz(oracle, points[10:15], gamma_prime[4:], steady.DENSITY),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    polars = reference_version()
    steady, shape, previous = polars.prerequisites()
    frozen = json.loads(PRIOR.read_bytes())
    oracle = steady.original_oracle()
    with localcontext() as context:
        context.prec = 75
        cases = []
        for wing, old, prior in zip(shape["cases"], previous["cases"], frozen["cases"], strict=True):
            detailed = polars.wing_data(steady, oracle, wing, old, D(shape["core_radius_m"]))
            for sample, published in zip(detailed["angles"], prior["angles"], strict=True):
                for key in ("force_N", "moment_Nm", "circulation_m2_s"):
                    assert all(
                        abs(value - D(expected)) <= D("1e-12")
                        for value, expected in zip(sample[key], published[key], strict=True)
                    )
            cases.append(observe(steady, oracle, wing, detailed, D(shape["core_radius_m"])))
    with args.output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump({"id": "finite-wing-span-lines-v1", "cases": cases}, handle, default=float, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
