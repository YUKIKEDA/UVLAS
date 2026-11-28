# SPDX-License-Identifier: Apache-2.0
"""Independent positive-core companion for the six frozen V-MULTI harmonic cases.

The ideal computation is checked against the archived Ptera Gauss48 matrices,
circulation, and common loads. Finite-core values have their own oracle and
partition convergence checks. They never inherit the ideal reference identity.
"""

import argparse
import hashlib
import io
import json
import math
import tarfile
from itertools import product
from pathlib import Path

import numpy as np

RECORDS = "rd/aero/multi-surface/records/"
REFINED = RECORDS + "multiple-surface-harmonic-refined-result.json"
SOURCES = {
    REFINED: "72427831e87dd0a25f98ebccd2d0e3dac6f3bcbeb3ffa4c7181f6ba5a5b537bf",
    RECORDS + "multiple-surface-harmonic-result.json":
        "824c8c956d72abcb76710aa97e86fe7f1cf5bbc2f8c254bf7a3aabec0151da20",
    RECORDS + "multiple-surface-harmonic-protocol.md":
        "e50d31db6dca90d30b11df6014fd3e53674073557138f93b55187dd3723285df",
    RECORDS + "multiple-surface-harmonic-quadrature-extension.md":
        "91cf81eefd7e73e3d0b5863b46372992995bb970e81543ed210dc16ee949f44f",
    "rd/aero/multi-surface/experiments/run_multiple_surface_harmonic.py":
        "5c2815bc91a897a5ed4b82605aed29e7d06411a1b0d85d1e92f6de33827456a9",
}
CASES = tuple((d, m, n) for d, (m, n) in product((0, 5), ((2, 8), (4, 16), (8, 32))))
DENSITY = 1.2
CORE_RADIUS = 1e-5
FLOW = np.array([20.0, 0.0, 0.0])
LENGTH = 128.0
ORIGIN = np.array([0.375, 0.0, 0.0])
OMEGA = 8.0 / 1.114285714285714
TERMINAL = LENGTH / FLOW[0]
MOTIONS = ["main-heave", "main-pitch", "tail-heave", "tail-pitch"]


def archived_ideal(root, bundle, cases):
    for name, expected in SOURCES.items():
        raw = (root / name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != expected:
            raise ValueError(f"Frozen source differs: {name}")
    record = json.loads((root / REFINED).read_bytes())
    evidence = json.loads((root / "rd/evidence-index.json").read_bytes())
    bundle_id = evidence["bundles"][bundle.name]
    with bundle.open("rb") as stream:
        if (hashlib.file_digest(stream, "sha256").hexdigest() != bundle_id["sha256"]
                or stream.tell() != bundle_id["bytes"]):
            raise ValueError("Frozen evidence bundle differs")
    names = {}
    for case in cases:
        d, m, n = case
        name = f"build/analysis/multiple-surface-harmonic-ideal-q48/reference-d{d}-m{m}-n{n}.npz"
        names[name] = case
    outputs = {}
    with tarfile.open(bundle, "r:") as archive:
        for member in archive:
            if member.name not in names:
                continue
            case = names[member.name]
            identity = evidence["files"][member.name]
            if (not member.isfile() or case in outputs or member.size != identity["bytes"]
                    or bundle.name != identity["bundle"]):
                raise ValueError(f"Unexpected archive member: {member.name}")
            stream = archive.extractfile(member)
            if stream is None:
                raise ValueError(f"Missing archive member: {member.name}")
            with stream:
                raw = stream.read()
            digest = hashlib.sha256(raw).hexdigest()
            if digest != identity["sha256"] or digest != record["sources"][member.name]:
                raise ValueError(f"Frozen reference differs: {member.name}")
            with np.load(io.BytesIO(raw), allow_pickle=False) as archive_data:
                arrays = {key: archive_data[key] for key in ("A", "gamma", "rhs", "dual")}
            _, m, n = case
            panel_count = 2 * m * n
            expected = {"A": (panel_count, panel_count), "gamma": (panel_count, 4),
                        "rhs": (panel_count, 4), "dual": (2, 6, 4)}
            for key, shape in expected.items():
                values = arrays[key]
                if (values.shape != shape or not np.iscomplexobj(values)
                        or not np.isfinite(values).all()):
                    raise ValueError(f"Invalid {key} in {member.name}")
            outputs[case] = (member.name, digest, arrays)
    if outputs.keys() != set(cases):
        raise ValueError("Incomplete frozen V-MULTI archive")
    return outputs


def wing(m, n, dihedral, *, tail):
    chord_nodes = (1 - np.cos(np.linspace(0.0, math.pi, m + 1))) / 2
    semi_span_nodes = (1 - np.cos(np.linspace(0.0, math.pi, n // 2 + 1))) / 2
    span_unit = np.r_[semi_span_nodes - 1, semi_span_nodes[1:]]
    ring_fraction = np.r_[0.75 * chord_nodes[:-1] + 0.25 * chord_nodes[1:], 1.0]
    query_fraction = 0.25 * chord_nodes[:-1] + 0.75 * chord_nodes[1:]

    def coordinates(fractions, span):
        chord, span = np.broadcast_arrays(fractions[:, None], span[None, :])
        radius = np.abs(span)
        if tail:
            x = 5.0 + 2.5 * radius * math.tan(math.radians(20)) + chord * (1.0 - 0.6 * radius)
            return np.stack((x, 2.5 * span, np.full_like(x, 0.6)), axis=-1)
        x = 4.0 * radius * math.tan(math.radians(35)) + chord * (1.5 - 0.9 * radius)
        return np.stack((x, 4.0 * span, 4.0 * radius * math.tan(math.radians(dihedral))), axis=-1)

    vertices = coordinates(ring_fraction, span_unit)
    query = coordinates(query_fraction, (span_unit[:-1] + span_unit[1:]) / 2).reshape(-1, 3)
    normal = np.cross(FLOW, np.diff(vertices[:-1], axis=1)).reshape(-1, 3)
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    return vertices, query, normal


def quads(vertices):
    return np.stack((vertices[:-1, :-1], vertices[:-1, 1:],
                     vertices[1:, 1:], vertices[1:, :-1]), axis=2).reshape(-1, 4, 3)


def segment_projection(points, normals, start, end, radius):
    """Line-integrated Biot-Savart velocity, dotted with each receiver normal."""
    edge = end - start
    length = np.linalg.norm(edge, axis=1)
    unit = edge / length[:, None]
    offset = points[:, None, :] - start
    head = np.einsum("qec,ec->qe", offset, unit)
    tail = head - length
    transverse = offset - head[:, :, None] * unit
    distance_sq = np.einsum("qec,qec->qe", transverse, transverse) + radius**2
    first = np.sqrt(head**2 + distance_sq)
    last = np.sqrt(tail**2 + distance_sq)
    same_side = head * tail > 0
    integral = np.empty_like(head)
    opposite = ~same_side
    integral[opposite] = (
        (head[opposite] / first[opposite] - tail[opposite] / last[opposite])
        / distance_sq[opposite]
    )
    head_same, tail_same = head[same_side], tail[same_side]
    first_same, last_same = first[same_side], last[same_side]
    integral[same_side] = (
        np.broadcast_to(length, head.shape)[same_side] * (head_same + tail_same)
        / (first_same * last_same * (head_same * last_same + tail_same * first_same))
    )
    return np.einsum("qec,qc->qe", np.cross(unit, offset), normals) * integral / (4 * math.pi)


def ring_projection(points, normals, corners, radius):
    return sum(segment_projection(points, normals, corners[:, edge],
                                  corners[:, (edge + 1) % 4], radius) for edge in range(4))


def age_partition(points, rosters):
    endpoints = np.concatenate([vertices[-1, :-1] for vertices in rosters])
    edges = np.concatenate([np.diff(vertices[-1], axis=0) for vertices in rosters])
    offset = points[:, None, :] - endpoints
    fractions = np.einsum("qec,ec->qe", offset, edges) / np.sum(edges**2, axis=1)
    distances = np.linalg.norm(offset - np.clip(fractions, 0.0, 1.0)[:, :, None] * edges, axis=2)
    interval = min(0.0005, float(distances.min()) / (2 * FLOW[0]))
    ages = [0.0]
    while ages[-1] < TERMINAL:
        next_age = min(TERMINAL, ages[-1] + min(interval, math.pi / OMEGA))
        if next_age <= ages[-1]:
            raise ValueError("Age partition cannot advance")
        ages.append(next_age)
        interval *= 1.5
    return np.array(ages)


def subdivide(boundaries, count):
    widths = np.diff(boundaries)
    offsets = np.arange(count) / count
    starts = (boundaries[:-1, None] + widths[:, None] * offsets).ravel()
    if not np.all(np.diff(starts) > 0) or starts[-1] >= math.nextafter(TERMINAL, 0.0):
        raise ValueError("Age starts are not admissible")
    return np.r_[starts, boundaries[-1]]


def additive_wake(points, normals, rosters, radius, boundaries, order):
    """Direct integration of (-iω)e^(-iωa) applied to the age-a-to-end ring."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    halves = np.diff(boundaries) / 2
    ages = (boundaries[:-1, None] + halves[:, None] * (nodes + 1)).ravel()
    weight = (halves[:, None] * weights).ravel() * (-1j * OMEGA) * np.exp(-1j * OMEGA * ages)
    left = np.concatenate([vertices[-1, :-1] for vertices in rosters])
    right = np.concatenate([vertices[-1, 1:] for vertices in rosters])
    tail_left = left + np.array([LENGTH, 0.0, 0.0])
    tail_right = right + np.array([LENGTH, 0.0, 0.0])
    count = len(left)
    update = np.zeros((len(points), count), dtype=complex)
    for first in range(0, len(ages), 14):
        age = ages[first:first + 14]
        shift = age[:, None, None] * FLOW
        corners = np.stack((left + shift, right + shift,
                            np.broadcast_to(tail_right, (len(age), count, 3)),
                            np.broadcast_to(tail_left, (len(age), count, 3))), axis=2)
        rings = ring_projection(points, normals, corners.reshape(-1, 4, 3), radius)
        update += np.einsum("qae,a->qe", rings.reshape(len(points), len(age), count),
                            weight[first:first + len(age)])
    return update


def loads(vertices, gamma):
    m, n = np.subtract(vertices.shape[:2], 1)
    amplitudes = gamma.reshape(m, n, len(MOTIONS))
    horse = np.diff(amplitudes, axis=0, prepend=np.zeros_like(amplitudes[:1])).reshape(gamma.shape)
    panels = quads(vertices)
    left, right = panels[:, 0], panels[:, 1]
    convective = DENSITY * np.cross(FLOW, right - left)[:, None, :] * horse[:, :, None]
    force = convective.sum(axis=0)
    moment = np.cross(((left + right) / 2 - ORIGIN)[:, None, :], convective).sum(axis=0)
    for ids in ((0, 2, 1), (0, 3, 2)):
        triangle = panels[:, ids, :]
        area = np.cross(triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0]) / 2
        derivative = 1j * OMEGA * DENSITY * area[:, None, :] * gamma[:, :, None]
        force += derivative.sum(axis=0)
        moment += np.cross((triangle.mean(axis=1) - ORIGIN)[:, None, :], derivative).sum(axis=0)
    return np.concatenate((force.T, moment.T))


def received_context(d, m, n):
    wings = [wing(m, n, d, tail=False), wing(m, n, 0, tail=True)]
    rosters, by_wing_query, by_wing_normal = zip(*wings, strict=True)
    points = np.concatenate(by_wing_query)
    normals = np.concatenate(by_wing_normal)
    rhs = np.zeros((len(points), len(MOTIONS)), dtype=complex)
    for index, (query, normal, pivot) in enumerate(zip(by_wing_query, by_wing_normal, (0.375, 5.25),
                                                       strict=True)):
        wing_slice = slice(index * m*n, (index + 1) * m*n)
        rhs[wing_slice, 2*index] = 1j * OMEGA * normal[:, 2]
        rhs[wing_slice, 2*index + 1] = -(FLOW[0] + 1j * OMEGA * (query[:, 0] - pivot)) * normal[:, 2]
    closed = [vertices.copy() for vertices in rosters]
    for vertices in closed:
        vertices[-1, :, 0] += LENGTH
    panels = np.concatenate([quads(vertices) for vertices in closed])
    return rosters, points, normals, rhs, panels


def evaluate(base, rosters, points, normals, rhs, radius, boundaries, order):
    a = base.astype(complex)
    update = additive_wake(points, normals, rosters, radius, boundaries, order)
    count = len(points) // 2
    span = rosters[0].shape[1] - 1
    a[:, count - span:count] += update[:, :span]
    a[:, 2*count - span:2*count] += update[:, span:]
    gamma = np.linalg.solve(a, rhs)
    dual = np.stack([loads(vertices, gamma[index*count:(index+1)*count])
                     for index, vertices in enumerate(rosters)])
    residual = np.max(np.linalg.norm(a @ gamma - rhs, axis=0) / np.linalg.norm(rhs, axis=0))
    if not residual < 1e-11:
        raise ValueError(f"Independent original residual exceeded budget: {residual}")
    return {"A": a, "gamma": gamma, "dual": dual, "rhs": rhs}, float(residual)


def compare(candidate, reference):
    """Include the frozen aggregate and isolate RHSs and both wing loads."""
    def worst(key, axis=None):
        actual, expected = candidate[key], reference[key]
        return float(np.max(np.linalg.norm(actual - expected, axis=axis)
                            / np.linalg.norm(expected, axis=axis)))

    return {"matrix": worst("A"), "circulation": worst("gamma"),
            "circulation_by_rhs": worst("gamma", 0), "dual": worst("dual"),
            "dual_by_wing_and_rhs": worst("dual", 1), "rhs_by_input": worst("rhs", 0)}


def verify_case(case, archive, include_fixture):
    rosters, points, normals, rhs, panels = received_context(*case)
    boundaries = age_partition(points, rosters)
    refined = subdivide(boundaries, 4)
    name, frozen_sha, frozen_values = archive
    ideal_base = ring_projection(points, normals, panels, 0)
    ideal, ideal_residual = evaluate(ideal_base, rosters, points, normals, rhs, 0, refined, 48)
    ideal_match = compare(ideal, frozen_values)
    if max(ideal_match.values()) >= 1e-9:
        raise ValueError(f"Independent ideal computation misses the frozen reference at {case}: {ideal_match}")
    finite_base = ring_projection(points, normals, panels, CORE_RADIUS)
    reference, finite_residual = evaluate(finite_base, rosters, points, normals, rhs,
                                          CORE_RADIUS, refined, 48)
    deltas = {}
    for count in (1, 2, 4):
        partition = subdivide(boundaries, count)
        rough, _ = evaluate(finite_base, rosters, points, normals, rhs, CORE_RADIUS, partition, 16)
        deltas[count] = compare(rough, reference)
    # Establish a reference at independently refined order and intervals.
    narrower, _ = evaluate(finite_base, rosters, points, normals, rhs, CORE_RADIUS,
                          subdivide(boundaries, 2), 48)
    reference_delta = compare(narrower, reference)
    if max(reference_delta.values()) >= 1e-9 or max(deltas[4].values()) >= 1e-9:
        raise ValueError(f"Finite-core quadrature fails to converge at {case}: {deltas}")
    result = {
        "case": case, "frozen_ideal_source": name, "frozen_ideal_sha256": frozen_sha,
        "base_age_panels": len(boundaries)-1, "ideal_against_archive": ideal_match,
        "finite_gauss16_against_48": deltas, "finite_gauss48_subdivision": reference_delta,
        "finite_against_ideal": compare(reference, ideal), "ideal_l2_residual": ideal_residual,
        "finite_l2_residual": finite_residual,
    }
    if include_fixture:
        result.update(
            vertices_m=[vertices.tolist() for vertices in rosters],
            collocation_m=points.tolist(), normals=normals.tolist(),
            panel_starts_seconds=refined[:-1].tolist(),
            **{key: np.stack((value.real, value.imag), axis=-1).tolist()
               for key, value in reference.items()},
        )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", type=int, nargs=3, choices=None, metavar=("DIHEDRAL", "M", "N"))
    parser.add_argument("--fixture", action="store_true", help="Include input and reference arrays")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    cases = CASES if args.case is None else (tuple(args.case),)
    if not set(cases) <= set(CASES):
        parser.error("Unsupported frozen case")
    archive = archived_ideal(args.root, args.bundle, cases)
    observations = []
    for case in cases:
        row = verify_case(case, archive[case], args.fixture)
        print(json.dumps({key: value for key, value in row.items() if key not in {
            "vertices_m", "collocation_m", "normals", "panel_starts_seconds", "A", "gamma", "dual", "rhs"}},
                         allow_nan=False), flush=True)
        observations.append(row)
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(
            {"schema": "uvlas.vmulti-harmonic-positive-core-companion/1", "sources": SOURCES,
             "core_radius_m": CORE_RADIUS, "wake_length_m": LENGTH, "omega_rad_s": OMEGA,
             "velocity_m_s": FLOW.tolist(), "density_kg_m3": DENSITY,
             "moment_origin_m": ORIGIN.tolist(), "motions": MOTIONS, "records": observations},
            stream, allow_nan=False, separators=(",", ":"),
        )
        stream.write("\n")


if __name__ == "__main__":
    main()
