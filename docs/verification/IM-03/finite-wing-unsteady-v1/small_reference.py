# SPDX-License-Identifier: Apache-2.0
"""Independent M2/N4 harmonic circulation and load observations."""

import argparse
import hashlib
import json
import math
from itertools import pairwise
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
SOURCE = "tests/aero/fixtures/finite-wing-drag-refinement-v1.json"
SOURCE_SHA = "29b7cf5fd5301a0b4ebfc684e81b9cf70a90fa7305f34451a88e9da5cbd15878"
LENGTH = 128.0
SPEED = 20.0
DENSITY = 1.2
CHORD = 1.114285714285714
ORIGIN = np.array([0.375, 0.0, 0.0])


def segment(points, a, b, radius):
    """Velocity along z and its x derivative, from the line antiderivative."""
    edge = b - a
    length = np.linalg.norm(edge)
    unit = edge / length
    r = points - a
    s = r @ unit
    t = s - length
    perpendicular = r - s[:, None] * unit
    h2 = np.sum(perpendicular**2, axis=1) + radius**2
    d = np.sqrt(s**2 + h2)
    e = np.sqrt(t**2 + h2)
    dd = r[:, 0] / d
    de = (r[:, 0] - edge[0]) / e
    integral = (s / d - t / e) / h2
    derivative = (
        unit[0] * (1 / d - 1 / e)
        - s * dd / d**2
        + t * de / e**2
        - 2 * integral * perpendicular[:, 0]
    ) / h2
    # Rationalize the antiderivative difference for two same-sign projections.
    mask = s * t > 0
    u, v, p, q = s[mask], t[mask], d[mask], e[mask]
    dp, dq = dd[mask], de[mask]
    denominator = u * q + v * p
    stable = length * (u + v) / (p * q * denominator)
    integral[mask] = stable
    derivative[mask] = stable * (
        2 * unit[0] / (u + v)
        - dp / p
        - dq / q
        - (unit[0] * (p + q) + u * dq + v * dp) / denominator
    )
    cross_z = np.cross(unit, r)[:, 2]
    return cross_z * integral / (4 * math.pi), (cross_z * derivative - unit[1] * integral) / (
        4 * math.pi
    )


def starts(points, trailing, omega):
    first, axes = trailing[:-1], np.diff(trailing, axis=0)
    r = points[:, None, :] - first
    t = np.sum(r * axes, axis=2) / np.sum(axes**2, axis=1)
    minimum = np.min(np.linalg.norm(r - np.clip(t, 0, 1)[:, :, None] * axes, axis=2))
    interval = min(0.0005, minimum / (2 * SPEED))
    boundaries = [0.0]
    while boundaries[-1] < LENGTH / SPEED:
        boundaries.append(min(LENGTH / SPEED, boundaries[-1] + min(interval, math.pi / omega)))
        interval *= 1.5
    return boundaries


def wake_update(points, trailing, radius, omega, boundaries, *, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    halves = np.diff(boundaries) / 2
    ages = (np.array(boundaries[:-1])[:, None] + halves[:, None] * (nodes + 1)).ravel()
    phase = omega * ages
    weighted = (halves[:, None] * weights).ravel() * (
        -2 * np.sin(phase / 2) ** 2 - 1j * np.sin(phase)
    )
    shifted = points[:, None, :] - ages[None, :, None] * np.array([SPEED, 0.0, 0.0])
    shifted = shifted.reshape(-1, 3)
    columns = []
    for first, last in pairwise(trailing):
        _, derivative = segment(shifted, first, last, radius)
        r_last, r_first = shifted - last, shifted - first
        last_leg = r_last[:, 1] / (4 * math.pi * (np.sum(r_last**2, axis=1) + radius**2) ** 1.5)
        first_leg = r_first[:, 1] / (4 * math.pi * (np.sum(r_first**2, axis=1) + radius**2) ** 1.5)
        age_derivative = SPEED * (derivative + last_leg - first_leg)
        columns.append(age_derivative.reshape(len(points), -1) @ weighted)
    return np.array(columns).T


def quads(vertices):
    return np.stack(
        [vertices[:-1, :-1], vertices[:-1, 1:], vertices[1:, 1:], vertices[1:, :-1]], axis=2
    ).reshape(-1, 4, 3)


def observations(vertices, gamma, omega):
    surfaces = quads(vertices)
    gamma_rows = gamma.reshape(vertices.shape[0] - 1, vertices.shape[1] - 1)
    differences = np.diff(gamma_rows, axis=0, prepend=np.zeros_like(gamma_rows[:1])).ravel()
    left, right = surfaces[:, 0], surfaces[:, 1]
    edge_forces = DENSITY * differences[:, None] * np.cross([SPEED, 0, 0], right - left)
    force = edge_forces.sum(axis=0)
    moment = np.cross((left + right) / 2 - ORIGIN, edge_forces).sum(axis=0)
    for indices in ([0, 2, 1], [0, 3, 2]):
        triangle = surfaces[:, indices]
        area = np.cross(triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0]) / 2
        facet_force = DENSITY * 1j * omega * gamma[:, None] * area
        force += facet_force.sum(axis=0)
        moment += np.cross(triangle.mean(axis=1) - ORIGIN, facet_force).sum(axis=0)
    return np.array([force[2], moment[1]])


def generate():
    contents = (ROOT / SOURCE).read_bytes()
    if hashlib.sha256(contents).hexdigest() != SOURCE_SHA:
        raise ValueError("Unexpected horizontal-wing source fixture")
    source = json.loads(contents)
    case = next(
        row
        for row in source["cases"]
        if (row["dihedral_degrees"], row["chord_panels"], row["span_panels"]) == (0, 2, 4)
    )
    radius = source["core_m"]
    vertices = np.array(case["vertices_m"])[:15].reshape(3, 5, 3)
    points = np.array(case["collocations_m"])
    closed = vertices.copy()
    closed[-1, :, 0] += LENGTH
    base = np.array(
        [
            sum(segment(points, panel[i], panel[(i + 1) % 4], radius)[0] for i in range(4))
            for panel in quads(closed)
        ]
    ).T
    result = []
    for frequency in (0.05, 0.2, 0.8):
        omega = 2 * frequency * SPEED / CHORD
        panels = starts(points, vertices[-1], omega)
        matrices = []
        for quadrature in (16, 24):
            matrix = base.astype(complex)
            matrix[:, -4:] += wake_update(
                points, vertices[-1], radius, omega, panels, order=quadrature
            )
            matrices.append(matrix)
        for motion, rhs in (
            ("heave", np.full(len(points), 1j * omega)),
            ("pitch", -(SPEED + 1j * omega * (points[:, 0] - ORIGIN[0]))),
        ):
            gamma = [np.linalg.solve(matrix, rhs) for matrix in matrices]
            loads = [observations(vertices, amplitudes, omega) for amplitudes in gamma]
            difference = np.max(abs(loads[0] / loads[1] - 1))
            if difference > 1e-11:
                raise ValueError(
                    f"Age quadrature did not converge for {frequency} {motion}: {difference}"
                )
            result.append(
                {
                    "k": frequency,
                    "motion": motion,
                    "omega": omega,
                    "panel_starts": panels[:-1],
                    "rhs": [[float(x.real), float(x.imag)] for x in rhs],
                    "gamma": [[float(x.real), float(x.imag)] for x in gamma[0]],
                    "loads": [[float(x.real), float(x.imag)] for x in loads[0]],
                    "quadrature_relative_difference": float(difference),
                }
            )
    return {
        "schema": "finite-wing-harmonic-m2n4/1",
        "source_sha256": SOURCE_SHA,
        "wake_length_m": LENGTH,
        "core_m": radius,
        "inflow_m_s": SPEED,
        "density_kg_m3": DENSITY,
        "cases": result,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = generate()
    with args.output.open("x", encoding="utf-8", newline="\n") as target:
        json.dump(result, target, allow_nan=False, separators=(",", ":"))
        target.write("\n")


if __name__ == "__main__":
    main()
