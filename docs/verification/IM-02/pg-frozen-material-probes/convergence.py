"""Replay frozen V-MAT observables and export smaller positive volume rules."""

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
import oracle
import transfer
from numpy.polynomial.legendre import legder, leggauss, legvander


def reference_data(archive):
    baseline, record = oracle.verified_baseline(archive)
    target, _, degree, frozen = transfer.source_displacement(baseline.shell, archive)
    return baseline, record["material_input"], target, degree, frozen


def kinematics(baseline, target, degree, points):
    _, J, _, H, dH = baseline.shell.geometry.reference(
        points[:, :2], points[:, 2], oracle.THICKNESS, oracle.PARAMETERS, True
    )
    g = baseline.shell.geometry.graph(points[:, :2], oracle.PARAMETERS)[4]
    dg = baseline.shell.geometry.metric_derivatives(points[:, :2], oracle.PARAMETERS)
    disp = np.empty((len(points), 3))
    grad = np.empty((len(points), 3, 3))
    for cx in range(2):
        for cy in range(2):
            use = ((points[:, 0] * 2).astype(int) == cx) & ((points[:, 1] * 2).astype(int) == cy)
            ps = points[use]
            if len(ps) == 0:
                continue
            a = legvander(4 * ps[:, 0] - (2 * cx + 1), degree)
            b = legvander(4 * ps[:, 1] - (2 * cy + 1), degree)
            coefficients = np.zeros((3, degree + 1, degree + 1, 3))
            for (field, x, y, i, j), values in target.items():
                if (cx, cy) == (x, y):
                    coefficients[field, i, j] = list(map(transfer.to_float, values))
            vals = np.einsum("na,nb,fabk->nfk", a, b, coefficients, optimize=True)
            deriv = [
                4
                * np.einsum(
                    "na,nb,fabk->nfk",
                    a[:, : der.shape[1]],
                    b[:, : der.shape[2]],
                    der,
                    optimize=True,
                )
                for der in (legder(coefficients, axis=1), legder(coefficients, axis=2))
            ]
            vals[:, 1] /= g[use, None] ** 1.5
            for j in range(2):
                deriv[j][:, 1] /= g[use, None] ** 1.5
                deriv[j][:, 1] -= 1.5 * vals[:, 1] * (dg[j][use] / g[use])[:, None]
            z = (H[use] * ps[:, 2] / 2)[:, None]
            through = vals[:, 1] + z * vals[:, 2]
            disp[use] = vals[:, 0] + z * vals[:, 1] + z * z * vals[:, 2] / 2
            coefs = [
                d[:, 0]
                + z * d[:, 1]
                + z * z * d[:, 2] / 2
                + through * (ps[:, 2] * dH[use, j] / 2)[:, None]
                for j, d in enumerate(deriv)
            ]
            coefs.append(through * (H[use] / 2)[:, None])
            grad[use] = np.stack(coefs, axis=-1) @ np.linalg.inv(J[use])
    return J, disp, grad


def eng(T):
    return np.stack(
        [T[:, 0, 0], T[:, 1, 1], T[:, 2, 2], 2 * T[:, 1, 2], 2 * T[:, 0, 2], 2 * T[:, 0, 1]],
        axis=-1,
    )


def stress_tensor(v):
    a, b, c, d, e, f = v.T
    return np.stack([a, f, e, f, b, d, e, d, c], axis=-1).reshape(len(v), 3, 3)


def actions(reference_model, rule):
    baseline, mats, target, degree = reference_model
    points, chart_weights, indices = rule
    J, ds, D = kinematics(baseline, target, degree, points)
    out = np.zeros(5)
    for index, m in enumerate(mats):
        idx = indices == index
        if not np.any(idx):
            continue
        request = m["material"]
        measure = chart_weights[idx] * np.linalg.det(J[idx])
        norm = np.cross(J[idx, :, 0], J[idx, :, 1])
        norm /= np.linalg.norm(norm, axis=1)[:, None]
        fibre = np.array(m["fiber_global"])
        fibre = fibre - (norm @ fibre)[:, None] * norm
        fibre /= np.linalg.norm(fibre, axis=1)[:, None]
        A = np.stack([fibre, np.cross(norm, fibre), norm], axis=-1)
        At = A.transpose(0, 2, 1)
        V = D[idx]
        F = np.eye(3) + V
        E = (F.transpose(0, 2, 1) @ F - np.eye(3)) / 2
        FV = F.transpose(0, 2, 1) @ V
        ED = (FV + FV.transpose(0, 2, 1)) / 2
        initial = np.array(request["e0"])
        initial[3:] /= np.sqrt(2)
        e = eng(At @ E @ A - oracle.tensor(initial))
        C = baseline.independent_Ceng(request)
        stress = e @ C.T
        S = A @ stress_tensor(stress) @ At
        SD = A @ stress_tensor(eng(At @ ED @ A) @ C.T) @ At
        P = F @ S
        PD = V @ S + F @ SD
        out += np.array(
            [
                measure @ (np.sum(e * stress, axis=1) / 2),
                measure @ np.sum(P * V, axis=(1, 2)),
                measure @ np.sum(PD * V, axis=(1, 2)),
                measure.sum() * request["rho"],
                measure @ np.sum(ds[idx] ** 2, axis=1) * request["rho"],
            ]
        )
    return out


def points_for(baseline, mats, order, depth, *, transverse=None):
    plane, wp = baseline.shell.geometry.quadrature(1, order)
    if transverse is not None:

        def interval_factors(n):
            nodes, weights = leggauss(n)
            return (
                np.concatenate([(cell + (nodes + 1) / 2) / 2 for cell in range(2)]),
                np.tile(weights / 4, 2),
            )

        s, ws = interval_factors(order)
        t, wt = interval_factors(transverse)
        ss, tt = np.meshgrid(s, t, indexing="ij")
        plane = np.column_stack((ss.ravel(), tt.ravel()))
        wp = np.multiply.outer(ws, wt).ravel()
    nodes, wz = leggauss(depth)
    points = []
    weights = []
    layer = []
    for i, m in enumerate(mats):
        lo, hi = float(m["lower"]), float(m["upper"])
        z = (lo + hi) / 2 + (hi - lo) / 2 * nodes
        points.append(np.column_stack([np.repeat(plane, len(z), axis=0), np.tile(z, len(plane))]))
        weights.append((wp[:, None] * wz * (hi - lo) / 2).ravel())
        layer.append(np.full(len(plane) * len(z), i))
    return np.concatenate(points), np.concatenate(weights), np.concatenate(layer)


def dyadic_rule(rule):
    points, weights, layers = rule
    points = np.rint(points * (1 << 28)) / (1 << 28)
    weights = np.rint(weights * (1 << 44)) / (1 << 44)
    if np.any(weights <= 0):
        raise ValueError("The exported finite rule must have positive weights")
    return points, weights, layers


def serialize_rule(rule):
    points, weights, layers = rule
    return [
        {
            "coordinates": [str(Fraction(float(value))) for value in coord],
            "layer": int(layer),
            "weight": str(Fraction(float(weight))),
        }
        for coord, layer, weight in zip(points, layers, weights, strict=True)
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("product_coefficients", type=Path)
    parser.add_argument(
        "--archive", type=Path, default=oracle.ROOT / "out" / "rd-evidence" / oracle.ARCHIVE
    )
    parser.add_argument("--mass-rule-output", type=Path)
    parser.add_argument("--elastic-rule-output", type=Path)
    args = parser.parse_args()
    baseline, materials, target, degree, frozen = reference_data(args.archive)
    fields = (
        "energy_joules",
        "displacement_internal_work_joules",
        "displacement_tangent_quadratic_joules",
        "mass_kg",
        "displacement_mass_quadratic_kg_metres_squared",
    )
    inputs = (baseline, materials, target, degree)

    samples = transfer.finite_volume_rule(args.product_coefficients)
    midpoint = actions(
        inputs,
        (
            np.array([list(map(float, row[1:4])) for row in samples]),
            np.array([float(row[4]) for row in samples]),
            np.array([row[0] for row in samples]),
        ),
    )
    direct = transfer.finite_actions(*inputs, samples)
    midpoint_reference = np.array([direct[field] for field in fields])
    if np.max(abs((midpoint - midpoint_reference) / midpoint_reference)) >= 1e-13:
        raise ValueError("Vector and scalar midpoint observations differ")
    reference = actions(inputs, points_for(baseline, materials, 45, 7))
    for i, field in enumerate(fields):
        if field in frozen and abs((reference[i] - frozen[field]) / frozen[field]) > 1e-8:
            raise ValueError(f"Frozen observable differs for {field}")
    print("45 x 45 x 7 replay", reference.tolist())
    print("order depth points", *fields)
    for order, depth in ((3, 3), (6, 3), (9, 1), (9, 2), (9, 3), (16, 2), (17, 2), (24, 3)):
        samples = points_for(baseline, materials, order, depth)
        values = actions(inputs, samples)
        print(order, depth, len(samples[0]), *abs((values - reference) / reference).tolist())

    for order, transverse in ((6, 9), (9, 17)):
        samples = points_for(baseline, materials, order, 2, transverse=transverse)
        values = actions(inputs, samples)
        print(
            f"{order} x {transverse}",
            2,
            len(samples[0]),
            *abs((values - reference) / reference).tolist(),
        )

    if args.mass_rule_output is not None:
        rule = dyadic_rule(points_for(baseline, materials, 9, 2))
        values = actions(inputs, rule)[3:]
        expected = np.array([reference[3], frozen[fields[4]]])
        if np.max(abs((values - expected) / expected)) > 1e-8:
            raise ValueError("The rounded finite rule failed a frozen mass observable")
        print("rounded 9 x 9 x 2", values.tolist(), abs((values - expected) / expected).tolist())
        result = {
            "points": serialize_rule(rule),
            "frozen_mass": float(reference[3]),
            "frozen_mass_work": frozen[fields[4]],
        }
        with args.mass_rule_output.open("x", encoding="utf-8", newline="\n") as output:
            json.dump(result, output)

    if args.elastic_rule_output is not None:
        rule = dyadic_rule(points_for(baseline, materials, 9, 2, transverse=17))
        values = actions(inputs, rule)[:3]
        expected = reference[:3]
        relative = abs((values - expected) / expected)
        if np.max(relative) > 1e-8:
            raise ValueError("The rounded finite rule failed a frozen elastic observable")
        print("rounded 9 x 17 x 2", values.tolist(), relative.tolist())
        result = {
            "schema": "uvlas.pg-frozen-material-elastic-rule/1",
            "archive_sha256": oracle.ARCHIVE_SHA256,
            "points": serialize_rule(rule),
            "reference": dict(zip(fields[:3], expected.tolist(), strict=True)),
            "finite": dict(zip(fields[:3], values.tolist(), strict=True)),
        }
        with args.elastic_rule_output.open("x", encoding="utf-8", newline="\n") as output:
            json.dump(result, output)


if __name__ == "__main__":
    main()
