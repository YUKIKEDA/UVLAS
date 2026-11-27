"""Transfer the frozen V-MAT displacement to the selected rational PG basis."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from functools import cache
from pathlib import Path

import oracle

# The oracle sets the replay's BLAS environment before either library loads NumPy.
# isort: split
import numpy as np
from flint import arb_mat, ctx, fmpq, fmpq_mat
from scipy.linalg import qr


def binary64(value):
    return fmpq(*float(value).as_integer_ratio())


def to_float(value):
    return int(value.numerator) / int(value.denominator)


def source_displacement(shell, archive: Path):
    with np.load(archive, allow_pickle=False) as frozen:
        displacement = frozen["u"]
        vector = displacement.ravel()
        actions = {
            "displacement_internal_work_joules": float(frozen["force"].ravel() @ vector),
            "displacement_tangent_quadratic_joules": float(vector @ frozen["shell_K"] @ vector),
            "displacement_mass_quadratic_kg_metres_squared": float(
                vector @ frozen["fullM"] @ vector
            ),
        }
    indices, basis, _, record = shell.source_basis(1, oracle.PARAMETERS, "extended")
    rank = basis.ncols()
    thickness_indices = np.arange(4)
    scalar_count = rank + len(thickness_indices)
    if displacement.shape != (scalar_count, 3):
        raise ValueError("Frozen displacement size does not match the reconstructed basis")
    # Recreate the pre-interpolation normalized-to-source transform at the
    # frozen precision, keeping exact dyadic coordinates for the transfer.
    ctx.prec = 512
    ctx.threads = 4
    mass, _ = shell.geometric_mass(
        1,
        oracle.PARAMETERS,
        oracle.THICKNESS,
        True,
        "extended",
        indices,
        basis,
        record["d_degree"],
        thickness_indices,
        record["d_degree"] + 2,
    )
    factor = shell.cholesky(mass)
    exact_state = fmpq_mat([[binary64(value) for value in row] for row in displacement])
    transformed = factor.transpose().solve(arb_mat(exact_state))
    coordinates = fmpq_mat(
        [[transformed[row, i].mid().fmpq() for i in range(3)] for row in range(scalar_count)]
    )

    cell_basis = {}
    for field in range(3):
        for cx in range(2):
            for cy in range(2):
                if field < 2:
                    coefficients = shell.field_coefficients(
                        1, indices, basis, record["d_degree"], cx, cy, field, rank
                    )
                    rows = {
                        key: values + [fmpq(0)] * len(thickness_indices)
                        for key, values in coefficients.items()
                    }
                else:
                    coefficients = shell.c_polynomials(1, cx, cy, thickness_indices)
                    rows = {
                        key: [fmpq(0)] * rank + [values.get(j, fmpq(0)) for j in thickness_indices]
                        for key, values in coefficients.items()
                    }
                keys = list(rows)
                cell = fmpq_mat(list(rows.values())) * coordinates
                for row, key in enumerate(keys):
                    cell_basis[(field, cx, cy, *key)] = [cell[row, i] for i in range(3)]
    return cell_basis, scalar_count, record["d_degree"], actions


def kinematics(shell, target, degree, points):
    reference, jacobian, _, height, thickness_derivatives = shell.geometry.reference(
        points[:, :2], points[:, 2], oracle.THICKNESS, oracle.PARAMETERS, True
    )
    metric = shell.geometry.graph(points[:, :2], oracle.PARAMETERS)[4]
    dg = shell.geometry.metric_derivatives(points[:, :2], oracle.PARAMETERS)
    displacements = []
    gradients = []
    for site, (s, t, zeta) in enumerate(points):
        cx, cy = int(2 * s), int(2 * t)
        u, v = 4 * s - (2 * cx + 1), 4 * t - (2 * cy + 1)
        polys = np.zeros((3, degree + 1, degree + 1, 3))
        for (field, x, y, a, b), values in target.items():
            if (x, y) == (cx, cy):
                polys[field, a, b, :] = list(map(to_float, values))

        def evaluate(coefficients, u, v):
            return np.array(
                [
                    [np.polynomial.legendre.legval2d(u, v, field[:, :, i]) for i in range(3)]
                    for field in coefficients
                ]
            )

        fields = evaluate(polys, u, v)
        derivatives = [
            4 * evaluate(np.polynomial.legendre.legder(polys, axis=j), u, v) for j in (1, 2)
        ]
        g = metric[site]
        fields[1] /= g**1.5
        for a in range(2):
            derivatives[a][1] = derivatives[a][1] / g**1.5 - 1.5 * fields[1] * dg[a][site] / g
        z = height[site] * zeta / 2
        through = fields[1] + z * fields[2]
        displacements.append(fields[0] + z * fields[1] + z * z * fields[2] / 2)
        chart_gradient = [
            d[0] + z * d[1] + z * z * d[2] / 2 + through * zeta * thickness_derivatives[site, a] / 2
            for a, d in enumerate(derivatives)
        ]
        chart_gradient.append(through * height[site] / 2)
        gradients.append(np.array(chart_gradient).T @ np.linalg.inv(jacobian[site]))
    return reference, jacobian, np.array(displacements), np.array(gradients)


def check_archive_probes(shell, target, degree):
    fixture = oracle.ROOT / "tests/structure/fixtures/pg-frozen-material-probes-v1.json"
    oracle.require_sha(fixture, "83ea890fbc83d651d63c51477c197b5e2ed8c167b08e6714758fa3715d51df97")
    probes = json.loads(fixture.read_text(encoding="utf-8"))["probes"]
    points = np.array([probe["coordinates"] for probe in probes])
    reference, _, displacement, gradients = kinematics(shell, target, degree, points)
    position = reference + displacement
    gradient = np.eye(3) + gradients
    worst_position = np.max(abs(position - [probe["current_position_metres"] for probe in probes]))
    worst_gradient = np.max(abs(gradient - [probe["deformation_gradient"] for probe in probes]))
    if worst_position > 1e-12 or worst_gradient > 1e-10:
        raise ValueError("Recreated source coefficients do not match the frozen kinematics")
    return {"position_metres": float(worst_position), "deformation_gradient": float(worst_gradient)}


def finite_volume_rule(path: Path):
    samples = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("PG_VOLUME "):
            _, layer, s, t, zeta, weight = line.split()
            samples.append((int(layer), *(Fraction(x) for x in (s, t, zeta, weight))))
    if len(samples) != len(set(samples)) or not samples or sum(row[4] for row in samples) != 2:
        raise ValueError("Midpoint rule must cover the unit square and the whole thickness")
    for layer, s, t, zeta, weight in samples:
        if not (0 <= layer < 3 and 0 < s < 1 and 0 < t < 1 and -1 < zeta < 1 and weight > 0):
            raise ValueError("Midpoint rule contains an invalid sample")
    return samples


def engineering(tensor):
    return np.array(
        [
            tensor[0, 0],
            tensor[1, 1],
            tensor[2, 2],
            2 * tensor[1, 2],
            2 * tensor[0, 2],
            2 * tensor[0, 1],
        ]
    )


def finite_actions(baseline, materials, target, degree, samples):
    points = np.array([list(map(float, row[1:4])) for row in samples])
    _, jacobians, displacement, gradients = kinematics(baseline.shell, target, degree, points)
    energy = force_work = tangent_work = mass = mass_work = 0.0
    for i, (layer, _, _, _, chart_weight) in enumerate(samples):
        material = materials[layer]
        request = material["material"]
        jacobian = jacobians[i]
        measure = float(chart_weight) * np.linalg.det(jacobian)
        normal = np.cross(jacobian[:, 0], jacobian[:, 1])
        normal /= np.linalg.norm(normal)
        fibre = np.array(material["fiber_global"])
        fibre -= (fibre @ normal) * normal
        fibre /= np.linalg.norm(fibre)
        axes = np.stack([fibre, np.cross(normal, fibre), normal], axis=1)
        variation = gradients[i]
        deformation = np.eye(3) + variation
        green = (deformation.T @ deformation - np.eye(3)) / 2
        mixed = deformation.T @ variation
        green_direction = (mixed + mixed.T) / 2
        initial = np.array(request["e0"])
        initial[3:] /= np.sqrt(2)
        strain = engineering(axes.T @ green @ axes - oracle.tensor(initial))
        stiffness = baseline.independent_Ceng(request)
        stress = stiffness @ strain
        second = axes @ oracle.tensor(stress) @ axes.T
        second_direction = (
            axes @ oracle.tensor(stiffness @ engineering(axes.T @ green_direction @ axes)) @ axes.T
        )
        first = deformation @ second
        first_direction = variation @ second + deformation @ second_direction
        energy += measure * (strain @ stress) / 2
        force_work += measure * np.sum(first * variation)
        tangent_work += measure * np.sum(first_direction * variation)
        density_measure = measure * request["rho"]
        mass += density_measure
        mass_work += density_measure * (displacement[i] @ displacement[i])
    return {
        "points": len(samples),
        "energy_joules": float(energy),
        "displacement_internal_work_joules": float(force_work),
        "displacement_tangent_quadratic_joules": float(tangent_work),
        "mass_kg": float(mass),
        "displacement_mass_quadratic_kg_metres_squared": float(mass_work),
    }


@cache
def local_power(degree, convert):
    # Both conventions use 4*n*s-(2*cx+1) as the local argument.
    coefficients = {(degree, 0): Fraction(1)}
    polys = convert(coefficients, degree)
    return {i: fmpq(value.numerator, value.denominator) for (i, _), value in polys.items()}


def product_basis(path: Path, convert):
    result = defaultdict(dict)
    ordinals = set()
    seen = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("PG_COEFFICIENT "):
            continue
        _, column, field, cx, cy, a, b, coefficient = line.split()
        column, field, cx, cy, a, b = map(int, (column, field, cx, cy, a, b))
        entry = column, field, cx, cy, a, b
        if entry in seen or min(entry) < 0 or field >= 3 or max(cx, cy) >= 2:
            raise ValueError(f"Invalid or repeated product coefficient: {entry}")
        seen.add(entry)
        ordinals.add(column)
        scale = fmpq(coefficient)
        for i, x in local_power(a, convert).items():
            for j, y in local_power(b, convert).items():
                row = result[field, cx, cy, i, j]
                row[column] = row.get(column, fmpq(0)) + scale * x * y
    if not ordinals or ordinals != set(range(len(ordinals))):
        raise ValueError("Missing product columns")
    return result, len(ordinals)


def compare(target, archive_count, product, count):
    if count != archive_count:
        raise ValueError(f"Product rank {count} differs from frozen rank {archive_count}")
    keys = sorted(set(target) | set(product))
    rows = [[product.get(key, {}).get(j, fmpq(0)) for j in range(count)] for key in keys]
    rhs = [target.get(key, [fmpq(0)] * 3) for key in keys]
    double = np.array([[to_float(value) for value in row] for row in rows])
    row_scale = np.max(abs(double), axis=1)
    scaled = double / np.where(row_scale > 0, row_scale, 1)[:, None]
    # Floating QR orders candidate equations. Exact elimination and replay of
    # every coefficient, including absent rows, establish the actual equality.
    _, _, pivot = qr(scaled.T, pivoting=True, mode="economic")
    selected = pivot[:count]
    matrix = fmpq_mat([rows[i] for i in selected])
    coefficients = matrix.solve(fmpq_mat([rhs[i] for i in selected]))
    rebuilt = fmpq_mat(rows) * coefficients
    error = rebuilt - fmpq_mat(rhs)
    for i, key in enumerate(keys):
        for j in range(3):
            if error[i, j]:
                raise ValueError(f"Nonzero exact reconstruction residual for {key}, component {j}")
    return {
        "rank": count,
        "verified_rows": len(keys),
        "displacement_metres": [[str(coefficients[i, j]) for j in range(3)] for i in range(count)],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("product_coefficients", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--archive", type=Path, default=oracle.ROOT / "out" / "rd-evidence" / oracle.ARCHIVE
    )
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    baseline, record = oracle.verified_baseline(args.archive)
    shell = baseline.shell
    product, count = product_basis(args.product_coefficients, shell.exact.p.to_legendre)
    samples = finite_volume_rule(args.product_coefficients)
    target, archive_count, degree, actions = source_displacement(shell, args.archive)
    errors = check_archive_probes(shell, target, degree)
    transfer = compare(target, archive_count, product, count)
    result = {
        "schema": "uvlas.pg-frozen-material-transfer/1",
        "archive_sha256": oracle.ARCHIVE_SHA256,
        "transfer": transfer,
        "rounded_archive_replay_maximum_absolute_error": errors,
        "frozen_actions": actions,
        "midpoint_actions": finite_actions(
            baseline, record["material_input"], target, degree, samples
        ),
    }
    with args.output.open("x", encoding="utf-8", newline="\n") as destination:
        json.dump(result, destination, indent=2, allow_nan=False)
        destination.write("\n")


if __name__ == "__main__":
    main()
