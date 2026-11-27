"""Extract physical material probes from the frozen V-MAT coefficient state."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
sys.dont_write_bytecode = True

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[4]
RECORD = "rd/structure/material/records/material-shell-integration-v2-result.json"
RECORD_SHA256 = "c3d928367419e9123cc72227d07f5b648964170751c82ab93b983109271f8914"
ARCHIVE = "build/analysis/material-shell-integration-v2-20261105-130237.npz"
ARCHIVE_SHA256 = "405a1b6d4116af6d9da00c60551544f24f38adc67a25dff9eb9d00e494a51001"
PARENT = "rd/structure/material/experiments/check_material_shell_integration_v2.py"
PARAMETERS = (0.4, 0.2, 0.07, 0.1, -0.06, 0.025)
THICKNESS = 0.05
TOLERANCE = 1e-8


def require_sha(path: Path, expected: str) -> None:
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError(f"{path.name}: expected {expected}, got {actual}")


def relative(actual: np.ndarray, expected: np.ndarray) -> float:
    if actual.shape != expected.shape or not np.isfinite(actual).all():
        raise ValueError(
            "Frozen state and reconstructed observations must have the same finite shape"
        )
    return float(np.linalg.norm(actual - expected) / np.linalg.norm(expected))


def tensor(components: np.ndarray) -> np.ndarray:
    xx, yy, zz, yz, xz, xy = components
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])


def probes(archive: Path) -> dict:
    require_sha(ROOT / RECORD, RECORD_SHA256)
    require_sha(archive, ARCHIVE_SHA256)
    record = json.loads((ROOT / RECORD).read_text(encoding="utf-8"))
    sources = record["sources"]
    for name, digest in sources.items():
        require_sha(ROOT / name, digest)
    sys.path.insert(0, str((ROOT / PARENT).parent))
    import check_material_shell_integration_v2 as baseline  # noqa: PLC0415

    shell = baseline.shell
    model = shell.prepare(1, PARAMETERS, THICKNESS, True, "extended")
    material_input = record["material_input"]
    layers = baseline.inp.prepare_stack(material_input)
    data = baseline.data_for(model, layers)
    with np.load(archive, allow_pickle=False) as frozen:
        displacement = frozen["u"]
        reconstructed_energy, force = baseline.independent_point(data, displacement, material_input)
        agreement = {
            "reference_coordinates": relative(model["reference_coordinates"], frozen["q0"]),
            "force": relative(force, frozen["force"]),
            "mass": relative(np.kron(data["M"], np.eye(3)), frozen["fullM"]),
        }
    if not all(np.isfinite(value) and value <= TOLERANCE for value in agreement.values()):
        raise ValueError(f"Frozen basis replay exceeds its tolerance: {agreement}")

    # The interface points and depths are the ones named in the frozen experiment.
    # Add one witness in each layer's interior to cover all three constitutive laws.
    selections = [
        (0.31, 0.47, -0.70, 0),
        (0.31, 0.47, -0.37, 0),
        (0.31, 0.47, -0.37, 1),
        (0.67, 0.73, -0.08, 1),
        (0.67, 0.73, 0.22, 1),
        (0.67, 0.73, 0.22, 2),
        (0.67, 0.73, 0.62, 2),
    ]
    coordinates = np.array([row[:3] for row in selections])
    shapes, gradients, _, _, position = shell.kinematics(model, coordinates, THICKNESS, True)
    _, jacobian, _, _, _ = shell.geometry.reference(
        coordinates[:, :2], coordinates[:, 2], THICKNESS, PARAMETERS, True
    )
    deformation = np.eye(3) + np.einsum("ia,pib->pab", displacement, gradients)
    position_now = position + shapes @ displacement

    samples = []
    for i, (s, t, zeta, layer) in enumerate(selections):
        material = material_input[layer]
        request = material["material"]
        a3 = np.cross(jacobian[i, :, 0], jacobian[i, :, 1])
        a3 /= np.linalg.norm(a3)
        fibre = np.array(material["fiber_global"])
        a1 = fibre - np.dot(fibre, a3) * a3
        a1 /= np.linalg.norm(a1)
        axes = np.stack((a1, np.cross(a3, a1), a3), axis=1)
        f = deformation[i]
        green_world = (f.T @ f - np.eye(3)) / 2
        initial = np.array(request["e0"])
        initial[3:] /= np.sqrt(2)
        green_elastic = axes.T @ green_world @ axes - tensor(initial)
        engineering = np.array(
            [
                green_elastic[0, 0],
                green_elastic[1, 1],
                green_elastic[2, 2],
                2 * green_elastic[1, 2],
                2 * green_elastic[0, 2],
                2 * green_elastic[0, 1],
            ]
        )
        stress = baseline.independent_Ceng(request) @ engineering
        second = axes @ tensor(stress) @ axes.T
        first = f @ second
        samples.append(
            {
                "coordinates": [s, t, zeta],
                "layer": layer,
                "reference_position_metres": position[i].tolist(),
                "current_position_metres": position_now[i].tolist(),
                "deformation_gradient": f.tolist(),
                "reference_axes": axes.tolist(),
                "energy_joules_per_cubic_metre": float(engineering @ stress / 2),
                "second_pk_world_pascal": second.tolist(),
                "first_pk_world_pascal": first.tolist(),
                "cauchy_world_pascal": (first @ f.T / np.linalg.det(f)).tolist(),
            }
        )

    return {
        "schema": "uvlas.pg-frozen-material-probes/1",
        "record": {"path": RECORD, "sha256": RECORD_SHA256},
        "archive": {"path": ARCHIVE, "sha256": ARCHIVE_SHA256},
        "sources": sources,
        "replay_relative_errors": agreement,
        "integrated_energy_joules": reconstructed_energy,
        "material_input": material_input,
        "height_coefficients": [
            [i, j, str(value.numerator), str(value.denominator)]
            for (i, j), value in sorted(shell.geometry.coefficients(PARAMETERS).items())
        ],
        "probes": samples,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--archive", type=Path, default=ROOT / "out" / "rd-evidence" / ARCHIVE)
    args = parser.parse_args()
    result = probes(args.archive)
    with args.output.open("x", encoding="utf-8", newline="\n") as destination:
        json.dump(result, destination, indent=2, allow_nan=False)
        destination.write("\n")


if __name__ == "__main__":
    main()
