"""Export one actual projection-failure witness from an isolated frozen experiment copy."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

PAYLOAD = r"""
import hashlib, json, sys
from pathlib import Path
import numpy as np
import scipy

sys.path.insert(0, str(Path.cwd() / "rd"))
import rd_paths
import check_uvlm_column_update as experiment

if "--probe" in sys.argv:
    print(json.dumps({"imported": True, "numpy": np.__version__, "scipy": scipy.__version__}))
    raise SystemExit(0)
m, k = 4096, 0.1
op = experiment.real_operator(m)
factor = experiment.lu_factor(np.array(op["A"], order="F"), overwrite_a=True, check_finite=False)
x = op["xc"]
v, integral_error = experiment.quad_vec(
    lambda t: experiment.O.D.density(1 + t, x, 5000.0) * np.expm1(-2j * k * t),
    0.0,
    160.0,
    epsrel=1e-11,
    epsabs=1e-12,
    limit=4096,
)
b = np.column_stack(
    [
        1j * 16 * k * mode["H"] - (8 + 1j * 16 * k * (x - mode["xp"])) * mode["alpha"]
        for mode in experiment.O.MODES
    ]
)
solve, denominator, ratio = experiment.update_solver(factor, v)
hi, lo, report = experiment.correction(op, v, b, solve)
if lo is None:
    raise RuntimeError("Frozen case did not produce TwoWord values")
r, _ = experiment.residual(op, v, hi + lo, b)
norms = np.linalg.norm(b, axis=0)
normalized = np.abs(r) / norms
row, column = np.unravel_index(np.argmax(normalized), normalized.shape)
print(
    json.dumps(
        {
            "worst_row": int(row),
            "column": int(column),
            "row_normalized": float(normalized[row, column]),
            "column_relative": report["float64_projection_compensated_relative"],
        }
    ),
    flush=True,
)
if normalized[row, column] <= 1e-11:
    raise RuntimeError("No single violating row; do not invent a synthetic witness")
a = experiment.block(op, v, int(row), int(row) + 1)[0]


def pairs(values):
    return [[float(z.real), float(z.imag)] for z in values]


original_record = rd_paths.RECORDS / "uvlm-column-update-result.json"
original = json.loads(original_record.read_text())
original_system = next(item for item in original["systems"] if item["m"] == m and item["k"] == k)
actual_hashes = [experiment.effective_hash(op, v), experiment.previous.hash_array(b)]
if actual_hashes != original_system["input_sha256"]:
    raise RuntimeError("Regenerated coefficients or rhs differ from the frozen original system")
sources = []
for module in list(sys.modules.values()):
    filename = getattr(module, "__file__", None)
    if (
        filename
        and Path(filename).suffix == ".py"
        and (Path.cwd() / "rd") in Path(filename).parents
    ):
        path = Path(filename)
        sources.append(
            dict(
                id=path.stem,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                hint=path.relative_to(Path.cwd()).as_posix(),
            )
        )
result = dict(
    m=m,
    k=k,
    row=int(row),
    column=int(column),
    mode=experiment.O.MODES[column]["name"],
    original_record=dict(
        id="uvlm-column-update-result",
        sha256=hashlib.sha256(original_record.read_bytes()).hexdigest(),
        hint=original_record.relative_to(Path.cwd()).as_posix(),
    ),
    sources=sorted(sources, key=lambda item: item["id"]),
    input_hashes_equal=True,
    projection_row_relative=float(normalized[row, column]),
    process_peak_bytes=experiment.O.P.memory.memory_bytes(),
    elapsed_seconds=experiment.time.monotonic() - experiment.START,
    coefficient_sha256=experiment.effective_hash(op, v),
    rhs_sha256=experiment.previous.hash_array(b),
    scale=float(norms[column]),
    tolerance=1e-11,
    a=pairs(a),
    b=pairs([b[row, column]]),
    hi=pairs(hi[:, column]),
    lo=pairs(lo[:, column]),
    selected_accuracy=report["selected_accuracy"],
    projection_relative=report["float64_projection_compensated_relative"],
    numpy=np.__version__,
    scipy=scipy.__version__,
    scope=("One original-system row is a rejection witness; "
           "not a product full-system acceptance result"),
)
Path("witness.json").write_text(
    json.dumps(result, indent=2, allow_nan=False) + "\n", newline="\n"
)
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not args.probe:
        if args.output is None:
            parser.error("--output is required unless --probe")
        if args.output.exists():
            raise FileExistsError("Refusing to overwrite existing witness")
    root = Path(__file__).resolve().parents[2]
    with tempfile.TemporaryDirectory(prefix="uvlas-precision-witness-") as temporary:
        work = Path(temporary)
        shutil.copy2(root / "AGENTS.md", work / "AGENTS.md")
        shutil.copytree(root / "rd", work / "rd", ignore=shutil.ignore_patterns("__pycache__"))
        command = [
            "uv",
            "run",
            "--no-project",
            "--python",
            "3.12",
            "--with",
            "numpy==2.3.5",
            "--with",
            "scipy==1.16.3",
            "--with",
            "scikit-fem==12.0.2",
            "python",
            "-c",
            PAYLOAD,
        ]
        if args.probe:
            command.append("--probe")
        subprocess.run(command, cwd=work, check=True, timeout=30 if args.probe else 320)
        if not args.probe:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with (work / "witness.json").open("rb") as source, args.output.open("xb") as target:
                shutil.copyfileobj(source, target)
            print(json.dumps({"sha256": hashlib.sha256(args.output.read_bytes()).hexdigest()}))


if __name__ == "__main__":
    main()
