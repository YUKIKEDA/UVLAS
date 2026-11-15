"""Import and verify the frozen 24-state PG reference; never run its R&D generator."""

from __future__ import annotations

import argparse
import ast
import hashlib
import io
import json
import math
import struct
import zipfile
from decimal import Context, Decimal, localcontext
from pathlib import Path
from typing import NoReturn, cast

from tools.rd import data_store

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = Path("tests/numerics/fixtures/j1-pg-free-history-24-v1")
AUDIT = Path("rd/foundation/numerics/records/numerical-foundation-adoption-audit.json")
SOURCES = (
    "rd/coupled/j1/experiments/export_j1_numerical_oracles.py",
    "rd/coupled/j1/experiments/pg_uvlm_free_step.py",
)
ARTIFACTS = ("pg-free-step.txt", "pg-free-oracle.npz", "oracle.json")
LOGICAL = "build/analysis/numeric-oracles/"
DIMENSION = 24
SCHEMA = "uvlas.foundation-reference/1"
REFERENCE_ID = "j1-pg-free-history-24-v1"


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def object_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_constant(value: str) -> NoReturn:
    raise ValueError(f"Nonstandard JSON numeric constant: {value}")


def json_object(data: bytes) -> dict[str, object]:
    value: object = json.loads(data, object_pairs_hook=object_pairs, parse_constant=reject_constant)
    if not isinstance(value, dict):
        raise ValueError("Expected a JSON object")
    return cast("dict[str, object]", value)


def mapping(value: object) -> dict[str, object]:
    if not isinstance(value, dict):
        raise ValueError("Expected a string-keyed map")
    raw = cast("dict[object, object]", value)
    if not all(isinstance(key, str) for key in raw):
        raise ValueError("Expected string map keys")
    return cast("dict[str, object]", raw)


def number(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("Expected a finite numeric scalar")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Nonfinite scalar")
    return result


def read_array(
    archive: zipfile.ZipFile, name: str, dtype: str, shape: tuple[int, ...]
) -> list[float]:
    """Read only the exact little-endian, row-major NPY-v1 types of this frozen archive."""
    size = math.prod(shape) * (16 if dtype == "<c16" else 8)
    info = archive.getinfo(name)
    if info.file_size > size + 4096:
        raise ValueError("Unexpected NPY size")
    raw = archive.read(name)
    if raw[:8] != b"\x93NUMPY\x01\x00" or len(raw) < 10:
        raise ValueError("Expected NPY v1")
    header_size = struct.unpack("<H", raw[8:10])[0]
    header: object = ast.literal_eval(raw[10 : 10 + header_size].decode("ascii").strip())
    expected = {"descr": dtype, "fortran_order": False, "shape": shape}
    if header != expected or len(raw) != 10 + header_size + size:
        raise ValueError("Unexpected NPY dtype, shape, order or payload extent")
    values = [item[0] for item in struct.iter_unpack("<d", raw[10 + header_size :])]
    if not all(math.isfinite(value) for value in values):
        raise ValueError("Nonfinite NPY input")
    return values


def solve_decimal(matrix: list[float], rhs: list[float], precision: int) -> list[str]:
    """Independent scalar elimination of the stored binary64 inputs, without Eigen/UVLAS."""
    n = len(rhs)
    if len(matrix) != n * n:
        raise ValueError("Wrong reference system shape")
    with localcontext(Context(prec=precision)):
        rows = [
            [Decimal(value) for value in matrix[i * n : (i + 1) * n]] + [Decimal(rhs[i])]
            for i in range(n)
        ]
        for column in range(n):
            pivot = max(range(column, n), key=lambda row: abs(rows[row][column]))
            rows[column], rows[pivot] = rows[pivot], rows[column]
            if rows[column][column] == 0:
                raise ValueError("Singular independent reference system")
            for row in range(column + 1, n):
                ratio = rows[row][column] / rows[column][column]
                rows[row][column] = Decimal(0)
                for index in range(column + 1, n + 1):
                    rows[row][index] -= ratio * rows[column][index]
        result = [Decimal(0)] * n
        for row in reversed(range(n)):
            value = rows[row][n]
            for column in range(row + 1, n):
                value -= rows[row][column] * result[column]
            result[row] = value / rows[row][row]
        return [str(value) for value in result]


def materialize(cache_root: Path | None = None) -> tuple[dict[str, object], dict[str, bytes]]:
    audit_raw = (ROOT / AUDIT).read_bytes()
    audit = json_object(audit_raw)
    indexed = mapping(data_store.load_index()["files"])
    hashes = mapping(audit["artifacts"])
    originals: dict[str, bytes] = {}
    origins: list[dict[str, object]] = []
    for name in ARTIFACTS:
        logical = LOGICAL + name
        entry = mapping(indexed[logical])
        path = data_store.cached(logical) if cache_root is None else cache_root / logical
        raw = path.read_bytes()
        if digest(raw) != entry["sha256"] or entry["sha256"] != hashes[logical]:
            raise ValueError(f"Frozen artifact hash mismatch: {logical}")
        originals[name] = raw
        origins.append(
            {
                "id": Path(name).stem,
                "hint": logical,
                "sha256": digest(raw),
                "bundle": entry.get("bundle"),
            }
        )
    text = originals["pg-free-step.txt"].decode("ascii").split()
    if len(text) != 1 + DIMENSION * DIMENSION or text[0] != str(DIMENSION):
        raise ValueError("Wrong stored full-history matrix extent")
    step = [number(float(value)) for value in text[1:]]
    with zipfile.ZipFile(io.BytesIO(originals["pg-free-oracle.npz"])) as archive:
        if sorted(archive.namelist()) != ["K.npy", "M.npy", "P.npy", "eigenvalues.npy"]:
            raise ValueError("Unexpected or duplicate archive members")
        binary_step = read_array(archive, "P.npy", "<f8", (DIMENSION, DIMENSION))
        eigenvalues = read_array(archive, "eigenvalues.npy", "<c16", (DIMENSION,))
        read_array(archive, "M.npy", "<f8", (2, 2))
        read_array(archive, "K.npy", "<f8", (2, 2))
    if b"".join(struct.pack("<d", x) for x in step) != b"".join(
        struct.pack("<d", x) for x in binary_step
    ):
        raise ValueError("Text and NPY matrices differ in binary64 bits")
    metadata = json_object(originals["oracle.json"])
    if metadata["dimension"] != DIMENSION or metadata["sha256"] != digest(originals[ARTIFACTS[0]]):
        raise ValueError("Original oracle metadata is inconsistent")
    matrix = [
        (1.0 if i == j else 0.0) - step[i * DIMENSION + j]
        for i in range(DIMENSION)
        for j in range(DIMENSION)
    ]
    rhs = [float(i + 1) / DIMENSION for i in range(DIMENSION)]
    # A second precision checks the frozen digits without relying on a UVLAS result.
    reference = solve_decimal(matrix, rhs, 120)
    wider = solve_decimal(matrix, rhs, 160)
    with localcontext(Context(prec=160)):
        if any(
            abs(Decimal(a) - Decimal(b)) > Decimal("1e-100")
            for a, b in zip(reference, wider, strict=True)
        ):
            raise ValueError("Independent reference precision did not stabilize")
    source_hashes = mapping(audit["sources"])
    sources = {name: digest((ROOT / name).read_bytes()) for name in SOURCES}
    if any(value != source_hashes[name] for name, value in sources.items()):
        raise ValueError("Frozen generator source changed")
    fixture: dict[str, object] = {
        "schema": SCHEMA,
        "id": REFERENCE_ID,
        "source_record": {"id": AUDIT.stem, "hint": AUDIT.as_posix(), "sha256": digest(audit_raw)},
        "source_artifacts": origins,
        "sources": sources,
        "artifacts": {
            (FIXTURE / "source" / name).as_posix(): digest(raw) for name, raw in originals.items()
        },
        "layout": {
            "rows": DIMENSION,
            "columns": DIMENSION,
            "order": "row-major",
            "state_order": "logical-index-0-23",
            "units": "abstract-dimensionless-replay-only",
        },
        "step_values": step,
        "operator": "rounded-identity-minus-stored-P",
        "matrix": matrix,
        "rhs": rhs,
        "reference_solution": reference,
        "reference_algorithm": (
            "stdlib Decimal scalar partial-pivot elimination; 120 digits checked at 160 digits"
        ),
        "eigenvalues": [[eigenvalues[2 * i], eigenvalues[2 * i + 1]] for i in range(DIMENSION)],
        "original_spectral_radius": number(metadata["spectral_radius"]),
        "relative_solution_tolerance": 1e-12,
        "normalized_residual_tolerance": 1e-12,
    }
    return fixture, originals


def encode(fixture: dict[str, object]) -> bytes:
    return (json.dumps(fixture, ensure_ascii=False, allow_nan=False, indent=2) + "\n").encode(
        "utf-8"
    )


def import_fixture(output: Path, cache_root: Path | None = None) -> None:
    # Refuse before reading/decompressing anything, and reserve the directory exclusively.
    output.mkdir(parents=True, exist_ok=False)
    fixture, originals = materialize(cache_root)
    (output / "source").mkdir()
    for name, raw in originals.items():
        with (output / "source" / name).open("xb") as stream:
            stream.write(raw)
    # The manifest is the completion marker; a failed import never looks complete.
    with (output / "fixture.json").open("xb") as stream:
        stream.write(encode(fixture))


def verify_fixture(directory: Path) -> None:
    fixture = json_object((directory / "fixture.json").read_bytes())
    if set(fixture) != {
        "schema",
        "id",
        "source_record",
        "source_artifacts",
        "sources",
        "artifacts",
        "layout",
        "step_values",
        "operator",
        "matrix",
        "rhs",
        "reference_solution",
        "reference_algorithm",
        "eigenvalues",
        "original_spectral_radius",
        "relative_solution_tolerance",
        "normalized_residual_tolerance",
    }:
        raise ValueError("Unknown or missing fixture fields")

    artifacts = mapping(fixture["artifacts"])
    expected = {(FIXTURE / "source" / name).as_posix() for name in ARTIFACTS}
    if set(artifacts) != expected:
        raise ValueError("Unexpected artifact map")
    for name, sha in artifacts.items():
        local = Path(name).relative_to(FIXTURE)
        if digest((directory / local).read_bytes()) != sha:
            raise ValueError("Fixture artifact changed")
    if fixture.get("schema") != SCHEMA or fixture.get("id") != REFERENCE_ID:
        raise ValueError("Unknown fixture schema or identity")
    layout = mapping(fixture["layout"])
    if layout != {
        "rows": DIMENSION,
        "columns": DIMENSION,
        "order": "row-major",
        "state_order": "logical-index-0-23",
        "units": "abstract-dimensionless-replay-only",
    }:
        raise ValueError("Unsupported layout or units")
    vectors: dict[str, list[float]] = {}
    for key, extent in (
        ("matrix", DIMENSION * DIMENSION),
        ("step_values", DIMENSION * DIMENSION),
        ("rhs", DIMENSION),
    ):
        values = fixture[key]
        if not isinstance(values, list):
            raise ValueError("Wrong fixture array")
        items = cast("list[object]", values)
        if len(items) != extent:
            raise ValueError("Wrong fixture extent")
        vectors[key] = [number(value) for value in items]
    matrix, step = vectors["matrix"], vectors["step_values"]
    tokens = (directory / "source" / "pg-free-step.txt").read_text(encoding="ascii").split()
    if len(tokens) != 1 + DIMENSION * DIMENSION or tokens[0] != str(DIMENSION):
        raise ValueError("Wrong original matrix shape")
    if any(
        struct.pack("<d", value) != struct.pack("<d", float(token))
        for value, token in zip(step, tokens[1:], strict=True)
    ):
        raise ValueError("JSON step differs from the frozen source")
    if (
        number(fixture["relative_solution_tolerance"]) != 1e-12
        or number(fixture["normalized_residual_tolerance"]) != 1e-12
    ):
        raise ValueError("Fixed tolerances changed")

    if fixture["operator"] != "rounded-identity-minus-stored-P" or any(
        value != (1.0 if i // DIMENSION == i % DIMENSION else 0.0) - step[i]
        for i, value in enumerate(matrix)
    ):
        raise ValueError("Operator adaptation changed")
    strings = fixture["reference_solution"]
    if not isinstance(strings, list):
        raise ValueError("Wrong independent reference array")
    items = cast("list[object]", strings)
    if len(items) != DIMENSION or not all(isinstance(value, str) for value in items):
        raise ValueError("Wrong independent reference shape")
    solution = [Decimal(value) for value in cast("list[str]", strings)]
    if not all(value.is_finite() for value in solution):
        raise ValueError("Nonfinite reference")
    # Validate the frozen oracle by its residual; do not regenerate or solve a new reference.
    with localcontext(Context(prec=160)):
        for row in range(DIMENSION):
            residual = Decimal(vectors["rhs"][row])
            for column in range(DIMENSION):
                residual -= Decimal(matrix[row * DIMENSION + column]) * solution[column]
            if abs(residual) > Decimal("1e-100"):
                raise ValueError("Frozen independent reference residual failed")
    with zipfile.ZipFile(directory / "source" / "pg-free-oracle.npz") as archive:
        frozen = read_array(archive, "eigenvalues.npy", "<c16", (DIMENSION,))
    pairs = fixture["eigenvalues"]
    if not isinstance(pairs, list):
        raise ValueError("Expected spectral pairs")
    raw_pairs = cast("list[object]", pairs)
    if len(raw_pairs) != DIMENSION:
        raise ValueError("Wrong spectral extent")
    eigenvalues: list[float] = []
    for pair in raw_pairs:
        if not isinstance(pair, list):
            raise ValueError("Expected a real/imaginary pair")
        components = cast("list[object]", pair)
        if len(components) != 2:
            raise ValueError("Wrong spectral pair extent")
        eigenvalues.extend(number(component) for component in components)
    if any(
        struct.pack("<d", a) != struct.pack("<d", b)
        for a, b in zip(eigenvalues, frozen, strict=True)
    ):
        raise ValueError("Frozen spectral oracle changed")
    metadata = json_object((directory / "source" / "oracle.json").read_bytes())
    if number(fixture["original_spectral_radius"]) != number(metadata["spectral_radius"]):
        raise ValueError("Original spectral metadata changed")
    audit = json_object((ROOT / AUDIT).read_bytes())
    indexed = mapping(data_store.load_index()["files"])
    origins: list[dict[str, object]] = []
    hashes = mapping(audit["artifacts"])
    for name in ARTIFACTS:
        logical = LOGICAL + name
        entry = mapping(indexed[logical])
        sha = digest((directory / "source" / name).read_bytes())
        if sha != entry["sha256"] or sha != hashes[logical]:
            raise ValueError("Packaged source does not match the frozen index")
        origins.append(
            {"id": Path(name).stem, "hint": logical, "sha256": sha, "bundle": entry.get("bundle")}
        )
    if fixture["source_artifacts"] != origins:
        raise ValueError("Source identities or provenance changed")
    provenance = mapping(fixture["source_record"])
    if provenance != {
        "id": AUDIT.stem,
        "hint": AUDIT.as_posix(),
        "sha256": digest((ROOT / AUDIT).read_bytes()),
    }:
        raise ValueError("Source record changed")
    sources = mapping(fixture["sources"])
    frozen_sources = mapping(audit["sources"])
    if set(sources) != set(SOURCES) or any(
        digest((ROOT / name).read_bytes()) != sha or sha != frozen_sources[name]
        for name, sha in sources.items()
    ):
        raise ValueError("Frozen source provenance changed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    create = commands.add_parser("import")
    create.add_argument("--output", type=Path, required=True)
    create.add_argument("--cache-root", type=Path)
    check = commands.add_parser("verify")
    check.add_argument("fixture", type=Path)
    args = parser.parse_args()
    if args.command == "import":
        import_fixture(args.output, args.cache_root)
    else:
        verify_fixture(args.fixture)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
