"""Reproduce the V-PG reference geometry and fixed support Cases from frozen R4/R5."""

from __future__ import annotations

import argparse
import ast
import json
from decimal import Decimal
from fractions import Fraction
from math import comb
from pathlib import Path
from typing import cast

from tools.qa.foundation_reference import (
    digest,
    json_object,
    mapping,
    object_pairs,
    reject_constant,
)

ROOT = Path(__file__).resolve().parents[2]
SHELL = Path("rd/structure/main-shell")
PAIRED = SHELL / "paired-geometry"
RECORDS = PAIRED / "records"
EXPERIMENTS = PAIRED / "experiments"
R4 = RECORDS / "paired-reference-r4-audit.json"
R5 = RECORDS / "paired-shell-r5-audit.json"
TEMPLATE = Path("docs/specs/foundation/case-schema-examples.json")
QUALITY = RECORDS / "mass-paired-quality-result.json"
EXPLICIT = RECORDS / "mass-paired-coefficient-input-result.json"
CUT = RECORDS / "paired-shell-r5-domain-v2-result.json"
CUT_SETUP = EXPERIMENTS / "check_paired_shell_r5_domain_v2.py"
INTERFACE = RECORDS / "mass-paired-interface-result.json"
INTERFACE_SETUP = EXPERIMENTS / "check_mass_paired_interface.py"
PINNED = {
    R4: "2af9e42f35928cb0d213f126fd481d03703ce6995b3c52e75ea1c8045cd9f0c9",
    R5: "48e4a9ecbcb3ea0ee6cff0269f71c0567032d9b0ddf6040421aba51f8e222167",
    TEMPLATE: "c957d97d4cf5eb2b919edd26cae004c85566f23b2d59dba306e3c01302e7dbfe",
}
COMMON = (
    SHELL / "polynomial-trace/experiments/polynomial_reference_surface.py",
    SHELL / "geometry-compatible/experiments/coefficient_compatible_finite.py",
    EXPERIMENTS / "exact_paired_finite.py",
)
SQUARE = tuple((Fraction(s), Fraction(t)) for s, t in ((0, 0), (1, 0), (1, 1), (0, 1)))
TAPER = (SQUARE[0], (Fraction(1), Fraction(3, 10)), (Fraction(1), Fraction(9, 10)), SQUARE[3])

Point = tuple[Fraction, Fraction]
Polynomial = dict[tuple[int, int], Fraction]


def array(value: object) -> list[object]:
    if not isinstance(value, list):
        raise ValueError("Expected an array")
    return cast("list[object]", value)


def exact(value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Decimal)):
        raise ValueError("Expected a finite JSON decimal")
    return Fraction(value)


def to_json(value: object) -> bytes:
    return (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode()


def rational(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def wire_point(point: Point) -> list[dict[str, str]]:
    return [rational(value) for value in point]


def decimal_number(value: Fraction) -> float:
    result = float(value)
    if Fraction(str(result)) != value:
        raise ValueError(f"Cannot serialize exact terminating Case decimal: {value}")
    return result


def pinned(path: Path, sha: str) -> bytes:
    data = (ROOT / path).read_bytes()
    if digest(data) != sha:
        raise ValueError(f"Frozen input SHA-256 mismatch: {path.as_posix()}")
    return data


def audited(audit: Path, paths: tuple[Path, ...], hashes: dict[str, str]) -> dict[Path, bytes]:
    sources = mapping(json_object(pinned(audit, PINNED[audit]))["sources"])
    result: dict[Path, bytes] = {}
    for path in paths:
        key = path.as_posix()
        sha = sources[key]
        if not isinstance(sha, str):
            raise ValueError(f"Missing source digest for {key}")
        result[path] = pinned(path, sha)
        hashes[key] = sha
    return result


def record(data: bytes) -> dict[str, object]:
    value: object = json.loads(
        data, parse_float=Decimal, object_pairs_hook=object_pairs, parse_constant=reject_constant
    )
    result = mapping(value)
    if result["passed"] is not True:
        raise ValueError("Frozen qualification did not pass")
    return result


def explicit_map(value: object) -> Polynomial:
    result: Polynomial = {}
    for key, scalar in mapping(value).items():
        powers: object = ast.literal_eval(key)
        if not isinstance(powers, tuple) or not isinstance(scalar, str):
            raise ValueError("Invalid explicit R4 coefficient")
        pair = cast("tuple[object, ...]", powers)
        if len(pair) != 2:
            raise ValueError("Invalid explicit R4 powers")
        i, j = pair
        if type(i) is not int or type(j) is not int:
            raise ValueError("Invalid explicit R4 powers")
        if min(i, j) < 0 or i + j > 4 or (i, j) in result:
            raise ValueError("Invalid R4 monomial")
        result[i, j] = Fraction(scalar)
    return result


def cut_parameters(source: bytes) -> list[Fraction]:
    # Read the pinned model constructor without importing the R&D program.
    matches = [
        node.value
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.keyword) and node.arg == "parameters"
    ]
    if len(matches) != 1 or not isinstance(matches[0], ast.Tuple):
        raise ValueError("Ambiguous swept-taper parameter tuple")
    return [Fraction(ast.unparse(value)) for value in matches[0].elts]


def graph(parameters: list[Fraction], chart: tuple[Fraction, Fraction, Fraction]) -> Polynomial:
    x0, width, skew = chart
    kx, ky, tau, cx, cy, mix = (parameters[:6] + [Fraction(0)] * 6)[:6]
    terms = (
        (2, 0, kx / 2),
        (0, 2, ky / 2),
        (1, 1, tau),
        (3, 0, cx / 3),
        (0, 3, cy / 3),
        (2, 2, mix),
    )
    result: Polynomial = {}
    for i, j, coefficient in terms:
        for s_power in range(i + 1):
            for t_power in range(i - s_power + 1):
                key = s_power, j + t_power
                remainder = i - s_power - t_power
                factor = comb(i, s_power) * comb(i - s_power, t_power)
                result[key] = result.get(key, Fraction(0)) + (
                    coefficient * factor * x0**remainder * width**s_power * skew**t_power
                )
    return {powers: value for powers, value in result.items() if value}


def make_case(
    template: bytes,
    *,
    n: int,
    h: Fraction,
    variable: bool,
    chart: tuple[Fraction, Fraction, Fraction],
    coefficients: Polynomial,
    domain: tuple[Point, ...],
) -> dict[str, object]:
    case = mapping(mapping(json_object(template)["schema_examples"])["structural_static"])
    model = mapping(case["model"])
    patch = mapping(array(model["geometry"])[0])
    x0, width, skew = chart
    patch.update(
        n=n,
        x0=rational(x0),
        width=rational(width),
        skew=rational(skew),
        degree=[max((powers[axis] for powers in coefficients), default=0) for axis in (0, 1)],
        coefficients=[
            {"powers": list(powers), "coefficient": rational(value)}
            for powers, value in sorted(coefficients.items())
        ],
        domain=[wire_point(point) for point in domain],
        boundaries=[
            {"name": "clamped_s0", "vertices": [wire_point(p) for p in (domain[0], domain[3])]}
        ],
    )
    mapping(patch["thickness"])["values"] = [
        decimal_number(h * (1 + s / 10 + t / 20) if variable else h) for s, t in SQUARE
    ]
    mapping(mapping(case["numerics"])["patch"])["maximum_partition"] = n
    # Frozen unit-density normalization and Lamé constants for E=1, nu=0.3.
    material = mapping(array(model["materials"])[0])
    material.update(
        poisson=0.3,
        provenance={
            "description": "凍結V-PGの単位密度と等方弾性定数",
            "references": [path.as_posix() for path in COMMON[1:]],
            "assumptions": [],
        },
    )
    return case


def interface_case(
    template: bytes, left_n: int, right_n: int, h: Fraction, parameters: list[Fraction]
) -> dict[str, object]:
    half = Fraction(1, 2)
    left_chart = (Fraction(0), half, Fraction(0))
    right_chart = (half, half, Fraction(0))
    case = make_case(
        template,
        n=left_n,
        h=h,
        variable=False,
        chart=left_chart,
        coefficients=graph(parameters, left_chart),
        domain=SQUARE,
    )
    right_case = make_case(
        template,
        n=right_n,
        h=h,
        variable=False,
        chart=right_chart,
        coefficients=graph(parameters, right_chart),
        domain=SQUARE,
    )
    model = mapping(case["model"])
    right_model = mapping(right_case["model"])
    geometry = array(model["geometry"])
    left = mapping(geometry[0])
    array(left["boundaries"]).append(
        {"name": "interface", "vertices": [wire_point(SQUARE[i]) for i in (1, 2)]}
    )
    right = mapping(array(right_model["geometry"])[0])
    right["id"] = "2"
    right_boundary = mapping(array(right["boundaries"])[0])
    right_boundary.update(name="interface", vertices=[wire_point(SQUARE[i]) for i in (3, 0)])
    geometry.append(right)
    right_shell = mapping(array(right_model["shells"])[0])
    right_shell.update(id="2", geometry="2")
    mapping(array(right_shell["layers"])[0])["id"] = "2"
    array(model["shells"]).append(right_shell)
    mapping(mapping(case["numerics"])["patch"])["maximum_partition"] = max(left_n, right_n)

    t = {"powers": [1, 0], "coefficient": rational(Fraction(1))}
    zeta = {"powers": [0, 1], "coefficient": rational(Fraction(1))}
    one = {"powers": [0, 0], "coefficient": rational(Fraction(1))}
    segments = 2 * max(left_n, right_n)
    maps: list[dict[str, object]] = []
    for index in range(segments):
        lo, hi = Fraction(index, segments), Fraction(index + 1, segments)
        maps.append(
            {
                "first_cell": [2 * left_n - 1, index * 2 * left_n // segments],
                "second_cell": [0, index * 2 * right_n // segments],
                "domain": [
                    wire_point((u, Fraction(z))) for u, z in ((lo, -1), (hi, -1), (hi, 1), (lo, 1))
                ],
                "first": [[one], [t], [zeta]],
                "second": [[], [t], [zeta]],
                "orientation": "reversed",
            }
        )
    model["joints"] = [
        {
            "kind": "adhesive",
            "id": "1",
            "first": {
                "target": {"kind": "shell", "id": "1"},
                "kind": "edge",
                "geometry": "1",
                "boundary": 1,
                "side": "interior",
            },
            "second": {
                "target": {"kind": "shell", "id": "2"},
                "kind": "edge",
                "geometry": "2",
                "boundary": 0,
                "side": "interior",
            },
            "mapping": maps,
            "split_at_attachment": False,
        }
    ]
    return case


def interface_input(source: bytes) -> tuple[list[Fraction], Fraction]:
    nodes = [
        node
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "prepare"
    ]
    if len(nodes) != 1 or len(nodes[0].args) != 5 or not isinstance(nodes[0].args[1], ast.Tuple):
        raise ValueError("Ambiguous frozen interface input")
    parts = nodes[0].args[1].elts
    if (
        len(parts) != 8
        or ast.unparse(parts[6]) != "offset"
        or Fraction(ast.unparse(parts[7])) != Fraction(1, 2)
    ):
        raise ValueError("Unexpected frozen interface chart")
    return [Fraction(ast.unparse(part)) for part in parts[:6]], Fraction(
        ast.unparse(nodes[0].args[2])
    )


def bundle() -> dict[Path, bytes]:
    hashes = {path.as_posix(): sha for path, sha in PINNED.items()}
    sources = audited(
        R4,
        (
            *COMMON,
            QUALITY,
            EXPLICIT,
            INTERFACE,
            INTERFACE_SETUP,
            EXPERIMENTS / "check_mass_paired_quality.py",
            EXPERIMENTS / "check_mass_paired_coefficient_input.py",
        ),
        hashes,
    )
    sources.update(audited(R5, (CUT, CUT_SETUP, EXPERIMENTS / "swept_taper_domain.py"), hashes))
    template = pinned(TEMPLATE, PINNED[TEMPLATE])
    explicit = record(sources[EXPLICIT])
    explicit_coefficients = explicit_map(explicit["coefficient_map"])
    swept_parameters = cut_parameters(sources[CUT_SETUP])
    files: dict[Path, bytes] = {}
    entries: list[dict[str, object]] = []
    for family, path, frozen, count in (
        ("r4-quality", QUALITY, record(sources[QUALITY]), 24),
        ("r4-explicit", EXPLICIT, explicit, 4),
        ("r5-taper", CUT, record(sources[CUT]), 4),
    ):
        rows = array(frozen["records"])
        if len(rows) != count:
            raise ValueError(f"Unexpected frozen series size: {family}")
        for index, item in enumerate(rows):
            row = mapping(item)
            n_value = row["n"]
            if type(n_value) is not int or n_value not in (1, 2, 4):
                raise ValueError("Invalid partition in frozen record")
            h = exact(row["h"])
            variable = row["variable"] if path == QUALITY else True
            if not isinstance(variable, bool):
                raise ValueError("Invalid thickness flag")
            parameters: list[Fraction] = []
            if path == QUALITY:
                parameters = [exact(x) for x in array(row["parameters"])]
            if path == CUT:
                parameters = swept_parameters
            chart = (Fraction(0), Fraction(1), Fraction(0))
            if len(parameters) == 9:
                chart = parameters[6], parameters[7], parameters[8]
            coefficients = explicit_coefficients if path == EXPLICIT else graph(parameters, chart)
            domain = TAPER if path == CUT else SQUARE
            target = Path(f"{family}-{index:02d}") / "case.json"
            data = to_json(
                make_case(
                    template,
                    n=n_value,
                    h=h,
                    variable=variable,
                    chart=chart,
                    coefficients=coefficients,
                    domain=domain,
                )
            )
            files[target] = data
            entries.append(
                {
                    "case": target.as_posix(),
                    "sha256": digest(data),
                    "result": path.as_posix(),
                    "record": index,
                    "n": n_value,
                    "h": rational(h),
                    "variable": variable,
                    "graph_parameters": [rational(x) for x in parameters],
                }
            )
    interface_parameters, thickness = interface_input(sources[INTERFACE_SETUP])
    interfaces = array(record(sources[INTERFACE])["records"])
    if len(interfaces) != 2:
        raise ValueError("Unexpected frozen interface series size")
    for index, item in enumerate(interfaces):
        partitions = array(mapping(item)["elements"])
        if partitions != [1, index + 1] or any(type(value) is not int for value in partitions):
            raise ValueError("Unexpected frozen interface partitions")
        left_n, right_n = cast("list[int]", partitions)
        target = Path(f"r4-interface-{index:02d}") / "case.json"
        data = to_json(interface_case(template, left_n, right_n, thickness, interface_parameters))
        files[target] = data
        entries.append(
            {
                "case": target.as_posix(),
                "sha256": digest(data),
                "result": INTERFACE.as_posix(),
                "record": index,
                "n": partitions,
                "h": rational(thickness),
                "variable": False,
                "graph_parameters": [rational(x) for x in interface_parameters],
            }
        )
    files[Path("manifest.json")] = to_json(
        {
            "schema": "uvlas.vpg-case-reference/1",
            "sources": dict(sorted(hashes.items())),
            "cases": entries,
        }
    )
    return files


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("generate", "verify"))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = cast("Path", args.output)
    files = bundle()
    if args.action == "generate":
        output.mkdir()
        for path, data in files.items():
            target = output / path
            target.parent.mkdir(exist_ok=True)
            with target.open("xb") as stream:
                stream.write(data)
    else:
        actual = {path.relative_to(output) for path in output.rglob("*") if path.is_file()}
        if actual != files.keys():
            raise ValueError("Unexpected fixture files")
        for path, data in files.items():
            if (output / path).read_bytes() != data:
                raise ValueError(f"Generated fixture differs: {path.as_posix()}")
    print(f"V-PG: {len(files) - 1} Case files, {digest(files[Path('manifest.json')])}")


if __name__ == "__main__":
    main()
