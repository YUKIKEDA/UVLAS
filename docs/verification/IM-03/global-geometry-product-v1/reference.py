"""Extract only frozen literal inputs; validate them with the independent Fraction oracle."""

import argparse
import ast
import hashlib
import importlib.util
import json
from fractions import Fraction
from pathlib import Path
from typing import Any


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Refusing to overwrite a fixture")
    root = Path(__file__).resolve().parents[4]
    source = root / "docs/verification/IM-03/global-geometry-te-start-v1/oracle.py"
    expected_hash = "854c1c13c02e1f15d53c08e57818328d5769fc9daf6dd1c1e9d827e16cde3f69"
    if hashlib.sha256(source.read_bytes()).hexdigest() != expected_hash:
        raise SystemExit("Frozen independent source hash mismatch")
    spec = importlib.util.spec_from_file_location("independent_geometry", source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    main_body = next(
        n.body
        for n in ast.parse(source.read_text(encoding="utf-8")).body
        if isinstance(n, ast.FunctionDef) and n.name == "main"
    )
    literals: list[ast.stmt] = [
        n
        for n in main_body
        if isinstance(n, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id in {"a", "cases"} for t in n.targets)
    ]
    assert len(literals) == 2
    namespace: dict[str, Any] = {"Q": Fraction}
    exec(compile(ast.Module(body=literals, type_ignores=[]), str(source), "exec"), namespace)
    a = namespace["a"]
    records: list[dict[str, object]] = []
    for name, b, expected, points, edge, accepted in namespace["cases"]:
        exact_b = tuple(tuple(map(Fraction, p)) for p in b)
        expected_vertices = {tuple(map(Fraction, p)) for p in expected}
        assert module.intersection(a, exact_b) == expected_vertices
        assert module.allowed(expected_vertices, points, edge) == accepted
        records.append(
            {
                "name": name,
                "a": a,
                "b": exact_b,
                "allowed_points": points,
                "allowed_edge": edge,
                "accepted": accepted,
                "vertices": sorted(expected_vertices),
            }
        )

    def encode(value: object) -> float:
        if isinstance(value, Fraction):
            result = float(value)
            assert Fraction(result) == value
            return result
        raise TypeError(type(value))

    payload = {
        "source": str(source.relative_to(root)).replace("\\", "/"),
        "source_sha256": expected_hash,
        "tolerance": "0",
        "cases": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, default=encode) + "\n", encoding="utf-8", newline="\n"
    )
    print("12 frozen literal cases independently rechecked; exact binary64 coordinates")


if __name__ == "__main__":
    main()
