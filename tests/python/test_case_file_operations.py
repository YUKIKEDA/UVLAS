"""IO4-01/02/05: verify parity of the product CLI and Python file entries."""

from __future__ import annotations

import copy
import hashlib
import importlib
import json
import runpy
import subprocess
import sys
import tempfile
import unittest
from decimal import Decimal
from pathlib import Path

sys.path.insert(0, str(Path(sys.argv.pop(1)).resolve()))
CLI = Path(sys.argv.pop(1)).resolve()
EXAMPLES = json.loads(Path(sys.argv.pop(1)).read_text(encoding="utf-8"))
uvlas = importlib.import_module("_uvlas")


def cli(*args, directory=None):
    return subprocess.run(
        [str(CLI), *map(str, args)], cwd=directory, capture_output=True, timeout=10, check=False
    )


def json_bytes(value):
    return json.dumps(value, ensure_ascii=False).encode("utf-8")


def expected_dictionary_source(value):
    if isinstance(value, dict):
        members = [
            json.dumps(key, ensure_ascii=False) + ":" + expected_dictionary_source(item)
            for key, item in value.items()
        ]
        return "{" + ",".join(members) + "}"
    if isinstance(value, list):
        return "[" + ",".join(map(expected_dictionary_source, value)) + "]"
    if isinstance(value, Decimal):
        return str(value)
    if isinstance(value, float):
        number = Decimal(value)
        sign, digits, exponent = number.as_tuple()
        if exponent >= 0:
            return format(number, ".1f")
        return ("-" if sign else "") + "".join(map(str, digits)) + "e" + str(exponent)
    return json.dumps(value, ensure_ascii=False)


class CaseFileOperations(unittest.TestCase):
    def assert_failure(self, result, code):
        self.assertFalse(result.ok)
        self.assertEqual(result.failure.code, code)
        with self.assertRaises(ValueError):
            _ = result.value

    def test_validation_and_owned_identity_agree_at_the_file_entry(self):
        data = json_bytes(EXAMPLES["schema_examples"]["structural_static"]) + b"\r\n"
        with tempfile.TemporaryDirectory() as workspace:
            path = Path(workspace) / "検証 case.json"
            path.write_bytes(data)
            output = cli("validate", path.name, directory=workspace)
            self.assertEqual(output.returncode, 0, output.stderr)
            self.assertIn(b"Case validation succeeded.", output.stdout)
            loaded = uvlas.load_case(path)
            self.assertTrue(loaded.ok)
            parent = loaded.value
            del loaded
            self.assertEqual(parent.source_sha256, hashlib.sha256(data).digest())
            self.assertEqual(parent.source_byte_count, len(data))
            self.assertEqual(path.read_bytes(), data)
            del parent

            for malformed in (b'{"schema":', b'{"schema":"uvlas.case/1","schema":"uvlas.case/1"}'):
                with self.subTest(malformed=malformed):
                    path.write_bytes(malformed)
                    self.assertEqual(cli("validate", path).returncode, 2)
                    self.assert_failure(uvlas.load_case(path), uvlas.FailureCode.InvalidInput)

            path.write_bytes(b" " * (8 * 1024 * 1024 + 1))
            self.assertEqual(cli("validate", path).returncode, 5)
            self.assert_failure(uvlas.load_case(path), uvlas.FailureCode.ResourceExceeded)

    def test_reference_digest_precedes_role_validation(self):
        case = copy.deepcopy(EXAMPLES["schema_examples"]["structural_static"])
        reference = copy.deepcopy(EXAMPLES["fragments"]["ComparisonReference.file"])
        reference["source"].update(
            source_id="900", path="sample.csv", sha256=hashlib.sha256(b"abc").hexdigest()
        )
        reference.update(columns=["x", "y", "z"], uncertainty=[0, 0, 0])
        case["observations"][0]["comparison"]["reference"] = reference
        with tempfile.TemporaryDirectory() as workspace:
            path = Path(workspace) / "case.json"
            external = path.with_name("sample.csv")
            path.write_bytes(json_bytes(case))
            external.write_bytes(b"xyz")
            failed_digest = cli("validate", path)
            self.assertEqual(failed_digest.returncode, 2)
            self.assertIn(b"SHA-256 mismatch", failed_digest.stderr)
            self.assert_failure(uvlas.load_case(path), uvlas.FailureCode.InvalidInput)

            external.write_bytes(b"abc")
            self.assertIn(b"dependency validation is unsupported", cli("validate", path).stderr)
            self.assert_failure(uvlas.load_case(path), uvlas.FailureCode.Unsupported)
            self.assert_failure(uvlas.case_from_dict(case), uvlas.FailureCode.InvalidInput)
            self.assert_failure(
                uvlas.case_from_dict(case, reference_root=path.parent),
                uvlas.FailureCode.Unsupported,
            )

    def test_sidecar_roles_are_read_at_the_shared_case_entry(self):
        case = copy.deepcopy(EXAMPLES["schema_examples"]["coupled_transient"])
        geometry = copy.deepcopy(EXAMPLES["fragments"]["Geometry.nastran"])
        deck = b"BEGIN BULK\nGRID,1,,3.0,4.0,5.0\nENDDATA\n"
        geometry["deck"].update(
            source_id="300", path="skin.bdf", sha256=hashlib.sha256(deck).hexdigest()
        )
        geometry["sidecar"].update(source_id="301", path="skin.json")
        case["model"]["geometry"].append(geometry)
        sidecar = {
            "schema": "uvlas.nastran-sidecar/1",
            "geometry": geometry["id"],
            "generator": "MSC Nastran 2022.1",
            "units": {"length": "m", "mass": "kg", "time": "s", "angle": "deg"},
            **{key: [] for key in ("shells", "materials", "beams", "masses", "supports", "joints")},
        }
        with tempfile.TemporaryDirectory() as workspace:
            root = Path(workspace)
            path = root / "case.json"
            (root / "skin.bdf").write_bytes(deck)
            cases = {
                "matching geometry": json_bytes(sidecar),
                "malformed sidecar": b"{}",
                "different geometry": json_bytes({**sidecar, "geometry": "999"}),
            }
            for scenario, data in cases.items():
                with self.subTest(scenario=scenario):
                    geometry["sidecar"]["sha256"] = hashlib.sha256(data).hexdigest()
                    (root / "skin.json").write_bytes(data)
                    path.write_bytes(json_bytes(case))
                    result = uvlas.load_case(path)
                    other = uvlas.case_from_dict(case, reference_root=root)
                    output = cli("validate", path)
                    if scenario == "matching geometry":
                        self.assertTrue(result.ok)
                        self.assertTrue(other.ok)
                        self.assertEqual(output.returncode, 0, output.stderr)
                    else:
                        self.assert_failure(result, uvlas.FailureCode.InvalidInput)
                        self.assert_failure(other, uvlas.FailureCode.InvalidInput)
                        self.assertEqual(output.returncode, 2, output.stderr)
                        self.assertIn(b"Nastran sidecar format", output.stderr)

    def test_static_observation_references_across_input_routes(self):
        verification = Path(__file__).resolve().parents[2] / "docs/verification/IM-06"
        witness = runpy.run_path(str(verification / "check-static-observation-example.py"))
        read_reference = witness["read_json"]
        condition_digest = witness["condition"]
        fixture = Path(__file__).with_name("fixtures") / "static-observation"
        case_bytes = (verification / "static-observation-case.json").read_bytes()
        reference_bytes = (verification / "parametric-observations.json").read_bytes()
        diagnostics = b"/observations/0/comparison/reference/source"

        with tempfile.TemporaryDirectory() as workspace:
            root = Path(workspace)
            path = root / "case.json"
            reference_path = root / "parametric-observations.json"

            def acquisition(filename, source, role):
                data = (fixture / filename).read_bytes()
                (root / filename).write_bytes(data)
                return {
                    "source_id": source,
                    "path": filename,
                    "sha256": hashlib.sha256(data).hexdigest(),
                    "role": role,
                }

            geometries = {
                "parametric": None,
                "gmsh": {
                    "kind": "gmsh",
                    "id": "2",
                    "mesh": acquisition("plate.msh", "1", "gmsh-msh4.1-ascii"),
                    "sidecar": acquisition("plate-gmsh.json", "2", "gmsh-sidecar"),
                    "reference": {
                        "kind": "fit",
                        "degree": 1,
                        "distance_budget": 1e-12,
                        "normal_budget": 1e-12,
                        "conversion_observations": ["1"],
                    },
                },
                "nastran": {
                    "kind": "nastran",
                    "id": "3",
                    "deck": acquisition("plate.bdf", "1", "nastran-deck"),
                    "sidecar": acquisition("plate-nastran.json", "2", "nastran-sidecar"),
                    "dialect": "msc-2022.1",
                    "mode": "structural-extraction",
                    "excluded_cards": [],
                },
            }

            def publish(case, reference):
                data = expected_dictionary_source(reference).encode("utf-8")
                reference_path.write_bytes(data)
                digest = hashlib.sha256(data).hexdigest()
                for observation in case["observations"]:
                    observation["comparison"]["reference"]["source"]["sha256"] = digest
                path.write_bytes(json_bytes(case))

            for route, external in geometries.items():
                with self.subTest(route=route):
                    case = json.loads(case_bytes)
                    if external is not None:
                        case["model"]["geometry"].append(external)
                    reference = read_reference(reference_bytes)
                    reference["condition_sha256"] = condition_digest(
                        read_reference(json_bytes(case))
                    )
                    reference["observations"].reverse()
                    for observation in reference["observations"]:
                        observation["components"].reverse()
                    publish(case, reference)
                    output = cli("validate", path)
                    self.assertEqual(output.returncode, 0, output.stderr)
                    self.assertTrue(uvlas.load_case(path).ok)

                    dictionary_condition = condition_digest(
                        read_reference(expected_dictionary_source(case))
                    )
                    if dictionary_condition != reference["condition_sha256"]:
                        self.assert_failure(
                            uvlas.case_from_dict(case, reference_root=root),
                            uvlas.FailureCode.InvalidInput,
                        )
                    reference["condition_sha256"] = dictionary_condition
                    publish(case, reference)
                    self.assertTrue(uvlas.case_from_dict(case, reference_root=root).ok)
                    reference_path.write_bytes(reference_path.read_bytes() + b" ")
                    self.assert_failure(
                        uvlas.case_from_dict(case, reference_root=root),
                        uvlas.FailureCode.InvalidInput,
                    )
                    self.assertIn(b"SHA-256 mismatch", cli("validate", path).stderr)

                    reference["condition_sha256"] = condition_digest(
                        read_reference(json_bytes(case))
                    )
                    point = reference["observations"][1]
                    reaction = reference["observations"][0]
                    self.assertEqual(point["id"], "1")
                    self.assertEqual(reaction["id"], "2")
                    controls = [
                        (reference, "condition_sha256", "0" * 64),
                        (point["observation"]["point"], "side", "upper"),
                        (point["observation"], "basis", {"kind": "material"}),
                        (point["comparison"], "time_s", 1),
                        (point["comparison"], "time_side", "left"),
                        (point, "components", point["components"][1:]),
                        (point, "components", [*point["components"], point["components"][0]]),
                        (point["components"][0], "unit", "mm"),
                        (point["components"][0], "uncertainty", 1),
                        (point["components"][0], "uncertainty", -1),
                        (reaction["observation"], "origin_m", [2, 1, 4]),
                        (reaction["observation"], "support", "999"),
                    ]
                    for fields, field, value in controls:
                        with self.subTest(field=field, value=value):
                            original = fields[field]
                            fields[field] = value
                            publish(case, reference)
                            self.assert_failure(
                                uvlas.load_case(path), uvlas.FailureCode.InvalidInput
                            )
                            output = cli("validate", path)
                            self.assertEqual(output.returncode, 2, output.stderr)
                            self.assertIn(diagnostics, output.stderr)
                            fields[field] = original
                    reference["schema"] = "uvlas.observation-reference/2"
                    publish(case, reference)
                    self.assert_failure(uvlas.load_case(path), uvlas.FailureCode.Unsupported)

            case = json.loads(case_bytes)
            template = read_reference(reference_bytes)
            condition = condition_digest(read_reference(json_bytes(case)))
            for index, (row, record) in enumerate(
                zip(case["observations"], template["observations"], strict=True)
            ):
                source = row["comparison"]["reference"]["source"]
                name = f"reference-{index}.json"
                data = expected_dictionary_source(
                    {
                        "schema": template["schema"],
                        "condition_sha256": condition,
                        "observations": [record],
                    }
                ).encode("utf-8")
                (root / name).write_bytes(data)
                source.update(
                    path=name, source_id=str(11 + index), sha256=hashlib.sha256(data).hexdigest()
                )
            path.write_bytes(json_bytes(case))
            self.assertTrue(uvlas.load_case(path).ok)
            self.assertEqual(cli("validate", path).returncode, 0)
            second = root / "reference-1.json"
            record = read_reference(second.read_bytes())
            record["observations"][0]["observation"]["origin_m"][0] = Decimal(5)
            data = expected_dictionary_source(record).encode("utf-8")
            second.write_bytes(data)
            case["observations"][1]["comparison"]["reference"]["source"]["sha256"] = hashlib.sha256(
                data
            ).hexdigest()
            path.write_bytes(json_bytes(case))
            rejected = cli("validate", path)
            self.assertEqual(rejected.returncode, 2, rejected.stderr)
            self.assertIn(b"/observations/1/comparison/reference/source", rejected.stderr)
            self.assertIn(b"reference node", rejected.stderr)

    def test_dictionary_checks_schema_and_owns_the_copied_input(self):
        case = copy.deepcopy(EXAMPLES["schema_examples"]["structural_static"])
        result = uvlas.case_from_dict(case)
        self.assertTrue(result.ok)
        checked = result.value
        expected = expected_dictionary_source(case).encode("utf-8")
        with tempfile.TemporaryDirectory() as workspace:
            path = Path(workspace) / "expected.json"
            path.write_bytes(expected)
            file_case = uvlas.load_case(path)
            self.assertTrue(file_case.ok)
            self.assertEqual(checked.source_sha256, file_case.value.source_sha256)
            self.assertEqual(checked.source_byte_count, file_case.value.source_byte_count)
            del file_case
        case["schema"] = "uvlas.invalid/1"
        self.assert_failure(uvlas.case_from_dict(case), uvlas.FailureCode.Unsupported)
        self.assertEqual(checked.source_sha256, hashlib.sha256(expected).digest())
        self.assertEqual(checked.source_byte_count, len(expected))

    def test_dictionary_preserves_binary64_and_ignores_the_python_digit_setting(self):
        integer = 2**63
        case = copy.deepcopy(EXAMPLES["schema_examples"]["structural_static"])
        comparison = case["observations"][0]["comparison"]
        comparison["window"] = [float(integer), integer]
        result = uvlas.case_from_dict(case)
        self.assertTrue(result.ok)
        expected = expected_dictionary_source(case).encode("utf-8")
        self.assertEqual(result.value.source_sha256, hashlib.sha256(expected).digest())

        comparison["window"] = [10**700] * 2
        with tempfile.TemporaryDirectory() as workspace:
            path = Path(workspace) / "large.json"
            path.write_bytes(json_bytes(case))
            self.assert_failure(uvlas.load_case(path), uvlas.FailureCode.Nonfinite)
            limit = sys.get_int_max_str_digits()
            try:
                sys.set_int_max_str_digits(640)
                self.assert_failure(uvlas.case_from_dict(case), uvlas.FailureCode.Nonfinite)
            finally:
                sys.set_int_max_str_digits(limit)

    def test_dictionary_rejects_non_json_values_and_over_budget_input(self):
        cycle = []
        cycle.append(cycle)

        class CustomDictionary(dict):
            pass

        for value in (
            None,
            CustomDictionary(),
            {1: True},
            {"x": cycle},
            {"x": (1,)},
            {"x": "\ud800"},
            {"x": float("nan")},
            {"x": float("inf")},
        ):
            with self.subTest(value=repr(value)):
                self.assert_failure(uvlas.case_from_dict(value), uvlas.FailureCode.InvalidInput)

        nested = []
        for _ in range(128):
            nested = [nested]
        for value in (nested, "あ" * 22_000, 1 << 4_097, [None] * 262_144, ["x" * 65_536] * 128):
            with self.subTest(kind=type(value).__name__):
                self.assert_failure(
                    uvlas.case_from_dict({"x": value}), uvlas.FailureCode.ResourceExceeded
                )

    def test_inspection_returns_owned_original_bytes_and_rejects_disagreement(self):
        with tempfile.TemporaryDirectory() as workspace:
            directory = Path(workspace) / "結果 case"
            directory.mkdir()
            resolved = directory / "resolved-case.json"
            resolved.write_bytes(b'{"sample":"verified"}\r\n')
            provenance = json_bytes(
                {
                    "build": json.loads(cli("--build-info").stdout),
                    "original_case": {
                        "source_id": "0",
                        "byte_count": "6",
                        "sha256": hashlib.sha256(b"source").hexdigest(),
                    },
                    "resolved_case_sha256": hashlib.sha256(resolved.read_bytes()).hexdigest(),
                }
            )
            (directory / "provenance.json").write_bytes(provenance)
            raw = (
                json_bytes(
                    {
                        "completed": True,
                        "exit_code": 0,
                        "failure": None,
                        "last_accepted": None,
                        "provenance_sha256": hashlib.sha256(provenance).hexdigest(),
                    }
                )
                + b"\r\n"
            )
            (directory / "status.json").write_bytes(raw)
            output = cli("inspect", directory)
            self.assertEqual(output.returncode, 0, output.stderr)
            inspected = uvlas.inspect_result(directory)
            self.assertTrue(inspected.ok)
            verified = inspected.value
            self.assertEqual(verified, raw)
            self.assertEqual(output.stdout, raw)

            resolved.write_bytes(b"{}")
            self.assertEqual(cli("inspect", directory).returncode, 6)
            self.assert_failure(uvlas.inspect_result(directory), uvlas.FailureCode.InvalidInput)
            self.assertEqual(verified, raw)


if __name__ == "__main__":
    unittest.main()
