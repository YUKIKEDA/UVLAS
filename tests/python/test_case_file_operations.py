"""IO4-01/02/05: verify parity of the product CLI and Python file entries."""

from __future__ import annotations

import copy
import hashlib
import importlib
import json
import subprocess
import sys
import tempfile
import unittest
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

    def test_dictionary_checks_schema_and_owns_the_copied_input(self):
        case = copy.deepcopy(EXAMPLES["schema_examples"]["structural_static"])
        result = uvlas.case_from_dict(case)
        self.assertTrue(result.ok)
        checked = result.value
        expected = json.dumps(case, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        with tempfile.TemporaryDirectory() as workspace:
            path = Path(workspace) / "expected.json"
            path.write_bytes(expected)
            file_case = uvlas.load_case(path)
            self.assertTrue(file_case.ok)
            self.assertEqual(checked.source_sha256, file_case.value.source_sha256)
            self.assertEqual(checked.source_byte_count, file_case.value.source_byte_count)
        case["schema"] = "uvlas.invalid/1"
        self.assert_failure(uvlas.case_from_dict(case), uvlas.FailureCode.InvalidInput)
        self.assertEqual(checked.source_sha256, hashlib.sha256(expected).digest())
        self.assertEqual(checked.source_byte_count, len(expected))

    def test_dictionary_rejects_non_json_values_and_over_budget_input(self):
        cycle = []
        cycle.append(cycle)

        class CustomDictionary(dict):
            pass

        for value in (None, CustomDictionary(), {1: True}, {"x": cycle}, {"x": (1,)},
                      {"x": "\ud800"}, {"x": float("nan")}, {"x": float("inf")}):
            with self.subTest(value=repr(value)):
                self.assert_failure(uvlas.case_from_dict(value), uvlas.FailureCode.InvalidInput)

        nested = []
        for _ in range(128):
            nested = [nested]
        for value in (nested, "あ" * 22_000, 1 << 4_097, [None] * 262_144,
                      ["x" * 65_536] * 128):
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
            provenance = json_bytes({
                "build": json.loads(cli("--build-info").stdout),
                "original_case": {
                    "source_id": "0",
                    "byte_count": "6",
                    "sha256": hashlib.sha256(b"source").hexdigest(),
                },
                "resolved_case_sha256": hashlib.sha256(resolved.read_bytes()).hexdigest(),
            })
            (directory / "provenance.json").write_bytes(provenance)
            raw = json_bytes({
                "completed": True, "exit_code": 0, "failure": None, "last_accepted": None,
                "provenance_sha256": hashlib.sha256(provenance).hexdigest(),
            }) + b"\r\n"
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
