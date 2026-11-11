from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr
import hashlib
import io
import json
import tempfile
import unittest
from check_evidence_manifest import check, main


class EvidenceChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="solver-evidence-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve() / "project"
        self.root.mkdir()

    def write(self, name, value):
        p = self.root / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(value) if isinstance(value, (dict, list)) else value,
                     encoding="utf-8")
        return hashlib.sha256(p.read_bytes()).hexdigest()

    def test_transitive_diamond_and_false_scientific_status(self):
        leaf = self.write("data.bin", "original")
        child = self.write("child.json", {"passed": False, "sources": {"data.bin": leaf}})
        self.write("root.json", {"sources": {"child.json": child, "data.bin": leaf}})
        r = check(self.root, ["root.json", "child.json"])
        self.assertTrue(r["integrity_passed"])
        self.assertEqual(r["checked_files"], 3)
        self.assertEqual(r["scientific_acceptance"], "not_evaluated")
        self.assertIn({"manifest": "child.json", "recorded_status": {"passed": False}},
                      r["manifest_statuses"])

    def test_modified_leaf_fails_even_when_result_says_passed(self):
        h = self.write("solver.cpp", "first")
        self.write("r.json", {"passed": True, "sources": {"solver.cpp": h}})
        self.write("solver.cpp", "changed")
        r = check(self.root, ["r.json"])
        self.assertFalse(r["integrity_passed"])
        self.assertIn("SHA256 mismatch: solver.cpp", r["errors"][0])

    def test_paths_outside_root_are_rejected(self):
        outside = self.root.parent / "outside.txt"
        outside.write_text("do not read", encoding="utf-8")
        self.write("r.json", {"sources": {"../outside.txt": "0" * 64}})
        r = check(self.root, ["r.json"])
        self.assertFalse(r["integrity_passed"])
        self.assertEqual(r["checked_files"], 1)
        self.assertIn("escapes root", r["errors"][0])

    def test_legacy_artifact_and_opaque_json_data(self):
        data = self.write("raw.json", {"sources": "not a manifest schema"})
        self.write("r.json", {"artifact": {"artifact": "raw.json", "artifact_sha256": data}})
        r = check(self.root, ["r.json"])
        self.assertTrue(r["integrity_passed"])
        self.assertEqual(r["inspected_manifests"], 1)

    def test_missing_invalid_hash_and_wrong_schema(self):
        for content in [{"sources": {"missing": "0" * 64}},
                        {"sources": {"a": "not a hash"}}, {"sources": []}, []]:
            with self.subTest(content=content):
                self.write("r.json", content)
                self.assertFalse(check(self.root, ["r.json"])["integrity_passed"])

    def test_limits_are_enforced(self):
        h = self.write("a", "a")
        self.write("r.json", {"sources": {"a": h}})
        self.assertFalse(check(self.root, ["r.json"], max_files=1)["integrity_passed"])
        self.assertFalse(check(self.root, ["r.json"], max_json_bytes=1)["integrity_passed"])

    def test_report_is_never_overwritten_and_can_not_escape_root(self):
        self.write("r.json", {"passed": False})
        argv = ["--root", str(self.root), "--manifest", "r.json", "--output", "out.json"]
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            self.assertEqual(main(argv), 0)
            old = (self.root / "out.json").read_bytes()
            self.assertEqual(main(argv), 2)
            self.assertEqual((self.root / "out.json").read_bytes(), old)
            self.assertEqual(main(argv[:-1] + ["../escape.json"]), 2)
        self.assertFalse((self.root.parent / "escape.json").exists())

    def test_duplicate_roots_are_checked_once(self):
        self.write("r.json", {"passed": True})
        r = check(self.root, ["r.json", "r.json"])
        self.assertTrue(r["integrity_passed"])
        self.assertEqual(r["inspected_manifests"], 1)


if __name__ == "__main__":
    unittest.main()
