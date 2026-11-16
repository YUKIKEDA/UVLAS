"""Standard-library guards for the isolated full precision-system exporter."""

from __future__ import annotations

import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import export_precision_system as exporter


class ExportSystemGuardTest(unittest.TestCase):
    def test_missing_output_never_runs_a_child(self) -> None:
        with (
            patch.object(sys, "argv", ["exporter"]),
            patch.object(exporter.subprocess, "run") as run,
            contextlib.redirect_stderr(io.StringIO()),
            self.assertRaises(SystemExit) as failure,
        ):
            exporter.main()
        self.assertEqual(failure.exception.code, 2)
        run.assert_not_called()

    def test_existing_output_is_immutable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "fixture.json"
            output.write_bytes(b"existing evidence")
            with (
                patch.object(exporter.subprocess, "run") as run,
                self.assertRaises(FileExistsError),
            ):
                exporter.run_export(root, output)
            self.assertEqual(output.read_bytes(), b"existing evidence")
            run.assert_not_called()

    def test_altered_original_rejected_before_copy_or_child(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            record = root / exporter.RECORD
            record.parent.mkdir(parents=True)
            record.write_bytes(b"a different original record")
            with (
                patch.object(exporter.shutil, "copytree") as copy,
                patch.object(exporter.subprocess, "run") as run,
                self.assertRaisesRegex(ValueError, "Frozen input hash mismatch"),
            ):
                exporter.run_export(root, root / "new.json")
            copy.assert_not_called()
            run.assert_not_called()
            self.assertEqual(record.read_bytes(), b"a different original record")

    def test_hash_guard_accepts_exact_bytes_and_refuses_changed_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "record"
            path.write_bytes(b"abc")
            expected = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
            exporter.require_hash(path, expected)
            path.write_bytes(b"abd")
            with self.assertRaises(ValueError):
                exporter.require_hash(path, expected)

    def test_raced_destination_is_not_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "raced.json"

            def race(*_args: object, **kwargs: object) -> None:
                work = kwargs["cwd"]
                assert isinstance(work, Path)
                (work / "system.json").write_bytes(b"new generated fixture")
                output.write_bytes(b"other evidence")

            with (
                patch.object(exporter, "require_hash"),
                patch.object(exporter, "resolve_record", return_value=root / "record"),
                patch.object(exporter, "require_record_sources"),
                patch.object(exporter, "source_hashes", return_value={}),
                patch.object(exporter.shutil, "copy2"),
                patch.object(exporter.shutil, "copytree"),
                patch.object(exporter.subprocess, "run", side_effect=race),
                self.assertRaises(FileExistsError),
            ):
                exporter.run_export(root, output)
            self.assertEqual(output.read_bytes(), b"other evidence")

    def test_original_mutation_is_detected_after_child_failure(self) -> None:
        with (
            tempfile.TemporaryDirectory() as temporary,
            patch.object(exporter, "require_hash"),
            patch.object(exporter, "resolve_record", return_value=Path(temporary) / "record"),
            patch.object(exporter, "require_record_sources"),
            patch.object(exporter, "source_hashes", side_effect=[{"x": "old"}, {"x": "new"}]),
            patch.object(exporter.shutil, "copy2"),
            patch.object(exporter.shutil, "copytree"),
            patch.object(exporter.subprocess, "run", side_effect=RuntimeError("child failed")),
            self.assertRaisesRegex(RuntimeError, "Original R&D files changed"),
        ):
            exporter.run_export(Path(temporary), Path(temporary) / "new.json")

    def test_frozen_source_content_cannot_change_with_record_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            record = root / exporter.RECORD
            record.parent.mkdir(parents=True)
            source = root / "rd/check_uvlm_column_update.py"
            source.write_bytes(b"abc")
            digest = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
            record.write_text(json.dumps({"source_sha256": digest, "dependency_sha256": {}}))
            exporter.require_record_sources(root)
            source.write_bytes(b"abd")
            with self.assertRaisesRegex(ValueError, "Frozen input hash mismatch"):
                exporter.require_record_sources(root)

    def test_missing_or_ambiguous_generator_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            record = root / exporter.RECORD
            record.parent.mkdir(parents=True)
            record.write_text(json.dumps({"source_sha256": "irrelevant", "dependency_sha256": {}}))
            with self.assertRaisesRegex(ValueError, "one location"):
                exporter.require_record_sources(root)
            for directory in ("a", "b"):
                source = root / "rd" / directory / "check_uvlm_column_update.py"
                source.parent.mkdir()
                source.write_bytes(b"abc")
            with self.assertRaisesRegex(ValueError, "one location"):
                exporter.require_record_sources(root)

    def test_relocated_record_is_resolved_and_still_hash_checked(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            relocated = root / "rd/another-family/records" / Path(exporter.RECORD).name
            relocated.parent.mkdir(parents=True)
            relocated.write_bytes(b"abc")
            resolved = exporter.resolve_record(root)
            self.assertEqual(resolved, relocated)
            exporter.require_hash(
                resolved, "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
            )
            relocated.write_bytes(b"abd")
            with self.assertRaisesRegex(ValueError, "Frozen input hash mismatch"):
                exporter.require_hash(
                    exporter.resolve_record(root),
                    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
                )

    def test_missing_or_duplicate_record_is_refused_before_subprocess(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (
                patch.object(exporter.subprocess, "run") as run,
                self.assertRaisesRegex(ValueError, "record must have one location"),
            ):
                exporter.run_export(root, root / "new.json")
            run.assert_not_called()
            for family in ("first", "second"):
                record = root / "rd" / family / Path(exporter.RECORD).name
                record.parent.mkdir(parents=True)
                record.write_bytes(b"abc")
            with (
                patch.object(exporter.subprocess, "run") as run,
                self.assertRaisesRegex(ValueError, "record must have one location"),
            ):
                exporter.run_export(root, root / "new.json")
            run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
