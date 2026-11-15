"""Cheap exporter guards; no numerical experiment or environment installation runs."""

from __future__ import annotations

import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import export_precision_witness as exporter


class ExportGuardTest(unittest.TestCase):
    def test_missing_output_rejected_before_running(self) -> None:
        with (
            patch.object(sys, "argv", ["exporter"]),
            patch.object(exporter.subprocess, "run") as run,
            contextlib.redirect_stderr(io.StringIO()),
        ):
            with self.assertRaises(SystemExit) as failure:
                exporter.main()
            self.assertEqual(failure.exception.code, 2)
            run.assert_not_called()

    def test_existing_output_rejected_before_running(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "existing.json"
            output.write_bytes(b"original")
            with (
                patch.object(sys, "argv", ["exporter", "--output", str(output)]),
                patch.object(exporter.subprocess, "run") as run,
            ):
                with self.assertRaises(FileExistsError):
                    exporter.main()
                run.assert_not_called()
            self.assertEqual(output.read_bytes(), b"original")

    def test_output_created_during_run_is_not_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "raced.json"

            def race(*_args: object, **kwargs: object) -> None:
                work = kwargs["cwd"]
                assert isinstance(work, Path)
                (work / "witness.json").write_bytes(b"new fixture")
                output.write_bytes(b"other evidence")

            with (
                patch.object(sys, "argv", ["exporter", "--output", str(output)]),
                patch.object(exporter.shutil, "copy2"),
                patch.object(exporter.shutil, "copytree"),
                patch.object(exporter.subprocess, "run", side_effect=race),
                self.assertRaises(FileExistsError),
            ):
                exporter.main()
            self.assertEqual(output.read_bytes(), b"other evidence")


if __name__ == "__main__":
    unittest.main()
