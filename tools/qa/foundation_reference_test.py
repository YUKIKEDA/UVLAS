"""Importer/fixture guards use the committed small originals, never rerun the R&D model."""

from __future__ import annotations

import copy
import json
import shutil
import tempfile
import unittest
from collections.abc import Callable
from pathlib import Path
from typing import cast
from unittest.mock import patch

from tools.qa import foundation_reference as reference

FIXTURE = reference.ROOT / reference.FIXTURE


class ReferenceTest(unittest.TestCase):
    def test_verify_reads_frozen_reference_without_regenerating(self) -> None:
        with patch.object(
            reference, "solve_decimal", side_effect=AssertionError("No regeneration")
        ):
            reference.verify_fixture(FIXTURE)

    def test_regeneration_is_byte_stable_and_sources_stay_read_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cache = root / "cache"
            source = cache / reference.LOGICAL
            source.mkdir(parents=True)
            originals = {
                name: (FIXTURE / "source" / name).read_bytes() for name in reference.ARTIFACTS
            }
            for name, data in originals.items():
                (source / name).write_bytes(data)
            output = root / "new"
            reference.import_fixture(output, cache)
            self.assertEqual(
                (output / "fixture.json").read_bytes(), (FIXTURE / "fixture.json").read_bytes()
            )
            self.assertNotIn(b"\r", (output / "fixture.json").read_bytes())
            for name, data in originals.items():
                self.assertEqual((source / name).read_bytes(), data)
                self.assertEqual((output / "source" / name).read_bytes(), data)
            with (
                self.assertRaises(FileExistsError),
                patch.object(reference, "materialize") as generate,
            ):
                reference.import_fixture(output, cache)
            generate.assert_not_called()
            reference.verify_fixture(output)

    def test_duplicate_json_keys_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            reference.json_object(b'{"schema":1,"schema":2}')

    def check_mutation(self, change: Callable[[dict[str, object]], None]) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "copy"
            shutil.copytree(FIXTURE, output)
            fixture = reference.json_object((output / "fixture.json").read_bytes())
            change(fixture)
            # Intentionally permit a malformed JSON constant to reach the reader.
            (output / "fixture.json").write_text(json.dumps(fixture), encoding="utf-8")
            with self.assertRaises(ValueError):
                reference.verify_fixture(output)

    def test_invalid_metadata_and_all_spectral_components_rejected(self) -> None:
        self.check_mutation(lambda item: item.update(schema="unknown"))
        self.check_mutation(lambda item: item.update(extra=True))
        self.check_mutation(lambda item: item.update(relative_solution_tolerance=1e-6))
        self.check_mutation(lambda item: item.update(matrix=[]))
        self.check_mutation(lambda item: item.update(reference_solution=["NaN"] * 24))
        self.check_mutation(lambda item: item.update(reference_solution=["1"] * 24))
        self.check_mutation(lambda item: item.update(source_artifacts=[]))
        for row in range(reference.DIMENSION):
            for component in (0, 1):
                for invalid in (False, True, float("inf"), "0", None):

                    def change(
                        item: dict[str, object],
                        pair_index: int = row,
                        index: int = component,
                        value: object = invalid,
                    ) -> None:
                        pairs = cast("list[list[object]]", copy.deepcopy(item["eigenvalues"]))
                        pairs[pair_index][index] = value
                        item["eigenvalues"] = pairs

                    with self.subTest(row=row, component=component, invalid=invalid):
                        self.check_mutation(change)
        self.check_mutation(lambda item: item.update(eigenvalues=[[0]] * 24))
        self.check_mutation(lambda item: item.update(eigenvalues=[]))

    def test_artifact_mutation_is_not_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "copy"
            shutil.copytree(FIXTURE, output)
            (output / "source" / "pg-free-step.txt").write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "artifact"):
                reference.verify_fixture(output)

    def test_source_and_manifest_cannot_be_restamped_together(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "fixture"
            shutil.copytree(FIXTURE, output)
            for name in (*reference.SOURCES, reference.AUDIT.as_posix()):
                destination = root / name
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes((reference.ROOT / name).read_bytes())
            source = root / reference.SOURCES[0]
            source.write_bytes(source.read_bytes() + b"\n# changed copied source\n")
            data = reference.json_object((output / "fixture.json").read_bytes())
            reference.mapping(data["sources"])[reference.SOURCES[0]] = reference.digest(
                source.read_bytes()
            )
            (output / "fixture.json").write_bytes(reference.encode(data))
            with (
                patch.object(reference, "ROOT", root),
                self.assertRaisesRegex(ValueError, "source provenance"),
            ):
                reference.verify_fixture(output)


if __name__ == "__main__":
    unittest.main()
