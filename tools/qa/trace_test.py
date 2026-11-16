"""Trace adopted IDs through formatted Catch2 declarations and fail on missing IDs."""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tools.qa import trace


class TraceTest(unittest.TestCase):
    def setUp(self) -> None:
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        for name, value in (
            ("REPO", self.root),
            ("MATRIX", self.root / "docs/quality/traceability.md"),
        ):
            patcher = patch.object(trace, name, value)
            patcher.start()
            self.addCleanup(patcher.stop)

    def write(self, name: str, text: str) -> None:
        target = self.root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")

    def test_eigen_definitions_and_boundaries_are_part_of_the_catalog(self) -> None:
        self.write("docs/specs/eigen.md", "## EIG-01：identity\n## EIG-14：acceptance\n")
        definitions, _, _ = trace.read_catalog()
        self.assertEqual(set(definitions), {"EIG-01", "EIG-14"})
        self.assertEqual(trace.kind_of("EIG-01"), "decision")
        self.assertEqual(trace.ANY.findall("xEIG-01 EIG-01x EIG-001 EIG-01"), ["EIG-01"])

    def test_formatted_macros_keep_tags_and_macro_line_numbers(self) -> None:
        self.write(
            "tests/sample.cpp",
            'TEST_CASE("one", "[EIG-01][N4-07]") {}\n'
            'TEST_CASE(\n    "wrapped (name)",\n    "[EIG-02]" "[negative]") {}\n'
            'SCENARIO("scenario",\n         "[EIG-03]") {}\n'
            '// TEST_CASE("comment", "[EIG-99]")\n'
            'const char* text = "[EIG-98]";\n',
        )
        tests, _, _ = trace.read_evidence()
        self.assertEqual(tests["EIG-01"], ["tests/sample.cpp:1"])
        self.assertEqual(tests["N4-07"], ["tests/sample.cpp:1"])
        self.assertEqual(tests["EIG-02"], ["tests/sample.cpp:2"])
        self.assertEqual(tests["EIG-03"], ["tests/sample.cpp:5"])
        self.assertNotIn("EIG-99", tests)
        self.assertNotIn("EIG-98", tests)

    def test_comments_raw_strings_and_concatenated_tags(self) -> None:
        text = (
            '/*\nTEST_CASE("comment", "[EIG-99]")\n*/\n'
            'const char* text = R"fixture(\nTEST_CASE("raw", "[EIG-98]")\n)fixture";\n'
            'TEST_CASE("actual", /* between arguments */ "[EIG-" "05]") {}\n'
        )
        self.assertEqual(trace.catch2_tags(text), [(7, ["EIG-05"])])

    def test_digit_separators_do_not_hide_the_following_declaration(self) -> None:
        text = (
            "const auto count = 1'000;\n"
            'TEST_CASE("actual", "[EIG-05]") {}\n'
            "const char value = 'x';\n"
        )
        self.assertEqual(trace.catch2_tags(text), [(2, ["EIG-05"])])

    def test_undefined_eigen_tags_and_records_fail_the_gate(self) -> None:
        self.write("docs/specs/eigen.md", "## EIG-01：identity\n")
        self.write("tests/sample.cpp", 'TEST_CASE("ok", "[EIG-01]") {}\n')
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(trace.main(["matrix"]), 0)
            self.assertEqual(trace.main(["check"]), 0)
        self.write("tests/sample.cpp", 'TEST_CASE(\n "missing",\n "[EIG-99]") {}\n')
        self.write("docs/verification/missing.json", '{"ids": ["EIG-98"]}')
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(trace.main(["matrix"]), 0)
            self.assertEqual(trace.main(["check"]), 1)
        self.assertIn("test names an undefined ID EIG-99", output.getvalue())
        self.assertIn("verification record names an undefined ID EIG-98", output.getvalue())


if __name__ == "__main__":
    unittest.main()
