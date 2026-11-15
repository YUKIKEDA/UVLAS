"""Regression checks for change selection and fail-closed CI aggregation."""

from __future__ import annotations

import subprocess
import tempfile
import unittest
from itertools import product
from pathlib import Path

from ci_scope import changed_paths, select, test_scope, verify_needs


class ScopeTest(unittest.TestCase):
    def test_path_groups(self) -> None:
        cases: dict[str, set[str]] = {
            "README.md": set(),
            ".cursor/rules/workflow.mdc": set(),
            ".agents/skills/solver-rd/SKILL.md": {"skills", "checks"},
            ".claude/skills/solver-rd/SKILL.md": {"skills", "checks"},
            "tools/dev/sync_skills.py": {"python", "skills", "checks"},
            "tools/rd/check_tree.py": {"python", "rd", "checks"},
            "tools/qa/trace.py": {"python", "trace", "checks"},
            "src/core/build_info.cpp": {"build", "benchmark", "format", "trace", "checks"},
            "tests/toolchain/arnoldi_test.cpp": {"build", "format", "trace", "checks"},
            "tests/toolchain/check_notices.cmake": {"build"},
            "docs/specs/foundation/contract.md": {"trace", "rd", "checks"},
            "rd/foundation/input.json": {"rd", "checks"},
            "ruff.toml": {"python", "checks"},
            "benchmarks/dimensions_bench.cpp": {"benchmark", "format", "checks"},
            "benchmarks/private.hpp": {"benchmark", "format", "checks"},
            "benchmarks/README.md": set(),
            "benchmarks/CMakeLists.txt": {"benchmark"},
            "tests/support/allocation_audit.hpp": {"build", "format", "trace", "checks"},
        }
        for path, expected in cases.items():
            with self.subTest(path=path):
                self.assertEqual(
                    {key for key, enabled in select([path]).items() if enabled}, expected
                )

    def test_shared_unknown_and_manual_full(self) -> None:
        for path in ("unknown-new-tool.sh", ".github/workflows/ci.yml", "tools/dev/ci_scope.py"):
            with self.subTest(path=path):
                self.assertTrue(all(select([path]).values()))
        self.assertTrue(all(select([], full=True).values()))
        self.assertFalse(any(select([]).values()))
        self.assertFalse(select(["benchmarks/dimensions_bench.cpp"])["build"])
        self.assertTrue(select(["benchmarks/dimensions_bench.cpp"])["benchmark"])
        self.assertFalse(select(["tests/core/instrumentation_test.cpp"])["benchmark"])
        self.assertTrue(select(["benchmarks/dimensions_bench.cpp", "src/core/source.cpp"])["build"])
        self.assertTrue(select(["README.md", "vcpkg.json"])["build"])

    def test_independent_tests_are_selected_without_product_tests(self) -> None:
        scoped = test_scope(["README.md", "tests/toolchain/arnoldi_test.cpp"])
        self.assertEqual(scoped["test_regex"], r"^(toolchain\.arnoldi)$")
        self.assertEqual(scoped["build_targets"], "uvlas_test_toolchain_arnoldi")
        self.assertEqual(scoped["cpp_files"], "tests/toolchain/arnoldi_test.cpp")
        self.assertEqual(test_scope(["src/core/build_info.cpp"])["test_regex"], ".*")
        self.assertEqual(
            test_scope(["tests/toolchain/arnoldi_test.cpp"], full=True)["test_regex"], ".*"
        )
        self.assertEqual(
            test_scope(["tests/toolchain/check_notices.cmake"])["build_targets"],
            "uvlas_git_revision",
        )

    def test_job_result_matrix(self) -> None:
        for build, benchmark in product((False, True), repeat=2):
            for checks in (False, True):
                needs: dict[str, object] = {
                    "changes": {
                        "result": "success",
                        "outputs": {
                            "build": str(build).lower(),
                            "benchmark": str(benchmark).lower(),
                            "checks": str(checks).lower(),
                        },
                    },
                    "windows": {"result": "success" if build else "skipped"},
                    "linux": {"result": "success" if build else "skipped"},
                    "benchmark": {"result": "success" if benchmark else "skipped"},
                    "checks": {"result": "success" if checks else "skipped"},
                }
                self.assertEqual(verify_needs(needs), [])
                for job in ("windows", "linux", "benchmark", "checks"):
                    original = needs[job]
                    for result in ("success", "skipped", "failure", "cancelled"):
                        if original == {"result": result}:
                            continue
                        with self.subTest(build=build, checks=checks, job=job, result=result):
                            needs[job] = {"result": result}
                            self.assertTrue(verify_needs(needs))
                    needs[job] = original
                needs["changes"] = {"result": "failure"}
                self.assertTrue(verify_needs(needs))

    def test_invalid_job_metadata_is_rejected(self) -> None:
        self.assertTrue(verify_needs({}))
        for outputs in ({}, {"build": "unknown", "checks": "false"}):
            needs: dict[str, object] = {
                job: {"result": "skipped"} for job in ("windows", "linux", "benchmark", "checks")
            }
            needs["changes"] = {"result": "success", "outputs": outputs}
            self.assertTrue(verify_needs(needs))
        with self.assertRaises(ValueError):
            verify_needs([])
        with self.assertRaises(ValueError):
            changed_paths("--invalid-ref", "0" * 40)

    def test_git_rename_and_delete_keep_original_scope(self) -> None:
        with tempfile.TemporaryDirectory(prefix="uvlas-ci-scope-") as folder:

            def git(*args: str) -> str:
                return subprocess.run(
                    ["git", "-C", folder, *args],
                    check=True,
                    capture_output=True,
                    text=True,
                ).stdout.strip()

            def commit() -> str:
                git("add", ".")
                git(
                    "-c",
                    "user.name=CI Test",
                    "-c",
                    "user.email=ci@example.invalid",
                    "commit",
                    "-qm",
                    "fixture",
                )
                return git("rev-parse", "HEAD")

            git("init", "-q")
            source = Path(folder) / "src/core/source.cpp"
            source.parent.mkdir(parents=True)
            source.write_text("int f() { return 0; }\n", encoding="utf-8")
            base = commit()
            git("mv", "src/core/source.cpp", "README.md")
            head = commit()
            paths = changed_paths(base, head, cwd=Path(folder))
            self.assertEqual(set(paths), {"src/core/source.cpp", "README.md"})
            self.assertTrue(select(paths)["build"])
            (Path(folder) / "README.md").unlink()
            last = commit()
            self.assertEqual(changed_paths(head, last, cwd=Path(folder)), ["README.md"])


if __name__ == "__main__":
    unittest.main()
