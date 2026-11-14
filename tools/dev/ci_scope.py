"""Select CI jobs from changed paths and verify the resulting required checks."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
from pathlib import Path
from typing import cast

FLAGS = ("build", "format", "python", "skills", "rd", "trace")
JOBS = {"changes", "windows", "linux", "checks"}
TOOLCHAIN_TESTS = {
    "tests/toolchain/pinned_dependencies_test.cpp": "toolchain.pinned_dependencies",
    "tests/toolchain/linear_algebra_test.cpp": "toolchain.linear_algebra",
    "tests/toolchain/arnoldi_test.cpp": "toolchain.arnoldi",
    "tests/toolchain/precise_arithmetic_test.cpp": "toolchain.precise_arithmetic",
    "tests/toolchain/check_notices.cmake": "toolchain.third_party_notices",
}


def test_scope(paths: list[str], *, full: bool = False) -> dict[str, str]:
    compiled = [path for path in paths if select([path])["build"]]
    if full or not compiled or any(path not in TOOLCHAIN_TESTS for path in compiled):
        return {"test_regex": ".*", "build_targets": "", "cpp_files": ""}
    tests = sorted({TOOLCHAIN_TESTS[path] for path in compiled})
    # The notice check has no executable target. This anchor builds no product kernel.
    targets = [
        "uvlas_test_" + name.replace(".", "_")
        for name in tests
        if name != "toolchain.third_party_notices"
    ]
    if not targets:
        targets = ["uvlas_git_revision"]
    return {
        "test_regex": "^(" + "|".join(re.escape(name) for name in tests) + ")$",
        "build_targets": " ".join(targets),
        "cpp_files": " ".join(path for path in compiled if path.endswith(".cpp")),
    }


def select(paths: list[str], *, full: bool = False) -> dict[str, bool]:
    selected: dict[str, bool] = dict.fromkeys(FLAGS, full)
    for path in paths:
        if path in {".github/workflows/ci.yml", "justfile"} or path in {
            "tools/dev/ci_scope.py",
            "tools/dev/ci_scope_test.py",
        }:
            selected = dict.fromkeys(FLAGS, True)
        elif path.startswith((".agents/", ".claude/")) or path == "skills-lock.json":
            selected["skills"] = True
        elif path.startswith(".cursor/") or path in {"LICENSE", "THIRD_PARTY_NOTICES.md"}:
            continue
        elif path == "tests/toolchain/check_notices.cmake":
            selected["build"] = True
        elif path.startswith(("src/", "include/", "tests/", "qa/", "cmake/")) or path in {
            "CMakeLists.txt",
            "CMakePresets.json",
            "vcpkg.json",
            "vcpkg-configuration.json",
            "tools/dev/msvc.ps1",
            ".clang-format",
            ".clang-tidy",
        }:
            selected["build"] = selected["format"] = selected["trace"] = True
        elif path.startswith(("rd/", "tools/rd/")):
            selected["rd"] = True
            selected["python"] |= path.startswith("tools/")
        elif path.startswith("tools/qa/"):
            # results.py is consumed by CTest; trace.py covers test tags and records.
            selected["python"] = selected["trace"] = True
            selected["build"] |= path == "tools/qa/results.py"
        elif path == "tools/dev/sync_skills.py":
            selected["python"] = selected["skills"] = True
        elif path in {"ruff.toml", "pyrightconfig.json"}:
            selected["python"] = True
        elif path.startswith("docs/"):
            # Adopted documents are inputs to trace and frozen evidence links.
            selected["trace"] = selected["rd"] = True
        elif ("/" not in path and path.endswith(".md")) or path.startswith(
            (".github/ISSUE_TEMPLATE/", ".github/pull_request_template")
        ):
            continue
        else:
            # Unknown dependencies must not silently disable required checks.
            selected = dict.fromkeys(FLAGS, True)
    selected["checks"] = any(selected[key] for key in FLAGS if key != "build")
    return selected


def changed_paths(base: str, head: str, *, cwd: Path | None = None) -> list[str]:
    for ref in (base, head):
        if re.fullmatch(r"[0-9a-fA-F]{40,64}", ref) is None:
            raise ValueError("Expected a commit SHA for each diff endpoint")
    result = subprocess.run(
        ["git", "diff", "--name-only", "--no-renames", "-z", f"{base}...{head}", "--"],
        check=True,
        capture_output=True,
        cwd=cwd,
    )
    # --no-renames preserves both the deleted and added path, including renames
    # from a compiled source to an otherwise documentation-only path.
    return [p.decode("utf-8") for p in result.stdout.split(b"\0") if p]


def object_map(value: object) -> dict[str, object]:
    if not isinstance(value, dict):
        raise ValueError("Expected an object with string keys")
    mapping = cast("dict[object, object]", value)
    if not all(isinstance(k, str) for k in mapping):
        raise ValueError("Expected an object with string keys")
    return cast("dict[str, object]", value)


def verify_needs(value: object) -> list[str]:
    needs = object_map(value)
    if needs.keys() != JOBS:
        return ["Unexpected job set; update the selector and aggregate together"]
    changes = object_map(needs["changes"])
    if changes.get("result") != "success":
        return ["Change selection did not succeed"]
    outputs = object_map(changes.get("outputs"))
    if any(outputs.get(flag) not in ("true", "false") for flag in ("build", "checks")):
        return ["Missing or invalid selection outputs"]
    errors: list[str] = []
    for job, flag in (("windows", "build"), ("linux", "build"), ("checks", "checks")):
        expected = "success" if outputs[flag] == "true" else "skipped"
        actual = object_map(needs[job]).get("result")
        if actual != expected:
            errors.append(f"{job}: expected {expected}, got {actual}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-needs", action="store_true")
    args = parser.parse_args()
    if args.verify_needs:
        value: object = json.loads(os.environ["CI_NEEDS_JSON"])
        errors = verify_needs(value)
        print(json.dumps({"errors": errors}))
        return int(bool(errors))
    base, head = os.environ.get("CI_BASE_SHA", ""), os.environ.get("CI_HEAD_SHA", "")
    full = os.environ.get("CI_FORCE_FULL") == "true" or not base or set(base) == {"0"}
    paths = [] if full else changed_paths(base, head)
    selected = select(paths, full=full)
    scope = test_scope(paths, full=full)
    print(json.dumps({"paths": paths, "selected": selected, "test_scope": scope}, indent=2))
    with Path(os.environ["GITHUB_OUTPUT"]).open("a", encoding="utf-8") as stream:
        stream.writelines(f"{key}={str(value).lower()}\n" for key, value in selected.items())
        stream.writelines(f"{key}={value}\n" for key, value in scope.items())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
