#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Exercise the full private runtime image workflow against a native Release build."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

INVENTORY = "inventory.json"
NOTICE = "share/uvlas/third-party-notices/boost-core.txt"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def files(directory):
    return {
        item.relative_to(directory).as_posix(): sha(item)
        for item in directory.rglob("*")
        if item.is_file()
    }


def record(result, expected, reason=""):
    require(
        result.returncode == expected, f"Unexpected status {result.returncode}\n{result.stderr}"
    )
    if expected != 0:
        require(not result.stdout and reason in result.stderr, f"Wrong failure:\n{result.stderr}")
    return result


def verify(args, base):
    work = base.resolve()
    child = [sys.executable, "-B", str(args.tool.resolve())]
    common = [
        "--commit",
        args.commit,
        "--numpy-wheel",
        str(args.numpy_wheel.resolve()),
        "--ownership-test",
        str(args.ownership_test.resolve()),
    ]
    staging = [
        "--build",
        str(args.build.resolve()),
        "--vcpkg-share",
        str(args.vcpkg_share.resolve()),
        "--cmake",
        args.cmake,
    ]

    def stage(image, *, flags=common, build_flags=staging):
        return [*child, "stage", "--image", str(image), *build_flags, *flags]

    def inspect(image):
        return [*child, "verify", "--image", str(image), *common]

    def run(command):
        return subprocess.run(
            command, cwd=work, capture_output=True, text=True, check=False, timeout=150
        )

    default = work / "default-install"
    record(
        run(
            [
                args.cmake,
                "--install",
                str(args.build.resolve()),
                "--config",
                "Release",
                "--prefix",
                str(default),
            ]
        ),
        0,
    )
    require(
        not default.exists() or not files(default), "Default install included the private image"
    )
    first = work / "first"
    second = work / "second"
    one = json.loads(record(run(stage(first)), 0).stdout)
    two = json.loads(record(run(stage(second)), 0).stdout)
    ready = files(first)
    require(ready == files(second) and one == two, "Repeat images differ")
    binary_source = args.build / "Release/uvlas.exe" if os.name == "nt" else args.build / "uvlas"
    module_source = args.build / ("src/python/Release" if os.name == "nt" else "src/python")
    extension = next((first / "python").iterdir())
    built_extension = module_source / extension.name
    require(sha(binary_source) == sha(first / "bin" / binary_source.name), "Staged CLI differs")
    require(sha(built_extension) == sha(extension), "Staged extension differs")
    require(one["files"] + 1 == len(ready), "Inventory length differs")
    require(
        one["normal_site_packages_import"] and one["ownership_cases"] == 10, "Missing ownership run"
    )
    require(one["cli_invalid_option_exit"] == 2, "Missing CLI rejection")

    transported = work / "transported"
    first.rename(transported)
    consumer = json.loads(record(run(inspect(transported)), 0).stdout)
    require(consumer == one and not first.exists(), "Transport did not resolve the same image")

    already = work / "already"
    already.mkdir()
    marker = already / "user-data"
    marker.write_bytes(b"retained")
    record(run(stage(already)), 1, "Image destination already exists")
    require(
        files(already) == {"user-data": hashlib.sha256(b"retained").hexdigest()},
        "Old output changed",
    )

    late = work / "bad-commit"
    wrong_commit = [*common]
    wrong_commit[1] = "a" * 40
    record(run(stage(late, flags=wrong_commit)), 1, "Expected the requested clean Release commit")
    require(not late.exists(), "Late failure exposed output")

    disabled = work / "disabled-component"
    other_build = [*staging]
    other_build[1] = str(args.no_component_build.resolve())
    record(
        run(stage(disabled, build_flags=other_build)), 1, "Build has no private Release component"
    )
    require(not disabled.exists(), "Disabled component exposed output")

    tampered = work / "tampered"
    shutil.copytree(second, tampered)
    with (tampered / "bin" / binary_source.name).open("ab") as modified:
        modified.write(b"changed")
    record(run(inspect(tampered)), 1, "Image file names or hashes differ")

    missing = work / "missing-transitive-notice"
    shutil.copytree(second, missing)
    (missing / NOTICE).unlink()
    record(run(inspect(missing)), 1, "Image file names or hashes differ")

    wrong_module = work / "incorrect-extension"
    shutil.copytree(second, wrong_module)
    broken = wrong_module / "python" / extension.name
    broken.write_bytes(b"not a native CPython module")
    inventory = wrong_module / INVENTORY
    coherent = json.loads(inventory.read_text(encoding="utf-8"))
    coherent["files"]["python/" + extension.name] = sha(broken)
    inventory.write_text(json.dumps(coherent), encoding="utf-8")
    record(run(inspect(wrong_module)), 1, "ImportError")

    staged_directory = work / "atomic-staged"
    staged_directory.mkdir()
    (staged_directory / "marker").write_bytes(b"staged source")
    empty_winner = work / "atomic-existing-empty"
    empty_winner.mkdir()
    atomic = """
from pathlib import Path
import runpy
import sys

tool, source, winner = map(Path, sys.argv[1:])
sys.path.insert(0, str(tool.parent))
publish = runpy.run_path(str(tool))["publish_directory"]
try:
    publish(source, winner)
except RuntimeError as error:
    if str(error) != "Image destination already exists":
        raise
else:
    raise RuntimeError("Replaced the concurrently created empty directory")
"""
    record(
        run(
            [
                sys.executable,
                "-B",
                "-c",
                atomic,
                str(args.tool.resolve()),
                str(staged_directory),
                str(empty_winner),
            ]
        ),
        0,
    )
    require(not files(empty_winner), "Concurrent empty directory changed")
    require(
        (staged_directory / "marker").read_bytes() == b"staged source", "Prepared image was lost"
    )

    race = work / "race"
    command = stage(race)
    processes = [
        subprocess.Popen(
            command, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        for _ in range(2)
    ]
    for process in processes:
        stdout, stderr = process.communicate(timeout=150)
        if process.returncode == 0:
            require(json.loads(stdout) == one, "Race winner produced different evidence")
        else:
            require(not stdout and "exist" in stderr.lower(), f"Unexpected race failure:\n{stderr}")
    require(
        sorted(process.returncode for process in processes) == [0, 1],
        "Concurrent outputs both succeeded or failed",
    )
    require(files(race) == ready, "Concurrent placement corrupted the result")
    require(not list(work.glob(".uvlas-image-*")), "Temporary placement directory survived")
    require(files(transported) == ready and files(second) == ready, "Verified inputs changed")
    return {
        "verified": one,
        "file_sha256": ready,
        "numpy_wheel_sha256": sha(args.numpy_wheel),
        "sha256_of_ownership_test": sha(args.ownership_test),
        "controls": {
            "default_install_private_files": 0,
            "preexisting_destination": 1,
            "bad_commit": 1,
            "disabled_component": 1,
            "tampered_cli": 1,
            "missing_transitive_notice": 1,
            "incorrect_native_extension": 1,
            "concurrent_exit_codes": [0, 1],
            "native_empty_destination_retained": True,
            "remaining_temporary_directories": 0,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for path in (
        "tool",
        "build",
        "no-component-build",
        "vcpkg-share",
        "numpy-wheel",
        "ownership-test",
    ):
        parser.add_argument("--" + path, type=Path, required=True)
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="uvlas-image-check-") as temporary:
        print(json.dumps(verify(args, Path(temporary)), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
