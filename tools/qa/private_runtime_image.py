#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Stage or rehearse a private native UVLAS Release image with trusted local inputs."""

from __future__ import annotations

import argparse
import ctypes
import errno
import json
import os
import re
import shutil
import subprocess
import sys
import sysconfig
import tempfile
from pathlib import Path

from runtime_artifact_audit import (
    REQUIRED_NOTICES,
    digest,
    document,
    isolated_environment,
    notice_hashes,
    output,
    require,
    run,
)

COMPONENT = "uvlas-private-runtime"
SCHEMA = "uvlas.private-runtime-image/1"
INVENTORY = "inventory.json"
BUILD_INFO = "build-info.json"
DESTINATION_EXISTS = "Image destination already exists"
NOTICES = Path("share/uvlas/third-party-notices")
LICENSES = (Path("share/uvlas/LICENSE"), Path("share/uvlas/THIRD_PARTY_NOTICES.md"))

SITE_PROBE = """
import json
from pathlib import Path
import sys
import sysconfig

site = Path(sysconfig.get_path("platlib")).resolve()
prefix = Path(sys.prefix).resolve()
if not sys.flags.isolated or sys.prefix == sys.base_prefix or not site.is_relative_to(prefix):
    raise RuntimeError("Expected an isolated venv site-packages directory")
print(json.dumps({"site": str(site), "extension_suffix": sysconfig.get_config_var("EXT_SUFFIX")}))
"""

IMPORT_PROBE = """
import _uvlas
import json
import numpy
from pathlib import Path
import sys

expected = Path(sys.argv[1])
if not sys.flags.isolated or Path(_uvlas.__file__).resolve() != expected:
    raise RuntimeError("Could not import the installed UVLAS extension")
if numpy.__version__ != "2.3.5":
    raise RuntimeError("Expected the audit's fixed NumPy version")
print(json.dumps({"cpython": sys.version.split()[0], "numpy": numpy.__version__}))
"""


def write_json(path: Path, value: dict[str, object]) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, sort_keys=True, indent=2)
        stream.write("\n")


def cli_path(image: Path) -> Path:
    return image / "bin" / ("uvlas.exe" if os.name == "nt" else "uvlas")


def extension_path(image: Path) -> Path:
    suffix = sysconfig.get_config_var("EXT_SUFFIX")
    require(isinstance(suffix, str), "Cannot identify the native CPython module suffix")
    return image / "python" / f"_uvlas{suffix}"


def inventory(image: Path) -> dict[str, str]:
    files: dict[str, str] = {}
    for item in sorted(image.rglob("*")):
        require(not item.is_symlink() and not item.is_junction(), "Image contains a linked entry")
        if item.is_dir():
            continue
        require(item.is_file(), "Image contains a non-regular file")
        relative = item.relative_to(image).as_posix()
        if relative != INVENTORY:
            files[relative] = digest(item)
    return files


def validate_layout(image: Path, files: dict[str, str]) -> None:
    required = {file.as_posix() for file in (*LICENSES, Path(BUILD_INFO))}
    required.update(
        path.relative_to(image).as_posix() for path in (cli_path(image), extension_path(image))
    )
    require(
        not required.difference(files),
        "Image is missing its native CLI, CPython extension or notices",
    )
    notices = {key for key in files if Path(key).parent == NOTICES and key.endswith(".txt")}
    for notice in REQUIRED_NOTICES:
        require((NOTICES / f"{notice}.txt").as_posix() in notices, f"Missing notice: {notice}")
    require(files.keys() == required | notices, "Unexpected files in private runtime image")


def observed_cli(image: Path, commit: str, work: Path) -> dict[str, object]:
    executable = str(cli_path(image))
    env = isolated_environment(work)
    banner = document(output([executable, "--build-info"], work, env))
    require(
        banner.get("schema") == "uvlas.build-info/1"
        and banner.get("commit") == commit
        and banner.get("configuration") == "Release"
        and banner.get("uncommitted_changes") is False,
        "Expected the requested clean Release commit",
    )
    platform_key = banner.get("platform_key")
    native = "windows-x64-msvc" if os.name == "nt" else "linux-x64-gnu14."
    require(
        isinstance(platform_key, str) and platform_key.startswith(native),
        "Unexpected native toolchain",
    )
    require(
        banner.get("dependency_repairs") == {"nlohmann_json_sax_allocation": 1},
        "Missing JSON repair",
    )
    result = run([executable, "--private-image-invalid"], work, env)
    require(
        result.returncode == 2 and "Unknown option:" in result.stderr,
        "CLI invalid option was accepted",
    )
    return banner


def validate_image(image: Path) -> dict[str, str]:
    files = inventory(image)
    validate_layout(image, files)
    manifest = document((image / INVENTORY).read_text(encoding="utf-8"))
    require(manifest == {"schema": SCHEMA, "files": files}, "Image file names or hashes differ")
    return files


def venv_python(venv: Path) -> Path:
    return venv / ("Scripts/python.exe" if os.name == "nt" else "bin/python")


def rehearse(
    image: Path, work: Path, numpy_wheel: Path, ownership_test: Path, commit: str
) -> dict[str, object]:
    # Hash before copying, then verify the copied image before loading a binary.
    files = validate_image(image)
    runtime = work / "relocated"
    shutil.copytree(image, runtime)
    require(files == validate_image(runtime), "Relocated image differs")
    banner = observed_cli(runtime, commit, work)
    require(
        banner == document((runtime / BUILD_INFO).read_text(encoding="utf-8")),
        "CLI and stored build info differ",
    )
    env = isolated_environment(work)
    venv = work / "venv"
    output([sys.executable, "-I", "-B", "-m", "venv", str(venv)], work, env)
    python = str(venv_python(venv))
    output(
        [
            python,
            "-I",
            "-B",
            "-m",
            "pip",
            "--disable-pip-version-check",
            "install",
            "--no-index",
            "--no-deps",
            "--only-binary",
            ":all:",
            str(numpy_wheel),
        ],
        work,
        env,
    )
    observed = document(output([python, "-I", "-B", "-c", SITE_PROBE], work, env))
    site = observed.get("site")
    if not isinstance(site, str):
        raise RuntimeError("Python did not report its venv site-packages")
    site_path = Path(site).resolve()
    require(
        site_path.is_relative_to(venv.resolve()),
        "Will only install the extension into the new venv",
    )
    extension = extension_path(runtime)
    suffix = observed.get("extension_suffix")
    require(
        isinstance(suffix, str) and extension.name == "_uvlas" + suffix,
        "Image extension ABI differs",
    )
    installed = site_path / extension.name
    with installed.open("xb") as destination, extension.open("rb") as source:
        shutil.copyfileobj(source, destination)
    require(digest(installed) == digest(extension), "Installed extension differs")
    runtime_info = document(
        output([python, "-I", "-B", "-c", IMPORT_PROBE, str(installed)], work, env)
    )
    script = work / "ownership.py"
    shutil.copyfile(ownership_test, script)
    # The suite takes a search path. Use an empty directory so import still uses the venv.
    empty = work / "empty"
    empty.mkdir()
    result = run([python, "-I", "-B", str(script), str(empty)], work, env)
    cases = re.search(r"\bRan ([1-9]\d*) tests? in ", result.stderr)
    if result.returncode != 0 or cases is None:
        raise RuntimeError(f"Installed ownership suite failed:\n{result.stderr}")
    return {
        "build_info": banner,
        "files": len(files),
        "inventory_sha256": digest(runtime / INVENTORY),
        "normal_site_packages_import": True,
        "extension_sha256": digest(installed),
        "cli_invalid_option_exit": 2,
        "ownership_cases": int(cases.group(1)),
        **runtime_info,
    }


def publish_directory(source: Path, destination: Path) -> None:
    if os.name == "nt":
        # On Windows rename never replaces an existing destination.
        try:
            source.rename(destination)
        except FileExistsError as error:
            raise RuntimeError(DESTINATION_EXISTS) from error
        return
    # On Linux an ordinary directory rename can replace an existing empty directory.
    # RENAME_NOREPLACE makes the commit safe even if another process wins the race.
    libc = ctypes.CDLL(None, use_errno=True)
    try:
        rename_noreplace = libc.renameat2
    except AttributeError as error:
        raise RuntimeError("Atomic no-replace directory rename is unavailable") from error
    rename_noreplace.argtypes = [
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    ]
    rename_noreplace.restype = ctypes.c_int
    at_fdcwd, rename_noreplace_flag = -100, 1
    status = rename_noreplace(
        at_fdcwd, os.fsencode(source), at_fdcwd, os.fsencode(destination), rename_noreplace_flag
    )
    if status != 0:
        error = ctypes.get_errno()
        if error == errno.EEXIST:
            raise RuntimeError(DESTINATION_EXISTS)
        raise OSError(error, os.strerror(error), destination)


def stage(
    build: Path,
    destination: Path,
    vcpkg_share: Path,
    cmake: str,
    *,
    numpy_wheel: Path,
    ownership_test: Path,
    commit: str,
) -> dict[str, object]:
    require(not os.path.lexists(destination), DESTINATION_EXISTS)
    require(not destination.is_relative_to(build), "Image destination must be outside the build")
    install_env = dict(os.environ)
    install_env.pop("DESTDIR", None)
    install_env.pop("CMAKE_INSTALL_MODE", None)
    with tempfile.TemporaryDirectory(prefix=".uvlas-image-", dir=destination.parent) as temporary:
        image = Path(temporary).resolve()
        # Remove every rehearsal output before committing the sole staged tree.
        with tempfile.TemporaryDirectory(prefix="uvlas-image-rehearsal-") as rehearsal:
            work = Path(rehearsal).resolve()
            output(
                [
                    cmake,
                    "--install",
                    str(build),
                    "--config",
                    "Release",
                    "--component",
                    COMPONENT,
                    "--prefix",
                    str(image),
                ],
                work,
                install_env,
            )
            require(
                cli_path(image).is_file(),
                "Build has no private Release component; enable it in CMake",
            )
            require(extension_path(image).is_file(), "Build has no native CPython 3.12 extension")
            notice_hashes(image / NOTICES, vcpkg_share)
            write_json(image / BUILD_INFO, observed_cli(image, commit, work))
            files = inventory(image)
            validate_layout(image, files)
            write_json(image / INVENTORY, {"schema": SCHEMA, "files": files})
            result = rehearse(image, work, numpy_wheel, ownership_test, commit)
        publish_directory(image, destination)
        return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="operation", required=True)
    create = commands.add_parser(
        "stage", help="install and rehearse a private image before publishing it locally"
    )
    create.add_argument("--build", required=True, type=Path)
    create.add_argument("--vcpkg-share", required=True, type=Path)
    create.add_argument("--cmake", default="cmake")
    verify = commands.add_parser(
        "verify", help="rehearse a previously staged and transported image"
    )
    for mode in (create, verify):
        mode.add_argument("--image", required=True, type=Path)
        mode.add_argument("--commit", required=True)
        mode.add_argument("--numpy-wheel", required=True, type=Path)
        mode.add_argument("--ownership-test", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        require(os.name == "nt" or sys.platform == "linux", "Expected Windows or Linux")
        require(
            sys.implementation.name == "cpython" and sys.version_info[:2] == (3, 12),
            "Expected CPython 3.12",
        )
        require(bool(re.fullmatch(r"[0-9a-f]{40}", args.commit)), "Use a full source commit")
        image = args.image.parent.resolve(strict=True) / args.image.name
        numpy_wheel = args.numpy_wheel.resolve(strict=True)
        ownership_test = args.ownership_test.resolve(strict=True)
        if args.operation == "stage":
            result = stage(
                args.build.resolve(strict=True),
                image,
                args.vcpkg_share.resolve(strict=True),
                args.cmake,
                numpy_wheel=numpy_wheel,
                ownership_test=ownership_test,
                commit=args.commit,
            )
        else:
            require(
                image.is_dir() and not image.is_symlink() and not image.is_junction(),
                "Image is not a directory",
            )
            with tempfile.TemporaryDirectory(prefix="uvlas-image-rehearsal-") as temporary:
                result = rehearse(
                    image, Path(temporary).resolve(), numpy_wheel, ownership_test, args.commit
                )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, RuntimeError, ValueError, subprocess.SubprocessError) as error:
        raise SystemExit(str(error)) from error


if __name__ == "__main__":
    main()
