# SPDX-License-Identifier: Apache-2.0
"""Reproduce the pinned R-06 rehearsal with native CPython 3.12."""

import base64
import copy
import csv
import ctypes
import hashlib
import importlib
import importlib.metadata
import importlib.util
import io
import json
import os
import platform
import shutil
import subprocess
import sys
import sysconfig
import zipfile
from pathlib import Path

COMMIT = "63df6cc00d7180725370201dca3454ad35073cfd"
DIST = "uvlas-private-rehearsal"
VERSION = f"0.1.0.dev0+g{COMMIT[:12]}"
INFO = f"uvlas_private_rehearsal-{VERSION}.dist-info"
WINDOWS = os.name == "nt"
TAG = "cp312-cp312-win_amd64" if WINDOWS else "cp312-cp312-linux_x86_64"
WHEEL = f"{INFO.removesuffix('.dist-info')}-{TAG}.whl"
EXTENSION = "_uvlas" + sysconfig.get_config_var("EXT_SUFFIX")
CLI = "uvlas.exe" if WINDOWS else "uvlas"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def digest(path):
    return sha(path.read_bytes())


def text_json(data):
    return json.dumps(data, ensure_ascii=False, indent=2) + "\n"


def producer(manifest):
    stamp = json.loads(manifest.read_bytes())
    info = stamp["build_info"]
    assert stamp["schema"] == "uvlas.runtime-artifact-producer/1"
    assert (info["schema"], info["commit"], info["configuration"]) == (
        "uvlas.build-info/1",
        COMMIT,
        "Release",
    )
    assert info["uncommitted_changes"] is False
    assert stamp["artifacts"]["cli"]["name"] == CLI
    assert stamp["artifacts"]["extension"]["name"] == EXTENSION
    return stamp


def build(stage, manifest, numpy_wheel, destination):
    stamp = producer(manifest)
    cli = stage / "bin" / CLI
    extension = stage / "python" / EXTENSION
    assert digest(cli) == stamp["artifacts"]["cli"]["sha256"]
    assert digest(extension) == stamp["artifacts"]["extension"]["sha256"]
    destination.mkdir()
    shutil.copytree(stage / "bin", destination / "cli" / "bin")
    notices = stage / "share" / "uvlas"
    shutil.copytree(notices, destination / "cli" / "share" / "uvlas")
    assert (notices / "LICENSE").is_file()
    assert (notices / "THIRD_PARTY_NOTICES.md").is_file()
    wheelhouse = destination / "wheelhouse"
    wheelhouse.mkdir()
    outputs = {
        EXTENSION: extension.read_bytes(),
        f"{INFO}/METADATA": (
            f"Metadata-Version: 2.4\nName: {DIST}\nVersion: {VERSION}\n"
            "Summary: Private offline runtime rehearsal\nRequires-Python: ==3.12.*\n"
            "Requires-Dist: numpy == 2.3.5\n"
        ).encode(),
        f"{INFO}/WHEEL": (
            "Wheel-Version: 1.0\nGenerator: uvlas-r06-rehearsal\n"
            f"Root-Is-Purelib: false\nTag: {TAG}\n"
        ).encode(),
        f"{INFO}/uvlas-runtime-artifacts.json": manifest.read_bytes(),
    }
    for path in sorted(notices.rglob("*")):
        if path.is_file():
            relative = path.relative_to(notices).as_posix()
            outputs[f"{INFO}/licenses/{relative}"] = path.read_bytes()
    records = io.StringIO(newline="")
    rows = csv.writer(records, lineterminator="\n")
    for name, data in outputs.items():
        checksum = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).decode().rstrip("=")
        rows.writerow((name, "sha256=" + checksum, str(len(data))))
    name = f"{INFO}/RECORD"
    rows.writerow((name, "", ""))
    outputs[name] = records.getvalue().encode()
    with zipfile.ZipFile(wheelhouse / WHEEL, "x", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, data in outputs.items():
            member = zipfile.ZipInfo(name)
            member.compress_type = zipfile.ZIP_DEFLATED
            member.create_system = 3
            member.external_attr = 0o100644 << 16
            archive.writestr(member, data)
    shutil.copy2(numpy_wheel, wheelhouse)
    shutil.copy2(manifest, destination / "runtime-artifacts.json")


def run(args, cwd, env=None, expected=0):
    result = subprocess.run(
        args,
        cwd=cwd,
        env=env,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        timeout=180,
        check=False,
    )
    assert result.returncode == expected, (args, result.returncode, result.stdout, result.stderr)
    return result


def runtime_env(root):
    if WINDOWS:
        system_root = Path(os.environ["SYSTEMROOT"])
        return {
            "SYSTEMROOT": str(system_root),
            "WINDIR": str(system_root),
            "PATH": str(system_root / "System32"),
            "TEMP": str(root),
            "TMP": str(root),
        }
    return {"PATH": "/usr/bin:/bin", "LANG": "C.UTF-8", "HOME": str(root), "TMPDIR": str(root)}


def probe(transport, fixtures, workspace):
    workspace.mkdir()
    env = runtime_env(workspace)
    tools = dict.fromkeys(("cl", "gcc", "g++", "cmake", "ninja", "git", "vcpkg", "dumpbin"))
    tools = {name: shutil.which(name, path=env["PATH"]) for name in tools}
    assert all(value is None for value in tools.values()), tools
    run([sys.executable, "-I", "-B", "-m", "venv", str(workspace / "venv")], workspace, env)
    python = workspace / ("venv/Scripts/python.exe" if WINDOWS else "venv/bin/python")
    run(
        [
            str(python),
            "-I",
            "-B",
            "-c",
            "import importlib.util; assert importlib.util.find_spec('_uvlas') is None",
        ],
        workspace,
        env,
    )
    wheelhouse = transport / "wheelhouse"
    wheel = wheelhouse / WHEEL
    pip = [str(python), "-I", "-B", "-m", "pip", "--isolated", "--disable-pip-version-check"]
    install = run(
        [
            *pip,
            "install",
            "--no-index",
            "--no-cache-dir",
            "--only-binary=:all:",
            "--find-links",
            str(wheelhouse),
            str(wheel),
        ],
        workspace,
        env,
    )
    check = run([*pip, "check"], workspace, env)
    result = run(
        [
            str(python),
            "-I",
            "-B",
            str(Path(__file__).resolve()),
            "installed",
            str(transport),
            str(fixtures),
            str(workspace),
        ],
        workspace,
        env,
    )
    evidence = json.loads(result.stdout)
    evidence.update(
        {
            "wheel_name": wheel.name,
            "wheel_sha256": digest(wheel),
            "numpy_wheel_sha256": digest(next(wheelhouse.glob("numpy-*.whl"))),
            "developer_tools": tools,
            "pip_install": install.stdout.splitlines()[-1],
            "pip_check": check.stdout.strip(),
        }
    )
    return evidence


def installed(transport, fixtures, workspace):
    numpy = importlib.import_module("numpy")
    uvlas = importlib.import_module("_uvlas")

    manifest = transport / "runtime-artifacts.json"
    stamp = producer(manifest)
    cli = transport / "cli" / "bin" / CLI
    extension = Path(uvlas.__file__).resolve()
    assert extension.parent == Path(sysconfig.get_path("platlib")).resolve()
    assert digest(cli) == stamp["artifacts"]["cli"]["sha256"]
    assert digest(extension) == stamp["artifacts"]["extension"]["sha256"]
    assert json.loads(run([str(cli), "--build-info"], workspace).stdout) == stamp["build_info"]
    dist = importlib.metadata.distribution(DIST)
    assert dist.version == VERSION and numpy.__version__ == "2.3.5"
    assert f"Tag: {TAG}" in dist.read_text("WHEEL")
    assert dist.files is not None
    licenses = [file for file in dist.files if str(file).startswith(f"{INFO}/licenses/")]
    notice_root = transport / "cli" / "share" / "uvlas"
    assert len(licenses) == sum(path.is_file() for path in notice_root.rglob("*"))
    for file in licenses:
        relative = str(file).removeprefix(f"{INFO}/licenses/")
        assert digest(dist.locate_file(file)) == digest(notice_root / relative)
    assert dist.read_text("uvlas-runtime-artifacts.json") == manifest.read_text(encoding="utf-8")

    examples = json.loads((fixtures / "case-schema-examples.json").read_bytes())
    case_file = workspace / "\u691c\u8a3c case.json"
    internal = text_json(examples["schema_examples"]["structural_static"]).encode()
    case_file.write_bytes(internal)
    passed = run([str(cli), "validate", case_file.name], workspace)
    assert "Case validation succeeded." in passed.stdout
    loaded = uvlas.load_case(case_file)
    assert loaded.ok
    owned = loaded.value
    del loaded
    if not WINDOWS:
        case_file.unlink()
    assert owned.source_sha256 == hashlib.sha256(internal).digest()
    assert owned.source_byte_count == len(internal)

    case_file = workspace / "external case.json"
    case = copy.deepcopy(examples["schema_examples"]["coupled_transient"])
    geometry = copy.deepcopy(examples["fragments"]["Geometry.nastran"])
    bdf = b"BEGIN BULK\nGRID,1,,3.0,4.0,5.0\nENDDATA\n"
    sidecar = b"{}"
    deck = workspace / "skin.bdf"
    deck.write_bytes(bdf)
    (workspace / "sidecar.json").write_bytes(sidecar)
    geometry["deck"].update(source_id="300", path="skin.bdf", sha256=sha(bdf))
    geometry["sidecar"].update(source_id="301", path="sidecar.json", sha256=sha(sidecar))
    assert len(case["model"]["geometry"]) == 1
    case["model"]["geometry"].append(geometry)
    case_file.write_text(text_json(case), encoding="utf-8")
    unsupported = run([str(cli), "validate", case_file.name], workspace, expected=2)
    assert (
        "external reference dependency validation is unsupported at /model/geometry/1/sidecar"
        in unsupported.stderr
    )
    no_case = uvlas.load_case(case_file)
    assert not no_case.ok and no_case.failure.code == uvlas.FailureCode.Unsupported
    deck.write_bytes(bdf + b"$ changed\n")
    mismatch = run([str(cli), "validate", case_file.name], workspace, expected=2)
    assert "reference SHA-256 mismatch at /model/geometry/1/deck" in mismatch.stderr
    bad_digest = uvlas.load_case(case_file)
    assert not bad_digest.ok and bad_digest.failure.code == uvlas.FailureCode.InvalidInput

    empty = workspace / "empty"
    empty.mkdir()
    ownership = run(
        [sys.executable, "-I", "-B", str(fixtures / "test_ownership.py"), str(empty)], workspace
    )
    assert "OK" in ownership.stderr, ownership.stderr
    assert owned.source_sha256 == hashlib.sha256(internal).digest()
    providers = {}
    if WINDOWS:
        win = ctypes.WinDLL("kernel32", use_last_error=True)
        win.GetModuleHandleW.argtypes = [ctypes.c_wchar_p]
        win.GetModuleHandleW.restype = ctypes.c_void_p
        win.GetModuleFileNameW.argtypes = [ctypes.c_void_p, ctypes.c_wchar_p, ctypes.c_uint]
        for name in (
            "python312.dll",
            "msvcp140.dll",
            "vcruntime140.dll",
            "vcruntime140_1.dll",
            "ucrtbase.dll",
        ):
            handle = win.GetModuleHandleW(name)
            if handle:
                output = ctypes.create_unicode_buffer(32768)
                assert win.GetModuleFileNameW(handle, output, len(output))
                providers[name] = digest(Path(output.value))
    else:
        for line in Path("/proc/self/maps").read_text().splitlines():
            fields = line.split(maxsplit=5)
            if len(fields) != 6:
                continue
            found = Path(fields[5])
            if found.is_absolute() and found.name.startswith(
                ("libstdc++.so.", "libgcc_s.so.", "libc.so.", "libm.so.", "ld-linux-")
            ):
                providers[found.name] = digest(found)
        assert any(name.startswith("libstdc++") for name in providers)

    if WINDOWS:
        path = workspace / "\u691c\u8a3c case.json"
        try:
            path.unlink()
        except PermissionError as error:
            assert error.winerror == 32, error
        else:
            raise AssertionError("Expected owning file handle to block deletion")
    owner = {"sha256": owned.source_sha256.hex(), "bytes": owned.source_byte_count}
    del owned
    if WINDOWS:
        path.unlink()

    python_notice = (
        Path(sys.base_prefix) / "LICENSE.txt"
        if WINDOWS
        else Path("/usr/share/doc/python3.12/copyright")
    )
    assert python_notice.is_file(), python_notice
    return {
        "build_info": stamp["build_info"],
        "manifest_sha256": digest(manifest),
        "artifact_hashes": {"cli": digest(cli), "extension": digest(extension)},
        "python": sys.version,
        "soabi": sysconfig.get_config_var("SOABI"),
        "machine": platform.machine(),
        "numpy": numpy.__version__,
        "wheel_license_count": len(licenses),
        "python_notice_sha256": digest(python_notice),
        "provider_hashes": providers,
        "ordinary_cli": {"exit": passed.returncode, "stdout": passed.stdout.strip()},
        "ordinary_owned_case": owner,
        "ordinary_manifest_remove": "after_owner_release" if WINDOWS else "while_owned",
        "matching_deck": {
            "exit": unsupported.returncode,
            "cli": unsupported.stderr.strip(),
            "python": "Unsupported",
        },
        "changed_deck": {
            "exit": mismatch.returncode,
            "cli": mismatch.stderr.strip(),
            "python": "InvalidInput",
        },
        "ownership": ownership.stderr.splitlines(),
    }


if __name__ == "__main__":
    assert sys.version_info[:2] == (3, 12) and sys.implementation.name == "cpython"
    assert platform.machine().lower() in {"amd64", "x86_64"}
    task, *arguments = sys.argv[1:]
    paths = [Path(path).resolve() for path in arguments]
    if task == "build":
        build(*paths)
    else:
        evaluators = {"probe": probe, "installed": installed}
        print(text_json(evaluators[task](*paths)), end="")
