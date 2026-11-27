#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Inspect and start locally built UVLAS artifacts from a private staging directory."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path, PurePosixPath
from typing import NamedTuple, cast

REQUIRED_NOTICES = (
    "boost-hash2",
    "boost-multiprecision",
    "eigen3",
    "nlohmann-json",
    "pybind11",
    "spectra",
)

# Run in the fresh interpreter, which supplies the actual extension loader and NumPy ABI.
PYTHON_PROBE = """
import ctypes
import hashlib
import importlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import sys
import sysconfig

target, cli_imports, extension_imports = sys.argv[1:]
module = Path(target).resolve()
sys.path.insert(0, str(module.parent))
loaded = importlib.import_module("_uvlas")
if Path(loaded.__file__).resolve() != module:
    raise RuntimeError("Resolved a different extension")

distribution = metadata.distribution("numpy")
notice_files = {}
bundled_runtimes = {}
wheel_tags = []
for file in distribution.files or ():
    relative = file.as_posix()
    if ".dist-info/" in relative:
        name = relative.split(".dist-info/", 1)[1]
        if name.lower().startswith(("license", "notice", "copying")):
            content = distribution.locate_file(file).read_bytes()
            notice_files[name] = hashlib.sha256(content).hexdigest()
        if name == "WHEEL":
            wheel_tags = [
                line.removeprefix("Tag: ")
                for line in distribution.locate_file(file).read_text(encoding="utf-8").splitlines()
                if line.startswith("Tag: ")
            ]
    elif file.parts[0] == "numpy.libs":
        content = distribution.locate_file(file).read_bytes()
        bundled_runtimes[file.name] = hashlib.sha256(content).hexdigest()
if not notice_files:
    raise RuntimeError("NumPy distribution did not report license files")

result = {
    "implementation": sys.implementation.name,
    "version": sys.version.split()[0],
    "soabi": sysconfig.get_config_var("SOABI"),
    "extension_suffix": sysconfig.get_config_var("EXT_SUFFIX"),
    "isolated": bool(sys.flags.isolated),
    "venv": sys.prefix != sys.base_prefix,
    "loaded_extension": module.name,
    "numpy": distribution.version,
    "numpy_license_expression": distribution.metadata.get("License-Expression"),
    "numpy_license_files": notice_files,
    "numpy_bundled_runtimes": bundled_runtimes,
    "numpy_wheel_tags": wheel_tags,
}
if os.name == "nt":
    python_license = Path(sys.base_prefix) / "LICENSE.txt"
    notice_path = "base_python/LICENSE.txt"
else:
    base = Path(sys.base_prefix).resolve()
    stdlib = Path(sysconfig.get_path("stdlib")).resolve()
    if not stdlib.is_relative_to(base):
        raise RuntimeError("CPython stdlib notice is outside the base interpreter")
    packaged_notice = Path("/usr/share/doc/python3.12/copyright")
    if base == Path("/usr") and packaged_notice.is_file():
        python_license = packaged_notice
        notice_path = str(python_license)
    else:
        python_license = stdlib / "LICENSE.txt"
        notice_path = "stdlib/LICENSE.txt"
try:
    notice_bytes = python_license.read_bytes()
except OSError as error:
    raise RuntimeError(f"Cannot read required CPython notice at {python_license}") from error
result["python_license"] = {
    "file": python_license.name,
    "path": notice_path,
    "sha256": hashlib.sha256(notice_bytes).hexdigest(),
}
if os.name == "nt":
    from ctypes import wintypes

    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    loaded_module = kernel.GetModuleHandleW
    loaded_module.argtypes = [wintypes.LPCWSTR]
    loaded_module.restype = wintypes.HMODULE
    module_name = kernel.GetModuleFileNameW
    module_name.argtypes = [wintypes.HMODULE, wintypes.LPWSTR, wintypes.DWORD]
    module_name.restype = wintypes.DWORD
    roots = {
        "windows": Path(os.environ["SYSTEMROOT"]).resolve(),
        "base_python": Path(sys.base_prefix).resolve(),
        "venv": Path(sys.prefix).resolve(),
        "stage": module.parent,
    }

    def origin(handle):
        path_buffer = ctypes.create_unicode_buffer(32768)
        size = module_name(handle, path_buffer, len(path_buffer))
        if size == 0 or size == len(path_buffer):
            raise ctypes.WinError(ctypes.get_last_error())
        path = Path(path_buffer.value).resolve()
        for name, root in roots.items():
            if path.is_relative_to(root):
                return name + "/" + path.relative_to(root).as_posix()
        raise RuntimeError("Loaded a runtime outside the known roots: " + path.name)

    extension_loaded = {}
    for name in json.loads(extension_imports):
        if name.lower().startswith(("api-ms-", "ext-ms-")):
            continue
        handle = loaded_module(name)
        if not handle:
            raise RuntimeError("Missing extension runtime module: " + name)
        extension_loaded[name] = origin(handle)
    result["extension_loaded_modules"] = extension_loaded

    # Independent loader resolution of the PE imports under the staged environment.
    handles = {}
    for name in set(json.loads(cli_imports)) | set(json.loads(extension_imports)):
        handles[name] = ctypes.WinDLL(name)
    result["windows_loader_probe"] = {
        name: origin(library._handle) for name, library in handles.items()
    }

print(json.dumps(result, sort_keys=True))
"""


PYTHON_EXCEPTION_PROBE = """
import importlib
from pathlib import Path
import sys

module, source = map(Path, sys.argv[1:])
sys.path.insert(0, str(module.parent))
uvlas = importlib.import_module("_uvlas")
if Path(uvlas.__file__).resolve() != module:
    raise RuntimeError("Resolved a different extension")
loaded = uvlas.load_case(source)
if loaded.ok or loaded.failure.code != uvlas.FailureCode.InvalidInput:
    raise RuntimeError("Malformed JSON did not return InvalidInput")
print("Malformed JSON rejected.")
"""

# GLIBC's binding report distinguishes the caller from the selected provider.
GNU_EXCEPTION_BINDING = re.compile(
    r"^\s*\d+:\s*binding file (.+?) \[[0-9]+\] to (.+?) \[[0-9]+\]: "
    r"normal symbol `_ZNSt13runtime_errorC1EPKc' \[GLIBCXX_3\.4\.21\]\s*$"
)


GNU_PARSER_HEADERS = frozenset(("stl_vector.h", "vector.tcc", "stl_bvector.h", "basic_string.tcc"))
GNU_PARSER_OBJECT = Path("CMakeFiles/uvlas_io.dir/src/io/input_document.cpp.o")

PYTHON_EXTERNAL_CASE_PROBE = """
import importlib
from pathlib import Path
import sys

module, source = map(Path, sys.argv[1:3])
sys.path.insert(0, str(module.parent))
uvlas = importlib.import_module("_uvlas")
if Path(uvlas.__file__).resolve() != module:
    raise RuntimeError("Resolved a different extension")
loaded = uvlas.load_case(source)
if not loaded.ok or loaded.value.source_sha256.hex() != sys.argv[3]:
    raise RuntimeError("Python did not accept the pinned external Case")
print("External Case accepted.")
"""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def digest(path: Path) -> str:
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def run(
    command: list[str], cwd: Path, env: dict[str, str] | None = None
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
        timeout=120,
    )


def output(command: list[str], cwd: Path, env: dict[str, str] | None = None) -> str:
    result = run(command, cwd, env)
    require(
        result.returncode == 0, f"{Path(command[0]).name}: {result.returncode}\n{result.stderr}"
    )
    return result.stdout


def object_members(value: object, description: str) -> dict[str, object]:
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object for {description}")
    return cast(dict[str, object], value)


def document(text: str) -> dict[str, object]:
    value: object = json.loads(text)
    return object_members(value, "document")


class ProducedArtifact(NamedTuple):
    name: str
    relative_path: PurePosixPath
    sha256: str


class RuntimeArtifactProducer(NamedTuple):
    build_info: dict[str, object]
    cli: ProducedArtifact
    extension: ProducedArtifact
    sha256: str


def producer_artifact(value: object) -> ProducedArtifact:
    record = object_members(value, "produced artifact")
    name, relative, expected = (record.get(key) for key in ("name", "relative_path", "sha256"))
    if not isinstance(name, str) or not isinstance(relative, str) or not isinstance(expected, str):
        raise ValueError("Incomplete produced artifact")
    path = PurePosixPath(relative)
    require(
        path.as_posix() == relative
        and path.name == name
        and all(part != ".." and re.fullmatch(r"[A-Za-z0-9_.+-]+", part) for part in path.parts),
        "Invalid producer artifact relative path",
    )
    require(bool(re.fullmatch(r"[0-9a-f]{64}", expected)), "Invalid producer artifact SHA-256")
    return ProducedArtifact(name, path, expected)


def read_producer_manifest(path: Path, commit: str) -> RuntimeArtifactProducer:
    require(bool(re.fullmatch(r"[0-9a-f]{40}", commit)), "Use the full source commit")
    data = path.read_bytes()
    manifest = document(data.decode("utf-8"))
    require(
        manifest.get("schema") == "uvlas.runtime-artifact-producer/1",
        "Unexpected runtime artifact producer manifest",
    )
    info = object_members(manifest.get("build_info"), "producer build info")
    require(
        info.get("schema") == "uvlas.build-info/1"
        and info.get("commit") == commit
        and info.get("configuration") == "Release"
        and info.get("uncommitted_changes") is False,
        "Producer manifest does not identify the requested clean Release commit",
    )
    artifacts = object_members(manifest.get("artifacts"), "producer artifacts")
    return RuntimeArtifactProducer(
        info,
        producer_artifact(artifacts.get("cli")),
        producer_artifact(artifacts.get("extension")),
        hashlib.sha256(data).hexdigest(),
    )


def verify_producer_artifacts(
    producer: RuntimeArtifactProducer, cli: Path, extension: Path
) -> None:
    for kind, path, expected in (
        ("CLI", cli, producer.cli),
        ("extension", extension, producer.extension),
    ):
        require(
            path.name == expected.name and digest(path) == expected.sha256,
            f"{kind} does not match the producer manifest: {path.name}",
        )


def isolated_environment(directory: Path) -> dict[str, str]:
    if os.name == "nt":
        root = Path(os.environ["SYSTEMROOT"])
        return {
            "SYSTEMROOT": str(root),
            "WINDIR": str(root),
            "PATH": str(root / "System32"),
            "TEMP": str(directory),
            "TMP": str(directory),
        }
    return {"PATH": "/usr/bin:/bin", "LANG": "C.UTF-8", "HOME": str(directory)}


def pe(target: Path, dumpbin: str, extension: bool) -> tuple[dict[str, object], list[str]]:
    headers = output([dumpbin, "/HEADERS", str(target)], target.parent)
    require("8664 machine (x64)" in headers, "Expected an AMD64 PE")
    require("20B magic # (PE32+)" in headers, "Expected PE32+")
    dependents = output([dumpbin, "/DEPENDENTS", str(target)], target.parent)
    imports = re.findall(r"(?im)^\s{4}([\w.-]+\.dll)\s*$", dependents)
    require(bool(imports), f"No PE imports observed in {target.name}")
    if extension:
        exports = output([dumpbin, "/EXPORTS", str(target)], target.parent)
        require(bool(re.search(r"\bPyInit__uvlas\b", exports)), "Missing Python init export")
    return {"format": "PE32+", "machine": "AMD64", "imports": imports}, imports


def elf(target: Path, extension: bool, env: dict[str, str]) -> dict[str, object]:
    headers = output(["readelf", "-hW", "-lW", str(target)], target.parent, env)
    require("ELF64" in headers and "X86-64" in headers, "Expected an x86_64 ELF")
    loader = re.search(r"Requesting program interpreter: ([^]]+)", headers)
    require((loader is None) == extension, "Unexpected ELF interpreter")
    dynamic = output(["readelf", "-dW", str(target)], target.parent, env)
    needed = re.findall(r"\(NEEDED\).*\[([^]]+)\]", dynamic)
    require(bool(needed), f"No ELF dependencies observed in {target.name}")
    search = re.findall(r"\((?:RUNPATH|RPATH)\).*\[([^]]*)\]", dynamic)
    require(not search, f"Unexpected embedded ELF library search path: {search}")
    versions = output(["readelf", "-VW", str(target)], target.parent, env)
    version_names = sorted(
        set(re.findall(r"\bName:\s+((?:GLIBCXX|GLIBC|CXXABI|GCC)_[\d.]+)", versions))
    )
    require(bool(version_names), "No ELF symbol version requirements observed")
    resolved = output(["ldd", str(target)], target.parent, env)
    require("not found" not in resolved, f"Missing ELF dynamic dependency in {target.name}")
    if extension:
        exports = output(["readelf", "--dyn-syms", "--wide", str(target)], target.parent, env)
        require(bool(re.search(r"\bPyInit__uvlas\b", exports)), "Missing Python init export")
    return {
        "format": "ELF64",
        "machine": "x86_64",
        "interpreter": loader.group(1) if loader else None,
        "needed": needed,
        "search_paths": search,
        "required_symbol_versions": version_names,
        "ldd": [line.split(" (0x", 1)[0].strip() for line in resolved.splitlines()],
    }


def notice_hashes(notices: Path, vcpkg_share: Path) -> dict[str, str]:
    installed = {original.parent.name: original for original in vcpkg_share.glob("*/copyright")}
    required_missing = sorted(set(REQUIRED_NOTICES) - installed.keys())
    require(not required_missing, f"Missing vcpkg source notices: {required_missing}")
    copied = {notice.stem: notice for notice in notices.glob("*.txt")}
    missing = sorted(installed.keys() - copied.keys())
    extra = sorted(copied.keys() - installed.keys())
    require(
        not missing and not extra, f"Notice inventory differs; missing={missing}, extra={extra}"
    )
    hashes: dict[str, str] = {}
    for package, original in sorted(installed.items()):
        value = digest(copied[package])
        require(value == digest(original), f"Notice bytes differ from vcpkg: {package}")
        hashes[package] = value
    return hashes


def verify_python(
    python: str,
    extension: Path,
    *,
    cli_imports: list[str],
    extension_imports: list[str],
    source_test: Path,
    stage: Path,
    env: dict[str, str],
) -> dict[str, object]:
    tested = document(
        output(
            [
                python,
                "-I",
                "-B",
                "-c",
                PYTHON_PROBE,
                str(extension),
                json.dumps(cli_imports),
                json.dumps(extension_imports),
            ],
            stage,
            env,
        )
    )
    version = tested.get("version")
    require(
        tested.get("implementation") == "cpython"
        and isinstance(version, str)
        and version.startswith("3.12."),
        "Expected CPython 3.12",
    )
    require(tested.get("numpy") == "2.3.5", "Unexpected NumPy version")
    require(tested.get("venv") is True and tested.get("isolated") is True, "Python not isolated")
    suffix = tested.get("extension_suffix")
    require(isinstance(suffix, str) and extension.name == "_uvlas" + suffix, "Wrong extension ABI")
    smoke = stage / "ownership.py"
    shutil.copyfile(source_test, smoke)
    owned = run([python, "-I", "-B", str(smoke), str(extension.parent)], stage, env)
    match = re.search(r"\bRan ([1-9]\d*) tests? in ", owned.stderr)
    if owned.returncode != 0 or match is None:
        raise RuntimeError(f"Ownership suite failed:\n{owned.stderr}")
    missing = stage / "empty"
    missing.mkdir()
    negative = run([python, "-I", "-B", str(smoke), str(missing)], stage, env)
    require(
        negative.returncode == 1 and "No module named '_uvlas'" in negative.stderr,
        f"Missing-extension probe did not reject the absent file:\n{negative.stderr}",
    )
    tested["ownership_test"] = {
        "script_sha256": digest(smoke),
        "cases": int(match.group(1)),
        "status": owned.returncode,
        "no_extension_status": negative.returncode,
        "no_extension_error": "ModuleNotFoundError: No module named '_uvlas'",
    }
    return tested


def exception_provider(trace: Path, caller: Path) -> Path:
    providers: set[Path] = set()
    with trace.open(encoding="utf-8", errors="replace") as lines:
        for line in lines:
            binding = GNU_EXCEPTION_BINDING.match(line)
            if binding is None or Path(binding[1]) != caller:
                continue
            provider = Path(binding[2])
            require(provider.is_absolute(), f"Non-absolute exception provider for {caller.name}")
            providers.add(provider.resolve(strict=True))
    if len(providers) != 1:
        raise RuntimeError(f"Expected one observed exception provider for {caller.name}")
    return providers.pop()


def audit_linux_exception_provider(
    python: Path, cli: Path, extension: Path, stage: Path, env: dict[str, str], *, expected: str
) -> dict[str, object]:
    # The byte offset in the normal CLI diagnostic establishes parser entry.
    source = stage / "invalid.json"
    source.write_bytes(b'{"schema":x}')
    runs = (
        ("cli", cli, [str(cli), "validate", str(source)], 2, "validate: JSON document at byte 10"),
        (
            "extension",
            extension,
            [str(python), "-I", "-B", "-c", PYTHON_EXCEPTION_PROBE, str(extension), str(source)],
            0,
            "Malformed JSON rejected.",
        ),
    )
    observations: dict[str, object] = {}
    for name, caller, command, status, message in runs:
        prefix = stage / f"binding-{name}"
        result = run(
            command, stage, {**env, "LD_DEBUG": "bindings", "LD_DEBUG_OUTPUT": str(prefix)}
        )
        report = result.stderr if name == "cli" else result.stdout
        require(
            result.returncode == status and report.strip() == message,
            f"{name} did not reach the malformed JSON path: {result.stderr.strip()}",
        )
        traces = list(stage.glob(f"{prefix.name}.*"))
        if len(traces) != 1:
            raise RuntimeError(f"Expected one glibc loader binding trace for {name}")
        provider = exception_provider(traces[0], caller)
        observed = digest(provider)
        require(
            observed == expected,
            f"{name} exception provider {provider.name} SHA-256 {observed} "
            f"does not match {expected}",
        )
        observations[name] = {"file": provider.name, "sha256": observed}
    return {"symbol": "std::runtime_error(char const*)@GLIBCXX_3.4.21", **observations}


class LinuxParserInputs(NamedTuple):
    build: Path
    expected_headers: dict[str, str]
    case: Path
    case_sha256: str


def make_dependency_paths(build: Path, target: Path) -> list[Path]:
    # GCC's -MD file is one Make rule; escaped line ends and names are significant.
    data = (build / target.with_suffix(".o.d")).read_text(encoding="utf-8")
    target_name, _, dependencies = data.replace("\\\n", " ").strip().partition(":")
    require(target_name == target.as_posix(), "Unexpected GNU Make parser dependency target")
    require("\n" not in dependencies, "Expected a single GNU Make parser dependency rule")
    names = re.findall(r"(?:\\.|[^\s\\])+", dependencies)
    return [build / re.sub(r"\\(.)", r"\1", name).replace("$$", "$") for name in names]


def verify_linux_parser_headers(
    inputs: LinuxParserInputs, producer: RuntimeArtifactProducer
) -> dict[str, object]:
    expected = inputs.expected_headers
    require(
        expected.keys() == GNU_PARSER_HEADERS, "Expected all four libstdc++ parser header hashes"
    )
    require(
        all(re.fullmatch(r"[0-9a-f]{64}", value) for value in expected.values()),
        "Invalid parser header SHA-256",
    )
    require(
        re.fullmatch(r"[0-9a-f]{64}", inputs.case_sha256) is not None,
        "Invalid external Case SHA-256",
    )
    build = inputs.build.resolve(strict=True)
    # Reuse the producer receipts against their original build outputs.
    verify_producer_artifacts(
        producer, build / producer.cli.relative_path, build / producer.extension.relative_path
    )
    require((build / GNU_PARSER_OBJECT).is_file(), "Missing built GNU parser object")
    observed: dict[str, str] = {}
    for path in make_dependency_paths(build, GNU_PARSER_OBJECT):
        if path.name not in expected:
            continue
        require(path.name not in observed, f"Multiple parser dependencies named {path.name}")
        observed[path.name] = digest(path)
        require(
            observed[path.name] == expected[path.name],
            f"Built parser header {path.name} SHA-256 {observed[path.name]} "
            f"does not match {expected[path.name]}",
        )
    require(observed.keys() == expected.keys(), "Missing a required GNU parser dependency")
    return {"object": GNU_PARSER_OBJECT.as_posix(), "header_sha256": observed}


def audit_linux_external_case(
    python: Path,
    cli: Path,
    extension: Path,
    stage: Path,
    env: dict[str, str],
    *,
    inputs: LinuxParserInputs,
) -> dict[str, object]:
    case = inputs.case.resolve(strict=True)
    require(digest(case) == inputs.case_sha256, "External Case SHA-256 differs before audit")
    result = output([str(cli), "validate", str(case)], stage, env)
    require(
        result.strip() == "Case validation succeeded.", "CLI did not validate the external Case"
    )
    result = output(
        [
            str(python),
            "-I",
            "-B",
            "-c",
            PYTHON_EXTERNAL_CASE_PROBE,
            str(extension),
            str(case),
            inputs.case_sha256,
        ],
        stage,
        env,
    )
    require(
        result.strip() == "External Case accepted.", "Python did not validate the external Case"
    )
    require(digest(case) == inputs.case_sha256, "External Case SHA-256 differs after audit")
    return {"case_sha256": inputs.case_sha256, "cli": "accepted", "python_file": "accepted"}


def audit(
    cli: Path,
    extension: Path,
    python: Path,
    source_test: Path,
    *,
    commit: str,
    producer_manifest: Path,
    dumpbin: str,
    linux_exception_provider_sha256: str | None = None,
    linux_parser_inputs: LinuxParserInputs | None = None,
) -> dict[str, object]:
    require(
        linux_parser_inputs is None or linux_exception_provider_sha256 is not None,
        "Combined audit requires an expected Linux exception provider",
    )
    if linux_exception_provider_sha256 is not None:
        require(sys.platform == "linux", "The exception provider audit requires Linux/glibc")
        require(
            re.fullmatch(r"[0-9a-f]{64}", linux_exception_provider_sha256) is not None,
            "Invalid expected exception provider SHA-256",
        )
    # The audit interpreter is outside the build; never resolve a venv's python symlink.
    require(python.is_absolute(), "Use an absolute path to the fresh venv interpreter")
    producer = read_producer_manifest(producer_manifest, commit)
    verify_producer_artifacts(producer, cli, extension)
    parser_headers = None
    if linux_parser_inputs is not None:
        parser_headers = verify_linux_parser_headers(linux_parser_inputs, producer)
    with tempfile.TemporaryDirectory(prefix="uvlas-runtime-audit-") as prepared:
        stage = Path(prepared).resolve()
        for source in (cli, extension, source_test):
            require(
                not stage.is_relative_to(source.resolve().parent),
                "Stage overlaps an input directory",
            )
        env = isolated_environment(stage)
        transport = stage / "artifacts"
        transport.mkdir()
        staged_cli = transport / cli.name
        staged_extension = transport / extension.name
        shutil.copy2(cli, staged_cli)
        shutil.copy2(extension, staged_extension)
        cli_hash = digest(staged_cli)
        extension_hash = digest(staged_extension)
        require(cli_hash == producer.cli.sha256, "CLI differs after copying")
        require(extension_hash == producer.extension.sha256, "Extension differs after copying")
        cli_imports: list[str] = []
        extension_imports: list[str] = []
        if os.name == "nt":
            cli_info, cli_imports = pe(staged_cli, dumpbin, extension=False)
            ext_info, extension_imports = pe(staged_extension, dumpbin, extension=True)
        elif sys.platform == "linux":
            cli_info = elf(staged_cli, extension=False, env=env)
            ext_info = elf(staged_extension, extension=True, env=env)
        else:
            raise RuntimeError("Expected Windows or Linux")
        banner = document(output([str(staged_cli), "--build-info"], stage, env))
        require(banner == producer.build_info, "CLI build info differs from the producer manifest")
        compiler = banner.get("compiler")
        required_prefix = "MSVC " if os.name == "nt" else "GNU 14."
        require(
            isinstance(compiler, str) and compiler.startswith(required_prefix),
            f"Expected {required_prefix.strip()} compiler evidence",
        )
        require(
            banner.get("dependency_repairs") == {"nlohmann_json_sax_allocation": 1},
            "Expected the current vcpkg JSON allocation repair",
        )
        negative = run([str(staged_cli), "--runtime-audit-invalid"], stage, env)
        require(
            negative.returncode == 2 and "Unknown option:" in negative.stderr,
            "CLI did not reject the invalid option",
        )
        python_info = verify_python(
            str(python),
            staged_extension,
            cli_imports=cli_imports,
            extension_imports=extension_imports,
            source_test=source_test,
            stage=stage,
            env=env,
        )
        report: dict[str, object] = {
            "platform": platform.system(),
            "producer_manifest_sha256": producer.sha256,
            "kernel_release": platform.release(),
            "platform_version": platform.version(),
            "isolated_environment_keys": sorted(env),
            "libc": platform.libc_ver() if sys.platform == "linux" else None,
            "cli": {"name": staged_cli.name, "sha256": cli_hash, **cli_info},
            "extension": {"name": staged_extension.name, "sha256": extension_hash, **ext_info},
            "relocated_cli": {"build_info": banner, "invalid_option_status": negative.returncode},
            "relocated_python": python_info,
        }
        if linux_exception_provider_sha256 is not None:
            report["linux_parser_exception_runtime"] = audit_linux_exception_provider(
                python,
                staged_cli,
                staged_extension,
                stage,
                env,
                expected=linux_exception_provider_sha256,
            )
        if linux_parser_inputs is not None:
            report["linux_parser_build"] = parser_headers
            report["linux_external_case"] = audit_linux_external_case(
                python, staged_cli, staged_extension, stage, env, inputs=linux_parser_inputs
            )
        return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("cli", "extension", "python", "ownership-test", "notices", "vcpkg-share"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--producer-manifest", required=True, type=Path)
    parser.add_argument("--dumpbin", default="dumpbin")
    parser.add_argument(
        "--linux-exception-provider-sha256",
        help="SHA-256 required for the libstdc++ provider of the JSON exception constructor",
    )
    parser.add_argument("--linux-parser-build", type=Path, help="GNU Make Release build directory")
    parser.add_argument(
        "--linux-parser-header",
        nargs=2,
        metavar=("NAME", "SHA256"),
        action="append",
        help="Expected header SHA-256; repeat for the four parser headers",
    )
    parser.add_argument("--linux-external-case", type=Path, help="Already prepared external Case")
    parser.add_argument("--linux-external-case-sha256")
    args = parser.parse_args()
    build: Path | None = args.linux_parser_build
    header_pairs: list[tuple[str, str]] | None = args.linux_parser_header
    case: Path | None = args.linux_external_case
    case_sha256: str | None = args.linux_external_case_sha256
    qualification = None
    if any(value is not None for value in (build, header_pairs, case, case_sha256)):
        if build is None or header_pairs is None or case is None or case_sha256 is None:
            parser.error("Supply build, four header hashes, external Case and its SHA-256 together")
        headers = dict(header_pairs)
        if len(headers) != len(header_pairs):
            parser.error("Supply each parser header only once")
        qualification = LinuxParserInputs(build, headers, case, case_sha256)
    result = audit(
        args.cli,
        args.extension,
        args.python,
        args.ownership_test,
        commit=args.commit,
        producer_manifest=args.producer_manifest,
        dumpbin=args.dumpbin,
        linux_exception_provider_sha256=args.linux_exception_provider_sha256,
        linux_parser_inputs=qualification,
    )
    result["vcpkg_notice_sha256"] = notice_hashes(args.notices, args.vcpkg_share)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
