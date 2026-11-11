"""Check recorded local SHA256 references, never scientific acceptance. Stdlib only."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

HASH_MAPS = ("sources", "artifacts", "immutable", "status_sync_source")


class EvidenceError(Exception):
    pass


def check(root: Path, manifests: list[str], *, max_files: int = 100000,
          max_json_bytes: int = 32 * 1024 * 1024) -> dict:
    root = root.resolve()
    errors: list[str] = []
    actual: dict[Path, str] = {}
    expected: dict[Path, str] = {}
    inspected: set[Path] = set()
    statuses: list[dict] = []
    queue: list[Path] = []

    def resolve(name: str) -> Path:
        if not isinstance(name, str) or not name or Path(name).is_absolute():
            raise EvidenceError(f"Expected a nonempty relative file path: {name!r}")
        p = (root / name).resolve()
        if not p.is_relative_to(root):
            raise EvidenceError(f"Path escapes root: {name}")
        return p

    def label(p: Path) -> str:
        return p.relative_to(root).as_posix()

    def digest(p: Path) -> str:
        if p in actual:
            return actual[p]
        if len(actual) >= max_files:
            raise EvidenceError(f"File count limit exceeded: {max_files}")
        if not p.is_file():
            raise EvidenceError(f"Missing file: {label(p)}")
        before = p.stat()
        h = hashlib.sha256()
        with p.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                h.update(chunk)
        after = p.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise EvidenceError(f"File changed during check: {label(p)}")
        actual[p] = h.hexdigest()
        return actual[p]

    def reference(name: str, h: str, recurse: bool) -> None:
        if not isinstance(h, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", h):
            raise EvidenceError(f"Invalid SHA256 for {name!r}")
        p = resolve(name)
        h = h.lower()
        if p in expected and expected[p] != h:
            raise EvidenceError(f"Conflicting recorded hashes: {label(p)}")
        expected[p] = h
        if digest(p) != h:
            raise EvidenceError(f"SHA256 mismatch: {label(p)}")
        if recurse and p.suffix.lower() == ".json":
            queue.append(p)

    if max_files <= 0 or max_json_bytes <= 0:
        errors.append("Limits must be positive")
    if not root.is_dir():
        errors.append(f"Root directory does not exist: {root}")
    for name in manifests:
        try:
            queue.append(resolve(name))
        except (EvidenceError, OSError, ValueError) as exc:
            errors.append(str(exc))
    if not manifests:
        errors.append("At least one manifest is required")
    while queue and not errors:
        p = queue.pop()
        if p in inspected:
            continue
        try:
            digest(p)
            if p.stat().st_size > max_json_bytes:
                raise EvidenceError(f"JSON size limit exceeded: {label(p)}")
            raw = p.read_bytes()
            if hashlib.sha256(raw).hexdigest() != actual[p]:
                raise EvidenceError(f"Manifest changed before parsing: {label(p)}")
            data = json.loads(raw.decode("utf-8-sig"))
            if not isinstance(data, dict):
                raise EvidenceError(f"Manifest must be an object: {label(p)}")
            inspected.add(p)
            status = {k: data[k] for k in ("passed", "all_checks_pass") if k in data}
            if status:
                statuses.append({"manifest": label(p), "recorded_status": status})
            for key in HASH_MAPS:
                if key not in data:
                    continue
                mapping = data[key]
                if not isinstance(mapping, dict):
                    raise EvidenceError(f"{label(p)}: {key} must map paths to hashes")
                for name, h in mapping.items():
                    reference(name, h, recurse=key != "artifacts")
            if isinstance(data.get("artifact"), dict):
                artifact = data["artifact"]
                if "artifact" in artifact or "artifact_sha256" in artifact:
                    if not {"artifact", "artifact_sha256"} <= artifact.keys():
                        raise EvidenceError(f"Incomplete legacy artifact in {label(p)}")
                    reference(artifact["artifact"], artifact["artifact_sha256"], False)
        except (EvidenceError, OSError, ValueError, TypeError, RecursionError) as exc:
            errors.append(f"{label(p)}: {exc}")
    return {
        "integrity_passed": not errors,
        "scientific_acceptance": "not_evaluated",
        "scope": "Recorded local bytes only; not dependency completeness or reproduction",
        "root": str(root),
        "checked_files": len(actual),
        "inspected_manifests": len(inspected),
        "manifest_statuses": sorted(statuses, key=lambda x: x["manifest"]),
        "errors": errors,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--manifest", action="append", required=True,
                        help="Root-relative JSON path; may repeat")
    parser.add_argument("--output", help="New root-relative report path; never overwritten")
    parser.add_argument("--max-files", type=int, default=100000)
    parser.add_argument("--max-json-bytes", type=int, default=32 * 1024 * 1024)
    args = parser.parse_args(argv)
    report = check(args.root, args.manifest, max_files=args.max_files,
                   max_json_bytes=args.max_json_bytes)
    if args.output:
        try:
            root = args.root.resolve()
            output = (root / args.output).resolve()
            if Path(args.output).is_absolute() or not output.is_relative_to(root):
                raise EvidenceError("Output must be a relative path inside root")
            output.parent.mkdir(parents=True, exist_ok=True)
            with output.open("x", encoding="utf-8", newline="\n") as stream:
                json.dump(report, stream, ensure_ascii=False, indent=2)
                stream.write("\n")
        except (EvidenceError, OSError, ValueError) as exc:
            print(f"Cannot save report: {exc}", file=sys.stderr)
            return 2
        print(json.dumps({k: report[k] for k in
                          ("integrity_passed", "scientific_acceptance", "checked_files",
                           "inspected_manifests", "errors")}, ensure_ascii=False))
    else:
        print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["integrity_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
