"""Static checks of the organized R&D tree. Standard library only.

    python tools/rd/check_tree.py

  links      every relative Markdown link in docs/, rd/ and the root files resolves
  scripts    every experiment script compiles, and every import of another
             experiment script resolves to exactly one file
  hashes     every "path": "sha256" pair in a record whose path is a file of
             the repository matches that file's bytes

Hash pairs whose file is outside git (build/...) are counted, not checked: the
large artifacts are not part of a clone. A pair that names one of the status
pages (README, roadmap and the like) is counted as historical: those pages are
edited after the record was written, by design.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from collections import Counter
from collections.abc import Callable
from pathlib import Path
from typing import cast

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import data_store  # noqa: E402
import relocation_plan as plan  # noqa: E402

LINK = re.compile(r"\]\(\s*<?([^)\s>]+)")
HEX = re.compile(r"^[0-9a-fA-F]{64}$")
REMOVED_LINKS: list[str] = []  # links to files that relocation_plan.REMOVED names
EVIDENCE = data_store.recorded_hashes()
EVIDENCE_LINKS: list[str] = []  # links to evidence files, which are not in the repository
IMPORT = re.compile(r"^\s*(?:from\s+([A-Za-z_]\w*)|import\s+([A-Za-z_][\w, ]*))", re.MULTILINE)


def tracked() -> list[Path]:
    out = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
        check=True,
        capture_output=True,
    ).stdout
    return [
        REPO / name for name in out.decode("utf-8").split("\0") if name and (REPO / name).is_file()
    ]


def check_links(files: list[Path]) -> list[str]:
    errors: list[str] = []
    for path in files:
        if path.suffix != ".md" or not path.is_file():
            continue
        text = re.sub(
            r"```.*?```", "", path.read_text(encoding="utf-8", errors="replace"), flags=re.DOTALL
        )
        for target in LINK.findall(text):
            if re.match(r"[a-zA-Z][a-zA-Z0-9+.\-]*:|#|//", target):
                continue
            file_part = target.split("#", 1)[0]
            if not file_part or (path.parent / file_part).exists():
                continue
            resolved = (path.parent / file_part).resolve()
            if (
                resolved.is_relative_to(REPO)
                and resolved.relative_to(REPO).as_posix() in plan.REMOVED
            ):
                REMOVED_LINKS.append(f"{path.relative_to(REPO).as_posix()}: {target}")
            elif (
                resolved.is_relative_to(REPO) and resolved.relative_to(REPO).as_posix() in EVIDENCE
            ):
                EVIDENCE_LINKS.append(f"{path.relative_to(REPO).as_posix()}: {target}")
            else:
                errors.append(f"{path.relative_to(REPO).as_posix()}: {target}")
    return errors


def check_scripts(files: list[Path]) -> list[str]:
    scripts = [p for p in files if p.suffix == ".py" and p.parent.name == "experiments"]
    modules = Counter(p.stem for p in scripts)
    errors = [f"module name used by {n} scripts: {name}" for name, n in modules.items() if n > 1]
    for path in scripts:
        source = path.read_text(encoding="utf-8", errors="replace")
        try:
            compile(source, str(path), "exec")
        except SyntaxError as exc:
            errors.append(f"{path.relative_to(REPO).as_posix()}: {exc}")
    return errors


def check_hashes(files: list[Path]) -> tuple[list[str], Counter[str]]:
    present = {p.relative_to(REPO).as_posix() for p in files}
    by_name: dict[str, list[str]] = {}
    for name in present:
        by_name.setdefault(name.rsplit("/", 1)[-1], []).append(name)
    registry = REPO / "rd" / "registry.json"
    baseline_hashes = {
        entry["baseline_sha256"]
        for entry in json.loads(registry.read_text(encoding="utf-8"))["files"].values()
    }
    files_moved = json.loads(registry.read_text(encoding="utf-8"))["files"]
    kept = {(files_moved[old]["path"], key) for old, key in plan.HASH_KEEP}
    evidence = data_store.recorded_hashes()
    digests: dict[str, str] = {}
    errors: list[str] = []
    counts: Counter[str] = Counter()
    status_pages: set[str] = set()

    def digest(name: str) -> str:
        if name not in digests:
            digests[name] = hashlib.sha256((REPO / name).read_bytes()).hexdigest()
        return digests[name]

    def walk(node: object, visit: Callable[[str, str], None]) -> None:
        if isinstance(node, dict):
            mapping = cast("dict[str, object]", node)
            for key, value in mapping.items():
                if isinstance(value, str) and HEX.match(value):
                    visit(key.replace("\\", "/"), value.lower())
                else:
                    walk(value, visit)
            entries = mapping.get("mutable_status_entries")
            if isinstance(entries, list):
                status_pages.update(e for e in cast("list[object]", entries) if isinstance(e, str))
        elif isinstance(node, list):
            for value in cast("list[object]", node):
                walk(value, visit)

    records = [p for p in files if p.suffix == ".json" and "records" in p.parts]
    parsed: list[tuple[Path, object]] = []
    for path in records:
        try:
            parsed.append((path, json.loads(path.read_bytes().decode("utf-8-sig"))))
        except ValueError as exc:
            errors.append(f"{path.relative_to(REPO).as_posix()}: not JSON: {exc}")
    for _, data in parsed:  # first pass: learn which pages are status pages
        walk(data, lambda _key, _value: None)
    for path, data in parsed:
        label = path.relative_to(REPO).as_posix()

        def visit(key: str, value: str, label: str = label) -> None:
            if key in present:
                if digest(key) == value:
                    counts["matches"] += 1
                elif key in status_pages:
                    counts["status page at an earlier state"] += 1
                else:
                    errors.append(f"{label}: {key} recorded {value[:12]}, is {digest(key)[:12]}")
            elif "/" not in key and len(by_name.get(key, ())) == 1:
                # Dependencies are often recorded by bare file name.
                if digest(by_name[key][0]) == value:
                    counts["matches by file name"] += 1
                elif (label, key) in kept:
                    counts["names a copy under build/, by plan"] += 1
                elif value not in baseline_hashes:
                    # Already different at the baseline: the record cites a version
                    # of the file that was later replaced.
                    counts["earlier version of a file, as at the baseline"] += 1
                else:
                    errors.append(
                        f"{label}: {key} recorded {value[:12]}, is {digest(by_name[key][0])[:12]}"
                    )
            elif key.startswith("build/"):
                indexed = evidence.get(key)
                if indexed is None:
                    errors.append(f"{label}: {key} is not in rd/evidence-index.json")
                elif indexed == value:
                    counts["large evidence file, matches the index"] += 1
                else:
                    errors.append(
                        f"{label}: {key} recorded {value[:12]}, the index has {indexed[:12]}"
                    )
            else:
                counts["name without a repository path, not checked"] += 1

        walk(data, visit)
    return errors, counts


def main() -> int:
    if sys.version_info < (3, 11):
        # Several scripts use star-unpacking inside subscripts, valid from Python 3.11.
        print("check_tree needs Python 3.11 or newer to compile the experiment scripts")
        return 2
    files = tracked()
    failed = False
    for title, errors in (("links", check_links(files)), ("scripts", check_scripts(files))):
        note = (
            f"; {len(REMOVED_LINKS)} links to files removed by plan, "
            f"{len(EVIDENCE_LINKS)} links to evidence files outside the repository"
            if title == "links"
            else ""
        )
        print(f"{title}: {len(errors)} problems{note}")
        for line in errors[:40]:
            print("   ", line)
        failed |= bool(errors)
    errors, counts = check_hashes(files)
    print(f"hashes: {len(errors)} problems; " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    for line in errors[:40]:
        print("   ", line)
    return 1 if failed or errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
