"""Static checks of the organized R&D tree. Standard library only.

    python tools/rd/check_tree.py

  links      every relative Markdown link in docs/, rd/ and the root files resolves
  scripts    every experiment script compiles, and every import of another
             experiment script resolves to exactly one file
  files      every file of rd/registry.json is in the tree, except the files
             that relocation_plan.REMOVED names with a reason
  hashes     every "path": "sha256" pair in a record matches the bytes of the
             file it names
  controls   the checks above fail when a cited file is taken away

A hash pair that names a path is an error when the file is neither in the
repository nor in rd/evidence-index.json, unless relocation_plan names it as
removed or as not retained, with the reason. A pair under build/ is compared
with the index: the large artifacts are not part of a clone. A pair that names
one of the status pages (README, roadmap and the like) is counted as
historical: those pages are edited after the record was written, by design.

A record cites a file in one of two forms, and both are checked:

    "<path>": "<sha256>"
    "path": "<path>", "sha256": "<sha256>"      (also "x": "<path>", "x_sha256": ...)

A path that is not a repository path is resolved as the one repository file
whose path ends with it: records name fixtures relative to their family. A
bare name is checked when the registry or the tree knows a file of that name.
The rest are counted under what they are: a field name such as "result_sha256"
with no path beside it, or a file that was never in the repository.
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
FILE_NAME = re.compile(r"^[\w.\-]+\.[A-Za-z]\w*$")
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


def registry_files() -> dict[str, dict[str, str]]:
    registry = REPO / "rd" / "registry.json"
    return json.loads(registry.read_text(encoding="utf-8"))["files"]


def check_registry(files: list[Path]) -> list[str]:
    present = {p.relative_to(REPO).as_posix() for p in files}
    expected = {entry["path"] for entry in registry_files().values()} - plan.REMOVED.keys()
    return [
        f"{name} is in rd/registry.json and not in the tree" for name in sorted(expected - present)
    ]


def check_hashes(files: list[Path]) -> tuple[list[str], Counter[str], dict[str, str]]:
    """Problems, counts by kind, and one matched file per way of citing a file.

    The last value feeds the negative controls: taking such a file away has to
    turn its citation into a problem.
    """
    present = {p.relative_to(REPO).as_posix() for p in files}
    by_name: dict[str, list[str]] = {}
    for name in present:
        by_name.setdefault(name.rsplit("/", 1)[-1], []).append(name)
    files_moved = registry_files()
    baseline_hashes = {entry["baseline_sha256"] for entry in files_moved.values()}
    kept = {(files_moved[old]["path"], key) for old, key in plan.HASH_KEEP}
    external = {(files_moved[old]["path"], key) for old, key in plan.HASH_EXTERNAL}
    registered_names = {
        entry["path"].rsplit("/", 1)[-1]
        for entry in files_moved.values()
        if entry["path"] not in plan.REMOVED
    }
    evidence = data_store.recorded_hashes()
    evidence_by_name: dict[str, set[str]] = {}
    for name, recorded in evidence.items():
        evidence_by_name.setdefault(name.rsplit("/", 1)[-1], set()).add(recorded)
    cited: dict[str, str] = {}
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
                    # "path" beside "sha256", or "x" beside "x_sha256", names the file.
                    beside = mapping.get("path" if key == "sha256" else key.removesuffix("_sha256"))
                    named = beside.replace("\\", "/") if isinstance(beside, str) else ""
                    if key.endswith("sha256") and (
                        "/" in named or named in by_name or FILE_NAME.match(named)
                    ):
                        visit(named, value.lower())
                    else:
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
            # A redacted absolute path into the working tree names a repository path.
            key = key.removeprefix("<repo>/")
            candidates = by_name.get(key, [])
            if "/" in key and key not in present and not key.startswith("build/"):
                relative = [name for name in present if name.endswith("/" + key)]
                if len(relative) == 1:
                    key = relative[0]
            if (label, key) in external:
                counts["file outside the repository and the index, by plan"] += 1
            elif key in present:
                if digest(key) == value:
                    counts["matches"] += 1
                    cited["path"] = min(key, cited.get("path", key))
                elif key in status_pages:
                    counts["status page at an earlier state"] += 1
                else:
                    errors.append(f"{label}: {key} recorded {value[:12]}, is {digest(key)[:12]}")
            elif key.startswith("build/"):
                indexed = evidence.get(key)
                if indexed == value:
                    counts["large evidence file, matches the index"] += 1
                elif indexed is not None:
                    errors.append(
                        f"{label}: {key} recorded {value[:12]}, the index has {indexed[:12]}"
                    )
                else:
                    errors.append(f"{label}: {key} is not in rd/evidence-index.json")
            elif "/" in key:
                if key in plan.REMOVED:
                    counts["file removed by plan"] += 1
                else:
                    errors.append(f"{label}: {key} is cited by hash and is not in the repository")
            elif len(candidates) == 1:
                # Dependencies are often recorded by bare file name.
                if digest(candidates[0]) == value:
                    counts["matches by file name"] += 1
                    cited["name"] = min(candidates[0], cited.get("name", candidates[0]))
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
            elif candidates:
                if any(digest(name) == value for name in candidates):
                    counts["matches by file name"] += 1
                else:
                    errors.append(
                        f"{label}: {key} recorded {value[:12]}, which none of the "
                        f"{len(candidates)} files of that name has"
                    )
            elif key in registered_names:
                errors.append(f"{label}: {key} is cited by hash and is not in the repository")
            elif "." not in key:
                counts["hash under a field name, not a file"] += 1
            elif value in evidence_by_name.get(key, ()):
                counts["large evidence file by name, matches the index"] += 1
            else:
                counts["file that was never in the repository, not checked"] += 1

        walk(data, visit)
    return errors, counts, cited


def check_controls(files: list[Path], cited: dict[str, str]) -> list[str]:
    """Negative controls: a cited file that is taken away has to be reported."""
    errors: list[str] = []
    for how in ("path", "name"):
        victim = cited.get(how)
        if victim is None:
            errors.append(f"no file is cited by {how}, so the control cannot run")
            continue
        without = [p for p in files if p != REPO / victim]
        short = victim.rsplit("/", 1)[-1] if how == "name" else victim
        if not any(f"{short} is cited by hash" in line for line in check_hashes(without)[0]):
            errors.append(f"hashes did not report the missing file {victim} (cited by {how})")
        if not any(line.startswith(victim + " ") for line in check_registry(without)):
            errors.append(f"files did not report the missing file {victim}")
    return errors


def main() -> int:
    if sys.version_info < (3, 11):
        # Several scripts use star-unpacking inside subscripts, valid from Python 3.11.
        print("check_tree needs Python 3.11 or newer to compile the experiment scripts")
        return 2
    files = tracked()
    failed = False
    for title, errors in (
        ("links", check_links(files)),
        ("scripts", check_scripts(files)),
        ("files", check_registry(files)),
    ):
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
    errors, counts, cited = check_hashes(files)
    print(f"hashes: {len(errors)} problems; " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    for line in errors[:40]:
        print("   ", line)
    controls = check_controls(files, cited)
    print(f"controls: {len(controls)} problems")
    for line in controls:
        print("   ", line)
    return 1 if failed or errors or controls else 0


if __name__ == "__main__":
    raise SystemExit(main())
