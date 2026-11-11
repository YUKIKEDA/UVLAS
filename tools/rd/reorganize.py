"""Reorganize the R&D baseline into the family tree, as a pure function of the baseline.

    python tools/rd/reorganize.py plan             # report only
    python tools/rd/reorganize.py reissue          # bring the working tree to the current plan
    python tools/rd/reorganize.py verify --rev REV # REV equals transform(baseline), byte for byte

To move a family or a file, change relocation_plan.py, run `reissue`, then
`verify --rev HEAD` after the commit. Pages listed as MAINTAINED in the plan are
edited by hand and are never written or compared; fix their links by hand.

The transform reads every file of the baseline tag from git, never from the
working tree. It moves each file as relocation_plan.py says and rewrites, in
text files only:

  * Markdown and HTML link targets, re-resolved from the new location;
  * repository-root-relative path tokens (docs/..., tools/..., output/...);
  * in experiment scripts, the root and directory expressions, which become
    lookups through rd/rd_paths.py so that scripts no longer depend on where
    they or the records are stored;
  * SHA-256 tokens that name a baseline file whose bytes changed above.

Nothing else changes. `verify` proves that by recomputing the transform.
Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import posixpath
import re
import subprocess
import sys
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import relocation_plan as plan

REPO = Path(__file__).resolve().parents[2]
REGISTRY = "rd/registry.json"
TEXT_SUFFIXES = {".md", ".json", ".py", ".sh", ".txt", ".html", ".cmake", ".cjs", ".yaml", ".csv"}
LINK_SUFFIXES = {".md", ".html"}
PROLOGUE = (
    "import sys,pathlib;sys.path.insert(0,str(next(p for p in "
    "pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths"
)
DIR_MOVES = {
    "docs/verification": "rd",
    "tools/design-experiments": "rd",
    "tools/skills": ".agents/skills",
    "output/pdf": "docs/papers/pdf",
    "output": "docs/papers",
    "docs": "docs",
    "tools": "tools",
}

TOKEN = re.compile(r"(?<![\w/.\-])((?:docs|tools|output)/[\w./\-]*[\w/\-])")
MD_LINK = re.compile(r"(\]\(\s*<?)([^)\s>]+)")
HTML_LINK = re.compile(r"""((?:href|src)\s*=\s*["'])([^"']+)""")
LINK = re.compile(MD_LINK.pattern + "|" + HTML_LINK.pattern)
HASH = re.compile(r"(?<![0-9a-fA-F])[0-9a-fA-F]{64}(?![0-9a-fA-F])")
KEYED_HASH = re.compile(r'"([^"\n]+)"\s*:\s*"([0-9a-fA-F]{64})"')


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_rev(rev: str) -> dict[str, bytes]:
    raw = subprocess.run(
        ["git", "-C", str(REPO), "archive", "--format=tar", rev], check=True, capture_output=True
    ).stdout
    files: dict[str, bytes] = {}
    with tarfile.open(fileobj=io.BytesIO(raw)) as tar:
        for member in tar:
            if member.isfile():
                source = tar.extractfile(member)
                assert source is not None
                files[member.name] = source.read()
    return files


class Transform:
    def __init__(self, baseline: dict[str, bytes]):
        self.baseline = baseline
        self.moves: dict[str, str] = {}
        self.notes: dict[str, Counter[str]] = defaultdict(Counter)
        unmapped: list[str] = []
        for old in baseline:
            new = plan.new_path(old)
            if new is None:
                unmapped.append(old)
            else:
                self.moves[old] = new
        if unmapped:
            raise SystemExit("No relocation rule for:\n  " + "\n  ".join(unmapped))
        targets = Counter(self.moves.values())
        clashes = [t for t, n in targets.items() if n > 1]
        if clashes:
            raise SystemExit("Two files map to: " + ", ".join(clashes))
        self.old_dirs = {posixpath.dirname(p) for p in baseline}
        self.base_sha = {old: sha256(data) for old, data in baseline.items()}
        self.by_sha: dict[str, list[str]] = defaultdict(list)
        for old, digest in self.base_sha.items():
            self.by_sha[digest].append(old)
        # The baseline was redacted before it was committed (tools/rd/redact.py).
        # The record gives, for every file that changed, the hash it had before,
        # so a hash a record wrote down then still finds its file.
        record: dict[str, dict[str, dict[str, str]]] = (
            json.loads(baseline[plan.REDACTION_RECORD]) if plan.REDACTION_RECORD in baseline else {}
        )
        for old, entry in record.get("files", {}).items():
            self.by_sha[entry["pristine_sha256"]].append(old)
        # Files that left the repository for the evidence bundles: old path -> cited name.
        self.external = {old: e["name"] for old, e in record.get("externalized", {}).items()}
        # Hash before redaction -> hash now, for files that are not in the tree.
        self.outside = {
            e["pristine_sha256"]: e["sha256"]
            for group in ("externalized", "evidence")
            for e in record.get(group, {}).values()
            if e["pristine_sha256"] != e["sha256"]
        }
        self.stage1: dict[str, bytes] = {}
        self.final: dict[str, bytes] = {}
        self._stack: list[str] = []

    # ---- path mapping ---------------------------------------------------
    def new_dir(self, old: str) -> str | None:
        old = old.rstrip("/")
        for source, target in plan.TREE_MOVES:
            if old == source or old.startswith(source + "/"):
                return target + old[len(source) :]
        if old in DIR_MOVES:
            return DIR_MOVES[old]
        return old if old in self.old_dirs and plan.new_path(old + "/x") == old + "/x" else None

    def map_token(self, token: str) -> str | None:
        """Root-relative old path (file or directory) -> new path, keeping a trailing slash."""
        core = token.rstrip("/")
        slash = token[len(core) :]
        if core in self.external:
            return self.external[core] + slash
        if core in self.moves:
            return self.moves[core] + slash
        mapped = self.new_dir(core)
        if mapped is None and re.search(r"\.\w{1,5}$", core) and "{" not in core:
            # A file the baseline names but does not hold (an output never kept):
            # it gets the place the rules give its name.
            mapped = plan.new_path(core)
            if mapped is not None:
                self.notes["path of a file absent from the baseline, mapped by rule"][core] += 1
        return None if mapped is None else mapped + slash

    def map_link(self, old_file: str, target: str) -> str | None:
        if re.match(r"[a-zA-Z][a-zA-Z0-9+.\-]*:|#|//", target):
            return None
        path, sep, fragment = target.partition("#")
        if not path:
            return None
        resolved = posixpath.normpath(posixpath.join(posixpath.dirname(old_file), path))
        if resolved.startswith(".."):
            return None
        mapped = self.map_token(resolved + ("/" if path.endswith("/") else ""))
        if mapped is None:
            self.notes["unresolved link"][f"{old_file} -> {target}"] += 1
            return None
        relative = posixpath.relpath(
            mapped.rstrip("/"), posixpath.dirname(self.moves[old_file]) or "."
        )
        return relative + ("/" if mapped.endswith("/") else "") + sep + fragment

    # ---- stage 1: locations --------------------------------------------
    def rewrite_tokens(self, old_file: str, text: str) -> str:
        def replace(match: re.Match[str]) -> str:
            token = match.group(1)
            mapped = self.map_token(token)
            if mapped is None:
                self.notes["unresolved path token"][token] += 1
                return token
            return mapped

        return TOKEN.sub(replace, text)

    def rewrite_links(self, old_file: str, text: str) -> str:
        """Retarget links; rewrite path tokens only in the text between link targets."""
        pieces: list[str] = []
        position = 0
        for match in LINK.finditer(text):
            lead, target = (
                (match.group(1), match.group(2))
                if match.group(1)
                else (match.group(3), match.group(4))
            )
            pieces.append(self.rewrite_tokens(old_file, text[position : match.start()]))
            mapped = self.map_link(old_file, target)
            pieces.append(lead + (target if mapped is None else mapped))
            position = match.end()
        pieces.append(self.rewrite_tokens(old_file, text[position:]))
        return "".join(pieces)

    def rewrite_script(self, old_file: str, text: str) -> str:
        def quoted(path: str) -> str:
            return rf"""(['"]){path}/?\1"""

        text = re.sub(r"Path\(__file__\)\.resolve\(\)\.parents\[2\]", "rd_paths.ROOT", text)
        # A script beside this one: it may now live with another family.
        text = re.sub(r"Path\(__file__\)(?:\.resolve\(\))?\.with_name\(", "rd_paths.sibling(", text)
        text = re.sub(
            r"Path\((\w+(?:\.\w+)*\.__file__)\)(?:\.resolve\(\))?\.with_name\(",
            "rd_paths.sibling(",
            text,
        )
        # "a module loaded from the directory of this script" meant "an experiment module".
        text = re.sub(
            r"Path\((\w+)\.__file__\)\.resolve\(\)\.parent\s*==\s*\w+\.parent\b",
            r"rd_paths.is_experiment(\1.__file__)",
            text,
        )
        text = re.sub(
            r"Path\(__file__\)(?:\.resolve\(\))?\.parent\s*/", "rd_paths.EXPERIMENTS/", text
        )
        for own in set(
            re.findall(r"\b(\w+)\s*=\s*Path\(__file__\)(?:\.resolve\(\))?\s*(?=[;\n\r])", text)
        ):
            text = re.sub(rf"\b{own}\.with_name\(", "rd_paths.sibling(", text)
            text = re.sub(rf"\b{own}\.parent\s*/", "rd_paths.EXPERIMENTS/", text)
            text = re.sub(rf"\b{own}\.parents\[2\]", "rd_paths.ROOT", text)
        text = re.sub(r"\b\w+\s*/\s*" + quoted("docs/verification"), "rd_paths.RECORDS", text)
        text = re.sub(
            r"\b\w+\s*/\s*" + quoted("tools/design-experiments"), "rd_paths.EXPERIMENTS", text
        )
        text = re.sub(r"\b\w+\s*/\s*" + quoted("docs") + r"(?=\s*/)", "rd_paths.DOCS", text)
        # 'tools/design-experiments/'+name  and  f'docs/verification/{stem}-result.json'
        text = re.sub(
            r"""(['"])tools/design-experiments/\1\s*\+\s*([\w.]+(?:\([^()]*\))?(?:\.\w+)*)""",
            r"rd_paths.experiment(\2)",
            text,
        )
        text = re.sub(r"""(['"])docs/\1\s*\+\s*(\w+)""", r"rd_paths.doc(\2)", text)
        text = re.sub(
            r"""\bf(['"])docs/verification/((?:(?!\1).)*\{(?:(?!\1).)*)\1""",
            r"rd_paths.record(f\1\2\1)",
            text,
        )
        for leftover in re.finditer(
            r"""['"](?:docs/verification|tools/design-experiments|docs)/?['"]""", text
        ):
            self.notes["script directory literal left as is"][
                f"{old_file}: {text[max(0, leftover.start() - 40) : leftover.end() + 30]!r}"
            ] += 1
        if "rd_paths" in text or re.search(r"^(?:import|from)\s", text, re.MULTILINE):
            first = re.search(r"^(?:import|from)\s(?!__future__)", text, re.MULTILINE)
            position = first.start() if first else 0
            newline = "\r\n" if "\r\n" in text else "\n"
            text = text[:position] + PROLOGUE + newline + text[position:]
        return text

    def locate(self, old: str) -> bytes:
        if old in self.stage1:
            return self.stage1[old]
        data = self.baseline[old]
        suffix = posixpath.splitext(old)[1].lower()
        text = None
        if suffix in TEXT_SUFFIXES or posixpath.basename(old) == "CMakeLists.txt":
            try:
                text = data.decode("utf-8")
            except UnicodeDecodeError:
                self.notes["not UTF-8, moved unchanged"][old] += 1
        if text is not None:
            for before, after in plan.SCRIPT_PATCHES.get(old, ()):
                if before not in text:
                    raise SystemExit(f"{old}: patch source not found: {before!r}")
                text = text.replace(before, after)
            if suffix == ".py" and old.startswith(plan.OLD_EXPERIMENTS + "/"):
                text = self.rewrite_script(old, text)
            if posixpath.basename(old) == "CMakeLists.txt" and old.startswith(plan.OLD_EXPERIMENTS):
                text = text.replace(
                    '"${CMAKE_CURRENT_SOURCE_DIR}/../../.."',
                    '"${CMAKE_CURRENT_SOURCE_DIR}/../../../../.."',
                )
            text = (
                self.rewrite_links(old, text)
                if suffix in LINK_SUFFIXES
                else self.rewrite_tokens(old, text)
            )
            data = text.encode("utf-8")
        self.stage1[old] = data
        return data

    # ---- stage 2: hashes that name rewritten files ----------------------
    def settle(self, old: str) -> bytes:
        if old in self.final:
            return self.final[old]
        data = self.locate(old)
        if b"\0" in data[:4096]:
            self.final[old] = data
            return data
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            self.final[old] = data
            return data
        if not any(
            m.group(0).lower() in self.by_sha or m.group(0).lower() in self.outside
            for m in HASH.finditer(text)
        ):
            self.final[old] = data
            return data
        self._stack.append(old)
        keyed = {m.start(2): m.group(1) for m in KEYED_HASH.finditer(text)}

        def replace(match: re.Match[str]) -> str:
            token = match.group(0)
            lowered = token.lower()
            candidates = self.by_sha.get(lowered)
            if not candidates:
                # An evidence file whose text was redacted has a new hash.
                return self.outside.get(lowered, token)
            key = keyed.get(match.start())
            if key is not None and (key.startswith("build/") or (old, key) in plan.HASH_KEEP):
                # The pair names an evidence file, not the repository file that had
                # the same content. Its hash changes only if that file was redacted.
                self.notes["hash of an evidence file, not of the repository file"][
                    f"{old}: {key}"
                ] += 1
                return self.outside.get(lowered, token)
            if key is not None and len(candidates) > 1:
                named = [c for c in candidates if key in (c, self.moves[c], posixpath.basename(c))]
                candidates = named or candidates
            results: set[str] = set()
            for candidate in candidates:
                if candidate in self._stack:
                    self.notes["hash cycle, token kept"][f"{old} <-> {candidate}"] += 1
                    return token
                results.add(sha256(self.settle(candidate)))
            if len(results) != 1:
                self.notes["hash names several different files, token kept"][
                    f"{old}: {token[:12]}"
                ] += 1
                return token
            return results.pop()

        text = HASH.sub(replace, text)
        self._stack.pop()
        self.final[old] = text.encode("utf-8")
        return self.final[old]

    def run(self) -> dict[str, bytes]:
        sys.setrecursionlimit(20000)
        tree = {self.moves[old]: self.settle(old) for old in sorted(self.baseline)}
        registry = {
            "schema": "uvlas.rd-registry/1",
            "baseline": plan.BASELINE,
            "note": "Baseline path -> current path. Regenerate with tools/rd/reorganize.py; "
            "sha256 is the content at reorganization, baseline_sha256 the content at the tag.",
            "files": {
                old: {
                    "path": self.moves[old],
                    "baseline_sha256": self.base_sha[old],
                    "sha256": sha256(tree[self.moves[old]]),
                }
                for old in sorted(self.baseline)
            },
        }
        tree[REGISTRY] = (json.dumps(registry, ensure_ascii=False, indent=1) + "\n").encode("utf-8")
        return tree


def summary(transform: Transform, tree: dict[str, bytes]) -> str:
    moved = sum(1 for old, new in transform.moves.items() if old != new)
    rewritten = sum(
        1 for old, new in transform.moves.items() if tree[new] != transform.baseline[old]
    )
    lines = [
        (
            f"baseline files {len(transform.baseline)}, moved {moved}, "
            f"content rewritten {rewritten}, byte-identical {len(transform.baseline) - rewritten}"
        )
    ]
    for kind, items in sorted(transform.notes.items()):
        lines.append(f"{kind}: {len(items)} distinct, {sum(items.values())} occurrences")
        lines += [f"    {n:4}  {text}" for text, n in items.most_common()]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("mode", choices=("plan", "reissue", "verify"))
    parser.add_argument(
        "--rev",
        default="HEAD",
        help="verify: commit to check. reissue: commit holding the previous placement",
    )
    parser.add_argument("--report", type=Path, help="Write the full report to this file")
    args = parser.parse_args(argv)
    transform = Transform(read_rev(plan.BASELINE))
    tree = transform.run()
    report = summary(transform, tree)
    if args.report:
        args.report.write_text(report + "\n", encoding="utf-8")
    print(report.splitlines()[0])
    for kind, items in sorted(transform.notes.items()):
        print(f"  {kind}: {len(items)} distinct")

    # Files in the claim "this tree is the transform of the baseline": everything the
    # transform writes, except product files that never moved, pages maintained by
    # hand since the reorganization, and files removed on purpose afterwards.
    scope = {
        path
        for path in tree
        if not path.startswith(plan.UNMOVED_PREFIXES)
        and path not in plan.MAINTAINED
        and path not in plan.REMOVED
    }

    if args.mode == "reissue":
        # Bring the working tree to the current plan without touching anything outside
        # the scope. The previous placement comes from the registry committed at --rev.
        previous = read_rev(args.rev)
        before = json.loads(previous[REGISTRY].decode("utf-8"))["files"]
        left = {entry["path"] for entry in before.values()} - set(tree)
        written = removed = 0
        for path in sorted(left | set(plan.REMOVED)):
            if path in plan.MAINTAINED or path.startswith(plan.UNMOVED_PREFIXES):
                continue
            target = REPO / path
            if target.is_file():
                target.unlink()
                removed += 1
                directory = target.parent
                while directory != REPO and not any(directory.iterdir()):
                    directory.rmdir()
                    directory = directory.parent
        for path in sorted(scope):
            target = REPO / path
            if target.is_file() and target.read_bytes() == tree[path]:
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(tree[path])
            written += 1
        print(f"reissue: {written} files written, {removed} removed, {len(scope)} in scope")
    elif args.mode == "verify":
        actual = read_rev(args.rev)
        wrong = [p for p in sorted(scope) if actual.get(p) != tree[p]]
        stale = [
            old
            for old, new in transform.moves.items()
            if old != new and old in actual and old not in plan.MAINTAINED
        ]
        present = [p for p in plan.REMOVED if p in actual]
        for path in wrong[:20]:
            print(f"  differs from transform(baseline): {path}")
        for path in stale[:20]:
            print(f"  baseline path still present: {path}")
        for path in present[:20]:
            print(f"  removed by plan but present: {path}")
        print(
            f"verify {args.rev}: {len(scope) - len(wrong)}/{len(scope)} files equal "
            f"transform({plan.BASELINE[:7]}), {len(stale)} stale, "
            f"{len(plan.MAINTAINED)} maintained by hand, {len(plan.REMOVED)} removed by plan"
        )
        return 1 if wrong or stale or present else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
