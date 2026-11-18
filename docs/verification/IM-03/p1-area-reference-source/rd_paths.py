"""Location lookups for the R&D experiment scripts.

The scripts were written against one flat ``docs/verification`` and one flat
``tools/design-experiments``. They now ask this module instead, so neither the
scripts nor the records depend on where they are stored:

    rd_paths.ROOT                      repository root
    rd_paths.RECORDS / "x-result.json" the record with that name, wherever it is
    rd_paths.EXPERIMENTS / "check_x.py"
    rd_paths.DOCS / "roadmap.md"
    rd_paths.record("x-result.json")   the same as a root-relative string

Importing it also puts every ``experiments`` directory on ``sys.path``, which
keeps the imports between scripts of different families working. Files are
found by name in the tree, with ``registry.json`` as the record of where each
baseline file went. Moving a family directory needs no edit here.
"""
from __future__ import annotations

import fnmatch
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "AGENTS.md").is_file())
_OLD = {"records": "docs/verification", "experiments": "tools/design-experiments", "docs": "docs"}


def _registry() -> dict[str, str]:
    path = HERE / "registry.json"
    if not path.is_file():
        return {}
    return {old: entry["path"] for old, entry in
            json.loads(path.read_text(encoding="utf-8"))["files"].items()}


_REGISTRY = _registry()


def _family_dirs(kind: str) -> list[Path]:
    return sorted(p for p in HERE.rglob(kind) if p.is_dir())


class _Tree:
    """Stands for one of the old flat directories."""

    def __init__(self, kind: str):
        self.kind = kind
        self.old = _OLD[kind]

    def _find(self, name: str) -> Path | None:
        registered = _REGISTRY.get(f"{self.old}/{name}")
        if registered and (ROOT / registered).exists():
            return ROOT / registered
        if self.kind == "docs":
            hits = [p for p in (ROOT / "docs").rglob(name) if p.is_file()]
            hits += [p for p in HERE.rglob(name) if p.is_file() and p.parent.name != "records"]
        else:
            hits = [d / name for d in _family_dirs(self.kind) if (d / name).exists()]
            # Decks and figures sit beside the records of their family.
            hits += [d.parent / name for d in _family_dirs(self.kind) if (d.parent / name).exists()
                     and "/" in name]
        if len(hits) > 1:
            raise LookupError(f"{name} exists in several places: {[str(h) for h in hits]}")
        return hits[0] if hits else None

    def __truediv__(self, name) -> Path:
        name = Path(name).as_posix()
        found = self._find(name)
        if found is not None:
            return found
        prefix = f"{self.old}/{name}/"
        inside = [new for old, new in _REGISTRY.items() if old.startswith(prefix)]
        if inside:  # a directory of the old tree
            return _common_parent(inside)
        return self._new_file(name)

    def _new_file(self, name: str) -> Path:
        """A file that does not exist yet goes to the family of the running script."""
        main = getattr(sys.modules.get("__main__"), "__file__", None)
        if main:
            script = Path(main).resolve()
            if script.parent.name == "experiments" and HERE in script.parents:
                family = script.parent.parent
                return family / name if self.kind == "docs" else family / self.kind / name
        return HERE / "unsorted" / self.kind / name

    def glob(self, pattern: str) -> list[Path]:
        if self.kind == "docs":
            return sorted(p for p in (ROOT / "docs").rglob(pattern) if p.is_file())
        return sorted(p for d in _family_dirs(self.kind) for p in d.glob(pattern))

    def iterdir(self) -> list[Path]:
        return self.glob("*")

    def __repr__(self) -> str:
        return f"<rd_paths {self.kind}: files formerly in {self.old}>"


def _common_parent(paths: list[str]) -> Path:
    parts = [Path(p).parts for p in paths]
    common = []
    for level in zip(*parts):
        if len(set(level)) != 1:
            break
        common.append(level[0])
    return ROOT.joinpath(*common)


RECORDS = _Tree("records")
EXPERIMENTS = _Tree("experiments")
DOCS = _Tree("docs")


def _relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def record(name: str) -> str:
    return _relative(RECORDS / name)


def experiment(name: str) -> str:
    return _relative(EXPERIMENTS / name)


def doc(name: str) -> str:
    return _relative(DOCS / name)


def matching(kind: str, pattern: str) -> list[str]:
    """Root-relative paths of registered baseline files of one kind whose name matches."""
    prefix = _OLD[kind] + "/"
    return sorted(new for old, new in _REGISTRY.items()
                  if old.startswith(prefix) and fnmatch.fnmatch(old[len(prefix):], pattern))


for _directory in _family_dirs("experiments"):
    if str(_directory) not in sys.path:
        sys.path.append(str(_directory))


def sibling(name) -> Path:
    """The experiment file with this name (formerly beside every other script)."""
    return EXPERIMENTS / name


def is_experiment(path) -> bool:
    """True for a file in one of the experiments directories."""
    resolved = Path(path).resolve()
    return resolved.parent.name == "experiments" and HERE in resolved.parents
