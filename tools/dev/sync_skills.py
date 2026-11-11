"""Copy .agents/skills to .claude/skills, or check that the copy is current.

.agents/skills is the source. Claude Code reads .claude/skills, so git tracks
both, and a commit that changes a skill changes both copies.

    python tools/dev/sync_skills.py          # copy
    python tools/dev/sync_skills.py --check  # exit 1 when the copies differ
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SOURCE = REPO / ".agents" / "skills"
COPY = REPO / ".claude" / "skills"


def files(root: Path) -> dict[str, bytes]:
    if not root.is_dir():
        return {}
    return {
        p.relative_to(root).as_posix(): p.read_bytes()
        for p in sorted(root.rglob("*"))
        if p.is_file() and "__pycache__" not in p.parts
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    source, copy = files(SOURCE), files(COPY)
    if args.check:
        different = sorted(
            name for name in source.keys() | copy.keys() if source.get(name) != copy.get(name)
        )
        for name in different:
            print(f"differs: {name}")
        print(f"{len(source)} skill files, {len(different)} differ")
        return 1 if different else 0
    if COPY.exists():
        shutil.rmtree(COPY)
    for name, data in source.items():
        target = COPY / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    print(f"copied {len(source)} files to {COPY.relative_to(REPO).as_posix()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
