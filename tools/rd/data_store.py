"""The large R&D evidence files: what the repository knows about them, and where they are.

The records in rd/ cite large files (mostly .npz arrays) as
"build/...": "<sha256>". The files are not in the repository. What is in the
repository is rd/evidence-index.json: every cited name with its SHA-256 and
size, and one of

    bundle         the asset of the GitHub release (`release` in the index) that holds it
    source         third-party material: where it is obtained from
    not_published  a build product: how it is made again

`tools/rd/evidence.py fetch` downloads the bundles into out/rd-evidence/ and
checks every file against the index. That directory is a cache: it can be
deleted and fetched again.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, cast

REPO = Path(__file__).resolve().parents[2]
INDEX = REPO / "rd" / "evidence-index.json"
CACHE = REPO / "out" / "rd-evidence"


def load_index() -> dict[str, Any]:
    return cast("dict[str, Any]", json.loads(INDEX.read_text(encoding="utf-8")))


def recorded_hashes() -> dict[str, str]:
    """Cited name -> SHA-256 the index records for it."""
    files = cast("dict[str, dict[str, Any]]", load_index()["files"])
    return {name: cast("str", entry["sha256"]) for name, entry in files.items()}


def cached(name: str) -> Path:
    """Where `fetch` puts the file for a cited name such as build/analysis/x.npz."""
    return CACHE / name
