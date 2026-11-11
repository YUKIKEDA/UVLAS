"""Fetch and check the large R&D evidence files (rd/evidence-index.json).

    python tools/rd/evidence.py summary          # what the index holds
    python tools/rd/evidence.py fetch            # download the bundles and check every file
    python tools/rd/evidence.py check            # check the files already in out/rd-evidence
    python tools/rd/evidence.py check-bundles DIR  # check bundle files in DIR against the index

`fetch` uses the GitHub CLI (`gh release download`). Every bundle is compared
with its SHA-256 before it is unpacked, and every unpacked file with its own.
Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
import tarfile
from collections import Counter
from pathlib import Path, PurePosixPath
from typing import Any, cast

sys.path.insert(0, str(Path(__file__).resolve().parent))
import data_store

REPOSITORY = "YUKIKEDA/UVLAS"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def summary(index: dict[str, Any]) -> int:
    files = cast("dict[str, dict[str, Any]]", index["files"])
    kinds: Counter[str] = Counter()
    sizes: Counter[str] = Counter()
    for entry in files.values():
        kind = next(k for k in ("bundle", "source", "not_published") if k in entry)
        kinds[kind] += 1
        sizes[kind] += cast("int", entry["bytes"])
    print(f"{len(files)} cited files; release {index['release']}")
    for kind in ("bundle", "source", "not_published"):
        print(f"  {kind:14} {kinds[kind]:5} files  {sizes[kind] / 2**30:6.2f} GiB")
    for name, info in cast("dict[str, dict[str, Any]]", index["bundles"]).items():
        print(f"  {name}: {info['files']} files, {info['bytes'] / 2**30:.2f} GiB")
    return 0


def check_bundles(index: dict[str, Any], directory: Path) -> int:
    failed = 0
    for name, info in cast("dict[str, dict[str, Any]]", index["bundles"]).items():
        path = directory / name
        if not path.is_file():
            print(f"missing bundle: {path}")
            failed += 1
        elif sha256(path) != info["sha256"]:
            print(f"bundle differs from the index: {path}")
            failed += 1
        else:
            print(f"bundle matches: {name}")
    return 1 if failed else 0


def unpack(index: dict[str, Any], directory: Path) -> None:
    files = cast("dict[str, dict[str, Any]]", index["files"])
    for bundle in cast("dict[str, Any]", index["bundles"]):
        with tarfile.open(directory / bundle) as tar:
            for member in tar:
                parts = PurePosixPath(member.name).parts
                if not member.isfile() or member.name not in files or ".." in parts:
                    raise SystemExit(f"{bundle}: unexpected member {member.name!r}")
                target = data_store.cached(member.name)
                target.parent.mkdir(parents=True, exist_ok=True)
                source = tar.extractfile(member)
                assert source is not None
                with source, target.open("wb") as out:
                    while chunk := source.read(1 << 20):
                        out.write(chunk)


def check_files(index: dict[str, Any]) -> int:
    files = cast("dict[str, dict[str, Any]]", index["files"])
    missing = wrong = good = 0
    for name, entry in files.items():
        if "bundle" not in entry:
            continue
        path = data_store.cached(name)
        if not path.is_file():
            missing += 1
        elif sha256(path) != entry["sha256"]:
            wrong += 1
            print(f"differs from the index: {name}")
        else:
            good += 1
    print(f"{good} files match the index, {missing} missing, {wrong} differ")
    return 1 if missing or wrong else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("mode", choices=("summary", "fetch", "check", "check-bundles"))
    parser.add_argument("directory", nargs="?", type=Path)
    args = parser.parse_args()
    index = data_store.load_index()
    if args.mode == "summary":
        return summary(index)
    if args.mode == "check":
        return check_files(index)
    if args.mode == "check-bundles":
        if args.directory is None:
            parser.error("check-bundles needs the directory that holds the bundles")
        return check_bundles(index, args.directory)
    downloads = data_store.CACHE.parent / "rd-evidence-download"
    downloads.mkdir(parents=True, exist_ok=True)
    done = subprocess.run(
        ["gh", "release", "download", str(index["release"]), "--repo", REPOSITORY,
         "--dir", str(downloads), "--pattern", "*.tar", "--skip-existing"],
        check=False,
    )  # fmt: skip
    if done.returncode != 0:
        return done.returncode
    if check_bundles(index, downloads) != 0:
        return 1
    unpack(index, downloads)
    return check_files(index)


if __name__ == "__main__":
    sys.exit(main())
