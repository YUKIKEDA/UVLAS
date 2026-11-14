"""Fetch and check the large R&D evidence files (rd/evidence-index.json).

    python tools/rd/evidence.py summary          # what the index holds
    python tools/rd/evidence.py fetch            # download the bundles and check every file
    python tools/rd/evidence.py check            # check the files already in out/rd-evidence
    python tools/rd/evidence.py check-bundles DIR  # check bundle files in DIR against the index
    python tools/rd/evidence.py pack BUNDLE LIST   # add cited files to the index, in a new bundle
    python tools/rd/evidence.py source NAME WHERE  # add a third-party file to the index
    python tools/rd/evidence.py not-published NAME HOW  # add a build product to the index

`fetch` uses the GitHub CLI (`gh release download`). Every bundle is compared
with its SHA-256 before it is unpacked, and every unpacked file with its own.

`pack`, `source` and `not-published` are for a cited file that tools/rd/check_tree.py reports as
missing from the index. They read the file from the working tree (build/...).
`pack` takes a JSON list of cited names, writes the bundle BUNDLE into
out/rd-evidence-bundles/ and adds the files and the bundle to the index. It
refuses a file that holds a user name or a machine path (tools/rd/redact.py).
The bundle is evidence only after its owner has uploaded it to the release and
`fetch` has passed. `source` records a third-party file that is not stored:
WHERE says where it is obtained from. `not-published` records a build product
that is not stored: HOW says what it is made from. Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tarfile
from collections import Counter
from pathlib import Path, PurePosixPath
from typing import Any, cast

sys.path.insert(0, str(Path(__file__).resolve().parent))
import data_store
import redact

REPOSITORY = "YUKIKEDA/UVLAS"
BUNDLES = data_store.CACHE.parent / "rd-evidence-bundles"


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


def write_index(index: dict[str, Any]) -> None:
    index["files"] = dict(sorted(cast("dict[str, Any]", index["files"]).items()))
    text = json.dumps(index, indent=1) + "\n"
    data_store.INDEX.write_text(text, encoding="utf-8", newline="\n")


def local_file(index: dict[str, Any], name: str) -> Path:
    """The working-tree file for a cited name that the index does not hold yet."""
    path = data_store.REPO / name
    if not name.startswith("build/") or ".." in PurePosixPath(name).parts:
        raise SystemExit(f"{name}: a cited name starts with build/")
    if name in index["files"]:
        raise SystemExit(f"{name}: already in the index")
    if not path.is_file():
        raise SystemExit(f"{name}: no such file in the working tree")
    return path


def pack(index: dict[str, Any], bundle: str, names: list[str]) -> int:
    target = BUNDLES / bundle
    if bundle in index["bundles"] or target.exists():
        raise SystemExit(
            f"{bundle}: a bundle of that name exists; a published bundle is not changed"
        )
    if not names or len(set(names)) != len(names):
        raise SystemExit("the list is empty or names a file twice")
    paths = {name: local_file(index, name) for name in sorted(names)}
    for name, path in paths.items():
        found = redact.leftovers(path.read_bytes())
        if found:
            raise SystemExit(f"{name}: holds a user name or a machine path: {found[0]!r}")
    BUNDLES.mkdir(parents=True, exist_ok=True)
    # The same files give the same bundle: sorted members, no owner, no time.
    with tarfile.open(target, "w", format=tarfile.PAX_FORMAT) as tar:
        for name, path in paths.items():
            member = tarfile.TarInfo(name)
            member.size = path.stat().st_size
            member.mode = 0o644
            with path.open("rb") as stream:
                tar.addfile(member, stream)
            index["files"][name] = {"sha256": sha256(path), "bytes": member.size, "bundle": bundle}
    index["bundles"][bundle] = {
        "sha256": sha256(target),
        "bytes": target.stat().st_size,
        "files": len(paths),
    }
    write_index(index)
    print(f"wrote {target.relative_to(data_store.REPO).as_posix()}: {len(paths)} files")
    print(f"upload it to the release {index['release']}, then run `evidence.py fetch`")
    return 0


def add_unstored(index: dict[str, Any], kind: str, name: str, text: str) -> int:
    """Record a file that is not stored: kind is `source` or `not_published`."""
    path = local_file(index, name)
    index["files"][name] = {"sha256": sha256(path), "bytes": path.stat().st_size, kind: text}
    write_index(index)
    print(f"added {name} as {kind}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument(
        "mode",
        choices=("summary", "fetch", "check", "check-bundles", "pack", "source", "not-published"),
    )
    parser.add_argument(
        "directory", nargs="?", help="check-bundles: DIR. pack: BUNDLE. source: NAME"
    )
    parser.add_argument("argument", nargs="?", help="pack: LIST. source: WHERE. not-published: HOW")
    args = parser.parse_args()
    index = data_store.load_index()
    if args.mode in ("pack", "source", "not-published"):
        if args.directory is None or args.argument is None:
            parser.error(f"{args.mode} needs two arguments")
        if args.mode != "pack":
            kind = args.mode.replace("-", "_")
            return add_unstored(index, kind, args.directory, args.argument)
        names = json.loads(Path(args.argument).read_text(encoding="utf-8"))
        return pack(index, args.directory, cast("list[str]", names))
    if args.mode == "summary":
        return summary(index)
    if args.mode == "check":
        return check_files(index)
    if args.mode == "check-bundles":
        if args.directory is None:
            parser.error("check-bundles needs the directory that holds the bundles")
        return check_bundles(index, Path(args.directory))
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
