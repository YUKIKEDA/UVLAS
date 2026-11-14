"""Traceability: where each requirement, decision and case is defined, tested and recorded.

    python tools/qa/trace.py show N4-06     # definition, tests, verification records, code
    python tools/qa/trace.py check          # gate
    python tools/qa/trace.py matrix         # rewrite docs/quality/traceability.md

The catalog is read from the documents; no list of IDs is kept by hand. In a
canonical document an ID is defined where a heading names it, where it starts
the first cell of a table row, or where it is followed by a full-width colon.
Evidence is found by the ID itself:

    tests         a Catch2 tag "[N4-06]" in tests/**/*.cpp
    records       "ids": [...] in docs/verification/**/*.json
    code          a mention in include/ or src/

`check` fails when a test or a record names an ID that no document defines,
when the canonical documents mention an ID that none of them defines, or when
docs/quality/traceability.md is not what `matrix` would write.
Standard library only.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import cast

REPO = Path(__file__).resolve().parents[2]
MATRIX = REPO / "docs" / "quality" / "traceability.md"

# Canonical documents, most authoritative first. The first definition of an ID wins;
# a later table row for the same ID is a restatement.
CANON = [
    "docs/requirements/product-requirements.md",
    "docs/requirements/acceptance-criteria.md",
    "docs/requirements/practical-acceptance-plan.md",
    "docs/requirements",
    "docs/specs",
    "docs/plans/implementation-start-review.md",
    "docs/quality",
]
KINDS = [
    ("requirement", r"(?:GEO|STR|MAT|PHY|SOL|RUN|QUAL|IO|API|ENV|DIST|OUT|BC|LOAD|PERF)-\d{2}"),
    ("decision", r"(?:SH-D|ST-|AE-|J1-D|MAT-D|N4-|IO4-|QA-)\d{2}"),
    ("acceptance", r"A5-\d{2}|PERF-FL-\d{2}|[VXP]-[A-Z][A-Z0-9]*(?:-[A-Z0-9]+)*"),
    ("package", r"IM-\d{2}"),
    ("risk", r"R-\d{2}"),
]
# IDs are ASCII; Japanese text may touch them on either side ("パッケージIM-00").
ANY = re.compile(r"(?<![A-Za-z0-9_-])(" + "|".join(p for _, p in KINDS) + r")(?![A-Za-z0-9_-])")
DEFINITION = re.compile(
    r"^\s*(?:#{1,6}\s+(?:\d+(?:\.\d+)*\.?\s+)?|\|\s*\**)("
    + "|".join(p for _, p in KINDS)
    + r")\**(?=[\s|:：、/])"
)
INLINE = re.compile(r"(?<![A-Za-z0-9_-])(" + "|".join(p for _, p in KINDS) + r")：")
TAG = re.compile(r"\[(" + "|".join(p for _, p in KINDS) + r")\]")
# Acceptance words that look like IDs but are ordinary text.
NOT_IDS = {"P-CURVED"}


def kind_of(identifier: str) -> str:
    for kind, pattern in KINDS:
        if re.fullmatch(pattern, identifier):
            return kind
    return "other"


def canon_files() -> list[Path]:
    seen: list[Path] = []
    for entry in CANON:
        path = REPO / entry
        files = [path] if path.is_file() else sorted(path.rglob("*.md")) if path.is_dir() else []
        for file in files:
            if file not in seen and file != MATRIX:
                seen.append(file)
    return seen


def read_catalog():
    definitions: dict[str, tuple[str, int, str]] = {}
    restated: dict[str, list[str]] = defaultdict(list)
    mentions: dict[str, set[str]] = defaultdict(set)
    for file in canon_files():
        label = file.relative_to(REPO).as_posix()
        for number, line in enumerate(
            file.read_text(encoding="utf-8", errors="replace").splitlines(), 1
        ):
            for match in ANY.finditer(line):
                mentions[match.group(1)].add(label)
            match = DEFINITION.match(line)
            # A heading defines every ID it names ("## 2. 最初のパッケージIM-00").
            named = (
                [m.group(1) for m in ANY.finditer(line)]
                if line.lstrip().startswith("#")
                else [match.group(1)]
                if match
                else []
            )
            # "MAT-D01：..." inside a paragraph is a definition as well.
            named += [m.group(1) for m in INLINE.finditer(line) if m.group(1) not in named]
            text = re.sub(r"\s+", " ", line.strip(" |#*")).strip()[:160]
            for identifier in named:
                if identifier in NOT_IDS:
                    continue
                if identifier in definitions:
                    restated[identifier].append(f"{label}:{number}")
                else:
                    definitions[identifier] = (label, number, text)
    for word in NOT_IDS:
        mentions.pop(word, None)
    return definitions, restated, mentions


def read_evidence():
    tests: dict[str, list[str]] = defaultdict(list)
    for file in sorted((REPO / "tests").rglob("*.cpp")):
        label = file.relative_to(REPO).as_posix()
        for number, line in enumerate(file.read_text(encoding="utf-8").splitlines(), 1):
            if "TEST_CASE" in line or "SCENARIO" in line:
                for tag in TAG.findall(line):
                    tests[tag].append(f"{label}:{number}")
    records: dict[str, list[str]] = defaultdict(list)
    root = REPO / "docs" / "verification"
    for file in sorted(root.rglob("*.json")) if root.is_dir() else []:
        label = file.relative_to(REPO).as_posix()
        try:
            data = json.loads(file.read_text(encoding="utf-8"))
        except ValueError:
            continue
        ids = cast("dict[str, list[str]]", data).get("ids", []) if isinstance(data, dict) else []
        for identifier in ids:
            records[str(identifier)].append(label)
    code: dict[str, list[str]] = defaultdict(list)
    for folder in ("include", "src"):
        for file in sorted((REPO / folder).rglob("*.[ch]pp")):
            label = file.relative_to(REPO).as_posix()
            for number, line in enumerate(file.read_text(encoding="utf-8").splitlines(), 1):
                for match in ANY.finditer(line):
                    code[match.group(1)].append(f"{label}:{number}")
    return tests, records, code


def history_definitions(missing: set[str]) -> dict[str, str]:
    """Where rd/ defines an ID that the canonical documents only mention."""
    found: dict[str, str] = {}
    for file in sorted((REPO / "rd").rglob("*.md")):
        if "records" in file.parts or not missing - found.keys():
            continue
        for number, line in enumerate(
            file.read_text(encoding="utf-8", errors="replace").splitlines(), 1
        ):
            match = DEFINITION.match(line)
            if match and match.group(1) in missing and match.group(1) not in found:
                found[match.group(1)] = f"{file.relative_to(REPO).as_posix()}:{number}"
    return found


def matrix_text(
    definitions: dict[str, tuple[str, int, str]],
    tests: dict[str, list[str]],
    records: dict[str, list[str]],
    code: dict[str, list[str]],
) -> str:
    lines = [
        "# 追跡表",
        "",
        (
            "`python tools/qa/trace.py matrix` が書く。手で直さない。IDの定義は文書から、"
            "試験はCatch2のタグから、検証記録は `docs/verification/` の `ids` から読む。"
            "一つのIDを調べるには `python tools/qa/trace.py show <ID>`。"
        ),
        "",
        (
            "「試験」「記録」の数が0であることは、そのIDがまだ実装も検証もされていないことを示す。"
            "数があることは、合格を意味しない。"
            "合否は検証記録と [品質保証計画](quality-plan.md) の状態で言う。"
        ),
        "",
    ]
    titles = {
        "requirement": "製品要求",
        "decision": "採用仕様の決定",
        "acceptance": "受入計画とケース",
        "package": "実装パッケージ",
        "risk": "残リスク",
    }
    for kind, title in titles.items():
        rows = sorted(i for i in definitions if kind_of(i) == kind)
        if not rows:
            continue
        lines += [
            f"## {title}",
            "",
            "| ID | 定義 | 試験 | 記録 | コード |",
            "|---|---|---:|---:|---:|",
        ]
        for identifier in rows:
            label, number, _ = definitions[identifier]
            link = "../" + label[len("docs/") :] if label.startswith("docs/") else "../../" + label
            lines.append(
                f"| {identifier} | [{Path(label).name}:{number}]({link}) | "
                f"{len(tests.get(identifier, []))} | {len(records.get(identifier, []))} | "
                f"{len(code.get(identifier, []))} |"
            )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("mode", choices=("show", "check", "matrix"))
    parser.add_argument("identifier", nargs="?")
    args = parser.parse_args(argv)
    definitions, restated, mentions = read_catalog()
    tests, records, code = read_evidence()

    if args.mode == "show":
        if not args.identifier:
            parser.error("show needs an ID")
        identifier = args.identifier
        if identifier not in definitions:
            where = history_definitions({identifier}).get(identifier)
            print(
                f"{identifier}: not defined in a canonical document"
                + (f"; rd/ defines it at {where}" if where else "")
            )
            return 1
        label, number, text = definitions[identifier]
        print(f"{identifier}  ({kind_of(identifier)})")
        print(f"  defined   {label}:{number}\n            {text}")
        for title, found in (
            ("restated", restated.get(identifier, [])),
            ("mentioned", sorted(mentions.get(identifier, set()) - {label})),
            ("tests", tests.get(identifier, [])),
            ("records", records.get(identifier, [])),
            ("code", code.get(identifier, [])),
        ):
            print(f"  {title:9} " + ("\n            ".join(found) if found else "none"))
        return 0

    text = matrix_text(definitions, tests, records, code)
    if args.mode == "matrix":
        MATRIX.parent.mkdir(parents=True, exist_ok=True)
        MATRIX.write_text(text, encoding="utf-8", newline="\n")
        print(f"wrote {MATRIX.relative_to(REPO).as_posix()}: {len(definitions)} IDs")
        return 0

    problems: list[str] = []
    status = (
        (MATRIX.parent / "case-status.md").read_text(encoding="utf-8")
        if (MATRIX.parent / "case-status.md").is_file()
        else ""
    )
    for identifier in sorted(definitions):
        if (
            kind_of(identifier) == "acceptance"
            and not identifier.startswith("A5-")
            and f"| {identifier} |" not in status
        ):
            problems.append(
                f"acceptance case without a row in docs/quality/case-status.md: {identifier}"
            )
    for source, found in (("test", tests), ("verification record", records)):
        for identifier, places in found.items():
            if identifier not in definitions:
                problems.append(f"{source} names an undefined ID {identifier}: {places[0]}")
    undefined = {i for i in mentions if i not in definitions}
    history = history_definitions(undefined)
    for identifier in sorted(undefined):
        place = min(mentions[identifier])
        hint = f" (rd/ defines it: {history[identifier]})" if identifier in history else ""
        problems.append(
            f"mentioned but not defined in a canonical document: {identifier} in {place}{hint}"
        )
    if not MATRIX.is_file() or MATRIX.read_text(encoding="utf-8") != text:
        problems.append("docs/quality/traceability.md is stale: run `trace.py matrix`")
    for line in problems:
        print(line)
    counts: defaultdict[str, int] = defaultdict(int)
    for identifier in definitions:
        counts[kind_of(identifier)] += 1
    print(
        f"{len(definitions)} IDs defined ({dict(counts)}), "
        f"{sum(1 for i in definitions if tests.get(i))} with tests, "
        f"{sum(1 for i in definitions if records.get(i))} with records, "
        f"{len(problems)} problems"
    )
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
