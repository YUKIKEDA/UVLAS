"""Replace the user names and machine paths that the R&D material recorded.

The R&D files named the directory they were produced in and the accounts that
ran them: in logs, in records and in a few scripts. Those strings identify a
person and a machine and say nothing about the result. Before the material
enters the repository they are replaced:

    <repo>   the repository root, in any of the ways it was written
    <home>   a home directory, which carries the account name

No account name is written in this file: a name is recognized by where it
stands in a path.

Only text is changed. A file that is not valid UTF-8 is left as it is and
reported. Every changed file is listed in rd/redaction-record.json with the
SHA-256 of its original bytes and of its new bytes, so a hash recorded before
the replacement can still be matched to its file. Standard library only.
"""

from __future__ import annotations

import re

# The separator between path components as it appears in text: /, \ or an escaped \\.
_SEP = r"(?:\\\\|\\|/)+"
_NAME = r"[A-Za-z0-9_.\-]+"

# The repository root is recognized by the name of its directory, not by the
# directories above it, so that this file does not hold a machine's layout.
_ROOT = "UVLAS"

_RULES: list[tuple[re.Pattern[str], str]] = [
    (re.compile(rf"[A-Za-z]:{_SEP}(?:{_NAME}{_SEP})*{_ROOT}(?![A-Za-z0-9_\-])"), "<repo>"),
    (re.compile(rf"/mnt/[a-z]/(?:{_NAME}/)*{_ROOT}(?![A-Za-z0-9_\-])"), "<repo>"),
    (re.compile(rf"[A-Za-z]:{_SEP}Users{_SEP}{_NAME}", re.IGNORECASE), "<home>"),
    (re.compile(rf"/home/{_NAME}"), "<home>"),
]


def redact_text(text: str) -> tuple[str, int]:
    """The text with the replacements made, and how many were made."""
    total = 0
    for pattern, placeholder in _RULES:
        text, count = pattern.subn(placeholder, text)
        total += count
    return text, total


def redact_bytes(data: bytes) -> tuple[bytes, int]:
    """Redact UTF-8 text. Other content is returned unchanged with a count of -1."""
    if b"\0" in data[:8192]:
        return data, -1
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return data, -1
    redacted, count = redact_text(text)
    return (redacted.encode("utf-8"), count) if count else (data, 0)


def leftovers(data: bytes) -> list[str]:
    """Strings a published file must not hold. Empty when the file is clean."""
    found: list[str] = []
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return found
    for pattern, _ in _RULES:
        found += [match.group(0) for match in pattern.finditer(text)]
    return found
