"""Prepare and exercise all frozen V-PG cases through the product test."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

from tools.qa.vpg_case_reference import bundle


def main() -> None:
    executable = Path(sys.argv[1]).resolve(strict=True)
    with TemporaryDirectory(prefix="uvlas-vpg-") as directory:
        root = Path(directory)
        for relative, data in bundle().items():
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        subprocess.run([str(executable), "[.frozen-vpg]"], cwd=root, check=True)


if __name__ == "__main__":
    main()
