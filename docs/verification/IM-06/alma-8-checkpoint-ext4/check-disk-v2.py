import hashlib
import json
from pathlib import Path
import re
import struct
import subprocess
import sys


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify(image):
    def read(command):
        return subprocess.run(
            ["debugfs", "-R", command, image], check=True, capture_output=True
        ).stdout

    base = "/recovery-proof"
    entries = read(f"ls -p {base}").decode()
    stage = re.findall(r"/([.]uvlas-checkpoint-[^/\n]+)/", entries)
    require(len(stage) == 1, f"{base}: expected one staging directory")
    fingerprints = dict.fromkeys(
        ("model", "basis", "material", "constraint", "algorithm", "precision"),
        "00" * 32,
    )
    fingerprints["version"] = "a7" + "00" * 31
    observations = []
    checkpoints = ((12, {"first": (42,)}), (13, {"first": (42,), "second": (100, 200)}))
    for number, values in checkpoints:
        parent = f"{base}/checkpoint-{number}"
        document = json.loads(read(f"cat {parent}/manifest.json"))
        require(
            document["schema"] == "uvlas.checkpoint/1" and document["byte_order"] == "little",
            f"{parent}: unexpected manifest schema or byte order",
        )
        require(document["fingerprints"] == fingerprints, f"{parent}: unexpected fingerprints")
        require(
            [item["name"] for item in document["chunks"]] == list(values),
            f"{parent}: unexpected chunk list",
        )
        for item in document["chunks"]:
            name = item["name"]
            expected = values[name]
            path = f"{parent}/chunk-{name}.bin"
            recorded = read(f"cat {path}")
            require(
                item["dtype"] == "u64" and item["shape"] == ([] if name == "first" else [2]),
                f"{path}: unexpected chunk type or shape",
            )
            require(item["bytes"] == len(recorded), f"{path}: wrong byte count")
            require(recorded == struct.pack(f"<{len(expected)}Q", *expected), f"{path}: wrong values")
            digest = hashlib.sha256(recorded).hexdigest()
            require(digest == item["sha256"], f"{path}: SHA-256 mismatch")
            observations.append(
                f"offline checkpoint-{number}/{name}: bytes={len(recorded)}, "
                f"values={expected}, sha256={digest}"
            )

    stage_path = f"{base}/{stage[0]}"
    stage_entries = read(f"ls -p {stage_path}").decode()
    require("manifest.json" not in stage_entries, f"{stage_path}: unexpected staging manifest")
    require(
        read(f"cat {stage_path}/chunk-first.bin") == struct.pack("<Q", 42),
        f"{stage_path}: wrong first staging chunk",
    )
    require(
        read(f"cat {stage_path}/chunk-second.bin") == struct.pack("<2Q", 100, 200)[:12],
        f"{stage_path}: wrong second staging chunk",
    )
    observations.append("offline staging: chunks=8/12 bytes, actual bytes match prefix, no manifest")
    return observations


if __name__ == "__main__":
    try:
        for line in verify(Path(sys.argv[1])):
            print(line)
    except (IndexError, KeyError, TypeError, ValueError, OSError, subprocess.CalledProcessError) as error:
        sys.exit(f"ext4 verification failed: {error}")
