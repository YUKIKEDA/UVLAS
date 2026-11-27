import hashlib
import json
from pathlib import Path
import re
import struct
import subprocess
import sys

image = Path(sys.argv[1])

def read(command):
    return subprocess.run(["debugfs", "-R", command, image], check=True, capture_output=True).stdout

entries = read("ls -p /recovery-proof").decode()
stage = re.findall(r"/([.]uvlas-checkpoint-[^/\n]+)/", entries)
assert len(stage) == 1
base = "/recovery-proof"

for number, values in ((12, {"first": (42,)}), (13, {"first": (42,), "second": (100, 200)})):
    parent = f"{base}/checkpoint-{number}"
    document = json.loads(read(f"cat {parent}/manifest.json"))
    assert document["schema"] == "uvlas.checkpoint/1" and document["byte_order"] == "little"
    assert set(document["fingerprints"]) == {"model", "basis", "material", "constraint", "algorithm", "precision", "version"}
    assert document["fingerprints"] == {
        name: "a7" + "00" * 31 if name == "version" else "00" * 32
        for name in document["fingerprints"]
    }
    assert [item["name"] for item in document["chunks"]] == list(values)
    for item in document["chunks"]:
        name = item["name"]
        expected = values[name]
        recorded = read(f"cat {parent}/chunk-{name}.bin")
        assert item["dtype"] == "u64" and item["shape"] == ([] if name == "first" else [2])
        assert item["bytes"] == len(recorded)
        assert recorded == struct.pack(f"<{len(expected)}Q", *expected)
        digest = hashlib.sha256(recorded).hexdigest()
        assert digest == item["sha256"]
        print(f"offline checkpoint-{number}/{name}: bytes={len(recorded)}, values={expected}, sha256={digest}")

stage_path = f"{base}/{stage[0]}"
stage_entries = read(f"ls -p {stage_path}").decode()
assert "manifest.json" not in stage_entries
assert read(f"cat {stage_path}/chunk-first.bin") == struct.pack("<Q", 42)
assert read(f"cat {stage_path}/chunk-second.bin") == struct.pack("<2Q", 100, 200)[:12]
print("offline staging: chunks=8/12 bytes, actual bytes match prefix, no manifest")
