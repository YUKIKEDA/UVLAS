"""Pin this launcher, verify the child mask and retain every measured sample."""
import ctypes
from ctypes import wintypes
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import statistics
import subprocess

output = Path("out/issue70-measure-cycles")
kernel = ctypes.WinDLL("kernel32", use_last_error=True)
kernel.GetCurrentProcess.restype = wintypes.HANDLE
kernel.GetProcessAffinityMask.argtypes = [wintypes.HANDLE, ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t)]
kernel.SetProcessAffinityMask.argtypes = [wintypes.HANDLE, ctypes.c_size_t]
process = kernel.GetCurrentProcess()
original, system = ctypes.c_size_t(), ctypes.c_size_t()
assert kernel.GetProcessAffinityMask(process, ctypes.byref(original), ctypes.byref(system))
selected = original.value & -original.value
assert selected and kernel.SetProcessAffinityMask(process, selected)
started = datetime.now(timezone.utc).isoformat()
try:
    run = subprocess.run([str(output / "probe.exe")], capture_output=True)
finally:
    assert kernel.SetProcessAffinityMask(process, original.value)
finished = datetime.now(timezone.utc).isoformat()
(output / "raw.log").write_bytes(run.stdout)
(output / "stderr.log").write_bytes(run.stderr)
report = {"started_utc": started, "finished_utc": finished, "exit": run.returncode,
          "launcher_affinity_mask": selected, "raw_sha256": hashlib.sha256(run.stdout).hexdigest()}
(output / "execution.json").write_text(json.dumps(report, indent=2) + "\n")
assert run.returncode == 0, run.stderr.decode(errors="replace")
text = run.stdout.decode("utf-8")
decoder = json.JSONDecoder()
rows = []
offset = 0
while offset < len(text):
    while offset < len(text) and text[offset].isspace():
        offset += 1
    if offset == len(text):
        break
    value, offset = decoder.raw_decode(text, offset)
    rows.append(value)
assert rows[0]["uncommitted_changes"] is False
assert rows[1]["process_affinity_mask"] == selected
samples = [row for row in rows if "case" in row]
expected = {"real8-pairs-both": 90, "real8-cluster-both": 9, "real96-pairs": 3, "complex6-pairs-both": 9}
summary = {"build_info": rows[0], "affinity": rows[1], "cases": {}}
for name, count in expected.items():
    case = [row for row in samples if row["case"] == name]
    assert len(case) == count and all(row["repeat_equal"] for row in case)
    assert len({(row["workspace_payload_bytes"], row["result_payload_bytes"], row["compared_words"]) for row in case}) == 1
    summary["cases"][name] = {
        "samples": count, "compared_words": case[0]["compared_words"],
        "workspace_payload_bytes": case[0]["workspace_payload_bytes"],
        "result_payload_bytes": case[0]["result_payload_bytes"],
        "median_stage_nanoseconds": [statistics.median(row["stages"][i]["nanoseconds"] for row in case) for i in range(4)],
        "minimum_stage_nanoseconds": [min(row["stages"][i]["nanoseconds"] for row in case) for i in range(4)],
        "maximum_stage_nanoseconds": [max(row["stages"][i]["nanoseconds"] for row in case) for i in range(4)],
    }
(output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps({"exit": run.returncode, "source_commit": rows[0]["commit"], "samples": len(samples), "cases": summary["cases"]}))
