"""Build only the private probe, borrowing the fixed Release library read-only."""
import argparse
import json
import shlex
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--link-only", action="store_true")
args = parser.parse_args()
root = Path.cwd()
output = root / "out/issue70-measure-cycles"
build = root / "out/issue70-release"
entry = next(e for e in json.loads((build / "compile_commands.json").read_text())
             if e["file"].replace("\\", "/").endswith("/qa/eigen_arnoldi/main.cpp"))
command = [a[1:-1] if a.startswith('"') and a.endswith('"') else a
           for a in shlex.split(entry["command"], posix=False)]
assert "/O2" in command and "/DNDEBUG" in command
command = [str(output / "probe.cpp") if a.replace("\\", "/").endswith("/qa/eigen_arnoldi/main.cpp")
           else "/Fo" + str(output / "probe.obj") if a.startswith("/Fo")
           else "/Fd" + str(output / "probe.pdb") if a.startswith("/Fd") else a for a in command]
link = ["link.exe", "/nologo", str(output / "probe.obj"), str(build / "uvlas_numerics.lib"),
        str(build / "uvlas_core.lib"), "psapi.lib", "/OUT:" + str(output / "probe.exe"),
        "/PDB:" + str(output / "link.pdb")]
(output / "commands.json").write_text(json.dumps({"compile": command, "link": link}, indent=2))
for name, invocation, cwd in ([] if args.link_only else [("compile", command, build)]) + [("link", link, output)]:
    result = subprocess.run(invocation, cwd=cwd, capture_output=True)
    (output / f"{name}.log").write_bytes(result.stdout + result.stderr)
    print(f"{name} exit={result.returncode}")
    if result.returncode:
        print((result.stdout + result.stderr).decode(errors="replace"))
        raise SystemExit(result.returncode)
