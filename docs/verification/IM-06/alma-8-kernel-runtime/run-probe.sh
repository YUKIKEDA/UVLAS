#!/bin/sh
set -u
export LC_ALL=C PYTHONUTF8=1 PATH=/bin HOME=/tmp
/bin/busybox uname -a
cd /case
/bin/busybox sha256sum /runtime/uvlas /work/qualified-tested/venv/lib/python3.12/site-packages/_uvlas*.so /usr/lib64/libc-2.28.so /usr/lib64/libstdc++.so.6.0.25 '/case/external case.json' /case/wing.bdf /case/wing.msh /case/gmsh.json /case/nastran.json
/runtime/uvlas --build-info
printf 'cli_build_info_exit=%s\n' "$?"
/runtime/uvlas validate '/case/external case.json'
printf 'cli_validate_exit=%s\n' "$?"
/work/qualified-tested/venv/bin/python -I -B /probe.py
printf 'python_exit=%s\n' "$?"
printf '%s\n' '{"schema":x}' > /tmp/malformed.json
LD_DEBUG=bindings /runtime/uvlas validate /tmp/malformed.json 2>/tmp/cli-bindings
printf 'cli_malformed_exit=%s\n' "$?"
LD_DEBUG=bindings /work/qualified-tested/venv/bin/python -I -B -c '
import _uvlas
r = _uvlas.load_case("/tmp/malformed.json")
assert not r.ok and r.failure.code == _uvlas.FailureCode.InvalidInput, r.failure
print("python_malformed_code=" + str(r.failure.code))
' 2>/tmp/python-bindings
printf 'python_malformed_probe_exit=%s\n' "$?"
/bin/busybox grep -E 'binding file (/runtime/uvlas|/work/qualified-tested/venv/.*/_uvlas[^ ]+).* to /lib64/libstdc[+][+][.]so[.]6 .*_ZNSt13runtime_errorC1EPKc' /tmp/cli-bindings /tmp/python-bindings
/bin/busybox realpath /lib64/libstdc++.so.6
/bin/busybox sha256sum /lib64/libstdc++.so.6
