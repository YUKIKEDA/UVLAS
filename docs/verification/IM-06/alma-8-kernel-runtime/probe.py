import ctypes
import json
import os
import pathlib
import sys

import _uvlas as uvlas

case = pathlib.Path("/case/external case.json")
libc = ctypes.CDLL(None, use_errno=True)
class OpenHow(ctypes.Structure):
    _fields_ = [(name, ctypes.c_uint64) for name in ("flags", "mode", "resolve")]

def openat2(label, parent, name, how):
    ctypes.set_errno(0)
    fd = libc.syscall(437, parent, name, ctypes.byref(how), ctypes.sizeof(how))
    errno = ctypes.get_errno()
    print(f"{label} openat2_return={fd} errno={errno}", flush=True)
    return fd

parent = openat2("native root", -100, b"/case", OpenHow(os.O_PATH | os.O_DIRECTORY | os.O_CLOEXEC, 0, 4))
opened = parent if parent >= 0 else os.open("/case", os.O_PATH | os.O_DIRECTORY)
child = openat2("native child" if parent >= 0 else "isolated child using ordinary directory fd", opened, b"wing.bdf", OpenHow(os.O_CLOEXEC | os.O_NONBLOCK | os.O_NOCTTY, 0, 12))
if child >= 0:
    os.close(child)
os.close(opened)

for kind, run in (
    ("file", lambda: uvlas.load_case(case)),
    ("dict", lambda: uvlas.case_from_dict(json.loads(case.read_text()), reference_root=case.parent)),
):
    try:
        result = run()
        print(f"python_{kind}={result.ok} failure={result.failure if not result.ok else None} code={result.failure.code if not result.ok else None}", flush=True)
    except Exception as error:
        print(f"python_{kind}={type(error).__name__}: {error}", flush=True)
        for name in ("code", "failure"):
            value = getattr(error, name, None)
            if value is not None:
                print(f"python_{kind}_{name}={value}", flush=True)
print("UVLAS_PYTHON_COMPLETE", flush=True)
