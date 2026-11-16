"""Regenerate the frozen 4096/.1 precision system in an isolated R&D copy."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import cast

RECORD = "rd/aero/wake-loads/records/uvlm-column-update-result.json"
RECORD_SHA256 = "ddd7649056e4253e5cafdf993322bbe1f8bf0c4ed82ec79c1ec7a9fe0ccc2fc9"
BASE_SHA256 = "c37f5c1a7d30c15777a84494917796d458b2f802a833de539be73a9f6ccfe5fb"
EFFECTIVE_SHA256 = "8ab8b7f97bc415b3f69ae4c6a70a4bfcd29a10931d441605302c68bed4988cbd"
RHS_SHA256 = "de478c4580657fb266760e41169a13270d7f799114ed58f1e7ffba69fa0aad6d"

PAYLOAD = r"""
import hashlib, json, sys
from pathlib import Path
import numpy as np
import scipy
sys.path.insert(0, str(Path.cwd() / "rd"))
import rd_paths
import check_uvlm_column_update as e
if "--probe" in sys.argv:
    print(json.dumps({"imported":True,"numpy":np.__version__,"scipy":scipy.__version__}),flush=True)
    raise SystemExit(0)
m,k=4096,0.1
op=e.real_operator(m)
base_hash=e.previous.hash_array(op["A"])
factor=e.lu_factor(np.array(op["A"],order="F"),overwrite_a=True,check_finite=False)
x=op["xc"]
v,integral_error=e.quad_vec(
    lambda t:e.O.D.density(1+t,x,5000.)*np.expm1(-2j*k*t),
    0.,160.,epsrel=1e-11,epsabs=1e-12,limit=4096)
b=np.column_stack([
    1j*16*k*mode["H"]-(8+1j*16*k*(x-mode["xp"]))*mode["alpha"] for mode in e.O.MODES])
actual=[base_hash,e.effective_hash(op,v),e.previous.hash_array(b)]
expected=json.loads(sys.argv[1])
if actual!=expected:
    raise RuntimeError("Regenerated full system differs from frozen base/effective/rhs hashes")
solve,denominator,ratio=e.update_solver(factor,v)
hi,lo,report=e.correction(op,v,b,solve)
if lo is None:
    raise RuntimeError("Frozen case did not produce TwoWord reference values")
rows=e.previous.evaluate(op,k,hi,lo)
record=rd_paths.RECORDS/"uvlm-column-update-result.json"
old=json.loads(record.read_text())
frozen=[r for r in old["continuous_cases"] if r["m"]==m and r["k"]==k]
if len(frozen)!=3 or [r["mode"] for r in frozen]!=[v["name"] for v in e.O.MODES]:
    raise RuntimeError("Frozen observation mode order changed")
def encoded(values):
    a=np.ascontiguousarray(values,dtype="<c16" if np.iscomplexobj(values) else "<f8")
    data=a.tobytes(order="C")
    if not np.all(np.isfinite(a)):
        raise RuntimeError("Nonfinite exported array")
    return dict(dtype="le-c128" if np.iscomplexobj(a) else "le-f64",
        shape=list(a.shape),order="C",sha256=hashlib.sha256(data).hexdigest(),
        payload_hex=data.hex())
arrays={name:encoded(op[name]) for name in ["xv","xc","fc","fu","mc0","mu0"]}
arrays.update({name:encoded(values) for name,values in [("v",v),("b",b),("hi",hi),("lo",lo)]})
observations=[]; differences=[]; modes=[]
for column,(mode,original,current) in enumerate(zip(e.O.MODES,frozen,rows)):
    scales=[float(s) for s in e.O.scales(mode)]
    modes.append(dict(**mode,scales=scales))
    if len(original["pressure_profile"])!=11 or len(current["pressure_profile"])!=11:
        raise RuntimeError("Frozen pressure window count changed")
    entries=[("force",None,None,original["loads"][0],current["loads"][0],scales[1]),("moment",None,None,original["loads"][1],current["loads"][1],scales[2])]
    for before,after in zip(original["pressure_profile"],current["pressure_profile"]):
        if [before["lower"],before["upper"]]!=[after["lower"],after["upper"]]:
            raise RuntimeError("Frozen window order or endpoints changed")
        entries.append(("pressure",before["lower"],before["upper"],before["pressure"],after["pressure"],scales[0]))
    for ordinal,(kind,lower,upper,value,generated,scale) in enumerate(entries):
        observation=dict(id=f"{mode['name']}/{ordinal}",mode=mode["name"],column=column,kind=kind,value=value,scale=scale)
        if kind=="pressure":observation.update(lower=lower,upper=upper)
        observations.append(observation)
        difference=abs(complex(*generated)-complex(*value))/max(scale,abs(complex(*value)))
        differences.append(dict(id=observation["id"],regenerated=generated,normalized_difference=float(difference)))
sources=[]
for module in list(sys.modules.values()):
    filename=getattr(module,"__file__",None)
    if filename and Path(filename).suffix==".py" and (Path.cwd()/"rd") in Path(filename).parents:
        path=Path(filename)
        sources.append(dict(id=path.stem,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),hint=path.relative_to(Path.cwd()).as_posix()))
e.budget()
result=dict(schema="uvlas.precision-system/1",id="precision-4096-system-v1",m=m,k=k,
    recipe=dict(signed_h="(sqrt((s-x)*(s-x)+5000*5000)/(s-x)-1)/(2*pi*5000)",
        pi_hex="0x1.921fb54442d18p+1",last_column_sources=["xv[-2]",161.0],
        update="A0 + v*e_last^T; rounded complex binary64 coefficients"),
    arrays=arrays,modes=modes,observations=observations,
    original=dict(id="uvlm-column-update-result",sha256=hashlib.sha256(record.read_bytes()).hexdigest(),hint=record.relative_to(Path.cwd()).as_posix(),base_sha256=actual[0],effective_sha256=actual[1],rhs_sha256=actual[2]),
    sources=sorted(sources,key=lambda item:item["id"]),
    generation=dict(id="precision-4096-system-v1-regenerated-reference",reference_planes_are_new=True,numpy=np.__version__,scipy=scipy.__version__,integral_error_estimate=float(integral_error),denominator=e.O.P.pair(denominator),relative_denominator=ratio,report=report,observations=differences,elapsed_seconds=e.time.monotonic()-e.START,peak_working_set_bytes=e.O.P.memory.memory_bytes(),budget=dict(seconds=300,peak_working_set_bytes=512*1024**2,outer_timeout_seconds=320)),
    qualification=dict(norm="L2",relative_residual_budget=1e-11,
        componentwise_backward_budget=1e-12,observation_budget=1e-8,
        scope="Reference recipe and regenerated planes; "
              "product full-system acceptance is a separate test"))
Path("system.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",newline="\n")
print(json.dumps(dict(stage="exported",
    array_bytes=sum(len(a["payload_hex"])//2 for a in arrays.values()),
    observations=len(observations),
    peak_working_set_bytes=result["generation"]["peak_working_set_bytes"],
    elapsed_seconds=result["generation"]["elapsed_seconds"])),flush=True)
"""


def require_hash(path: Path, expected: str) -> None:
    """Reject an altered frozen input before starting the numerical subprocess."""
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"Frozen input hash mismatch: {path.name}")


def source_hashes(root: Path) -> dict[str, str]:
    """Capture all frozen R&D files so the caller can verify the source copy is read-only."""
    return {
        p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted((root / "rd").rglob("*"))
        if p.is_file() and "__pycache__" not in p.parts
    }


def require_record_sources(root: Path) -> None:
    """Bind the generator and its frozen dependencies to the immutable record."""
    raw: object = json.loads((root / RECORD).read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("Frozen record must be an object")
    record = cast("dict[str, object]", raw)
    source = record.get("source_sha256")
    dependencies = record.get("dependency_sha256")
    if not isinstance(source, str) or not isinstance(dependencies, dict):
        raise ValueError("Frozen record source hashes are missing")
    expected = {"check_uvlm_column_update.py": source}
    for name, digest in cast("dict[str, object]", dependencies).items():
        if not isinstance(digest, str) or Path(name).name != name:
            raise ValueError("Invalid frozen dependency hash")
        expected[name] = digest
    for name, digest in expected.items():
        matches = list((root / "rd").rglob(name))
        if len(matches) != 1:
            raise ValueError(f"Frozen source must have one location: {name}")
        require_hash(matches[0], digest)


def run_export(root: Path, output: Path | None, *, probe: bool = False) -> None:
    """Run only the requested case and atomically refuse any existing destination."""
    if not probe and output is None:
        raise ValueError("An output path is required")
    if output is not None and output.exists():
        raise FileExistsError("Refusing to overwrite existing system fixture")
    require_hash(root / RECORD, RECORD_SHA256)
    require_record_sources(root)
    before = source_hashes(root)
    try:
        with tempfile.TemporaryDirectory(prefix="uvlas-precision-system-") as temporary:
            work = Path(temporary)
            shutil.copy2(root / "AGENTS.md", work / "AGENTS.md")
            shutil.copytree(root / "rd", work / "rd", ignore=shutil.ignore_patterns("__pycache__"))
            require_hash(work / RECORD, RECORD_SHA256)
            require_record_sources(work)
            command = [
                "uv",
                "run",
                "--no-project",
                "--python",
                "3.12",
                "--with",
                "numpy==2.3.5",
                "--with",
                "scipy==1.16.3",
                "--with",
                "scikit-fem==12.0.2",
                "python",
                "-c",
                PAYLOAD,
                json.dumps([BASE_SHA256, EFFECTIVE_SHA256, RHS_SHA256]),
            ]
            if probe:
                command.append("--probe")
            subprocess.run(command, cwd=work, check=True, timeout=10 if probe else 320)
            if not probe:
                if output is None:
                    raise ValueError("An output path is required")
                output.parent.mkdir(parents=True, exist_ok=True)
                with (work / "system.json").open("rb") as source, output.open("xb") as target:
                    shutil.copyfileobj(source, target)
                print(json.dumps({"sha256": hashlib.sha256(output.read_bytes()).hexdigest()}))
    finally:
        if before != source_hashes(root):
            raise RuntimeError("Original R&D files changed during isolated export")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not args.probe and args.output is None:
        parser.error("--output is required unless --probe")
    run_export(Path(__file__).resolve().parents[2], args.output, probe=args.probe)


if __name__ == "__main__":
    main()
