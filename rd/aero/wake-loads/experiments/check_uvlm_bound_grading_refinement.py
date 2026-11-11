"""Predeclared 128-to-256 refinement after the retained 64-to-128 failure."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import hashlib
import json
import time
from pathlib import Path
import check_uvlm_bound_grading as grading

ROOT=rd_paths.ROOT
START=time.monotonic()


def main():
    diagnostic=grading.diagnostic
    previous=json.loads((ROOT/'rd/aero/wake-loads/records/uvlm-bound-grading-result.json').read_text())
    mesh=grading.install_grid(256,True);rows=[];changes=[];sensitivities=[]
    for k in [.1,.5,1.]:
        kw,error=diagnostic.continuous(256,k,160.)
        row,_,_=diagnostic.solve(256,k,kw)
        row['mesh']='cosine-physical-panels-v1';row['integration_error_estimate']=error
        rows.append(row)
        coarse=next(r for r in previous['harmonic_cases'] if r['m']==128 and r['k']==k)
        changes.append(dict(k=k,relative=diagnostic.previous.relative(coarse,row)))
        tighter,estimate=diagnostic.continuous(256,k,160.,tight=True)
        tightrow,_,_=diagnostic.solve(256,k,tighter)
        sensitivities.append(dict(k=k,relative=diagnostic.previous.relative(row,tightrow),tight_error_estimate=estimate))
    checks=dict(positive_cells=mesh['positive_cells'],
        complex_accuracy=all(r['complex_relative']<=.02 for r in rows),
        amplitude_accuracy=all(r['amplitude_relative']<=.02 for r in rows),
        phase_accuracy=all(abs(r['phase_degrees'])<=1 for r in rows),
        spatial_change=all(c['relative']<=.005 for c in changes),
        nonpenetration=all(r['nonpenetration_relative']<=1e-11 for r in rows),
        integral_sensitivity=all(c['relative']<=1e-8 for c in sensitivities))
    source=Path(__file__);names=['check_uvlm_bound_grading.py',*previous['dependency_sha256']]
    result=dict(protocol='rd/aero/wake-loads/records/uvlm-trailing-edge-diagnostic.md#7',
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={name:hashlib.sha256(rd_paths.sibling(name).read_bytes()).hexdigest() for name in names},
        baseline_result_sha256=hashlib.sha256((ROOT/'rd/aero/wake-loads/records/uvlm-bound-grading-result.json').read_bytes()).hexdigest(),
        geometry_id=diagnostic.previous.geometry.GEOMETRY,mesh=mesh,harmonic_cases=rows,
        refinement_changes=changes,integral_sensitivity=sensitivities,checks=checks,
        all_criteria_pass=bool(all(checks.values())),numpy=diagnostic.np.__version__,scipy=diagnostic.scipy.__version__,
        U=8.,rho=1.2,H=.001,c=1.,span=10000.,length=160.,
        seconds=time.monotonic()-START,peak_process_working_bytes=diagnostic.previous.memory.memory_bytes(),
        limits='frequency-only continuous wake and exact derivative; no complete transient method or Kutta pressure proof')
    (ROOT/'rd/aero/wake-loads/records/uvlm-bound-grading-refinement-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=diagnostic.previous.core.json_scalar)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,changes=changes,rows=rows,seconds=result['seconds'])),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':main()
