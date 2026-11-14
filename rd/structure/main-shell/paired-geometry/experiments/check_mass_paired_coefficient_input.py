"""Extra arbitrary coefficient-map probe, independently of the six-term family.

The process-local adapter is for this design experiment only. Production
reference data must own the immutable map, rather than use a module binding.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from fractions import Fraction as F
from types import SimpleNamespace
import hashlib,json,time,traceback
import numpy as np
import mass_paired_ancf as shell
import general_paired_fields as target
import exact_paired_finite as engine
from check_mass_paired_quality import inclusion
from paired_reference_admissibility import certify
ROOT=rd_paths.ROOT
COEFFICIENTS={(2,0):F(1,5),(0,2):F(1,10),(1,1):F(7,100),(3,0):F(1,30),(0,3):F(-1,50),
    (2,2):F(1,40),(4,0):F(7,1000),(0,4):F(-1,100),(3,1):F(3,250),(1,3):F(-1,125),(2,1):F(9,1000)}

def main():
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/mass-paired-coefficient-input-result.json';assert not dest.exists();records=[];stops=[];start=time.monotonic()
    original=shell.geometry.coefficients
    # A unique reference identity; no cache may cross a change of coefficient data.
    pars=(.3141592653589793,);shell.source_basis.cache_clear();shell.geometry.coefficients=lambda parameters:dict(COEFFICIENTS)
    engine.shell=shell;engine.coarse=target.Coarse
    engine.target=SimpleNamespace(fields=lambda p,h,parameters,variable:target.fields(p,h,parameters,variable,rates=(3.,5.,7.)))
    try:
        for n in [1,2]:
            for h in [.01,.0001]:
                tick=time.monotonic();model=shell.prepare(n,pars,h,True);proof=certify(pars,h,True);inc=inclusion(model)
                sol=engine.solve(model,h,True,order=model['degree']+3);obs=engine.observe(model,h,sol,True,order=model['degree']+5)
                checks=dict(stationarity=max(obs['stationarity'])<=1e-8,inclusion=max(max(a.values()) for a in inc)<=1e-8,
                    rcond=sol['rcond']>=1e-10,mass=sol['mass_min']>1e-10,resource=time.monotonic()-tick<=300.)
                if n==2:checks.update(energy=max(obs['error'])<=.01,position=max(obs['position'])<=.01)
                record=dict(n=n,h=h,observation=obs,inclusion=inc,rcond=sol['rcond'],mass_min=sol['mass_min'],stiffness_ranks=sol['ranks'],
                    mass_root_ranks=sol['mass_root_ranks'],checks=checks,certificate=proof,seconds=time.monotonic()-tick)
                records.append(record);print('COEFFICIENT_CASE',record,flush=True)
    except Exception as exc:stops.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print(stops[-1],flush=True)
    finally:shell.geometry.coefficients=original;shell.source_basis.cache_clear()
    comparisons=[]
    if len(records)==4:
        for n in [1,2]:
            a,b=[r for r in records if r['n']==n]
            for quantity in ['error','position']:
                ratio=np.maximum(b['observation'][quantity],1e-12)/np.maximum(a['observation'][quantity],1e-12)
                comparisons.append(dict(kind='thin_amplification',n=n,quantity=quantity,values=ratio.tolist(),passed=bool(np.all(ratio<=2))))
        for h in [.01,.0001]:
            a,b=[r for r in records if r['h']==h]
            for quantity in ['error','position']:
                comparisons.append(dict(kind='spatial_decrease',h=h,quantity=quantity,
                    passed=bool(np.all(np.array(b['observation'][quantity])<=np.maximum(a['observation'][quantity],1e-12)))))
    paths=[Path(__file__),Path(shell.__file__),Path(target.__file__),Path(engine.__file__),ROOT/'rd/structure/main-shell/paired-geometry/experiments/check_mass_paired_quality.py',ROOT/'rd/structure/main-shell/paired-geometry/experiments/paired_reference_admissibility.py']
    result=dict(passed=len(records)==4 and not stops and all(all(r['checks'].values()) for r in records) and all(c['passed'] for c in comparisons),
        coefficient_map={str(k):str(v) for k,v in COEFFICIENTS.items()},records=records,comparisons=comparisons,stops=stops,
        seconds=time.monotonic()-start,R4_complete=False,adapter_is_product_runtime=False,
        sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('COEFFICIENT_RESULT',result['passed'],flush=True)

if __name__=='__main__':main()
