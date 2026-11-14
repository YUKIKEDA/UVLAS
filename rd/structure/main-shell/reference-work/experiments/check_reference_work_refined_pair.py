"""Eight fixed shared-trace configurations, RW-3 against immutable strong C0 data."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,itertools,time
import numpy as np
import reference_work_refined_fit as fit
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def case(k,n,amp,h,out):
    tick=time.monotonic();a=fit.assemble(k,n,amp,h);sol=fit.solve(a);o=fit.observe(a,sol);ref=[fit.m.reference_h1(t,k,40) for t in fit.KIND_TARGETS]
    checks=dict(local_condition=sol['fit']['local_rcond']>=1e-10,global_condition=sol['fit']['global_rcond']>=1e-10,local_ranks=all(c['local_ranks']==[c['n']-6]*3 for c in a['local']),
        linear=max(sol['fit']['linear_residual'])<=1e-8,kkt=max(sol['fit']['kkt'])<=1e-8,direct_kkt=max(sol['fit']['direct_kkt'])<=1e-8,quadrature=max(abs(x-y) for x,y in zip(o['error'],sol['fit']['error']))<=1e-8,
        H1=all(x<=2*y*(1+1e-8) for x,y in zip(o['h1'],ref)),resources=time.monotonic()-tick<=180 and memory.memory_bytes()<=2*1024**3)
    path=out/f'k{k:g}-a{amp:g}-n{n}-h{h:g}.npz';np.savez_compressed(path,local=sol['local'],trace=sol['trace'],ref=a['original']['ref'],cells=a['original']['cells'])
    return dict(kind='RW3',kappa=k,n=n,amplitude=amp,h=h,fit=sol['fit'],observed=o,checks={k:bool(v) for k,v in checks.items()},seconds=time.monotonic()-tick,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path))
def main():
    dest=ROOT/'rd/structure/main-shell/reference-work/records/reference-work-refined-pair-result.json';assert not dest.exists();start=time.monotonic();records=[];stop=None
    assert json.loads((ROOT/'rd/structure/main-shell/reference-work/records/reference-work-v3-local-result.json').read_text())['all_local_checks']
    assert json.loads((ROOT/'rd/structure/main-shell/reference-work/records/reference-work-v3-local-audit.json').read_text())['all_reproduction_checks']
    out=ROOT/'build/analysis'/('reference-work-refined-pair-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    try:
        for k,amp,h in itertools.product([0.],[0.],[.01,.0001]):
            r=case(k,2,amp,h,out);records.append(r);(out/'checkpoint.json').write_text(json.dumps(dict(records=records),indent=2),encoding='utf-8');print(k,amp,h,r['observed'],r['checks'],r['seconds'],flush=True)
            if not all(r['checks'].values()) or time.monotonic()-start>900:raise RuntimeError('shared_necessary_gate')
    except Exception as exc:stop=repr(exc)
    baseline=json.loads((ROOT/'rd/structure/main-shell/hybrid-gradient/records/hybrid-gradient-shared-result.json').read_text());comparisons=[]
    for r in records:
        b=next(x for x in baseline['records'] if x['kind']=='baseline' and (x['kappa'],x['amplitude'],x['h'])==(r['kappa'],r['amplitude'],r['h']))
        comparisons.append(dict(kappa=r['kappa'],amplitude=r['amplitude'],h=r['h'],ratios=[x/y for x,y in zip(r['observed']['error'],b['observed']['error'])]))
    numeric=len(records)==2 and stop is None;improve=numeric and all(all(x<.99 for x in c['ratios']) for c in comparisons if c['kappa'])
    ps=[Path(__file__),Path(fit.__file__),ROOT/'rd/structure/main-shell/reference-work/experiments/reference_work_ancf_v3.py',ROOT/'rd/structure/main-shell/reference-work/records/reference-work-v3-protocol.md',ROOT/'rd/structure/main-shell/reference-work/records/reference-work-refined-preparation.json']
    dest.write_text(json.dumps(dict(records=records,stop_reason=stop,all_global_checks=numeric,all_curved_improved=bool(improve),comparisons=comparisons,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in ps},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
