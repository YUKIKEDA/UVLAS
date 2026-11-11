"""Conditional curved n2 then n4 finite-thickness projection acceptance."""
from pathlib import Path
import gc,hashlib,json,time
import numpy as np
import mixed_mean_stress_finite as fit
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    dest=ROOT/'docs/verification/mixed-mean-stress-finite-curved-result.json';assert not dest.exists()
    parent=ROOT/'docs/verification/mixed-mean-stress-finite-pair-disposition.json'
    assert json.loads(parent.read_text())['initial_package_completed']
    start=time.monotonic();records=[];stops=[]
    directory=ROOT/'build/analysis'/('mixed-mean-stress-finite-curved-'+time.strftime('%Y%m%d-%H%M%S'));directory.mkdir()
    for n in [2,4]:
        for amp,h in [(0.,.01),(0.,.0001),(.2,.01),(.2,.0001)]:
            try:
                tick=time.monotonic();print('START',n,amp,h,flush=True)
                a=fit.assemble(.4,n,amp,h);sol=fit.solve(a);obs=fit.observe(a,sol)
                checks=dict(local_ranks=all(c['local_ranks']==[552]*3 for c in a['local']),local_condition=sol['fit']['local_rcond']>=1e-10,global_condition=sol['fit']['global_rcond']>=1e-10,linear=max(sol['fit']['linear_residual'])<=1e-8,compressed_stationarity=max(sol['fit']['compressed_kkt'])<=1e-8,original_stationarity=max(sol['fit']['direct_kkt'])<=1e-8,independent=max(abs(x-y) for x,y in zip(sol['fit']['error'],obs['error']))<=1e-8,recovery_action=not obs['recovery_global_identity_allocated'],resources=time.monotonic()-tick<=180 and memory.memory_bytes()<=2*1024**3)
                path=directory/f'n{n}-a{amp:g}-h{h:g}.npz'
                np.savez_compressed(path,local=sol['local'],trace=sol['trace'],ref=a['original']['ref'],cells=a['original']['cells'])
                r=dict(kappa=.4,n=n,amplitude=amp,h=h,fit=sol['fit'],observed=obs,local_ranks=[c['local_ranks'] for c in a['local']],checks={k:bool(v) for k,v in checks.items()},seconds=time.monotonic()-tick,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path))
                records.append(r);print('CURVED',n,amp,h,r['fit'],obs,r['checks'],r['seconds'],flush=True)
                (directory/'checkpoint.json').write_text(json.dumps(dict(records=records,stops=stops),indent=2),encoding='utf-8')
                del a,sol;gc.collect()
                if not all(checks.values()):raise RuntimeError('MMS1_finite_curved_numerical_gate')
                if time.monotonic()-start>900:raise RuntimeError('MMS1_finite_curved_series_resources')
            except Exception as exc:stops.append(dict(n=n,amplitude=amp,h=h,reason=repr(exc)));print('STOP',repr(exc),flush=True);break
        if stops:break
    complete=len(records)==8 and not stops;accuracy=[]
    if complete:
        lookup={(r['n'],r['amplitude'],r['h']):r for r in records}
        for amp in [0.,.2]:
            for field in ['error','recovered_position_error']:
                for target in [1,2]:
                    errors={h:[lookup[(n,amp,h)]['observed'][field][target] for n in [2,4]] for h in [.01,.0001]}
                    ratio=errors[.0001][1]/max(errors[.01][1],1e-30)
                    accuracy.append(dict(amplitude=amp,field=field,target=fit.KINDS[target],errors=errors,thin_amplification=ratio,convergence=all(v[1]<v[0] for v in errors.values()),finest_accuracy=all(v[1]<=.01 for v in errors.values()),thin_uniformity=ratio<=2))
    paths=[Path(__file__),Path(fit.__file__),ROOT/'tools/design-experiments/mixed_mean_stress_element_v2.py',ROOT/'docs/verification/mixed-mean-stress-finite-curved-protocol.md',parent]
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_numerical_checks=complete,accuracy=accuracy,finite_projection_accuracy_gate=bool(complete and all(x['convergence'] and x['finest_accuracy'] and x['thin_uniformity'] for x in accuracy)),main_shell_adopted=False,physical_PDE_acceptance=False,independent_global_resolve=False,kinds=fit.KINDS,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),sources={p.relative_to(ROOT).as_posix():digest(p) for p in paths}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
