"""First two full-DOF MMS-1 finite-thickness condensed solves."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
import mixed_mean_stress_finite as fit
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    out=ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-finite-pair-result.json';assert not out.exists();assert json.loads((ROOT/'rd/structure/main-shell/moment-enriched/records/moment-position-recovery-action-result.json').read_text())['all_checks_pass'];start=time.monotonic();directory=ROOT/'build/analysis'/('mixed-mean-stress-finite-pair-'+time.strftime('%Y%m%d-%H%M%S'));directory.mkdir();records=[];stops=[]
    for h in [.01,.0001]:
        try:
            tick=time.monotonic();a=fit.assemble(0.,2,0.,h);sol=fit.solve(a);obs=fit.observe(a,sol);checks=dict(local_ranks=all(c['local_ranks']==[552]*3 for c in a['local']),local_condition=sol['fit']['local_rcond']>=1e-10,global_condition=sol['fit']['global_rcond']>=1e-10,linear=max(sol['fit']['linear_residual'])<=1e-8,compressed_stationarity=max(sol['fit']['compressed_kkt'])<=1e-8,original_stationarity=max(sol['fit']['direct_kkt'])<=1e-8,independent=max(abs(x-y) for x,y in zip(sol['fit']['error'],obs['error']))<=1e-8,recovery_action=not obs['recovery_global_identity_allocated'],resources=time.monotonic()-tick<=180 and memory.memory_bytes()<=2*1024**3)
            path=directory/f'h{h:g}.npz';np.savez_compressed(path,local=sol['local'],trace=sol['trace'],ref=a['original']['ref'],cells=a['original']['cells']);r=dict(kappa=0.,n=2,amplitude=0.,h=h,fit=sol['fit'],observed=obs,local_ranks=[c['local_ranks'] for c in a['local']],checks={k:bool(v) for k,v in checks.items()},seconds=time.monotonic()-tick,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path));records.append(r);print('FINITE',h,r['fit'],r['observed'],r['checks'],r['seconds'],flush=True)
            (directory/'checkpoint.json').write_text(json.dumps(dict(records=records,stops=stops),indent=2),encoding='utf-8')
            if not all(r['checks'].values()):raise RuntimeError('MMS1_finite_pair_gate')
            if time.monotonic()-start>900:raise RuntimeError('MMS1_finite_pair_series_resources')
        except Exception as exc:stops.append(dict(h=h,reason=repr(exc)));print('STOP',repr(exc),flush=True);break
    paths=[Path(__file__),Path(fit.__file__),ROOT/'rd/structure/main-shell/mixed-mean-stress/experiments/mixed_mean_stress_element_v2.py',ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-finite-pair-protocol.md',ROOT/'rd/structure/main-shell/mixed-mean-stress/mixed-mean-stress-finite-formulation.md',ROOT/'rd/structure/main-shell/moment-enriched/records/moment-position-recovery-action-result.json']
    out.write_text(json.dumps(dict(records=records,stops=stops,all_checks_pass=len(records)==2 and not stops,kinds=fit.KINDS,main_shell_adopted=False,physical_PDE_acceptance=False,thin_accuracy_gate_evaluated=False,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),sources={p.relative_to(ROOT).as_posix():digest(p) for p in paths}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
