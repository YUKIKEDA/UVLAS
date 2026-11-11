"""Independent uncondensed reference for the two thin curved n2 cases."""
from pathlib import Path
import gc,hashlib,json,time
import numpy as np
import check_mixed_mean_stress_finite_full_v2 as reference
fit=reference.fit;memory=reference.memory;ROOT=reference.ROOT
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    source=ROOT/'docs/verification/mixed-mean-stress-finite-curved-result.json';parent=json.loads(source.read_text())
    coarse=[r for r in parent['records'] if r['n']==2];assert len(coarse)==4 and all(all(r['checks'].values()) for r in coarse)
    selected=[r for r in coarse if r['h']==.0001];assert len(selected)==2
    dest=ROOT/'docs/verification/mixed-mean-stress-curved-full-result.json';assert not dest.exists()
    start=time.monotonic();records=[];stops=[];out=ROOT/'build/analysis'/('mixed-mean-stress-curved-full-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    for r in selected:
        try:
            tick=time.monotonic();print('START_FULL',r['amplitude'],flush=True)
            a=fit.assemble(r['kappa'],r['n'],r['amplitude'],r['h']);sol,metrics=reference.solve_full(a);obs=fit.observe(a,sol)
            with np.load(ROOT/r['artifact']) as old:field=np.linalg.norm(sol['local']-old['local'],axis=(0,1))/np.linalg.norm(old['local'],axis=(0,1))
            error=max(abs(x-y) for x,y in zip(obs['error'],r['observed']['error']));position=max(abs(x-y) for x,y in zip(obs['recovered_position_error'],r['observed']['recovered_position_error']))
            checks=dict(original_stationarity=max(metrics['direct_kkt'])<=1e-8,linear=max(metrics['linear_residual'])<=1e-8,preconditioned_condition=metrics['preconditioned_full_rcond']>=1e-10,basis_map_backward_error=max(metrics['basis_map_backward_error'])<=1e-12,fields=max(field)<=1e-8,error=error<=1e-8,position=position<=1e-8,resources=time.monotonic()-tick<=180 and memory.memory_bytes()<=2*1024**3)
            path=out/f"a{r['amplitude']:g}.npz";np.savez_compressed(path,local=sol['local'],trace=sol['trace'])
            record=dict(kappa=r['kappa'],n=2,amplitude=r['amplitude'],h=r['h'],metrics=metrics,observed=obs,field_difference=field.tolist(),error_difference=error,position_difference=position,checks={k:bool(v) for k,v in checks.items()},seconds=time.monotonic()-tick,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path));records.append(record)
            print('CURVED_FULL',record,flush=True);del a,sol;gc.collect()
            if not all(checks.values()):raise RuntimeError('MMS1_curved_full_reference_gate')
            if time.monotonic()-start>900:raise RuntimeError('MMS1_curved_full_reference_series_resources')
        except Exception as exc:stops.append(dict(amplitude=r['amplitude'],h=r['h'],reason=repr(exc)));print('STOP',repr(exc),flush=True);break
    paths=[Path(__file__),Path(reference.__file__),Path(fit.__file__),source,ROOT/'docs/verification/mixed-mean-stress-curved-full-protocol.md']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_checks_pass=len(records)==2 and not stops,independent_global_resolve=True,n4_independent_global_resolve=False,main_shell_adopted=False,physical_PDE_acceptance=False,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),sources={p.relative_to(ROOT).as_posix():digest(p) for p in paths}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
