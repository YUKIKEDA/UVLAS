"""R4.1 whole-field, reaction and work checks for independent physical forcing."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import gc,hashlib,json,time
import numpy as np
import pg_physical_pde_observe as obs
import check_q9_plate_convergence as memory
pde=obs.pde;ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-physical-pde-result.json';assert not dest.exists();start=time.monotonic();records=[];stops=[];geos={};sources=[];prior={}
    for name in ['paired-geometry-ancf-candidate-result.json','paired-geometry-ancf-spatial-result.json']:
        path=rd_paths.RECORDS/name;sources.append(path)
        for r in json.loads(path.read_text())['records']:prior[r.get('mode','rational'),r['amplitude'],r['n'],r['h']]=r
    out=ROOT/'build/analysis'/('pg-physical-pde-full-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    def forbidden(*args,**kwargs):raise RuntimeError('manufactured_field_called_during_solve')
    for amp,n,h in [(amp,n,h) for amp in [0.,.2] for n in [2,4] for h in [.01,.0001]]:
        tick=time.monotonic()
        try:
            pde.f.old.reference=forbidden;a=pde.prepare(.4,n,amp,h);print('PREPARED',amp,n,h,a['rcond'],flush=True)
            for mode in ['rational','standard']:
                pde.f.old.reference=forbidden;s=pde.solve(a,mode);s['coarse_load']=pde.external(a,mode)[2]
                old=prior[mode,amp,n,h];oldpath=ROOT/old['artifact'];sources.append(oldpath)
                with np.load(oldpath) as data:
                    previous=dict(coef=data['coef'].copy(),beta=data['beta'].copy());geometry_equal=all(np.array_equal(a['geo']['enrich'][0][key],data[f'enrich0_{key}']) for key in ['coeff','V','s'])
                metrics,reaction=obs.observe(a,s,mode,previous);other=pde.external(a,mode,13,18);load_diff=float(np.linalg.norm(other[0]-s['F'])/np.linalg.norm(s['F']));checks=dict(geometry=geometry_equal,condition=a['rcond']>=1e-10,linear=s['linear']<=1e-8,equilibrium=max(metrics['stationarity'])<=1e-8,projection=max(metrics['projection_field_energy_difference']+metrics['projection_field_position_difference'])<=1e-8,load=load_diff<=1e-8,work=max(metrics['relative_work_defect'])<=1e-8,balance=max(metrics['relative_force_moment_balance'])<=1e-8,reaction=max(metrics['relative_reaction_resultant_difference'])<=1e-8,boundary=metrics['fixed_boundary']<=1e-9,resources=time.monotonic()-tick<=300 and memory.memory_bytes()<=2*1024**3)
                artifact=out/f'{mode}-a{amp:g}-n{n}-h{h:g}.npz';np.savez_compressed(artifact,coef=s['coef'],beta_global=s['beta'],F=s['F'],fixed_reference=s['fixed'],reaction=reaction,Zs=a['Zs'],coarse_load=s['coarse_load']);record=dict(mode=mode,amplitude=amp,n=n,h=h,metrics=metrics,load_difference=load_diff,rcond=a['rcond'],linear=s['linear'],checks={k:bool(v) for k,v in checks.items()},configuration_elapsed=time.monotonic()-tick,artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=digest(artifact));records.append(record);(out/'checkpoint.json').write_text(json.dumps(dict(records=records,stops=stops),indent=2),encoding='utf-8');print('CHECK',mode,amp,n,h,checks,metrics['error'],metrics['projection_field_position_difference'],flush=True);del s;gc.collect()
                if not all(checks.values()):raise RuntimeError('physical_PDE_gate')
            del a;gc.collect()
            if time.monotonic()-start>1200:raise RuntimeError('physical_PDE_series_resources')
        except Exception as exc:stops.append(dict(amplitude=amp,n=n,h=h,reason=repr(exc)));print('STOP',stops[-1],flush=True);break
    accuracy=[]
    if len(records)==16 and not stops:
        lookup={(r['mode'],r['amplitude'],r['n'],r['h']):r for r in records}
        for mode in ['rational','standard']:
            for amp in [0.,.2]:
                for field in ['error','position']:
                    for target in range(3):
                        for n in [2,4]:
                            thick=lookup[mode,amp,n,.01]['metrics'][field][target];thin=lookup[mode,amp,n,.0001]['metrics'][field][target];floor=max(thick,thin)<1e-12;ratio=thin/max(thick,1e-300);accuracy.append(dict(mode=mode,amplitude=amp,n=n,field=field,target=target,thick=thick,thin=thin,numerical_floor=floor,ratio=ratio,passed=(floor or ratio<=2) and (n!=4 or max(thick,thin)<=.01)))
    paths=[Path(__file__),Path(obs.__file__),Path(pde.__file__),Path(pde.loads.__file__),ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-physical-pde-full-protocol.md',ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-physical-pde-screen.json',ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-load-moment-audit.json',ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-candidate-audit.json']+sources
    result=dict(passed=len(records)==16 and not stops and all(a['passed'] for a in accuracy),records=records,accuracy=accuracy,stops=stops,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),manufactured_fields_forbidden_during_solve=True,main_shell_adopted=False,R4_complete=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in paths});dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('PASSED',result['passed'],flush=True)

if __name__=='__main__':main()
