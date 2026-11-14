"""Finish the same four physical cases with a six-RHS shared geometry pass."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import gc,hashlib,json,time
import numpy as np
import pg_physical_pde_observe_batch as obs
import check_q9_plate_convergence as memory
pde=obs.pde;ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-physical-pde-batch-result.json';assert not dest.exists();start=time.monotonic();parent_path=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-physical-pde-result.json';parent=json.loads(parent_path.read_text());records=[dict(r,execution_route='original_passed') for r in parent['records'] if all(r['checks'].values())];assert len(records)==12;sources=[parent_path];prior={};stops=[]
    for name in ['paired-geometry-ancf-candidate-result.json','paired-geometry-ancf-spatial-result.json']:
        path=rd_paths.RECORDS/name;sources.append(path)
        for r in json.loads(path.read_text())['records']:prior[r.get('mode','rational'),r['amplitude'],r['n'],r['h']]=r
    out=ROOT/'build/analysis'/('pg-physical-pde-batch-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    def forbidden(*args,**kwargs):raise RuntimeError('manufactured_field_called_during_solve')
    for h in [.01,.0001]:
        try:
            n=4;amp=.2;pde.f.old.reference=forbidden;a=pde.prepare(.4,n,amp,h);tick=time.monotonic()-a['prepare_seconds'];sources.append(ROOT/a['preparation_manifest']);print('PREPARED',h,a['prepare_seconds'],flush=True)
            stage=time.monotonic();s=pde.solve(a,'combined');s['coarse_load']=pde.external(a,'combined')[2];print('SOLVED',h,time.monotonic()-stage,flush=True);previous={};geometry_equal=True
            for mode in ['rational','standard']:
                path=ROOT/prior[mode,amp,n,h]['artifact'];sources.append(path)
                with np.load(path) as d:
                    previous[mode]=dict(coef=d['coef'].copy(),beta=d['beta'].copy());geometry_equal=geometry_equal and all(np.array_equal(a['geo']['enrich'][0][key],d[f'enrich0_{key}']) for key in ['coeff','V','s'])
            previous={key:np.column_stack((previous['rational'][key],previous['standard'][key])) for key in ['coef','beta']};stage=time.monotonic();metrics,reaction=obs.observe(a,s,'combined',previous);print('OBSERVED',h,time.monotonic()-stage,flush=True);other=pde.external(a,'combined',13,18);elapsed=time.monotonic()-tick
            for j,mode in enumerate(['rational','standard']):
                sl=slice(3*j,3*j+3);part={key:(value[sl] if isinstance(value,list) else value) for key,value in metrics.items()};load_diff=float(np.linalg.norm(other[0][:,sl]-s['F'][:,sl])/np.linalg.norm(s['F'][:,sl]));checks=dict(geometry=geometry_equal,condition=a['rcond']>=1e-10,linear=s['linear']<=1e-8,equilibrium=max(part['stationarity'])<=1e-8,projection=max(part['projection_field_energy_difference']+part['projection_field_position_difference'])<=1e-8,load=load_diff<=1e-8,work=max(part['relative_work_defect'])<=1e-8,balance=max(part['relative_force_moment_balance'])<=1e-8,reaction=max(part['relative_reaction_resultant_difference'])<=1e-8,boundary=part['fixed_boundary']<=1e-9,resources=elapsed<=300 and memory.memory_bytes()<=2*1024**3)
                artifact=out/f'{mode}-a0.2-n4-h{h:g}.npz';np.savez_compressed(artifact,coef=s['coef'][:,sl],beta_global=s['beta'][:,sl],F=s['F'][:,sl],fixed_reference=s['fixed'][:,sl],reaction=reaction[:,sl],Zs=a['Zs'],coarse_load=s['coarse_load'][:,sl]);record=dict(execution_route='six_rhs_shared_geometry',mode=mode,amplitude=amp,n=n,h=h,metrics=part,load_difference=load_diff,rcond=a['rcond'],linear=s['linear'],checks={key:bool(value) for key,value in checks.items()},configuration_elapsed=elapsed,commutation_error=a['commutation_error'],load_value_error=a['load_value_error'],basis_cache_error=a['basis_cache_error'],plane_kernel_error=a['plane_kernel_error'],preparation_manifest=a['preparation_manifest'],artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=digest(artifact));records.append(record);print('CHECK',mode,h,checks,elapsed,flush=True)
            (out/'checkpoint.json').write_text(json.dumps(dict(records=records,stops=stops),indent=2),encoding='utf-8');del a,s;gc.collect()
            if not all(all(r['checks'].values()) for r in records):raise RuntimeError('physical_PDE_gate')
        except Exception as exc:stops.append(dict(h=h,reason=repr(exc)));print('STOP',stops[-1],flush=True);break
    accuracy=[]
    if len(records)==16 and not stops:
        lookup={(r['mode'],r['amplitude'],r['n'],r['h']):r for r in records}
        for mode in ['rational','standard']:
            for amp in [0.,.2]:
                for field in ['error','position']:
                    for target in range(3):
                        for n in [2,4]:
                            thick=lookup[mode,amp,n,.01]['metrics'][field][target];thin=lookup[mode,amp,n,.0001]['metrics'][field][target];floor=max(thick,thin)<1e-12;ratio=thin/max(thick,1e-300);accuracy.append(dict(mode=mode,amplitude=amp,n=n,field=field,target=target,thick=thick,thin=thin,numerical_floor=floor,ratio=ratio,passed=(floor or ratio<=2) and (n!=4 or max(thick,thin)<=.01)))
    paths=[Path(__file__),Path(obs.__file__),Path(obs.batch.__file__),ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-physical-pde-batch-protocol.md',ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-physical-pde-plane-result.json']+sources
    result=dict(passed=len(records)==16 and not stops and all(r['passed'] for r in accuracy) and time.monotonic()-start<=1200,records=records,accuracy=accuracy,stops=stops,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),manufactured_fields_forbidden_during_solve=True,main_shell_adopted=False,R4_complete=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in paths});dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('PASSED',result['passed'],flush=True)

if __name__=='__main__':main()
