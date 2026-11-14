"""Bounded qualification of GC-ANCF; preserve fields, references and failures."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import gc,hashlib,json,time
import numpy as np
from scipy.linalg import lstsq
import geometry_compatible_ancf_validation as v
import check_q9_plate_convergence as memory
f=v.f;ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(path,a,s):
    g=a['geo'];items=dict(coef=s['coef'],beta=s['beta'],Zs=a['Zs'],refs=np.array(g['refs']),geo_coeff=g['coeff'],geo_V=g['V'],geo_s=g['s'],metric=np.array(g['metric']),sizes=np.array([a['nr'],a['nd'],a['nc']]))
    for j,e in enumerate(g['enrich']):
        for key in ['coeff','V','s']:items[f'enrich{j}_{key}']=e[key]
    np.savez_compressed(path,**items)
def main():
    dest=ROOT/'rd/structure/main-shell/geometry-compatible/records/geometry-compatible-ancf-candidate-result.json';assert not dest.exists();start=time.monotonic();records=[];stops=[];geos={};out=ROOT/'build/analysis'/('geometry-compatible-ancf-qualified-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    cases=[('standard',amp,n,h) for amp in [0.,.2] for n in [2,4] for h in [.01,.0001]]+[('rational',amp,4,h) for amp in [0.,.2] for h in [.01,.0001]]
    for mode,amp,n,h in cases:
        try:
            tick=time.monotonic();v.select(mode);print('START',mode,amp,n,h,flush=True);a=f.initialize(.4,n,amp,h,13,7,geos.get((n,amp)));geos[n,amp]=a['geo'];print('SOLVE',mode,amp,n,h,a['ns'],flush=True)
            free=mode=='standard' and h==.01;s=f.solve(a,free_rank=free)
            if not s['passed_rank']:raise ValueError('GC_fixed_rank:'+str(s))
            primary=f.observe(a,s,13,7);ind=f.observe(a,s,15,9);traces=v.interface(a,s);anchor=v.anchor_consistency(a,s,13,7)
            mechanics=v.mechanics(a,13,7) if mode=='standard' and n==2 and h==.01 else None;svd_record=None
            if mode=='standard' and n==2 and h==.0001:
                scale=1/np.linalg.norm(a['A']@a['Z'],axis=0);x,_,rank,sv=lstsq((a['A']@a['Z'])*scale,a['b'],cond=1e-12,lapack_driver='gelsd',check_finite=False);other=dict(coef=a['Z']@(scale[:,None]*x),beta=s['beta']);oi=f.observe(a,other,15,9)
                svd_record=dict(rank=int(rank),columns=len(scale),error_difference=float(max(abs(x-y) for x,y in zip(oi['error'],ind['error']))),position_error_difference=float(max(abs(x-y) for x,y in zip(oi['position'],ind['position']))),original_kkt=oi['original_kkt'])
            checks=dict(fixed_rank=s['passed_rank'],condition=s['rcond']>=1e-10,linear=max(s['linear'])<=1e-8,original_stationarity=max(primary['original_kkt'])<=1e-8,compression=max(abs(x-y) for x,y in zip(primary['error'],s['error']))<=1e-8,independent=max(abs(x-y) for x,y in zip(primary['error'],ind['error']))<=1e-8,mass=a['mass_scaled_min']>1e-10,trace=max(traces.values())<=1e-9,anchor=max(anchor)<=1e-8,free_rigid=(not free or s['free_nullities']==[6]*3),mechanics=(mechanics is None or max(mechanics.values())<=1e-9),svd=(svd_record is None or (svd_record['rank']==svd_record['columns'] and svd_record['error_difference']<=1e-8 and svd_record['position_error_difference']<=1e-8 and max(svd_record['original_kkt'])<=1e-8)),resources=time.monotonic()-tick<=300 and memory.memory_bytes()<=2*1024**3)
            path=out/f'{mode}-a{amp:g}-n{n}-h{h:g}.npz';save(path,a,s);record=dict(mode=mode,amplitude=amp,n=n,h=h,fit={key:value for key,value in s.items() if key not in ['coef','beta']},primary=primary,independent=ind,traces=traces,anchor_consistency=anchor,mechanics=mechanics,svd=svd_record,mass_scaled_min=a['mass_scaled_min'],anchor_ranks=a['anchor_ranks'],anchor_columns=a['anchor_columns'],boundary_ranks=a['boundary_ranks'],sizes=[a['nr'],a['nd'],a['nc']],checks={k:bool(x) for k,x in checks.items()},seconds=time.monotonic()-tick,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path));records.append(record)
            (out/'checkpoint.json').write_text(json.dumps(dict(records=records,stops=stops),indent=2),encoding='utf-8');print('QUALIFY',mode,amp,n,h,ind,checks,record['seconds'],flush=True);del a,s;gc.collect()
            if not all(checks.values()):raise RuntimeError('GC_candidate_numerical_gate')
            if time.monotonic()-start>1800:raise RuntimeError('GC_candidate_series_resources')
        except Exception as exc:stops.append(dict(mode=mode,amplitude=amp,n=n,h=h,reason=repr(exc)));print('STOP',repr(exc),flush=True);break
    accuracy=[];complete=len(records)==len(cases) and not stops
    if complete:
        lookup={(r['mode'],r['amplitude'],r['n'],r['h']):r for r in records}
        for mode in ['standard','rational']:
            for amp in [0.,.2]:
                for field in ['error','position']:
                    for target in ([1,2] if mode=='standard' else [0,1,2]):
                        coarse=lookup[(mode,amp,4,.01)]['independent'][field][target];thin=lookup[(mode,amp,4,.0001)]['independent'][field][target];ratio=thin/max(coarse,1e-300)
                        convergence=True if mode=='rational' else all(lookup[(mode,amp,4,h)]['independent'][field][target]<lookup[(mode,amp,2,h)]['independent'][field][target] for h in [.01,.0001])
                        accuracy.append(dict(mode=mode,amplitude=amp,field=field,target=target,coarse_thickness_error=coarse,thin_error=thin,ratio=ratio,finest=max(coarse,thin)<=.01,uniform=ratio<=2,convergence=convergence,numerical_floor=max(coarse,thin)<1e-12))
    qualified=bool(complete and all(r['finest'] and r['uniform'] and (r['convergence'] or r['numerical_floor']) for r in accuracy))
    paths=[Path(__file__),Path(v.__file__),Path(f.__file__),Path(f.c.__file__),ROOT/'rd/structure/main-shell/reference-work/experiments/stable_curvature_profiles.py',ROOT/'rd/structure/main-shell/reference-work/experiments/curvature_work_profiles.py',ROOT/'rd/structure/main-shell/geometry-compatible/records/geometry-compatible-ancf-qualification-protocol.md']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_numerical_checks=complete,accuracy=accuracy,adoption_candidate_qualified=qualified,main_shell_adopted=False,physical_PDE_acceptance=False,human_independent_review=False,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),sources={p.relative_to(ROOT).as_posix():digest(p) for p in paths}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
