"""Finite-thickness quality of the exact graph PG space; no thresholds fitted."""
from pathlib import Path
from types import SimpleNamespace
import gc,hashlib,json,time,traceback
import numpy as np
import mass_paired_ancf as shell
import general_paired_fields as target
import exact_paired_finite as engine
import check_exact_paired_primary as prior
ROOT=Path(__file__).resolve().parents[2]
engine.shell=shell;engine.coarse=target.Coarse;prior.shell=shell

def inclusion(model):
    p=np.random.default_rng(20261105).uniform(.001,.999,(91,2));F=shell.fields(model,p)
    v,d,ds=target.jets(p,model['parameters'],'polynomial');Q=model['coarse_coordinates'];out=[]
    for a in range(3):
        r=(F[a][0]@Q).reshape(len(p),3,8);dr=(F[a][1]@Q).reshape(len(p),3,8)
        ref=d if a==0 else ds[a-1]
        out.append(dict(r=float(np.linalg.norm(r-v[a])/max(1.,np.linalg.norm(v[a]))),
            d=float(np.linalg.norm(dr-ref)/max(1.,np.linalg.norm(ref)))))
    return out

def main():
    dest=ROOT/'docs/verification/mass-paired-quality-result.json';assert not dest.exists()
    out=ROOT/'build/analysis'/('mass-paired-quality-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    jobs=[('primary',p,n,v,(2.,4.,6.)) for p,v in [((.4,0.,0.),False),((.4,.2,0.),False),((.4,0.,.07),False),((.4,-.2,.07),True)] for n in [2,4]]
    jobs += [('nonquadratic',p,n,True,(3.,5.,7.)) for p in [(.4,.2,.07,.1,-.06,.025),(.4,.2,.07,.1,-.06,.025,0.,1.,.2)] for n in [1,2]]
    records=[];failures=[];start=time.monotonic()
    for family,pars,n,var,rates in jobs:
        engine.target=SimpleNamespace(fields=lambda p,h,parameters,variable:target.fields(p,h,parameters,variable,rates=rates))
        for h in [.01,.0001]:
            try:
                tick=time.monotonic();miss=shell.source_basis.cache_info().misses;model=shell.prepare(n,pars,h,var)
                prep=time.monotonic()-tick;reuse=shell.source_basis.cache_info().misses==miss
                cold=prep+(model['exact_source']['source_seconds'] if reuse else 0.)
                inc=inclusion(model);equiv=prior.local_equivalence(model,h,var)
                sol=engine.solve(model,h,var,order=model['degree']+3)
                obs=engine.observe(model,h,sol,var,order=model['degree']+5)
                face=np.array([[0.,y,z] for y in np.linspace(0,1,71) for z in [-1.,-.31,0.,.47,1.]])
                N,*_=shell.kinematics(model,face,h,var)
                trace=float(np.max(abs(N)));seconds=time.monotonic()-tick+(model['exact_source']['source_seconds'] if reuse else 0.)
                checks=dict(inclusion=max(max(a.values()) for a in inc)<=1e-8,equivalence=max(equiv.values())<=1e-12,
                    stationarity=max(obs['stationarity'])<=1e-8,rcond=sol['rcond']>=1e-10,mass=sol['mass_min']>1e-10,
                    trace=trace<=1e-9,resource=seconds<=300.)
                artifact=out/f'case-{len(records):02d}.npz'
                np.savez_compressed(artifact,**{k:sol[k] for k in ['coef','beta','scale','root','massroot']},
                    **{k:model[k] for k in ['Rcoeff','Dcoeff','Ccoeff','coarse_coordinates']},
                    parameters=np.array(pars),h=h,n=n,variable=var,degrees=np.array([model['r_degree'],model['degree']]))
                r=dict(family=family,parameters=pars,n=n,h=h,variable=var,rates=rates,ns=model['ns'],observation=obs,
                    inclusion=inc,equivalence=equiv,trace=trace,rcond=sol['rcond'],mass_min=sol['mass_min'],
                    stiffness_ranks=sol['ranks'],mass_root_ranks=sol['mass_root_ranks'],checks=checks,
                    preparation=prep,cold_preparation=cold,cold_case_seconds=seconds,exact_source=model['exact_source'],normalization=model['normalization'],
                    artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(artifact.read_bytes()).hexdigest())
                records.append(r);print('CASE',len(records),family,pars,n,h,obs,checks,seconds,flush=True)
                del sol,model;gc.collect()
            except Exception as exc:
                failures.append(dict(family=family,parameters=pars,n=n,h=h,reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
            (out/'checkpoint.json').write_text(json.dumps(dict(records=records,failures=failures),indent=2),encoding='utf-8')
    comparisons=[]
    for family,pars,n,var,rates in jobs:
        pair=[next((r for r in records if r['parameters']==pars and r['n']==n and r['h']==h),None) for h in [.01,.0001]]
        if any(r is None for r in pair):continue
        for quantity in ['error','position']:
            ratio=np.maximum(pair[1]['observation'][quantity],1e-12)/np.maximum(pair[0]['observation'][quantity],1e-12)
            comparisons.append(dict(kind='thin_amplification',parameters=pars,n=n,quantity=quantity,values=ratio.tolist(),passed=bool(np.all(ratio<=2))))
        fine=4 if family=='primary' else 2
        if n!=fine:continue
        for r in pair:
            lower=next((s for s in records if s['parameters']==pars and s['n']==fine//2 and s['h']==r['h']),None)
            for quantity in ['error','position']:
                v=np.array(r['observation'][quantity]);checks=dict(accuracy=bool(np.all(v<=.01)))
                if lower:checks['decrease']=bool(np.all(v<=np.maximum(lower['observation'][quantity],1e-12)))
                comparisons.append(dict(kind='finest_quality',parameters=pars,h=r['h'],quantity=quantity,values=v.tolist(),**checks,passed=all(checks.values())))
    paths=[Path(__file__),ROOT/'docs/paired-reference-coefficient-definition.md']
    paths += [ROOT/'tools/design-experiments'/s for s in ['mass_paired_ancf.py','general_paired_polynomials.py','general_paired_ancf.py','polynomial_reference_surface.py','general_paired_fields.py','exact_paired_finite.py','paired_exact_polynomials.py','coefficient_compatible_ancf.py','exact_paired_ancf.py','check_exact_paired_primary.py']]
    result=dict(passed=len(records)==24 and not failures and all(all(r['checks'].values()) for r in records) and all(c['passed'] for c in comparisons),
        records=records,failures=failures,comparisons=comparisons,seconds=time.monotonic()-start,R4_complete=False,
        sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('QUALITY_RESULT',result['passed'],flush=True)

if __name__=='__main__':main()
