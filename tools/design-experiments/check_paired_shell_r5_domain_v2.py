"""Reuse immutable space/roots; split tip traction at actual background knots."""
from pathlib import Path
from types import SimpleNamespace
import hashlib,json,time
import numpy as np
import check_paired_shell_r5_domain as parent
ROOT=Path(__file__).resolve().parents[2]
shell=parent.shell;engine=parent.engine;body=parent.body;domain=parent.domain

def forces(model,order):
    # Parent integrates each s-face after rescaling uniform edge subintervals;
    # these are not the native t knots when the face length changes. Correct
    # only those face contributions, keeping body/oblique/top/bottom unchanged.
    F,weak,loads=parent.forces(model,order)
    n=model['n'];h=model['h'];q,w=np.polynomial.legendre.leggauss(order);z,wz=np.polynomial.legendre.leggauss(7)
    oldedge=np.concatenate([(j+(q+1)/2)/(2*n) for j in range(2*n)]);oldweight=np.tile(w/(4*n),2*n)
    for sign in [-1,1]:
        s=(sign+1)/2;lo=domain.lower(s);hi=domain.upper(s)
        knots=sorted({lo,hi}|{j/(2*n) for j in range(2*n+1) if lo<j/(2*n)<hi})
        tt=np.concatenate([(a+b)/2+(b-a)*q/2 for a,b in zip(knots[:-1],knots[1:])])
        tw=np.concatenate([w*(b-a)/2 for a,b in zip(knots[:-1],knots[1:])])
        for factor,t,weight in [(-1.,lo+(hi-lo)*oldedge,oldweight*(hi-lo)),(1.,tt,tw)]:
            p=np.column_stack((np.full(len(t)*len(z),s),np.repeat(t,len(z)),np.tile(z,len(t))))
            ww=(weight[:,None]*wz).ravel();N,G,B,det,X=shell.kinematics(model,p,h,True)
            J=shell.geometry.reference(p[:,:2],p[:,2],h,model['parameters'],True)[1]
            area=sign*body.cofactor(J,0);sigma=2*body.amplitude*X[:,0,None]*np.array([body.lam+2*body.mu,body.lam,body.lam]);traction=sigma*area
            F+=factor*(N.T@(ww[:,None]*traction))
    return F,weak,loads

def main():
    dest=ROOT/'docs/verification/paired-shell-r5-domain-v2-result.json';assert not dest.exists();start=time.monotonic()
    source=ROOT/'docs/verification/paired-shell-r5-domain-result.json';old=json.loads(source.read_text());assert len(old['records'])==4
    out=ROOT/'build/analysis'/('paired-shell-r5-domain-v2-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    records=[];paths=[Path(__file__),Path(parent.__file__),source,Path(domain.__file__),Path(shell.__file__),Path(body.__file__),Path(engine.__file__)]
    original=shell.geometry.quadrature;shell.geometry.quadrature=domain.quadrature;engine.target=SimpleNamespace(fields=body.exact)
    try:
        for r in old['records']:
            path=ROOT/r['artifact'];assert hashlib.sha256(path.read_bytes()).hexdigest()==r['artifact_sha256'];paths.append(path);a=np.load(path)
            model=dict(n=r['n'],h=r['h'],ns=r['ns'],variable=True,boundary='clamped',parameters=(.4,.2,0.,.1,-.06,.025),r_degree=12,degree=21,
                Rcoeff=a['Rcoeff'],Dcoeff=a['Dcoeff'],Ccoeff=a['Ccoeff'])
            sol=dict(root=a['root'],scale=a['scale']);F,weak,loads=forces(model,28);coef=parent.solve_force(sol,F,r['h']);cw=parent.solve_force(sol,weak,r['h'])
            physical=dict(sol,coef=np.repeat(coef,3,axis=1),beta=np.zeros((8,3)))
            obs=engine.observe(model,r['h'],physical,True,order=30)
            delta=float(np.linalg.norm(sol['root']@((coef-cw)/sol['scale'][:,None]))/np.linalg.norm(sol['root']@(cw/sol['scale'][:,None])))
            workerr=abs(float(coef[:,0]@F.ravel())-loads['exact_external_work'])/abs(loads['exact_external_work'])
            checks=dict(r['checks']);checks.update(physical_energy=max(obs['error'])<=.01,physical_position=max(obs['position'])<=.01,
                physical_stationarity=max(obs['stationarity'])<=1e-8,independent_physical_solution=delta<=1e-8,work=workerr<=.01)
            artifact=out/f"n{r['n']}-h{r['h']:g}.npz";np.savez_compressed(artifact,force=F,weak=weak,physical_coef=coef,weak_coef=cw)
            record=dict(r,physical=obs,physical_vs_weak=delta,work_error=workerr,checks=checks,
                parent_artifact=r['artifact'],parent_artifact_sha256=r['artifact_sha256'],artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(artifact.read_bytes()).hexdigest())
            records.append(record);print('CORRECTED_DOMAIN',record,flush=True)
    finally:shell.geometry.quadrature=original
    result=dict(passed=all(all(r['checks'].values()) for r in records) and all(c['passed'] for c in old['comparisons']),
        records=records,comparisons=old['comparisons'],seconds=time.monotonic()-start,
        correction='s=1 face split at actual native t knots; unchanged N, space, material, roots, and thresholds',
        sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('DOMAIN_V2_RESULT',result['passed'],flush=True)

if __name__=='__main__':main()
