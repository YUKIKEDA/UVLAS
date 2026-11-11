"""Same R4 space on swept tapered footprint: quality and independent PDE."""
from pathlib import Path
from types import SimpleNamespace
import hashlib,json,time,traceback
import numpy as np
from scipy.linalg import solve_triangular
import mass_paired_ancf as shell
import general_paired_fields as target
import exact_paired_finite as engine
import check_mass_paired_body as body
import swept_taper_domain as domain
ROOT=Path(__file__).resolve().parents[2]
engine.shell=shell;engine.coarse=target.Coarse

def forces(model,order):
    n=model['n'];h=model['h'];ns=model['ns'];force=np.zeros((ns,3));weak=np.zeros_like(force);total=np.zeros(3);work=0.
    b=np.array([-2*body.amplitude*(body.lam+2*body.mu),0.,0.])
    plane,wp=domain.quadrature(n,order);z,wz=np.polynomial.legendre.leggauss(7)
    def stress(X):return 2*body.amplitude*X[:,0,None]*np.array([body.lam+2*body.mu,body.lam,body.lam])
    for begin in range(0,len(plane),128):
        xy=plane[begin:begin+128];p=engine.points(xy,z);N,G,B,det,X=shell.kinematics(model,p,h,True);ww=(wp[begin:begin+128,None]*wz).ravel()*det
        force+=(N.T@ww)[:,None]*b;total+=np.sum(ww)*b;work+=float(np.sum(ww*body.amplitude*X[:,0]**2*b[0]))
        weak+=np.einsum('p,pna,pa->na',ww,G,stress(X))
    q,w=np.polynomial.legendre.leggauss(order);edge=np.concatenate([(j+(q+1)/2)/(2*n) for j in range(2*n)]);ew=np.tile(w/(4*n),2*n)
    faces=[]
    for sign in [-1,1]:
        faces.append((np.column_stack((plane,np.full(len(plane),sign))),wp,np.array([0.,0.,sign])))
        s=(sign+1)/2;t=domain.lower(s)+(domain.upper(s)-domain.lower(s))*edge
        pp=np.column_stack((np.full(len(t)*len(z),s),np.repeat(t,len(z)),np.tile(z,len(t))))
        faces.append((pp,(ew[:,None]*wz*(domain.upper(s)-domain.lower(s))).ravel(),np.array([sign,0.,0.])))
    for function,covector in [(domain.lower,np.array([.3,-1.,0.])),(domain.upper,np.array([.1,1.,0.]))]:
        # Split the oblique face wherever it crosses a background cell.
        breaks={0.,1.};N=2*n
        for j in range(N+1):
            breaks.add(j/N)
            root=(j/N)/.3 if function is domain.lower else (1-j/N)/.1
            if 0<root<1:breaks.add(root)
        ss=[];ww=[];breaks=sorted(breaks)
        for l,r in zip(breaks[:-1],breaks[1:]):ss.extend((l+r)/2+(r-l)*q/2);ww.extend(w*(r-l)/2)
        ss=np.array(ss);ww=np.array(ww);p=np.column_stack((np.repeat(ss,len(z)),np.repeat(function(ss),len(z)),np.tile(z,len(ss))))
        faces.append((p,(ww[:,None]*wz).ravel(),covector))
    for points,weights,covector in faces:
        for begin in range(0,len(points),384):
            p=points[begin:begin+384];ww=weights[begin:begin+384];N,G,B,det,X=shell.kinematics(model,p,h,True)
            J=shell.geometry.reference(p[:,:2],p[:,2],h,model['parameters'],True)[1]
            area=sum(covector[a]*body.cofactor(J,a) for a in range(3));traction=stress(X)*area
            force+=N.T@(ww[:,None]*traction);total+=np.sum(ww[:,None]*traction,axis=0);work+=float(np.sum(ww*body.amplitude*X[:,0]**2*traction[:,0]))
    return force,weak,dict(resultant=total.tolist(),exact_external_work=work)

def solve_force(sol,F,h):
    y=solve_triangular(sol['root'].T,sol['scale'][:,None]*F.ravel()[:,None]/h**3,lower=True)
    return sol['scale'][:,None]*solve_triangular(sol['root'],y)

def main():
    dest=ROOT/'docs/verification/paired-shell-r5-domain-result.json';assert not dest.exists();records=[];failures=[];start=time.monotonic()
    out=ROOT/'build/analysis'/('paired-shell-r5-domain-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    old=shell.geometry.quadrature;shell.geometry.quadrature=domain.quadrature
    try:
        for n in [1,2]:
            p,w=domain.quadrature(n,23);moments=np.array([sum(w),w@p[:,0],w@p[:,1]])
            assert np.max(abs(moments-np.array([.8,11/30,131/300])))<1e-13
            for h in [.01,.0001]:
                try:
                    tick=time.monotonic();miss=shell.source_basis.cache_info().misses;model=shell.prepare(n,(.4,.2,0.,.1,-.06,.025),h,True)
                    extra=model['exact_source']['source_seconds'] if shell.source_basis.cache_info().misses==miss else 0.
                    engine.target=SimpleNamespace(fields=lambda p,h,parameters,variable:target.fields(p,h,parameters,variable,rates=(3.,5.,7.)))
                    sol=engine.solve(model,h,True,order=model['degree']+3);obs=engine.observe(model,h,sol,True,order=model['degree']+5)
                    quality_seconds=time.monotonic()-tick+extra
                    F,weak,loads=forces(model,model['degree']+7);coef=solve_force(sol,F,h);weakcoef=solve_force(sol,weak,h)
                    engine.target=SimpleNamespace(fields=body.exact);physical=dict(sol,coef=np.repeat(coef,3,axis=1),beta=np.zeros((8,3)))
                    physicalobs=engine.observe(model,h,physical,True,order=model['degree']+9)
                    delta=float(np.linalg.norm(sol['root']@((coef-weakcoef)/sol['scale'][:,None]))/np.linalg.norm(sol['root']@(weakcoef/sol['scale'][:,None])))
                    workerr=abs(float(coef[:,0]@F.ravel())-loads['exact_external_work'])/abs(loads['exact_external_work'])
                    checks=dict(original_stationarity=max(obs['stationarity'])<=1e-8,rcond=sol['rcond']>=1e-10,mass=sol['mass_min']>1e-10,resource=quality_seconds<=300.,
                        physical_energy=max(physicalobs['error'])<=.01,physical_position=max(physicalobs['position'])<=.01,physical_stationarity=max(physicalobs['stationarity'])<=1e-8,
                        independent_physical_solution=delta<=1e-8,work=workerr<=.01,force_balance=bool(np.linalg.norm(loads['resultant'])<=1e-10))
                    artifact=out/f'n{n}-h{h:g}.npz';np.savez_compressed(artifact,coef=sol['coef'],beta=sol['beta'],root=sol['root'],scale=sol['scale'],force=F,weak=weak,physical_coef=coef,
                        Rcoeff=model['Rcoeff'],Dcoeff=model['Dcoeff'],Ccoeff=model['Ccoeff'])
                    record=dict(n=n,h=h,ns=model['ns'],observation=obs,physical=physicalobs,physical_vs_weak=delta,work_error=workerr,loads=loads,
                        stiffness_ranks=sol['ranks'],mass_root_ranks=sol['mass_root_ranks'],rcond=sol['rcond'],mass_min=sol['mass_min'],checks=checks,
                        footprint_moments=moments.tolist(),cold_quality_seconds=quality_seconds,artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(artifact.read_bytes()).hexdigest())
                    records.append(record);print('DOMAIN',record,flush=True)
                except Exception as exc:failures.append(dict(n=n,h=h,reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    finally:shell.geometry.quadrature=old
    comparisons=[]
    if len(records)==4:
        for n in [1,2]:
            a,b=[r for r in records if r['n']==n]
            for key in ['error','position']:
                ratio=np.maximum(b['observation'][key],1e-12)/np.maximum(a['observation'][key],1e-12)
                comparisons.append(dict(kind='thin_amplification',n=n,quantity=key,values=ratio.tolist(),passed=bool(np.all(ratio<=2))))
        for h in [.01,.0001]:
            a,b=[r for r in records if r['h']==h]
            for key in ['error','position']:
                values=np.array(b['observation'][key]);comparisons.append(dict(kind='finest_quality',h=h,quantity=key,values=values.tolist(),
                    passed=bool(np.all(values<=.01) and np.all(values<=np.maximum(a['observation'][key],1e-12)))))
    paths=[Path(__file__),Path(domain.__file__),Path(shell.__file__),Path(target.__file__),Path(engine.__file__),Path(body.__file__),ROOT/'docs/verification/paired-shell-r5-domain-protocol.md']
    result=dict(passed=len(records)==4 and not failures and all(all(r['checks'].values()) for r in records) and all(c['passed'] for c in comparisons),
        records=records,comparisons=comparisons,failures=failures,seconds=time.monotonic()-start,
        sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('DOMAIN_RESULT',result['passed'],flush=True)

if __name__=='__main__':main()
