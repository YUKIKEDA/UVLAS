"""Independent nonaffine 3D body/traction PDE on the same fixed PG space."""
from pathlib import Path
from types import SimpleNamespace
import hashlib,json,time,traceback
import numpy as np
from scipy.linalg import solve_triangular
import mass_paired_ancf as shell
import exact_paired_finite as engine
import general_paired_fields as coarse
ROOT=Path(__file__).resolve().parents[2]
mu=1/2.6;lam=.3/(1.3*.4);amplitude=.01
engine.shell=shell;engine.coarse=coarse.Coarse

def exact(p,h,parameters,variable=False):
    X=shell.geometry.reference(p[:,:2],p[:,2],h,parameters,variable)[0]
    U=np.zeros((len(p),3,3));E=np.zeros((len(p),6,3))
    U[:,0,:]=amplitude*X[:,0,None]**2;E[:,0,:]=2*amplitude*X[:,0,None]
    return E,U
engine.target=SimpleNamespace(fields=exact)

def cofactor(J,axis):
    return np.cross(J[:,:,(axis+1)%3],J[:,:,(axis+2)%3])

def external(model,order):
    h=model['h'];var=model['variable'];n=model['n'];force=np.zeros((model['ns'],3));resultant=np.zeros(3);work=0.
    body=np.array([-2*amplitude*(lam+2*mu),0.,0.])
    plane,wp=shell.geometry.quadrature(n,order);z,wz=np.polynomial.legendre.leggauss(6)
    for begin in range(0,len(plane),128):
        xy=plane[begin:begin+128];p=engine.points(xy,z);N,G,B,det,X=shell.kinematics(model,p,h,var)
        weight=(wp[begin:begin+128,None]*wz).ravel()*det
        force+=(N.T@weight)[:,None]*body;resultant+=sum(weight)*body
        work+=float(np.sum(weight*amplitude*X[:,0]**2*body[0]))
    q,w=np.polynomial.legendre.leggauss(order);edge=np.concatenate([(j+(q+1)/2)/(2*n) for j in range(2*n)]);ew=np.tile(w/(4*n),2*n)
    for axis in range(3):
        for sign in [-1,1]:
            if axis==2:p=np.column_stack((plane,np.full(len(plane),sign)));weight=wp
            else:
                p=np.zeros((len(edge)*len(z),3));p[:,axis]=(sign+1)/2;p[:,1-axis]=np.repeat(edge,len(z));p[:,2]=np.tile(z,len(edge));weight=(ew[:,None]*wz).ravel()
            for begin in range(0,len(p),512):
                pp=p[begin:begin+512];ww=weight[begin:begin+512];N,G,B,det,X=shell.kinematics(model,pp,h,var)
                J=shell.geometry.reference(pp[:,:2],pp[:,2],h,model['parameters'],var)[1];area=sign*cofactor(J,axis)
                diagonal=2*amplitude*X[:,0,None]*np.array([lam+2*mu,lam,lam]);traction=diagonal*area
                force+=N.T@(ww[:,None]*traction);resultant+=np.sum(ww[:,None]*traction,axis=0)
                work+=float(np.sum(ww*amplitude*X[:,0]**2*traction[:,0]))
    return force,dict(resultant=resultant.tolist(),exact_external_work=work)

def main():
    dest=ROOT/'docs/verification/mass-paired-body-result.json';assert not dest.exists();records=[];stops=[];start=time.monotonic()
    out=ROOT/'build/analysis'/('mass-paired-body-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    for n in [1,2]:
        for h in [.01,.0001]:
            try:
                model=shell.prepare(n,(.4,.2,0.,.1,-.06,.025),h,True)
                sol=engine.solve(model,h,True,order=model['degree']+3,beta_value=np.zeros((8,3)))
                force,loads=external(model,model['degree']+7)
                rhs=force.ravel()[:,None]/h**3
                y=solve_triangular(sol['root'].T,sol['scale'][:,None]*rhs,lower=True)
                coef=sol['scale'][:,None]*solve_triangular(sol['root'],y)
                physical=dict(sol,coef=np.repeat(coef,3,axis=1));obs=engine.observe(model,h,physical,True,order=model['degree']+9)
                observation=engine.observe(model,h,sol,True,order=model['degree']+9)
                difference=physical['coef']-sol['coef'];discrete=np.linalg.norm(sol['root']@(difference/sol['scale'][:,None]))/max(1e-30,np.linalg.norm(sol['root']@(sol['coef']/sol['scale'][:,None])))
                work=float(coef[:,0]@force.ravel());loadwork=abs(work-loads['exact_external_work'])/abs(loads['exact_external_work'])
                checks=dict(energy=max(obs['error'])<=.01,position=max(obs['position'])<=.01,
                    stationary=max(obs['stationarity'])<=1e-8,physical_vs_projection=discrete<=1e-8,
                    external_work=loadwork<=.01,force_balance=np.linalg.norm(loads['resultant'])<=1e-10)
                artifact=out/f'n{n}-h{h:g}.npz';np.savez_compressed(artifact,coef=coef,force=force,projection=sol['coef'])
                record=dict(n=n,h=h,observation=obs,projection=observation,physical_vs_projection=float(discrete),external_work_error=loadwork,
                    loads=loads,checks=checks,artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(artifact.read_bytes()).hexdigest())
                records.append(record);print('BODY',record,flush=True)
            except Exception as exc:stops.append(dict(n=n,h=h,reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',stops[-1],flush=True)
    paths=[Path(__file__),Path(shell.__file__),Path(engine.__file__),Path(coarse.__file__)]
    result=dict(passed=len(records)==4 and not stops and all(all(r['checks'].values()) for r in records),records=records,stops=stops,
        independent_load_definition='u_x=a X^2; sigma=diag(2aX(lambda+2mu),2aXlambda,2aXlambda); b=-div(sigma); t=sigma*n',
        seconds=time.monotonic()-start,R4_complete=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2,default=lambda value:value.item()),encoding='utf-8');print('BODY_RESULT',result['passed'],flush=True)

if __name__=='__main__':main()
