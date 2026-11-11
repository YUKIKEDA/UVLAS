"""Decisive raw-3D energy tests on the new general reference/unit-thickness space."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time,gc
import numpy as np
from scipy.linalg import svd,qr,solve_triangular
import pg_ritz_shell_v2 as shell
import pg_ritz_targets as target
ROOT=rd_paths.ROOT
mu=1/2.6;lam=.3/(1.3*.4);C=np.diag([2*mu]*6);C[:3,:3]+=lam;CH=np.linalg.cholesky(C).T

def support(model):
    p=np.column_stack((np.zeros(227),np.linspace(0,1,227)));components=shell.fields(model,p)[0];A=np.vstack(components);_,s,V=svd(A,full_matrices=True,check_finite=False);r=[int(sum(s>tol*max(1.,np.linalg.norm(A)))) for tol in [1e-9,1e-10,1e-11]]
    if len(set(r))!=1:raise ValueError('full_trace_support_rank:'+str(r))
    return V[r[0]:].T,r

def solve(model,h,variable=False,order=13):
    Zs,ranks=support(model);Z=np.kron(Zs,np.eye(3));plane,wp=shell.lift.geometry.quadrature(model['n'],order);z,wz=np.polynomial.legendre.leggauss(7);p=np.column_stack((np.repeat(plane,7,axis=0),np.tile(z,len(plane))));weights=(wp[:,None]*wz).ravel();R=None;mass=np.zeros((model['ns'],model['ns']));reference=np.zeros(3)
    for start in range(0,len(p),1024):
        pp=p[start:start+1024];w=weights[start:start+1024];N,G,B,det,X=shell.kinematics(model,pp,h,variable);E,U=target.fields(pp,h,model['parameters'],variable);Gz=np.einsum('pni,nj->pji',G,Zs,optimize=True);MB=np.zeros((6,3,3));MB[0,0,0]=MB[1,1,1]=MB[2,2,2]=1
        for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):MB[k,i,j]=MB[k,j,i]=1/np.sqrt(2)
        Bz=np.einsum('kij,pnj->pkni',MB,Gz).reshape(len(pp),6,3*Zs.shape[1]);weight=np.sqrt(w*det)[:,None,None]/h**1.5;A=(CH@Bz*weight).reshape(-1,Z.shape[1]);b=(CH@E*weight).reshape(-1,3);block=np.column_stack((A,b));_,rr=qr(block if R is None else np.vstack((R,block)),mode='economic',check_finite=False);R=rr;mass+=N.T@((w*det)[:,None]*N);reference+=np.einsum('p,pkt,kl,plt->t',w*det/h**3,E,C,E)
    n=Z.shape[1];scale=1/np.linalg.norm(R[:,:n],axis=0);Q,T=qr(R[:,:n]*scale,mode='economic',check_finite=False);sv=svd(T,compute_uv=False,check_finite=False);r=[int(sum(sv>tol*sv[0])) for tol in [1e-9,1e-10,1e-11]]
    if r!=[n]*3:raise ValueError('general_stiffness_rank:'+str(r)+'/'+str(n))
    y=solve_triangular(T,Q.T@R[:,n:],check_finite=False);coef=Z@(scale[:,None]*y);ms=1/np.sqrt(np.diag(mass));mmin=float(np.linalg.eigvalsh(mass*ms[:,None]*ms[None,:])[0]);return dict(coef=coef,Z=Z,Zs=Zs,scale=scale,rcond=float(sv[-1]/sv[0]),mass_min=mmin,support_ranks=ranks,reference=reference)

def observe(model,h,sol,variable=False,order=15):
    plane,wp=shell.lift.geometry.quadrature(model['n'],order);z,wz=np.polynomial.legendre.leggauss(9);p=np.column_stack((np.repeat(plane,9,axis=0),np.tile(z,len(plane))));weights=(wp[:,None]*wz).ravel();error=np.zeros(3);reference=np.zeros(3);pe=np.zeros(3);pn=np.zeros(3);gradient=np.zeros_like(sol['coef'])
    for start in range(0,len(p),1024):
        pp=p[start:start+1024];w=weights[start:start+1024];N,G,B,det,X=shell.kinematics(model,pp,h,variable);E,U=target.fields(pp,h,model['parameters'],variable);e=B@sol['coef']-E;u=(N@sol['coef'].reshape(model['ns'],9)).reshape(len(pp),3,3);weight=w*det/h**3;error+=np.einsum('p,pkt,kl,plt->t',weight,e,C,e);reference+=np.einsum('p,pkt,kl,plt->t',weight,E,C,E);pe+=np.einsum('p,pat,pat->t',w*det,u-U,u-U);pn+=np.einsum('p,pat,pat->t',w*det,U,U);gradient+=np.einsum('p,pkn,pkt->nt',weight,B,C@e,optimize=True)
    residual=np.linalg.norm(sol['scale'][:,None]*(sol['Z'].T@gradient),axis=0)/np.sqrt(sol['Z'].shape[1]*reference);return dict(error=np.sqrt(error/reference).tolist(),position=np.sqrt(pe/pn).tolist(),stationarity=residual.tolist())

def main():
    path=ROOT/'rd/structure/main-shell/paired-geometry/records/pg-ritz-finite-result.json';assert not path.exists();records=[];stops=[];models={};start=time.monotonic()
    for parameters,n,variable in [((.4,0.,0.),2,False),((.4,0.,0.),4,False),((.4,.2,0.),2,False),((.4,.2,0.),4,False),((.4,0.,.07),2,False),((.4,0.,.07),4,False),((.4,-.2,.07),2,True),((.4,-.2,.07),4,True)]:
        try:
            model=shell.prepare(n,parameters);models[parameters,n]=model;print('PREPARED',parameters,n,model['ns'],flush=True)
            for h in [.01,.0001]:
                tick=time.monotonic();sol=solve(model,h,variable);o=observe(model,h,sol,variable);record=dict(parameters=parameters,n=n,h=h,variable=variable,observation=o,rcond=sol['rcond'],mass_min=sol['mass_min'],seconds=time.monotonic()-tick);records.append(record);print(record,flush=True)
                if max(o['stationarity'])>1e-8:raise RuntimeError('original_stationarity')
        except Exception as exc:stops.append(dict(parameters=parameters,n=n,reason=repr(exc)));print('STOP',stops[-1],flush=True);break
    paths=[Path(__file__),Path(shell.__file__),Path(shell.lift.__file__),Path(shell.lift.particular.__file__),Path(shell.lift.geometry.__file__),Path(target.__file__),ROOT/'rd/structure/main-shell/paired-geometry/records/pg-ritz-general-protocol.md',ROOT/'rd/structure/main-shell/paired-geometry/records/pg-ritz-shell-v1-rank-stop.json'];result=dict(records=records,stops=stops,all_cases=len(records)==16 and not stops,main_shell_adopted=False,R4_complete=False,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});path.write_text(json.dumps(result,indent=2),encoding='utf-8')

if __name__=='__main__':main()
