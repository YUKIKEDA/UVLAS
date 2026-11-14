"""Decisive raw-3D energy tests on the new general reference/unit-thickness space."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time,gc
import numpy as np
from scipy.linalg import svd,qr,solve_triangular
import coefficient_compatible_ancf as shell
import pg_ritz_targets as target
import pg_ritz_trace as trace
import pg_ritz_coarse_v2 as coarse
ROOT=rd_paths.ROOT
mu=1/2.6;lam=.3/(1.3*.4);C=np.diag([2*mu]*6);C[:3,:3]+=lam;CH=np.linalg.cholesky(C).T

def support(model):
    return np.eye(model['ns']),[0,0,0]

def solve(model,h,variable=False,order=15):
    plane0,weight0=shell.lift.geometry.quadrature(model['n'],13);zz,wzz=np.polynomial.legendre.leggauss(7);pp0=np.column_stack((np.repeat(plane0,7,axis=0),np.tile(zz,len(plane0))));ww0=(weight0[:,None]*wzz).ravel();Rb=None
    for start in range(0,len(pp0),1024):
        pp=pp0[start:start+1024];w=ww0[start:start+1024];Eg,Ug=coarse.fields(model,shell,pp,h,variable);E,U=target.fields(pp,h,model['parameters'],variable);det=shell.lift.geometry.reference(pp[:,:2],pp[:,2],h,model['parameters'],variable)[2];weight=np.sqrt(w*det)[:,None,None]/h**1.5;block=np.column_stack(((CH@Eg*weight).reshape(-1,8),(CH@E*weight).reshape(-1,3)));_,rr=qr(block if Rb is None else np.vstack((Rb,block)),mode='economic',check_finite=False);Rb=rr
    beta=solve_triangular(Rb[:8,:8],Rb[:8,8:],check_finite=False)
    Zs,ranks=support(model);Z=np.kron(Zs,np.eye(3));plane,wp=shell.lift.geometry.quadrature(model['n'],order);z,wz=np.polynomial.legendre.leggauss(7);p=np.column_stack((np.repeat(plane,7,axis=0),np.tile(z,len(plane))));weights=(wp[:,None]*wz).ravel();R=None;mass=np.zeros((model['ns'],model['ns']));reference=np.zeros(3)
    for start in range(0,len(p),1024):
        pp=p[start:start+1024];w=weights[start:start+1024];N,G,B,det,X=shell.kinematics(model,pp,h,variable);E,U=target.fields(pp,h,model['parameters'],variable);Eg,Ug=coarse.fields(model,shell,pp,h,variable);E=E-Eg@beta;Bz=B;weight=np.sqrt(w*det)[:,None,None]/h**1.5;A=(CH@Bz*weight).reshape(-1,Z.shape[1]);b=(CH@E*weight).reshape(-1,3);block=np.column_stack((A,b));_,rr=qr(block if R is None else np.vstack((R,block)),mode='economic',check_finite=False);R=rr;mass+=N.T@((w*det)[:,None]*N);reference+=np.einsum('p,pkt,kl,plt->t',w*det/h**3,E,C,E)
    n=Z.shape[1];scale=1/np.linalg.norm(R[:,:n],axis=0);Q,T=qr(R[:,:n]*scale,mode='economic',check_finite=False);sv=svd(T,compute_uv=False,check_finite=False);r=[int(sum(sv>tol*sv[0])) for tol in [1e-9,1e-10,1e-11]]
    if r!=[n]*3:raise ValueError('general_stiffness_rank:'+str(r)+'/'+str(n))
    y=solve_triangular(T,Q.T@R[:,n:],check_finite=False);coef=Z@(scale[:,None]*y);ms=1/np.sqrt(np.diag(mass));mmin=float(np.linalg.eigvalsh(mass*ms[:,None]*ms[None,:])[0]);return dict(coef=coef,beta=beta,Z=Z,Zs=Zs,scale=scale,rcond=float(sv[-1]/sv[0]),mass_min=mmin,support_ranks=ranks,reference=reference)

def observe(model,h,sol,variable=False,order=17):
    plane,wp=shell.lift.geometry.quadrature(model['n'],order);z,wz=np.polynomial.legendre.leggauss(9);p=np.column_stack((np.repeat(plane,9,axis=0),np.tile(z,len(plane))));weights=(wp[:,None]*wz).ravel();error=np.zeros(3);reference=np.zeros(3);pe=np.zeros(3);pn=np.zeros(3);gradient=np.zeros_like(sol['coef'])
    for start in range(0,len(p),1024):
        pp=p[start:start+1024];w=weights[start:start+1024];N,G,B,det,X=shell.kinematics(model,pp,h,variable);E,U=target.fields(pp,h,model['parameters'],variable);Eg,Ug=coarse.fields(model,shell,pp,h,variable);e=B@sol['coef']+(Eg@sol['beta']-E);u=(N@sol['coef'].reshape(model['ns'],9)).reshape(len(pp),3,3)+Ug@sol['beta'];weight=w*det/h**3;error+=np.einsum('p,pkt,kl,plt->t',weight,e,C,e);reference+=np.einsum('p,pkt,kl,plt->t',weight,E,C,E);pe+=np.einsum('p,pat,pat->t',w*det,u-U,u-U);pn+=np.einsum('p,pat,pat->t',w*det,U,U);gradient+=np.einsum('p,pkn,pkt->nt',weight,B,C@e,optimize=True)
    residual=np.linalg.norm(sol['scale'][:,None]*(sol['Z'].T@gradient),axis=0)/np.sqrt(sol['Z'].shape[1]*reference);return dict(error=np.sqrt(error/reference).tolist(),position=np.sqrt(pe/pn).tolist(),stationarity=residual.tolist())

