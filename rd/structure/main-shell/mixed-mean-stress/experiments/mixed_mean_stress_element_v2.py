"""MMS-1 mixed stress elimination; conditional analytic-cylinder diagnostic."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.linalg import cholesky,solve_triangular
import moment_enriched_ancf as base
import mixed_mean_stress_space_v4 as space
m=base.m;c=base.c;ranks=base.ranks;NV=base.NV;PID=base.PID;VID=base.VID;VI=base.VI;basis=base.basis;reduce=base.reduce;boundary_map=base.boundary_map
def average(s,E):return np.einsum('z,pzk...->pk...',s['wz']/2,E.reshape(s['nplane'],s['zorder'],6,*E.shape[2:]),optimize=True)
def project(s,E):
    mean=average(s,E);rhs=np.einsum('p,pkr,pk...->r...',s['mean_weight'],s['stress_basis'],mean,optimize=True)
    flat=rhs.reshape(rhs.shape[0],-1);low=solve_triangular(s['stress_chol'],flat,lower=True,check_finite=False);lam=solve_triangular(s['stress_chol'].T,low,check_finite=False).reshape(rhs.shape)
    projected=np.einsum('pkr,r...->pk...',s['stress_basis'],lam,optimize=True);factor=s['mean_factor'].reshape(s['nplane'],s['zorder'],*([1]*(E.ndim-1)))
    tiled=np.broadcast_to(mean[:,None],(s['nplane'],s['zorder'],*mean.shape[1:]));high=(E.reshape(tiled.shape)-factor*tiled).reshape(E.shape)
    correction=(factor*projected[:,None]).reshape(E.shape);S=np.einsum('ij,pj...->pi...',m.C,high,optimize=True)+correction
    return S,high,low.reshape(rhs.shape)
def prepare(q,h,kappa,order=13,zorder=7,independent=False):
    s=base.prepare(q,h,order=order,zorder=zorder,independent=independent);plane,wp=m.quadrature(order,2);z,wz=np.polynomial.legendre.leggauss(zorder);L,meta=space.evaluate(q,kappa,plane,independent)
    det=s['det'].reshape(len(plane),zorder);gamma=(1/det)@(wz/2);mean_weight=2*wp/gamma
    gram=np.einsum('p,pir,ij,pjs->rs',mean_weight,L,np.linalg.inv(m.C),L,optimize=True);chol=cholesky(gram,lower=True,check_finite=False)
    s.update(kappa=kappa,nplane=len(plane),zorder=zorder,wz=wz,mean_weight=mean_weight,mean_factor=(1/(gamma[:,None]*det)).reshape(-1),stress_basis=L,stress_chol=chol,stress_space=meta,stress_gram_condition=float(np.linalg.cond(gram)))
    _,high,low=project(s,s['B']);elastic=np.vstack(((np.einsum('ij,pjn->pin',m.CH,high)*np.sqrt(s['weight'])[:,None,None]).reshape(-1,3*s['nscalar']),low))
    stab=s['root'][s['elastic_rows']:];s.update(root=np.vstack((elastic,stab)),elastic_rows=len(elastic),kind='MMS1')
    return s
def nonlinear_increment(s,increment,v1,v2):
    ns=s['nscalar'];u=np.asarray(increment).reshape(ns,3);Gw=s['Gw'];w=s['weight'];F=np.eye(3)+np.einsum('na,pjn->paj',u,Gw)
    E=np.einsum('kij,pij->pk',m.MB,(F.transpose(0,2,1)@F-np.eye(3))/2);S,high,low=project(s,E);St=np.einsum('kij,pk->pij',m.MB,S);P=F@St
    B=np.einsum('kij,pai,pjn->pkna',m.MB,F,Gw).reshape(-1,6,3*ns);e1=B@v1;e2=B@v2;S2,_,_=project(s,e2)
    d1=np.einsum('na,pjn->paj',v1.reshape(ns,3),Gw);d2=np.einsum('na,pjn->paj',v2.reshape(ns,3),Gw);sec=np.einsum('kij,pij->pk',m.MB,(d1.transpose(0,2,1)@d2+d2.transpose(0,2,1)@d1)/2)
    U=(np.einsum('p,pi,ij,pj->',w,high,m.C,high)+np.sum(low*low))/2;force=np.einsum('p,pjn,paj->na',w,Gw,P);K=np.einsum('p,pk,pk->',w,e1,S2)+np.einsum('p,pk,pk->',w,S,sec)
    physical=np.einsum('pjr,pjk->pkr',s['Sbody'],s['invJ'],optimize=True);bulk=np.einsum('p,pjr,paj->ar',w,physical,P,optimize=True);eta=np.zeros((3,s['H'].shape[0]),dtype=F.dtype);jumps=[]
    for face in s['faces']:
        jump=face['jump']@u;jumps.append(jump);ww=s['tau']*face['measure'];U+=np.einsum('p,pa,pa->',ww,jump,jump)/2;force+=np.einsum('p,pn,pa->na',ww,face['jump'],jump)
        K+=np.einsum('p,pa,pa->',ww,face['jump']@v1.reshape(ns,3),face['jump']@v2.reshape(ns,3));eta+=np.einsum('p,pr,pa->ar',ww,face['E'],jump)
    dual=(bulk+eta@s['H'])/8;fluxes=[face['sign']*np.einsum('pr,ar->pa',face['S'][:,face['axis'],:],dual)/face['area_norm'][:,None]-s['tau']*jump for face,jump in zip(s['faces'],jumps)]
    return U,force.reshape(-1),K,dict(min_F=float(np.min(np.linalg.det(F.real))),fluxes=fluxes)

def nonlinear(s,current,v1,v2):
    return nonlinear_increment(s,np.asarray(current).reshape(s['nscalar'],3)-s['X0'],v1,v2)
