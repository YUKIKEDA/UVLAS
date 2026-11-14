"""Stable geometry-derived bending anchor plus a full-space correction."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import solve_triangular
import geometry_compatible_ancf_v4 as c
old=c.old;m=c.m;f=c.f
def bending(q,p,h,k,n):
    N,_=m.shape(p,0);x=(N@q.reshape(75,3)[:25])[:,0];z=h*p[:,2]/2
    bpp,I=c.curvature.evaluate(x,n,0,[k,0.]);i0,i1,i2=I;b=x[:,None]*i0-i1;bp=i0;primitive=(x*x)[:,None]*i0/2-i2/2
    a=np.column_stack((np.ones(len(x)),np.zeros(len(x)),k*x));tensor=-z[:,None,None,None]*a[:,:,None,None]*a[:,None,:,None]*bpp[:,None,None,:]/(1+k*k*x*x-k*z)[:,None,None,None]
    E=np.einsum('kij,pijn->pkn',m.MB,tensor);u=np.zeros((len(x),3,b.shape[1]));u[:,0]=-k*primitive-z[:,None]*bp;u[:,2]=b-k*x[:,None]*z[:,None]*bp
    return E,u
def anchor(a):
    rows=[]
    for q in a['geo']['refs']:
        p,w=c.volume_rule(q,a['n'],a['order'],a['zorder']);_,_,det,_,_=f.base.geometry(q,a['h'],p);E,_=bending(q,p,a['h'],a['k'],a['n']);target,_=old.reference(q,p,a['h'],a['k']);weight=np.sqrt(w*det)[:,None,None]/a['h']**1.5
        A=(np.einsum('ij,pjn->pin',m.CH,E)*weight).reshape(-1,E.shape[-1]);b=(np.einsum('ij,pjt->pit',m.CH,target)*weight).reshape(-1,3);rows.append(old.qrtools.rfactor(np.asfortranarray(np.column_stack((A,b)))))
    combined=np.vstack(rows);n=combined.shape[1]-3;scale=1/np.linalg.norm(combined[:,:n],axis=0);combined[:,:n]*=scale;R=old.qrtools.rfactor(np.asfortranarray(combined));value=scale[:,None]*solve_triangular(R[:n,:n],R[:n,n:],check_finite=False)
    return value
def residual_system(a,beta):
    rows=[];rhs=[]
    for cell,q in enumerate(a['geo']['refs']):
        p,w=c.volume_rule(q,a['n'],a['order'],a['zorder']);N,G,B,det,_=c.kinematics(a['n'],cell,q,a['h'],p,a['geo']);E,_=bending(q,p,a['h'],a['k'],a['n']);target,_=old.reference(q,p,a['h'],a['k']);target-=E@beta;weight=np.sqrt(w*det)[:,None,None]/a['h']**1.5
        A=(np.einsum('ij,pjn->pin',m.CH,B)*weight).reshape(-1,3*a['ns']);bb=(np.einsum('ij,pjt->pit',m.CH,target)*weight).reshape(-1,3);R=old.qrtools.rfactor(np.asfortranarray(np.column_stack((A,bb))));rows.append(R[:,:3*a['ns']]);rhs.append(R[:,3*a['ns']:])
    out=a.copy();out['A']=np.vstack(rows);out['b']=np.vstack(rhs);return out
def solve(a,free_rank=False):
    beta=anchor(a);residual=residual_system(a,beta);s=c.solve(residual,free_rank)
    if s['passed_rank']:s['beta']=beta
    return s
def observe(a,sol,order,zorder):
    error=np.zeros(3);norm=np.zeros(3);pe=np.zeros(3);pn=np.zeros(3);gradient=np.zeros((3*a['ns'],3));coef=sol['coef'];beta=sol['beta']
    for cell,q in enumerate(a['geo']['refs']):
        p,w=c.volume_rule(q,a['n'],order,zorder);N,G,B,det,_=c.kinematics(a['n'],cell,q,a['h'],p,a['geo']);E,U=bending(q,p,a['h'],a['k'],a['n']);target,pos=old.reference(q,p,a['h'],a['k']);e=B@coef+(E@beta-target);stress=np.einsum('ij,pjt->pit',m.C,e)
        error+=np.einsum('p,pkt,pkt->t',w*det,e,stress)/a['h']**3;norm+=np.einsum('p,pkt,kl,plt->t',w*det,target,m.C,target)/a['h']**3
        gradient+=np.einsum('p,pkn,pkt->nt',w*det,B,stress,optimize=True)/a['h']**3
        value=(N@coef.reshape(a['ns'],9)).reshape(len(p),3,3)+U@beta;pe+=np.einsum('p,pat,pat->t',w*det,value-pos,value-pos);pn+=np.einsum('p,pat,pat->t',w*det,pos,pos)
    scale=1/np.linalg.norm(a['A']@a['Z'],axis=0);g=scale[:,None]*(a['Z'].T@gradient);kkt=np.linalg.norm(g,axis=0)/np.maximum(1e-30,np.sqrt(len(scale))*np.sqrt(norm))
    return dict(error=np.sqrt(error/norm).tolist(),position=np.sqrt(pe/pn).tolist(),original_kkt=kkt.tolist())
