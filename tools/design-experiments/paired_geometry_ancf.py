"""Joint scalar position enrichment, with Cartesian rotation closure.

phi1=-k*C-z*I0 and phi3=b-k*x*z*I0 are scalar shape functions,
each multiplying three independent Cartesian coefficients. Their r/d pieces
are never independently enriched. Baseline r,d,c remain independent.
"""
import numpy as np
from scipy.linalg import svd
import geometry_compatible_ancf_stream as f
import paired_curvature_profiles as cp
base=f.c.base;m=f.m

def pool(q,p,n,metric):
    N,dN=m.shape(np.column_stack((p,np.zeros(len(p)))),0);X=N@q.reshape(75,3)[:25];J=np.einsum('pna,nc->pca',dN,q.reshape(75,3)[:25]);values=[];gx=[];gy=[]
    for axis in [0,1]:
        v,d=cp.profiles(X[:,axis],n,axis,0,metric);values.append(v);gx.append(d*J[None,:,axis,0,None]);gy.append(d*J[None,:,axis,1,None])
    return np.concatenate(values,axis=2),np.concatenate(gx,axis=2),np.concatenate(gy,axis=2)

def prepare(k,n,amp,order):
    geo=base.geometry_space(k,n,amp,order);geo.update(n=n,metric=[k,0.]);As=[];Ps=[]
    for cell,q in enumerate(geo['refs']):
        p,w=f.c.plane_rule(q,n,order);R=base.surface(n,cell,p,geo)[0][0];D=base.surface(n,cell,p,geo)[1][0];v,_,_=pool(q,p,n,geo['metric']);nr=R.shape[1];nd=D.shape[1];A=np.zeros((2*len(p),nr+nd));A[:len(p),:nr]=R;A[len(p):,nr:]=D
        weight=np.tile(np.sqrt(w),2);As.append(A*weight[:,None]);Ps.append(np.vstack(v)*weight[:,None])
    A=np.vstack(As);P=np.vstack(Ps);coef=np.linalg.lstsq(A,P,rcond=1e-13)[0];U,s,V=svd(P-A@coef,full_matrices=False,check_finite=False);ranks=[int(sum(s>tol*max(1.,np.linalg.norm(P)))) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError('paired_enrichment_rank:'+str(ranks))
    rank=ranks[0];geo['enrich']=[dict(coeff=coef,V=V[:rank].T,s=s[:rank],rank=rank,ranks=ranks)];return geo

def extra(n,cell,q,h,p,geo):
    fs=base.surface(n,cell,p[:,:2],geo);v,gx,gy=pool(q,p[:,:2],n,geo['metric']);e=geo['enrich'][0];nr=fs[0][0].shape[1];coef=e['coeff'];T=e['V']/e['s'];out=[]
    for deriv,raw in enumerate([v,gx,gy]):
        r=(raw[0]-fs[0][deriv]@coef[:nr])@T;d=(raw[1]-fs[1][deriv]@coef[nr:])@T;out.append((r,d))
    z=h*p[:,2]/2;N=out[0][0]+z[:,None]*out[0][1];G=np.stack((out[1][0]+z[:,None]*out[1][1],out[2][0]+z[:,None]*out[2][1],h/2*out[0][1]),axis=2)
    return N,G

def kinematics(n,cell,q,h,p,geo):
    N,G,_,det,X=base.kinematics(n,cell,q,h,p,geo);E,dE=extra(n,cell,q,h,p,geo);_,_,_,inv,_=f.f.base.geometry(q,h,p);G=np.concatenate((G,np.einsum('pni,pij->pnj',dE,inv,optimize=True)),axis=1);N=np.column_stack((N,E));B=np.einsum('kij,pnj->pkni',m.MB,G).reshape(len(p),6,3*N.shape[1]);return N,G,B,det,X

def shapes(a,cell,q,p):
    fs=base.surface(a['n'],cell,p[:,:2],a['geo']);z=a['h']*p[:,2]/2;N=np.column_stack((fs[0][0],z[:,None]*fs[1][0],z[:,None]**2*fs[2][0]/2));return np.column_stack((N,extra(a['n'],cell,q,a['h'],p,a['geo'])[0]))

def sizes(n,cell,q,p,geo):
    # Parent initializer only uses these shapes to determine column counts.
    fs=base.surface(n,cell,p,geo);rank=geo['enrich'][0]['rank'];fs[0]=tuple(np.column_stack((v,np.zeros((len(p),rank)))) for v in fs[0]);return fs

_initialize=f.initialize
def initialize(*args,**kwargs):
    a=_initialize(*args,**kwargs);a['paired_columns']=a['geo']['enrich'][0]['rank'];a['nr']-=a['paired_columns'];return a

f.c.prepare_geometry=prepare;f.c.fields=sizes;f.c.kinematics=kinematics;f.shape_values=shapes;f.initialize=initialize
