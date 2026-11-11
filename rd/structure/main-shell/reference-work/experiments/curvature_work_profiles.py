"""Fixed reference-metric curvature space and its first two primitives."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
def evaluate(x,n,axis,metric):
    lo=0. if axis==0 else -.5;u=np.clip(np.asarray(x)-lo,0.,1.);unique,back=np.unique(u,return_inverse=True);nb=16*n
    vals=np.zeros((len(unique),nb));I=np.zeros((3,len(unique),nb));g,w=np.polynomial.legendre.leggauss(32);curv=metric[axis]
    for cell in range(2*n):
        left=cell/(2*n);right=(cell+1)/(2*n);upper=np.clip(unique,left,right);length=upper-left;quad=left+length[:,None]*(g+1)/2;local=2*(quad-left)/(right-left)-1
        pol=np.polynomial.legendre.legvander(local,7)/(1+curv*curv*quad*quad)[:,:,None]**2
        for power in range(3):I[power,:,8*cell:8*(cell+1)]=np.einsum('pq,pqk,q,p->pk',quad**power,pol,w,length/2)
        inside=(unique>=left)&((unique<right) if cell<2*n-1 else (unique<=right))
        xx=unique[inside];vals[inside,8*cell:8*(cell+1)]=np.polynomial.legendre.legvander(2*(xx-left)/(right-left)-1,7)/(1+curv*curv*xx*xx)[:,None]**2
    return vals[back],I[:,back]
def profiles(x,n,axis,kind,metric):
    u=np.asarray(x)-(0. if axis==0 else -.5);k,I=evaluate(x,n,axis,metric);i0,i1,i2=I
    if kind==0:return np.column_stack((u[:,None]*i0-i1,(u*u)[:,None]*i0/2-i2/2)),np.column_stack((i0,u[:,None]*i0))
    return np.column_stack((i0,u[:,None]*i0)),np.column_stack((k,i0+u[:,None]*k))
