"""Finite closed straight wake specialization, exclusively for reference trials."""
import numpy as np
from scipy.linalg import solve
from check_vortex_core_load import kernel

def geometry(m,n,dihedral=0.):
    if n%2:raise ValueError('even span count')
    xi=(1-np.cos(np.linspace(0,np.pi,m+1)))/2
    t=(1-np.cos(np.linspace(0,np.pi,n//2+1)))/2
    y=np.r_[4*(t-1),4*t[1:]]
    xv=np.r_[.75*xi[:-1]+.25*xi[1:],1.]
    xc=.25*xi[:-1]+.75*xi[1:]
    def points(x,y):
        x,y=np.broadcast_arrays(x,y);c=1.5-.9*abs(y)/4
        return np.stack([abs(y)*np.tan(np.deg2rad(35))+x*c,y,abs(y)*np.tan(np.deg2rad(dihedral))],axis=-1)
    X=points(xv[:,None],y[None,:]);P=points(xc[:,None],((y[:-1]+y[1:])/2)[None,:])
    left=X[:-1,:-1].reshape(-1,3);right=X[:-1,1:].reshape(-1,3)
    normal=np.cross(np.broadcast_to([1.,0,0],left.shape),right-left);normal/=np.linalg.norm(normal,axis=1)[:,None]
    return dict(m=m,n=n,X=X,P=P.reshape(-1,3),left=left,right=right,normal=normal,y=y)

def endpoints(g,length):
    X=g['X'].copy();X[-1,:,0]+=length
    a=np.stack([X[:-1,:-1],X[:-1,1:],X[1:,1:],X[1:,:-1]],axis=2).reshape(-1,4,3)
    # Last span leg is the finite downstream closure; adjacent x-directed legs have merged.
    return a,np.roll(a,-1,axis=1)

def matrix(g,eps,length,block=32):
    a,b=endpoints(g,length);p=g['P'];N=len(p);A=np.empty((N,N))
    for lo in range(0,N,block):
        hi=min(lo+block,N);K=kernel(p[lo:hi],a.reshape(-1,3),b.reshape(-1,3),eps).reshape(hi-lo,N,4,3).sum(axis=2)
        A[lo:hi]=np.einsum('pij,pj->pi',K,g['normal'][lo:hi])
    return A

def linear_load(g,horse,velocity=20.,rho=1.2):
    left,right=g['left'],g['right'];F=rho*horse[:,None]*np.cross([velocity,0,0],right-left)
    mid=(left+right)/2;M=np.cross(mid-[.375,0,0],F);bands=np.zeros(20)
    bounds=np.linspace(-4,4,21)
    for i,(a,b) in enumerate(zip(left,right)):
        for k,(lo,hi) in enumerate(zip(bounds[:-1],bounds[1:])):
            fraction=max(0.,min(b[1],hi)-max(a[1],lo))/(b[1]-a[1]);bands[k]+=fraction*F[i,2]
    return F.sum(axis=0),M.sum(axis=0),bands

def run(m,n,dihedral=0.,eps=1e-5,length=1024.):
    g=geometry(m,n,dihedral);A=matrix(g,eps,length);rhs=-20*g['normal'][:,2];gamma=solve(A,rhs,assume_a='gen',check_finite=True)
    residual=np.linalg.norm(A@gamma-rhs)/np.linalg.norm(rhs)
    G=gamma.reshape(m,n);horse=(G-np.vstack([np.zeros((1,n)),G[:-1]])).ravel()
    F,M,bands=linear_load(g,horse)
    return dict(g=g,A=A,gamma=gamma,horse=horse,F=F,M=M,bands=bands,residual=residual)
