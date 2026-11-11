"""Steady ring specialization with exact material self-line exclusion before summing."""
import numpy as np
from numba import njit,prange
from uvlm_frequency_fast import segment
@njit(parallel=True)
def _matrix(P,N,X,eps,own):
    m=X.shape[0]-1;n=X.shape[1]-1;C=np.empty((len(P),m*n))
    for k in prange(len(P)):
        for i in range(m):
            for j in range(n):
                col=i*n+j;a=X[i,j];b=X[i,j+1];c=X[i+1,j+1];d=X[i+1,j];v=segment(P[k],b,c,eps,N[k])[0]+segment(P[k],d,a,eps,N[k])[0]
                if col!=own[k]:v+=segment(P[k],a,b,eps,N[k])[0]
                if own[k]<0 or col!=own[k]-n:v+=segment(P[k],c,d,eps,N[k])[0]
                C[k,col]=v
    return C
def matrix(P,N,X,eps,own):
    P=np.asarray(P,float);N=np.asarray(N,float);X=np.asarray(X,float);own=np.asarray(own,int);m=X.shape[0]-1;n=X.shape[1]-1
    if eps<0 or own.shape!=(len(P),) or N.shape!=P.shape or np.any((own<-1)|(own>=m*n)):raise ValueError('invalid_force_kernel_input')
    if not all(np.isfinite(v).all() for v in [P,N,X]):raise ValueError('nonfinite_force_geometry')
    use=own>=0;i=own[use]//n;j=own[use]%n;a=X[i,j];b=X[i,j+1];target=(a+b)/2;scale=np.maximum(np.maximum(np.linalg.norm(a,axis=1),np.linalg.norm(b,axis=1)),np.linalg.norm(b-a,axis=1))
    if np.any(np.linalg.norm(P[use]-target,axis=1)>64*np.finfo(float).eps*scale):raise ValueError('false_material_self_line')
    C=_matrix(P,N,X,eps,own)
    if not np.isfinite(C).all():raise ValueError('unresolved_force_kernel')
    return C
