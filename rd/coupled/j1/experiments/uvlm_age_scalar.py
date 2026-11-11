"""Allocation-free scalar evaluation of the already qualified point derivative."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import math,numpy as np
from numba import njit,prange
@njit
def point_derivative(rx,ry,rz,ex,ey,ez,L,eps,nx,ny,nz):
    s=rx*ex+ry*ey+rz*ez;t=s-L;px=rx-s*ex;py=ry-s*ey;pz=rz-s*ez;k2=px*px+py*py+pz*pz+eps*eps;aa=math.sqrt(rx*rx+ry*ry+rz*rz+eps*eps);bx=rx-L*ex;by=ry-L*ey;bz=rz-L*ez;bb=math.sqrt(bx*bx+by*by+bz*bz+eps*eps)
    if aa==0 or bb==0:return np.nan
    ax=rx/aa;bbx=bx/bb
    if s*t>0:
        D=s*bb+t*aa;den=aa*bb*D
        if not math.isfinite(den) or den==0:return np.nan
        I=L*(s+t)/den;Ix=I*(2*ex/(s+t)-ax/aa-bbx/bb-(ex*(aa+bb)+s*bbx+t*ax)/D)
    else:
        if k2==0:return np.nan
        I=(s/aa-t/bb)/k2;Ix=(ex*(1/aa-1/bb)-s*ax/(aa*aa)+t*bbx/(bb*bb))/k2-I*2*px/k2
    C=nx*(ey*rz-ez*ry)+ny*(ez*rx-ex*rz)+nz*(ex*ry-ey*rx);Cx=ny*ez-nz*ey;right=(-ny*bz+nz*by)/(bb*bb*bb);left=(-ny*rz+nz*ry)/(aa*aa*aa)
    return (Cx*I+C*Ix+right-left)/(4*np.pi)
@njit(parallel=True)
def _integral(P,N,te,eps,ages,weights):
    W=np.empty((len(P),len(te)-1),np.complex128)
    for k in prange(len(P)):
        for j in range(len(te)-1):
            dx=te[j+1,0]-te[j,0];dy=te[j+1,1]-te[j,1];dz=te[j+1,2]-te[j,2];L=math.sqrt(dx*dx+dy*dy+dz*dz);ex=dx/L;ey=dy/L;ez=dz/L;rx=P[k,0]-te[j,0];ry=P[k,1]-te[j,1];rz=P[k,2]-te[j,2];value=0j;comp=0j
            for q in range(len(ages)):
                v=20*weights[q]*point_derivative(rx-20*ages[q],ry,rz,ex,ey,ez,L,eps,N[k,0],N[k,1],N[k,2]);y=v-comp;z=value+y;comp=(z-value)-y;value=z
            W[k,j]=value
    return W
def age_integral(P,N,te,eps,ages,weights):
    if eps<0 or np.asarray(P).shape!=np.asarray(N).shape or len(ages)!=len(weights) or np.any(np.linalg.norm(np.diff(te,axis=0),axis=1)==0):raise ValueError('invalid_age_geometry')
    if not all(np.isfinite(a).all() for a in [P,N,te,ages,weights]):raise ValueError('nonfinite_age_input')
    W=_integral(P,N,te,eps,ages,weights)
    if not np.isfinite(W).all():raise ValueError('age_derivative_precision_unresolved')
    return W
