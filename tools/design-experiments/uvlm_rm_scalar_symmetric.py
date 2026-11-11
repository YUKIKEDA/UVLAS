import math,numpy as np
from numba import njit,prange
@njit
def line(P,N,a,b,eps):
    dx=b[0]-a[0];dy=b[1]-a[1];dz=b[2]-a[2];L=math.sqrt(dx*dx+dy*dy+dz*dz);ex=dx/L;ey=dy/L;ez=dz/L;rx=P[0]-a[0];ry=P[1]-a[1];rz=P[2]-a[2];s=rx*ex+ry*ey+rz*ez;t=s-L;px=rx-s*ex;py=ry-s*ey;pz=rz-s*ez;k2=px*px+py*py+pz*pz+eps*eps;aa=math.sqrt(s*s+k2);bb=math.sqrt(t*t+k2)
    if s*t>0:I=L*(s+t)/(aa*bb*(s*bb+t*aa))
    else:I=(s/aa-t/bb)/k2
    C=N[0]*(ey*pz-ez*py)+N[1]*(ez*px-ex*pz)+N[2]*(ex*py-ey*px);return C*I/(4*np.pi)
@njit
def ring(P,N,X,i,j,eps):
    a=X[i,j];b=X[i,j+1];c=X[i+1,j+1];d=X[i+1,j];return line(P,N,a,b,eps)+line(P,N,b,c,eps)+line(P,N,c,d,eps)+line(P,N,d,a,eps)
@njit(parallel=True)
def condensed(P,N,X,eps):
    m=X.shape[0]-1;n=X.shape[1]-1;nh=n//2;A=np.empty((len(P),m*nh))
    for k in prange(len(P)):
        for i in range(m):
            for j in range(nh):A[k,i*nh+j]=ring(P[k],N[k],X,i,nh+j,eps)+ring(P[k],N[k],X,i,nh-1-j,eps)
    return A
