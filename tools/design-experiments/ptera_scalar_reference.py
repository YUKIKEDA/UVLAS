"""Independent Ptera endpoint formula projected to scalars, zero-age core only."""
import math,numpy as np
from numba import njit,prange
from numpy.polynomial.legendre import leggauss
@njit
def line(px,py,pz,ax,ay,az,bx,by,bz,eps,nx,ny,nz):
    ux=ax-px;uy=ay-py;uz=az-pz;vx=bx-px;vy=by-py;vz=bz-pz;dx=bx-ax;dy=by-ay;dz=bz-az;L2=dx*dx+dy*dy+dz*dz;L=math.sqrt(L2)
    if L<2.220446049250313e-16:return 0.
    ra=math.sqrt(ux*ux+uy*uy+uz*uz);rb=math.sqrt(vx*vx+vy*vy+vz*vz)
    if min(ra,rb)<L*1e-10:return 0.
    cx=uy*vz-uz*vy;cy=uz*vx-ux*vz;cz=ux*vy-uy*vx;c2=cx*cx+cy*cy+cz*cz;rr=ra*rb
    if math.sqrt(c2)<1e-10*rr:return 0.
    dot=ux*vx+uy*vy+uz*vz;f=(ra+rb)*(rr-dot)/(rr*(c2+L2*eps*eps))/(4*np.pi);return f*(nx*cx+ny*cy+nz*cz)
@njit
def semi(rx,ry,rz,eps,ny,nz):
    perp=ry*ry+rz*rz;r=math.sqrt(rx*rx+perp)
    if perp==0:return 0.
    if rx<0:f=perp/(r*(r-rx)*(perp+eps*eps))
    else:f=(1+rx/r)/(perp+eps*eps)
    return (-ny*rz+nz*ry)*f/(4*np.pi)
@njit
def horseshoe(P,N,a,b,eps,shift):
    ax=a[0]+shift;bx=b[0]+shift
    return line(P[0],P[1],P[2],ax,a[1],a[2],bx,b[1],b[2],eps,N[0],N[1],N[2])+semi(P[0]-bx,P[1]-b[1],P[2]-b[2],eps,N[1],N[2])-semi(P[0]-ax,P[1]-a[1],P[2]-a[2],eps,N[1],N[2])
@njit(parallel=True)
def matrix(P,N,X,eps):
    m=X.shape[0]-1;n=X.shape[1]-1;A=np.empty((len(P),m*n))
    for k in prange(len(P)):
        for i in range(m):
            for j in range(n):
                a=X[i,j];b=X[i,j+1];c=X[i+1,j+1];d=X[i+1,j];v=0.
                for aa,bb in [(a,b),(b,c),(c,d),(d,a)]:v+=line(P[k,0],P[k,1],P[k,2],aa[0],aa[1],aa[2],bb[0],bb[1],bb[2],eps,N[k,0],N[k,1],N[k,2])
                A[k,i*n+j]=v
    return A
@njit(parallel=True)
def integral(P,N,te,eps,length,ages,weights):
    W=np.empty((len(P),len(te)-1),np.complex128)
    for k in prange(len(P)):
        for j in range(len(te)-1):
            tail=horseshoe(P[k],N[k],te[j],te[j+1],eps,length);value=0j;comp=0j
            for i in range(len(ages)):
                term=weights[i]*(horseshoe(P[k],N[k],te[j],te[j+1],eps,20*ages[i])-tail);y=term-comp;v=value+y;comp=(v-value)-y;value=v
            W[k,j]=value
    return W
def baseline(g,side,eps,length):
    X=g['X'].copy();X[-1,:,0]+=length;return matrix(g['P'],g['normal'],X,eps)
def delta(g,side,omega,eps,length,order):
    te=g['X'][-1];e=np.diff(te,axis=0);r=g['P'][:,None]-te[:-1][None];s=np.clip(np.einsum('pij,ij->pi',r,e)/np.sum(e*e,axis=1)[None],0,1);distance=np.sqrt(np.min(np.sum((r-s[:,:,None]*e[None])**2,axis=-1)));step=min(.0005,distance/40);T=length/20;nodes=[0.];a=0.
    while a<T:a=min(T,a+min(step,np.pi/max(omega,1e-10)));nodes.append(a);step*=1.4
    z,w=leggauss(order);lo=np.array(nodes[:-1]);hi=np.array(nodes[1:]);ages=((lo+hi)[:,None]/2+(hi-lo)[:,None]*z/2).ravel();weights=((hi-lo)[:,None]*w/2).ravel()*(-1j*omega)*np.exp(-1j*omega*ages);W=integral(g['P'],g['normal'],te,eps,length,ages,weights)
    if not np.isfinite(W).all():raise ValueError('native_scalar_reference_unresolved')
    return W
