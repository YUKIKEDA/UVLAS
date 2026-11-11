"""Experimental acceleration of the same RM kernel and its exact point derivative."""
import math
import numpy as np
from numba import njit,prange
from scipy.linalg import solve
from numpy.polynomial.legendre import leggauss
from uvlm_steady_reference import geometry
@njit
def segment(p,a,b,eps,n):
    dl=b-a;L=np.sqrt(np.dot(dl,dl));e=dl/L;r=p-a;s=np.dot(r,e);t=s-L;perp=r-s*e;k2=np.dot(perp,perp)+eps*eps;aa=np.sqrt(s*s+k2);bb=np.sqrt(t*t+k2);ax=r[0]/aa;bx=(r[0]-L*e[0])/bb
    if s*t>0:
        D=s*bb+t*aa;I=L*(s+t)/(aa*bb*D);Ix=I*(2*e[0]/(s+t)-ax/aa-bx/bb-(e[0]*(aa+bb)+s*bx+t*ax)/D)
    else:
        I=(s/aa-t/bb)/k2;Ix=(e[0]*(1/aa-1/bb)-s*ax/(aa*aa)+t*bx/(bb*bb))/k2-I*2*perp[0]/k2
    C=np.cross(e,perp);Cx=np.array([0.,e[2],-e[1]])
    return np.dot(n,C)*I/(4*np.pi),(np.dot(n,Cx)*I+np.dot(n,C)*Ix)/(4*np.pi)
@njit(parallel=True)
def bound(P,N,X,eps):
    m=X.shape[0]-1;n=X.shape[1]-1;A=np.empty((len(P),m*n))
    for k in prange(len(P)):
        for i in range(m):
            for j in range(n):
                a=X[i,j];b=X[i,j+1];c=X[i+1,j+1];d=X[i+1,j]
                A[k,i*n+j]=segment(P[k],a,b,eps,N[k])[0]+segment(P[k],b,c,eps,N[k])[0]+segment(P[k],c,d,eps,N[k])[0]+segment(P[k],d,a,eps,N[k])[0]
    return A
@njit(parallel=True)
def age_integral(P,N,te,eps,ages,weights):
    n=len(te)-1;W=np.zeros((len(P),n),np.complex128)
    for k in prange(len(P)):
        for j in range(n):
            value=0j
            for q in range(len(ages)):
                p=P[k].copy();p[0]-=20*ages[q];_,dv=segment(p,te[j],te[j+1],eps,N[k]);rr=p-te[j+1];rl=p-te[j];right=(N[k,2]*rr[1]-N[k,1]*rr[2])/(np.dot(rr,rr)+eps*eps)**1.5;left=(N[k,2]*rl[1]-N[k,1]*rl[2])/(np.dot(rl,rl)+eps*eps)**1.5;value+=weights[q]*20*(dv+(right-left)/(4*np.pi))
            W[k,j]=value
    return W
def age_nodes(g,omega,length,order):
    te=g['X'][-1];a=te[:-1];e=te[1:]-a;r=g['P'][:,None]-a[None];t=np.clip(np.einsum('pij,ij->pi',r,e)/np.sum(e*e,axis=1)[None],0,1);distance=np.sqrt(np.min(np.sum((r-t[:,:,None]*e[None])**2,axis=-1)));step=min(.0005,distance/20/2);T=length/20;nodes=[0.];x=0.
    while x<T:x=min(T,x+min(step,np.pi/max(abs(omega),1e-8)));nodes.append(x);step*=1.5
    z,w=leggauss(order);lo=np.array(nodes[:-1]);hi=np.array(nodes[1:]);ages=((lo+hi)[:,None]/2+(hi-lo)[:,None]*z/2).ravel();weights=((hi-lo)[:,None]*w/2).ravel()*np.exp(-1j*omega*ages);return ages,weights
def run(m,n,k=.2,eps=1e-5,length=128.,order=12,motion='heave',axis=.375):
    g=geometry(m,n);omega=2*k*20/1.114285714285714;ages,weights=age_nodes(g,omega,length,order);W=age_integral(g['P'],g['normal'],g['X'][-1],eps,ages,weights);A=bound(g['P'],g['normal'],g['X'],eps).astype(complex);A[:,-n:]+=W;rhs=1j*omega*g['normal'][:,2] if motion=='heave' else -(20+1j*omega*(g['P'][:,0]-axis))*g['normal'][:,2];gamma=solve(A,rhs);G=gamma.reshape(m,n);horse=(G-np.vstack([np.zeros((1,n)),G[:-1]])).ravel();edgeF=1.2*horse[:,None]*np.cross([20.,0,0],g['right']-g['left']);F=edgeF.sum(axis=0);M=np.cross((g['left']+g['right'])/2-[axis,0,0],edgeF).sum(axis=0);X=g['X']
    for i in range(m):
        for j in range(n):
            q=X[[i,i,i+1,i+1],[j,j+1,j+1,j]]
            for ids in [[0,2,1],[0,3,2]]:
                tri=q[ids];f=1.2j*omega*gamma[i*n+j]*np.cross(tri[1]-tri[0],tri[2]-tri[0])/2;F+=f;M+=np.cross(tri.mean(axis=0)-[axis,0,0],f)
    return dict(g=g,gamma=gamma,F=F,M=M,W=W,residual=float(np.linalg.norm(A@gamma-rhs)/np.linalg.norm(rhs)))
