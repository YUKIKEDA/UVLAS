"""Exact constant/delta wake change of basis; one real LU shared by frequencies."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.linalg import lu_factor,lu_solve
from numpy.polynomial.legendre import leggauss
from frequency_motion_columns import observe

def baseline(g,side,eps,length):
    X=g['X'].copy();X[-1,:,0]+=length
    if side=='candidate':
        from uvlm_frequency_fast import bound
        return bound(g['P'],g['normal'],X,eps)
    from ptera_frequency_reference import project
    return project(g,X[:-1,:-1].reshape(-1,3),X[:-1,1:].reshape(-1,3),X[1:,1:].reshape(-1,3),X[1:,:-1].reshape(-1,3),eps)

def delta(g,side,omega,eps,length,order):
    if side=='candidate':
        from uvlm_frequency_fast import age_nodes
        from uvlm_age_scalar import age_integral
        ages,w=age_nodes(g,omega,length,order);phase=omega*ages
        weights=abs(w)*(-2*np.sin(phase/2)**2-1j*np.sin(phase))
        return age_integral(g['P'],g['normal'],g['X'][-1],eps,ages,weights)
    from ptera_frequency_reference import project
    te=g['X'][-1];e=np.diff(te,axis=0);r=g['P'][:,None]-te[:-1][None];s=np.clip(np.einsum('pij,ij->pi',r,e)/np.sum(e*e,axis=1)[None],0,1);distance=np.sqrt(np.min(np.sum((r-s[:,:,None]*e[None])**2,axis=-1)));step=min(.0005,distance/40);T=length/20;nodes=[0.];a=0.
    while a<T:
        a=min(T,a+min(step,np.pi/max(omega,1e-10)));nodes.append(a);step*=1.4
    z,w=leggauss(order);lo=np.array(nodes[:-1]);hi=np.array(nodes[1:]);ages=((lo+hi)[:,None]/2+(hi-lo)[:,None]*z/2).ravel();weights=((hi-lo)[:,None]*w/2).ravel()*(-1j*omega)*np.exp(-1j*omega*ages);W=np.zeros((len(g['P']),g['n']),complex);tail=te+[length,0,0]
    for i in range(0,len(ages),8):
        tau=ages[i:i+8];front=te[None]+np.column_stack([20*tau,np.zeros((len(tau),2))])[:,None,:];back=np.broadcast_to(tail,front.shape);B=project(g,front[:,:-1].reshape(-1,3),front[:,1:].reshape(-1,3),back[:,1:].reshape(-1,3),back[:,:-1].reshape(-1,3),eps).reshape(len(g['P']),len(tau),g['n']);W+=np.einsum('q,pqn->pn',weights[i:i+8],B)
    return W

class Solver:
    def __init__(self,g,side,eps=1e-5,length=128.):
        self.g=g;self.side=side;self.eps=eps;self.length=length;self.A=baseline(g,side,eps,length);self.factor=lu_factor(self.A)
    def apply(self,rhs):return lu_solve(self.factor,rhs.real)+1j*lu_solve(self.factor,rhs.imag)
    def run(self,k,order=16):
        g=self.g;n=g['n'];omega=40*k/1.114285714285714;W=delta(g,self.side,omega,self.eps,self.length,order);nz=g['normal'][:,2];rhs=np.column_stack([1j*omega*nz,-(20+1j*omega*(g['P'][:,0]-.375))*nz]);Y=self.apply(W);small=lu_factor(np.eye(n)+Y[-n:])
        def invert(b):
            v=self.apply(b);return v-Y@lu_solve(small,v[-n:])
        gamma=invert(rhs);res=[]
        for iteration in range(4):
            R=rhs-self.A@gamma-W@gamma[-n:];res=np.linalg.norm(R,axis=0)/np.linalg.norm(rhs,axis=0)
            if max(res)<1e-11:break
            if iteration==3:raise ValueError(('per_rhs_residual_unresolved',res.tolist()))
            gamma+=invert(R)
        return dict(gamma=gamma,transfer=observe(g,gamma,omega,self.side),W=W,residuals=res.tolist(),refinement_steps=iteration)
