"""Linear homogeneous-C EAS21 prototype; not a finite-strain formulation."""
import itertools
import numpy as np
from scipy.linalg import qr, solve_triangular
import check_q9_curved_strain_screen as raw
import check_q9_independent_load as common

MB=raw.EB.copy(); MB[3:]*=np.sqrt(2.)
S=np.diag([1.,1.,1.,np.sqrt(2.),np.sqrt(2.),np.sqrt(2.)])
C=S@common.elastic@S
CH=np.linalg.cholesky(C).T

def quadrature(order):
    p,w=np.polynomial.legendre.leggauss(order)
    return np.array(list(itertools.product(p,p,p))),np.array([a*b*c for a,b,c in itertools.product(w,w,w)])

def modes(points):
    s,t=points[:,:2].T
    def cardinal(x):
        return np.array([x*(x-1)/2,1-x*x,x*(x+1)/2]).T,np.array([x-.5,-2*x,x+.5]).T
    ls,ds=cardinal(s);lt,dt=cardinal(t);gs=s*(s*s-1);gt=t*(t*t-1);dgs=3*s*s-1;dgt=3*t*t-1
    f=np.column_stack((gs[:,None]*lt,ls*gt[:,None],gs*gt))
    g=np.zeros((len(s),7,3))
    g[:,:3,0]=dgs[:,None]*lt;g[:,:3,1]=gs[:,None]*dt
    g[:,3:6,0]=ds*gt[:,None];g[:,3:6,1]=ls*dgt[:,None]
    g[:,6,0]=dgs*gt;g[:,6,1]=gs*dgt
    return f,g

def source(q,h,points):
    J,B=raw.kinematics(q,h,0.,points);B=B.copy();B[:,3:]*=1/np.sqrt(2.)
    _,g=modes(points);grad=np.einsum('pna,paj->pnj',g,np.linalg.inv(J))
    M=np.einsum('kij,pnj->pkni',MB,grad).reshape(len(points),6,21)
    return J,B,M

def prepare(q,h,order=9):
    p,w=quadrature(order);J,B,M=source(q,h,p);weight=w*np.linalg.det(J);volume=weight.sum()
    mean=np.einsum('p,pki->ki',weight,M)/volume;M=M-mean
    A=(np.einsum('ij,pjk->pik',CH,M)*np.sqrt(weight)[:,None,None]).reshape(-1,21)
    scale=1/np.linalg.norm(A,axis=0);Q,R=qr(A*scale,mode='economic')
    W=(np.einsum('ij,pjk->pik',CH,B)*np.sqrt(weight)[:,None,None]).reshape(-1,81)
    coefficient=scale[:,None]*solve_triangular(R,Q.T@W)
    return dict(q=q,h=h,mean=mean,coefficient=coefficient,mode_matrix=A,volume=volume)

def evaluate(state,points):
    J,B,M=source(state['q'],state['h'],points);M=M-state['mean']
    return J,B,B-np.einsum('pki,in->pkn',M,state['coefficient']),M
