"""Research-only full-3D geometry-based stress projection; no product API."""
import itertools
import numpy as np
from scipy.linalg import qr, solve, solve_triangular, svdvals
import mq3_factor_family as reference

BASE=reference.base
C=BASE.C
CI=np.linalg.inv(C)
CIH=np.linalg.cholesky(CI).T
CH=BASE.CH
MB=BASE.MB


class PreparationFailure(ValueError):
    def __init__(self,report,witness):
        super().__init__('stress_generator_preparation_failed:'+str(report))
        self.report=report;self.witness=witness


def ranks(s):
    return [int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]]


def kinematics(q,h,points):
    N,d=BASE.geometry.shape(points);z=h*points[:,2]/2
    g=np.zeros((len(points),48,3),dtype=np.result_type(h,float))
    for power,factor in [(0,1.),(1,1.),(2,.5)]:
        g[:,power::3,:2]=d*factor*z[:,None,None]**power
        if power:g[:,power::3,2]=N*factor*power*z[:,None]**(power-1)
    J=np.einsum('na,pnb->pab',q.reshape(48,3),g)
    return N,J,g@np.linalg.inv(J)


def linear_raw(q,h,points):
    _,J,grad=kinematics(q,h,points)
    B=np.einsum('kij,pnj->pkni',MB,grad).reshape(-1,6,144)
    return J,np.concatenate((B,np.broadcast_to(np.eye(6),(len(points),6,6))),axis=2)


def green(state,q,h,points):
    _,J,grad=kinematics(state['q'],h,points)
    F=np.eye(3)+np.einsum('na,pnb->pab',(q-state['q']).reshape(48,3),grad)
    E=(F.transpose(0,2,1)@F-np.eye(3))/2
    B=np.einsum('kij,pai,pnj->pkna',MB,F,grad).reshape(-1,6,144)
    return np.einsum('kij,pij->pk',MB,E),B,F,grad,J


def geometry_state(q):
    N,J,_=kinematics(q,0.,np.zeros((1,3)))
    scale=np.linalg.norm(J[0],axis=0)
    if np.min(scale)<=0:raise ValueError('degenerate_midplane_frame')
    return dict(q=q,frame_scale=scale,length=float(np.sqrt(scale[0]*scale[1])),center=N@q.reshape(16,3,3)[:,0])


def stress_generator(state,points):
    N,J,_=kinematics(state['q'],0.,points)
    A=J/state['frame_scale'][None,None,:];det=np.linalg.det(A)
    if np.min(det)<=0:raise ValueError('nonpositive_midplane_frame')
    polynomials=[np.polynomial.legendre.legvander(points[:,axis],3)*np.sqrt((2*np.arange(4)+1)/2) for axis in range(3)]
    hat=np.zeros((len(points),6,203));col=0
    for component,imax,jmax,kmax in [(0,2,3,2),(1,3,2,2),(2,3,3,1),(3,3,2,2),(4,2,3,2),(5,2,2,2)]:
        for i,j,k in itertools.product(range(imax+1),range(jmax+1),range(kmax+1)):
            hat[:,component,col]=polynomials[0][:,i]*polynomials[1][:,j]*polynomials[2][:,k];col+=1
    assert col==203
    tensor=np.einsum('kij,pkn->pijn',MB,hat)
    mapped=np.einsum('pai,pijn,pbj->pabn',A,tensor,A,optimize=True)/det[:,None,None,None]
    mapped=np.einsum('kab,pabn->pkn',MB,mapped)
    coords=state['q'].reshape(16,3,3)
    phi=np.column_stack((np.ones(len(points)),(N@coords[:,0]-state['center'])/state['length'],
        points[:,2,None]*(N@coords[:,1]),state['length']*points[:,2,None]**2*(N@coords[:,2])))
    physical=(phi[:,None,:,None]*np.eye(6)[None,:,None,:]).reshape(len(points),6,60)
    return np.concatenate((mapped,physical),axis=2)


def measure(state,h,points,w):
    _,J,_=kinematics(state['q'],h,points)
    weight=w*np.linalg.det(J)
    if not np.iscomplexobj(h) and np.min(np.linalg.det(J))<=0:raise ValueError('nonpositive_reference_volume')
    return weight/weight.sum(),weight.sum()*h/2


def prepare(q):
    state=geometry_state(q);p,w=BASE.quadrature(9)
    S=stress_generator(state,p);weight,_=measure(state,0.,p,w)
    A=(np.einsum('ij,pjn->pin',CIH,S)*np.sqrt(weight)[:,None,None]).reshape(-1,263)
    norms=np.linalg.norm(A,axis=0);active=np.flatnonzero(norms>1e-14*norms.max())
    normalized=A[:,active]/norms[active]
    singular=svdvals(normalized,check_finite=False)
    _,R,pivot=qr(normalized,mode='economic',pivoting=True,check_finite=False)
    sr,rr=ranks(singular),ranks(abs(np.diag(R)))
    report=dict(stress_ranks_svd=sr,stress_ranks_qr=rr,active_columns=len(active),generated_columns=263)
    if sr!=rr or sr!=[sr[1]]*3:
        raise PreparationFailure(report,dict(q=q,normalized_generator=normalized,svd_singular=singular,qr_diagonal=np.diag(R),active=active,column_norms=norms))
    rank=sr[1];selected=active[pivot[:rank]]
    T=np.zeros((263,rank))
    T[selected]=solve_triangular(R[:rank,:rank],np.eye(rank),check_finite=False)/norms[selected,None]
    state.update(T=T,active=active,column_norms=norms,stress_ranks=sr,svd_singular=singular,qr_diagonal=np.diag(R))
    return state,report


def system(state,h,order=9):
    p,w=BASE.quadrature(order);weight,volume=measure(state,h,p,w)
    S=stress_generator(state,p)@state['T']
    H=np.einsum('p,pki,kl,plj->ij',weight,S,CI,S,optimize=True)
    return p,weight,volume,S,H


def project(state,h,E,order=9):
    p,weight,volume,S,H=system(state,h,order)
    rhs=np.einsum('p,pki,pkn->in',weight,S,E,optimize=True)
    a=solve(H,rhs,assume_a='sym',check_finite=False)
    sigma=S@a
    return dict(a=a,sigma=sigma,R=np.einsum('ij,pjn->pin',CI,sigma),H=H,rhs=rhs,
        S=S,points=p,weight=weight,volume=volume)


def adjoint(P,tau):
    rhs=np.einsum('p,pki,kl,pln->in',P['weight'],P['S'],CI,tau,optimize=True)
    return P['S']@solve(P['H'],rhs,assume_a='sym',check_finite=False)


def energy_force(state,q,h):
    p,w=BASE.quadrature(9)
    E,B,F,grad,J=green(state,q,h,p)
    P=project(state,h,E[:,:,None])
    sigma=P['sigma'][:,:,0]
    energy=(P['rhs'][:,0].T@P['a'][:,0])/2
    force=np.einsum('p,pki,pk->i',P['weight'],B,sigma,optimize=True)
    return energy,force,P,B,F,grad


def variations(state,q,h,v1,v2):
    energy,force,P,B,F,grad=energy_force(state,q,h)
    d1=np.einsum('na,pnb->pab',v1.reshape(48,3),grad)
    d2=np.einsum('na,pnb->pab',v2.reshape(48,3),grad)
    l1=np.einsum('p,pki,pk->i',P['weight'],P['S'],B@v1,optimize=True)
    l2=np.einsum('p,pki,pk->i',P['weight'],P['S'],B@v2,optimize=True)
    material=l1.T@solve(P['H'],l2,assume_a='sym',check_finite=False)
    d2E=(d1.transpose(0,2,1)@d2+d2.transpose(0,2,1)@d1)/2
    geometric=np.einsum('p,pk,kij,pij->',P['weight'],P['sigma'][:,:,0],MB,d2E,optimize=True)
    return energy,force,material,geometric
