"""Physical-z M-Q3 family; fixed-basis KKT differentiation, no singular J_zeta."""
import numpy as np
from scipy.linalg import qr,solve_triangular,solve
import mq3_independent as base
DIM=117
def raw(q,h,points):
    N,d=base.geometry.shape(points);z=h*points[:,2]/2;dtype=np.result_type(h,float);g=np.zeros((len(z),48,3),dtype=dtype)
    for node in range(16):
        for power,factor in [(0,1.),(1,1.),(2,.5)]:
            col=3*node+power;g[:,col,:2]=d[:,node]*factor*z[:,None]**power
            if power:g[:,col,2]=N[:,node]*factor*power*z**(power-1)
    J=np.einsum('na,pnb->pab',q.reshape(48,3),g);grad=g@np.linalg.inv(J);B=np.einsum('kij,pnj->pkni',base.MB,grad).reshape(-1,6,144);E=np.concatenate((B,np.broadcast_to(np.eye(6),(len(z),6,6))),axis=2);return J,E
def data(q,h,points):
    points=np.asarray(points);J,E=raw(q,h,points);J0,E0=raw(q,0.,points);cov=base.covariant(J,E);covm=base.covariant(J0,E0);m=np.stack((covm[:,0,0],covm[:,1,1],covm[:,0,1]),axis=1);g=np.stack((cov[:,0,2],cov[:,1,2]),axis=1);M=np.zeros((len(points),3,DIM),dtype=E.dtype);S=np.zeros((len(points),2,DIM),dtype=E.dtype);s,t,zeta=points.T;col=0
    for component,imax,jmax in [(0,2,3),(1,3,2),(2,2,2)]:
        for j in range(jmax+1):
            for i in range(imax+1):M[:,component,col]=s**i*t**j;col+=1
    M[:,:,33:39]=m[:,:,144:];col=39
    for component,imax,jmax in [(0,2,3),(1,3,2)]:
        for k in range(3):
            for j in range(jmax+1):
                for i in range(imax+1):S[:,component,col]=s**i*t**j*zeta**k;col+=1
    S[:,:,111:]=g[:,:,144:];covW=np.zeros((len(points),DIM,3,3),dtype=E.dtype);covF=np.zeros((len(points),150,3,3),dtype=E.dtype)
    for dst,mem,sh in [(covW,M,S),(covF,m,g)]:
        dst[:,:,0,0]=mem[:,0];dst[:,:,1,1]=mem[:,1];dst[:,:,0,1]=dst[:,:,1,0]=mem[:,2];dst[:,:,0,2]=dst[:,:,2,0]=sh[:,0];dst[:,:,1,2]=dst[:,:,2,1]=sh[:,1]
    return dict(J=J,E=E,W=base.push(J,covW.transpose(0,2,3,1)),F=base.push(J,covF.transpose(0,2,3,1)),M=M,S=S,m=m,g=g)
def moments(q,h,d,weight):
    p,w=np.polynomial.legendre.leggauss(9);F=[];rhs=[]
    for tangent in range(2):
        for sign in [-1.,1.]:
            pts=np.zeros((9,3));pts[:,tangent]=p;pts[:,1-tangent]=sign;edge=data(q,h,pts)
            for degree in range(3):F.append((w*p**degree)@edge['M'][:,tangent]);rhs.append((w*p**degree)@edge['m'][:,tangent])
            face=np.zeros((81,3));face[:,tangent]=np.repeat(p,9);face[:,1-tangent]=sign;face[:,2]=np.tile(p,9);side=data(q,h,face)
            for i in range(3):
                for j in range(3):
                    wf=np.outer(w*p**i,w*p**j).reshape(-1);F.append(wf@side['S'][:,tangent]);rhs.append(wf@side['g'][:,tangent])
    F.extend(np.einsum('p,pkn->kn',weight,d['W']));rhs.extend(np.einsum('p,pkn->kn',weight,d['F']));return np.array(F),np.array(rhs)
def volume(q,h,order=9):
    p,w=base.quadrature(order);d=data(q,h,p);weight=w*np.linalg.det(d['J']);weight/=weight.sum();return p,d,weight
def choose(q):
    _,d,w=volume(q,0.);A=(np.einsum('kl,pln->pkn',base.CH,d['W'])*np.sqrt(w)[:,None,None]).reshape(-1,DIM);norm=np.linalg.norm(A,axis=0);active=np.flatnonzero(norm>1e-14*norm.max());Q,R,pivot=qr(A[:,active]/norm[active],mode='economic',pivoting=True,check_finite=False);ranks=base.ranks(np.diag(R));rank=ranks[1];selected=active[pivot[:rank]];T=np.zeros((DIM,rank));T[selected]=solve_triangular(R[:rank,:rank],np.eye(rank),check_finite=False)/norm[selected,None];F,g=moments(q,0.,d,w);D=F@T;rownorm=np.linalg.norm(D,axis=1);active_rows=np.flatnonzero(rownorm>1e-13*rownorm.max());Fn=D[active_rows]/rownorm[active_rows,None];Qf,Rf,piv=qr(Fn.T,mode='economic',pivoting=True,check_finite=False);cranks=base.ranks(np.diag(Rf));r=cranks[1];rows=active_rows[piv[:r]];scale=1/rownorm[rows];Qd,Rd=qr((D[rows]*scale[:,None]).T,mode='economic',check_finite=False);U=solve_triangular(Rd.T,np.eye(r),lower=True,check_finite=False)*scale[None,:]
    return dict(q=q,T=T,rows=rows,U=U,space_ranks=ranks,constraint_ranks=cranks)
def system(state,h):
    _,d,w=volume(state['q'],h);W=d['W']@state['T'];H=np.einsum('p,pki,kl,plj->ij',w,W,base.C,W,optimize=True);l=np.einsum('p,pki,kl,pln->in',w,W,base.C,d['F'],optimize=True);F,g=moments(state['q'],h,d,w);D=state['U']@(F[state['rows']]@state['T']);rhs=state['U']@g[state['rows']];n=len(H);m=len(D);S=np.zeros((n+m,n+m),dtype=H.dtype);S[:n,:n]=H;S[:n,n:]=D.T;S[n:,:n]=D;return S,np.vstack((l,rhs)),F,g
def coefficients(state,h):
    S,b,F,g=system(state,h);X=solve(S,b,assume_a='sym',check_finite=False);return X,S,b,F,g
def evaluate(state,h,X,points):
    d=data(state['q'],h,points);c=state['T']@X[:state['T'].shape[1]];return d['E']-d['F']+d['W']@c
