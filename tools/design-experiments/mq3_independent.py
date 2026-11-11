"""Monomial/Vandermonde and QR-nullspace implementation of M-Q3."""
import itertools
import numpy as np
from scipy.linalg import qr,solve_triangular
import q3_curved_independent as geometry
DIM=117
MB=np.zeros((6,3,3))
for i in range(3):MB[i,i,i]=1
for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):MB[k,i,j]=MB[k,j,i]=1/np.sqrt(2.)
MU=1/2.6;LAM=.3/(1.3*.4);C=2*MU*np.eye(6);C[:3,:3]+=LAM;CH=np.linalg.cholesky(C).T
def quadrature(order):
    p,w=np.polynomial.legendre.leggauss(order);return np.array(list(itertools.product(p,p,p))),np.array([a*b*c for a,b,c in itertools.product(w,w,w)])
def strain_matrix(grad):
    B=np.zeros((len(grad),6,144))
    for i in range(3):B[:,i,i::3]=grad[:,:,i]
    for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):B[:,k,i::3]=grad[:,:,j]/np.sqrt(2.);B[:,k,j::3]=grad[:,:,i]/np.sqrt(2.)
    return B
def covariant(J,E):return np.einsum('pia,kij,pjb,pkn->pabn',J,MB,J,E,optimize=True)
def push(J,cov):
    inv=np.linalg.inv(J);t=inv.transpose(0,2,1)[:,None]@cov.transpose(0,3,1,2)@inv[:,None]
    return np.einsum('kij,pnij->pkn',MB,t,optimize=True)
def data(q,h,points):
    points=np.asarray(points);_,J,g=geometry.data(q,h,points);B=strain_matrix(g);J=J.copy();J[:,:,2]*=2/h;E=np.concatenate((B,np.broadcast_to(np.eye(6),(len(points),6,6))),axis=2);mid=points.copy();mid[:,2]=0;_,Jm,gm=geometry.data(q,h,mid);Jm=Jm.copy();Jm[:,:,2]*=2/h;Em=np.concatenate((strain_matrix(gm),np.broadcast_to(np.eye(6),(len(points),6,6))),axis=2);cov=covariant(J,E);covm=covariant(Jm,Em);m=np.stack((covm[:,0,0],covm[:,1,1],covm[:,0,1]),axis=1);shear=np.stack((cov[:,0,2],cov[:,1,2]),axis=1);Gm=m[:,:,144:];Gs=shear[:,:,144:];M=np.zeros((len(points),3,DIM));S=np.zeros((len(points),2,DIM));x,y,z=points.T;column=0
    for component,imax,jmax in [(0,2,3),(1,3,2),(2,2,2)]:
        for j,i in itertools.product(range(jmax+1),range(imax+1)):M[:,component,column]=x**i*y**j;column+=1
    M[:,:,33:39]=Gm;column=39
    for component,imax,jmax in [(0,2,3),(1,3,2)]:
        for k,j,i in itertools.product(range(3),range(jmax+1),range(imax+1)):S[:,component,column]=x**i*y**j*z**k;column+=1
    S[:,:,111:117]=Gs;Wcov=np.zeros((len(points),DIM,3,3));Fcov=np.zeros((len(points),150,3,3))
    for dst,mem,sh in [(Wcov,M,S),(Fcov,m,shear)]:
        dst[:,:,0,0]=mem[:,0];dst[:,:,1,1]=mem[:,1];dst[:,:,0,1]=dst[:,:,1,0]=mem[:,2];dst[:,:,0,2]=dst[:,:,2,0]=sh[:,0];dst[:,:,1,2]=dst[:,:,2,1]=sh[:,1]
    return dict(J=J,E=E,W=push(J,Wcov.transpose(0,2,3,1)),F=push(J,Fcov.transpose(0,2,3,1)),M=M,S=S,m=m,g=shear)
def constraints(q,h,d,weight):
    p,w=np.polynomial.legendre.leggauss(9);F=[];rhs=[]
    for axis in range(2):
        for sign in [-1.,1.]:
            pts=np.zeros((9,3));pts[:,axis]=p;pts[:,1-axis]=sign;edge=data(q,h,pts)
            for power in range(3):
                weights=w*p**power;F.append(weights@edge['M'][:,axis]);rhs.append(weights@edge['m'][:,axis])
            face=np.zeros((81,3));face[:,axis]=np.repeat(p,9);face[:,1-axis]=sign;face[:,2]=np.tile(p,9);side=data(q,h,face)
            for a,b in itertools.product(range(3),range(3)):
                weights=np.outer(w*p**a,w*p**b).reshape(-1);F.append(weights@side['S'][:,axis]);rhs.append(weights@side['g'][:,axis])
    F.extend(np.einsum('p,pkn->kn',weight,d['W']));rhs.extend(np.einsum('p,pkn->kn',weight,d['F']));return np.array(F),np.array(rhs)
def ranks(diagonal):
    s=abs(diagonal);return [int(np.count_nonzero(s>v*s[0])) for v in [1e-9,1e-10,1e-11]]
def prepare(q,h):
    p,w=quadrature(9);d=data(q,h,p);weight=w*np.linalg.det(d['J']);weight/=weight.sum();A=(np.einsum('ij,pjn->pin',CH,d['W'])*np.sqrt(weight)[:,None,None]).reshape(-1,DIM);Y=(np.einsum('ij,pjn->pin',CH,d['F'])*np.sqrt(weight)[:,None,None]).reshape(-1,150);norm=np.linalg.norm(A,axis=0);active=np.flatnonzero(norm>1e-14*norm.max());Q,R,piv=qr(A[:,active]/norm[active],mode='economic',pivoting=True,check_finite=False);br=ranks(np.diag(R));rank=br[1];selected=active[piv[:rank]];T=np.zeros((DIM,rank));T[selected]=solve_triangular(R[:rank,:rank],np.eye(rank),check_finite=False)/norm[selected,None];coef0=Q[:,:rank].T@Y;F,g=constraints(q,h,d,weight);F=F@T;rownorm=np.linalg.norm(F,axis=1);rowactive=np.flatnonzero(rownorm>1e-13*rownorm.max());Fn=F[rowactive]/rownorm[rowactive,None];gn=g[rowactive]/rownorm[rowactive,None];Qp,Rp,indices=qr(Fn.T,mode='economic',pivoting=True,check_finite=False);cr=ranks(np.diag(Rp));nconstraints=cr[1];independent=indices[:nconstraints];Frows=Fn[independent];grows=gn[independent];Qc,Rc=qr(Frows.T,mode='full',check_finite=False);part=Qc[:,:nconstraints]@solve_triangular(Rc[:nconstraints,:nconstraints].T,grows,lower=True,check_finite=False);null=Qc[:,nconstraints:];coef=part+null@(null.T@coef0);coeff=T@coef;residual=float(np.linalg.norm(F@coef-g)/max(1.,np.linalg.norm(g)));zero=float(np.linalg.norm(g[np.setdiff1d(np.arange(len(g)),rowactive)])/max(1.,np.linalg.norm(g)))
    return dict(q=q,h=h,coeff=coeff,basis_ranks=br,constraint_ranks=cr,constraint_residual=residual,zero_rhs=zero,basis_diagonal=np.diag(R),constraint_diagonal=np.diag(Rp))
def evaluate(state,points):
    d=data(state['q'],state['h'],points);return d['E']-d['F']+d['W']@state['coeff']
