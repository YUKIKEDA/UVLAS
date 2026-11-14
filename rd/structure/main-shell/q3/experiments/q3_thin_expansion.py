"""Mid-thickness compatible strain and its physical-z derivative."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.linalg import svd,solve_triangular
from scipy.optimize import brentq
import q3_curved_element as element
import q3_quartic_approximation as h1tool
RD=np.array([9*a+b for a in range(16) for b in range(6)])
CC=np.array([9*a+b for a in range(16) for b in range(6,9)])
def expansion(q,points):
    N,g=element.surface(points);coords=q.reshape(16,3,3);J0=np.zeros((len(points),3,3));J1=np.zeros_like(J0);J0[:,:,:2]=np.einsum('na,pnb->pab',coords[:,0],g);J0[:,:,2]=N@coords[:,1];J1[:,:,:2]=np.einsum('na,pnb->pab',coords[:,1],g);Ji=np.linalg.inv(J0);g0=np.zeros((len(points),48,3));g1=np.zeros_like(g0);g0[:,0::3,:2]=g;g0[:,1::3,2]=N;g1[:,1::3,:2]=g;g1[:,2::3,2]=N;grad0=g0@Ji;grad1=g1@Ji-grad0@J1@Ji;B0=np.einsum('kij,pnj->pkni',element.MB,grad0).reshape(-1,6,144);B1=np.einsum('kij,pnj->pkni',element.MB,grad1).reshape(-1,6,144);return J0,B0,B1
def frame_map(k,x):
    D=1+k*k*x*x;Q=np.zeros((len(x),3,3));Q[:,0,0]=Q[:,2,2]=1/np.sqrt(D);Q[:,1,1]=1;Q[:,2,0]=k*x/np.sqrt(D);Q[:,0,2]=-k*x/np.sqrt(D);t=np.einsum('pia,kab,pjb->pkij',Q,element.MB,Q);return np.einsum('pkij,lij->pkl',t,element.MB)
def target(kind,k,x):
    a=np.column_stack((np.ones_like(x),np.zeros_like(x),k*x));v=np.ones_like(x) if kind=='quadratic' else x*x;E=-v[:,None,None]*a[:,:,None]*a[:,None,:]/(1+k*k*x*x)[:,None,None];return np.einsum('kij,pij->pk',element.MB,E)
def field(kind,k,xy):
    x=xy[:,0];u=np.zeros((len(x),9))
    if kind=='quadratic':u[:,:3]=np.column_stack((-k*x**3/3,np.zeros(len(x)),x*x/2));u[:,3:6]=np.column_stack((-x,np.zeros(len(x)),-k*x*x))
    else:u[:,:3]=np.column_stack((-k*x**5/15,np.zeros(len(x)),x**4/12));u[:,3:6]=np.column_stack((-x**3/3,np.zeros(len(x)),-k*x**4/3))
    return u
def reference_h1(kind,k,order):
    p,w=np.polynomial.legendre.leggauss(order);x=(p+1)/2;w=w/2;D=1+k*k*x*x
    value=(x*x+k*k*x**4+(1+4*k*k*x*x)/D) if kind=='quadratic' else ((x**6+k*k*x**8)/9+(x**4+16*k*k*x**6/9)/D)
    return float(np.sqrt((w@(value*np.sqrt(D)))/(w@np.sqrt(D))))
def assemble(k,n,amplitude,order):
    xy,cells,coords,_=element.mesh(k,n,amplitude);free=np.flatnonzero(np.arange(len(xy))%(3*n+1));lookup={int(v):i for i,v in enumerate(free)};p,w=np.polynomial.legendre.leggauss(order);points=np.array(list(__import__('itertools').product(p,p)));ws=np.array([a*b for a,b in __import__('itertools').product(w,w)]);rows=len(cells)*len(points)*6;A0=np.zeros((rows,len(free)*6));C0=np.zeros_like(A0);A1=np.zeros_like(A0);D1=np.zeros((rows,len(free)*3));targets={kind:np.zeros(rows) for kind in ['quadratic','quartic']};zero_c=0.
    for ci,(_,_,ids) in enumerate(cells):
        J,B0,B1=expansion(coords[ids].reshape(144),points);det=np.linalg.det(J)
        if np.min(det)<=0:raise ValueError('nonpositive_reference_measure')
        weight=np.sqrt(ws*det);base,_=element.bilinear(points,xy[ids][[0,3,15,12]]);H=frame_map(k,base[:,0]);L0=np.einsum('pkl,pli->pki',H,B0);W1=np.einsum('kl,pli->pki',element.CH,B1);sl=slice(ci*len(points)*6,(ci+1)*len(points)*6);zero_c=max(zero_c,float(np.linalg.norm(B0[:,:,CC])))
        for local,node in enumerate(ids):
            if int(node) not in lookup:continue
            index=lookup[int(node)];A0[sl,6*index:6*index+6]=(B0[:,:,9*local:9*local+6]*weight[:,None,None]).reshape(-1,6);C0[sl,6*index:6*index+6]=(L0[:,:,9*local:9*local+6]*weight[:,None,None]).reshape(-1,6);A1[sl,6*index:6*index+6]=(W1[:,:,9*local:9*local+6]*weight[:,None,None]).reshape(-1,6);D1[sl,3*index:3*index+3]=(W1[:,:,9*local+6:9*local+9]*weight[:,None,None]).reshape(-1,3)
        for kind in targets:targets[kind][sl]=(target(kind,k,base[:,0])@element.CH.T*weight[:,None]).reshape(-1)
    M,area=h1tool.h1_matrix(xy,cells,coords,free,order)
    return dict(xy=xy,cells=cells,coords=coords,free=free,A0=A0,C0=C0,A1=A1,D1=D1,targets=targets,M=M,area=area,zero_c=zero_c)
def induced_metric(N,M):
    nf=len(M);nd=N.shape[1];D=N.reshape(nf,6,nd)[:,3:6];H=np.zeros((nd+nf*3,nd+nf*3))
    for c in range(3):H[:nd,:nd]+=D[:,c,:].T@M@D[:,c,:]
    H[nd:,nd:]=np.kron(M,np.eye(3));return H
def solve_ball(T,b,H,radius):
    L=np.linalg.cholesky(H);A=solve_triangular(L,T.T,lower=True,check_finite=False).T;U,s,Vh=svd(A,full_matrices=False,check_finite=False);ranks=[int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]]
    if ranks!=[A.shape[1]]*3:raise ValueError('first_order_rank_'+str(ranks))
    rhs=U.T@b
    def z(lam):return Vh.T@(s*rhs/(s*s+lam))
    uncon=z(0.);lam=0.
    if np.linalg.norm(uncon)>radius:
        high=1.
        while np.linalg.norm(z(high))>radius:high*=10
        lam=brentq(lambda v:np.linalg.norm(z(v))/radius-1,0.,high,xtol=1e-16,rtol=1e-13)
    zfit=z(lam);q=solve_triangular(L.T,zfit,check_finite=False);r=A@zfit-b;g=A.T@r+lam*zfit;kkt=float(np.linalg.norm(g)/(np.linalg.norm(A)*np.linalg.norm(b)+lam*np.linalg.norm(zfit)));return q,dict(error=float(np.linalg.norm(r)/np.linalg.norm(b)),h1=float(np.linalg.norm(zfit)),unconstrained_h1=float(np.linalg.norm(uncon)),multiplier=float(lam),kkt=kkt,ranks=ranks)
def full_field(a,N,q):
    u=np.zeros((len(a['xy']),9));nd=N.shape[1];u[a['free'],:6]=(N@q[:nd]).reshape(-1,6);u[a['free'],6:9]=q[nd:].reshape(-1,3);return u
