"""External body/traction moments for a parabolic-cylinder linear MMS.

No FE stiffness, displacement reference or strain-reference RHS is used.
Combine body and top/bottom polynomial coefficients before conversion to a
small resultant: their O(h) terms cancel identically, not numerically.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np

def curvature(x,mode):
    if mode=='standard':return np.column_stack((x*x,20*x**3,np.exp(x))),np.column_stack((2*x,60*x*x,np.exp(x)))
    a=np.array([2.,4.,6.]);g=1+x[:,None]**2*a*a
    return 1/g,-2*x[:,None]*a*a/g**2

def moment(h,p):
    return 0. if p%2 else h**(p+1)/(2**p*(p+1))

def stress_coefficients(x,k,v,vp,count):
    nu=.3;mu=1/(2*(1+nu));lam=nu/((1+nu)*(1-2*nu));g=1+k*k*x*x;gp=2*k*k*x
    a=np.column_stack((np.ones(len(x)),np.zeros(len(x)),k*x));ap=np.broadcast_to([0.,0.,k],a.shape)
    T=lam*g[:,None,None]*np.eye(3)+2*mu*a[:,:,None]*a[:,None,:]
    Tx=lam*gp[:,None,None]*np.eye(3)+2*mu*(ap[:,:,None]*a[:,None,:]+a[:,:,None]*ap[:,None,:])
    for m in range(1,count+1):
        scale=k**(m-1)/g**m
        S=-scale[:,None,None,None]*T[:,:,:,None]*v[:,None,None,:]
        Sx=-scale[:,None,None,None]*(Tx[:,:,:,None]*v[:,None,None,:]+T[:,:,:,None]*(vp-m*gp[:,None]/g[:,None]*v)[:,None,None,:])
        yield m,S,Sx

def moments(x,k,h,v,vp,count=14):
    """Return combined volume+top/bottom and signed-positive lateral moments.

    Scalar N=R+z D+z² C/2; returned component l already includes 1/l!.
    Lateral normals are +x and +y; caller supplies boundary signs.
    """
    a=np.column_stack((np.ones(len(x)),np.zeros(len(x)),k*x));d=np.column_stack((-k*x,np.zeros(len(x)),np.ones(len(x))));g=1+k*k*x*x
    W=np.zeros((3,len(x),3,v.shape[1]));Fx=np.zeros_like(W);Fy=np.zeros_like(W)
    for m,S,Sx in stress_coefficients(x,k,v,vp,count):
        Sa=np.einsum('pijt,pj->pit',S,a);Sxa=np.einsum('pijt,pj->pit',Sx,a);Sd=np.einsum('pijt,pj->pit',S,d);Sz=S[:,:,2,:]
        for l in range(3):
            fac=2 if l==2 else 1;mm=moment(h,m+l)
            W[l]+=(-Sxa*mm-k*(l+1)*Sz*mm+(l*Sd*moment(h,m+l-1) if l else 0))/fac
            Fx[l]+=Sa*mm/fac
            Fy[l]+=S[:,:,1,:]*(g*mm-k*moment(h,m+l+1))[:,None,None]/fac
    return W,Fx,Fy

def global_basis(x,k,order=32):
    """Eight smooth curvature modes, selected before any load evaluation."""
    g=1+k*k*x*x;t=2*x-1;P=np.polynomial.legendre.legvander(t,7)
    dP=np.column_stack([np.polynomial.legendre.legval(t,np.polynomial.legendre.legder(np.eye(8)[j]))*2 for j in range(8)])
    K=P/g[:,None]**2;Kx=dP/g[:,None]**2-4*k*k*x[:,None]*P/g[:,None]**3
    q,w=np.polynomial.legendre.leggauss(order);s=x[:,None]*(q+1)/2;B=np.polynomial.legendre.legvander(2*s-1,7)/(1+k*k*s*s)[:,:,None]**2
    I=[np.einsum('pq,pqn,q,p->pn',s**l,B,w,x/2) for l in range(3)];i0,i1,i2=I;b=x[:,None]*i0-i1;C=(x[:,None]**2*i0-i2)/2
    U=np.zeros((3,len(x),3,8));U[0,:,0]=-k*C;U[0,:,2]=b;U[1,:,0]=-i0;U[1,:,2]=-k*x[:,None]*i0
    return K,Kx,U

def global_fields(q,p,h,k,m):
    N,_=m.shape(p,0);x=(N@q.reshape(75,3)[:25])[:,0];z=h*p[:,2]/2;K,_,U=global_basis(x,k)
    a=np.column_stack((np.ones(len(x)),np.zeros(len(x)),k*x));T=-z[:,None,None,None]*a[:,:,None,None]*a[:,None,:,None]*K[:,None,None,:]/(1+k*k*x*x-k*z)[:,None,None,None]
    E=np.einsum('kij,pijn->pkn',m.MB,T)
    return E,U[0]+z[:,None,None]*U[1]
