"""Monomial/Vandermonde shape differentiation and tensor-based work path."""
import itertools
import numpy as np
NODES=np.array([-1.,-1/3,1/3,1.]);COEFFICIENT=np.linalg.inv(np.vander(NODES,4,increasing=True))
def shape(points):
    s,t=points[:,:2].T
    def rows(x):
        return np.column_stack((np.ones(len(x)),x,x*x,x**3))@COEFFICIENT,np.column_stack((np.zeros(len(x)),np.ones(len(x)),2*x,3*x*x))@COEFFICIENT
    ls,ds=rows(s);lt,dt=rows(t);N=[];gs=[];gt=[]
    for j,i in itertools.product(range(4),range(4)):N.append(ls[:,i]*lt[:,j]);gs.append(ds[:,i]*lt[:,j]);gt.append(ls[:,i]*dt[:,j])
    return np.array(N).T,np.stack((np.array(gs).T,np.array(gt).T),axis=2)
def data(q,h,points):
    N,deriv=shape(np.asarray(points));z=h*np.asarray(points)[:,2]/2;phi=np.zeros((len(z),48));g=np.zeros((len(z),48,3))
    for a in range(16):
        for power,factor in [(0,1.),(1,1.),(2,.5)]:
            col=3*a+power;phi[:,col]=N[:,a]*factor*z**power;g[:,col,:2]=deriv[:,a,:]*factor*z[:,None]**power
            if power:g[:,col,2]=N[:,a]*factor*power*z**(power-1)*h/2
    J=np.einsum('na,pnb->pab',q.reshape(48,3),g);return phi,J,g@np.linalg.inv(J)
def boundary(q,h,faces=None):
    if faces is None:faces=list(itertools.product(range(3),[-1.,1.]))
    p,w=np.polynomial.legendre.leggauss(9);force=np.zeros((48,3,6));stress=np.zeros((6,3,3))
    for i in range(3):stress[i,i,i]=1
    for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):stress[k,i,j]=stress[k,j,i]=1/np.sqrt(2.)
    for axis,sign in faces:
        points=np.zeros((81,3));points[:,axis]=sign;points[:,[a for a in range(3) if a!=axis]]=np.array(list(itertools.product(p,p)));weight=np.array([a*b for a,b in itertools.product(w,w)]);phi,J,_=data(q,h,points);area=sign*np.cross(J[:,:,(axis+1)%3],J[:,:,(axis+2)%3]);traction=np.einsum('kab,pb->pak',stress,area);force+=np.einsum('p,pn,pak->nak',weight,phi,traction)
    return force.reshape(144,6)
