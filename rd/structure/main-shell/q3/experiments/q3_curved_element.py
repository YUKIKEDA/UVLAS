"""Q3 x thickness-P2 compatible linear solid-shell research element."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from functools import lru_cache
import itertools
import numpy as np
import check_q9_independent_load as material
import check_q9_distorted_plate as old_grid
NODES=np.linspace(-1.,1.,4)
MB=np.zeros((6,3,3))
for i in range(3):MB[i,i,i]=1.
for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):MB[k,i,j]=MB[k,j,i]=1/np.sqrt(2.)
S=np.diag([1.,1.,1.,np.sqrt(2.),np.sqrt(2.),np.sqrt(2.)]);C=S@material.elastic@S;CH=np.linalg.cholesky(C).T
def cardinal(x):
    x=np.asarray(x);values=[];derivatives=[]
    for i in range(4):
        others=[j for j in range(4) if j!=i];den=np.prod([NODES[i]-NODES[j] for j in others]);v=np.ones_like(x);d=np.zeros_like(x)
        for j in others:v*=x-NODES[j]
        for j in others:
            term=np.ones_like(x)
            for k in others:
                if k!=j:term*=x-NODES[k]
            d+=term
        values.append(v/den);derivatives.append(d/den)
    return np.array(values).T,np.array(derivatives).T
def surface(points):
    ls,ds=cardinal(points[:,0]);lt,dt=cardinal(points[:,1]);N=np.einsum('pi,pj->pij',lt,ls).reshape(-1,16);gs=np.einsum('pi,pj->pij',lt,ds).reshape(-1,16);gt=np.einsum('pi,pj->pij',dt,ls).reshape(-1,16)
    return N,np.stack((gs,gt),axis=2)
@lru_cache(None)
def quadrature(order):
    p,w=np.polynomial.legendre.leggauss(order);return np.array(list(itertools.product(p,p,p))),np.array([a*b*c for a,b,c in itertools.product(w,w,w)])
def tables(points,h):
    N,g=surface(points);z=h*points[:,2]/2;powers=np.column_stack((np.ones(len(z)),z,z*z/2));dz=np.column_stack((np.zeros(len(z)),np.ones(len(z)),z));phi=(N[:,:,None]*powers[:,None,:]).reshape(-1,48);gradient=np.zeros((len(z),48,3));gradient[:,:,:2]=(g[:,:,None,:]*powers[:,None,:,None]).reshape(-1,48,2);gradient[:,:,2]=(N[:,:,None]*dz[:,None,:]).reshape(-1,48)*h/2
    return phi,gradient
def evaluate(q,h,points):
    phi,g=tables(points,h);J=np.einsum('ia,pib->pab',q.reshape(48,3),g)
    if np.any(np.linalg.det(J)<=0):raise ValueError('nonpositive_reference_Jacobian')
    grad=g@np.linalg.inv(J);B=np.einsum('kij,pnj->pkni',MB,grad).reshape(-1,6,144)
    return phi,J,B
def bilinear(points,corners):
    s,t=points[:,:2].T;N=np.column_stack(((1-s)*(1-t),(1+s)*(1-t),(1+s)*(1+t),(1-s)*(1+t)))/4
    ds=np.column_stack((-(1-t),1-t,1+t,-(1+t)))/4;dt=np.column_stack((-(1-s),-(1+s),1+s,1-s))/4
    return N@corners,np.stack((ds@corners,dt@corners),axis=2)
def mesh(kappa,n,amplitude):
    vertices,_,_=old_grid.mesh(n,amplitude);vertices=vertices.copy();vertices[:,:,1]-=.5;side=3*n+1;xy=np.full((side*side,2),np.nan);cells=[]
    for j,i in itertools.product(range(n),range(n)):
        ids=np.array([(3*j+b)*side+3*i+a for b,a in itertools.product(range(4),range(4))]);corners=vertices[[j,j,j+1,j+1],[i,i+1,i+1,i]];points=np.array([[s,t] for t,s in itertools.product(NODES,NODES)]);local,_=bilinear(points,corners)
        for index,value in zip(ids,local):
            if np.all(np.isfinite(xy[index])):assert np.max(abs(xy[index]-value))<1e-14
            xy[index]=value
        cells.append((i,j,ids))
    x,y=xy.T;coords=np.zeros((len(x),3,3));coords[:,0]=np.column_stack((x,y,kappa*x*x/2));coords[:,1]=np.column_stack((-kappa*x,np.zeros(len(x)),np.ones(len(x))));u=np.zeros((len(x),9));u[:,:3]=np.column_stack((-kappa*x**3/3,np.zeros(len(x)),x*x/2));u[:,3:6]=np.column_stack((-x,np.zeros(len(x)),-kappa*x*x));return xy,cells,coords,u
def target(kappa,x,z):
    a=np.column_stack((np.ones_like(x),np.zeros_like(x),kappa*x));j=1+kappa*kappa*x*x-kappa*z;return np.einsum('kij,pij->pk',MB,-z[:,None,None]*a[:,:,None]*a[:,None,:]/j[:,None,None])
