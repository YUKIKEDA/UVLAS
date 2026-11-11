"""Reference-independent polynomial closure, with quadratic cross enrichment.

For a quadratic graph the exact companion numerator has total degree <=14
on a high-degree arm, transverse degree <=5; the native Q2 contribution is
in the conforming degree-six serendipity space. This dictionary therefore
contains all director companions without geometry-small singular directions.
"""
import numpy as np
from scipy.interpolate import BSpline
from numpy.polynomial import legendre as L
import coefficient_compatible_ancf as base
geometry=base.geometry
lift=base.lift
fields=base.fields
kinematics=base.kinematics
compatible_particular=base.compatible_particular


def rraw(p,n,clamped):
    def trial(p):
        native=geometry.trial(p,n)
        count=(n+2)**2
        jets=[a[:,:count] for a in native]
        for axis in [0,1]:
            U=geometry.univariate(p[:,axis],n);other=p[:,1-axis]
            for power in range(3):
                for j,(dx,dy) in enumerate([(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]):
                    along=dx if axis==0 else dy;across=dy if axis==0 else dx
                    factor=other**power if across==0 else power*other**(power-1) if across==1 and power>=1 else power*(power-1)*other**(power-2) if across==2 and power>=2 else np.zeros(len(p))
                    jets[j]=np.column_stack((jets[j],U[along]*factor[:,None]))
        return jets
    jets=trial(p)
    if not clamped:return jets
    edge=p.copy();edge[:,0]=0;boundary=trial(edge);x=p[:,0,None]
    return [jets[0]-(1-x)*boundary[0],jets[1]+boundary[0],jets[2]-(1-x)*boundary[2],jets[3],jets[4]+boundary[2],jets[5]-(1-x)*boundary[5]]


def c0(x,cells,degree):
    knots=np.r_[np.zeros(degree+1),np.repeat(np.arange(1,cells)/cells,degree),np.ones(degree+1)]
    return BSpline(knots,np.eye(cells*degree+1),degree)(x)


def serendipity(p,n):
    x,y=p.T;cx=np.minimum((x*n).astype(int),n-1);cy=np.minimum((y*n).astype(int),n-1)
    tx=x*n-cx;ty=y*n-cy
    X=c0(x,n,1);Y=c0(y,n,1)
    parts=[np.einsum('pi,pj->pij',X,Y).reshape(len(p),-1)]
    for axis in [0,1]:
        cell=cx if axis==0 else cy;t=tx if axis==0 else ty;transverse=Y if axis==0 else X
        polynomial=L.legvander(2*t-1,3)*t[:,None]*(1-t[:,None])
        functions=np.zeros((len(p),n,n+1,4))
        functions[np.arange(len(p)),cell]=transverse[:,:,None]*polynomial[:,None,:]
        parts.append(functions.reshape(len(p),-1))
    interior=[]
    for i,j in [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]:
        value=tx*(1-tx)*ty*(1-ty)*L.legval(2*tx-1,np.eye(3)[i])*L.legval(2*ty-1,np.eye(3)[j])
        f=np.zeros((len(p),n,n));f[np.arange(len(p)),cx,cy]=value;interior.append(f.reshape(len(p),-1))
    return np.column_stack(parts+interior)


def numerator_closure(p,n,clamped):
    def raw(p):
        x,y=p.T;parts=[]
        for axis in [0,1]:
            along=x if axis==0 else y;other=y if axis==0 else x
            Q=L.legvander(2*other-1,5)
            for j in range(6):parts.append(c0(along,2*n,14-j)*Q[:,j,None])
        return np.column_stack(parts+[serendipity(p,n)])
    result=raw(p)
    if clamped:
        edge=p.copy();edge[:,0]=0;result-=(1-p[:,0,None])*raw(edge)
    return result


def prepare(n,parameters,order=16,clamped=True):
    p,w=geometry.quadrature(n,order)
    T,rr,rc=base.independent(rraw(p,n,clamped)[0],np.sqrt(w))
    hx,ht=base.hermite_nodes(n);X,Y=np.meshgrid(hx,hx,indexing='ij')
    dx,dy=np.meshgrid(ht,ht,indexing='ij')
    jet=np.where(dx.ravel(),np.where(dy.ravel(),4,1),np.where(dy.ravel(),2,0))
    hp=np.column_stack((X.ravel(),Y.ravel()));jets=rraw(hp,n,clamped);control=np.zeros_like(jets[0])
    for j in [0,1,2,4]:control[jet==j]=jets[j][jet==j]
    control*=((1/(4*n))**(dx+dy)).ravel()[:,None]
    Rcoeff=(control@T).reshape(len(hx),len(hx),-1)
    raw=numerator_closure(p,n,clamped)
    D,dr,dc=base.independent(raw,np.sqrt(w))
    degree=14;grid=base.nodes(n,degree);X,Y=np.meshgrid(grid,grid,indexing='ij')
    points=np.column_stack((X.ravel(),Y.ravel()))
    coefficients=(numerator_closure(points,n,clamped)@D).reshape(len(grid),len(grid),-1)
    C=base.native(p,n)[2][0];keep=np.arange(C.shape[1])
    if clamped:keep=keep[n+1:]
    CT,cr,cc=base.independent(C[:,keep],np.sqrt(w))
    interp=base.interpolation(p,n,degree,coefficients)[0]
    error=float(np.linalg.norm(interp-raw@D)/np.linalg.norm(raw@D))
    if error>1e-9:raise ValueError('polynomial_closure_interpolation:'+str(error))
    return dict(n=n,parameters=parameters,clamped=clamped,Tr=T,Rcoeff=Rcoeff,Dcoeff=coefficients,degree=degree,
                Cindices=keep,Tc=CT,nr=T.shape[1],nd=D.shape[1],nc=CT.shape[1],ns=T.shape[1]+D.shape[1]+CT.shape[1],
                lift=dict(n=n,parameters=parameters),ranks=dict(r=rr,d=dr,c=cr),rcond=dict(r=rc,d=dc,c=cc),numerator_interpolation=error)
