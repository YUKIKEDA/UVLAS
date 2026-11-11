"""Polynomial graph reference patches with true unit-normal thickness.

Parameters: kx, ky, tau, cubic_x, cubic_y, quartic_xy. The represented surface
is exact input geometry; no pointwise curvature substitution is performed.
"""
from fractions import Fraction as F
from math import comb
import numpy as np
from pg_ritz_geometry import quadrature


def coefficients(parameters):
    values=tuple(parameters);kx,ky,tau,cx,cy,mixed=(values[:6]+(0.,)*6)[:6]
    offset,width=chart(parameters);skew=shear(parameters)
    values={(2,0):F(str(kx))/2,(0,2):F(str(ky))/2,(1,1):F(str(tau)),
            (3,0):F(str(cx))/3,(0,3):F(str(cy))/3,(2,2):F(str(mixed))}
    result={}
    for (i,j),a in values.items():
        for k in range(i+1):
            for l in range(i-k+1):
                result[(k,j+l)]=result.get((k,j+l),F(0))+a*comb(i,k)*comb(i-k,l)*offset**(i-k-l)*width**k*skew**l
    return {k:v for k,v in result.items() if v}


def chart(parameters):
    if len(parameters)>6:
        offset,width=F(str(parameters[6])),F(str(parameters[7]))
        if width<=0:raise ValueError('nonpositive_chart_width')
        return offset,width
    return F(0),F(1)


def shear(parameters):return F(str(parameters[8])) if len(parameters)>8 else F(0)


def values(p,parameters,dx=0,dy=0):
    out=np.zeros(len(p),dtype=np.result_type(p,float));x,y=p.T
    for (a,b),v in coefficients(parameters).items():
        if a<dx or b<dy:continue
        factor=float(v)
        for j in range(dx):factor*=a-j
        for j in range(dy):factor*=b-j
        out+=factor*x**(a-dx)*y**(b-dy)
    return out


def graph(p,parameters):
    x,y=p.T;offset,width=map(float,chart(parameters));skew=float(shear(parameters))
    fx,fy,fxs,fys,fxt,fyt=slope_jets(p,parameters)
    g=1+fx*fx+fy*fy;raw=np.column_stack((-fx,-fy,np.ones(len(p))));normal=raw/np.sqrt(g)[:,None]
    gx=2*(fx*fxs+fy*fys);gy=2*(fx*fxt+fy*fyt)
    rx=np.column_stack((-fxs,-fys,np.zeros(len(p))));ry=np.column_stack((-fxt,-fyt,np.zeros(len(p))))
    nx=rx/np.sqrt(g)[:,None]-raw*gx[:,None]/(2*g[:,None]**1.5)
    ny=ry/np.sqrt(g)[:,None]-raw*gy[:,None]/(2*g[:,None]**1.5)
    r=np.column_stack((offset+width*x+skew*y,y,values(p,parameters)))
    a1=np.column_stack((np.full(len(p),width),np.zeros(len(p)),values(p,parameters,1,0)))
    a2=np.column_stack((np.full(len(p),skew),np.ones(len(p)),values(p,parameters,0,1)))
    return r,np.stack((a1,a2),axis=2),normal,np.stack((nx,ny),axis=2),g


def metric_derivatives(p,parameters):
    fx,fy,fxs,fys,fxt,fyt=slope_jets(p,parameters)
    return [2*(fx*fxs+fy*fys),2*(fx*fxt+fy*fyt)]


def slope_jets(p,parameters):
    width=float(chart(parameters)[1]);skew=float(shear(parameters))
    fx=values(p,parameters,1,0)/width;fy=values(p,parameters,0,1)-skew*fx
    fxs=values(p,parameters,2,0)/width;fxt=values(p,parameters,1,1)/width
    fys=values(p,parameters,1,1)-skew*fxs;fyt=values(p,parameters,0,2)-skew*fxt
    return fx,fy,fxs,fys,fxt,fyt


def reference(p,zeta,h,parameters,variable=False):
    r,A,n,dn,g=graph(p,parameters)
    H=h*(1+.1*p[:,0]+.05*p[:,1]) if variable else np.full(len(p),h)
    Ha=np.tile([.1*h,.05*h],(len(p),1)) if variable else np.zeros((len(p),2))
    z=H*zeta/2;J=np.empty((len(p),3,3),dtype=np.result_type(p,zeta,float))
    J[:,:,:2]=A+z[:,None,None]*dn+n[:,:,None]*(zeta[:,None]*Ha/2)[:,None,:]
    J[:,:,2]=H[:,None]*n/2
    return r+z[:,None]*n,J,np.linalg.det(J),H,Ha
