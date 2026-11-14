"""Exact global P7 curvature fields, used as protected reference-only coordinates."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from fractions import Fraction as F
from math import comb
import numpy as np
import paired_exact_polynomials as p


def compose_univariate(poly,offset,slope):
    out={}
    for (a,b),v in poly.items():
        for i in range(a+1):out[(i,0)]=out.get((i,0),F(0))+v*comb(a,i)*offset**(a-i)*slope**i
    return {k:v for k,v in out.items() if v}


def build(n,parameters):
    result=np.zeros((2,2*n,2*n,14,14,24));kx,ky,tau=[F(str(v)) for v in parameters]
    global_fields=[p.primitive(0,0,k,F(1,2)) for k in range(8)]
    for cx in range(2*n):
        for cy in range(2*n):
            one={(0,0):F(1)};x={(0,0):F(2*cx+1,4*n),(1,0):F(1,4*n)};y={(0,0):F(2*cy+1,4*n),(0,1):F(1,4*n)}
            fx=p.add(p.scale(x,kx),p.scale(y,tau));fy=p.add(p.scale(y,ky),p.scale(x,tau))
            xx=p.mul(fx,fx);yy=p.mul(fy,fy);xy=p.mul(fx,fy)
            for k,(bg,Cg) in enumerate(global_fields):
                b=compose_univariate(bg,F(2*cx+1,2*n)-1,F(1,2*n));C=compose_univariate(Cg,F(2*cx+1,2*n)-1,F(1,2*n))
                v=[p.add(p.scale(C,-kx),p.scale(p.mul(y,b),-tau)),p.add(p.scale(p.mul(y,b),-ky),p.scale(p.mul(x,b),tau),p.scale(C,-2*tau)),b]
                vx=[p.derivative(a,0,n) for a in v];vy=[p.derivative(a,1,n) for a in v]
                t1=p.add(p.scale(p.mul(fx,vx[0]),-1),p.scale(p.mul(fy,vx[1]),-1),vx[2])
                t2=p.add(p.scale(p.mul(fx,vy[0]),-1),p.scale(p.mul(fy,vy[1]),-1),vy[2])
                d=[p.add(p.scale(p.mul(p.add(one,yy),t1),-1),p.mul(xy,t2)),p.add(p.mul(xy,t1),p.scale(p.mul(p.add(one,xx),t2),-1)),p.scale(p.add(p.mul(fx,t1),p.mul(fy,t2)),-1)]
                for field,polys in enumerate([v,d]):
                    for component,a in enumerate(polys):
                        for (i,j),value in p.to_legendre(a).items():result[field,cx,cy,i,j,component*8+k]=float(value)
    return result
