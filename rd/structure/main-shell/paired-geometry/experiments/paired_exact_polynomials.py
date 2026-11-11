"""Exact rational polynomial coefficients of the particular PG dictionary.

No sampled fit is used to decide algebraic dependencies. Coordinates in each
2n cell are t_x,t_y in [-1,1]. Geometry coefficients are exact decimal inputs.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from scipy.interpolate import BSpline


def add(*polynomials):
    out={}
    for p in polynomials:
        for k,v in p.items():out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}


def scale(p,a):return {k:v*a for k,v in p.items() if v*a}


def mul(p,q):
    out={}
    for (a,b),v in p.items():
        for (c,d),w in q.items():
            k=(a+c,b+d);out[k]=out.get(k,F(0))+v*w
    return {k:v for k,v in out.items() if v}


def derivative(p,axis,n):
    out={}
    for k,v in p.items():
        if k[axis]:
            q=list(k);q[axis]-=1;out[tuple(q)]=v*k[axis]*(4*n)
    return out


def edge(p,cx):
    t=F(-2*cx-1);out={}
    for (a,b),v in p.items():out[(0,b)]=out.get((0,b),F(0))+v*t**a
    return {k:v for k,v in out.items() if v}


def clamp(p,x,cx):return add(p,scale(mul(add({(0,0):F(1)},scale(x,-1)),edge(p,cx)),-1))


@lru_cache(None)
def legendre(k):
    if k==0:return {(0,0):F(1)}
    if k==1:return {(1,0):F(1)}
    return scale(add(scale(mul({(1,0):F(1)},legendre(k-1)),2*k-1),scale(legendre(k-2),1-k)),F(1,k))


def integrate_t(p,n):
    out={(a+1,0):v/F((a+1)*4*n) for (a,b),v in p.items()}
    value=sum(v*(-1)**a for (a,b),v in out.items())
    return add(out,{(0,0):-value})


def primitive(cell,origin,k,n):
    x={(0,0):F(2*cell+1,4*n),(1,0):F(1,4*n)}
    if cell<origin:return {},{}
    if cell==origin:
        P=legendre(k);I0=integrate_t(P,n);I1=integrate_t(mul(x,P),n);I2=integrate_t(mul(mul(x,x),P),n)
    else:
        width=F(1,2*n);center=F(2*origin+1,4*n)
        moments=[[width,0,0],[width*center,width**2/6,0],
                 [width*(center**2+width**2/12),center*width**2/3,width**3/30]]
        I0,I1,I2=[{(0,0):F(row[k])} if k<3 and row[k] else {} for row in moments]
    b=add(mul(x,I0),scale(I1,-1));C=scale(add(mul(mul(x,x),I0),scale(I2,-1)),F(1,2))
    return b,C


def transpose(p):return {(b,a):v for (a,b),v in p.items()}


def bspline_polynomials(n,cell,c0=False):
    knots=[F(0)]*3+[F(j,n) for j in range(1,n) for _ in range(2 if c0 else 1)]+[F(1)]*3
    count=len(knots)-3
    def evaluate(x):
        if x==1:return [F(int(j==count-1)) for j in range(count)]
        values=[F(int(knots[j]<=x<knots[j+1])) for j in range(len(knots)-1)]
        for degree in [1,2]:
            values=[((x-knots[j])*values[j]/(knots[j+degree]-knots[j]) if knots[j+degree]!=knots[j] else F(0))
                    +((knots[j+degree+1]-x)*values[j+1]/(knots[j+degree+1]-knots[j+1]) if knots[j+degree+1]!=knots[j+1] else F(0)) for j in range(len(values)-1)]
        return values
    values=[evaluate(F(2*cell+j,4*n)) for j in range(3)]
    result=[]
    for j in range(count):
        lo,mid,hi=[row[j] for row in values]
        result.append({k:v for k,v in {(0,0):mid,(1,0):(hi-lo)/2,(2,0):(hi+lo)/2-mid}.items() if v})
    return result


def cell_dictionary(n,cx,cy,parameters):
    kx,ky,tau=[F(str(v)) for v in parameters]
    one={(0,0):F(1)};x={(0,0):F(2*cx+1,4*n),(1,0):F(1,4*n)}
    y={(0,0):F(2*cy+1,4*n),(0,1):F(1,4*n)}
    fx=add(scale(x,kx),scale(y,tau));fy=add(scale(y,ky),scale(x,tau))
    xx=mul(fx,fx);yy=mul(fy,fy);xy=mul(fx,fy);g=add(one,xx,yy)
    nr=(n+2)**2;nd=(2*n+1)**2+3;npair=96*n
    R=[{} for _ in range(nr+nd+npair)];D=[{} for _ in R]
    BX=bspline_polynomials(n,cx);BY=[transpose(p) for p in bspline_polynomials(n,cy)]
    for i,a in enumerate(BX):
        for j,b in enumerate(BY):
            p=mul(a,b);boundary=b if i==0 else {}
            R[i*(n+2)+j]=add(p,scale(mul(add(one,scale(x,-1)),boundary),-1))
    BX=bspline_polynomials(n,cx,True);BY=[transpose(p) for p in bspline_polynomials(n,cy,True)]
    for i,a in enumerate(BX):
        for j,b in enumerate(BY):
            p=mul(a,b);boundary=b if i==0 else {}
            D[nr+i*(2*n+1)+j]=add(p,scale(mul(add(one,scale(x,-1)),boundary),-1))
    for j,p in enumerate([scale(mul(g,fx),-1),scale(mul(g,fy),-1),g]):D[nr+(2*n+1)**2+j]=clamp(p,x,cx)
    for axis in [0,1]:
        along,other=(x,y) if axis==0 else (y,x);ka,kb=(kx,ky) if axis==0 else (ky,kx)
        for origin in range(2*n):
            for k in range(8):
                b,C=primitive(cx if axis==0 else cy,origin,k,n)
                if axis:b,C=transpose(b),transpose(C)
                v=[add(scale(C,-ka),scale(mul(other,b),-tau)),
                   add(scale(mul(other,b),-kb),scale(mul(along,b),tau),scale(C,-2*tau)),b]
                if axis:v=v[1],v[0],v[2]
                vx=[derivative(p,0,n) for p in v];vy=[derivative(p,1,n) for p in v]
                t1=add(scale(mul(fx,vx[0]),-1),scale(mul(fy,vx[1]),-1),vx[2])
                t2=add(scale(mul(fx,vy[0]),-1),scale(mul(fy,vy[1]),-1),vy[2])
                num=[add(scale(mul(add(one,yy),t1),-1),mul(xy,t2)),
                     add(mul(xy,t1),scale(mul(add(one,xx),t2),-1)),
                     scale(add(mul(fx,t1),mul(fy,t2)),-1)]
                for j in range(3):
                    index=nr+nd+j*(32*n)+axis*(16*n)+origin*8+k
                    R[index]=v[j] if axis==0 else clamp(v[j],x,cx)
                    D[index]=num[j] if axis==0 else clamp(num[j],x,cx)
    return R,D


@lru_cache(None)
def power_to_legendre(degree):
    matrix={}
    for p in range(degree+1):
        matrix[p]=[]
        for j in range(p+1):
            value=sum(c*F(2,a+p+1) for (a,b),c in legendre(j).items() if (a+p)%2==0)*F(2*j+1,2)
            if value:matrix[p].append((j,value))
    return matrix


def to_legendre(p,degree=13):
    M=power_to_legendre(degree);intermediate={};out={}
    for (a,b),v in p.items():
        for i,c in M[a]:intermediate[(i,b)]=intermediate.get((i,b),F(0))+v*c
    for (i,b),v in intermediate.items():
        if v:
            for j,c in M[b]:out[(i,j)]=out.get((i,j),F(0))+v*c
    return {k:v for k,v in out.items() if v}


def build(n,parameters,basis='legendre'):
    degree=13;columns=(n+2)**2+(2*n+1)**2+3+96*n
    coefficients=np.zeros((2,2*n,2*n,degree+1,degree+1,columns))
    for cx in range(2*n):
        for cy in range(2*n):
            R,D=cell_dictionary(n,cx,cy,parameters)
            for field,polys in enumerate([R,D]):
                for j,p in enumerate(polys):
                    if basis=='legendre':p=to_legendre(p,degree)
                    for (a,b),v in p.items():coefficients[field,cx,cy,a,b,j]=float(v)
    return coefficients


def rational_matrix(n,parameters,rows,columns):
    """Selected coefficient matrix over Q, for exact rank and basis construction."""
    from flint import fmpq_mat,fmpq
    degree=13;nr=(n+2)**2+(2*n+1)**2+3+96*n
    rowmap={int(v):i for i,v in enumerate(rows)};colmap={int(v):j for j,v in enumerate(columns)}
    matrix=fmpq_mat(len(rows),len(columns))
    for cx in range(2*n):
        for cy in range(2*n):
            R,D=cell_dictionary(n,cx,cy,parameters)
            for field,polys in enumerate([R,D]):
                for j in columns:
                    p=to_legendre(polys[int(j)],degree)
                    for (a,b),v in p.items():
                        index=(((field*(2*n)+cx)*(2*n)+cy)*(degree+1)+a)*(degree+1)+b
                        if index in rowmap:matrix[rowmap[index],colmap[int(j)]]=fmpq(v.numerator,v.denominator)
    return matrix
