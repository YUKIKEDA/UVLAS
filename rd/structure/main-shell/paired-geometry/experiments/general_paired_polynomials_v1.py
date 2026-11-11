"""Direct reference lifting for arbitrary polynomial graph patches, over Q.

For w''=k, w(0)=w'(0)=0:
 v1=-f I0+integral(f k),
 v2=-x integral(f_y k)+integral(s f_y k), v3=x I0-I1.
This gives eta_xx=eta_xy=0 exactly, without fixing the physical solution to
a membrane-Ritz minimizer. The remaining eta_yy is retained in full energy.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
import paired_exact_polynomials as p
from paired_coarse_polynomials import compose_univariate
import polynomial_reference_surface as surface


def specification(parameters):
    coeffs=surface.coefficients(parameters)
    q=max((sum(k) for k in coeffs),default=0)
    r=max(9,q+8);d=max(2,r+3*max(0,q-1))
    return coeffs,r,d


def power(a,k):
    out={(0,0):F(1)}
    for _ in range(k):out=p.mul(out,a)
    return out


def graph_polynomial(n,cx,cy,coeffs):
    x={(0,0):F(2*cx+1,4*n),(1,0):F(1,4*n)}
    y={(0,0):F(2*cy+1,4*n),(0,1):F(1,4*n)}
    f={}
    for (i,j),a in coeffs.items():f=p.add(f,p.scale(p.mul(power(x,i),power(y,j)),a))
    return x,y,f


@lru_cache(None)
def moments(cell,origin,k,n,maximum):
    if cell<origin:return tuple({} for _ in range(maximum+1))
    seed=origin if cell>origin else cell
    x={(0,0):F(2*seed+1,4*n),(1,0):F(1,4*n)};P=p.legendre(k);out=[]
    for m in range(maximum+1):
        value=p.integrate_t(p.mul(power(x,m),P),n)
        if cell>origin:
            endpoint=sum(value.values());value={(0,0):endpoint} if endpoint else {}
        out.append(value)
    return tuple(out)


def lifting(n,cx,cy,coeffs,origin,k,axis=0,global_seed=False,width=F(1)):
    if axis:
        local=lifting(n,cy,cx,{(j,i):a for (i,j),a in coeffs.items()},origin,k,0,global_seed)
        return [p.scale(p.transpose(local[1]),1/width),p.transpose(local[0]),p.transpose(local[2])]
    x,y,f=graph_polynomial(n,cx,cy,coeffs)
    maximum=max((i for i,j in coeffs),default=0)+1
    if global_seed:
        I=[compose_univariate(a,F(2*cx+1,2*n)-1,F(1,2*n)) for a in moments(0,0,k,F(1,2),maximum)]
    else:I=moments(cx,origin,k,n,maximum)
    J={};Jy={};Jsy={}
    for (i,j),a in coeffs.items():
        J=p.add(J,p.scale(p.mul(power(y,j),I[i]),a))
        if j:
            factor=p.scale(power(y,j-1),a*j)
            Jy=p.add(Jy,p.mul(factor,I[i]));Jsy=p.add(Jsy,p.mul(factor,I[i+1]))
    return [p.scale(p.add(p.scale(p.mul(f,I[0]),-1),J),1/width),p.add(p.scale(p.mul(x,Jy),-1),Jsy),p.add(p.mul(x,I[0]),p.scale(I[1],-1))]


def director_numerator(v,fx,fy,n,width=F(1)):
    one={(0,0):F(1)};xx=p.mul(fx,fx);yy=p.mul(fy,fy);xy=p.mul(fx,fy)
    vx=[p.scale(p.derivative(a,0,n),1/width) for a in v];vy=[p.derivative(a,1,n) for a in v]
    t1=p.add(p.scale(p.mul(fx,vx[0]),-1),p.scale(p.mul(fy,vx[1]),-1),vx[2])
    t2=p.add(p.scale(p.mul(fx,vy[0]),-1),p.scale(p.mul(fy,vy[1]),-1),vy[2])
    return [p.add(p.scale(p.mul(p.add(one,yy),t1),-1),p.mul(xy,t2)),
            p.add(p.mul(xy,t1),p.scale(p.mul(p.add(one,xx),t2),-1)),
            p.scale(p.add(p.mul(fx,t1),p.mul(fy,t2)),-1)]


def cell_raw(n,cx,cy,parameters):
    coeffs,rd,dd=specification(parameters);x,y,f=graph_polynomial(n,cx,cy,coeffs)
    offset,width=surface.chart(parameters);fx=p.scale(p.derivative(f,0,n),1/width);fy=p.derivative(f,1,n);one={(0,0):F(1)}
    g=p.add(one,p.mul(fx,fx),p.mul(fy,fy))
    nr=(n+2)**2+1;nd=(2*n+1)**2+3;count=nr+nd+96*n
    R=[{} for _ in range(count)];D=[{} for _ in R];ER=[{} for _ in R];ED=[{} for _ in R]
    BX=p.bspline_polynomials(n,cx);BY=[p.transpose(a) for a in p.bspline_polynomials(n,cy)]
    for i,a in enumerate(BX):
        for j,b in enumerate(BY):
            index=i*(n+2)+j;R[index]=p.mul(a,b);ER[index]=b if i==0 else {}
    R[nr-1]=f;ER[nr-1]=p.edge(f,cx)
    BX=p.bspline_polynomials(n,cx,True);BY=[p.transpose(a) for a in p.bspline_polynomials(n,cy,True)]
    for i,a in enumerate(BX):
        for j,b in enumerate(BY):
            index=nr+i*(2*n+1)+j;D[index]=p.mul(a,b);ED[index]=b if i==0 else {}
    for j,a in enumerate([p.scale(p.mul(g,fx),-1),p.scale(p.mul(g,fy),-1),g]):
        D[nr+(2*n+1)**2+j]=a;ED[nr+(2*n+1)**2+j]=p.edge(a,cx)
    for axis in [0,1]:
        for origin in range(2*n):
            for k in range(8):
                v=lifting(n,cx,cy,coeffs,origin,k,axis,width=width);d=director_numerator(v,fx,fy,n,width)
                for j in range(3):
                    index=nr+nd+j*(32*n)+axis*(16*n)+origin*8+k
                    R[index]=v[j];D[index]=d[j]
                    if axis:ER[index]=p.edge(v[j],cx);ED[index]=p.edge(d[j],cx)
    return R,D,ER,ED,x


@lru_cache(None)
def trace_indices(n,parameters):
    keepR=set();keepD=set()
    for cy in range(2*n):
        R,D,ER,ED,x=cell_raw(n,0,cy,parameters)
        keepR.update(j for j,a in enumerate(ER) if a);keepD.update(j for j,a in enumerate(ED) if a)
    return tuple(sorted(keepR)),tuple(sorted(keepD))


def cell_dictionary(n,cx,cy,parameters,boundary='clamped'):
    R,D,ER,ED,x=cell_raw(n,cx,cy,parameters);psi=p.add({(0,0):F(1)},p.scale(x,-1))
    if boundary=='clamped':
        return [p.add(a,p.scale(p.mul(psi,b),-1)) for a,b in zip(R,ER)],[p.add(a,p.scale(p.mul(psi,b),-1)) for a,b in zip(D,ED)]
    if boundary!='extended':raise ValueError('boundary')
    ir,id=trace_indices(n,tuple(parameters));zero={}
    return R+[p.mul(psi,ER[j]) for j in ir]+[zero for j in id],D+[zero for j in ir]+[p.mul(psi,ED[j]) for j in id]


def dimensions(n,parameters,boundary):
    coeffs,rd,dd=specification(parameters);raw=(n+2)**2+1+(2*n+1)**2+3+96*n
    if boundary=='extended':raw+=sum(len(a) for a in trace_indices(n,tuple(parameters)))
    return raw,rd,dd


def build(n,parameters,boundary='clamped'):
    columns,rd,degree=dimensions(n,parameters,boundary)
    out=np.zeros((2,2*n,2*n,degree+1,degree+1,columns))
    for cx in range(2*n):
        for cy in range(2*n):
            R,D=cell_dictionary(n,cx,cy,parameters,boundary)
            for field,polys in enumerate([R,D]):
                for j,a in enumerate(polys):
                    for (i,k),v in p.to_legendre(a,degree).items():out[field,cx,cy,i,k,j]=float(v)
    return out


def rational_matrix(n,parameters,rows,columns,boundary='clamped'):
    from flint import fmpq_mat,fmpq
    total,rd,degree=dimensions(n,parameters,boundary)
    rowmap={int(v):i for i,v in enumerate(rows)};colmap={int(v):j for j,v in enumerate(columns)}
    matrix=fmpq_mat(len(rows),len(columns))
    for cx in range(2*n):
        for cy in range(2*n):
            R,D=cell_dictionary(n,cx,cy,parameters,boundary)
            for field,polys in enumerate([R,D]):
                for j in columns:
                    for (a,b),v in p.to_legendre(polys[int(j)],degree).items():
                        index=(((field*(2*n)+cx)*(2*n)+cy)*(degree+1)+a)*(degree+1)+b
                        if index in rowmap:matrix[rowmap[index],colmap[int(j)]]=fmpq(v.numerator,v.denominator)
    return matrix


def coarse(n,parameters,boundary='clamped'):
    coeffs,rd,degree=specification(parameters);out=np.zeros((2,2*n,2*n,degree+1,degree+1,27))
    for cx in range(2*n):
        for cy in range(2*n):
            x,y,f=graph_polynomial(n,cx,cy,coeffs);offset,width=surface.chart(parameters);fx=p.scale(p.derivative(f,0,n),1/width);fy=p.derivative(f,1,n)
            for j,a in enumerate([p.add({(0,0):offset},p.scale(x,width)),y,f]):
                if boundary=='clamped':a=p.clamp(a,x,cx)
                for (i,k),v in p.to_legendre(a,degree).items():out[0,cx,cy,i,k,j]=float(v)
            for k in range(8):
                v=lifting(n,cx,cy,coeffs,0,k,0,True,width);d=director_numerator(v,fx,fy,n,width)
                for field,polys in enumerate([v,d]):
                    for j,a in enumerate(polys):
                        for (i,l),value in p.to_legendre(a,degree).items():out[field,cx,cy,i,l,3+j*8+k]=float(value)
    return out
