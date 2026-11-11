import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from decimal import Decimal,localcontext
import numpy as np
def field(g,X,gamma,P,N,eps):
    if np.any(X[:,:,2]!=0) or np.any(P[:,2]!=0) or np.any(N!=[0,0,1]):raise ValueError('planar_decimal_domain')
    m,n=g['m'],g['n'];G=gamma.reshape(m,n,-1);raw=[]
    for i in range(m+1):
        for j in range(n):raw.append((X[i,j],X[i,j+1],(G[i,j] if i<m else 0)-(G[i-1,j] if i else 0)))
    for i in range(m):
        for j in range(n+1):raw.append((X[i,j],X[i+1,j],(G[i,j-1] if j else 0)-(G[i,j] if j<n else 0)))
    out=np.zeros((len(P),gamma.shape[1]),complex);D=Decimal.from_float
    with localcontext() as ctx:
        ctx.prec=70;core=D(float(eps))**2;pi=Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406');sources=[]
        for a,b,s in raw:
            ax,ay=D(float(a[0])),D(float(a[1]));dx=D(float(b[0]))-ax;dy=D(float(b[1]))-ay;L=(dx*dx+dy*dy).sqrt();sources.append((ax,ay,dx/L,dy/L,L,[D(float(v.real)) for v in s],[D(float(v.imag)) for v in s]))
        for row,p in enumerate(P):
            px,py=D(float(p[0])),D(float(p[1]));rr=[Decimal(0) for _ in range(gamma.shape[1])];ii=rr.copy()
            for ax,ay,ex,ey,L,sr,si in sources:
                rx=px-ax;ry=py-ay;projection=rx*ex+ry*ey;t=projection-L;k2=rx*rx+ry*ry-projection*projection+core;ra=(rx*rx+ry*ry+core).sqrt();bx=rx-L*ex;by=ry-L*ey;rb=(bx*bx+by*by+core).sqrt();v=(ex*ry-ey*rx)*(projection/ra-t/rb)/k2/(4*pi)
                for j in range(len(rr)):rr[j]+=v*sr[j];ii[j]+=v*si[j]
            out[row]=[complex(float(a),float(b)) for a,b in zip(rr,ii)]
    return out
def complex_dot(A,x):
    out=np.empty((len(A),x.shape[1]),complex);D=Decimal.from_float
    with localcontext() as ctx:
        ctx.prec=70
        for i in range(len(A)):
            for j in range(x.shape[1]):
                rr=Decimal(0);ii=Decimal(0)
                for a,b in zip(A[i],x[:,j]):
                    ar,ai,br,bi=[D(float(v)) for v in [a.real,a.imag,b.real,b.imag]];rr+=ar*br-ai*bi;ii+=ar*bi+ai*br
                out[i,j]=complex(float(rr),float(ii))
    return out
