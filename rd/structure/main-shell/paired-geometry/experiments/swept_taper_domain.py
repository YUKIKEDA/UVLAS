"""Nonrectangular footprint; split at every native cell/edge intersection."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
def lower(s):return .3*s
def upper(s):return 1-.1*s
def quadrature(n,order):
    q,w=np.polynomial.legendre.leggauss(order);points=[];weights=[];N=2*n
    for i in range(N):
        a,b=i/N,(i+1)/N;cuts={a,b}
        for j in range(N+1):
            edge=j/N
            for root in [edge/.3,(1-edge)/.1]:
                if a<root<b:cuts.add(root)
        cuts=sorted(cuts)
        for l,r in zip(cuts[:-1],cuts[1:]):
            s=(l+r)/2+(r-l)*q/2;ws=w*(r-l)/2
            for j in range(N):
                lo=np.maximum(j/N,lower(s));hi=np.minimum((j+1)/N,upper(s));active=hi>lo
                if not np.any(active):continue
                t=(lo[active,None]+hi[active,None])/2+(hi[active,None]-lo[active,None])*q/2
                points.append(np.column_stack((np.repeat(s[active],order),t.ravel())))
                weights.append((ws[active,None]*w*(hi[active,None]-lo[active,None])/2).ravel())
    return np.vstack(points),np.concatenate(weights)
