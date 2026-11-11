"""Legendre primitives with exact endpoint moment identities for polynomial seeds."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from numpy.polynomial import legendre as L


def evaluate(x,n,*unused):
    count=16*n;k=np.zeros((len(x),count));I=[np.zeros_like(k) for _ in range(3)]
    width=1/(2*n)
    for cell in range(2*n):
        left=cell*width;center=left+width/2
        t=np.clip(2*(x-left)/width-1,-1,1)
        partial=(x>left)&(x<left+width);full=x>=left+width
        inside=(x>=left)&((x<left+width) if cell<2*n-1 else x<=1)
        for degree in range(8):
            column=8*cell+degree;P=np.eye(8)[degree]
            k[inside,column]=L.legval(t[inside],P)
            power=np.array([1.])
            for moment in range(3):
                integral=L.legint(L.legmul(P,power))*width/2
                I[moment][partial,column]=L.legval(t[partial],integral)-L.legval(-1.,integral)
                # Analytic full-interval moments eliminate artificial long tails.
                values=[[width,0.,0.],
                        [width*center,width**2/6,0.],
                        [width*(center**2+width**2/12),center*width**2/3,width**3/30]]
                I[moment][full,column]=values[moment][degree] if degree<3 else 0.
                power=L.legmul(power,[center,width/2])
    return k,I
