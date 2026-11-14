"""High-precision row-space conditioning of complete trace coefficients.

No tiny nonzero coefficient constraint is deleted by a float64 value Gram.
The prepared coefficient matrix is treated as fixed input; all nonzero rows
are retained to 70-digit relative precision before downconverting a null basis.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from decimal import Decimal as D,localcontext
import numpy as np
from scipy.linalg import qr,svd

def nullspace(matrix):
    rows,cols=matrix.shape
    with localcontext() as ctx:
        ctx.prec=80;A=[[D.from_float(float(x)) for x in row] for row in matrix];permutation=list(range(cols));largest=max(abs(x) for row in A for x in row);rank=0;last=None
        for k in range(min(rows,cols)):
            maximum=D(0);pivot=None
            for i in range(k,rows):
                for j in range(k,cols):
                    if abs(A[i][j])>maximum:maximum=abs(A[i][j]);pivot=(i,j)
            if maximum<=largest*D('1e-70'):break
            i,j=pivot;A[k],A[i]=A[i],A[k]
            if j!=k:
                for row in A:row[k],row[j]=row[j],row[k]
                permutation[k],permutation[j]=permutation[j],permutation[k]
            value=A[k][k];last=value;A[k]=[x/value for x in A[k]]
            for i in range(rows):
                if i==k:continue
                value=A[i][k]
                if value:
                    A[i]=[x-value*y for x,y in zip(A[i],A[k])];A[i][k]=D(0)
            rank+=1
        P=np.array([[float(x) for x in row[rank:]] for row in A[:rank]]);C=np.column_stack((np.eye(rank),P));N=np.vstack((-P,np.eye(cols-rank)));Z=np.empty_like(N);Z[np.array(permutation)]=N;Z=qr(Z,mode='economic',check_finite=False)[0]
    s=svd(C,compute_uv=False,check_finite=False);ranks=[int(sum(s>tol*max(1.,np.linalg.norm(C)))) for tol in [1e-9,1e-10,1e-11]]
    if ranks!=[rank]*3:raise ValueError('conditioned_trace_rank:'+str(ranks))
    residual=float(np.linalg.norm(matrix@Z)/max(1.,np.linalg.norm(matrix)));assert residual<=1e-12
    return Z,ranks,dict(decimal_rank=rank,last_pivot=str(last),relative_original_constraint_residual=residual,conditioned_rcond=float(s[-1]/s[0]))
