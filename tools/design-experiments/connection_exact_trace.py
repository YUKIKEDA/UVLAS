"""Exact flat-PG trace restrictions, rank quotient and verified row scaling."""
from fractions import Fraction as F
import numpy as np
from flint import fmpq,fmpq_mat,arb_mat,ctx
import mass_paired_ancf as shell
def rational(x):return fmpq(x.numerator,x.denominator)
def legendre(x,degree):
    values=[F(1),x]
    for j in range(1,degree):values.append(((2*j+1)*x*values[-1]-j*values[-2])/(j+1))
    return values[:degree+1]
def evaluation(model,points):
    n=model['n'];rows,B,Z,record=shell.source_basis(n,tuple(model['parameters']),'extended');degree=record['d_degree'];k=B.ncols();ns=model['ns'];N=2*n;size=degree+1
    decoded=[]
    for index in rows:
        b=int(index)%size;q=int(index)//size;a=q%size;q//=size;cy=q%N;q//=N;cx=q%N;kind=q//N;decoded.append((kind,cx,cy,a,b))
    W=fmpq_mat(len(points),len(rows));C=fmpq_mat(len(points),ns-k)
    for p,(x,y,zeta) in enumerate(points):
        cx=min(int(x*N),N-1);cy=min(int(y*N),N-1);xx=legendre(2*(x*N-cx)-1,degree);yy=legendre(2*(y*N-cy)-1,degree);z=F(str(model['h']))*zeta/2
        for i,(kind,i0,j0,a,b) in enumerate(decoded):
            if i0==cx and j0==cy:W[p,i]=rational(xx[a]*yy[b]*(1 if kind==0 else z))
        ix=min(int(x*n),n-1);iy=min(int(y*n),n-1);u=x*n-ix;v=y*n-iy
        for a in [0,1]:
            for b in [0,1]:C[p,(ix+a)*(n+1)+iy+b]=rational(z*z/2*(u if a else 1-u)*(v if b else 1-v))
    A=W*B;out=fmpq_mat(len(points),ns)
    for i in range(len(points)):
        for j in range(k):out[i,j]=A[i,j]
        for j in range(ns-k):out[i,k+j]=C[i,j]
    return out
def assemble(models,blocks):
    total=sum(m['ns'] for m in models);nrows=sum(len(maps[0][1]) for maps in blocks);A=fmpq_mat(nrows,total);offset=np.cumsum([0]+[m['ns'] for m in models]);start=0
    for maps in blocks:
        count=len(maps[0][1])
        for index,points,sign in maps:
            E=evaluation(models[index],points)
            for i in range(count):
                for j in range(E.ncols()):A[start+i,int(offset[index])+j]+=sign*E[i,j]
        start+=count
    reduced,rank=A.rref();reduced=fmpq_mat([[reduced[i,j] for j in range(total)] for i in range(rank)])
    ctx.prec=512;T=arb_mat(total,total)
    for index,model in enumerate(models):
        t=model['source_to_normalized'];off=int(offset[index])
        for i in range(t.nrows()):
            for j in range(t.ncols()):T[off+i,off+j]=t[i,j]
    C=arb_mat(reduced)*T;Gram=C*C.transpose();L=shell.cholesky(Gram);W=L.inv();normalized=W*C;check=W*Gram*W.transpose();bound=0.
    for i in range(rank):
        row=0.
        for j in range(rank):
            v=check[i,j]-(1 if i==j else 0);row+=float(abs(v).upper())
        bound=max(bound,row)
    if bound>1e-20:raise ValueError('unverified_constraint_row_preconditioner')
    matrix=np.array([[float(normalized[i,j]) for j in range(total)] for i in range(rank)])
    numeric=np.linalg.svd(matrix,compute_uv=False);ranks=[int(sum(numeric>tol*numeric[0])) for tol in [1e-9,1e-10,1e-11]]
    if ranks!=[rank]*3:raise ValueError('completed_constraint_rank')
    return matrix,dict(raw_rows=nrows,exact_rank=int(rank),completed_ranks=ranks,verified_row_gram_defect=bound)
