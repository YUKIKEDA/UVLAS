"""Exact y-invariant subspace of the adopted flat PG space; no new element."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from flint import fmpq_mat,arb_mat,ctx
import connection_shell as pg
import connection_exact_trace as exact

def kernel(A):
    R,rank=A.rref();piv=[next(j for j in range(R.ncols()) if R[i,j]!=0) for i in range(rank)];free=[j for j in range(R.ncols()) if j not in piv];S=fmpq_mat(R.ncols(),len(free))
    for j,c in enumerate(free):
        S[c,j]=1
        for i,k in enumerate(piv):S[k,j]=-R[i,c]
    if A*S!=fmpq_mat(A.nrows(),len(free)):raise ValueError('exact_plane_kernel')
    return S

def floats(A):return np.array([[float(A[i,j]) for j in range(A.ncols())] for i in range(A.nrows())])

def plane(model):
    n=model['n'];rows,B,_,record=pg.shell.source_basis(n,tuple(model['parameters']),'extended');dd=record['d_degree'];N=2*n;size=dd+1;ns=model['ns'];k=B.ncols();coeff={}
    for i,index in enumerate(rows):
        b=int(index)%size;q=int(index)//size;a=q%size;q//=size;cy=q%N;q//=N;cx=q%N;kind=q//N;coeff[kind,cx,cy,a,b]=[B[i,j] for j in range(k)]+[0]*(ns-k)
    H=[];zero=[0]*ns
    for (kind,cx,cy,a,b),v in coeff.items():
        if b:H.append(v)
        elif cy:H.append([x-y for x,y in zip(v,coeff.get((kind,cx,0,a,0),zero))])
    # Missing coefficient rows also impose equality to cell cy=0.
    for kind,cx,cy,a,b in list(coeff):
        if b==0 and cy==0:
            for j in range(1,N):
                if (kind,cx,j,a,0) not in coeff:H.append(coeff[kind,cx,0,a,0])
    for ix in range(n+1):
        for iy in range(1,n+1):
            v=[0]*ns;v[k+ix*(n+1)+iy]=1;v[k+ix*(n+1)]=-1;H.append(v)
    S=kernel(fmpq_mat(H))
    # Exact closure under averaging in y proves that the plane-invariant solve
    # is stationary in the FULL 3D PG space for these nu=0, uniform-y cases.
    # Jensen's inequality then identifies the same unique clamped minimizer.
    keys=set(coeff)
    for kind,cx,cy,a,b in list(coeff):
        for j in range(N):keys.add((kind,cx,j,a,0))
    keys=sorted(keys);raw=fmpq_mat(len(keys)+ns-k,ns);averaged=fmpq_mat(raw.nrows(),ns)
    for i,key in enumerate(keys):
        kind,cx,cy,a,b=key;values=coeff.get(key,zero)
        for j in range(ns):raw[i,j]=values[j]
        if b==0:
            for j in range(ns):averaged[i,j]=sum((coeff.get((kind,cx,l,a,0),zero)[j] for l in range(N)),pg.shell.fmpq(0))/N
    for ix in range(n+1):
        for iy in range(n+1):
            row=len(keys)+ix*(n+1)+iy;raw[row,k+ix*(n+1)+iy]=1
            for j in range(n+1):averaged[row,k+ix*(n+1)+j]=pg.shell.fmpq(1 if j in [0,n] else 2,2*n)
    rawplane=raw*S;rr,rrank=rawplane.transpose().rref();piv=[next(j for j in range(rr.ncols()) if rr[i,j]!=0) for i in range(rrank)];D=fmpq_mat([[rawplane[i,j] for j in range(S.ncols())] for i in piv]).inv()*fmpq_mat([[averaged[i,j] for j in range(ns)] for i in piv])
    if rawplane*D!=averaged:raise ValueError('full_PG_space_not_closed_under_y_average')
    ctx.prec=512;V=model['source_to_normalized'].inv()*arb_mat(S);Gram=V.transpose()*V;L=pg.shell.cholesky(Gram);A=L.inv().transpose();P=V*A
    check=P.transpose()*P;bound=max(sum(float(abs(check[i,j]-(1 if i==j else 0)).upper()) for j in range(check.ncols())) for i in range(check.nrows()))
    if bound>1e-20:raise ValueError('plane_normalization')
    return dict(S=S,A=A,P=floats(P),orthogonality_bound=bound,dimension=S.ncols(),exact_full_space_average_closure=True)

def constraints(models,planes,blocks):
    offsets=np.cumsum([0]+[p['dimension'] for p in planes]);total=int(offsets[-1]);nr=sum(len(maps[0][1]) for maps in blocks);C=fmpq_mat(nr,total);start=0
    for maps in blocks:
        count=len(maps[0][1])
        for i,points,sign in maps:
            E=exact.evaluation(models[i],points)*planes[i]['S']
            for a in range(count):
                for b in range(E.ncols()):C[start+a,int(offsets[i])+b]+=sign*E[a,b]
        start+=count
    null=kernel(C);ctx.prec=512;A=arb_mat(total,total)
    for i,p in enumerate(planes):
        for a in range(p['dimension']):
            for b in range(p['dimension']):A[int(offsets[i])+a,int(offsets[i])+b]=p['A'][a,b]
    V=A.inv()*arb_mat(null);L=pg.shell.cholesky(V.transpose()*V);Z=V*L.inv().transpose();defect=arb_mat(C)*A*Z
    bound=max(sum(float(abs(defect[i,j]).upper()) for j in range(defect.ncols())) for i in range(defect.nrows()))
    if bound>1e-20:raise ValueError('plane_joint_exact_kernel')
    return floats(Z),dict(raw_rows=nr,exact_rank=total-null.ncols(),plane_dofs=total,free_dofs=null.ncols(),verified_kernel_bound=bound)

def data(model,plane,origin,Q,E=1000.):
    xx,ww=np.polynomial.legendre.leggauss(model['degree']+1);n=model['n'];x=np.concatenate([(i+(xx+1)/2)/(2*n) for i in range(2*n)]);wx=np.tile(ww/(4*n),2*n);xy=np.column_stack((x,np.full(len(x),.5)));zz,wz=np.polynomial.legendre.leggauss(4);pp,N,G,det,X=pg.engine.kinematics_from_fields(model,xy,zz,model['h'],False,pg.shell.fields(model,xy));P=plane['P'];scale=model['physical_scale'];N=N@P;G=np.einsum('pib,ij->pjb',G,P)@Q.T/scale;weights=(wx[:,None]*wz).ravel()*det*scale**3
    if np.max(abs(G[:,:,1]))>1e-8:raise ValueError('plane_invariance_physical_gradient')
    out=dict(G=G,weights=weights,E=E);K=pg.linear_stiffness(out);indices=np.array([3*i+j for i in range(P.shape[1]) for j in [0,2]]);return K[np.ix_(indices,indices)],dict(N=N,G=G,weights=weights,P=P,indices=indices)
