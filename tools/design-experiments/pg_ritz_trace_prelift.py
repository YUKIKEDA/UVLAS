"""Pre-kernel coordinates; physical trace acceptance applies AFTER exact lifting."""
import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import svd
from pg_trace_nullspace import nullspace
from scipy.linalg import svd,lstsq,qr

def basis(y,segments,degree,continuity):
    multiplicity=degree-continuity;knots=np.r_[np.zeros(degree+1),np.repeat(np.arange(1,segments)/segments,multiplicity),np.ones(degree+1)];count=degree+1+(segments-1)*multiplicity;return BSpline(knots,np.eye(count),degree)(y)

def coefficients(model,shell,x=0.):
    n=model['n'];q,w=np.polynomial.legendre.leggauss(18);y=np.concatenate([(j+(q+1)/2)/(2*n) for j in range(2*n)]);p=np.column_stack((np.full(len(y),x),y));native=shell.native(p,n);R,D,C=[a[0] for a in native];v,d,ds=shell.lift.evaluate(model['lift'],p);r,A,normal,dn,g=shell.lift.geometry.graph(p,model['parameters']);fx=A[:,2,0];fy=A[:,2,1];t1=-fx[:,None]*v[1][:,0]-fy[:,None]*v[1][:,1]+v[1][:,2];t2=-fx[:,None]*v[2][:,0]-fy[:,None]*v[2][:,1]+v[2][:,2];numerator=np.stack((-(1+fy[:,None]**2)*t1+fx[:,None]*fy[:,None]*t2,fx[:,None]*fy[:,None]*t1-(1+fx[:,None]**2)*t2,-fx[:,None]*t1-fy[:,None]*t2),axis=1).reshape(len(y),-1);normal_num=normal*g[:,None]**1.5;rawr=v[0].reshape(len(y),-1);e=model['completion'];T=e['V']/e['s'];ng=model['normal_completion'];G=ng['V']/ng['s'];nd0=D.shape[1];poly_map=np.column_stack((np.eye(nd0),-ng['coeff']@G));weight_map=np.column_stack((np.zeros((3,nd0)),G));ns=model['ns'];nr,nd,nc=model['nr'],model['nd'],model['nc'];offset=nr+nd+nc;fit_error=0.
    def fit(B,V):
        nonlocal fit_error
        a=np.linalg.lstsq(B,V,rcond=1e-13)[0];fit_error=max(fit_error,float(np.linalg.norm(B@a-V)/max(1.,np.linalg.norm(V))));return a
    BR=basis(y,2*n,10,1);cr=fit(BR,R);vr=fit(BR,rawr);Cr=np.zeros((BR.shape[1],ns));Cr[:,:nr]=cr;Cr[:,offset:]=(vr-cr@e['coeff'][:nr])@T
    constant=np.max(abs(g-g[0]))<1e-14
    if constant:
        BD=basis(y,2*n,8,0);pd=fit(BD,D);nw=fit(BD,normal_num/g[0]**1.5);dw=fit(BD,numerator/g[0]**1.5);extended=pd@poly_map+nw@weight_map;Cd=np.zeros((BD.shape[1],ns));Cd[:,nr:nr+nd]=extended;Cd[:,offset:]=(dw-extended@e['coeff'][nr:])@T;blocks=[Cr,Cd]
    else:
        BD=basis(y,n,2,0);BW=basis(y,2*n,13,0);dp=fit(BD,D)@poly_map;dw=fit(BW,normal_num)@weight_map;rawdw=fit(BW,numerator);Cp=np.zeros((BD.shape[1],ns));Cw=np.zeros((BW.shape[1],ns));Cp[:,nr:nr+nd]=dp;Cw[:,nr:nr+nd]=dw;Cp[:,offset:]=-dp@e['coeff'][nr:]@T;Cw[:,offset:]=(rawdw-dw@e['coeff'][nr:])@T;blocks=[Cr,Cp,Cw]
    BC=basis(y,n,1,0);cc=fit(BC,C);Cc=np.zeros((BC.shape[1],ns));Cc[:,nr+nd:offset]=cc;blocks.append(Cc);matrix=np.vstack(blocks);Z,ranks,conditioning=protected_nullspace(model,shell,matrix)
    if fit_error>1e-9:raise ValueError('trace_space_inclusion:'+str(fit_error))
    return Z,ranks,dict(fit_error=fit_error,rows=matrix.shape[0],constant_weight=bool(constant),matrix=matrix,conditioning=conditioning)


def protected_nullspace(model,shell,matrix):
    """Keep analytically homogeneous x-seed and native fields in the trace kernel.

    Coefficient fits are not allowed to turn roundoff in an exact identity into
    an additional essential boundary condition. All other trace coefficients
    remain constrained, without a weakened physical trace tolerance.
    """
    n=model['n'];p,w=shell.lift.geometry.quadrature(n,13)
    F=shell.fields(model,p)[0].transpose(1,0,2).reshape(-1,model['ns'])
    v,d,_=shell.lift.evaluate(model['lift'],p)
    indices=np.arange(32*n,32*n+8);width=3*len(indices)
    raw=np.stack((v[0][:,:,indices].reshape(len(p),-1),d[:,:,indices].reshape(len(p),-1),np.zeros((len(p),width))),axis=1).reshape(3*len(p),width)
    weight=np.repeat(np.sqrt(w),3);A=F*weight[:,None];B=raw*weight[:,None]
    norms=np.linalg.norm(B,axis=0);keep=norms>1e-11*max(norms);B=B[:,keep]/norms[keep]
    q=lstsq(A,B,cond=1e-13,lapack_driver='gelsy')[0]
    inclusion=float(np.linalg.norm(A@q-B)/max(1.,np.linalg.norm(B)))
    edge=np.column_stack((np.zeros(227),np.linspace(0,1,227)))
    native=model['nr']+model['nd']+model['nc'];values=shell.fields(model,edge)[0][:,:,:native].reshape(-1,native)
    _,ss,vh=svd(values,full_matrices=True,check_finite=False)
    rr=[int(sum(ss>tol*ss[0])) for tol in [1e-9,1e-10,1e-11]]
    if len(set(rr))!=1:raise ValueError('protected_native_trace_rank:'+str(rr))
    N=np.zeros((model['ns'],native-rr[0]));N[:native]=vh[rr[0]:].T
    q=np.column_stack((q,N));norms=np.linalg.norm(q,axis=0);q=q/np.maximum(norms,1e-300)
    U,ss,_=svd(q,full_matrices=True,check_finite=False);ranks=[int(sum(ss>tol*ss[0])) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError('protected_homogeneous_rank:'+str(ranks))
    H=U[:,:ranks[0]];W=U[:,ranks[0]:];V,cr,conditioning=nullspace(matrix@W)
    Z=np.column_stack((H,W@V));defect=float(np.linalg.norm(shell.fields(model,edge)[0].reshape(-1,model['ns'])@Z)/max(1.,np.linalg.norm(shell.fields(model,edge)[0])))
    if inclusion>1e-9:raise ValueError('protected_trace_identity:'+str((inclusion,defect)))
    return Z,[model['ns']-Z.shape[1]]*3,dict(remainder=conditioning,known_homogeneous=ranks[0],native_rank=rr,inclusion=inclusion,physical_trace=defect)
