"""Complete trace coefficients in their known polynomial/rational function spaces."""
import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import svd

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
    BC=basis(y,n,1,0);cc=fit(BC,C);Cc=np.zeros((BC.shape[1],ns));Cc[:,nr+nd:offset]=cc;blocks.append(Cc);matrix=np.vstack(blocks);_,s,V=svd(matrix,full_matrices=True,check_finite=False);ranks=[int(sum(s>tol*max(1.,np.linalg.norm(matrix)))) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError('complete_trace_coefficient_rank:'+str(ranks))
    if fit_error>1e-9:raise ValueError('trace_space_inclusion:'+str(fit_error))
    return V[ranks[0]:].T,ranks,dict(fit_error=fit_error,rows=matrix.shape[0],constant_weight=bool(constant),matrix=matrix)
