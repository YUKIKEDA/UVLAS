"""Pure shear from target-field support before constraint row whitening."""
import numpy as np
from scipy.linalg import cholesky, solve, solve_triangular, svdvals
import mq3_interior_stationarity as original


def prepare(state):
    family=original.method.family
    d=original.data(state,9)
    H=d['H']
    L=cholesky(H,lower=True,check_finite=False)
    K=solve_triangular(L.T,np.eye(len(L)),check_finite=False)
    points,field,weight=family.volume(state['q'],0.)
    selected=np.flatnonzero(np.linalg.norm(state['T'],axis=1)>0.)
    nM=np.linalg.norm(field['M'].reshape(-1,117),axis=0)
    nS=np.linalg.norm(field['S'].reshape(-1,117),axis=0)
    pure=np.array([j for j in selected if nM[j]==0. and nS[j]>0.])
    remaining=np.array([j for j in selected if nS[j]==0. and nM[j]>0.])
    closed=len(pure)+len(remaining)==len(selected)==len(H)
    E=np.eye(117)[selected][:,pure]
    P=solve(state['T'][selected],E,check_finite=False)
    Hs=P.T@H@P
    Ls=cholesky(Hs,lower=True,check_finite=False)
    Ns=solve_triangular(Ls.T,np.eye(len(Ls)),check_finite=False)
    F,g=family.moments(state['q'],0.,field,weight)
    Dpure=F[:,pure]@Ns
    Z,ranks,s=original.null(Dpure)
    Qold=P@Ns@Z
    Zs=solve(K,Qold,check_finite=False)
    metrics=dict(pure_target_dimension=len(pure),constraint_ranks=ranks,
        directions=Zs.shape[1],basis_closure=closed,
        constraint_residual=float(np.linalg.norm(Dpure@Z)/max(1.,np.linalg.norm(Dpure))),
        orthogonality=float(np.linalg.norm(Zs.T@Zs-np.eye(Zs.shape[1]))),
        field_purity=float(np.linalg.norm(field['M']@state['T']@Qold)/max(1.,np.linalg.norm(field['M']@state['T']))))
    checks=dict(closure=closed,rank=ranks==[ranks[1]]*3,
        residual=metrics['constraint_residual']<=1e-10,
        orthogonality=metrics['orthogonality']<=1e-10,purity=metrics['field_purity']<=1e-10)
    return dict(K=K,Zs=Zs,checks={k:bool(v) for k,v in checks.items()},metrics=metrics,singular=s)
