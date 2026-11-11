"""Shared kinematics, independently supplied source kernels and load definitions."""
import numpy as np
from scipy.linalg import lu_factor,lu_solve,solve
def schur_solve(A0,W,rhs,n):
    factor=lu_factor(A0);Y=lu_solve(factor,W.real)+1j*lu_solve(factor,W.imag);G=lu_solve(factor,rhs.real)+1j*lu_solve(factor,rhs.imag);gamma=G-Y@solve(np.eye(n)+Y[-n:],G[-n:]);R=A0@gamma+W@gamma[-n:]-rhs;res=np.linalg.norm(R)/np.linalg.norm(rhs);mode='real_LU_TE_Schur'
    if res>=1e-11:
        A=A0.astype(complex);A[:,-n:]+=W;gamma=solve(A,rhs);res=np.linalg.norm(A@gamma-rhs)/np.linalg.norm(rhs);mode='full_complex_after_residual_rejection'
    if res>=1e-11:raise ValueError('frequency_original_residual_unresolved')
    return gamma,float(res),mode
def observe(g,gamma,omega,kind):
    m,n=g['m'],g['n'];G=gamma.reshape(m,n,-1);horse=(G-np.concatenate([np.zeros((1,n,G.shape[-1])),G[:-1]],axis=0)).reshape(m*n,-1);dl=g['right']-g['left'];edge=1.2*horse[:,None,:]*np.cross([20.,0,0],dl)[:,:,None];F=edge.sum(axis=0);M=np.cross(((g['left']+g['right'])/2-[.375,0,0])[:,:,None],edge,axisa=1,axisb=1,axisc=1).sum(axis=0)
    if kind=='reference':
        f=1.2j*omega*gamma[:,None,:]*g['area'][:,None,None]*g['normal'][:,:,None];F+=f.sum(axis=0);M+=np.cross((g['P']-[.375,0,0])[:,:,None],f,axisa=1,axisb=1,axisc=1).sum(axis=0)
    else:
        X=g['X'];quads=np.stack([X[:-1,:-1],X[:-1,1:],X[1:,1:],X[1:,:-1]],axis=2).reshape(-1,4,3)
        for ids in [[0,2,1],[0,3,2]]:
            tri=quads[:,ids];area=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])/2;f=1.2j*omega*gamma[:,None,:]*area[:,:,None];F+=f.sum(axis=0);M+=np.cross((tri.mean(axis=1)-[.375,0,0])[:,:,None],f,axisa=1,axisb=1,axisc=1).sum(axis=0)
    return np.stack([F[2],M[1]],axis=0).T
