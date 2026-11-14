"""Use qualified age kernels; preserve complex loads without float band storage."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.linalg import solve
from uvlm_harmonic_reference import wake_matrix,H,finite,semi
from uvlm_steady_reference import geometry,matrix
def run(m,n,k=.2,eps=1e-5,length=256.,order=12,motion='heave',axis=.375):
    g=geometry(m,n);omega=2*k*20/1.114285714285714;A=matrix(g,eps,0.).astype(complex);W=wake_matrix(g,omega,eps,length,order);A[:,-n:]+=W
    rhs=1j*omega*g['normal'][:,2] if motion=='heave' else -(20+1j*omega*(g['P'][:,0]-axis))*g['normal'][:,2]
    gamma=solve(A,rhs);G=gamma.reshape(m,n);horse=(G-np.vstack([np.zeros((1,n)),G[:-1]])).ravel();left,right=g['left'],g['right'];edgeF=1.2*horse[:,None]*np.cross([20.,0,0],right-left);F=edgeF.sum(axis=0);M=np.cross((left+right)/2-[axis,0,0],edgeF).sum(axis=0);Fc=F.copy();Mc=M.copy();X=g['X']
    for i in range(m):
        for j in range(n):
            q=X[[i,i,i+1,i+1],[j,j+1,j+1,j]]
            for ids in [[0,2,1],[0,3,2]]:
                tri=q[ids];area=np.cross(tri[1]-tri[0],tri[2]-tri[0])/2;f=1.2j*omega*gamma[i*n+j]*area;F+=f;M+=np.cross(tri.mean(axis=0)-[axis,0,0],f)
    return dict(g=g,omega=omega,gamma=gamma,F=F,M=M,circulatory_F=Fc,circulatory_M=Mc,residual=float(np.linalg.norm(A@gamma-rhs)/np.linalg.norm(rhs)),W=W)
