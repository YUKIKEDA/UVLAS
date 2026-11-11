"""Same proven symmetry, with a per-RHS residual and no dense complex fallback."""
import sys,numpy as np
from scipy.linalg import lu_factor,lu_solve
import run_high_frequency_final_symmetric as experiment
def bounded_schur(A,W,rhs,n):
    factor=lu_factor(A);Y=lu_solve(factor,W.real)+1j*lu_solve(factor,W.imag);small=lu_factor(np.eye(n)+Y[-n:])
    def inverse(b):
        g=lu_solve(factor,b.real)+1j*lu_solve(factor,b.imag);return g-Y@lu_solve(small,g[-n:])
    gamma=inverse(rhs)
    for i in range(4):
        R=rhs-A@gamma-W@gamma[-n:];errors=np.linalg.norm(R,axis=0)/np.linalg.norm(rhs,axis=0)
        if max(errors)<1e-11:return gamma,float(max(errors)),'bounded_real_LU_TE_Schur_per_rhs'
        if i<3:gamma+=inverse(R)
    raise ValueError(('no_dense_complex_fallback',errors.tolist()))
experiment.schur_solve=bounded_schur
if __name__=='__main__':experiment.main()
