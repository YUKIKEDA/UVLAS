"""Same square-root solve, with one Fortran QR buffer for the n8 trace system."""
import numpy as np
from scipy.linalg import solve_triangular
from polynomial_trace_fit import target,h1,make_local,assemble,observe,rigid_graph,KIND_TARGETS,BASE_PERM,INS,BND,m,f

def solve(a):
    A=a['A'];b=a['b'];scale=1/np.linalg.norm(A,axis=0);n=A.shape[1]
    combined=np.empty((A.shape[0],n+2),order='F');np.multiply(A,scale[None,:],out=combined[:,:n]);combined[:,n:]=b
    a['A']=None;del A
    R=f.qrtools.rfactor(combined);del combined
    T=R[:n,:n];rhs=R[:n,n:n+2];y=solve_triangular(T,rhs,check_finite=False);u=scale[:,None]*y
    count=len(a['original']['ref']) if a['baseline'] else len(a['grid']['trace_ref']);trace=np.zeros((count,3,2));trace[a['free']]=u.reshape(-1,3,2)
    local_solutions=[];error=np.zeros(2);gradient=np.zeros((count,3,2))
    for cell in a['local']:
        ni,nlocal=cell['ni'],cell['n'];Rl=cell['R'];s=cell['scale'];ub=trace[cell['trace_ids']].reshape(-1,2);yb=ub/s[ni:,None]
        yi=solve_triangular(Rl[:ni,:ni],Rl[:ni,nlocal:nlocal+2]-Rl[:ni,ni:nlocal]@yb,check_finite=False)
        value=np.vstack((s[:ni,None]*yi,ub));original=np.zeros_like(value);original[cell['perm']]=value
        residual=Rl[:,:nlocal]@(value/s[:,None])-Rl[:,nlocal:nlocal+2];error+=np.sum(residual*residual,axis=0);local_solutions.append(original)
        condensed=cell['tail']@ub-cell['b'];g=cell['tail'].T@condensed;gradient[cell['trace_ids']]+=g.reshape(-1,3,2)
    kkt=np.linalg.norm(scale[:,None]*gradient[a['free']].reshape(-1,2),axis=0)/np.maximum(1e-30,np.sqrt(n)*np.linalg.norm(b,axis=0))
    linear=np.linalg.norm(T@y-rhs,axis=0)/np.maximum(1e-30,np.linalg.norm(rhs,axis=0))
    return dict(trace=trace,local=np.array(local_solutions),fit=dict(error=np.sqrt(error/a['reference']).tolist(),kkt=kkt.tolist(),linear_residual=linear.tolist(),local_rcond=a['min_local'],global_rcond=f.qrtools.rcond(T)))
