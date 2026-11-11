"""Mean-constrained covariant tying; reference fixed oblique projection."""
import itertools
import numpy as np
from scipy.linalg import svd,qr,solve_triangular
import qem_el2_q4 as base

SPACES=[(0,0,3,4),(1,1,4,3),(0,1,3,3),(0,2,3,4),(1,2,4,3),(2,2,8,8)]

def interpolation(nodes,points,independent):
    if independent:
        V=np.polynomial.legendre.legvander(nodes,len(nodes)-1)
        return np.linalg.solve(V.T,np.polynomial.legendre.legvander(points,len(nodes)-1).T).T
    result=np.ones((len(points),len(nodes)))
    for j,x in enumerate(nodes):
        for k,y in enumerate(nodes):
            if j!=k:result[:,j]*=(points-y)/(x-y)
    return result

def tying_matrix(q,plane,J,independent=False):
    order=int(np.sqrt(len(plane)));nodes=np.polynomial.legendre.leggauss(order)[0]
    blocks=[]
    for a,b,ds,dt in SPACES:
        xs=np.polynomial.legendre.leggauss(ds+1)[0];ys=np.polynomial.legendre.leggauss(dt+1)[0]
        ties=np.array([(x,y) for y,x in itertools.product(ys,xs)])
        Is=interpolation(nodes,ties[:,0],independent);It=interpolation(nodes,ties[:,1],independent)
        cardinal=np.einsum('pi,pj->pij',It,Is).reshape(len(ties),len(plane))
        # Interpolate the covariant polynomial, not rational physical strain.
        pull=np.einsum('pi,kij,pj->pk',J[:,:,a],base.MB,J[:,:,b])
        _,Jt,_=base.raw(q,ties,0.,independent)
        covnorm=np.linalg.norm(np.einsum('pi,kij,pj->pk',Jt[:,:,a],base.MB,Jt[:,:,b]),axis=1)
        block=cardinal[:,:,None]*pull[None,:,:]/covnorm[:,None,None]
        blocks.append(block.reshape(len(ties),-1))
    return np.concatenate(blocks)

def construct(state,plane_mu,independent=False):
    D=state['S'].reshape(-1,state['rank']);L=state['tying'];p=len(plane_mu)
    G=(plane_mu[:,None,None]*np.broadcast_to(np.eye(6),(p,6,6))).transpose(1,0,2).reshape(6,6*p)
    F=G@D
    U,s,Vh=svd(F,full_matrices=True,lapack_driver='gesvd' if independent else 'gesdd',check_finite=False)
    if s[-1]<1e-10*s[0]:raise ValueError('physical_mean_rank_failed')
    right=(Vh[:6].T/s)@U.T;Z=Vh[6:].T
    A=L@D@Z
    singular=svd(A,compute_uv=False,lapack_driver='gesvd' if independent else 'gesdd',check_finite=False)
    if singular[-1]<1e-10*singular[0]:raise ValueError('constrained_tying_rank_failed')
    Q,R=qr(A,mode='economic',check_finite=False)
    affine=right@G
    coeff=affine+Z@solve_triangular(R,Q.T@(L-L@D@affine),check_finite=False)
    T=D@coeff
    return T,dict(mean_condition=float(s[0]/s[-1]),tying_condition=float(singular[0]/singular[-1]),
                  idempotence=float(np.linalg.norm(T@T-T)/np.linalg.norm(T)),mean_residual=float(np.linalg.norm(G@T-G)/np.linalg.norm(G)))
