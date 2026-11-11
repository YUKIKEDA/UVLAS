"""Algebraic probe of a minimal next space change; not a new solver acceptance."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
import numpy as np
from scipy.linalg import svd
import reference_work_thin_limit as t
import reference_work_exact_isometry as exact
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def legendre(n):
    a=[F(1)];b=[F(0),F(1)]
    if n==0:return a
    for k in range(1,n):
        c=[F(0)]*(k+2)
        for i,v in enumerate(b):c[i+1]+=(2*k+1)*v/(k+1)
        for i,v in enumerate(a):c[i]-=k*v/(k+1)
        a,b=b,c
    return b
def main():
    dest=ROOT/'docs/verification/reference-work-missing-moment-result.json';assert not dest.exists()
    p4,p5,p6=legendre(4),legendre(5),legendre(6);b=[(p6[i]-(p4[i] if i<len(p4) else 0))/11 for i in range(7)]
    derivative=[(i+1)*b[i+1] for i in range(6)];endpoints=[sum(v*F(sign)**i for i,v in enumerate(b)) for sign in [-1,1]]
    exact_checks=dict(divergence_corner=derivative==p5,zero_normal_trace=endpoints==[0,0],parent_axis_exchange_pair=True)
    grid=t.m.mesh(.4,2,.2);q=grid['ref'][grid['cells'][0]].reshape(-1);s=t.f.prepare(q,.01);ids=[(i,j,0) for i in range(6) for j in range(6)]
    phi,grad=t.f.basis(s['p'],ids);mom=t.f.reduce(s,lambda p:t.f.basis(p,ids)[0]);extra=np.polynomial.legendre.legval(s['p'][:,0],[0,0,0,0,0,1])*np.polynomial.legendre.legval(s['p'][:,1],[0,0,0,0,0,1])
    newrow=(s['mu']*extra)@phi;current=t.rank_info(svd(mom,compute_uv=False));added=t.rank_info(svd(np.vstack((mom,newrow)),compute_uv=False))
    beta=np.zeros((len(s['p']),3,2));bc=np.array([0.,0.,0.,0.,-1/11,0.,1/11])
    beta[:,0,0]=np.polynomial.legendre.legval(s['p'][:,0],bc)*np.polynomial.legendre.legval(s['p'][:,1],[0,0,0,0,0,1]);beta[:,1,1]=np.polynomial.legendre.legval(s['p'][:,0],[0,0,0,0,0,1])*np.polynomial.legendre.legval(s['p'][:,1],bc)
    A=np.einsum('p,pjr,pjn->rn',s['mu'],s['S'],grad,optimize=True);bubble=np.einsum('p,pjr,pjn->rn',s['mu'],beta,grad,optimize=True)
    corner=ids.index((5,5,0));desired=np.zeros_like(bubble);desired[:,corner]=-s['mu']@(extra*extra)
    projector_current=t.rank_info(svd(A,compute_uv=False));projector_added=t.rank_info(svd(np.vstack((A,bubble)),compute_uv=False))
    x=exact.bilinear(q.reshape(75,3)[[0,4,20,24],0]);delta=x.get((1,1),F(0));missing=-F(2,5)*delta**5/15
    result=dict(exact_checks=exact_checks,current_mid_observation_ranks=current,with_one_internal_moment_ranks=added,current_projected_gradient_ranks=projector_current,with_axis_symmetric_flux_pair_ranks=projector_added,
        bubble_work_error=float(np.linalg.norm(bubble-desired)),bilinear_cross_coefficient=str(delta),missing_quartic_u0x_monomial_coefficient=str(missing),missing_quartic_coefficient_nonzero=bool(missing),
        internal_position_coordinates_added=3,global_trace_coordinates_added=0,new_local_stability_tested=False,nonlinear_tested=False,practical_geometry_mapping_closed=False,solver_adopted=False)
    assert all(exact_checks.values()) and current==[35]*3 and added==[36]*3 and projector_current==[34]*3 and projector_added==[35]*3 and result['bubble_work_error']<=1e-12 and missing
    ps=[Path(__file__),ROOT/'docs/reference-work-missing-moment-proposal.md'];result['sources']={p.relative_to(ROOT).as_posix():digest(p) for p in ps}
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
