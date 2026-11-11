from pathlib import Path
import json,hashlib,numpy as np
from scipy.integrate import quad_vec
from probe_finite_wing_drag import far_drag
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/finite-wing-2-3'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    te=np.array([[0.,-.8,0.],[.2,-.1,.08],[.4,.4,.1],[.6,.9,0.]]);mu=np.array([.4,.7,.3]);kap=np.r_[0,mu]-np.r_[mu,0];eps=.03;total=0.
    for j in range(3):
        a,b=te[j,1:],te[j+1,1:];dl=b-a
        def fun(t):
            r=a+t*dl-te[:,1:];return np.sum(kap*(r@dl)/(np.sum(r*r,axis=1)+eps*eps))/(2*np.pi)
        total+=mu[j]*quad_vec(fun,0.,1.,epsabs=1e-12,epsrel=1e-12)[0]
    ref=-.6*total;candidate=far_drag(te,mu,eps);error=float(abs(candidate/ref-1));assert error<1e-10
    records=[dict(case='exact_log_Trefftz_vs_independent_line_integral',relative_error=error)];artifacts={}
    for p in sorted(OUT.glob('drag-*.json')):
        r=json.loads(p.read_text());records.append(dict(case=p.stem,result=r));artifacts[p.relative_to(ROOT).as_posix()]=sha(p)
    sources={p.relative_to(ROOT).as_posix():sha(p) for p in [Path(__file__),ROOT/'tools/design-experiments/probe_finite_wing_drag.py',ROOT/'tools/design-experiments/check_drag_quadrature.py',ROOT/'docs/verification/finite-wing-drag-diagnostic-protocol.md',ROOT/'docs/verification/edge-load-integral-result.json']}
    result=dict(passed=False,operator_log_integral_verified=True,near_force_independent_comparison_passed=True,finite_grid_near_far_identity_passed=False,physical_drag_continuum_qualified=False,failures=['near midpoint force and exact regularized far-energy functional differ at finite grid; no common continuum error budget qualified'],records=records,sources=sources,artifacts=artifacts);dest=ROOT/'docs/verification/finite-wing-drag-diagnostic-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(passed=False,log_formula_error=error,physical_drag_continuum_qualified=False),indent=2))
if __name__=='__main__':main()
