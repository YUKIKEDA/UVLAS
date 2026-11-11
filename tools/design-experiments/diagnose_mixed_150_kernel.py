"""Fixed saved-kernel test and exact bound for a specified bending sector."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
import numpy as np
from scipy.linalg import qr,lstsq,svdvals
import mixed_order_150 as m
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1.,np.linalg.norm(b)))
def main():
    dest=ROOT/'docs/verification/mixed-150-root-diagnostic.json';assert not dest.exists()
    parent=ROOT/'docs/verification/mixed-order-150-result.json';data=json.loads(parent.read_text());records=[]
    for item in [r for r in data['records'] if r['kappa']]:
        p=ROOT/item['artifact'];assert digest(p)==item['artifact_sha256'];saved=np.load(p)
        mesh=m.mesh(item['kappa'],2,item['amplitude']);assert np.array_equal(mesh['ref'],saved['ref']) and np.array_equal(mesh['cells'],saved['cells'])
        a=m.assemble(mesh,item['kappa'],9);b=m.assemble(mesh,item['kappa'],11)
        N=saved['kernel'];scale=saved['scale'];V=N/scale[:,None]
        modes=np.column_stack([m.displacement(mesh,kind,item['kappa'])[a['free_rd']].reshape(-1) for kind in ['quadratic','cubic']])
        Q,_=qr(modes/scale[:,None],mode='economic',check_finite=False)
        contain=rel((np.eye(len(Q))-V@V.T)@Q,Q)
        projector=float(np.linalg.norm(V@V.T-Q@Q.T))
        T=np.column_stack((a['A1']@modes,a['D1']));sc=1/np.linalg.norm(T,axis=0);s=svdvals(T*sc,check_finite=False)
        ranks=[int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]]
        coef,_,rank,_=lstsq(T*sc,a['targets']['quartic'],lapack_driver='gelsy',cond=1e-10,check_finite=False);coef*=sc
        error=float(np.linalg.norm(b['A1']@(modes@coef[:2])+b['D1']@coef[2:]-b['targets']['quartic'])/np.linalg.norm(b['targets']['quartic']))
        u=np.zeros_like(mesh['ref']);u[a['free_rd']]=(modes@coef[:2]).reshape(-1,3);u[a['free_c']]=coef[2:].reshape(-1,3)
        dc=u[np.r_[mesh['free'][1],mesh['free'][2]]];h1=float(np.sqrt(np.sum(dc*(b['M']@dc))))
        checks=dict(modes_in_kernel=contain<=1e-10,restricted_rank=ranks==[T.shape[1]]*3 and rank==T.shape[1],
            sphere=h1<=2*m.reference_h1('quartic',item['kappa'],40)*(1+1e-8))
        if item['amplitude']:
            checks['whole_kernel_two_modes']=N.shape[1]==2 and projector<=1e-7
            checks['restricted_reproduces_full_fit']=abs(error-item['solutions']['quartic']['observed_error'])<=1e-8
        records.append(dict(kappa=item['kappa'],amplitude=item['amplitude'],kernel_dimension=N.shape[1],two_mode_containment=contain,
            projector_difference=projector,restricted_columns=T.shape[1],restricted_ranks=ranks,restricted_error=error,
            full_kernel_error=item['solutions']['quartic']['observed_error'],h1=h1,checks={k:bool(v) for k,v in checks.items()},
            artifact=item['artifact'],artifact_sha256=item['artifact_sha256']))
    k=F(2,5);nu=F(3,10)
    moment=lambda j:F(1,j+1)+k*k/F(j+3)
    A,B,C=moment(0),moment(1),moment(2);r0,r1=moment(2),moment(3);det=A*C-B*B
    alpha=(C*r0-B*r1)/det;beta=(A*r1-B*r0)/det
    residual=moment(4)-r0*alpha-r1*beta
    bound2=(1-(nu/(1-nu))**2)*residual/moment(4)
    exact=dict(kappa=str(k),nu=str(nu),alpha=str(alpha),beta=str(beta),
        relaxed_relative_error_squared=str(bound2),relaxed_relative_error=float(np.sqrt(float(bound2))),
        scope='global_quadratic_cubic_bending_sector_plus_arbitrary_normal_shear_c_field',
        full_discrete_kernel_equivalence='numerical_2x2_distorted_only_not_an_exact_global_kernel_certificate',
        finite_h_lower_bound=False)
    checks=dict(local_sector_bound_positive=bound2>0,
        bound_below_restricted_fits=all(np.sqrt(float(bound2))<=r['restricted_error']+1e-10 for r in records))
    sources=[Path(__file__),Path(m.__file__),parent,ROOT/'docs/verification/mixed-150-root-protocol.md',ROOT/'docs/verification/mixed-150-root-preparation.json']
    result=dict(records=records,exact_sector_bound=exact,checks={k:bool(v) for k,v in checks.items()},
        all_checks_pass=all(all(r['checks'].values()) for r in records) and all(checks.values()),
        candidate_changed=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
