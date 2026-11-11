from pathlib import Path
import numpy as np,json,time,hashlib
from trefftz_observations import midpoint
from check_vortex_core_load import reference_segment
ROOT=Path(__file__).resolve().parents[2]
def main():
    start=time.monotonic();rows=[]
    for n in [32,64,128,256,512,1024,2048]:
        y=-np.cos(np.linspace(0,np.pi,n+1));mu=np.sqrt(1-((y[:-1]+y[1:])/2)**2);te=np.column_stack([np.zeros(n+1),y,np.zeros(n+1)]);D=midpoint(te,mu);exact=np.pi*1.2/8;rows.append(dict(case='prescribed_elliptic_circulation',n=n,drag=float(D),exact=float(exact),relative_error=float(abs(D/exact-1))))
    assert rows[-1]['relative_error']<.001 and all(rows[i+1]['relative_error']<rows[i]['relative_error'] for i in range(6))
    te=np.array([[0.,-.8,0.],[.2,-.1,.08],[.4,.4,.1],[.6,.9,0.]]);mu=np.array([.4,.7,.3]);eps=.03;value=midpoint(te,mu,eps=eps);kap=np.r_[0.,mu]-np.r_[mu,0.];total=0.
    for j in range(3):
        p=(te[j]+te[j+1])/2;p[0]=0;velocity=np.zeros(3)
        for node,k in zip(te,kap):
            a=node.copy();b=node.copy();a[0]=-1024.;b[0]=1024.;velocity+=k*reference_segment(p,a,b,eps)
        total+=mu[j]*velocity@np.cross([1.,0,0],te[j+1]-te[j])
    independent=-.6*total;error=float(abs(value/independent-1));assert error<1e-6;rows.append(dict(case='independent_full_3D_filament_limit',relative_error=error))
    angle=.41;R=np.array([[1,0,0],[0,np.cos(angle),-np.sin(angle)],[0,np.sin(angle),np.cos(angle)]]);cov=[]
    for lam in [.001,1.,1000.]:
        v=midpoint(lam*(te@R.T)+[12.,-.5,.3],mu,eps=lam*eps);cov.append(float(abs(v/value-1)))
    reversal=float(abs(midpoint(te[::-1],mu[::-1],eps=eps)/value-1));assert max(cov+[reversal])<1e-10;rows.append(dict(case='transverse_rotation_scale_translation_reversal',errors=cov,reversal=reversal))
    sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'tools/design-experiments/trefftz_observations.py',ROOT/'docs/verification/induced-drag-common-limit-protocol.md']};result=dict(passed=True,scope='discrete midpoint far observer and prescribed circulation continuum, not finite-wing accuracy',records=rows,sources=sources,seconds=time.monotonic()-start);dest=ROOT/'docs/verification/trefftz-midpoint-operator-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
