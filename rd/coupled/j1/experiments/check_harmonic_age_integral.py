"""Independent full triangle W-P1 check of the age-derivative representation."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np,json,time,hashlib
from uvlm_harmonic_reference import H,wake_matrix,finite
from uvlm_steady_reference import geometry
from check_vortex_core_load import kernel
import uvlm_analytic_induction_v2 as full
ROOT=rd_paths.ROOT
def main():
    start=time.monotonic();g=geometry(2,4);omega=3.;eps=.02;length=4.;h=length/20/8;W=wake_matrix(g,omega,eps,length,16,h);independent=np.zeros_like(W);te=g['X'][-1];ages=np.linspace(0,length/20,9)
    for r in range(8):
        a=te+[20*ages[r],0,0];b=te+[20*ages[r+1],0,0]
        for j in range(4):
            q=np.array([a[j],a[j+1],b[j+1],b[j]]);mu=np.exp(-1j*omega*ages[[r,r,r+1,r+1]])
            for tri in [[0,2,1],[0,3,2]]:
                for k,(p,n) in enumerate(zip(g['P'],g['normal'])):
                    B,_,_=full.triangle(p,q[tri],eps,False);independent[k,j]+=n@(mu[tri]@B)
    error=float(np.linalg.norm(W-independent)/np.linalg.norm(independent));assert error<1e-8
    p=g['P'];a=te[:-1];b=te[1:];lineerror=float(np.linalg.norm(finite(p,a,b,eps)-kernel(p,a,b,eps))/np.linalg.norm(kernel(p,a,b,eps)));assert lineerror<1e-12
    sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'rd/aero/unsteady-reference/experiments/uvlm_harmonic_reference.py',ROOT/'rd/aero/unsteady-reference/experiments/uvlm_analytic_induction_v2.py']}
    result=dict(passed=True,records=[dict(case='WP1_triangle_vs_age_derivative',relative_error=error),dict(case='complex_capable_finite_segment',relative_error=lineerror)],sources=sources,seconds=time.monotonic()-start)
    (ROOT/'rd/coupled/j1/records/harmonic-age-integral-result.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
