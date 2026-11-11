"""Per-column acceptance: never hide a weak basis in a matrix norm."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,time,traceback
import numpy as np
import uvlm_analytic_induction_v2 as new
import uvlm_analytic_induction as previous
import uvlm_induction_reference as ref
import check_uvlm_wp1_three_dimensional as old
ROOT=rd_paths.ROOT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(a,b):
    absolute=float(np.linalg.norm(a-b));relative=absolute/max(float(np.linalg.norm(b)),1e-30)
    if not (absolute<=1e-10 or relative<=1e-8):raise AssertionError(dict(absolute=absolute,relative=relative))
    return dict(absolute=absolute,relative=relative)
def main():
    start=time.monotonic();dest=ROOT/'rd/aero/kernel/records/uvlm-induction-column-result.json';assert not dest.exists();records=[];failures=[]
    try:
        values=[];derivatives=[]
        for shape,x in old.SHAPES.items():
            for eps in [0.,.02,.08]:
                for p in [old.POINTS[0],old.POINTS[1]]:
                    v,_=new.patch(p,x,eps)
                    for k,mu in enumerate([[1.,0.],[0.,1.]]):
                        r=old.area_reference(x,mu,p[None,:],eps,96)[0];low=old.area_reference(x,mu,p[None,:],eps,64)[0];check(low,r);values.append(dict(group='patch',shape=shape,basis=k,eps=eps,**check(v[k],r)))
        x=np.array([[0.,0,0],[1.,.1,.07],[.2,1.,-.05]]);n=np.cross(x[1]-x[0],x[2]-x[0]);n/=np.linalg.norm(n)
        for beta,label in [(np.array([.3,.4,.3]),'interior'),(np.array([.5,.5,0.]),'edge'),(np.array([1.,0,0.]),'vertex')]:
            for eps,d in [(0.,.01),(0.,.001),(0.,.00001),(.02,0.),(.08,0.)]:
                p=beta@x+d*n;v,_,_=new.triangle(p,x,eps,d==0);lo=ref.integrate(p,x,eps,48,True);hi=ref.integrate(p,x,eps,64,True)
                for k in range(3):check(lo[k],hi[k]);values.append(dict(group='near',location=label,basis=k,eps=eps,distance=d,**check(v[k],hi[k])))
        for beta,d in [(np.array([.3,.4,.3]),.04),(np.array([.3,.4,.3]),0.),(np.array([1.,0,0.]),0.)]:
            p=beta@x+d*n;_,J,_=new.triangle(p,x,.08,d==0);lo=ref.jacobian(p,x,.08,128);hi=ref.jacobian(p,x,.08,192)
            for k in range(3):
                for j in range(13):check(lo[k,:,j],hi[k,:,j]);derivatives.append(dict(basis=k,variable=j,distance=d,**check(J[k,:,j],hi[k,:,j])))
        assert len(values)==81 and len(derivatives)==117;records.append(dict(case='individual_basis_and_derivative_acceptance',velocity_columns=values,geometric_derivative_columns=derivatives))
        with np.errstate(all='ignore'):
            bad=previous.triangle(np.array([1e100]*3),x,.02);unqualified=not (np.isfinite(bad[0]).all() and np.isfinite(bad[1]).all())
            try:new.triangle(np.array([1e100]*3),x,.02)
            except ValueError as exc:reason=str(exc)
            else:raise AssertionError('nonfinite far derivative accepted')
        assert unqualified and reason=='far_induction_precision_unresolved';records.append(dict(case='far_nonfinite_guard_regression',previous_returned_unqualified_result=unqualified,current_rejection=reason,no_tolerance_change=True))
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    paths=[Path(__file__),Path(new.__file__),Path(previous.__file__),Path(ref.__file__),Path(old.__file__),ROOT/'rd/aero/kernel/records/uvlm-induction-adoption-result.json',ROOT/'rd/aero/kernel/records/uvlm-induction-adoption-v2-result.json'];result=dict(passed=len(records)==2 and not failures,records=records,failures=failures,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():sha(p) for p in paths});dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('UVLM_COLUMN_RESULT',result['passed'],flush=True)
if __name__=='__main__':main()
