"""Observe full-position roundtrip defect without altering acceptance thresholds."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
import numpy as np
import mixed_mean_stress_element_v2 as f
ROOT=rd_paths.ROOT
def main():
    out=ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-affine-state-diagnostic.json';assert not out.exists();grid=f.m.mesh(.4,2,0.);q=grid['ref'][grid['cells'][0]].reshape(-1);records=[];F=np.array([[1.003,.002,0.],[0.,1.003,.001],[0.,0.,1.003]])
    for h in [.01,.001,.0001]:
        s=f.prepare(q,h,.4);zero=np.zeros(3*s['nscalar']);a=s['X0']@(F-np.eye(3)).T;current=s['X0']@F.T;b=current-s['X0'];delta=b-a;_,fa,_,_=f.nonlinear_increment(s,a,zero,zero);_,fb,_,_=f.nonlinear_increment(s,b,zero,zero);_,fc,_,_=f.nonlinear_increment(s,a.astype(complex)+1e-25j*delta,zero,zero);_,fw,_,_=f.nonlinear(s,current,zero,zero)
        E=np.einsum('kij,ij->k',f.m.MB,(F.T@F-np.eye(3))/2);P=F@np.einsum('kij,k->ij',f.m.MB,E@f.m.C.T);expected=(f.boundary_map(s)@P.T).reshape(-1);change=fb-fa;prediction=fc.imag/1e-25
        records.append(dict(h=h,position_roundtrip_increment_defect=float(np.linalg.norm(delta)),direct_position_force_error=float(np.linalg.norm(fb-expected)/np.linalg.norm(expected)),preserved_increment_force_error=float(np.linalg.norm(fa-expected)/np.linalg.norm(expected)),force_change_prediction_relative_residual=float(np.linalg.norm(change-prediction)/max(1e-30,np.linalg.norm(change))),wrapper_increment_difference=float(np.linalg.norm(fw-fb))))
    paths=[Path(__file__),Path(f.__file__),ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-increment-contract.md',ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-local-result.json',ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-local-v2-result.json']
    out.write_text(json.dumps(dict(records=records,diagnostic_only=True,original_absolute_state_failure_retained=True,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8');print(records)
if __name__=='__main__':main()
