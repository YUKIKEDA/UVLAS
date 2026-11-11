"""Finite-h linear energy vs the derived MMS-1 t^3 coefficient."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
import mixed_mean_stress_bending as t
import mixed_mean_stress_element_v2 as e
ROOT=rd_paths.ROOT
def main():
    out=ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-bending-coefficient-result.json';assert not out.exists();records=[];stops=[];start=time.monotonic();rng=np.random.default_rng(20261104)
    for k,a in [(0.,0.),(.4,0.),(.4,.2)]:
        try:
            grid=t.m.mesh(k,2,a);q=grid['ref'][grid['cells'][0]].reshape(-1);cell=t.mixed.local(q,k);op=t.operators(q);roots,_=t.leading_rows(op,k);v=cell['Z']@rng.normal(size=(cell['ny'],2));nc=cell['nc'];q0,q1=v[:nc],v[nc:];leading=np.sum((roots[0]@q0+roots[1]@q1)**2,axis=0);opi=t.operators(q,15,9,True);ri,_=t.leading_rows(opi,k);li=np.sum((ri[0]@q0+ri[1]@q1)**2,axis=0);errors=[]
            for h in [.04,.02,.01,.005]:
                state=e.prepare(q,h,k);finite=np.sum((state['root']@(q0+h*q1))**2,axis=0)/h**3;errors.append((abs(finite-leading)/leading).tolist())
            errors=np.array(errors);orders=np.log2(errors[:-1]/errors[1:]);exact=bool(np.max(errors)<=1e-8);trend=exact or bool(np.max(errors[-1])<=1e-3 and np.min(orders)>=1.8 and np.max(orders)<=2.2);ind=float(np.max(abs(li-leading)/leading));r=dict(kappa=k,amplitude=a,relative_errors=errors.tolist(),observed_orders=orders.tolist(),roundoff_identity=exact,independent_difference=ind,passed=bool(trend and ind<=1e-8));records.append(r);print('COEFFICIENT',r,flush=True)
            if not r['passed']:break
        except Exception as exc:stops.append(dict(kappa=k,amplitude=a,reason=repr(exc)));print('STOP',repr(exc),flush=True);break
    paths=[Path(__file__),Path(t.__file__),Path(e.__file__),ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-bending-protocol.md',ROOT/'rd/structure/main-shell/mixed-mean-stress/mixed-mean-stress-leading-energy.md']
    out.write_text(json.dumps(dict(records=records,stops=stops,all_checks_pass=len(records)==3 and not stops and all(r['passed'] for r in records),seconds=time.monotonic()-start,main_shell_adopted=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
