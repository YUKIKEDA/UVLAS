import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import time,json,hashlib,traceback,numpy as np
from scipy.linalg import solve,eig
from pg_uvlm_coupled_step import Problem
from run_frequency_delta import peak_bytes
ROOT=rd_paths.ROOT
def main():
    start=time.monotonic();records=[];failures=[]
    try:
        p=Problem();s=p.sstar;x=s[:4];R,out,J=p.evaluate(x,s,True);assert np.linalg.norm(R)<1e-10 and np.linalg.norm(out-s)<1e-10;records.append(dict(case='deformed_PG_finite_wake_fixed_point',residual=float(np.linalg.norm(R)),step_difference=float(np.linalg.norm(out-s)),polynomial_projection_error=p.projection_error))
        rng=np.random.default_rng(5103);direction=rng.normal(size=4);hh=1e-4;fd=(-p.evaluate(x+2*hh*direction,s)[0]+8*p.evaluate(x+hh*direction,s)[0]-8*p.evaluate(x-hh*direction,s)[0]+p.evaluate(x-2*hh*direction,s)[0])/(12*hh);jerr=np.linalg.norm(fd-J@direction)/np.linalg.norm(fd);assert jerr<1e-6;solved,it,res=p.advance(s,x+np.array([.003,-.002,.01,-.005]));assert np.linalg.norm(solved-s)<1e-8;records.append(dict(case='actual_PG_UVLM_strong_step_and_JVP',relative_JVP_error=float(jerr),iterations=it,residual=res))
        d=p.dimension;Rp=np.empty((4,d));Cp=np.empty((d,d));Cx=np.empty((d,4));dh=2e-6
        for j in range(d):
            v=np.eye(d)[j];a,ca,_=p.evaluate(x,s+dh*v);b,cb,_=p.evaluate(x,s-dh*v);Rp[:,j]=(a-b)/(2*dh);Cp[:,j]=(ca-cb)/(2*dh)
        for j in range(4):
            v=np.eye(4)[j];a=p.evaluate(x+dh*v,s)[1];b=p.evaluate(x-dh*v,s)[1];Cx[:,j]=(a-b)/(2*dh)
        P=Cp-Cx@solve(J,Rp);v=rng.normal(size=d);v/=np.linalg.norm(v);a=p.advance(s+dh*v)[0];b=p.advance(s-dh*v)[0];fd=(a-b)/(2*dh);err=np.linalg.norm(fd-P@v)/np.linalg.norm(fd);assert err<1e-4;records.append(dict(case='full_history_implicit_step_linearization',dimension=d,relative_directional_error=float(err)))
        ev,V=eig(P);eligible=[i for i,z in enumerate(ev) if z.imag>.01 and np.linalg.norm(V[:2,i])>1e-5];assert eligible;pick=max(eligible,key=lambda i:abs(ev[i]));z=ev[pick];mode=V[:,pick];amp=1e-4;plus=s+amp*mode.real;minus=s-amp*mode.real
        for j in range(4):plus=p.advance(plus)[0];minus=p.advance(minus)[0]
        response=(plus-minus)/(2*amp);prediction=(z**4*mode).real;terr=np.linalg.norm(response-prediction)/np.linalg.norm(prediction);assert terr<1e-3;records.append(dict(case='deformed_PG_mode_and_nonlinear_time_response',eigenvalue=[float(z.real),float(z.imag)],relative_time_mode_error=float(terr),spectral_radius=float(max(abs(ev)))))
    except Exception as err:failures.append(dict(reason=repr(err),traceback=traceback.format_exc()));print(failures[-1],flush=True)
    peak=peak_bytes();seconds=time.monotonic()-start;passed=len(records)==4 and not failures and peak<2*2**30 and seconds<300;paths=[Path(__file__),ROOT/'rd/coupled/j1/experiments/pg_uvlm_coupled_step.py',ROOT/'rd/coupled/j1/records/pg-uvlm-coupled-step-protocol.md',ROOT/'docs/specs/coupled/coupled-time-stability-specification.md'];result=dict(passed=passed,scope='two actual PG Ritz directions, nonlinear deformed equilibrium, prescribed finite wake, full step history; not full-shell accuracy/free-wake adoption',records=records,failures=failures,seconds=seconds,peak_bytes=peak,sources={q.relative_to(ROOT).as_posix():hashlib.sha256(q.read_bytes()).hexdigest() for q in paths});dest=ROOT/'rd/coupled/j1/records/pg-uvlm-coupled-step-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2));assert passed
if __name__=='__main__':main()
