from pathlib import Path
import time,json,hashlib,traceback,numpy as np
from scipy.linalg import solve,eig
from pg_uvlm_free_step import FreeProblem
from run_frequency_delta import peak_bytes
ROOT=Path(__file__).resolve().parents[2]
def main():
    start=time.monotonic();records=[];failures=[]
    try:
        p=FreeProblem();s=p.sstar;x=p.guess(s);R,out,J=p.evaluate(x,s,True);assert np.linalg.norm(R)<1e-10 and np.linalg.norm(out-s)<1e-10;records.append(dict(case='same_free_discrete_equilibrium_and_step',age_Newton_residuals=p.mean_records,full_step_residual=float(np.linalg.norm(R)),fixed_point_difference=float(np.linalg.norm(out-s))))
        rng=np.random.default_rng(5106);v=rng.normal(size=16);v/=np.linalg.norm(v);hh=1e-5;fd=(-p.evaluate(x+2*hh*v,s)[0]+8*p.evaluate(x+hh*v,s)[0]-8*p.evaluate(x-hh*v,s)[0]+p.evaluate(x-2*hh*v,s)[0])/(12*hh);err=np.linalg.norm(fd-J@v)/np.linalg.norm(fd);assert err<1e-6;records.append(dict(case='full_free_step_all_unknown_JVP',relative_error=float(err)))
        d=p.dimension;Rp=np.empty((16,d));Cp=np.empty((d,d));Cx=np.empty((d,16));dh=2e-6
        for j in range(d):
            v=np.eye(d)[j];a,ca,_=p.evaluate(x,s+dh*v);b,cb,_=p.evaluate(x,s-dh*v);Rp[:,j]=(a-b)/(2*dh);Cp[:,j]=(ca-cb)/(2*dh)
        for j in range(16):
            v=np.eye(16)[j];Cx[:,j]=(p.evaluate(x+dh*v,s)[1]-p.evaluate(x-dh*v,s)[1])/(2*dh)
        P=Cp-Cx@solve(J,Rp);v=rng.normal(size=d);v/=np.linalg.norm(v);a=p.advance(s+dh*v)[0];b=p.advance(s-dh*v)[0];fd=(a-b)/(2*dh);err=np.linalg.norm(fd-P@v)/np.linalg.norm(fd);assert err<1e-4;records.append(dict(case='all_free_history_step_linearization',dimension=d,relative_directional_error=float(err)))
        ev,V=eig(P);eligible=[i for i,z in enumerate(ev) if z.imag>.01 and np.linalg.norm(V[:2,i])>1e-5];assert eligible;pick=max(eligible,key=lambda i:abs(ev[i]));z=ev[pick];mode=V[:,pick];amp=1e-4;plus=s+amp*mode.real;minus=s-amp*mode.real
        for j in range(4):plus=p.advance(plus)[0];minus=p.advance(minus)[0]
        response=(plus-minus)/(2*amp);prediction=(z**4*mode).real;err=np.linalg.norm(response-prediction)/np.linalg.norm(prediction);assert err<1e-3;records.append(dict(case='free_deformed_PG_mode_and_time_response',eigenvalue=[float(z.real),float(z.imag)],relative_time_mode_error=float(err),spectral_radius=float(max(abs(ev)))))
    except Exception as err:failures.append(dict(reason=repr(err),traceback=traceback.format_exc()));print(failures[-1],flush=True)
    peak=peak_bytes();seconds=time.monotonic()-start;passed=len(records)==4 and not failures and peak<2*2**30 and seconds<300;paths=[Path(__file__),ROOT/'tools/design-experiments/pg_uvlm_free_step.py',ROOT/'docs/verification/pg-uvlm-free-step-protocol.md',ROOT/'docs/verification/pg-uvlm-coupled-step-result.json'];result=dict(passed=passed,scope='actual two-direction PG, nonlinear deformed equilibrium, all free wake/source couplings and history; not full structural accuracy/practical flutter validation',records=records,failures=failures,seconds=seconds,peak_bytes=peak,sources={q.relative_to(ROOT).as_posix():hashlib.sha256(q.read_bytes()).hexdigest() for q in paths});dest=ROOT/'docs/verification/pg-uvlm-free-step-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2));assert passed
if __name__=='__main__':main()
