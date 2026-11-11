import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np,json,time,hashlib
from uvlm_age_scalar import point_derivative,age_integral
from uvlm_harmonic_reference import H
from uvlm_frequency_fast import geometry,age_nodes,age_integral as previous
ROOT=rd_paths.ROOT
def main():
    rng=np.random.default_rng(82);errors=[]
    for i in range(100):
        a=rng.normal(size=3);b=a+rng.normal(size=3);p=rng.normal(size=3);n=rng.normal(size=3);eps=[0.,.01,.1][i%3];dl=b-a;L=np.linalg.norm(dl);e=dl/L;r=p-a;v=point_derivative(*r,*e,L,eps,*n);pc=p.astype(complex);pc[0]+=1e-20j;ref=H(pc[None],np.array([a,b]),eps)[0,0].imag@n/1e-20;errors.append(float(abs(v-ref)/max(abs(ref),1e-12)))
    assert max(errors)<1e-8
    g=geometry(8,48);omega=8/1.114285714285714;ages,weights=age_nodes(g,omega,128.,12);start=time.monotonic();a=age_integral(g['P'],g['normal'],g['X'][-1],1e-5,ages,weights);elapsed=time.monotonic()-start;start=time.monotonic();b=previous(g['P'],g['normal'],g['X'][-1],1e-5,ages,weights);oldtime=time.monotonic()-start;err=float(np.linalg.norm(a-b)/np.linalg.norm(b));assert err<1e-10
    try:point_derivative(1e200,0.,0.,1.,0.,0.,1.,.01,0.,0.,1.)
    except (OverflowError,ZeroDivisionError):pass
    else:assert not np.isfinite(point_derivative(1e200,0.,0.,1.,0.,0.,1.,.01,0.,0.,1.))
    sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'rd/coupled/j1/experiments/uvlm_age_scalar.py',ROOT/'rd/aero/unsteady-reference/experiments/uvlm_frequency_fast.py']};result=dict(passed=True,max_individual_point_derivative_relative_error=max(errors),field_relative_error=err,scalar_seconds_including_first_compile=elapsed,previous_seconds=oldtime,sources=sources);dest=ROOT/'rd/coupled/j1/records/age-scalar-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
