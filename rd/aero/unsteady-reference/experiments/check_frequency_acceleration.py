import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np,json,time,hashlib
from uvlm_frequency_fast import segment,bound,run
from uvlm_harmonic_reference import finite
from uvlm_harmonic_reference_v2 import run as previous
from uvlm_steady_reference import geometry,matrix
ROOT=rd_paths.ROOT
def main():
    start=time.monotonic();rng=np.random.default_rng(710);errors=[]
    for _ in range(50):
        a=rng.normal(size=3);b=a+rng.normal(size=3);p=rng.normal(size=3);n=rng.normal(size=3);v,d=segment(p,a,b,.02,n);pc=p.astype(complex);pc[0]+=1e-20j;ref=finite(pc[None],a[None],b[None],.02)[0,0]@n;errors.append(abs(d-ref.imag/1e-20)/max(abs(d),1e-12));assert abs(v-ref.real)<1e-12
    g=geometry(2,8);A=bound(g['P'],g['normal'],g['X'],1e-5);B=matrix(g,1e-5,0);err=float(np.linalg.norm(A-B)/np.linalg.norm(B));assert max(errors)<1e-9 and err<1e-12
    x=run(8,48,length=256.);y=previous(8,48,length=256.);field=float(np.linalg.norm(x['W']-y['W'])/np.linalg.norm(x['W']));load=float(np.linalg.norm(x['F']-y['F'])/np.linalg.norm(x['F']));assert max(field,load)<1e-8
    sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'rd/aero/unsteady-reference/experiments/uvlm_frequency_fast.py',ROOT/'rd/aero/unsteady-reference/experiments/uvlm_harmonic_reference_v2.py']}
    result=dict(passed=True,records=[dict(case='point_derivative_50_directions',max_relative_error=max(errors)),dict(case='bound_matrix',relative_error=err),dict(case='frequency_field_and_load',field_error=field,load_error=load)],sources=sources,seconds=time.monotonic()-start);dest=ROOT/'rd/aero/unsteady-reference/records/frequency-acceleration-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
