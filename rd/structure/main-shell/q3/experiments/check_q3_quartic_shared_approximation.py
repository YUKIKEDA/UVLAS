"""Selected constrained-by-feasibility best approximation; no load equilibrium."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
import q3_quartic_approximation as method
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    target=ROOT/'rd/structure/main-shell/q3/records/q3-quartic-shared-approximation-result.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();parent=ROOT/'rd/structure/main-shell/q3/records/developable-general-reference-v3-result.json';prior=json.loads(parent.read_text());assert prior['all_checks_pass'];references={r['h']:2*float(r['decimal_energy']) for r in prior['records'] if r['kind']=='quartic' and r['kappa']==.4};r32,a32=method.reference_h1(32);rh,area=method.reference_h1(40);assert max(abs(r32/rh-1),abs(a32/area-1))<=1e-12;out=ROOT/'build/analysis'/('q3-quartic-shared-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();records=[];stop=None
    try:
        for n,h in [(4,.01),(4,.0001),(8,.001),(8,.0001)]:
            tick=time.monotonic();a=method.assemble(n,h,np.sqrt(references[h]));u,fit,factors=method.solve(a);obs=method.observe(a,h,u);artifact=out/f'n{n}-h{h:g}.npz';data=dict(u=u,xy=a['xy'],coords=a['coords'],ids=np.array([cell['ids'] for cell in a['local']]),free=a['free'],boundary=a['boundary'],local_R=np.array([cell['R'] for cell in a['local']]),local_scale=np.array([cell['scale'] for cell in a['local']]),normalization=a['normalization'],**factors)
            if n==4:data.update(A=a['A'],b=a['b'])
            np.savez_compressed(artifact,**data);checks=dict(ball=obs['h1']<=2*rh*(1+1e-8),kkt=fit['kkt']<=1e-8,local_rank=a['local_rcond']>=1e-11,global_rank=fit['global_rcond']>=1e-11,linear_solve=fit['linear_residual']<=1e-10,quadrature=abs(obs['error']-fit['error'])<=1e-8,reference=abs(obs['reference_norm']/references[h]-1)<=1e-9,h1=abs(obs['h1']/fit['h1']-1)<=1e-8,resources=memory.memory_bytes()<=4*1024**3 and time.monotonic()-tick<=300 and time.monotonic()-start<=900)
            records.append(dict(n=n,h=h,free_dofs=9*len(a['free']),condensed_dofs=9*len(a['boundary']),elements=len(a['cells']),reference_h1=rh,h1_bound=2*rh,local_rcond=a['local_rcond'],fit=fit,observation=obs,interpolation=a['interpolation'],reference_ratio9=a['reference_ratio'],checks={key:bool(v) for key,v in checks.items()},seconds=time.monotonic()-tick,peak_bytes=memory.memory_bytes(),artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=digest(artifact),artifact_bytes=artifact.stat().st_size));print(n,h,'error',obs['error'],'H1',obs['h1'],'checks',checks,'seconds',time.monotonic()-tick,flush=True)
            del a,u,data,factors
            if not all(checks.values()):stop='required_numerical_or_ball_gate_failed';break
    except Exception as exc:stop=repr(exc)
    gate=None
    if len(records)==4 and stop is None:
        amp=records[1]['observation']['error']/records[0]['observation']['error'];gate=dict(absolute_1pct=all(r['observation']['upper']<=.01 for r in records[2:]),thin_amplification=amp,thin_uniformity=amp<=2.)
    sources=[Path(__file__),Path(method.__file__),Path(method.element.__file__),Path(method.reference.__file__),Path(method.rounding.__file__),parent,ROOT/'rd/structure/main-shell/q3/records/q3-quartic-shared-approximation-experiment.md']
    target.write_text(json.dumps(dict(records=records,all_numerical_checks_pass=len(records)==4 and stop is None,stop_reason=stop,gate=gate,physical_adoption=False,rigorous_optimality_lower_bound=False,reference_h1_32_40_difference=abs(r32/rh-1),sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
