"""First two load-driven solves; preserve numerical stops without changing gates."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time,gc
import numpy as np
import pg_physical_pde as pde
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT

def main():
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-physical-pde-screen.json';assert not dest.exists();records=[];stops=[];start=time.monotonic();out=ROOT/'build/analysis'/('pg-physical-pde-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    def forbidden(*args,**kwargs):raise RuntimeError('manufactured_strain_called_during_solve')
    for h in [.01,.0001]:
        try:
            tick=time.monotonic();pde.f.old.reference=forbidden;a=pde.prepare(.4,2,0.,h);print('PREPARED',h,a['ns'],a['rcond'],flush=True);s=pde.solve(a,'rational');print('SOLVED',h,flush=True)
            obs,reaction=pde.observe(a,s,'rational');other=pde.external(a,'rational',13,18);load_diff=float(np.linalg.norm(other[0]-s['F'])/np.linalg.norm(s['F']));checks=dict(condition=a['rcond']>=1e-10,linear=s['linear']<=1e-8,equilibrium=max(obs['stationarity'])<=1e-8,load_independent=load_diff<=1e-8,resources=time.monotonic()-tick<=180 and memory.memory_bytes()<=2*1024**3)
            artifact=out/f'n2-h{h:g}.npz';np.savez_compressed(artifact,coef=s['coef'],beta_global=s['beta'],F=s['F'],fixed=s['fixed'],reaction=reaction,Zs=a['Zs'])
            record=dict(h=h,n=2,amplitude=0.,mode='rational',observation=obs,load_difference=load_diff,rcond=a['rcond'],linear=s['linear'],checks={k:bool(v) for k,v in checks.items()},seconds=time.monotonic()-tick,artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(artifact.read_bytes()).hexdigest());records.append(record);print(record,flush=True);del a,s;gc.collect()
            if not all(checks.values()):raise RuntimeError('physical_PDE_numerical_gate')
        except Exception as exc:stops.append(dict(h=h,reason=repr(exc)));print('STOP',stops[-1],flush=True);break
    paths=[Path(__file__),Path(pde.__file__),Path(pde.loads.__file__),ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-physical-pde-protocol.md',ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-candidate-audit.json']
    result=dict(records=records,stops=stops,passed=len(records)==2 and not stops,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),main_shell_adopted=False,sources={q.relative_to(ROOT).as_posix():hashlib.sha256(q.read_bytes()).hexdigest() for q in paths});dest.write_text(json.dumps(result,indent=2),encoding='utf-8')

if __name__=='__main__':main()
