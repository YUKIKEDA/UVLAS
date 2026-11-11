"""Stop-on-failure local joint moment-space feasibility screening."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time,itertools
import numpy as np
import mq3_joint_moments as method
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    target=ROOT/'rd/structure/main-shell/mq3/records/mq3-local-feasibility-result.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();out=ROOT/'build/analysis'/('mq3-local-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();records=[];stop=None
    try:
        for k,(n,amp),h in itertools.product([0.,.4],[(1,0.),(2,.2)],[.01,.0001]):
            xy,cells,coords,u=method.element.mesh(k,n,amp)
            for cell,(_,_,ids) in enumerate(cells):
                state=method.prepare(coords[ids].reshape(144),h);state['log']['checks']['resources']=time.monotonic()-start<=180 and memory.memory_bytes()<=2*1024**3;artifact=out/f'k{k:g}-n{n}-h{h:g}-c{cell}.npz';np.savez_compressed(artifact,**{key:val for key,val in state.items() if key!='log'});records.append(dict(kappa=k,n=n,amplitude=amp,h=h,cell=cell,**state['log'],artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=digest(artifact)));print(k,n,h,cell,state['log'],flush=True)
                if not all(state['log']['checks'].values()):stop='required_local_gate_failed';break
            if stop:break
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(method.__file__),Path(method.element.__file__),ROOT/'rd/structure/main-shell/mq3/curved-shell-moment-candidate.md',ROOT/'rd/structure/main-shell/mq3/records/mq3-local-feasibility-experiment.md']
    target.write_text(json.dumps(dict(records=records,all_checks_pass=len(records)==20 and stop is None,stop_reason=stop,physical_adoption=False,shared_kernel_evaluated=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
