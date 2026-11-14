"""Reproduce recovery checks with the missing peak-memory observation."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import check_moment_position_recovery as p
import check_q9_plate_convergence as memory
R=p.ROOT
def main():
    out=R/'rd/structure/main-shell/moment-enriched/records/moment-recovery-resource-result.json';assert not out.exists();records=[];stops=[]
    for name,k,n,amp in [('flat',0.,2,0.),('curved-distorted',.4,2,.2),('quartic-reference',None,1,0.)]:
        try:
            if k is None:qs=[p.local.prior.generic()]
            else:
                grid=p.f.m.mesh(k,n,amp);qs=[grid['ref'][ids].reshape(-1) for ids in grid['cells']]
            rec=p.run(name,n,qs,.001);peak=memory.memory_bytes();records.append(dict(name=name,seconds=rec['seconds'],peak_bytes=peak,checks_reproduced=rec['passed'],passed=rec['passed'] and rec['seconds']<=900 and peak<=2*1024**3));print(records[-1],flush=True)
        except Exception as exc:stops.append(dict(name=name,reason=repr(exc)))
    paths=[Path(__file__),Path(p.__file__),Path(p.r.__file__),R/'rd/structure/main-shell/moment-enriched/records/moment-position-recovery-result.json']
    out.write_text(json.dumps(dict(records=records,stops=stops,passed=len(records)==3 and not stops and all(x['passed'] for x in records),sources={x.relative_to(R).as_posix():hashlib.sha256(x.read_bytes()).hexdigest() for x in paths}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
