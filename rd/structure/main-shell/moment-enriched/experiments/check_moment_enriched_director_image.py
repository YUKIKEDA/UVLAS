"""Basic ME-1 normal-gradient image on the full TB-stabilization kernel."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.linalg import svd
import moment_enriched_thin as t
ROOT=rd_paths.ROOT
def main():
    dest=ROOT/'rd/structure/main-shell/moment-enriched/records/moment-enriched-director-image-result.json';assert not dest.exists();records=[];stops=[]
    target=[t.c.LOOK[i,j,k] for i in range(6) for j in range(6) if (i,j)!=(5,5) for k in range(2)]
    for k,amp in [(0.,0.),(.4,0.),(.4,.2)]:
        grid=t.m.mesh(k,2,amp);q=grid['ref'][grid['cells'][0]].reshape(-1)
        for independent in [False,True]:
            try:
                s=t.f.prepare(q,.001,order=15 if independent else 13,zorder=9 if independent else 7,independent=independent)
                J=np.vstack([x['hat'].T@((x['w']/4)[:,None]*x['jump']) for x in s['faces'][:2]]);Z,info=t.nullspace(J,independent=independent)
                G=(s['bodycoeff'][:,2,:]@s['D'])@Z;U,sv,V=svd(G,full_matrices=False,lapack_driver='gesvd' if independent else 'gesdd');ranks=t.rank_info(sv)
                outside=G.copy();outside[target]=0;outside_error=float(np.linalg.norm(outside)/max(1.,np.linalg.norm(G)))
                basis=U[:,:ranks[1]];T=np.eye(len(t.c.IDS))[:,target];span=float(np.linalg.norm(T-basis@(basis.T@T))/np.linalg.norm(T))
                records.append(dict(kappa=k,amplitude=amp,independent=independent,TB_kernel=info,ranks=ranks,outside=outside_error,span=span,passed=ranks==[70]*3 and outside_error<=1e-9 and span<=1e-9))
            except Exception as exc:stops.append(dict(kappa=k,amplitude=amp,independent=independent,reason=repr(exc)))
    paths=[Path(__file__),ROOT/'rd/structure/main-shell/moment-enriched/experiments/moment_enriched_thin.py',ROOT/'rd/structure/main-shell/moment-enriched/experiments/moment_enriched_ancf.py',ROOT/'rd/structure/main-shell/moment-enriched/records/moment-enriched-thin-protocol.md']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_checks_pass=len(records)==6 and not stops and all(x['passed'] for x in records),sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8')
    print('director',records,'stops',stops,flush=True)
if __name__=='__main__':main()
