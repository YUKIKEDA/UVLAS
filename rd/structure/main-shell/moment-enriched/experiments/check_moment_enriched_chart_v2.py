"""Proper chart exchange (xi,eta,zeta) -> (eta,xi,-zeta), all DOFs."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
import moment_enriched_ancf as f
ROOT=rd_paths.ROOT
def rel(a,b):return float(np.linalg.norm(a-b)/max(1e-30,np.linalg.norm(b)))
def main():
    out=ROOT/'rd/structure/main-shell/moment-enriched/records/moment-enriched-chart-v2-result.json';assert not out.exists()
    src=ROOT/'rd/structure/main-shell/moment-enriched/records/moment-enriched-local-v2-result.json';parent=json.loads(src.read_text());records=[];stops=[]
    for r in parent['records']:
        try:
            d=np.load(ROOT/r['artifact']);q=d['q'].reshape(-1,3);qnew=np.concatenate([q[a:b].reshape(n,n,3).transpose(1,0,2).reshape(-1,3)*sign for a,b,n,sign in [(0,25,5,1),(25,50,5,-1),(50,75,5,1)]]).reshape(-1)
            s=f.prepare(d['q'],r['h']);t=f.prepare(qnew,r['h']);T=np.zeros((t['nscalar'],s['nscalar']))
            for i,(x,y,z) in enumerate(f.VID):T[i,f.VID.index((y,x,z))]=(-1)**z
            for face in t['faces']:
                axis=[1,0,2][face['axis']];sign=face['sign']*(-1 if axis==2 else 1)
                old=next(v for v in s['faces'] if v['axis']==axis and v['sign']==sign)
                p=face['points'][:,[1,0,2]].copy();p[:,2]*=-1;M=f.old.face_values(s,old,p)
                T[face['offset']:face['offset']+face['n'],old['offset']:old['offset']+old['n']]=face['hat'].T@((face['w']/4)[:,None]*M)
            U,force,K,_=f.nonlinear(t,T@d['current'],(T@d['v1'].reshape(-1,3)).reshape(-1),(T@d['v2'].reshape(-1,3)).reshape(-1))
            e=dict(energy=rel(U-d['energy'],d['energy']),force=rel((T.T@force.reshape(-1,3)).reshape(-1)-d['force'],d['force']),tangent=rel(K-d['tangent'],d['tangent']),reference=rel(T@s['X0']-t['X0'],t['X0']),orthogonal=rel(T.T@T-np.eye(len(T)),np.eye(len(T))))
            records.append(dict(geometry=r['geometry'],h=r['h'],metrics=e,passed=max(e.values())<=1e-8))
        except Exception as exc:stops.append(dict(geometry=r['geometry'],h=r['h'],reason=repr(exc)))
    ps=[Path(__file__),Path(f.__file__),src];out.write_text(json.dumps(dict(records=records,stops=stops,passed=len(records)==9 and not stops and all(r['passed'] for r in records),physical_adoption=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in ps}),indent=2),encoding='utf-8')
    print('chart',len(records),'stops',stops,'failed',[r for r in records if not r['passed']],flush=True)
if __name__=='__main__':main()
