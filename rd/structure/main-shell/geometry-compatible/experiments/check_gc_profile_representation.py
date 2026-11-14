"""Independent float64 checks of the preserved reference-only profile conversion."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.linalg import lstsq
import stable_curvature_profiles as stable
import curvature_work_profiles as raw

ROOT=rd_paths.ROOT

def main():
    records=[]
    g,_=np.polynomial.legendre.leggauss(29)
    for n in [2,4]:
        x=np.concatenate([(j+(g+1)/2)/(2*n) for j in range(2*n)])
        for k in [0.,.4]:
            for kind in [0,1]:
                c=stable.get(n,k,kind)
                A,dA=stable.profiles(x,n,0,kind,[k,0.])
                B,dB=raw.profiles(x,n,0,kind,[k,0.])
                T,_,rank,_=lstsq(A,B,cond=1e-12)
                value=float(np.linalg.norm(A@T-B)/np.linalg.norm(B))
                derivative=float(np.linalg.norm(dA@T-dB)/np.linalg.norm(dB))
                jump=slope=0.
                for j in range(2*n-1):
                    jump=max(jump,float(np.max(abs(np.polynomial.chebyshev.chebval(1.,c[j])-np.polynomial.chebyshev.chebval(-1.,c[j+1])))))
                    if kind==0:
                        dc=np.polynomial.chebyshev.chebder(c,axis=1)
                        slope=max(slope,float(np.max(abs(np.polynomial.chebyshev.chebval(1.,dc[j])-np.polynomial.chebyshev.chebval(-1.,dc[j+1]))))*(4*n))
                C=c.reshape(-1,c.shape[-1]);orth=float(np.max(abs(C.T@C-np.eye(C.shape[1]))))
                stem=ROOT/'build/analysis/stable-curvature-profiles'/f'n{n}-k{k:g}-f{kind}-d40-p100'
                meta=json.loads(Path(str(stem)+'.json').read_text())
                sha=hashlib.sha256(Path(str(stem)+'.npz').read_bytes()).hexdigest()
                checks=dict(span=rank==A.shape[1],value=value<1e-10,derivative=derivative<1e-10,continuity=jump<1e-9 and slope<1e-8,orthogonality=orth<1e-12,hash=sha==meta['artifact_sha256'])
                records.append(dict(n=n,k=k,kind=kind,value=value,derivative=derivative,jump=jump,slope=slope,orthogonality=orth,rank=rank,tail=float(np.max(abs(c[:,-5:]))),checks=checks))
    sources=[Path(__file__),Path(stable.__file__),Path(raw.__file__)]
    sources+=list((ROOT/'build/analysis/stable-curvature-profiles').glob('*'))
    result=dict(passed=all(all(r['checks'].values()) for r in records),records=records,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources if p.is_file()})
    path=ROOT/'rd/structure/main-shell/geometry-compatible/records/geometry-compatible-ancf-profile-audit.json'
    assert not path.exists();path.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result['records'],indent=2));print('PASSED',result['passed'])

if __name__=='__main__':main()
