"""Raw primitive, continuity, and native-space inclusion checks, independent of loads."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.linalg import lstsq
import paired_geometry_ancf as paired
import curvature_work_profiles as raw
cp=paired.cp;f=paired.f;ROOT=rd_paths.ROOT

def conversion(n,k):
    g,_=np.polynomial.legendre.leggauss(29);x=np.concatenate([(j+(g+1)/2)/(2*n) for j in range(2*n)])
    A,dA=cp.profiles(x,n,0,0,[k,0.]);curv,I=raw.evaluate(x,n,0,[k,0.]);i0,i1,i2=I
    B=np.stack((np.column_stack((-k*(x[:,None]**2*i0-i2)/2,x[:,None]*i0-i1)),np.column_stack((-i0,-k*x[:,None]*i0))))
    dB=np.stack((np.column_stack((-k*x[:,None]*i0,i0)),np.column_stack((-curv,-k*(i0+x[:,None]*curv)))))
    T,_,rank,_=lstsq(np.vstack(A),np.vstack(B),cond=1e-12)
    value=float(np.linalg.norm(np.vstack(A)@T-np.vstack(B))/np.linalg.norm(B));derivative=float(np.linalg.norm(np.vstack(dA)@T-np.vstack(dB))/np.linalg.norm(dB))
    c=cp.get(n,k,0);C=c.reshape(-1,c.shape[-1]);orth=float(np.max(abs(C.T@C-np.eye(C.shape[1]))));jump=slope=0.
    for j in range(2*n-1):
        for part in [0,1]:jump=max(jump,float(np.max(abs(np.polynomial.chebyshev.chebval(1.,c[part,j])-np.polynomial.chebyshev.chebval(-1.,c[part,j+1])))))
        dc=np.polynomial.chebyshev.chebder(c[0],axis=1);slope=max(slope,float(np.max(abs(np.polynomial.chebyshev.chebval(1.,dc[j])-np.polynomial.chebyshev.chebval(-1.,dc[j+1]))))*(4*n))
    checks=dict(span=rank==A.shape[-1],value=value<1e-10,derivative=derivative<1e-10,continuity=jump<1e-9 and slope<1e-8,orthogonality=orth<1e-12)
    return T,dict(n=n,k=k,rank=rank,value=value,derivative=derivative,jump=jump,slope=slope,orthogonality=orth,checks=checks)

def main():
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-basis-audit.json';assert not dest.exists();records=[];conversions={};inclusion=[]
    for n in [2,4]:
        for k in [0.,.4]:conversions[n,k],r=conversion(n,k);records.append(r)
    for n in [2,4]:
        for amp in [0.,.2]:
            geo=f.c.prepare_geometry(.4,n,amp,13);r0=paired.base.surface(n,0,np.zeros((1,2)),geo);nr,nd,nc=[v.shape[1] for v,dx,dy in r0];e=geo['enrich'][0];ns=nr+nd+nc+e['rank'];T=conversions[n,.4];m=16*n;beta=np.sin(np.arange(m*3).reshape(m,3)*.37+.19)
            Q=np.zeros((e['coeff'].shape[1],3));Q[:T.shape[0],0]=T[:,:m]@beta[:,0];Q[:T.shape[0],2]=T[:,m:]@beta[:,0]
            # Test three different full bending fields, each with x/z components.
            coeff=np.zeros((ns,3,3))
            for j in range(3):
                Q[:]=0;Q[:T.shape[0],0]=T[:,:m]@beta[:,j];Q[:T.shape[0],2]=T[:,m:]@beta[:,j];coeff[:nr+nd,:,j]=e['coeff']@Q;coeff[nr+nd+nc:,:,j]=e['s'][:,None]*(e['V'].T@Q)
            for h in [.01,.0001]:
                ve=vn=se=sn=0.
                p=np.array([[x,y,z] for x in np.linspace(-.96,.96,7) for y in np.linspace(-.91,.91,5) for z in [-1.,-.27,.43,1.]])
                for cell,q in enumerate(geo['refs']):
                    N,G,B,det,X=f.c.kinematics(n,cell,q,h,p,geo);E,U=f.bending(q,p,h,.4,n);actual=(N@coeff.reshape(ns,9)).reshape(len(p),3,3);target=U@beta;ds=B@coeff.reshape(3*ns,3)-E@beta
                    ve+=np.sum((actual-target)**2);vn+=np.sum(target**2);se+=np.sum(ds**2);sn+=np.sum((E@beta)**2)
                pos=float(np.sqrt(ve/vn));strain=float(np.sqrt(se/sn));inclusion.append(dict(n=n,amplitude=amp,h=h,position=pos,strain=strain,passed=max(pos,strain)<=1e-8));print('INCLUSION',inclusion[-1],flush=True)
    sources=[Path(__file__),Path(paired.__file__),Path(cp.__file__),Path(raw.__file__)]
    for path in (ROOT/'build/analysis/paired-curvature-profiles').glob('*.npz'):
        meta=path.with_suffix('.json');record=json.loads(meta.read_text());assert hashlib.sha256(path.read_bytes()).hexdigest()==record['artifact_sha256'];sources.extend([path,meta])
    passed=all(all(r['checks'].values()) for r in records) and all(r['passed'] for r in inclusion)
    result=dict(passed=passed,representation=records,native_inclusion=inclusion,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources});dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('PASSED',passed,flush=True)

if __name__=='__main__':main()
