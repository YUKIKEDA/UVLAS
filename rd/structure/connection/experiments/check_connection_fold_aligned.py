"""Same folded continuum, with declared joint footprint chart boundaries."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from fractions import Fraction as F
import json,time,traceback,hashlib,argparse
import numpy as np
from scipy.linalg import null_space,block_diag,solve
import connection_shell as pg
import connection_exact_trace as exact
from check_connection_fold import reference,nodes,Q
ROOT=rd_paths.ROOT

def candidate(h,n):
    start=time.monotonic();hf=F(str(h));models=[pg.flat_model(n,h,1-h,.2),pg.flat_model(n,h,h,.2),pg.flat_model(n,h,.5,.2)];origins=[np.zeros(3),np.array([1-h,0.,0.]),np.array([1-h/2,0.,h/2])];rot=[np.eye(3),np.eye(3),Q];ys=nodes(n,models[0]['degree']);zs=[F(-1),F(0),F(1)]
    root=[(F(0),y,z) for y in ys for z in zs];right=[(F(1),y,z) for y in ys for z in zs];left=root
    # The pad's surface coordinate x=(1-zeta_rib)/2 spans its complete face.
    # Split at every polynomial-cell boundary on either face.
    zz=nodes(n,models[0]['degree']);a=[(x,y,F(1)) for x in zz for y in ys];b=[(F(0),y,1-2*x) for x in zz for y in ys]
    C,record=exact.assemble(models,[[(0,root,1)],[(0,right,1),(1,left,-1)],[(1,a,1),(2,b,-1)]]);Z0=null_space(C);Z=np.kron(Z0,np.eye(3));print('ALIGNED_TRACE',n,h,record,Z.shape,flush=True)
    data=[pg.data(m,origin=o,Q=q,order=m['degree']+1) for m,o,q in zip(models,origins,rot)];K=block_diag(*(pg.linear_stiffness(d) for d in data));sizes=[3*m['ns'] for m in models];offsets=np.cumsum([0]+sizes);f=np.zeros(sum(sizes));p,w,area=pg.face(models[2],axis=0,side=1.);N,_,_=pg.point(models[2],p);f[offsets[2]::3]=N.T@(w*np.linalg.norm(area,axis=1)/(h*.2));Kr=Z.T@K@Z;v=solve(Kr,Z.T@f,assume_a='pos');u=Z@v
    gap=0.;rng=np.random.default_rng(744);y=rng.uniform(size=43);z=rng.uniform(-1,1,43)
    for i,j,pa,pb in [(0,1,np.column_stack((np.ones(43),y,z)),np.column_stack((np.zeros(43),y,z))),(1,2,np.column_stack(((1-z)/2,y,np.ones(43))),np.column_stack((np.zeros(43),y,z)))]:
        Ni,_,Xi=pg.point(models[i],pa,origins[i],rot[i]);Nj,_,Xj=pg.point(models[j],pb,origins[j],rot[j]);assert np.max(abs(Xi-Xj))<1e-12;gap=max(gap,float(np.max(abs(Ni@u[offsets[i]:offsets[i+1]].reshape(-1,3)-Nj@u[offsets[j]:offsets[j+1]].reshape(-1,3)))))
    row=dict(h=h,n=n,tip=float(f@u),free_residual=float(np.linalg.norm(Z.T@(K@u-f))/np.linalg.norm(f)),constraint=float(np.linalg.norm(np.kron(C,np.eye(3))@u)),holdout_gap=gap,trace=record,seconds=time.monotonic()-start,dofs=len(u),reduced_dofs=len(v));print('ALIGNED_FOLD',json.dumps(row),flush=True);return row

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--h',type=float,nargs='*',default=[.1,.05]);parser.add_argument('--n',type=int,nargs='*',default=[1,2,4]);parser.add_argument('--publish',action='store_true');args=parser.parse_args();records=[];failures=[];start=time.monotonic()
    try:
        for h in args.h:
            refs=[reference(h,d) for d in [2,4,8]];print('REFERENCE',refs,flush=True);cc=[candidate(h,n) for n in args.n];records.append(dict(h=h,reference=refs,candidates=cc))
            for c in cc:c['tip_error']=abs(c['tip']/refs[-1]['tip']-1)
            if args.publish:
                assert [c['n'] for c in cc]==[1,2,4];assert abs(refs[-1]['tip']/refs[-2]['tip']-1)<.002
                assert cc[-1]['tip_error']<.01 and cc[-1]['tip_error']<cc[0]['tip_error'];assert all(c['free_residual']<1e-8 and c['constraint']<1e-9 and c['holdout_gap']<1e-9 for c in cc)
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    if args.publish:
        paths=[Path(__file__),Path(pg.__file__),Path(exact.__file__),ROOT/'rd/structure/connection/experiments/check_connection_fold.py',ROOT/'rd/structure/connection/records/connection-adoption-protocol.md'];dest=ROOT/'rd/structure/connection/records/connection-fold-aligned-result.json';assert not dest.exists();dest.write_text(json.dumps(dict(passed=not failures,records=records,failures=failures,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8');print('CONNECTION_ALIGNED_FOLD_RESULT',not failures,flush=True)
    else:print(json.dumps(dict(records=records,failures=failures),indent=2),flush=True)
if __name__=='__main__':main()
