"""Thickness refinement of folded joint using existing PG child patches."""
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.linalg import block_diag,solve
import time,json,argparse,traceback,hashlib
import connection_shell as pg
import connection_plane_subspace as plane
from check_connection_fold import nodes,Q,reference
ROOT=Path(__file__).resolve().parents[2]

def candidate(h,n,layers):
    start=time.monotonic();hf=F(str(h));models=[];origins=[];rot=[];blocks=[];ys=[F(1,3)];zz=[F(-1),F(0),F(1)];xs=nodes(n,9)
    for part in [0,1,2]:
        for j in range(layers):
            if part<2:
                length=1-h if part==0 else h;origin=np.array([0 if part==0 else 1-h,0.,-h/2+(j+.5)*h/layers]);q=np.eye(3)
            else:length=.5;origin=np.array([1-h+(j+.5)*h/layers,0.,h/2]);q=Q
            models.append(pg.flat_model(n,h/layers,length,.2));origins.append(origin);rot.append(q)
    for j in range(layers):
        root=[(F(0),y,z) for y in ys for z in zz];right=[(F(1),y,z) for y in ys for z in zz];blocks.extend([[ (j,root,1)],[(j,right,1),(layers+j,root,-1)]])
    # Adjacent thickness layers share complete physical faces, with orientation
    # Q reversing the rib's local thickness direction.
    for part in [0,1,2]:
        for j in range(layers-1):
            a=part*layers+j;b=a+1;sa=F(1) if part<2 else F(-1);sb=-sa;blocks.append([(a,[(x,y,sa) for x in xs for y in ys],1),(b,[(x,y,sb) for x in xs for y in ys],-1)])
    for j in range(layers):
        # Pad surface x in [j/layers,(j+1)/layers], rib j bottom all thickness.
        # Split the restricted pad surface at every source polynomial cell.
        cuts=sorted({F(j,layers),F(j+1,layers)}|{F(k,2*n) for k in range(2*n+1) if F(j,layers)<F(k,2*n)<F(j+1,layers)})
        xj=[a+(b-a)*F(2*k+1,20) for a,b in zip(cuts[:-1],cuts[1:]) for k in range(10)];blocks.append([(2*layers-1,[(x,y,F(1)) for x in xj for y in ys],1),(2*layers+j,[(F(0),y,1-2*(layers*x-j)) for x in xj for y in ys],-1)])
    cache={}
    for m in models:
        if id(m) not in cache:cache[id(m)]=plane.plane(m)
    planes=[cache[id(m)] for m in models];Z0,record=plane.constraints(models,planes,blocks);Z=np.kron(Z0,np.eye(2));data=[plane.data(m,p,o,q) for m,p,o,q in zip(models,planes,origins,rot)];K=block_diag(*(d[0] for d in data));sizes=[len(d[0]) for d in data];offset=np.cumsum([0]+sizes);f=np.zeros(len(K))
    for j in range(layers):
        i=2*layers+j;p,w,area=pg.face(models[i],0,1.);N,_,_=pg.point(models[i],p);load=(N@planes[i]['P']).T@(w*np.linalg.norm(area,axis=1)/(h*.2));f[offset[i]:offset[i+1]:2]=load
    v=solve(Z.T@K@Z,Z.T@f,assume_a='pos');u=Z@v;res=float(np.linalg.norm(Z.T@(K@u-f))/np.linalg.norm(f));tip=float(f@u)
    rng=np.random.default_rng(149);gap=0.
    for maps in blocks:
        indices=rng.integers(len(maps[0][1]),size=(23,2));t=rng.uniform(size=(23,1));values=[]
        for i,points,sign in maps:
            points=np.array(points,float);pp=t*points[indices[:,0]]+(1-t)*points[indices[:,1]];N,_,_=pg.point(models[i],pp);v2=u[offset[i]:offset[i+1]].reshape(-1,2);q2=planes[i]['P']@v2;values.append(sign*N@q2)
        gap=max(gap,float(np.max(abs(sum(values)))))
    row=dict(h=h,n=n,layers=layers,tip=tip,free_residual=res,holdout_gap=gap,trace=record,plane_dimensions=[p['dimension'] for p in planes],exact_full_space_average_closure=all(p['exact_full_space_average_closure'] for p in planes),seconds=time.monotonic()-start);print('THICKNESS_FOLD',json.dumps(row),flush=True);return row

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--h',type=float,nargs='*',default=[.1,.05]);parser.add_argument('--n',type=int,nargs='*',default=[1,2,4]);parser.add_argument('--layers',type=int,nargs='*',default=[1,2,4]);parser.add_argument('--publish',action='store_true');args=parser.parse_args();rows=[];failures=[];references=[];start=time.monotonic()
    try:
        for h in args.h:
            refs=[reference(h,d) for d in [2,4,8]];references.append(dict(h=h,records=refs));print('REFERENCE',refs,flush=True)
            for layers in args.layers:
                for n in args.n:
                    r=candidate(h,n,layers);r['tip_error']=abs(r['tip']/refs[-1]['tip']-1);rows.append(r)
        if args.publish:
            assert args.h==[.1,.05] and args.n==[1,2,4] and args.layers==[1,2,4] and len(rows)==18
            for h in args.h:
                r=[x for x in rows if x['h']==h];fine=[x for x in r if x['n']==4 and x['layers']==4][0];coarse=[x for x in r if x['n']==1 and x['layers']==4][0];unresolved=[x for x in r if x['n']==4 and x['layers']==1][0];refs=next(x['records'] for x in references if x['h']==h)
                assert fine['tip_error']<.01 and fine['tip_error']<coarse['tip_error'] and fine['tip_error']<unresolved['tip_error'];assert abs(refs[-1]['tip']/refs[-2]['tip']-1)<.002
            assert all(r['free_residual']<1e-8 and r['holdout_gap']<1e-9 and r['exact_full_space_average_closure'] for r in rows)
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    if args.publish:
        paths=[Path(__file__),Path(plane.__file__),Path(pg.__file__),ROOT/'tools/design-experiments/connection_exact_trace.py',ROOT/'tools/design-experiments/check_connection_fold.py',ROOT/'docs/verification/connection-adoption-protocol.md'];dest=ROOT/'docs/verification/connection-fold-thickness-result.json';assert not dest.exists();dest.write_text(json.dumps(dict(passed=not failures and len(rows)==18,records=rows,references=references,failures=failures,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8');print('CONNECTION_FOLD_THICKNESS_RESULT',not failures and len(rows)==18,flush=True)
    else:print(json.dumps(rows,indent=2),flush=True)
if __name__=='__main__':main()
