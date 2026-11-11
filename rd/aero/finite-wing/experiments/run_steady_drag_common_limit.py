import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,time,gc,numpy as np
from scipy.linalg import solve
from uvlm_steady_reference import geometry
from uvlm_frequency_fast import bound,segment
from check_q9_plate_convergence import memory_bytes
from trefftz_observations import midpoint
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/drag-common-limit';OUT.mkdir(exist_ok=True)
OLD=ROOT/'build/analysis/finite-wing-2-3'
def calculate(m,n,d,force_solve=False):
    start=time.monotonic();g=geometry(m,n,d);X=g['X'];N=g['normal'];idx=np.arange(m*n).reshape(m,n);right=idx[:,n//2:].ravel();left=idx[:,n//2-1::-1].ravel();eps=1.114285714285714e-5;tail=X.copy();tail[-1,:,0]+=1024.
    source=OLD/f'candidate-d{d}-m{m}-n{n}.npz'
    if source.exists() and not force_solve:gamma=np.load(source)['gamma'];res=None;reused=True
    else:
        # Enforce the exact known reflection before symmetry condensation.
        # Independently sampled cosine coordinates differ by up to ~1e-15 m.
        nh=n//2;X[:,:nh+1]=X[:,n:nh-1:-1]*[1,-1,1];pp=g['P'].reshape(m,n,3);pp[:,:nh]=pp[:,n-1:nh-1:-1]*[1,-1,1];nn=N.reshape(m,n,3);nn[:,:nh]=nn[:,n-1:nh-1:-1]*[1,-1,1];g['left']=X[:-1,:-1].reshape(-1,3);g['right']=X[:-1,1:].reshape(-1,3);tail=X.copy();tail[-1,:,0]+=1024.
        B=bound(g['P'][right],N[right],tail,eps);A=B[:,right]+B[:,left];del B;rhs=-20*N[right,2];reduced=solve(A,rhs);res=float(np.linalg.norm(A@reduced-rhs)/np.linalg.norm(rhs));G=np.empty((m,n));G[:,n//2:]=reduced.reshape(m,n//2);G[:,:n//2]=G[:,n//2:][:,::-1];gamma=G.ravel();assert res<1e-11;reused=False;del A
    G=gamma.reshape(m,n);horse=(G-np.vstack([np.zeros((1,n)),G[:-1]])).ravel();p=(g['left'][right]+g['right'][right])/2;dl=g['right'][right]-g['left'][right];normal=np.column_stack([np.zeros(len(dl)),-dl[:,2],dl[:,1]]);C=bound(p,normal,tail,eps)
    for row,col in enumerate(right):C[row,col]-=segment(p[row],g['left'][col],g['right'][col],eps,normal[row])[0]
    near=-2*1.2*horse[right]@(C@gamma);far=midpoint(X[-1],gamma[-n:]);del C
    chosen=np.unique(np.r_[right[np.linspace(0,len(right)-1,8,dtype=int)],left[np.linspace(0,len(left)-1,8,dtype=int)]]);R=bound(g['P'][chosen],N[chosen],tail,eps)@gamma+20*N[chosen,2];fullres=float(np.linalg.norm(R)/np.linalg.norm(20*N[chosen,2]));assert fullres<1e-11,('uncondensed_residual',fullres)
    peak=memory_bytes();elapsed=time.monotonic()-start;assert elapsed<1200 and peak<24*2**30
    return dict(near_drag=float(near),trefftz_midpoint_drag=float(far),near_far_difference=float(near/far-1),residual=res,uncondensed_sample_residual=fullres,reused_circulation=reused,seconds=elapsed,peak_bytes=peak),gamma
def main():
    from uvlm_steady_reference import run as independent
    for d in [0,5]:
        _,gamma=calculate(4,12,d,True);full=independent(4,12,d,1.114285714285714e-5,1024.);error=float(np.linalg.norm(gamma-full['gamma'])/np.linalg.norm(full['gamma']));assert error<1e-10;print('SYMMETRY_PROOF',d,error,flush=True)
    for d in [0,5]:
        for m,n in [(16,96),(32,96),(64,96),(16,192),(32,192),(64,192)]:
            p=OUT/f'candidate-d{d}-m{m}-n{n}.json'
            if p.exists():continue
            print('START',p.stem,flush=True);r,gamma=calculate(m,n,d);r.update(m=m,n=n,dihedral=d);np.savez_compressed(p.with_suffix('.npz'),gamma=gamma);p.write_text(json.dumps(r,indent=2),encoding='utf-8');print('DONE',p.stem,json.dumps(r),flush=True);gc.collect()
if __name__=='__main__':main()
