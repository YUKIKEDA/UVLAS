"""Full-system fallback when floating-point symmetry condensation fails raw residual."""
from pathlib import Path
import json,time,gc,numpy as np
from scipy.linalg import solve
from uvlm_steady_reference import geometry,run as independent
from uvlm_frequency_fast import bound,segment
from check_q9_plate_convergence import memory_bytes
from trefftz_observations import midpoint
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/drag-common-limit-v2';OUT.mkdir(exist_ok=True);OLD=ROOT/'build/analysis/finite-wing-2-3'
def calculate(m,n,d,force_solve=False):
    start=time.monotonic();g=geometry(m,n,d);X=g['X'];N=g['normal'];idx=np.arange(m*n).reshape(m,n);right=idx[:,n//2:].ravel();left=idx[:,n//2-1::-1].ravel();eps=1.114285714285714e-5;source=OLD/f'candidate-d{d}-m{m}-n{n}.npz';rejected=None
    if source.exists() and not force_solve:gamma=np.load(source)['gamma'];mode='reused_full';res=None
    else:
        nh=n//2;X[:,:nh+1]=X[:,n:nh-1:-1]*[1,-1,1];pp=g['P'].reshape(m,n,3);pp[:,:nh]=pp[:,n-1:nh-1:-1]*[1,-1,1];nn=N.reshape(m,n,3);nn[:,:nh]=nn[:,n-1:nh-1:-1]*[1,-1,1];g['left']=X[:-1,:-1].reshape(-1,3);g['right']=X[:-1,1:].reshape(-1,3)
        tail=X.copy();tail[-1,:,0]+=1024.;B=bound(g['P'][right],N[right],tail,eps);A=B[:,right]+B[:,left];del B;rhs=-20*N[right,2];reduced=solve(A,rhs);res=float(np.linalg.norm(A@reduced-rhs)/np.linalg.norm(rhs));G=np.empty((m,n));G[:,n//2:]=reduced.reshape(m,n//2);G[:,:n//2]=G[:,n//2:][:,::-1];gamma=G.ravel();mode='symmetric';del A
    tail=X.copy();tail[-1,:,0]+=1024.;chosen=np.unique(np.r_[right[np.linspace(0,len(right)-1,8,dtype=int)],left[np.linspace(0,len(left)-1,8,dtype=int)]])
    def raw_residual():return float(np.linalg.norm(bound(g['P'][chosen],N[chosen],tail,eps)@gamma+20*N[chosen,2])/np.linalg.norm(20*N[chosen,2]))
    fullres=raw_residual()
    if fullres>=1e-11 and mode=='symmetric':
        rejected=fullres;print('REJECT_CONDENSATION',m,n,d,rejected,flush=True);A=bound(g['P'],N,tail,eps);rhs=-20*N[:,2];gamma=solve(A,rhs);res=float(np.linalg.norm(A@gamma-rhs)/np.linalg.norm(rhs));mode='full_after_symmetry_rejection';del A;fullres=raw_residual()
    assert fullres<1e-11 and (res is None or res<1e-11),(mode,res,fullres)
    G=gamma.reshape(m,n);horse=(G-np.vstack([np.zeros((1,n)),G[:-1]])).ravel();p=(g['left']+g['right'])/2;dl=g['right']-g['left'];normal=np.column_stack([np.zeros(len(dl)),-dl[:,2],dl[:,1]]);C=bound(p,normal,tail,eps)
    for row in range(m*n):C[row,row]-=segment(p[row],g['left'][row],g['right'][row],eps,normal[row])[0]
    near=-1.2*horse@(C@gamma);far=midpoint(X[-1],gamma[-n:]);del C;elapsed=time.monotonic()-start;peak=memory_bytes();assert elapsed<1200 and peak<24*2**30
    return dict(near_drag=float(near),trefftz_midpoint_drag=float(far),near_far_difference=float(near/far-1),residual=res,uncondensed_sample_residual=fullres,algebraic_mode=mode,rejected_condensation_residual=rejected,seconds=elapsed,peak_bytes=peak),gamma
def main():
    for d in [0,5]:
        _,gamma=calculate(4,12,d,True);full=independent(4,12,d,1.114285714285714e-5,1024.);error=float(np.linalg.norm(gamma-full['gamma'])/np.linalg.norm(full['gamma']));assert error<1e-10;print('COARSE_PROOF',d,error,flush=True)
    for d in [0,5]:
        for m,n in [(16,96),(32,96),(64,96),(16,192),(32,192),(64,192)]:
            p=OUT/f'candidate-d{d}-m{m}-n{n}.json'
            if p.exists():continue
            print('START',p.stem,flush=True);r,gamma=calculate(m,n,d);r.update(m=m,n=n,dihedral=d);np.savez_compressed(p.with_suffix('.npz'),gamma=gamma);p.write_text(json.dumps(r,indent=2),encoding='utf-8');print('DONE',p.stem,json.dumps(r),flush=True);gc.collect()
if __name__=='__main__':main()
