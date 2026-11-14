import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np,json,time,gc
from scipy.linalg import solve
from uvlm_steady_reference import geometry
from uvlm_frequency_fast import bound
from uvlm_drag_force_kernel import matrix
from trefftz_observations import midpoint
ROOT=rd_paths.ROOT;OLD=ROOT/'build/analysis/finite-wing-2-3';OUT=ROOT/'build/analysis/drag-core-resolution';OUT.mkdir(exist_ok=True)
def main():
    for m,n in [(16,96),(32,96),(64,96)]:
        p=OUT/f'm{m}-n{n}.json'
        if p.exists():continue
        print('START',p.stem,flush=True);start=time.monotonic();g=geometry(m,n);X=g['X'].copy();X[-1,:,0]+=1024.;A=bound(g['P'],g['normal'],X,0.);rhs=-20*g['normal'][:,2];ideal=solve(A,rhs);res=float(np.linalg.norm(A@ideal-rhs)/np.linalg.norm(rhs));assert res<1e-11;del A
        source=OLD/f'candidate-d0-m{m}-n{n}.npz'
        if source.exists():regular=np.load(source)['gamma']
        else:regular=np.load(ROOT/f'build/analysis/drag-common-limit-v2/candidate-d0-m{m}-n{n}.npz')['gamma']
        P=(g['left']+g['right'])/2;dl=g['right']-g['left'];N=np.column_stack([np.zeros(m*n),-dl[:,2],dl[:,1]]);results=[]
        for eps in [0.,1.114285714285714e-5]:
            C=matrix(P,N,X,eps,np.arange(m*n))
            for label,Gamma in [('ideal',ideal),('regularized',regular)]:
                G=Gamma.reshape(m,n);h=(G-np.vstack([np.zeros((1,n)),G[:-1]])).ravel();D=-1.2*h@(C@Gamma);far=midpoint(g['X'][-1],Gamma[-n:]);results.append(dict(circulation=label,field_eps=eps,near_drag=float(D),trefftz_midpoint_drag=float(far),relative_difference=float(D/far-1)))
            del C
        r=dict(m=m,n=n,records=results,ideal_residual=res,seconds=time.monotonic()-start);np.savez_compressed(p.with_suffix('.npz'),gamma_ideal=ideal,gamma_regularized=regular);p.write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r),flush=True);gc.collect()
if __name__=='__main__':main()
