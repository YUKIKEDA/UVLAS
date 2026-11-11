"""The single predeclared M128/N192 extension, with hard per-case timeout."""
from pathlib import Path
import json,time,sys,subprocess,numpy as np
from scipy.linalg import solve
from uvlm_steady_reference import geometry
from uvlm_frequency_fast import bound,segment
from check_q9_plate_convergence import memory_bytes
from trefftz_observations import midpoint
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/drag-common-limit-v2'
def worker(d):
    start=time.monotonic();m,n=128,192;g=geometry(m,n,d);tail=g['X'].copy();tail[-1,:,0]+=1024.;eps=1.114285714285714e-5;A=bound(g['P'],g['normal'],tail,eps);rhs=-20*g['normal'][:,2];gamma=solve(A,rhs);res=float(np.linalg.norm(A@gamma-rhs)/np.linalg.norm(rhs));assert res<1e-11;del A
    chosen=np.linspace(0,m*n-1,32,dtype=int);R=bound(g['P'][chosen],g['normal'][chosen],tail,eps)@gamma-rhs[chosen];check=float(np.linalg.norm(R)/np.linalg.norm(rhs[chosen]));assert check<1e-11
    G=gamma.reshape(m,n);horse=(G-np.vstack([np.zeros((1,n)),G[:-1]])).ravel();p=(g['left']+g['right'])/2;dl=g['right']-g['left'];normal=np.column_stack([np.zeros(len(dl)),-dl[:,2],dl[:,1]]);C=bound(p,normal,tail,eps)
    for j in range(m*n):C[j,j]-=segment(p[j],g['left'][j],g['right'][j],eps,normal[j])[0]
    D=-1.2*horse@(C@gamma);far=midpoint(g['X'][-1],gamma[-n:]);peak=memory_bytes();seconds=time.monotonic()-start;assert peak<24*2**30;result=dict(m=m,n=n,dihedral=d,near_drag=float(D),trefftz_midpoint_drag=float(far),near_far_difference=float(D/far-1),residual=res,uncondensed_sample_residual=check,seconds=seconds,peak_bytes=peak,algebraic_mode='full_direct_after_previous_condensation_rejections');p=OUT/f'candidate-d{d}-m128-n192.json';np.savez_compressed(p.with_suffix('.npz'),gamma=gamma);p.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result),flush=True)
def main():
    if len(sys.argv)>1:worker(int(sys.argv[1]));return
    for d in [0,5]:
        p=OUT/f'candidate-d{d}-m128-n192.json'
        if p.exists():continue
        print('START',p.stem,flush=True)
        try:subprocess.run([sys.executable,str(Path(__file__).resolve()),str(d)],check=True,timeout=1200,cwd=ROOT)
        except subprocess.TimeoutExpired:
            (OUT/f'candidate-d{d}-m128-n192-timeout.json').write_text(json.dumps(dict(passed=False,reason='predeclared_1200_second_case_limit',dihedral=d)),encoding='utf-8');raise
if __name__=='__main__':main()
