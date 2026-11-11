from pathlib import Path
import sys,time,json,numpy as np
from scipy.linalg import lu_factor,lu_solve
from uvlm_frequency_fast import geometry,age_nodes,bound,run as full_run
from uvlm_age_scalar import age_integral
from uvlm_rm_scalar_symmetric import condensed
from frequency_motion_columns import schur_solve,observe
from run_frequency_delta import peak_bytes
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/high-frequency-final'
def run(m,n,eps):
    start=time.monotonic();g=geometry(m,n);nh=n//2;X=g['X'];P=g['P'].reshape(m,n,3);oldX=X.copy();oldP=P.copy();X[:,:nh+1]=X[:,n:nh-1:-1]*[1,-1,1];P[:,:nh]=P[:,n-1:nh-1:-1]*[1,-1,1];g['P']=P.reshape(-1,3);g['left']=X[:-1,:-1].reshape(-1,3);g['right']=X[:-1,1:].reshape(-1,3);rounding=max(np.max(abs(X-oldX)),np.max(abs(P-oldP)));assert rounding<1e-14;idx=np.arange(m*n).reshape(m,n);right=idx[:,nh:].ravel();left=idx[:,nh-1::-1].ravel();pp=g['P'][right];nn=g['normal'][right];tail=X.copy();tail[-1,:,0]+=128.;A=condensed(pp,nn,tail,eps);omega=32/1.114285714285714;ages,w=age_nodes(g,omega,128.,16);phase=omega*ages;weights=abs(w)*(-2*np.sin(phase/2)**2-1j*np.sin(phase));W=age_integral(pp,nn,X[-1],eps,ages,weights);W=W[:,nh:]+W[:,nh-1::-1];rhs=np.column_stack([1j*omega*nn[:,2],-(20+1j*omega*(pp[:,0]-.375))*nn[:,2]]);gamma,res,mode=schur_solve(A,W,rhs,nh);GG=np.empty((m,n,2),complex);GG[:,nh:]=gamma.reshape(m,nh,2);GG[:,:nh]=GG[:,nh:][:,::-1];full=GG.reshape(m*n,2);selected=np.unique(np.r_[right[np.linspace(0,len(right)-1,16,dtype=int)],left[np.linspace(0,len(left)-1,16,dtype=int)]]);B=bound(g['P'][selected],g['normal'][selected],tail,eps);WW=age_integral(g['P'][selected],g['normal'][selected],X[-1],eps,ages,weights);rr=np.column_stack([1j*omega*g['normal'][selected,2],-(20+1j*omega*(g['P'][selected,0]-.375))*g['normal'][selected,2]]);sample=np.linalg.norm(B@full+WW@full[-n:]-rr,axis=0)/np.linalg.norm(rr,axis=0);assert max(sample)<1e-11;transfer=observe(g,full,omega,'candidate');return dict(gamma=full,transfer=transfer,residual=res,uncondensed_sample_residuals=sample.tolist(),coordinate_rounding=float(rounding),seconds=time.monotonic()-start,peak_bytes=peak_bytes(),algebraic_mode=mode)
def main():
    OUT.mkdir(exist_ok=True)
    if sys.argv[1]=='proof':
        records=[]
        for m,n in [(2,8),(4,24)]:
            for eps in [0.,1e-6,1e-5]:
                r=run(m,n,eps);columns=[full_run(m,n,k=.8,eps=eps,length=128.,order=16,motion=motion) for motion in ['heave','pitch']];G=np.column_stack([v['gamma'] for v in columns]);tr=np.array([[v['F'][2],v['M'][1]] for v in columns]);errors=dict(gamma=float(np.linalg.norm(G-r['gamma'])/np.linalg.norm(G)),transfer=float(np.max(abs((tr-r['transfer'])/tr))));assert max(errors.values())<1e-10;records.append(dict(m=m,n=n,core=eps,errors=errors))
        import hashlib
        paths=[Path(__file__),ROOT/'tools/design-experiments/uvlm_rm_scalar_symmetric.py',ROOT/'docs/verification/high-frequency-final-space-core-protocol.md'];dest=ROOT/'docs/verification/high-frequency-final-symmetry-proof.json';assert not dest.exists();dest.write_text(json.dumps(dict(passed=True,records=records,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8');print(json.dumps(records),flush=True);return
    m=int(sys.argv[1]);r=run(m,192,1e-6);assert r['peak_bytes']<24*2**30 and r['seconds']<1200;p=OUT/f'candidate-m{m}-n192-core1e-6.json';assert not p.exists();np.savez_compressed(p.with_suffix('.npz'),gamma=r.pop('gamma'),transfer=r.pop('transfer'));r.update(m=m,n=192,k=.8,core=1e-6,physical_unknowns=m*192,algebraic_unknowns=m*96);p.write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r),flush=True)
if __name__=='__main__':main()
