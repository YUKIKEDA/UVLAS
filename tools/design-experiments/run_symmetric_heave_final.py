from pathlib import Path
import json,time,hashlib,numpy as np
from scipy.linalg import solve
from uvlm_frequency_fast import geometry,bound,age_nodes,age_integral,run as full_run
from check_q9_plate_convergence import memory_bytes
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/unsteady-2-3'
def run(m,n):
    start=time.monotonic();g=geometry(m,n);nh=n//2;X=g['X'];P=g['P'].reshape(m,n,3);origX=X.copy();origP=P.copy();X[:,:nh+1]=X[:,n:nh-1:-1]*[1,-1,1];P[:,:nh]=P[:,n-1:nh-1:-1]*[1,-1,1];g['P']=P.reshape(-1,3);g['left']=X[:-1,:-1].reshape(-1,3);g['right']=X[:-1,1:].reshape(-1,3);rounding=max(np.max(abs(X-origX)),np.max(abs(P-origP)));assert rounding<1e-14
    indices=np.arange(m*n).reshape(m,n);right=indices[:,nh:].ravel();left=indices[:,nh-1::-1].ravel();pp=g['P'][right];nn=g['normal'][right];eps=1e-5;omega=8/1.114285714285714;ages,weights=age_nodes(g,omega,128.,12);W=age_integral(pp,nn,X[-1],eps,ages,weights);B=bound(pp,nn,X,eps);A=np.array(B[:,right]+B[:,left],dtype=complex,order='F');del B;A[:,-nh:]+=W[:,nh:]+W[:,nh-1::-1];rhs=1j*omega*nn[:,2];gamma=solve(A,rhs);res=float(np.linalg.norm(A@gamma-rhs)/np.linalg.norm(rhs));G=np.empty((m,n),complex);G[:,nh:]=gamma.reshape(m,nh);G[:,:nh]=G[:,nh:][:,::-1];GG=G.ravel()
    # Check uncondensed rows on both halves, without allocating the full system.
    chosen=np.unique(np.r_[right[np.linspace(0,len(right)-1,12,dtype=int)],left[np.linspace(0,len(left)-1,12,dtype=int)]]);C=bound(g['P'][chosen],g['normal'][chosen],X,eps);WC=age_integral(g['P'][chosen],g['normal'][chosen],X[-1],eps,ages,weights);check=C@GG+WC@G[-1]-1j*omega*g['normal'][chosen,2];fullres=float(np.linalg.norm(check)/np.linalg.norm(1j*omega*g['normal'][chosen,2]));assert max(res,fullres)<1e-11
    horse=(G-np.vstack([np.zeros((1,n)),G[:-1]])).ravel();edgeF=1.2*horse[:,None]*np.cross([20.,0,0],g['right']-g['left']);F=edgeF.sum(axis=0);M=np.cross((g['left']+g['right'])/2-[.375,0,0],edgeF).sum(axis=0);quads=np.stack([X[:-1,:-1],X[:-1,1:],X[1:,1:],X[1:,:-1]],axis=2).reshape(-1,4,3)
    for ids in [[0,2,1],[0,3,2]]:
        tri=quads[:,ids];area=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])/2;f=1.2j*omega*GG[:,None]*area;F+=f.sum(axis=0);M+=np.cross(tri.mean(axis=1)-[.375,0,0],f).sum(axis=0)
    return dict(gamma=GG,F=F,M=M,residual=res,uncondensed_sample_residual=fullres,coordinate_rounding=float(rounding),peak_bytes=memory_bytes(),seconds=time.monotonic()-start)
def main():
    proof=[]
    for m,n in [(2,8),(4,24)]:
        a=run(m,n);b=full_run(m,n,length=128.);errs={key:float(np.linalg.norm(a[key]-b[key])/np.linalg.norm(b[key])) for key in ['gamma','F','M']};assert max(errs.values())<1e-10;proof.append(dict(m=m,n=n,errors=errs));print('PROOF',proof[-1],flush=True)
    print('START symmetric-m128-n192',flush=True);r=run(128,192);assert r['seconds']<1200 and r['peak_bytes']<24*2**30;coef=np.array([r['F'][2],r['M'][1]]);np.savez_compressed(OUT/'candidate-symmetric-m128-n192.npz',gamma=r['gamma'],F=r['F'],M=r['M']);meta={k:v for k,v in r.items() if k not in ['gamma','F','M']};meta.update(m=128,n=192,physical_unknowns=24576,algebraic_unknowns=12288,transfer=[[float(z.real),float(z.imag)] for z in coef],symmetry_proofs=proof);(OUT/'candidate-symmetric-m128-n192.json').write_text(json.dumps(meta,indent=2),encoding='utf-8');print(json.dumps(meta,indent=2),flush=True)
if __name__=='__main__':main()
