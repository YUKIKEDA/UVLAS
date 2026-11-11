from pathlib import Path
import argparse,json,time,gc,numpy as np
from scipy.linalg import solve
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/finite-wing-2-3'
def main():
    p=argparse.ArgumentParser();p.add_argument('side');args=p.parse_args()
    if args.side=='reference':from probe_aerosandbox_reference import run
    else:
        from uvlm_frequency_fast import bound
        from uvlm_steady_reference import geometry,linear_load
    for d in [0,5]:
        for n in [192,384,768]:
            m=8;label=f'{args.side}-d{d}-band-n{n}';p=OUT/(label+'.json')
            if p.exists():continue
            start=time.monotonic();print('START',label,flush=True)
            if args.side=='candidate':
                g=geometry(m,n,d);X=g['X'].copy();X[-1,:,0]+=1024.;A=bound(g['P'],g['normal'],X,0.);rhs=-20*g['normal'][:,2];gamma=solve(A,rhs);res=float(np.linalg.norm(A@gamma-rhs)/np.linalg.norm(rhs));G=gamma.reshape(m,n);horse=(G-np.vstack([np.zeros((1,n)),G[:-1]])).ravel();F,M,bands=linear_load(g,horse);r=dict(F=F.tolist(),M=M.tolist(),bands=bands.tolist(),residual=res,eps=0.,length=1024.);np.savez_compressed(p.with_suffix('.npz'),gamma=gamma,horse=horse,F=F,M=M,bands=bands)
            else:
                a,r0,U,elapsed=run(m,n,d);horse=a.vortex_strengths/np.sin(np.deg2rad(1.));F=1.2*horse[:,None]*np.cross([20.,0,0],a.right_vortex_vertices-a.left_vortex_vertices);M=np.cross((a.left_vortex_vertices+a.right_vortex_vertices)/2-[.375,0,0],F);bands=[]
                for lo in np.linspace(-4,3.6,20):
                    f=np.maximum(0,np.minimum(a.right_vortex_vertices[:,1],lo+.4)-np.maximum(a.left_vortex_vertices[:,1],lo))/(a.right_vortex_vertices[:,1]-a.left_vortex_vertices[:,1]);bands.append(float(f@F[:,2]))
                r=dict(F=F.sum(axis=0).tolist(),M=M.sum(axis=0).tolist(),bands=bands,core=1e-8);np.savez_compressed(p.with_suffix('.npz'),horse=horse,F=F.sum(axis=0),M=M.sum(axis=0),bands=bands);del a
            r.update(m=m,n=n,dihedral=d,seconds=time.monotonic()-start);p.write_text(json.dumps(r,indent=2),encoding='utf-8');print('DONE',label,json.dumps(r),flush=True);gc.collect()
if __name__=='__main__':main()
