from pathlib import Path
import argparse,json,time,gc,numpy as np
from frequency_motion_columns import schur_solve,observe
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/pitch-frequency-2-3';OUT.mkdir(exist_ok=True)
def main():
    p=argparse.ArgumentParser();p.add_argument('side');args=p.parse_args()
    if args.side=='candidate':from uvlm_frequency_fast import geometry,bound,age_nodes,age_integral
    else:from ptera_frequency_reference import geometry,project,run
    for m,n in [(2,8),(4,24),(8,48),(16,96),(32,96)]:
        g=geometry(m,n);X=g['X']
        if args.side=='candidate':A0=bound(g['P'],g['normal'],X,1e-5)
        else:A0=project(g,X[:-1,:-1].reshape(-1,3),X[:-1,1:].reshape(-1,3),X[1:,1:].reshape(-1,3),X[1:,:-1].reshape(-1,3),1e-5)
        for k in ([.2] if m==2 else [.05,.2,.8]):
            p=OUT/f'{args.side}-m{m}-n{n}-k{k}.json'
            if p.exists():continue
            print('START',p.stem,flush=True);start=time.monotonic();omega=40*k/1.114285714285714
            if args.side=='candidate':ages,weights=age_nodes(g,omega,128.,12);W=age_integral(g['P'],g['normal'],X[-1],1e-5,ages,weights)
            else:baseline=run(m,n,k,1e-5,128.,12,True);W=baseline['W']
            nz=g['normal'][:,2];rhs=np.column_stack([1j*omega*nz,-(20+1j*omega*(g['P'][:,0]-.375))*nz]);gamma,res,mode=schur_solve(A0,W,rhs,n);transfer=observe(g,gamma,omega,args.side);direct=None
            if m<=4:
                from scipy.linalg import solve
                A=A0.astype(complex);A[:,-n:]+=W;d=solve(A,rhs);direct=float(np.linalg.norm(d-gamma)/np.linalg.norm(d));assert direct<1e-10
            if args.side=='reference':assert np.linalg.norm(gamma[:,0]-baseline['gamma'])/np.linalg.norm(gamma[:,0])<1e-10
            meta=dict(side=args.side,m=m,n=n,k=k,motions=['heave_m','nose_up_pitch_rad'],transfer=[[[float(z.real),float(z.imag)] for z in motion] for motion in transfer],residual=res,algebraic_mode=mode,direct_gamma_error=direct,seconds=time.monotonic()-start);np.savez_compressed(p.with_suffix('.npz'),gamma=gamma,transfer=transfer,W=W);p.write_text(json.dumps(meta,indent=2),encoding='utf-8');print('DONE',p.stem,json.dumps(meta),flush=True);gc.collect()
if __name__=='__main__':main()
