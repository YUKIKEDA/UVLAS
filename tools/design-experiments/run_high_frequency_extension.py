from pathlib import Path
import argparse,json,time,numpy as np
from frequency_motion_columns import schur_solve,observe
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/pitch-frequency-2-3'
def main():
    p=argparse.ArgumentParser();p.add_argument('side');side=p.parse_args().side;m,n,k=64,96,.8;dest=OUT/f'{side}-m{m}-n{n}-k{k}-q16.json'
    if dest.exists():return
    print('START',dest.stem,flush=True);start=time.monotonic();omega=40*k/1.114285714285714
    if side=='candidate':
        from uvlm_frequency_fast import geometry,bound,age_nodes,age_integral
        g=geometry(m,n);A0=bound(g['P'],g['normal'],g['X'],1e-5);ages,weights=age_nodes(g,omega,128.,16);W=age_integral(g['P'],g['normal'],g['X'][-1],1e-5,ages,weights)
    else:
        from ptera_frequency_reference import run,project
        r=run(m,n,k,1e-5,128.,16,True);g=r['g'];W=r['W'];X=g['X'];A0=project(g,X[:-1,:-1].reshape(-1,3),X[:-1,1:].reshape(-1,3),X[1:,1:].reshape(-1,3),X[1:,:-1].reshape(-1,3),1e-5)
    nz=g['normal'][:,2];rhs=np.column_stack([1j*omega*nz,-(20+1j*omega*(g['P'][:,0]-.375))*nz]);gamma,res,mode=schur_solve(A0,W,rhs,n);transfer=observe(g,gamma,omega,side);np.savez_compressed(dest.with_suffix('.npz'),gamma=gamma,transfer=transfer,W=W);meta=dict(side=side,m=m,n=n,k=k,order=16,transfer=[[[float(z.real),float(z.imag)] for z in motion] for motion in transfer],residual=res,algebraic_mode=mode,seconds=time.monotonic()-start);dest.write_text(json.dumps(meta,indent=2),encoding='utf-8');print(json.dumps(meta),flush=True)
if __name__=='__main__':main()
