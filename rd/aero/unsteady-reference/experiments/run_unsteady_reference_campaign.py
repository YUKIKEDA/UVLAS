"""Finite Ptera campaign and retained original/audit-core parity."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,time,gc
import numpy as np
from probe_ptera_reference import run
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/unsteady-2-3';OUT.mkdir(exist_ok=True)
def main():
    tasks=[('original',2,8,64,3,1e-4,None,1e-12),('original-core-parity',2,8,64,3,1e-4,.0315,1e-12)]
    tasks += [(f'time{s}',2,8,s,3,1e-4,1e-5,1e-12) for s in [64,128,256]]
    tasks += [('cycles6',2,8,256,6,1e-4,1e-5,1e-12),('core1e-3',2,8,256,3,1e-4,1e-3,1e-12),('core1e-4',2,8,256,3,1e-4,1e-4,1e-12),('amplitude1e-5',2,8,256,3,1e-5,1e-5,1e-12),('nu1e-13',2,8,256,3,1e-4,1e-5,1e-13),('space4x24',4,24,256,3,1e-4,1e-5,1e-12),('space8x48',8,48,256,3,1e-4,1e-5,1e-12)]
    for label,m,n,s,cycles,amp,core,nu in tasks:
        dest=OUT/(label+'.json')
        if dest.exists():continue
        print('START',label,flush=True);solver,problem,r=run(m,n,s,cycles,amp,core=core,nu=nu);last=r['t']>r['t'][-1]-2*np.pi/r['omega']-1e-9;D=np.column_stack([np.ones(last.sum()),np.cos(r['omega']*r['t'][last]),np.sin(r['omega']*r['t'][last])]);L=-r['F'][:,2];My=r['M'][:,1]-.375*L;fit=np.linalg.lstsq(D,np.column_stack([L[last],My[last]]),rcond=None)[0];coef=(fit[1]-1j*fit[2])/amp
        np.savez_compressed(dest.with_suffix('.npz'),**r,transfer=coef)
        result=dict(label=label,m=m,n=n,samples=s,cycles=cycles,amplitude=amp,core=core,nu=nu,seconds=r['seconds'],transfer=[[float(x.real),float(x.imag)] for x in coef]);dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('DONE',label,json.dumps(result),flush=True);del solver,problem,r;gc.collect()
if __name__=='__main__':main()
