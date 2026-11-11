"""Ideal-core reference using the separately qualified native scalar projection."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import sys,json,time,numpy as np
import uvlm_frequency_delta as common
import ptera_scalar_reference as native
from ptera_frequency_reference import geometry
from run_frequency_delta import peak_bytes
common.baseline=native.baseline;common.delta=native.delta
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/frequency-ideal-reference'
def main():
    m,n=map(int,sys.argv[1:3]);OUT.mkdir(exist_ok=True);start=time.monotonic();g=geometry(m,n);solver=common.Solver(g,'reference',eps=0.)
    for k in [.05,.2,.8]:
        p=OUT/f'reference-m{m}-n{n}-k{k}.json';assert not p.exists();r=solver.run(k);peak=peak_bytes();assert peak<24*2**30;np.savez_compressed(p.with_suffix('.npz'),gamma=r['gamma'],transfer=r['transfer'],W=r['W']);meta=dict(side='reference',core=0.,m=m,n=n,k=k,residuals=r['residuals'],refinement_steps=r['refinement_steps'],seconds=time.monotonic()-start,peak_bytes=peak,transfer=[[[float(z.real),float(z.imag)] for z in v] for v in r['transfer']]);p.write_text(json.dumps(meta,indent=2),encoding='utf-8');print(json.dumps(meta),flush=True)
if __name__=='__main__':main()
