import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import argparse,json,time,gc,numpy as np
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/unsteady-2-3'
def main():
    p=argparse.ArgumentParser();p.add_argument('side');p.add_argument('--intersection',action='store_true');args=p.parse_args()
    if args.side=='candidate':from uvlm_frequency_fast import run
    else:from ptera_frequency_reference import run
    grids=[(32,96),(64,96),(32,192),(16,384)] if not args.intersection else [(64,192)]
    for m,n in grids:
        label=f'{args.side}-extension-m{m}-n{n}';p=OUT/(label+'.json')
        if p.exists():continue
        print('START',label,flush=True);start=time.monotonic();r=run(m,n,length=128.,order=12);coef=[r['F'][2],r['M'][1]];result=dict(label=label,m=m,n=n,length=128.,order=12,eps=1e-5,transfer=[[float(z.real),float(z.imag)] for z in coef],residual=r['residual'],seconds=time.monotonic()-start);np.savez_compressed(p.with_suffix('.npz'),gamma=r['gamma'],F=r['F'],M=r['M'],W=r['W']);p.write_text(json.dumps(result,indent=2),encoding='utf-8');print('DONE',label,json.dumps(result),flush=True);del r;gc.collect()
if __name__=='__main__':main()
