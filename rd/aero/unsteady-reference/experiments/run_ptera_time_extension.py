import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,numpy as np,gc
from probe_ptera_reference import run
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/unsteady-2-3'
def main():
    for samples in [512,1024]:
        label=f'time{samples}';p=OUT/(label+'.json')
        if p.exists():continue
        print('START',label,flush=True);s,prob,r=run(2,8,samples,3,1e-4,core=1e-5);np.savez_compressed(p.with_suffix('.npz'),**r);p.write_text(json.dumps(dict(label=label,m=2,n=8,samples=samples,cycles=3,amplitude=1e-4,core=1e-5,nu=1e-12,seconds=r['seconds']),indent=2),encoding='utf-8');print('DONE',label,r['seconds'],flush=True);del s,prob,r;gc.collect()
if __name__=='__main__':main()
