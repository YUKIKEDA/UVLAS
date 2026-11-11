from pathlib import Path
import numpy as np,json,gc,time
from probe_aerosandbox_reference import run
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/drag-common-limit';OUT.mkdir(exist_ok=True)
def main():
    for d in [0,5]:
        for core in [1e-9,1e-10]:
            p=OUT/f'reference-d{d}-m8-n768-core{core}.npz'
            if p.exists():continue
            print('START',p.stem,flush=True);a,_,_,elapsed=run(8,768,d,core=core);np.savez_compressed(p,horse=a.vortex_strengths/np.sin(np.deg2rad(1.)),left=a.left_vortex_vertices,right=a.right_vortex_vertices,collocation=a.collocation_points,core=core);print('DONE',p.stem,elapsed,flush=True);del a;gc.collect()
if __name__=='__main__':main()
