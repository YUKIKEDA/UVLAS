from pathlib import Path
import numpy as np,json,gc
from probe_aerosandbox_reference import run
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/drag-common-limit'
def main():
    for d in [0,5]:
        for m,n in [(8,96),(8,192),(8,384),(16,384)]:
            p=OUT/f'reference-d{d}-m{m}-n{n}-core1e-10.npz'
            if p.exists():continue
            print('START',p.stem,flush=True);a,_,_,elapsed=run(m,n,d,core=1e-10);np.savez_compressed(p,horse=a.vortex_strengths/np.sin(np.deg2rad(1.)),left=a.left_vortex_vertices,right=a.right_vortex_vertices,collocation=a.collocation_points,core=1e-10);print('DONE',p.stem,elapsed,flush=True);del a;gc.collect()
if __name__=='__main__':main()
