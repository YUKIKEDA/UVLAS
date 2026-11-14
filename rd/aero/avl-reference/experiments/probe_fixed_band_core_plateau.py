import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,numpy as np,gc
from probe_aerosandbox_reference import run
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/finite-wing-2-3'
def main():
    for d in [0,5]:
        dest=OUT/f'reference-d{d}-band768-core1e-10.json'
        if dest.exists():continue
        print('START',dest.stem,flush=True);a,_,_,elapsed=run(8,768,d,core=1e-10);horse=a.vortex_strengths/np.sin(np.deg2rad(1.));F=1.2*horse[:,None]*np.cross([20.,0,0],a.right_vortex_vertices-a.left_vortex_vertices);M=np.cross((a.left_vortex_vertices+a.right_vortex_vertices)/2-[.375,0,0],F);bands=[]
        for lo in np.linspace(-4,3.6,20):
            f=np.maximum(0,np.minimum(a.right_vortex_vertices[:,1],lo+.4)-np.maximum(a.left_vortex_vertices[:,1],lo))/(a.right_vortex_vertices[:,1]-a.left_vortex_vertices[:,1]);bands.append(float(f@F[:,2]))
        r=dict(F=F.sum(axis=0).tolist(),M=M.sum(axis=0).tolist(),bands=bands,core=1e-10,seconds=elapsed);dest.write_text(json.dumps(r,indent=2),encoding='utf-8');print('DONE',dest.stem,flush=True);del a;gc.collect()
if __name__=='__main__':main()
