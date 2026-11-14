"""Bounded 2.3 campaign. Run in either independent environment; cache raw evidence."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import argparse,json,time,hashlib,gc
import numpy as np
ROOT=rd_paths.ROOT
OUT=ROOT/'build/analysis/finite-wing-2-3';OUT.mkdir(exist_ok=True)
def main():
    p=argparse.ArgumentParser();p.add_argument('side',choices=['candidate','reference']);p.add_argument('--max-panels',type=int,default=6144);args=p.parse_args()
    if args.side=='reference':
        from probe_aerosandbox_reference import run
    else:from uvlm_steady_reference import run
    grids=[(4,12),(4,24),(8,24),(8,48),(16,48),(16,96),(32,96),(16,192),(32,192)]
    for d in [0.,5.]:
        for m,n in grids:
            if m*n>args.max_panels:continue
            name=f'{args.side}-d{int(d)}-m{m}-n{n}';dest=OUT/(name+'.json')
            if dest.exists():print('CACHE',name,flush=True);continue
            start=time.monotonic();print('START',name,flush=True)
            if args.side=='candidate':
                a=run(m,n,d,eps=1.114285714285714e-5,length=1024.)
                np.savez_compressed(OUT/(name+'.npz'),gamma=a['gamma'],horse=a['horse'],left=a['g']['left'],right=a['g']['right'],collocation=a['g']['P'],normal=a['g']['normal'],F=a['F'],M=a['M'],bands=a['bands'])
                result=dict(F=a['F'].tolist(),M=a['M'].tolist(),bands=a['bands'].tolist(),residual=float(a['residual']),eps=1.114285714285714e-5,length=1024.)
            else:
                a,r,U,elapsed=run(m,n,d);h=a.vortex_strengths/np.sin(np.deg2rad(1.));left=a.left_vortex_vertices;right=a.right_vortex_vertices;F=1.2*h[:,None]*np.cross([20.,0,0],right-left);M=np.cross((left+right)/2-[.375,0,0],F);bands=np.zeros(20);bounds=np.linspace(-4,4,21)
                for k,(lo,hi) in enumerate(zip(bounds[:-1],bounds[1:])):
                    fraction=np.maximum(0,np.minimum(right[:,1],hi)-np.maximum(left[:,1],lo))/(right[:,1]-left[:,1]);bands[k]=fraction@F[:,2]
                np.savez_compressed(OUT/(name+'.npz'),horse=h,left=left,right=right,collocation=a.collocation_points,normal=a.normal_directions,F=F.sum(axis=0),M=M.sum(axis=0),bands=bands,raw_F=a.forces_geometry,raw_M=a.moments_geometry,raw_U=U)
                result=dict(F=F.sum(axis=0).tolist(),M=M.sum(axis=0).tolist(),bands=bands.tolist(),raw_F=np.asarray(r['F_g']).tolist(),raw_M=np.asarray(r['M_g']).tolist(),core=1e-8)
            result.update(side=args.side,m=m,n=n,dihedral=d,seconds=time.monotonic()-start,snapshot_sha256=hashlib.sha256((OUT/(name+'.npz')).read_bytes()).hexdigest())
            dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('DONE',name,json.dumps(result),flush=True);del a;gc.collect()
if __name__=='__main__':main()
