"""Exploratory geometry/circulation equivalence; not a continuum acceptance."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np,json
from scipy.spatial import cKDTree
from uvlm_steady_reference import run,linear_load
ROOT=rd_paths.ROOT
def main():
    r=run(4,12,eps=1e-8,length=4096.);a=np.load(ROOT/'build/reference-aerosandbox/probe.npz');g=r['g'];dist,index=cKDTree(a['collocation']).query(g['P']);assert max(dist)<1e-12
    scale=np.sin(np.deg2rad(1.));hr=a['gamma'][index]/scale
    F,M,bands=linear_load(g,hr)
    print(json.dumps(dict(geometry_error=float(max(dist)),left_error=float(np.max(abs(g['left']-a['left'][index]))),right_error=float(np.max(abs(g['right']-a['right'][index]))),circulation_error=float(np.linalg.norm(r['horse']-hr)/np.linalg.norm(hr)),force_slope=r['F'].tolist(),ref_force_slope=F.tolist(),moment_slope=r['M'].tolist(),ref_moment_slope=M.tolist(),span_difference=float(np.max(abs(r['bands']-bands)/abs(bands))),residual=float(r['residual'])),indent=2))
if __name__=='__main__':main()
