"""Diagnostic surface/panel normal with a directional differential."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
def normal(t1,t2,dt1=None,dt2=None):
    area=np.cross(t1,t2);length=np.linalg.norm(area,axis=-1)
    scale=np.linalg.norm(t1,axis=-1)*np.linalg.norm(t2,axis=-1)
    if np.any(~np.isfinite(length)) or np.any(scale==0) or np.any(length<=1e-12*scale):raise ValueError('degenerate_surface_or_panel')
    n=area/length[...,None]
    if dt1 is None:return n
    da=np.cross(dt1,t2)+np.cross(t1,dt2)
    dn=(da-n*np.sum(n*da,axis=-1)[...,None])/length[...,None]
    return n,dn
def panel(vertices,dvertices=None):
    # Two diagonals, oriented by the stated cyclic four-vertex connectivity.
    t1=vertices[...,2,:]-vertices[...,0,:];t2=vertices[...,3,:]-vertices[...,1,:]
    if dvertices is None:return normal(t1,t2)
    return normal(t1,t2,dvertices[...,2,:]-dvertices[...,0,:],dvertices[...,3,:]-dvertices[...,1,:])
