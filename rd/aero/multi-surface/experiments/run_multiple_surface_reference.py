import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,time,numpy as np
from probe_aerosandbox_reference import FixedDensity
import aerosandbox as asb
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/multiple-surface-2-3';OUT.mkdir(exist_ok=True)
CONFIG=[(4.,1.5,.6,35.,[0.,0.,0.]),(2.5,1.,.4,20.,[5.,0.,.6])]
def main():
    for m,n in [(2,12),(4,24),(8,48)]:
        wings=[]
        for s,cr,ct,sweep,offset in CONFIG:
            sections=[asb.WingXSec(xyz_le=(np.array(offset)+[abs(y)*np.tan(np.deg2rad(sweep)),y,0.]).tolist(),chord=cr+(ct-cr)*abs(y)/s,airfoil=asb.Airfoil('naca0012')) for y in [-s,0.,s]];wings.append(asb.Wing(xsecs=sections,symmetric=False))
        for label,selected in [('coupled',wings),('main',[wings[0]]),('tail',[wings[1]])]:
            p=OUT/f'reference-{label}-m{m}-n{n}.npz'
            if p.exists():continue
            start=time.monotonic();plane=asb.Airplane(wings=selected,xyz_ref=[0,0,0],s_ref=8.4,b_ref=8,c_ref=1.114285714285714);op=asb.OperatingPoint(atmosphere=FixedDensity(),velocity=20,alpha=1,beta=0);a=asb.VortexLatticeMethod(plane,op,spanwise_resolution=n//2,chordwise_resolution=m,vortex_core_radius=1e-9,align_trailing_vortices_with_wind=False);a.run();np.savez_compressed(p,left=a.left_vortex_vertices,right=a.right_vortex_vertices,collocation=a.collocation_points,horse=a.vortex_strengths/np.sin(np.deg2rad(1.)));print(label,m,n,time.monotonic()-start,flush=True)
if __name__=='__main__':main()
