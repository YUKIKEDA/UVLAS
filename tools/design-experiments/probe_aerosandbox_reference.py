"""Read-only independent code probe; not an accepted physical reference yet."""
from pathlib import Path
import json,time,hashlib,sys
import numpy as np
import aerosandbox as asb
ROOT=Path(__file__).resolve().parents[2]
class FixedDensity(asb.Atmosphere):
    def density(self):return 1.2
def run(m=4,n=12,dihedral=0.,alpha=1.,core=1e-8):
    sections=[]
    for y in [-4.,0.,4.]:
        eta=abs(y)/4;sections.append(asb.WingXSec(xyz_le=[abs(y)*np.tan(np.deg2rad(35)),y,abs(y)*np.tan(np.deg2rad(dihedral))],chord=1.5-.9*eta,twist=0.,airfoil=asb.Airfoil('naca0012')))
    airplane=asb.Airplane(wings=[asb.Wing(name='wing',xsecs=sections,symmetric=False)],xyz_ref=[.375,0.,0.],s_ref=8.4,b_ref=8.,c_ref=1.114285714285714);op=asb.OperatingPoint(atmosphere=FixedDensity(),velocity=20.,alpha=alpha,beta=0.,p=0.,q=0.,r=0.);analysis=asb.VortexLatticeMethod(airplane=airplane,op_point=op,spanwise_resolution=n//2,chordwise_resolution=m,vortex_core_radius=core,align_trailing_vortices_with_wind=False,run_symmetric_if_possible=False);start=time.monotonic();aero=analysis.run();elapsed=time.monotonic()-start;U=op.compute_freestream_velocity_geometry_axes();return analysis,aero,U,elapsed
def main():
    a,result,U,seconds=run();out=ROOT/'build/reference-aerosandbox/probe.npz';np.savez_compressed(out,front_left=a.front_left_vertices,front_right=a.front_right_vertices,back_left=a.back_left_vertices,back_right=a.back_right_vertices,left=a.left_vortex_vertices,right=a.right_vortex_vertices,collocation=a.collocation_points,normals=a.normal_directions,gamma=a.vortex_strengths,F=a.forces_geometry,M=a.moments_geometry,U=U)
    print(json.dumps(dict(version=asb.__version__,numpy=np.__version__,panels=len(a.vortex_strengths),U=np.asarray(U).tolist(),force=np.asarray(result['F_g']).tolist(),moment=np.asarray(result['M_g']).tolist(),CL=float(result['CL']),seconds=seconds,snapshot=out.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(out.read_bytes()).hexdigest()),indent=2))
if __name__=='__main__':main()
