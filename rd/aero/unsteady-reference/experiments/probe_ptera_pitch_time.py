import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np,json,time,gc
import pterasoftware as ps
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/pitch-frequency-2-3';OUT.mkdir(exist_ok=True)
def run(samples=256,cycles=3):
    foil=ps.geometry.airfoil.Airfoil('NACA0012');sections=[ps.geometry.wing_cross_section.WingCrossSection(foil,4,chord=1.5,spanwise_spacing='cosine',control_surface_symmetry_type='symmetric'),ps.geometry.wing_cross_section.WingCrossSection(foil,None,chord=.6,Lp_Wcsp_Lpp=[4*np.tan(np.deg2rad(35)),4,0],control_surface_symmetry_type='symmetric')];wing=ps.geometry.wing.Wing(sections,num_chordwise_panels=2,symmetric=True,symmetryNormal_G=[0,1,0],symmetryPoint_G_Cg=[0,0,0]);plane=ps.geometry.airplane.Airplane([wing],s_ref=8.4,c_ref=1.114285714285714,b_ref=8.);omega=8/1.114285714285714;period=2*np.pi/omega;sm=[ps.movements.wing_cross_section_movement.WingCrossSectionMovement(s) for s in sections];wm=ps.movements.wing_movement.WingMovement(wing,sm,ampAngles_Gs_to_Wn_ixyz=[0,.005,0],periodAngles_Gs_to_Wn_ixyz=[0,period,0],phaseAngles_Gs_to_Wn_ixyz=[0,90,0],rotationPointOffset_Gs_Ler=[.375,0,0]);pm=ps.movements.airplane_movement.AirplaneMovement(plane,[wm]);op=ps.operating_point.OperatingPoint(rho=1.2,vCg__E=20,alpha=0,nu=1e-12);om=ps.movements.operating_point_movement.OperatingPointMovement(op);movement=ps.movements.movement.Movement([pm],om,delta_time=period/samples,num_steps=samples*cycles+1);problem=ps.problems.UnsteadyProblem(movement,only_final_results=False)
    class Solver(ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver):
        def _collapse_geometry(self):
            super()._collapse_geometry();self._currentStackBoundRc0s[:]=1e-5;self._currentStackWakeRc0s[:]=1e-5
    solver=Solver(problem);start=time.monotonic();solver.run(prescribed_wake=True,calculate_streamlines=False,show_progress=False);seconds=time.monotonic()-start;F=np.array([p.airplanes[0].forces_W for p in problem.steady_problems]);M=np.array([p.airplanes[0].moments_W_CgP1 for p in problem.steady_problems]);normal=problem.steady_problems[0].airplanes[0].wings[0].panels[0,0].unitNormal_G;theta=float(np.arctan2(normal[0],normal[2]));assert abs(abs(theta)-np.deg2rad(.005))<1e-10;t=np.arange(len(F))*period/samples;last=t>=t[-1]-period-1e-10;D=np.column_stack([np.ones(last.sum()),np.cos(omega*t[last]),np.sin(omega*t[last])]);L=-F[:,2];My=M[:,1]+.375*L;fit=np.linalg.lstsq(D,np.column_stack([L[last],My[last]]),rcond=None)[0];transfer=(fit[1]-1j*fit[2])/theta
    return dict(F=F,M=M,t=t,omega=omega,signed_pitch_amplitude=theta,transfer=transfer,seconds=seconds)
def main():
    for samples,cycles in [(256,3),(512,3),(1024,3),(512,6)]:
        p=OUT/f'pitch-time-s{samples}-c{cycles}.json'
        if p.exists():continue
        print('START',p.stem,flush=True);r=run(samples,cycles);np.savez_compressed(p.with_suffix('.npz'),**r);meta=dict(samples=samples,cycles=cycles,seconds=r['seconds'],signed_pitch_amplitude=r['signed_pitch_amplitude'],transfer=[[float(z.real),float(z.imag)] for z in r['transfer']]);p.write_text(json.dumps(meta,indent=2),encoding='utf-8');print('DONE',p.stem,json.dumps(meta),flush=True);gc.collect()
if __name__=='__main__':main()
