"""Unmodified PteraSoftware probe; no reference qualification is implied."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np,json,time,importlib.metadata as md
import pterasoftware as ps
ROOT=rd_paths.ROOT
def run(m=2,n=8,steps_per_cycle=64,cycles=3,amplitude=1e-4,k=.2,core=None,nu=1e-12):
    airfoil=ps.geometry.airfoil.Airfoil(name='NACA0012')
    sections=[ps.geometry.wing_cross_section.WingCrossSection(airfoil=airfoil,num_spanwise_panels=n//2,chord=1.5,spanwise_spacing='cosine',control_surface_symmetry_type='symmetric'),ps.geometry.wing_cross_section.WingCrossSection(airfoil=airfoil,num_spanwise_panels=None,chord=.6,Lp_Wcsp_Lpp=[4*np.tan(np.deg2rad(35)),4,0],control_surface_symmetry_type='symmetric')]
    wing=ps.geometry.wing.Wing(wing_cross_sections=sections,num_chordwise_panels=m,symmetric=True,symmetryNormal_G=[0,1,0],symmetryPoint_G_Cg=[0,0,0])
    airplane=ps.geometry.airplane.Airplane(wings=[wing],s_ref=8.4,c_ref=1.114285714285714,b_ref=8.)
    omega=2*k*20/1.114285714285714;period=2*np.pi/omega
    section_moves=[ps.movements.wing_cross_section_movement.WingCrossSectionMovement(s) for s in sections]
    wing_move=ps.movements.wing_movement.WingMovement(wing,section_moves,ampLer_Gs_Cgs=[0,0,amplitude],periodLer_Gs_Cgs=[0,0,period],phaseLer_Gs_Cgs=[0,0,90])
    plane_move=ps.movements.airplane_movement.AirplaneMovement(airplane,[wing_move]);op=ps.operating_point.OperatingPoint(rho=1.2,vCg__E=20.,alpha=0.,nu=nu);op_move=ps.movements.operating_point_movement.OperatingPointMovement(op)
    movement=ps.movements.movement.Movement([plane_move],op_move,delta_time=period/steps_per_cycle,num_steps=cycles*steps_per_cycle+1)
    problem=ps.problems.UnsteadyProblem(movement,only_final_results=False)
    cls=ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    if core is not None:
        class PrecisionCoreSolver(cls):
            def _collapse_geometry(self):
                super()._collapse_geometry()
                self._currentStackBoundRc0s[:]=core;self._currentStackWakeRc0s[:]=core
        cls=PrecisionCoreSolver
    solver=cls(problem);start=time.monotonic();solver.run(prescribed_wake=True,calculate_streamlines=False,show_progress=False);elapsed=time.monotonic()-start
    F=np.array([p.airplanes[0].forces_W for p in problem.steady_problems]);M=np.array([p.airplanes[0].moments_W_CgP1 for p in problem.steady_problems]);t=np.arange(len(F))*movement.delta_time
    return solver,problem,dict(F=F,M=M,t=t,omega=omega,amplitude=amplitude,seconds=elapsed,dt=movement.delta_time)
def main():
    solver,problem,r=run();out=ROOT/'build/reference-ptera/probe.npz';np.savez_compressed(out,**r)
    last=r['t']>=r['t'][-1]-2*np.pi/r['omega'];D=np.column_stack([np.ones(last.sum()),np.cos(r['omega']*r['t'][last]),np.sin(r['omega']*r['t'][last])]);fit=np.linalg.lstsq(D,-r['F'][last,2],rcond=None)[0];coef=(fit[1]-1j*fit[2])/r['amplitude']
    print(json.dumps(dict(version=md.version('PteraSoftware'),steps=len(r['t']),seconds=r['seconds'],dt=r['dt'],lift_transfer=[coef.real,coef.imag],initial_force=r['F'][0].tolist(),final_force=r['F'][-1].tolist(),snapshot=out.relative_to(ROOT).as_posix()),indent=2))
if __name__=='__main__':main()
