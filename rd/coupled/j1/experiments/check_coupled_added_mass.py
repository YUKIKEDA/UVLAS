"""Actual adopted spatial aero operators; constrained linear initialization witness."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,numpy as np
from scipy.linalg import solve
from uvlm_frequency_fast import geometry,bound
ROOT=rd_paths.ROOT
def main():
    g=geometry(4,24);X=g['X'];tail=X.copy();tail[-1,:,0]+=128.;A=bound(g['P'],g['normal'],tail,1e-5);quads=np.stack([X[:-1,:-1],X[:-1,1:],X[1:,1:],X[1:,:-1]],axis=2).reshape(-1,4,3);H=np.zeros((2,len(g['P'])))
    for ids in [[0,2,1],[0,3,2]]:
        tri=quads[:,ids];area=1.2*np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])/2;H[0]+=area[:,2];H[1]+=np.cross(tri.mean(axis=1)-[.375,0,0],area)[:,1]
    # +theta around +y gives vz=-(x-axis)*theta_dot. Same rigid velocity map in Ra.
    V=np.column_stack([np.ones(len(g['P'])),-(g['P'][:,0]-.375)]);B=-g['normal'][:,2,None]*V;M=np.diag([40.,8.]);force=np.array([1.3,-.4]);d=np.linspace(-.01,.02,len(A));G=np.array([[.2,1.]]);bc=np.array([.003]);z=np.zeros((len(A),1));K=np.block([[M,-H,G.T],[B,A,z],[G,np.zeros((1,len(A))),np.zeros((1,1))]]);rhs=np.r_[force,-d,bc];solution=solve(K,rhs);a=solution[:2];gd=solution[2:-1];lam=solution[-1:];Meff=M+H@solve(A,B);red=np.block([[Meff,G.T],[G,np.zeros((1,1))]]);xr=solve(red,np.r_[force-H@solve(A,d),bc]);difference=np.linalg.norm(solution[:2]-xr[:2]);res=np.linalg.norm(K@solution-rhs)/np.linalg.norm(rhs);assert difference<1e-12 and res<1e-11
    added_heave=float((Meff-M)[0,0]);assert added_heave>0;wrong=M-H@solve(A,B);wrong_solution=solve(np.block([[wrong,G.T],[G,np.zeros((1,1))]]),np.r_[force+H@solve(A,d),bc]);wrong_difference=np.linalg.norm(wrong_solution[:2]-a)/np.linalg.norm(a);assert wrong_difference>.05
    # Independently assemble physical force from triangle centroids for solved Gamma-dot.
    wrench=np.zeros(2)
    for p,q in enumerate(quads):
        for ids in [[0,2,1],[0,3,2]]:
            tri=q[ids];f=1.2*gd[p]*np.cross(tri[1]-tri[0],tri[2]-tri[0])/2;wrench+=np.array([f[2],np.cross(tri.mean(axis=0)-[.375,0,0],f)[1]])
    force_error=np.linalg.norm(wrench-H@gd);assert force_error<1e-12;constraint=float(np.linalg.norm(G@a-bc));assert constraint<1e-12
    result=dict(passed=True,scope='actual spatial heave/pitch operators and a linear acceleration constraint; not full transient/free-wake adoption',full_block_residual=float(res),schur_acceleration_difference=float(difference),added_heave_mass=added_heave,physical_noncirculatory_wrench_difference=float(force_error),acceleration_constraint_residual=constraint,wrong_sign_negative_control_relative=float(wrong_difference),sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'docs/specs/coupled/coupled-initialization-residual-contract.md',ROOT/'rd/aero/unsteady-reference/experiments/uvlm_frequency_fast.py']});dest=ROOT/'rd/coupled/j1/records/coupled-added-mass-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
