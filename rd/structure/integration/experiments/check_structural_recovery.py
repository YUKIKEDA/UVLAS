import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np
import time,json,hashlib,traceback,sys
import structural_recovery as recovery
import material_equivalent_input as material
from check_material_shell_integration import requests,independent_Ceng
import check_beam_se3_physical_v4 as physical
import connection_shell as pg
ROOT=rd_paths.ROOT
def rel(a,b):return float(np.linalg.norm(a)/max(np.linalg.norm(b),1e-30))
def main():
    start=time.monotonic();rows=[];failures=[]
    try:
        stack=material.prepare_stack(requests());F=np.array([[1.07,.03,-.01],[.02,.96,.04],[.01,-.02,1.03]]);E=(F.T@F-np.eye(3))/2;errors=[];values=[]
        for layer,req in zip(stack,requests()):
            angle=np.arctan2(req['fiber_global'][1],req['fiber_global'][0]);c=np.cos(angle);s=np.sin(angle);A=np.array([[c,-s,0],[s,c,0],[0,0,1.]]);T=material.transform(A);C=T@layer.material.C@T.T;e0=T@layer.material.e0;out=recovery.stress(F,C,e0);Em=A.T@E@A-material.tensor(layer.material.e0);ee=np.array([Em[0,0],Em[1,1],Em[2,2],2*Em[1,2],2*Em[0,2],2*Em[0,1]]);ss=independent_Ceng(req['material'])@ee;Sm=np.array([[ss[0],ss[5],ss[4]],[ss[5],ss[1],ss[3]],[ss[4],ss[3],ss[2]]]);Sref=A@Sm@A.T;sigref=F@Sref@F.T/np.linalg.det(F);err=max(rel(out['second_PK']-Sref,Sref),rel(out['cauchy']-sigref,sigref));errors.append(err);values.append(out['second_PK']);assert err<1e-8
        assert np.linalg.norm(values[0]-values[1])>1e-4;rows.append(dict(case='layerwise_stress_measure_axes_and_sides',max_independent_stress_error=max(errors),interface_stresses_not_averaged=True))
        # Independent physical traction transform and moment transport on a
        # curved cut, using a constant-F exact affine displacement field.
        model=pg.flat_model(1,.1,1.,.2);p,w,area=pg.face(model,0,1.);N,G,X=pg.point(model,p);xx=X@F.T+np.array([.1,-.2,.3]);C=np.eye(6)*1000;out=recovery.stress(np.broadcast_to(F,(len(X),3,3)),C,np.zeros(6));origin=np.array([.2,.3,-.1]);cut=recovery.cut_wrench(xx,np.broadcast_to(F,(len(X),3,3)),out['first_PK'],area,w,origin);spatial=np.einsum('pab,pb->pa',out['cauchy'],cut['current_area']);nanson=rel(spatial-cut['reference_traction_density'],spatial);area_sum=np.sum(w[:,None]*area,axis=0);expectedF=out['first_PK'][0]@area_sum;centroid=np.average(xx,weights=w,axis=0);expectedM=np.cross(centroid-origin,expectedF);offset=np.array([.3,-.1,.2]);moved=recovery.cut_wrench(xx,np.broadcast_to(F,(len(X),3,3)),out['first_PK'],area,w,origin+offset);transport=rel(moved['moment']-cut['moment']+np.cross(offset,cut['force']),cut['moment']);err=max(nanson,rel(cut['force']-expectedF,expectedF),rel(cut['moment']-expectedM,expectedM),transport);assert err<1e-8;rows.append(dict(case='physical_cut_traction_and_wrench_transport',max_error=err))
        K,_,_=physical.section(.1);P=.36*K[5,5];cases=[]
        for n in [8,16,32]:
            r,R,history=physical.nonlinear(n,P,K);forces=[];moments=[];reference_force=np.array([0.,P,0.]);reference_moment=np.cross(r[-1],reference_force)
            for i in range(n):
                for t in [.23,.51,.79]:
                    cut=recovery.beam_cut(r[i:i+2],R[i:i+2],1/n,K,np.array([1.,0,0,0,0,0]),np.zeros(6),t,np.zeros(3));forces.append(rel(cut['force']-reference_force,reference_force));moments.append(rel(cut['moment']-reference_moment,reference_moment))
            row=dict(n=n,force_error=max(forces),moment_error=max(moments),original_residual=history[-1]);cases.append(row);print('CUT_CONVERGENCE',row,flush=True)
        assert max(cases[-1]['force_error'],cases[-1]['moment_error'])<.01;assert cases[-1]['force_error']<cases[0]['force_error'] and cases[-1]['moment_error']<cases[0]['moment_error'];rows.append(dict(case='nonlinear_beam_cut_wrench_vs_boundary_balance',cases=cases,no_equilibrium_correction=True))
        invalid=F.copy();invalid[:,0]*=-1
        try:recovery.stress(invalid,C,np.zeros(6))
        except ValueError:pass
        else:raise AssertionError('inverted recovery accepted')
        rows.append(dict(case='recovery_admissibility',inverted_state_rejected=True))
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    print(json.dumps(rows,indent=2),flush=True)
    if '--publish' in sys.argv:
        dest=ROOT/'rd/structure/integration/records/structural-recovery-result.json';assert not dest.exists();paths=[Path(__file__),Path(recovery.__file__),Path(material.__file__),Path(physical.__file__),Path(pg.__file__),ROOT/'rd/structure/material/experiments/check_material_shell_integration.py',ROOT/'rd/structure/integration/records/structural-final-integration-protocol.md'];dest.write_text(json.dumps(dict(passed=len(rows)==4 and not failures,records=rows,failures=failures,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8');print('STRUCTURAL_RECOVERY_RESULT',len(rows)==4 and not failures,flush=True)
if __name__=='__main__':main()
