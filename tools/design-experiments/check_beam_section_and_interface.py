"""Section validation and body-moment/point-force shell transfer contract."""
from pathlib import Path
import copy,hashlib,json,time,traceback
import numpy as np
from scipy.linalg import expm
import beam_se3_v2 as b
import beam_se3_high_precision as hp
import beam_decimal_reference as independent
import beam_section_input as section
from check_beam_se3_physical_v4 import section as iso
ROOT=Path(__file__).resolve().parents[2]
def main():
    dest=ROOT/'docs/verification/beam-section-interface-result.json';assert not dest.exists();records=[];failures=[];start=time.monotonic()
    try:
        K,M,props=iso(.1);m=M[0,0];J=M[3:,3:];c=np.array([0.,.035,-.021]);C=np.array(b.skew(c),float);S=np.block([[np.eye(3),-C],[np.zeros((3,3)),np.eye(3)]])
        req=dict(units='SI',order=section.ORDER,source_id='synthetic-rectangle-free-SV',reference_scale=1.,warping='saint_venant_free',K=(S.T@K@S).tolist(),
            initial_generalized_strain=[.001,-.002,.001,.003,-.001,.002],mass_per_length=m,centroid=c.tolist(),centroid_inertia_per_length=J.tolist())
        before=copy.deepcopy(req);p=section.prepare(req);assert req==before
        rejected=[]
        for key,value in [('units','mm_N'),('order','unknown'),('warping','restrained'),('restrained_warping',True),('source_id',''),('mass_per_length',0.),('centroid_inertia_per_length',np.diag([3.,1.,1.]).tolist()),('initial_generalized_strain',[0.]*5),('K',np.diag([1,1,-1,1,1,1]).tolist())]:
            r=copy.deepcopy(req);r[key]=value
            try:section.prepare(r)
            except ValueError:rejected.append(key)
            else:raise AssertionError('accepted '+key)
        assert len(rejected)==9
        records.append(dict(case='explicit_section_mass_and_warping_controls',negative_controls=rejected,prepared_immutable=not p['K'].flags.writeable))
        rng=np.random.default_rng(20261105);ref=np.array([1.,.02,-.01,.1,.2,-.15]);L=.8;d=L*(ref+np.array([.02,-.01,.01,.07,-.03,.02]));x,R=b.pose_at(np.zeros(3),np.eye(3),d,1.);r=np.array([[0.,0,0],b.float_array(x)]);RR=np.array([np.eye(3),b.float_array(R)])
        e=b.evaluate(r,RR,L,p['K'],ref,p['e0']);high=hp.evaluate(r,RR,L,p['K'],ref,p['e0']);direc=rng.normal(size=12);direc/=np.linalg.norm(direc)
        fd,sd=independent.directional(r,RR,L,p['K'],ref,p['e0'],direc)
        metrics=dict(force=abs(np.asarray(high['f'],float)@direc-fd)/abs(fd),tangent=abs(direc@np.asarray(high['K'],float)@direc-sd)/abs(sd))
        assert max(metrics.values())<=1e-8;records.append(dict(case='precision_backend_general_coupled_section_derivatives',metrics=metrics))
        point,Q=b.point_jet(r,RR,.43,[0.,.03,-.02]);value=b.float_array(point);Rpoint=b.float_array(Q);Jx=np.array([a.g for a in point]);JR=np.array([[a.g for a in row] for row in Q]);W=np.array([b.vee(Rpoint.T@JR[:,:,j]) for j in range(12)],float).T
        force=np.array([.7,-.4,.2]);moment=np.array([.3,.5,-.2]);nodal=Jx.T@force+W.T@moment;v=rng.normal(size=12)
        error=abs(nodal@v-force@(Jx@v)-moment@(W@v));assert error<=1e-10
        rigid=np.zeros((12,6))
        for i in range(2):rigid[6*i:6*i+3,:3]=np.eye(3);rigid[6*i:6*i+3,3:]=-np.array(b.skew(r[i]),float);rigid[6*i+3:6*i+6,3:]=RR[i].T
        expected=np.r_[force,np.cross(value,force)+Rpoint@moment];wrench_error=np.linalg.norm(rigid.T@nodal-expected)/np.linalg.norm(expected);assert wrench_error<=1e-10
        records.append(dict(case='shell_point_and_frame_dual_work',virtual_work_error=float(error),force_moment_wrench_error=float(wrench_error),position_H_available=True,frame_H_available=True))
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    paths=[Path(__file__),Path(b.__file__),Path(hp.__file__),Path(hp.ad.__file__),Path(section.__file__),Path(independent.__file__),ROOT/'tools/design-experiments/check_beam_se3_physical_v4.py']
    result=dict(passed=len(records)==3 and not failures,records=records,failures=failures,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result['records'],indent=2));print('SECTION_INTERFACE_RESULT',result['passed'])
if __name__=='__main__':main()
