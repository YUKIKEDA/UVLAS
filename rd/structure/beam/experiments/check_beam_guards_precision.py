import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
from flint import arb_mat,arb,ctx
import beam_admissibility as guards
import check_beam_se3_physical_v4 as physics
ROOT=rd_paths.ROOT
def main():
    dest=ROOT/'rd/structure/beam/records/beam-guards-precision-result.json';assert not dest.exists();start=time.monotonic();K,M,props=physics.section(.1)
    req=dict(units='SI',order=guards.inputs.ORDER,source_id='synthetic',reference_scale=1.,warping='saint_venant_free',K=K.tolist(),initial_generalized_strain=[0.]*6,
        mass_per_length=M[0,0],centroid=[0.,0.,0.],centroid_inertia_per_length=M[3:,3:].tolist(),section_envelope=[.05,.1])
    section=guards.prepare(req);r=np.array([[0.,0,0],[1.,0,0]]);R=np.array([np.eye(3)]*2);proof=guards.regularity(r,R,1.,section);rejected=[]
    for name,rr,RR in [('collapse',np.zeros((2,3)),R),('inverted',np.array([[0.,0,0],[-1.,0,0]]),R)]:
        try:guards.regularity(rr,RR,1.,section)
        except ValueError:rejected.append(name)
        else:raise AssertionError(name)
    x,Q=physics.b.pose_at(np.zeros(3),np.eye(3),np.array([1.,0,0,0,0,30.]),.05);curved_r=np.array([np.zeros(3),physics.b.float_array(x)]);curved_R=np.array([np.eye(3),physics.b.float_array(Q)])
    try:guards.regularity(curved_r,curved_R,.05,section)
    except ValueError:rejected.append('tube')
    else:raise AssertionError('accepted folded/uncertified tube')
    for bounds in [None,[0.,.1],[-.1,.1],[float('nan'),.1]]:
        bad=dict(req,section_envelope=bounds)
        try:guards.prepare(bad)
        except ValueError:rejected.append('bounds_'+str(bounds))
        else:raise AssertionError('bounds')
    assert len(rejected)==7
    parent=ROOT/'rd/structure/beam/records/beam-se3-physical-v4-result.json';report=json.loads(parent.read_text());assert report['passed']
    linear=next(r for r in report['records'] if r['case']=='Timoshenko_tip_force_convergence')['cases'];replayed=[]
    for case in linear:
        state=case['state_precision']
        if state['precision_bits']==53:continue
        ctx.prec=256;root,_,_=physics.linear_operators(case['n'],case['h']);K,_,_=physics.section(case['h']);F=np.zeros(6*case['n']);F[-5]=K[5,5]*1e-6
        A=arb_mat(root.tolist());f=arb_mat([[float(x)] for x in F]);u=arb_mat([[arb(x)] for x in state['coordinates']]);res=A.transpose()*(A*u)-f
        retained=max(float(abs(res[i,0]).upper()) for i in range(len(F)))/max(abs(F));projected=arb_mat([[float(u[i,0])] for i in range(len(F))]);pr=A.transpose()*(A*projected)-f
        projection=max(float(abs(pr[i,0]).upper()) for i in range(len(F)))/max(abs(F));assert retained<=1e-8
        replayed.append(dict(n=case['n'],h=case['h'],retained_residual=retained,binary64_projection_residual=projection,projection_qualified=projection<=1e-8))
    assert len(replayed)==3 and any(not r['projection_qualified'] for r in replayed)
    paths=[Path(__file__),Path(guards.__file__),Path(guards.inputs.__file__),Path(physics.__file__),parent]
    result=dict(passed=True,records=[dict(case='reference_and_current_section_tube',certificate=proof,negative_controls=rejected),dict(case='retained_state_and_projection_precision',cases=replayed)],
        seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2,default=lambda x:x.item()),encoding='utf-8');print(json.dumps(result['records'],indent=2,default=lambda x:x.item()))
if __name__=='__main__':main()
