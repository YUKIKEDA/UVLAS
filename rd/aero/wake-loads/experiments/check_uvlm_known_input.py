import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,numpy as np
from uvlm_known_input import KnownFlow
import uvlm_wake_assembly_v2 as aero
import uvlm_surface_bindings as binding
ROOT=rd_paths.ROOT
def main():
    times=[0,.4,1.,2.];values=np.array([[20.,0,.1,.03],[21,0,.2,.06],[19,0,.4,.08],[20,0,.3,.1]]);flow=KnownFlow(times,values,wave=[.4,.2,0],direction=[0,0,1],omega=.7,phase=.2);rng=np.random.default_rng(110507);p=rng.normal(size=(8,3));dp=rng.normal(size=(8,3));t=.73;U,J,dt=flow.evaluate(p,t);h=1e-5;fd=(flow.value(p+h*dp,t)-flow.value(p-h*dp,t))/(2*h);jerr=np.linalg.norm(fd-flow.jvp(p,t,dp))/np.linalg.norm(fd);ft=(flow.value(p,t+h)-flow.value(p,t-h))/(2*h);terr=np.linalg.norm(ft-dt)/np.linalg.norm(ft);assert jerr<1e-7 and terr<1e-7;assert np.max(abs(np.trace(J,axis1=-2,axis2=-1)))<1e-15;records=[dict(case='space_time_derivatives_and_divergence',spatial_relative=float(jerr),temporal_relative=float(terr))]
    for knot in times[1:-1]:
        a=flow.evaluate(p,knot-1e-8);b=flow.evaluate(p,knot+1e-8);assert np.max(abs(a[0]-b[0]))<1e-6 and np.max(abs(a[2]-b[2]))<1e-5
    for j in range(len(times)-1):
        sample=flow.polynomial(np.linspace(times[j],times[j+1],51));assert np.all(sample>=np.minimum(values[j],values[j+1])-1e-14) and np.all(sample<=np.maximum(values[j],values[j+1])+1e-14)
    records.append(dict(case='C1_knots_and_component_monotonicity'))
    encoded=flow.encode();restored=KnownFlow.restore(encoded);assert restored.identity==flow.identity and restored.encode()==encoded
    for x,y in zip(flow.evaluate(p,t),restored.evaluate(p,t)):assert np.array_equal(x,y)
    records.append(dict(case='coefficient_restart_bitwise'))
    rejected=[]
    badsaved=json.loads(encoded);badsaved['coefficients'][3][1][0]+=.1
    for label,call in [('range',lambda:flow.value(p,3)),('duplicate',lambda:KnownFlow([0,0],values[:2])),('missing',lambda:KnownFlow([0,1],[[1,0,0,0],[2,0,float('nan'),0]])),('gust_divergence',lambda:KnownFlow(times,values,wave=[0,0,1])),('nonfinite_time',lambda:flow.value(p,float('nan'))),('saved_discontinuity',lambda:KnownFlow.restore(json.dumps(badsaved)))]:
        try:call()
        except ValueError:rejected.append(label)
        else:raise AssertionError(label)
    try:flow.omega=1
    except AttributeError:rejected.append('mutation')
    else:raise AssertionError('mutation')
    records.append(dict(case='invalid_input_rejections',rejected=rejected))
    b=binding.q1('wing',1,1,[[0,-.5,0],[0,.5,0],[1,.5,0],[1,-.5,0]]);model=aero.Model([b],.025);s=aero.cold_state(model,b['q0']);trial=aero.begin_trial(model,s,s.q,s.v,.1,flow,mode='prescribed');aero.solve_trial(model,s,trial,flow,aero.Workspace());state=trial['state'];saved=aero.checkpoint(model,state);r=aero.restore(model,saved,restored);out=[]
    for st,fl in [(state,flow),(r,restored)]:
        tr=aero.begin_trial(model,st,st.q,st.v,.1,fl);aero.solve_trial(model,st,tr,fl,aero.Workspace());out.append(tr)
    assert out[0]['state'].fingerprint()==out[1]['state'].fingerprint() and np.array_equal(out[0]['loads']['Q'],out[1]['loads']['Q']);changed=KnownFlow(times,values,wave=[.4,.2,0],direction=[0,0,1],omega=.701,phase=.2)
    try:aero.restore(model,saved,changed)
    except ValueError:pass
    else:raise AssertionError('changed_input_restart')
    records.append(dict(case='adopted_assembly_input_restart_and_identity',next_state_and_load_bitwise=True))
    paths=[Path(__file__),ROOT/'rd/aero/wake-loads/experiments/uvlm_known_input.py',ROOT/'docs/specs/aero/uvlm-integrated-specification.md',ROOT/'rd/aero/wake-loads/records/uvlm-wake-policy-result.json'];dest=ROOT/'rd/aero/wake-loads/records/uvlm-known-input-result.json';assert not dest.exists();result=dict(passed=True,groups=len(records),records=records,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
