from pathlib import Path
import numpy as np,json,hashlib,time
import uvlm_wake_assembly_v2 as aero
import uvlm_surface_bindings as binding
from uvlm_edge_load_integral import loads
ROOT=Path(__file__).resolve().parents[2]
def main():
    start=time.monotonic();b=binding.q1('wing',1,1,[[0.,-.4,0.],[.1,.4,.02],[1.1,.4,.02],[1.,-.4,0.]]);model=aero.Model([b],.08);q=b['q0'];v=np.zeros_like(q);flow=aero.Flow([1.,0,.05]);parent=aero.cold_state(model,q);t=aero.begin_trial(model,parent,q,v,.2,flow,mode='prescribed');aero.solve_trial(model,parent,t,flow);g=model.geometry(q,v);a=loads(model,g,t['wakes'],t['gamma'],t['gamma_dot'],flow,.2,32);fine=loads(model,g,t['wakes'],t['gamma'],t['gamma_dot'],flow,.2,64);err=float(np.linalg.norm(a['Q']-fine['Q'])/np.linalg.norm(fine['Q']));assert err<1e-8
    rng=np.random.default_rng(222);dq=rng.normal(size=q.shape)*.01;dv=rng.normal(size=q.shape)*.01;dU=np.array([.01,-.002,.005]);base=aero.trial_jvp(model,parent,t,flow,dq,dv,dU);dw=[np.zeros_like(w.vertices) for w in t['wakes']];dw[0][0]=model.Nv[model.surface_info[0]['te']]@dq;dw[0][1:]=.2*dU;dmu=np.zeros(t['wakes'][0].mu.size);dmu[:model.nte]=base['dgamma'][model.surface_info[0]['gamma_te']];dg=model.geometry(q,v,dq,dv);der=loads(model,dg,t['wakes'],t['gamma'],t['gamma_dot'],flow,.2,32,base['dgamma'],base['dgamma']/.2,dw,dmu,dU)
    def response(s):
        fl=aero.Flow(flow.U+s*dU);trial=aero.begin_trial(model,parent,q+s*dq,v+s*dv,.2,fl,mode='prescribed');aero.solve_trial(model,parent,trial,fl);gg=model.geometry(q+s*dq,v+s*dv);return loads(model,gg,trial['wakes'],trial['gamma'],trial['gamma_dot'],fl,.2,32)['Q']
    h=.001;numeric=(-response(2*h)+8*response(h)-8*response(-h)+response(-2*h))/(12*h);derr=float(np.linalg.norm(numeric-der['dQ'])/np.linalg.norm(numeric));assert derr<1e-8
    work=float(np.sum(a['Q']*dq));nodalwork=float(np.sum(a['Fv']*(model.Nv@dq)));assert abs(work-nodalwork)<1e-12;assert np.linalg.norm(a['F']-a['Fv'].sum(axis=0))<1e-12;assert np.linalg.norm(a['M']-np.cross(g['X'],a['Fv']).sum(axis=0))<1e-12
    sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'tools/design-experiments/uvlm_edge_load_integral.py',ROOT/'docs/verification/finite-wing-drag-diagnostic-protocol.md']};r=dict(passed=True,quadrature_32_64_error=err,full_implicit_load_tangent_error=derr,work_force_moment_preserved=True,adopted=False,sources=sources,seconds=time.monotonic()-start);dest=ROOT/'docs/verification/edge-load-integral-result.json';assert not dest.exists();dest.write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r,indent=2))
if __name__=='__main__':main()
