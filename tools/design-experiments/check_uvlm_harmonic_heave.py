"""Linear B1 heave trial; prewritten protocol in docs/verification/."""
import hashlib
import json
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import lu_factor, lu_solve
from scipy.special import hankel2, jv, yv
import check_aerogrid_geometry as geometry
import check_vortex_core_load as core
import check_q9_plate_convergence as memory

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()
U=8.;RHO=1.2;H=.001
CACHE={}


def budget():
    if time.monotonic()-START>180 or memory.memory_bytes()>1024**3:
        raise RuntimeError('experiment_budget_exceeded')


def pair(z):return [float(np.real(z)),float(np.imag(z))]


def reference(k):
    omega=2*k*U
    c=hankel2(1,k)/(hankel2(1,k)+1j*hankel2(0,k))
    j0,j1,y0,y1=jv(0,k),jv(1,k),yv(0,k),yv(1,k)
    den=(j1+y0)**2+(y1-j0)**2
    independent=(j1*(j1+y0)+y1*(y1-j0)-1j*(y1*y0+j1*j0))/den
    return RHO*np.pi*.25*omega**2*H-1j*omega*RHO*np.pi*U*c*H,c,abs(c-independent)


def operators(m,delta,length,span=10000.,eps=0.):
    key=(m,delta,length,span,eps)
    if key in CACHE:return CACHE[key]
    budget()
    grid=geometry.make_grid(m,1,span,wake=False)
    a=geometry.aic(grid,eps)
    nw=int(np.ceil(length/delta))
    wake=np.empty((m,nw))
    for start in range(0,nw,256):
        budget();ids=np.arange(start,min(start+256,nw))
        x=1.+ids*delta
        p=np.zeros((len(ids),4,3))
        p[:,:,0]=np.stack([x,x,x+delta,x+delta],axis=1)
        p[:,:,1]=np.array([-span/2,span/2,span/2,-span/2])[None,:]
        aa=p.reshape(-1,3);bb=np.roll(p,-1,axis=1).reshape(-1,3)
        v=core.kernel(grid['collocation'],aa,bb,eps)
        wake[:,ids]=v[:,:,2].reshape(m,len(ids),4).sum(axis=2)
    # Build both force operators from the same finite edges and triangles.
    ends=np.array([[grid['points'][i],grid['points'][j]] for i,j in grid['edges']])
    dl=ends[:,1]-ends[:,0];ids=grid['body_ids']
    circulating=RHO*U*(dl[ids,1]@grid['body'][ids])/span
    p=grid['points'][grid['triangles']]
    area=.5*np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])[:,2]
    unsteady=np.bincount(grid['tri_owner'],weights=RHO*area,minlength=m)/span
    # Bound memory while retaining reusable main comparison operators.
    result=(a,wake,circulating,unsteady)
    CACHE[key]=result
    if len(CACHE)>5:CACHE.pop(next(iter(CACHE)))
    return result


def harmonic(m,delta,length,k,span=10000.,eps=0.,handoff='old',wake_enabled=True):
    a,wake,fc,fu=operators(m,delta,length,span,eps)
    dt=delta/U;omega=2*k*U;z=np.exp(1j*omega*dt)
    phase=z**(-(np.arange(wake.shape[1])+(1 if handoff=='old' else 0)))
    kw=wake@phase if wake_enabled else np.zeros(m,complex)
    matrix=a.astype(complex);matrix[:,-1]+=kw
    rhs=np.full(m,1j*omega*H)
    gamma=np.linalg.solve(matrix,rhs)
    derivative=(1-1/z)/dt
    circ=fc@gamma;noncirc=derivative*(fu@gamma);force=circ+noncirc
    ref,c,reference_difference=reference(k)
    row=dict(m=m,delta=delta,length=length,span=span,eps=eps,k=k,nwake=wake.shape[1],
        force=pair(force),circulating=pair(circ),unsteady=pair(noncirc),reference=pair(ref),
        theodorsen=pair(c),reference_arithmetic_difference=float(reference_difference),
        complex_relative=float(abs(force-ref)/abs(ref)),amplitude_relative=float(abs(abs(force)/abs(ref)-1)),
        phase_degrees=float(np.angle(force/ref,deg=True)),
        nonpenetration_relative=float(np.linalg.norm(matrix@gamma-rhs)/np.linalg.norm(rhs)))
    return row,gamma,force


def relative(row_a,row_b):
    a=complex(*row_a['force']);b=complex(*row_b['force'])
    return float(abs(a-b)/abs(b))


def time_march():
    m=32;k=.5;steps=512;delta=np.pi/(k*steps);length=40.
    row,gamma,force=harmonic(m,delta,length,k)
    a,wake,fc,fu=operators(m,delta,length)
    factor=lu_factor(a);dt=delta/U;omega=2*k*U;z=np.exp(1j*omega*dt)
    # At t=0 the newest existing wake carries gamma_TE at t=-dt.
    periodic_w=(z**(-(np.arange(wake.shape[1])+1))*gamma[-1]).real
    outputs=[]
    for kind in ['periodic','cold']:
        old=gamma.real.copy() if kind=='periodic' else np.zeros(m)
        w=periodic_w.copy() if kind=='periodic' else np.zeros(wake.shape[1])
        errors=[];trace=[];residual=0.
        for n in range(1,2*steps+1):
            if n%128==0:budget()
            w[1:]=w[:-1].copy();w[0]=old[-1]
            velocity=float(np.real(1j*omega*H*z**n))
            rhs=np.full(m,velocity)-wake@w
            current=lu_solve(factor,rhs)
            dot=(current-old)/dt
            load=float(fc@current+fu@dot)
            exact=float(np.real(force*z**n))
            errors.append(abs(load-exact)/abs(force))
            residual=max(residual,float(np.linalg.norm(a@current-rhs)/max(np.linalg.norm(rhs),1e-15)))
            if n%64==0:trace.append(dict(step=n,force=load,periodic_target=exact))
            old=current
        sample=np.array([entry['force'] for entry in trace])
        outputs.append(dict(initialization=kind,max_periodic_difference=max(errors),
            last_period_max_difference=max(errors[-steps:]),nonpenetration_relative=residual,
            finite=bool(np.all(np.isfinite(sample))),trace=trace))
    return dict(m=m,k=k,steps_per_period=steps,delta=delta,length=length,
        nwake=wake.shape[1],retired_cells_per_run=2*steps,cycles=2,
        periodic_frequency_solution=row,runs=outputs)


def main():
    rows=[];differences=[];negative=[]
    for k in [.1,.5,1.]:
        cases={}
        for m,delta,length,span,eps in [
            (16,.025,80.,10000.,0.),(32,.025,80.,10000.,0.),(64,.025,80.,10000.,0.),
            (64,.1,80.,10000.,0.),(64,.05,80.,10000.,0.),(64,.0125,80.,10000.,0.),
            (64,.0125,160.,10000.,0.),(64,.0125,160.,20000.,0.),
            (64,.0125,160.,10000.,.001/64)]:
            row,_,_=harmonic(m,delta,length,k,span,eps)
            rows.append(row);cases[(m,delta,length,span,eps)]=row
            print(json.dumps(dict(stage='harmonic',**row)),flush=True)
        get=lambda m,d,l:cases[(m,d,l,10000.,0.)]
        fine=get(64,.0125,160.)
        differences.append(dict(k=k,spatial=relative(get(32,.025,80.),get(64,.025,80.)),
            temporal=relative(get(64,.025,80.),get(64,.0125,80.)),
            tail=relative(get(64,.0125,80.),fine),
            span=relative(fine,cases[(64,.0125,160.,20000.,0.)]),
            core=relative(fine,cases[(64,.0125,160.,10000.,.001/64)])))
        wrong,_,_=harmonic(64,.0125,160.,k,handoff='current')
        missing_wake,_,_=harmonic(64,.0125,160.,k,wake_enabled=False)
        force=complex(*fine['force']);circ=complex(*fine['circulating'])
        negative.append(dict(k=k,current_handoff=relative(wrong,fine),
            missing_unsteady=float(abs(circ-force)/abs(force)),missing_wake=relative(missing_wake,fine)))
    march=time_march();fine=[r for r in rows if r['m']==64 and r['delta']==.0125 and r['length']==160 and r['span']==10000 and r['eps']==0]
    checks=dict(complex_accuracy=all(r['complex_relative']<=.02 for r in fine),
        amplitude_accuracy=all(r['amplitude_relative']<=.02 for r in fine),
        phase_accuracy=all(abs(r['phase_degrees'])<=1 for r in fine),
        spatial_change=all(d['spatial']<=.005 for d in differences),
        temporal_change=all(d['temporal']<=.01 for d in differences),
        tail_change=all(d['tail']<=.002 for d in differences),
        span_change=all(d['span']<=.0005 for d in differences),
        core_change=all(d['core']<=.001 for d in differences),
        nonpenetration=all(r['nonpenetration_relative']<=1e-11 for r in rows),
        reference_arithmetic=all(r['reference_arithmetic_difference']<=1e-12 for r in rows),
        periodic_time_match=march['runs'][0]['max_periodic_difference']<=1e-10,
        negative_current=max(r['current_handoff'] for r in negative)>.001,
        negative_unsteady=max(r['missing_unsteady'] for r in negative)>.001,
        negative_wake=max(r['missing_wake'] for r in negative)>.001)
    source=Path(__file__)
    dependencies=['check_aerogrid_geometry.py','check_vortex_core_load.py','check_q9_plate_convergence.py',
        'check_q9_interface_reparameterization.py','check_q9_curved_load.py','check_q9_independent_load.py']
    result=dict(protocol='docs/verification/uvlm-harmonic-heave-experiment.md',geometry_id=geometry.GEOMETRY,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={name:hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() for name in dependencies},
        numpy=np.__version__,scipy=scipy.__version__,U=U,rho=RHO,H=H,c=1.,
        harmonic_cases=rows,refinement_changes=differences,negative_controls=negative,time_march=march,
        checks=checks,all_criteria_pass=bool(all(checks.values())),seconds=time.monotonic()-START,
        peak_process_working_bytes=memory.memory_bytes(),
        limits='linear heave, prescribed advection, finite closed wake; no free wake, ANCF coupling or moment accuracy')
    output=ROOT/'docs/verification/uvlm-harmonic-heave-result.json'
    output.write_text(json.dumps(result,indent=2,allow_nan=False,default=core.json_scalar)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,refinement=differences,time_march=march['runs'],seconds=result['seconds']),default=core.json_scalar),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':main()
