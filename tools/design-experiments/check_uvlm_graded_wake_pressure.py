"""Predeclared graded W-P1 transient and conservative TE-pressure observations."""
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import lu_factor,lu_solve
import check_uvlm_bound_grading as grading

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()
D=grading.diagnostic
P=D.previous
MODELS={}
MESHES={}
ROWS=[]


def budget():
    if time.monotonic()-START>180 or P.memory.memory_bytes()>512*1024**2:
        raise RuntimeError('experiment_budget_exceeded')


def activate(m):
    if m not in MODELS:
        MESHES[m]=grading.install_grid(m,True)
        MODELS[m]=D.BASE[m]
    D.BASE.clear();D.BASE[m]=MODELS[m]
    return MESHES[m]


def weak_pressure(gamma,derivative,m,a):
    mesh=activate(m);x=np.array(mesh['vortex_x'])
    overlap=np.maximum(0.,x[1:]-np.maximum(x[:-1],a))
    cell=int(np.searchsorted(x,a,side='right')-1)
    left=gamma[cell] if cell>=0 else 0.
    integral=overlap@gamma
    temporal=1.2*derivative*integral/(1-a)
    flux=1.2*8*(gamma[-1]-left)/(1-a)
    return temporal+flux,temporal,flux,dict(a=float(a),width=float(1-a),left_cell=cell,
        overlapping_cells=int(np.count_nonzero(overlap)),left_trace='right',TE_trace='left')


def observation(gamma,derivative,m,k,a):
    total,temporal,flux,metadata=weak_pressure(gamma,derivative,m,a)
    scale=1.2*8*(16*k)*.001
    return dict(**metadata,pressure=P.pair(total),temporal=P.pair(temporal),flux=P.pair(flux),
        reference_pressure=scale,normalized_amplitude=float(abs(total)/scale))


def preliminary():
    q,w=leggauss(64);q=(q+1)/2;w=w/2;quadrature=[];manufactured=[]
    for m in [64,128,256]:
        mesh=activate(m);xc=mesh['collocation_x'][-1]
        for s in [1.,81.]:
            for delta in [.025,.003125]:
                r0=s-xc;ell=np.log1p(delta/r0);r=r0*np.exp(ell*q)
                tau=(r-r0)/delta
                fds=5000./(2*np.pi*r*r*np.sqrt(r*r+5000.**2))*r*ell
                a,b=D.weights(s,delta,xc,5000.)
                aq=np.dot(w,fds*(1-tau));bq=np.dot(w,fds*tau)
                eligible=lambda v,ref:bool(abs(v-ref)<=1e-10 or abs(v-ref)/abs(ref)<=1e-9)
                quadrature.append(dict(m=m,s=s,delta=delta,left_error=float(a-aq),right_error=float(b-bq),
                    eligible=eligible(a,aq) and eligible(b,bq)))
        derivative=(1-np.exp(-1j*8*(.00625/8)))/(.00625/8);a=1-1/np.sqrt(m)
        constant=np.full(m,1+2j);total,_,_,_=weak_pressure(constant,derivative,m,a)
        exact=1.2*derivative*(1+2j)
        split=m-4;xsplit=mesh['vortex_x'][split]
        jump=np.ones(m,complex);jump[split:]=2.
        jtotal,jtime,jflux,_=weak_pressure(jump,derivative,m,a)
        expected=1.2*(derivative*((xsplit-a)+2*(1-xsplit))+8)/(1-a)
        manufactured.append(dict(m=m,constant_relative=float(abs(total-exact)/abs(exact)),
            jump_relative=float(abs(jtotal-expected)/abs(expected)),
            missing_flux_difference=float(abs(jtotal-jtime)/abs(jtotal)),
            forced_zero_difference=float(abs(0.-exact)/abs(exact)),
            constant_pressure=P.pair(total),jump_pressure=P.pair(jtotal)))
    return quadrature,manufactured


def frequency(m,k,delta,length):
    budget();activate(m)
    ops=D.discrete(m,k,delta,length)
    row,gamma,force=D.solve(m,k,ops['p1'],delta,'graded-W-P1','BDF1',ops['effective_length'])
    row['requested_length']=length
    derivative=(1-np.exp(-2j*k*delta))/(delta/8)
    row['weak_pressure']=observation(gamma,derivative,m,k,1-1/np.sqrt(m))
    return row,gamma,force


def march(steps_per_period):
    m=128;k=.5;delta=np.pi/(k*steps_per_period);activate(m)
    ops=D.discrete(m,k,delta,40.,keep_nodes=True)
    row,gamma,force=D.solve(m,k,ops['p1'],delta,'graded-W-P1','BDF1',ops['effective_length'])
    a,fc,fu,_=D.base(m);matrix=a.copy();matrix[:,-1]+=ops['left0']
    factor=lu_factor(matrix);g=ops['nodes'];z=np.exp(2j*k*delta)
    history=(z**(-np.arange(ops['n']+1))*gamma[-1]).real
    old=gamma.real.copy();dt=delta/8;window=1-1/np.sqrt(m)
    derivative=(1-1/z)/dt;weak=weak_pressure(gamma,derivative,m,window)[0]
    pscale=1.2*8*8*.001;ref=P.reference(k)[0]
    max_discrete=0.;max_physical=0.;max_pressure=0.;max_residual=0.;trace=[];loads=[]
    for n in range(1,2*steps_per_period+1):
        if n%128==0:budget()
        trial=np.empty_like(history);trial[1:]=history[:-1]
        velocity=float(np.real(1j*8*.001*z**n))
        rhs=np.full(m,velocity)-g[:,1:]@trial[1:]
        current=lu_solve(factor,rhs);trial[0]=current[-1]
        dot=(current-old)/dt;load=float(fc@current+fu@dot)
        pressure=float(weak_pressure(dot,1.,m,window)[1]+weak_pressure(current,0.,m,window)[2])
        expected=float(np.real(force*z**n));ptarget=float(np.real(weak*z**n))
        max_discrete=max(max_discrete,abs(load-expected)/abs(force))
        max_physical=max(max_physical,abs(load-float(np.real(ref*z**n)))/abs(ref))
        max_pressure=max(max_pressure,abs(pressure-ptarget)/pscale)
        residual=a@current+g@trial-velocity
        max_residual=max(max_residual,float(np.linalg.norm(residual)/max(np.linalg.norm(rhs),1e-12)))
        loads.append(load)
        if n%(steps_per_period//8)==0:trace.append(dict(step=n,force=load,target=expected,
            pressure=pressure,pressure_target=ptarget))
        history=trial;old=current
    result=dict(m=m,k=k,steps_per_period=steps_per_period,periods=2,delta=delta,
        requested_length=40.,effective_length=ops['effective_length'],nwake=ops['n'],
        initialization='periodic nodal history',retired_cells=2*steps_per_period,
        max_discrete_force_difference=max_discrete,max_theodorsen_difference=max_physical,
        max_discrete_pressure_difference=max_pressure,max_nonpenetration_relative=max_residual,
        trace=trace,frequency_case=row,frequency_pressure=observation(gamma,derivative,m,k,window))
    return result,np.array(loads)


def main():
    quadrature,manufactured=preliminary();baseline=[];continuous_pressure=[];fixed_windows=[]
    initial=json.loads((ROOT/'docs/verification/uvlm-bound-grading-result.json').read_text())
    refined=json.loads((ROOT/'docs/verification/uvlm-bound-grading-refinement-result.json').read_text())
    saved=initial['harmonic_cases']+refined['harmonic_cases'];changes=[]
    for k in [.1,.5,1.]:
        cases={}
        for m in [64,128,256]:
            activate(m);kw,_=D.continuous(m,k,160.)
            row,gamma,_=D.solve(m,k,kw)
            old=next(v for v in saved if v['m']==m and v['k']==k)
            baseline.append(dict(m=m,k=k,relative=P.relative(row,old)))
            continuous_pressure.append(dict(m=m,k=k,observation=observation(gamma,1j*16*k,m,k,1-1/np.sqrt(m))))
            current,_,_=frequency(m,k,.003125,160.);ROWS.append(current);cases[(m,.003125,160.)]=current
        for delta in [.025,.0125,.00625]:
            current,_,_=frequency(256,k,delta,160.);ROWS.append(current);cases[(256,delta,160.)]=current
        for length in [80.,320.]:
            current,_,_=frequency(256,k,.003125,length);ROWS.append(current);cases[(256,.003125,length)]=current
        fine,gamma,_=frequency(256,k,.003125,160.)
        for width in [.1,.05,.025]:
            derivative=(1-np.exp(-2j*k*.003125))/(.003125/8)
            fixed_windows.append(dict(k=k,observation=observation(gamma,derivative,256,k,1-width)))
        changes.append(dict(k=k,spatial=P.relative(cases[(128,.003125,160.)],fine),
            temporal=P.relative(cases[(256,.00625,160.)],fine),
            tail_80_160=P.relative(cases[(256,.003125,80.)],fine),
            tail_160_320=P.relative(fine,cases[(256,.003125,320.)])))
        print(json.dumps(dict(stage='frequency_done',k=k,changes=changes[-1])),flush=True)
    runs=[];loads={}
    for steps in [256,512,1024]:
        result,values=march(steps);runs.append(result);loads[steps]=values
        print(json.dumps(dict(stage='time_done',steps=steps,discrete=result['max_discrete_force_difference'],
            physical=result['max_theodorsen_difference'])),flush=True)
    time_change=float(np.max(np.abs(loads[512]-loads[1024][1::2]))/abs(complex(*runs[-1]['frequency_case']['force'])))
    fine=[r for r in ROWS if r['m']==256 and r['delta']==.003125 and r['requested_length']==160.]
    checks=dict(near_TE_weights=all(r['eligible'] for r in quadrature),
        manufactured_pressure=all(r['constant_relative']<=1e-11 and r['jump_relative']<=1e-11 for r in manufactured),
        negative_missing_flux=all(r['missing_flux_difference']>.001 for r in manufactured),
        negative_forced_zero=all(r['forced_zero_difference']>.001 for r in manufactured),
        continuous_saved_match=all(r['relative']<=1e-10 for r in baseline),
        complex_accuracy=all(r['complex_relative']<=.02 for r in fine),
        amplitude_accuracy=all(r['amplitude_relative']<=.02 for r in fine),
        phase_accuracy=all(abs(r['phase_degrees'])<=1 for r in fine),
        spatial_change=all(r['spatial']<=.005 for r in changes),
        temporal_change=all(r['temporal']<=.01 for r in changes),
        tail_change=all(r['tail_160_320']<=.002 for r in changes),
        frequency_nonpenetration=all(r['nonpenetration_relative']<=1e-11 for r in ROWS),
        time_force_match=all(r['max_discrete_force_difference']<=1e-10 for r in runs),
        time_pressure_match=all(r['max_discrete_pressure_difference']<=1e-10 for r in runs),
        time_nonpenetration=all(r['max_nonpenetration_relative']<=1e-11 for r in runs),
        time_physical_accuracy=runs[-1]['max_theodorsen_difference']<=.02,time_refinement=time_change<=.01)
    source=Path(__file__);deps=['check_uvlm_bound_grading.py',*initial['dependency_sha256']]
    result=dict(protocol='docs/verification/uvlm-graded-wake-pressure-experiment.md',geometry_id=P.geometry.GEOMETRY,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={name:hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() for name in deps},
        baseline_result_sha256={name:hashlib.sha256((ROOT/'docs/verification'/name).read_bytes()).hexdigest() for name in ['uvlm-bound-grading-result.json','uvlm-bound-grading-refinement-result.json']},
        U=8.,rho=1.2,H=.001,c=1.,span=10000.,numpy=np.__version__,scipy=D.scipy.__version__,
        meshes=MESHES,quadrature_cases=quadrature,manufactured_pressure_cases=manufactured,
        continuous_baseline=baseline,continuous_weak_pressure=continuous_pressure,
        harmonic_cases=ROWS,refinement_changes=changes,fixed_pressure_windows=fixed_windows,
        time_runs=runs,time_refinement_difference=time_change,checks=checks,
        all_criteria_pass=bool(all(checks.values())),
        pressure_physical_status='diagnostic only; no pointwise Kutta or local pressure acceptance',
        seconds=time.monotonic()-START,peak_process_working_bytes=P.memory.memory_bytes(),
        limits='linear periodic heave, one span, prescribed wake, weak pressure diagnostic; no pitch/moment/free wake/ANCF coupling')
    (ROOT/'docs/verification/uvlm-graded-wake-pressure-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=P.core.json_scalar)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,changes=changes,time_change=time_change,seconds=result['seconds']),default=P.core.json_scalar),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':
    try:main()
    except Exception as error:
        (ROOT/'build/analysis/uvlm-graded-wake-pressure-stopped.json').write_text(json.dumps(dict(
            error=repr(error),partial_harmonic_cases=ROWS,seconds=time.monotonic()-START),indent=2,allow_nan=False,default=P.core.json_scalar)+'\n',encoding='utf-8')
        raise
