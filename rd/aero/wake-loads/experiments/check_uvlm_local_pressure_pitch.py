"""Independent acceleration-potential pressure; graded UVLM heave/pitch trial."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.special import kv
from scipy.linalg import lu_factor,lu_solve
from numpy.polynomial.legendre import leggauss
import check_uvlm_graded_wake_pressure as existing

ROOT=rd_paths.ROOT
START=time.monotonic()
D=existing.D
P=existing.P
MODES=[dict(name='heave-quarter',H=.001,alpha=0.,xp=.25),
       dict(name='pitch-quarter',H=0.,alpha=.001,xp=.25),
       dict(name='pitch-mid',H=0.,alpha=.001,xp=.5)]
MODELS={}
CONTINUOUS=[]
WAVE=[]


def budget():
    if time.monotonic()-START>180 or P.memory.memory_bytes()>512*1024**2:
        raise RuntimeError('experiment_budget_exceeded')


def scales(mode):
    p=1.2*8**2*(abs(mode['H'])/.5+abs(mode['alpha']))
    return p,.5*p,.25*p


def coefficients(k,mode):
    c=1.+0j if k==0 else kv(1,1j*k)/(kv(0,1j*k)+kv(1,1j*k))
    eta=mode['H']/.5;alpha=mode['alpha'];a=(mode['xp']-.5)/.5
    b0=-2j*k*c*eta+(2*c*(1+1j*k*(.5-a))-1j*k)*alpha
    b1=k*k*eta+(2j*k+k*k*a)*alpha;b2=-k*k*alpha/4
    return np.array([b0,b1,b2]),c


def reference_load(k,mode):
    b,c=coefficients(k,mode);a=(mode['xp']-.5)/.5
    return np.array([1.2*np.pi*8**2*.5*(b[0]+b[1]),
        1.2*np.pi*8**2*.5**2*((a+.5)*b[0]+a*b[1]-b[2]/2)])


def expanded_load(k,mode):
    c=coefficients(k,mode)[1];omega=16*k;b=.5;a=(mode['xp']-b)/b
    h=mode['H'];alpha=mode['alpha'];w=-1j*omega*h+(8+1j*omega*b*(.5-a))*alpha
    force=1.2*np.pi*b*b*(omega**2*h+(1j*omega*8+b*a*omega**2)*alpha)+2*np.pi*1.2*8*b*c*w
    moment=1.2*np.pi*b*b*(a*b*omega**2*h+
        (-1j*omega*8*b*(.5-a)+b*b*(.125+a*a)*omega**2)*alpha)+2*np.pi*1.2*8*b*b*(a+.5)*c*w
    return np.array([force,moment])


def pressure_integral(k,mode,lower,upper):
    b,_=coefficients(k,mode)
    left=np.arccos(2*lower-1);right=np.arccos(2*upper-1)
    r=.5*((left-np.sin(left))-(right-np.sin(right)))
    s=.5*((left/2-np.sin(2*left)/4)-(right/2-np.sin(2*right)/4))
    xs=.5*(np.sin(left)**3-np.sin(right)**3)/3
    return 1.2*8**2*(b[0]*r+2*b[1]*s+4*b[2]*xs)


def reference_tests():
    q,w=leggauss(80);theta=(q+1)*np.pi/2;wt=w*np.pi/2;rows=[]
    for k in [.1,.5,1.]:
        for mode in MODES:
            b,c=coefficients(k,mode);x=(1+np.cos(theta))/2
            pdx=1.2*8**2*.5*(b[0]*(1-np.cos(theta))+2*b[1]*np.sin(theta)**2+
                4*b[2]*np.cos(theta)*np.sin(theta)**2)
            integral=np.array([wt@pdx,wt@(-(x-mode['xp'])*pdx)])
            closed=reference_load(k,mode);expanded=expanded_load(k,mode)
            _,fs,ms=scales(mode);den=np.maximum(np.abs(closed),.01*np.array([fs,ms]))
            limit=(2/3)*1.2*8**2*(b[0]+4*b[1]+8*b[2]);limits=[]
            for width in [.01,.001,.0001,.00001]:
                average=pressure_integral(k,mode,1-width,1)/width
                limits.append(dict(width=width,average=P.pair(average),
                    divided_by_sqrt=P.pair(average/np.sqrt(width))))
            last=complex(*limits[-1]['divided_by_sqrt'])
            rows.append(dict(k=k,mode=mode['name'],coefficients=[P.pair(v) for v in b],
                C=P.pair(c),C_hankel_difference=float(abs(c-P.reference(k)[1])),
                loads=[P.pair(v) for v in closed],quadrature_error=float(np.max(np.abs(integral-closed)/den)),
                expanded_error=float(np.max(np.abs(expanded-closed)/den)),
                trailing_limit=P.pair(limit),limit_relative=float(abs(last-limit)/max(abs(limit),scales(mode)[0])),
                trailing_means=limits,TE_point_limit=[0.,0.]))
    return rows


def model(m):
    mesh=existing.activate(m)
    if m not in MODELS:
        grid=P.geometry.make_grid(m,1,10000.,wake=False)
        grid['points'][:,0]=np.repeat(np.array(mesh['vortex_x']),2)
        ends=np.array([[grid['points'][i],grid['points'][j]] for i,j in grid['edges']])
        dl=ends[:,1]-ends[:,0];centers=ends.mean(axis=1);ids=grid['body_ids']
        mc0=-1.2*8*((centers[ids,0]*dl[ids,1])@grid['body'][ids])/10000.
        tri=grid['points'][grid['triangles']]
        area=.5*np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])[:,2]
        mu0=np.bincount(grid['tri_owner'],weights=-1.2*area*tri.mean(axis=1)[:,0],minlength=m)/10000.
        MODELS[m]=(mc0,mu0)
    return mesh,D.base(m),MODELS[m]


def weak_average(gamma,derivative,m,lower,upper):
    mesh=existing.activate(m);x=np.array(mesh['vortex_x'])
    overlap=np.maximum(0.,np.minimum(x[1:],upper)-np.maximum(x[:-1],lower))
    li=int(np.searchsorted(x,lower,side='right')-1)
    ui=int(np.searchsorted(x,upper,side='left')-1)
    left=gamma[li] if li>=0 else 0.;right=gamma[ui] if ui>=0 else 0.
    return 1.2*(derivative*(overlap@gamma)+8*(right-left))/(upper-lower)


def pressure_profile(gamma,derivative,m,k,mode):
    edges=[.05,.15,.3,.5,.7,.85,.95,1.]
    intervals=list(zip(edges[:-1],edges[1:]))+[(1-width,1.) for width in [.1,.05,.025,.0125]]
    rows=[];ps=scales(mode)[0]
    for lower,upper in intervals:
        value=weak_average(gamma,derivative,m,lower,upper)
        ref=pressure_integral(k,mode,lower,upper)/(upper-lower)
        rows.append(dict(lower=lower,upper=upper,pressure=P.pair(value),reference=P.pair(ref),
            normalized_error=float(abs(value-ref)/max(ps,abs(ref)))))
    width=1/np.sqrt(m);value=weak_average(gamma,derivative,m,1-width,1.)
    b,_=coefficients(k,mode);limit=(2/3)*1.2*8**2*(b[0]+4*b[1]+8*b[2])
    finite=pressure_integral(k,mode,1-width,1)/width
    return rows,dict(width=width,pressure=P.pair(value),reference_mean=P.pair(finite),
        divided_by_sqrt=P.pair(value/np.sqrt(width)),asymptotic=P.pair(limit),
        normalized_asymptotic_difference=float(abs(value/np.sqrt(width)-limit)/max(ps,abs(limit))),
        normalized_finite_window_difference=float(abs(value-finite)/max(ps,abs(finite))))


def solve_cases(m,k,kw,delta=None,modes=MODES):
    budget();mesh,(a,fc,fu,xc),(mc0,mu0)=model(m)
    matrix=a.astype(complex);matrix[:,-1]+=kw
    rhs=np.column_stack([1j*16*k*mode['H']-(8+1j*16*k*(xc-mode['xp']))*mode['alpha'] for mode in modes])
    gamma=np.linalg.solve(matrix,rhs)
    derivative=1j*16*k if delta is None else (1-np.exp(-2j*k*delta))/(delta/8)
    rows=[]
    for j,mode in enumerate(modes):
        g=gamma[:,j];xp=mode['xp'];mc=mc0+xp*fc;mu=mu0+xp*fu
        load=np.array([fc@g+derivative*(fu@g),mc@g+derivative*(mu@g)])
        ref=reference_load(k,mode);_,fs,ms=scales(mode);qscale=np.array([fs,ms])
        den=np.maximum(np.abs(ref),.01*qscale)
        x=np.array(mesh['vortex_x']);width=np.diff(x)
        first=.5*((x[1:]-xp)**2-(x[:-1]-xp)**2)
        conservative=-1.2*(derivative*(first@g)+8*((1-xp)*g[-1]-width@g))
        shifted=(mc+.1*fc)@g+derivative*((mu+.1*fu)@g)
        profile,near=pressure_profile(g,derivative,m,k,mode)
        phase=[float(np.angle(load[l]/ref[l],deg=True)) if abs(ref[l])>1e-3*qscale[l] else None for l in range(2)]
        rows.append(dict(m=m,k=k,mode=mode['name'],H=mode['H'],alpha=mode['alpha'],xp=xp,delta=delta,
            loads=[P.pair(v) for v in load],reference=[P.pair(v) for v in ref],
            complex_error=(np.abs(load-ref)/den).tolist(),amplitude_error=(np.abs(np.abs(load)-np.abs(ref))/den).tolist(),
            phase_degrees=phase,nonpenetration_relative=float(np.linalg.norm(matrix@g-rhs[:,j])/max(np.linalg.norm(rhs[:,j]),1e-15)),
            conservative_moment_error=float(abs(load[1]-conservative)/max(ms,abs(load[1]))),
            shifted_moment_error=float(abs(shifted-load[1]-.1*load[0])/max(ms,abs(load[1]))),
            wrong_reference_point_difference=float(abs(.1*load[0])/den[1]),
            pressure_profile=profile,trailing_pressure=near))
    return rows,gamma,rhs


def time_response():
    m=128;k=.5;steps=1024;delta=np.pi/(k*steps);modes=MODES[:2]
    model(m);ops=D.discrete(m,k,delta,40.,keep_nodes=True)
    frequency,gamma,rhs_hat=solve_cases(m,k,ops['p1'],delta,modes)
    _,(a,fc,fu,_),(mc0,mu0)=model(m);xp=.25;mc=mc0+xp*fc;mu=mu0+xp*fu
    matrix=a.copy();matrix[:,-1]+=ops['left0'];factor=lu_factor(matrix)
    z=np.exp(2j*k*delta);history=(z**(-np.arange(ops['n']+1))[:,None]*gamma[-1][None,:]).real
    old=gamma.real.copy();dt=delta/8;g=ops['nodes']
    predicted=np.array([[complex(*v) for v in r['loads']] for r in frequency]).T
    reference=np.column_stack([reference_load(k,mode) for mode in modes])
    den=np.column_stack([np.maximum(np.abs(reference_load(k,mode)),.01*np.array(scales(mode)[1:])) for mode in modes])
    numerical=np.zeros((2,2));physical=np.zeros((2,2));residual=0.;trace=[]
    for n in range(1,2*steps+1):
        if n%128==0:budget()
        trial=np.empty_like(history);trial[1:]=history[:-1]
        target_rhs=np.real(rhs_hat*z**n)
        rhs=target_rhs-g[:,1:]@trial[1:]
        current=lu_solve(factor,rhs);trial[0]=current[-1]
        dot=(current-old)/dt
        load=np.array([fc@current+fu@dot,mc@current+mu@dot])
        numerical=np.maximum(numerical,np.abs(load-np.real(predicted*z**n))/den)
        physical=np.maximum(physical,np.abs(load-np.real(reference*z**n))/den)
        residual=max(residual,float(np.linalg.norm(a@current+g@trial-target_rhs)/max(np.linalg.norm(rhs),1e-12)))
        if n%128==0:trace.append(dict(step=n,loads=load.tolist(),frequency_target=np.real(predicted*z**n).tolist()))
        history=trial;old=current
    return dict(m=m,k=k,xp=xp,steps_per_period=steps,periods=2,delta=delta,
        requested_length=40.,effective_length=ops['effective_length'],initialization='periodic nodal history',
        output_rows=['Fz','My'],input_columns=['heave','pitch'],max_discrete_difference=numerical.tolist(),
        max_physical_difference=physical.tolist(),max_nonpenetration_relative=residual,trace=trace,frequency_cases=frequency)


def main():
    tests=reference_tests()
    for m in [128,256,512]:
        for k in [.1,.5,1.]:
            model(m);kw,_=D.continuous(m,k,160.)
            rows,_,_=solve_cases(m,k,kw);CONTINUOUS.extend(rows)
            print(json.dumps(dict(stage='continuous',m=m,k=k,max_pressure_error=max(p['normalized_error'] for r in rows for p in r['pressure_profile']))),flush=True)
    for k in [.1,.5,1.]:
        model(256);ops=D.discrete(256,k,.003125,160.)
        rows,_,_=solve_cases(256,k,ops['p1'],.003125);WAVE.extend(rows)
    pressure_changes=[]
    for k in [.1,.5,1.]:
        for mode in MODES:
            coarse=next(r for r in CONTINUOUS if r['m']==256 and r['k']==k and r['mode']==mode['name'])
            fine=next(r for r in CONTINUOUS if r['m']==512 and r['k']==k and r['mode']==mode['name'])
            ps=scales(mode)[0]
            values=[abs(complex(*c['pressure'])-complex(*f['pressure']))/max(ps,abs(complex(*f['reference']))) for c,f in zip(coarse['pressure_profile'],fine['pressure_profile'])]
            pressure_changes.append(dict(k=k,mode=mode['name'],maximum=float(max(values))))
    model(256);static_ops=D.discrete(256,0.,.003125,160.)
    static,gamma,rhs=solve_cases(256,0.,static_ops['p1'],.003125,[MODES[1]])
    _,(a,fc,fu,_),_=model(256);matrix=a.copy();matrix[:,-1]+=static_ops['p1'].real
    wrong=np.linalg.solve(matrix,np.zeros(256));wrong_force=fc@wrong
    static_row=static[0];ref_static=reference_load(0.,MODES[1])[0]
    static_result=dict(case=static_row,missing_normal_force=P.pair(wrong_force),
        missing_normal_relative=float(abs(wrong_force-ref_static)/abs(ref_static)),
        normalized_quarter_moment=float(abs(complex(*static_row['loads'][1]))/scales(MODES[1])[2]))
    march=time_response();fine=[r for r in CONTINUOUS if r['m']==512]
    checks=dict(C_functions=all(r['C_hankel_difference']<=1e-12 for r in tests),
        reference_integrals=all(r['quadrature_error']<=1e-10 and r['expanded_error']<=1e-10 for r in tests),
        reference_TE_asymptotic=all(r['limit_relative']<=.001 for r in tests),
        local_pressure=all(p['normalized_error']<=.05 for r in fine for p in r['pressure_profile']),
        local_space_change=all(r['maximum']<=.02 for r in pressure_changes),
        W_P1_complex=all(max(r['complex_error'])<=.02 for r in WAVE),
        W_P1_amplitude=all(max(r['amplitude_error'])<=.02 for r in WAVE),
        W_P1_phase=all(all(v is None or abs(v)<=1 for v in r['phase_degrees']) for r in WAVE),
        conservative_moment=all(r['conservative_moment_error']<=1e-11 for r in CONTINUOUS+WAVE),
        reference_point_shift=all(r['shifted_moment_error']<=1e-11 for r in CONTINUOUS+WAVE),
        static_lift=static_row['complex_error'][0]<=.02,static_quarter_moment=static_result['normalized_quarter_moment']<=.001,
        negative_normal=static_result['missing_normal_relative']>.1,
        negative_reference_point=max(r['wrong_reference_point_difference'] for r in WAVE)>.001,
        nonpenetration=all(r['nonpenetration_relative']<=1e-11 for r in CONTINUOUS+WAVE),
        time_discrete_match=bool(np.max(march['max_discrete_difference'])<=1e-10),
        time_nonpenetration=march['max_nonpenetration_relative']<=1e-11,
        time_physical_accuracy=bool(np.max(march['max_physical_difference'])<=.02))
    source=Path(__file__);old=json.loads((ROOT/'rd/aero/wake-loads/records/uvlm-graded-wake-pressure-result.json').read_text())
    deps=['check_uvlm_graded_wake_pressure.py',*old['dependency_sha256']]
    result=dict(protocol='rd/aero/wake-loads/records/uvlm-local-pressure-pitch-experiment.md',geometry_id=P.geometry.GEOMETRY,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={name:hashlib.sha256(rd_paths.sibling(name).read_bytes()).hexdigest() for name in deps},
        U=8.,rho=1.2,c=1.,span=10000.,length=160.,modes=MODES,numpy=np.__version__,scipy=D.scipy.__version__,
        reference_qualification=tests,continuous_cases=CONTINUOUS,W_P1_cases=WAVE,
        pressure_refinement=pressure_changes,static_case=static_result,time_response=march,checks=checks,
        all_criteria_pass=bool(all(checks.values())),seconds=time.monotonic()-START,peak_process_working_bytes=P.memory.memory_bytes(),
        limits='linear 2D reference versus wide finite sheet; local means not P0 point pressure; no free wake or ANCF coupling')
    (ROOT/'rd/aero/wake-loads/records/uvlm-local-pressure-pitch-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=P.core.json_scalar)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,time_errors=march['max_physical_difference'],seconds=result['seconds']),default=P.core.json_scalar),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':
    try:main()
    except Exception as error:
        (ROOT/'build/analysis/uvlm-local-pressure-pitch-stopped.json').write_text(json.dumps(dict(
            error=repr(error),continuous_cases=CONTINUOUS,W_P1_cases=WAVE,seconds=time.monotonic()-START),indent=2,allow_nan=False,default=P.core.json_scalar)+'\n',encoding='utf-8')
        raise
