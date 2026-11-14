"""Limited comparison. Protocol: rd/structure/main-shell/q9-dynamics/records/q9-integrator-candidate-comparison.md."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from numpy.linalg import eigvals
import check_q9_transient_accuracy as diagnostic
import check_q9_plate_convergence as plate
import check_q9_material_tangent as material

ROOT=rd_paths.ROOT
START=time.monotonic()


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parameters(rho):
    tm=(2-rho)/(1+rho);tf=1/(1+rho)
    return tm,tf,.25*(1+tm-tf)**2,.5+tm-tf


def amplification(rho,eta):
    tm,tf,beta,gamma=parameters(rho)
    ex=np.array([1.,0.,0.]);ev=np.array([0.,1.,0.]);ea=np.array([0.,0.,1.])
    pred=ex+ev+(.5-beta)*ea
    acc=-((1-tm)*ea+eta**2*((1-tf)*ex+tf*pred))/(tm+tf*beta*eta**2)
    return np.array([pred+beta*acc,ev+(1-gamma)*ea+gamma*acc,acc])


def scalar_step(state,rho,eta):
    tm,tf,beta,gamma=parameters(rho);x,v,a=state
    xp=x+v+(.5-beta)*a
    # Scalar stage equilibrium solved directly; state uses h=1.
    anew=-((1-tm)*a+eta**2*((1-tf)*x+tf*xp))/(tm+tf*beta*eta**2)
    return np.array([xp+beta*anew,v+(1-gamma)*a+gamma*anew,anew])


def linear_checks():
    records=[];rng=np.random.default_rng(20261030)
    for rho in [1.,.8,.5,0.]:
        state_error=0.;spectral_max=0.
        for eta in np.logspace(-3,4,71):
            a=amplification(rho,eta);state=rng.standard_normal(3)
            state_error=max(state_error,float(np.max(np.abs(a@state-scalar_step(state,rho,eta)))))
            spectral_max=max(spectral_max,float(np.max(np.abs(eigvals(a)))))
        high=float(np.max(np.abs(eigvals(amplification(rho,1e7)))))
        errors=[]
        for n in [40,80,160,320]:
            h=2*np.pi/n;a=amplification(rho,h)
            state=np.array([1.,0.,-h*h])
            for step in range(n):state=a@state
            errors.append(float(np.hypot(state[0]-1,state[1]/h)))
        order=float(np.log2(errors[-2]/errors[-1]))
        records.append(dict(rho=rho,coefficients=parameters(rho),matrix_scalar_error=state_error,
            spectral_max=spectral_max,high_frequency_radius=high,one_period_phase_errors=errors,order=order,
            eligible=bool(state_error<=1e-12 and spectral_max<=1+1e-9 and abs(high-rho)<=1e-4 and 1.7<=order<=2.3)))
    return records


def step(model,y,v,alg,h,rho):
    if time.monotonic()-START>240 or plate.memory_bytes()>4*1024**3:raise TimeoutError('experiment budget')
    tm,tf,beta,gamma=parameters(rho)
    pred=y+h*v+h*h*(.5-beta)*alg
    d=beta*h*h*alg
    for iteration in range(13):
        stage=(1-tf)*y+tf*(pred+d)
        r,k,det=model.evaluate(stage)
        if not np.all(np.isfinite(r)) or det<=0:raise ValueError('invalid_trial')
        residual=tm*d+beta*h*h*((1-tm)*alg+r)
        delta=np.linalg.solve(tm*np.eye(69)+tf*beta*h*h*k,-residual)
        relative=float(np.linalg.norm(residual)/model.amplitude)
        correction=float(np.max(np.abs(model.physical(delta))))
        if relative<=1e-11 and correction<=1e-12:break
        if iteration==12:raise ValueError(dict(reason='newton',rho=rho,relative=relative))
        d+=delta
    anew=d/(beta*h*h)
    ynew=pred+d;vnew=v+h*((1-gamma)*alg+gamma*anew)
    fend=model.fstar+np.einsum('paJk,k->paJ',model.d,ynew)
    enddet=float(np.linalg.det(fend).min())
    if enddet<=0 or not all(np.all(np.isfinite(x)) for x in [ynew,vnew,anew]):raise ValueError('invalid_endpoint')
    return ynew,vnew,anew,dict(residual=relative,iterations=iteration,det_min=min(det,enddet))


def path(model,duration,n,rho,start_state=None,start_number=0):
    h=duration/n
    if start_state is None:
        y=model.y0.copy();v=np.zeros(69);alg=-model.evaluate(y,False)
    else:y,v,alg=[np.array(a) for a in start_state]
    states=[(y.copy(),v.copy(),alg.copy())];ys=[y.copy()];vs=[v.copy()]
    maxres=0.;mindet=1e30;maxiteration=0
    for number in range(start_number,n):
        y,v,alg,metadata=step(model,y,v,alg,h,rho)
        states.append((y.copy(),v.copy(),alg.copy()))
        if (number+1)%(n//160)==0:ys.append(y.copy());vs.append(v.copy())
        maxres=max(maxres,metadata['residual']);mindet=min(mindet,metadata['det_min']);maxiteration=max(maxiteration,metadata['iterations'])
    return dict(y=np.array(ys),v=np.array(vs),states=states,metadata=dict(rho=rho,n=n,residual_max=maxres,
        sampled_det_f_min=mindet,max_iterations=maxiteration))


def stress(model,ys):
    result=[]
    for y in ys:
        h=np.einsum('paJk,k->paJ',model.d,y)
        de=.5*(np.einsum('paI,paJ->pIJ',model.fstar,h)+np.einsum('paI,paJ->pIJ',h,model.fstar)+np.einsum('paI,paJ->pIJ',h,h))
        result.append(material.stress_tensor(material.engineering(de)@model.c.T))
    return np.array(result)


def observables(model,case,reference,ref_acc,ref_stress):
    error=diagnostic.errors(case,reference,model)
    acc=np.array([-model.evaluate(y,False) for y in case['y']])
    error['physical_acceleration_relative']=float(np.linalg.norm(acc-ref_acc,axis=1).max()/np.linalg.norm(ref_acc,axis=1).max())
    ds=stress(model,case['y'])-ref_stress
    norm=np.sqrt(np.einsum('p,tpij,tpij->t',model.w,ds,ds))
    scale=np.sqrt(np.einsum('p,tpij,tpij->t',model.w,ref_stress,ref_stress)).max()
    error['pk2_increment_relative']=float(norm.max()/scale)
    error['low4_phase_space_max']=float(np.sqrt(np.sum(((case['y']-reference['y'])[:,:4]/model.amplitude)**2+
        ((case['v']-reference['v'])[:,:4]/model.initial_speed_scale)**2,axis=1)).max())
    grades={key:bool(error[key]<=1e-3) for key in ['displacement_mass_max','velocity_scaled_mass_max',
        'physical_acceleration_relative','pk2_increment_relative','low4_phase_space_max']}
    return dict(errors=error,observable_eligible=grades)


def restart_check(model,duration,normal):
    n=640;index=237
    encoded=json.dumps(dict(n=index,t=index*duration/n,y=normal['states'][index][0].tolist(),
        v=normal['states'][index][1].tolist(),algorithmic_acceleration=normal['states'][index][2].tolist()),allow_nan=False)
    restored=json.loads(encoded)
    state=[restored[k] for k in ['y','v','algorithmic_acceleration']]
    resumed=path(model,duration,n,.5,state,index)
    equal=all(np.array_equal(a,b) for sa,sb in zip(normal['states'][index:],resumed['states']) for a,b in zip(sa,sb))
    physical=-model.evaluate(np.array(state[0]),False)
    reset=[state[0],state[1],physical.tolist()]
    wrong=path(model,duration,n,.5,reset,index)
    end_ref=dict(y=np.array([normal['states'][-1][0]]),v=np.array([normal['states'][-1][1]]))
    end_bad=dict(y=np.array([wrong['states'][-1][0]]),v=np.array([wrong['states'][-1][1]]))
    difference=diagnostic.errors(end_bad,end_ref,model)
    later_alg=float(np.linalg.norm(wrong['states'][1][2]-resumed['states'][1][2])/(model.w1**2*model.amplitude))
    return dict(checkpoint_step=index,all_endpoint_states_bit_equal=equal,reset_algorithmic_to_physical_difference=difference,
        next_algorithmic_acceleration_scaled_difference=later_alg,
        negative_detected=bool(difference['phase_space_max']>1e-8 or later_alg>1e-8))


def main():
    checks=linear_checks()
    print(json.dumps(dict(stage='linear_checks',cases=checks)),flush=True)
    if not all(c['eligible'] for c in checks):raise ValueError('linear verification failed')
    previous=json.loads((ROOT/'rd/structure/main-shell/q9-dynamics/records/q9-transient-accuracy-result.json').read_text(encoding='utf-8'))
    reference_path=ROOT/previous['trajectory_file']
    if sha(reference_path)!=previous['trajectory_sha256']:raise ValueError('reference hash')
    archive=np.load(reference_path)
    model=diagnostic.Model()
    if not np.array_equal(archive['omega'],model.omega) or not np.array_equal(archive['y0'],model.y0):raise ValueError('reference basis environment')
    reference=dict(y=archive['reference_y'],v=archive['reference_v'])
    duration=previous['short_window_seconds']
    ref_acc=np.array([-model.evaluate(y,False) for y in reference['y']]);ref_stress=stress(model,reference['y'])
    cases=[];normal=None;midpoint_match=None
    for rho in [1.,.8,.5,0.]:
        for n in ([640] if rho==1 else [320,640,1280]):
            case=path(model,duration,n,rho)
            row=dict(**case['metadata'],**observables(model,case,reference,ref_acc,ref_stress))
            cases.append(row)
            if rho==1:
                old=dict(y=archive['short_y_640'],v=archive['short_v_640'])
                midpoint_match=diagnostic.errors(case,old,model)
            if rho==.5 and n==640:normal=case
            print(json.dumps(dict(stage='nonlinear',**row)),flush=True)
    restart=restart_check(model,duration,normal)
    verification=dict(linear_cases=all(c['eligible'] for c in checks),rho1_midpoint=midpoint_match['phase_space_max']<=1e-7,
        algorithmic_history_restart=restart['all_endpoint_states_bit_equal'],incorrect_history_detected=restart['negative_detected'])
    source=Path(__file__)
    result=dict(protocol='rd/structure/main-shell/q9-dynamics/records/q9-integrator-candidate-comparison.md',source_sha256=sha(source),
        dependency_sha256={name:sha(rd_paths.sibling(name)) for name in ['check_q9_transient_accuracy.py','check_q9_material_tangent.py']},
        reference_result_sha256=sha(ROOT/'rd/structure/main-shell/q9-dynamics/records/q9-transient-accuracy-result.json'),reference_trajectory_sha256=sha(reference_path),
        numpy=np.__version__,short_window_seconds=duration,linear_cases=checks,nonlinear_cases=cases,
        rho1_midpoint_difference=midpoint_match,restart=restart,verification_checks=verification,
        all_verification_checks_pass=all(verification.values()),seconds=time.monotonic()-START,
        peak_process_working_bytes=plate.memory_bytes(),
        limitations='same isotropic fixed-support short IVP, all 69 coordinates; candidate observable failures are retained; no full-wing, aero, damping identification or product method adoption')
    (ROOT/'rd/structure/main-shell/q9-dynamics/records/q9-integrator-candidate-result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(verification=verification,seconds=result['seconds'],restart=restart)),flush=True)
    if not result['all_verification_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()
