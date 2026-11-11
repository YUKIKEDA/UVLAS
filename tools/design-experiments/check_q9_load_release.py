"""Limited INIT-RELEASE trial. Protocol: docs/verification/q9-load-release-experiment.md."""
import copy
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
import check_q9_transient_accuracy as diagnostic
import check_q9_integrator_candidates as integrators
import check_q9_material_tangent as material
import check_q9_plate_convergence as plate

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()
PROGRESS=ROOT/'build/analysis/q9-load-release-progress.json'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def progress(data):
    PROGRESS.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def initialize(model):
    original=model.y0.copy();g=model.w1**2*original[0]
    # Direct volume traction-free body load, independent of M times a vector.
    u=model.f.s@model.basis[:,0].reshape(27,3)
    density=diagnostic.previous.full.RHO
    body=density*g*u
    hold=(model.f.s.T@(model.w[:,None]*body)).reshape(81)
    modal=model.basis.T@hold;expected=np.zeros(69);expected[0]=g
    load_error=float(np.linalg.norm(modal-expected)/abs(g))
    y=original.copy();history=[]
    for iteration in range(21):
        r,k,det=model.evaluate(y);residual=r-modal
        delta=np.linalg.solve(k,-residual)
        rel=float(np.linalg.norm(residual)/abs(g))
        correction=float(np.max(np.abs(model.physical(delta))))
        history.append(dict(iteration=iteration,relative_residual=rel,physical_correction_m=correction,sampled_det_f_min=det))
        if rel<=1e-10 and correction<=1e-12:break
        if iteration==20:raise ValueError('holding equilibrium did not converge')
        y+=delta
    q=model.f.qstar+model.basis@y
    cert=diagnostic.previous.full.determinant_certificate(q)
    if not cert['exact_positive']:raise ValueError('held global geometry not certified')
    force=material.analytic_evaluation(q,model.f.ctx)['force']
    full_residual=force-model.f.load-hold
    reaction=np.zeros(81);reaction[diagnostic.previous.full.FIXED]=full_residual[diagnostic.previous.full.FIXED]
    rigid=diagnostic.previous.full.curved.rigid_modes(q,np.zeros(3))
    wrench=rigid.T@(reaction+model.f.load+hold)
    scale=max(np.linalg.norm(rigid.T@model.f.load),np.linalg.norm(rigid.T@hold))
    wrench_error=float(np.linalg.norm(wrench)/scale)
    before=-(model.evaluate(y,False)-modal);after=-model.evaluate(y,False)
    before_error=float(np.linalg.norm(before)/abs(g))
    after_error=float(np.linalg.norm(after+modal)/abs(g))
    initial=dict(kind='INIT-RELEASE',release_time=0.,release_completed=True,
        model_source_sha256=sha(Path(diagnostic.__file__)),persistent_load_id='right-dead-prestress-v1',
        holding_load_id='reference-body-lowest-mode-v1',active_loads=['right-dead-prestress-v1'],
        holding_force_sha256=hashlib.sha256(hold.tobytes()).hexdigest(),
        y=y.tolist(),v=np.zeros(69).tolist(),physical_acceleration=after.tolist(),
        algorithmic_acceleration=after.tolist(),initial_strain=model.e0.tolist())
    restored=json.loads(json.dumps(initial,allow_nan=False))
    roundtrip=(restored==initial and restored['release_completed'] and restored['active_loads']==['right-dead-prestress-v1'])
    ctx=copy.deepcopy(model.f.ctx);ctx['e0'][:]=0.
    wrong_force=material.analytic_evaluation(q,ctx)['force']-model.f.load-hold
    strain_error=float(np.linalg.norm(model.basis.T@wrong_force)/abs(g))
    negatives=[dict(name='holding_force_still_active_after_release',mismatch=float(np.linalg.norm(after-before)/abs(g))),
        dict(name='zero_acceleration_after_release',mismatch=float(np.linalg.norm(after)/abs(g))),
        dict(name='initial_strain_reset_in_holding_equilibrium',mismatch=strain_error)]
    for negative in negatives:negative['detected']=bool(negative['mismatch']>1e-10)
    # Change only the initial condition for this separate physical preparation.
    model.y0=y.copy();model.amplitude=float(np.linalg.norm(y));model.initial_speed_scale=model.w1*model.amplitude
    energy=diagnostic.previous.prestress.potential_increment(model.f.qstar,q,model.f.ctx,model.f.load)
    data=dict(load_coefficient=g,modal_load_relative_error=load_error,newton=history,geometry_certificate=cert,
        holding_body_force_resultant=(np.sum(model.w[:,None]*body,axis=0)).tolist(),
        holding_load_current_wrench=(rigid.T@hold).tolist(),reaction_current_wrench=(rigid.T@reaction).tolist(),
        total_wrench_relative=wrench_error,held_amplitude_m=float(np.max(np.abs(model.physical(y)))),
        previous_ivp_difference_mass_relative=float(np.linalg.norm(y-original)/np.linalg.norm(original)),
        acceleration_before_relative=before_error,acceleration_after_equation_relative=after_error,
        high65_initial_acceleration_relative=float(np.linalg.norm(after[4:])/abs(g)),
        initial_persistent_potential_increment_J=float(energy),potential_definition_switch_J=float(modal@y),
        event_state_roundtrip=roundtrip,negative_controls=negatives,initial_state=initial,
        eligible=bool(load_error<=1e-10 and rel<=1e-10 and correction<=1e-12 and wrench_error<=1e-8 and
            before_error<=1e-10 and after_error<=1e-10 and roundtrip and all(n['detected'] for n in negatives)))
    return data


def observe(model,case,reference):
    acc=np.array([-model.evaluate(y,False) for y in reference['y']])
    stress=integrators.stress(model,reference['y'])
    return integrators.observables(model,case,reference,acc,stress)


def midpoint(model,duration,n):
    h=duration/n;y=model.y0.copy();v=np.zeros(69);ys=[y.copy()];vs=[v.copy()]
    maxres=0.;mindet=1e30;maxit=0
    for number in range(n):
        if time.monotonic()-START>240 or plate.memory_bytes()>4*1024**3:raise TimeoutError('experiment budget')
        z=v.copy()
        for it in range(13):
            r,k,det=model.evaluate(y+.5*h*z)
            residual=z-v+.5*h*r
            delta=np.linalg.solve(np.eye(69)+.25*h*h*k,-residual)
            rel=float(np.linalg.norm(residual)/model.initial_speed_scale)
            corr=float(np.max(np.abs(model.physical(h*delta))))
            if rel<=1e-11 and corr<=1e-12:break
            if it==12:raise ValueError('midpoint Newton limit')
            z+=delta
        y+=h*z;v=2*z-v
        fend=model.fstar+np.einsum('paJk,k->paJ',model.d,y)
        enddet=float(np.linalg.det(fend).min())
        if min(det,enddet)<=0 or not all(np.all(np.isfinite(a)) for a in [y,v]):raise ValueError('invalid geometry or state')
        maxres=max(maxres,rel);mindet=min(mindet,det,enddet);maxit=max(maxit,it)
        if (number+1)%(n//160)==0:ys.append(y.copy());vs.append(v.copy())
    return dict(y=np.array(ys),v=np.array(vs),metadata=dict(n=n,residual_max=maxres,sampled_det_f_min=mindet,max_iterations=maxit))


def dynamic_alpha(model,duration,n):
    tm,tf,beta,gamma=integrators.parameters(.8)
    h=duration/n;y=model.y0.copy();v=np.zeros(69);alg=-model.evaluate(y,False)
    ys=[y.copy()];vs=[v.copy()];maxforce=0.;maxvelocity=0.;mindet=1e30;maxit=0
    for number in range(n):
        if time.monotonic()-START>240:raise TimeoutError('experiment budget')
        pred=y+h*v+h*h*(.5-beta)*alg;d=beta*h*h*alg
        for it in range(13):
            r,k,det=model.evaluate((1-tf)*y+tf*(pred+d))
            residual=tm*d+beta*h*h*((1-tm)*alg+r)
            delta=np.linalg.solve(tm*np.eye(69)+tf*beta*h*h*k,-residual)
            force=float(np.linalg.norm(residual)/(beta*h*h*model.w1**2*model.amplitude))
            speed=float(np.linalg.norm(gamma*delta/(beta*h))/model.initial_speed_scale)
            correction=float(np.max(np.abs(model.physical(delta))))
            if force<=1e-10 and speed<=1e-10 and correction<=1e-12:break
            if it==12:raise ValueError('dynamic alpha Newton limit')
            d+=delta
        anew=d/(beta*h*h);y=pred+d;v+=h*((1-gamma)*alg+gamma*anew);alg=anew
        fend=model.fstar+np.einsum('paJk,k->paJ',model.d,y)
        enddet=float(np.linalg.det(fend).min())
        if min(det,enddet)<=0 or not all(np.all(np.isfinite(a)) for a in [y,v,alg]):raise ValueError('invalid state')
        maxforce=max(maxforce,force);maxvelocity=max(maxvelocity,speed);mindet=min(mindet,det,enddet);maxit=max(maxit,it)
        if (number+1)%(n//160)==0:ys.append(y.copy());vs.append(v.copy())
    return dict(y=np.array(ys),v=np.array(vs),metadata=dict(n=n,rho=.8,residual_scale='dynamic-acceleration',
        force_residual_max=maxforce,velocity_correction_max=maxvelocity,sampled_det_f_min=mindet,max_iterations=maxit))


def main():
    dynamic='--dynamic-residual' in sys.argv
    model=diagnostic.Model();init=initialize(model)
    progress(dict(stage='initialized',initialization=init))
    print(json.dumps(dict(stage='initialized',eligible=init['eligible'],newton=len(init['newton'])-1,
        amplitude=init['held_amplitude_m'],initial_high_acceleration=init['high65_initial_acceleration_relative'])),flush=True)
    if not init['eligible']:raise ValueError('initialization checks failed')
    duration=2*np.pi/model.w1/20
    reused=None;previous_archive=None
    if dynamic:
        reused=json.loads((ROOT/'docs/verification/q9-load-release-initial-result.json').read_text(encoding='utf-8'))
        previous_archive=np.load(ROOT/reused['trajectory_file'])
        if sha(ROOT/reused['trajectory_file'])!=reused['trajectory_sha256']:raise ValueError('reference checksum')
        if not all(np.array_equal(previous_archive[key],value) for key,value in [('omega',model.omega),('basis',model.basis),('y0',model.y0)]):raise ValueError('reused model/basis mismatch')
        reference=dict(metadata=reused['references'][0])
    else:reference=diagnostic.oracle(model,duration,1e-11,1e-14)
    print(json.dumps(dict(stage='reference',**reference['metadata'])),flush=True)
    progress(dict(stage='first_reference',initialization=init,reference=reference['metadata']))
    if dynamic:
        strict=dict(y=previous_archive['reference_y'],v=previous_archive['reference_v'],metadata=reused['references'][1])
        sensitivity=reused['reference_sensitivity']
    else:
        strict=diagnostic.oracle(model,duration,1e-13,1e-16)
        sensitivity=observe(model,reference,strict)
    print(json.dumps(dict(stage='reference_sensitivity',**sensitivity)),flush=True)
    cases=[];trajectories={}
    for method in ['midpoint','generalized-alpha-0.8']:
        for n in [160,320,640,1280]:
            if dynamic and method=='midpoint':
                case=dict(y=previous_archive[f'{method}_{n}_y'],v=previous_archive[f'{method}_{n}_v'],
                    metadata={k:v for k,v in next(r for r in reused['cases'] if r['method']==method and r['n']==n).items()
                        if k not in ['method','errors','observable_eligible','high65_linear_energy_fraction_max']})
            elif method=='midpoint':case=midpoint(model,duration,n)
            else:case=dynamic_alpha(model,duration,n) if dynamic else integrators.path(model,duration,n,.8)
            metrics=observe(model,case,strict)
            energies=.5*np.sum(case['v'][:,4:]**2+(case['y'][:,4:]*model.omega[4:])**2,axis=1)
            initial=.5*np.sum(model.y0**2*model.omega**2)
            row=dict(method=method,**case['metadata'],**metrics,high65_linear_energy_fraction_max=float(energies.max()/initial))
            cases.append(row)
            trajectories[f'{method}_{n}_y']=case['y'];trajectories[f'{method}_{n}_v']=case['v']
            print(json.dumps(dict(stage='time_case',**row)),flush=True)
            progress(dict(stage='time_cases',initialization=init,reference_sensitivity=sensitivity,cases=cases))
    orders=[]
    for method in ['midpoint','generalized-alpha-0.8']:
        rows=[r for r in cases if r['method']==method]
        orders.append(dict(method=method,last_orders={key:float(np.log2(rows[-2]['errors'][key]/rows[-1]['errors'][key]))
            for key in ['displacement_mass_max','velocity_scaled_mass_max','physical_acceleration_relative','pk2_increment_relative']}))
    checks=dict(initialization=init['eligible'],reference_displacement_velocity=sensitivity['errors']['phase_space_max']<=1e-7,
        reference_acceleration=sensitivity['errors']['physical_acceleration_relative']<=1e-5,
        reference_stress=sensitivity['errors']['pk2_increment_relative']<=1e-5)
    for method in ['midpoint','generalized-alpha-0.8']:
        fine=[r for r in cases if r['method']==method][-1]
        checks[method+'_finest_observables']=all(fine['observable_eligible'].values())
    archive=ROOT/('build/analysis/q9-load-release-dynamic-trajectories.npz' if dynamic else 'build/analysis/q9-load-release-trajectories.npz')
    np.savez_compressed(archive,omega=model.omega,basis=model.basis,y0=model.y0,reference_y=strict['y'],reference_v=strict['v'],**trajectories)
    source=Path(__file__)
    result=dict(protocol='docs/verification/q9-load-release-experiment.md',dynamic_residual=dynamic,source_sha256=sha(source),
        reused_initial_result_sha256=sha(ROOT/'docs/verification/q9-load-release-initial-result.json') if dynamic else None,
        dependency_sha256={name:sha(source.with_name(name)) for name in ['check_q9_transient_accuracy.py','check_q9_integrator_candidates.py','check_q9_material_tangent.py']},
        initialization=init,duration_seconds=duration,references=[reference['metadata'],strict['metadata']],
        reference_sensitivity=sensitivity,cases=cases,orders=orders,checks=checks,all_criteria_pass=all(checks.values()),
        trajectory_file=str(archive.relative_to(ROOT)),trajectory_sha256=sha(archive),seconds=time.monotonic()-START,
        peak_process_working_bytes=plate.memory_bytes(),limitations='different prepared IVP, fixed 69-DOF one-element isotropic model, sampled transient geometry; no general release/wake event or product integrator adoption')
    (ROOT/('docs/verification/q9-load-release-dynamic-result.json' if dynamic else 'docs/verification/q9-load-release-result.json')).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    progress(dict(stage='completed',checks=checks))
    print(json.dumps(dict(stage='completed',checks=checks,seconds=result['seconds'],orders=orders)),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':
    try:main()
    except Exception as exc:
        failure=dict(status='stopped',reason=repr(exc),source_sha256=sha(Path(__file__)),seconds=time.monotonic()-START)
        (ROOT/'build/analysis/q9-load-release-stopped.json').write_text(json.dumps(failure,indent=2)+'\n',encoding='utf-8')
        raise
