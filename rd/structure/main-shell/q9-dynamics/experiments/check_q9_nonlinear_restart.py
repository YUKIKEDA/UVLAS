"""Limited full-coordinate structural midpoint/restart experiment; not product code.
Prewritten protocol: rd/structure/main-shell/q9-dynamics/records/q9-nonlinear-restart-experiment.md.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import copy
import hashlib
import json
import os
import sys
import time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
import check_q9_full_space_equilibrium as full
import check_q9_material_tangent as material
import check_q9_prestress_equilibrium as prestress
import check_q9_plate_convergence as plate

ROOT=rd_paths.ROOT
START=time.monotonic()


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def seal(payload):
    return dict(payload=payload,sha256=digest(payload))


def clone(state):
    return json.loads(json.dumps(state,allow_nan=False))


class Fixture:
    def __init__(self,manufactured=False):
        self.manufactured=manufactured
        self.ctx,self.mass,self.s=full.make_context(full.Q0,np.eye(3))
        self.p=np.eye(81)[:,full.FREE]*full.W[full.FREE]
        self.qstar=full.exact_q(np.eye(3),np.zeros(3))
        self.load=full.face_load(full.Q0,np.eye(3),1)
        self.m=self.p.T@self.mass@self.p
        initial=material.analytic_evaluation(self.qstar,self.ctx)
        self.scale=float(np.linalg.norm(self.p.T@self.load))
        self.equilibrium_residual=float(np.linalg.norm(self.p.T@(initial['force']-self.load))/self.scale)
        values,vectors=eigh(self.p.T@initial['k']@self.p,self.m)
        self.period=2*np.pi/np.sqrt(values[0])
        phi=vectors[:,0]
        self.x0=phi*(2e-4/np.max(np.abs(self.s@(self.p@phi).reshape(27,3))))
        self.model=dict(q0=full.Q0.tolist(),qstar=self.qstar.tolist(),p=self.p.tolist(),
            mass=self.mass.tolist(),load=self.load.tolist(),e0=self.ctx['e0'].tolist(),
            material=self.ctx['c'].tolist(),axes=self.ctx['axes'].tolist(),
            gradients=self.ctx['g'].tolist(),weights=self.ctx['w'].tolist())
        self.model['forcing']='manufactured-cosine-v1' if manufactured else 'dead-only'
        if manufactured:
            self.model['forcing_parameters']=dict(target=self.x0.tolist(),omega=float(2*np.pi/self.period))
        self.model_id=digest(self.model)

    def evaluate(self,x):
        result=material.analytic_evaluation(self.qstar+self.p@x,self.ctx)
        return self.p.T@(result['force']-self.load),self.p.T@result['k']@self.p,result['sampled_det_f_min']

    def energy(self,state):
        x=np.array(state['x']);v=np.array(state['v'])
        return float(.5*v@self.m@v+prestress.potential_increment(self.qstar,self.qstar+self.p@x,self.ctx,self.load))

    def forcing(self,t):
        if not self.manufactured:return np.zeros(69)
        omega=2*np.pi/self.period
        target=self.x0*np.cos(omega*t)
        return self.evaluate(target)[0]-omega**2*(self.m@target)


def initial_state(f):
    return dict(accepted=True,n=0,t=0.,last_dt=None,x=f.x0.tolist(),v=np.zeros(69).tolist(),output_cursor=0)


def policy(schedule):
    return dict(algorithm='structural-midpoint-v1',schedule=schedule,residual_tolerance=1e-9,
        physical_correction_tolerance=1e-11,max_newton=12,dtype='float64',coordinates=69)


def try_step(f,state,dt,max_newton=12):
    # Never assign through an input array or mutate the accepted dictionary.
    x=np.array(state['x']);v=np.array(state['v']);y=x+dt*v
    forcing=f.forcing(state['t']+dt/2)
    for iteration in range(max_newton+1):
        if time.monotonic()-START>240 or plate.memory_bytes()>4*1024**3:
            raise TimeoutError('prewritten experiment budget')
        r,k,det=f.evaluate((x+y)/2)
        residual=2*f.m@(y-x-dt*v)/dt**2+r-forcing
        delta=np.linalg.solve(2*f.m/dt**2+.5*k,-residual)
        relative=float(np.linalg.norm(residual)/f.scale)
        correction=float(np.max(np.abs(f.s@(f.p@delta).reshape(27,3))))
        if relative<=1e-9 and correction<=1e-11:
            _,_,enddet=f.evaluate(y)
            result=dict(accepted=True,n=state['n']+1,t=state['t']+dt,last_dt=dt,
                x=y.tolist(),v=(2*(y-x)/dt-v).tolist(),output_cursor=state['n']+1)
            return dict(parent=digest(state),candidate=result,residual=relative,
                correction=correction,det_min=min(det,enddet),iterations=iteration)
        if iteration==max_newton:
            return dict(parent=digest(state),rejected='iteration_limit',residual=relative)
        y+=delta
    raise AssertionError('unreachable')


def commit(state,trial):
    if trial['parent']!=digest(state):raise ValueError('stale_parent_or_double_commit')
    if 'rejected' in trial:raise ValueError('rejected_trial')
    return clone(trial['candidate'])


def validate(envelope,f,settings):
    payload=envelope['payload']
    if envelope['sha256']!=digest(payload):raise ValueError('checksum')
    if payload['schema']!='q9-structural-checkpoint-v1':raise ValueError('schema')
    if payload['model_id']!=f.model_id:raise ValueError('model_identity')
    if payload['policy_id']!=digest(settings):raise ValueError('algorithm_identity')
    state=payload['state']
    if state['accepted'] is not True:raise ValueError('trial_state')
    n=state['n']
    if type(n) is not int or n<0 or n>len(settings['schedule']):raise ValueError('step_number')
    if state['output_cursor']!=n:raise ValueError('output_cursor')
    if state['t']!=sum(settings['schedule'][:n]):raise ValueError('time')
    if state['last_dt']!=(settings['schedule'][n-1] if n else None):raise ValueError('last_dt')
    for key in ['x','v']:
        value=np.asarray(state[key],dtype=float)
        if value.shape!=(69,) or not np.all(np.isfinite(value)):raise ValueError('state_array')
    # Rebuild the equations and check current geometry; no factorization is read.
    f.evaluate(np.array(state['x']))
    return clone(state)


def save(path,f,settings,state):
    envelope=seal(dict(schema='q9-structural-checkpoint-v1',model_id=f.model_id,
        policy_id=digest(settings),state=clone(state)))
    validate(envelope,f,settings)
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(envelope,allow_nan=False),encoding='utf-8')
    validate(json.loads(temporary.read_text(encoding='utf-8')),f,settings)
    os.replace(temporary,path)
    return envelope


def run(f,schedule,restart=False):
    settings=policy(schedule);state=initial_state(f);history=[];keys=[0]
    e0=f.energy(state);maxenergy=0.;maxres=0.;mindet=1e30;iterations=0;work=0.
    rollback=True;guards=True
    for step,dt in enumerate(schedule):
        if restart and step==17:
            before=digest(state);failed=try_step(f,state,dt,max_newton=0)
            rollback=bool('rejected' in failed and digest(state)==before)
            try:commit(state,failed);guards=False
            except ValueError:pass
        trial=try_step(f,state,dt)
        if 'rejected' in trial:raise RuntimeError(dict(step=step,trial=trial))
        work+=float(f.forcing(state['t']+dt/2)@(np.array(trial['candidate']['x'])-np.array(state['x'])))
        state=commit(state,trial)
        if restart and step==17:
            try:commit(state,trial);guards=False
            except ValueError:pass
        if restart and state['n'] in [17,53]:
            path=ROOT/'build/analysis/q9-structural-checkpoint.json'
            save(path,f,settings,state)
            # Fresh fixture reconstructs M/K/model context, never reuses a cache.
            fresh=Fixture(f.manufactured)
            state=validate(json.loads(path.read_text(encoding='utf-8')),fresh,settings)
            f=fresh
        maxenergy=max(maxenergy,abs((f.energy(state)-e0-work)/e0))
        maxres=max(maxres,trial['residual']);mindet=min(mindet,trial['det_min'])
        iterations=max(iterations,trial['iterations'])
        history.append(digest(state));keys.append(state['output_cursor'])
    return state,dict(steps=len(schedule),energy_relative_max=maxenergy,residual_max=maxres,
        sampled_det_f_min=mindet,max_newton_iterations=iterations,rollback_unchanged=rollback,
        commit_guards=guards,trajectory_hash=digest(history),output_keys=keys,
        history=history)


def negatives(f,schedule):
    settings=policy(schedule);state=initial_state(f)
    envelope=seal(dict(schema='q9-structural-checkpoint-v1',model_id=f.model_id,
        policy_id=digest(settings),state=state))
    records=[]
    def check(name,modified,model=f,config=settings):
        try:validate(modified,model,config)
        except (ValueError,KeyError,TypeError):records.append(dict(name=name,detected=True));return
        records.append(dict(name=name,detected=False))
    bad=copy.deepcopy(envelope);bad['sha256']='0'*64;check('checksum',bad)
    for name,key,value in [('schema','schema','future'),('model','model_id','other'),
        ('algorithm','policy_id','other')]:
        payload=copy.deepcopy(envelope['payload']);payload[key]=value;check(name,seal(payload))
    for name,key,value in [('trial','accepted',False),('step','n',-1),('time','t',1.),
        ('cursor','output_cursor',1),('last_dt','last_dt',1.),('array_size','x',[0.])]:
        payload=copy.deepcopy(envelope['payload']);payload['state'][key]=value;check(name,seal(payload))
    payload=copy.deepcopy(envelope['payload']);del payload['state']['v'];check('missing_velocity',seal(payload))
    # Nonfinite JSON is also refused before constructing a valid envelope.
    payload=copy.deepcopy(envelope['payload']);payload['state']['v'][0]=float('nan')
    try:seal(payload)
    except ValueError:records.append(dict(name='nonfinite',detected=True))
    else:records.append(dict(name='nonfinite',detected=False))
    altered=Fixture(f.manufactured);altered.model['e0'][0][0][0]+=1e-6;altered.model_id=digest(altered.model)
    check('changed_initial_strain',envelope,model=altered)
    for name,text in [('truncated',json.dumps(envelope)[:100]),('corrupt','{x}')]:
        try:validate(json.loads(text),f,settings)
        except (ValueError,KeyError,TypeError):records.append(dict(name=name,detected=True))
        else:records.append(dict(name=name,detected=False))
    return records


def main():
    manufactured='--manufactured' in sys.argv
    f=Fixture(manufactured);cases=[];states=[]
    refined='--refine' in sys.argv
    for n in ([20,40,80,160,320,640,1280,2560] if refined else [20,40,80,160]):
        state,record=run(f,[f.period/(2*n)]*n)
        states.append(state);cases.append(record)
        print(json.dumps(dict(stage='time_convergence',**{k:v for k,v in record.items() if k not in ['history','output_keys']})),flush=True)
    norm0=np.sqrt(f.x0@f.m@f.x0)
    differences=[]
    for a,b in zip(states[:-1],states[1:]):
        d=np.array(a['x'])-np.array(b['x'])
        differences.append(float(np.sqrt(d@f.m@d)/norm0))
    exact_errors=[]
    if manufactured:
        omega=2*np.pi/f.period
        for state in states:
            dx=np.array(state['x'])-f.x0*np.cos(omega*state['t'])
            dv=np.array(state['v'])/omega+f.x0*np.sin(omega*state['t'])
            exact_errors.append(float(np.sqrt(dx@f.m@dx+dv@f.m@dv)/norm0))
    order=float(np.log2(exact_errors[-2]/exact_errors[-1] if manufactured else differences[-2]/differences[-1]))
    schedule=[f.period/160]*80
    resumed,record=run(f,schedule,True)
    equal=bool(digest(resumed)==digest(states[2]) and record['history']==cases[2]['history'] and record['output_keys']==cases[2]['output_keys'])
    variable=schedule[:40]+[f.period/320]*80
    normal,normal_record=run(f,variable)
    split,split_record=run(f,variable,True)
    variable_equal=bool(digest(normal)==digest(split) and normal_record['history']==split_record['history'] and normal_record['output_keys']==split_record['output_keys'])
    controls=negatives(f,schedule)
    checks=dict(equilibrium=f.equilibrium_residual<=1e-9,time_difference=(exact_errors[-1] if manufactured else differences[-1])<=1e-3,
        order=1.7<=order<=2.3,energy=cases[-1]['energy_relative_max']<=1e-3,
        step_residual=max(c['residual_max'] for c in cases+[record,normal_record,split_record])<=1e-9,
        restart_bit_equal=equal,variable_step_restart_bit_equal=variable_equal,
        rollback=record['rollback_unchanged'] and split_record['rollback_unchanged'],
        commit_guards=record['commit_guards'] and split_record['commit_guards'],
        invalid_checkpoints=all(c['detected'] for c in controls))
    for case in cases+[record,normal_record,split_record]:case.pop('history',None)
    result=dict(protocol='rd/structure/main-shell/q9-dynamics/records/q9-nonlinear-restart-experiment.md',refined=refined,manufactured=manufactured,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        equilibrium_residual=f.equilibrium_residual,period=f.period,initial_physical_amplitude=2e-4,
        time_cases=cases,successive_endpoint_mass_relative=differences,exact_phase_space_errors=exact_errors,last_order=order,
        restart=record,variable_step=normal_record,variable_step_restart=split_record,
        negative_controls=controls,checks=checks,all_criteria_pass=all(checks.values()),
        seconds=time.monotonic()-START,peak_process_working_bytes=plate.memory_bytes(),
        limitations='one element, fixed linear support, isotropic material, nonzero eigenstrain, dead load; sampled geometry only; no aero/wake, OS-crossing or power-loss test')
    path=ROOT/('rd/structure/main-shell/q9-dynamics/records/q9-nonlinear-restart-manufactured-result.json' if manufactured else 'rd/structure/main-shell/q9-dynamics/records/q9-nonlinear-restart-result.json')
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],last_order=order,differences=differences)),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':main()
