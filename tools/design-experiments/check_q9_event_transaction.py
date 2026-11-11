"""Limited event transaction/restart trial; no production time-type adoption.
Protocol: docs/verification/q9-event-transaction-experiment.md.
"""
import copy
import hashlib
import json
import os
import time
from pathlib import Path
import numpy as np
import check_q9_load_release as release
import check_q9_transient_accuracy as diagnostic
import check_q9_integrator_candidates as candidates
import check_q9_plate_convergence as plate

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()
EVENT=37
END=357
HOLD=['persistent-right','holding-body']
BASE=['persistent-right']


class AcceptedState(dict):
    """Internal API marker; deliberately absent from uncommitted candidates."""
    pass


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def clone(value):return json.loads(json.dumps(value,allow_nan=False))


class Fixture:
    def __init__(self):
        self.model=diagnostic.Model();self.init=release.initialize(self.model)
        if not self.init['eligible']:raise ValueError('holding initialization')
        self.unit=(2*np.pi/self.model.w1/20)/320
        u=self.model.f.s@self.model.basis[:,0].reshape(27,3)
        body=diagnostic.previous.full.RHO*self.init['load_coefficient']*u
        full=(self.model.f.s.T@(self.model.w[:,None]*body)).reshape(81)
        self.hold=self.model.basis.T@full
        self.identity=digest(dict(model=self.model.f.model,basis=self.model.basis.tolist(),hold=self.hold.tolist(),
            initial_y=self.model.y0.tolist(),velocity_factor=.02,unit=self.unit,event_tick=EVENT,end_tick=END,
            integrator='GA08-dynamic-v1',force_tol=1e-10,velocity_tol=1e-10,physical_tol=1e-12,
            event_history='physical-acceleration-restart-v1'))

    def acceleration(self,y,cursor):
        return -self.model.evaluate(np.asarray(y),False)+(self.hold if cursor==0 else 0.)

    def initial(self):
        y=self.model.y0.copy();v=.02*self.model.w1*y;a=self.acceleration(y,0)
        return AcceptedState(tick=0,t=0.,n=0,revision=0,output_cursor=0,event_cursor=0,history_epoch=0,
            phase='holding',loads=HOLD.copy(),y=y.tolist(),v=v.tolist(),algorithmic_acceleration=a.tolist(),
            physical_acceleration=a.tolist())

    def try_step(self,state,ticks,max_newton=12):
        parent=digest(state)
        def rejected(reason):return dict(parent=parent,rejected=reason)
        if type(ticks) is not int or ticks<=0:return rejected('invalid_tick_step')
        if state['tick']==EVENT and state['event_cursor']==0:return rejected('pending_event')
        if state['tick']<EVENT<state['tick']+ticks:return rejected('event_crossing')
        if state['tick']+ticks>END:return rejected('past_end')
        h=ticks*self.unit;model=self.model
        tm,tf,beta,gamma=candidates.parameters(.8)
        y=np.array(state['y']);v=np.array(state['v']);alg=np.array(state['algorithmic_acceleration'])
        applied=self.hold if state['event_cursor']==0 else np.zeros(69)
        pred=y+h*v+h*h*(.5-beta)*alg;d=beta*h*h*alg
        try:
            for it in range(max_newton+1):
                if time.monotonic()-START>240 or plate.memory_bytes()>4*1024**3:raise TimeoutError('experiment budget')
                r,k,det=model.evaluate((1-tf)*y+tf*(pred+d));r-=applied
                if det<=0 or not np.all(np.isfinite(r)):return rejected('geometry_or_nonfinite')
                residual=tm*d+beta*h*h*((1-tm)*alg+r)
                delta=np.linalg.solve(tm*np.eye(69)+tf*beta*h*h*k,-residual)
                force=float(np.linalg.norm(residual)/(beta*h*h*model.w1**2*model.amplitude))
                speed=float(np.linalg.norm(gamma*delta/(beta*h))/model.initial_speed_scale)
                physical=float(np.max(np.abs(model.physical(delta))))
                if force<=1e-10 and speed<=1e-10 and physical<=1e-12:break
                if it==max_newton:return rejected('newton_budget')
                d+=delta
        except (ValueError,np.linalg.LinAlgError):return rejected('evaluation_or_linear_solve')
        anew=d/(beta*h*h);yn=pred+d;vn=v+h*((1-gamma)*alg+gamma*anew)
        f=model.fstar+np.einsum('paJk,k->paJ',model.d,yn)
        if np.linalg.det(f).min()<=0 or not all(np.all(np.isfinite(a)) for a in [yn,vn,anew]):return rejected('endpoint_geometry')
        candidate=clone(state);tick=state['tick']+ticks
        candidate.update(tick=tick,t=tick*self.unit,n=state['n']+1,revision=state['revision']+1,
            output_cursor=state['output_cursor']+1,phase='event-left' if tick==EVENT else ('released' if state['event_cursor'] else 'holding'),
            y=yn.tolist(),v=vn.tolist(),algorithmic_acceleration=anew.tolist(),physical_acceleration=self.acceleration(yn,state['event_cursor']).tolist())
        return dict(parent=parent,candidate=candidate,residual=force,iterations=it)

    def try_event(self,state):
        parent=digest(state)
        if state['event_cursor']!=0:return dict(parent=parent,rejected='already_applied')
        if state['tick']!=EVENT or state['phase']!='event-left':return dict(parent=parent,rejected='not_at_left_boundary')
        candidate=clone(state);a=self.acceleration(state['y'],1).tolist()
        candidate.update(event_cursor=1,history_epoch=1,revision=state['revision']+1,output_cursor=state['output_cursor']+1,
            phase='event-right',loads=BASE.copy(),algorithmic_acceleration=a,physical_acceleration=a)
        return dict(parent=parent,candidate=candidate)


def commit(state,trial):
    if trial['parent']!=digest(state):raise ValueError('stale_parent_or_double_commit')
    if 'rejected' in trial:raise ValueError(trial['rejected'])
    return AcceptedState(clone(trial['candidate']))


def validate(envelope,fixture):
    payload=envelope['payload']
    if envelope['checksum']!=digest(payload):raise ValueError('checksum')
    if payload['schema']!='q9-event-checkpoint-v1' or payload['identity']!=fixture.identity:raise ValueError('identity_or_schema')
    state=payload['state'];tick=state['tick'];cursor=state['event_cursor']
    if type(tick) is not int or tick<0 or tick>END:raise ValueError('clock')
    if type(cursor) is not int or cursor not in [0,1] or (tick<EVENT and cursor!=0) or (tick>EVENT and cursor!=1):raise ValueError('event_cursor')
    if tick<EVENT:
        if tick%2:raise ValueError('clock_schedule')
        n=tick//2;phase='holding'
    else:
        if (tick-EVENT)%2:raise ValueError('clock_schedule')
        n=19+(tick-EVENT)//2
        phase=('event-right' if cursor else 'event-left') if tick==EVENT else 'released'
    if state['n']!=n or state['revision']!=n+cursor or state['output_cursor']!=n+cursor:raise ValueError('generation_or_output')
    if state['t']!=tick*fixture.unit or state['history_epoch']!=cursor or state['phase']!=phase:raise ValueError('phase_time_history')
    if state['loads']!=(HOLD if cursor==0 else BASE):raise ValueError('loads')
    for key in ['y','v','algorithmic_acceleration','physical_acceleration']:
        value=np.asarray(state[key])
        if value.shape!=(69,) or not np.all(np.isfinite(value)):raise ValueError('state_array')
    a=fixture.acceleration(state['y'],cursor)
    if np.linalg.norm(a-np.asarray(state['physical_acceleration']))>1e-10*fixture.model.w1**2*fixture.model.amplitude:raise ValueError('physical_acceleration')
    if phase=='event-right' and not np.array_equal(a,np.array(state['algorithmic_acceleration'])):raise ValueError('right_history')
    f=fixture.model.fstar+np.einsum('paJk,k->paJ',fixture.model.d,np.array(state['y']))
    if np.linalg.det(f).min()<=0:raise ValueError('geometry')
    return AcceptedState(clone(state))


def envelope(fixture,state):
    payload=dict(schema='q9-event-checkpoint-v1',identity=fixture.identity,state=clone(state))
    return dict(payload=payload,checksum=digest(payload))


def save_load(fixture,state,label):
    if not isinstance(state,AcceptedState):raise ValueError('uncommitted_candidate_cannot_be_saved')
    path=ROOT/f'build/analysis/q9-event-{label}.json';temporary=path.with_suffix('.tmp')
    data=envelope(fixture,state);validate(data,fixture)
    temporary.write_text(json.dumps(data,allow_nan=False),encoding='utf-8')
    validate(json.loads(temporary.read_text(encoding='utf-8')),fixture);os.replace(temporary,path)
    fresh=Fixture()
    return fresh,validate(json.loads(path.read_text(encoding='utf-8')),fresh)


def run(fixture,start=None,inject=False,breakpoint=None):
    state=fixture.initial() if start is None else AcceptedState(clone(start))
    states=[clone(state)];rollback=True;double=True;resumed=False
    while state['tick']<END:
        if breakpoint and not resumed and (state['tick'],state['phase'])==breakpoint:
            before=digest(state);fixture,state=save_load(fixture,state,f'{state["tick"]}-{state["phase"]}')
            if before!=digest(state):raise ValueError('checkpoint roundtrip')
            resumed=True
        if state['tick']==EVENT and state['event_cursor']==0:
            trial=fixture.try_event(state)
            if inject:
                before=digest(state);discarded=fixture.try_event(state)
                rollback=rollback and digest(state)==before and 'candidate' in discarded
            state=commit(state,trial)
            if inject:
                try:commit(state,trial);double=False
                except ValueError:pass
        else:
            ticks=min(2,END-state['tick'])
            if state['tick']<EVENT: ticks=min(ticks,EVENT-state['tick'])
            if inject and state['tick']==EVENT+4:
                before=digest(state);trial=fixture.try_step(state,ticks,max_newton=0)
                rollback=rollback and 'rejected' in trial and trial['rejected']=='newton_budget' and digest(state)==before
                try:commit(state,trial);rollback=False
                except ValueError:pass
            trial=fixture.try_step(state,ticks)
            state=commit(state,trial)
        states.append(clone(state))
    return states,dict(rollback=rollback,double_commit_refused=double,resumed=resumed)


def segmented_reference(fixture):
    state=fixture.initial();states=[clone(state)]
    for ticks in [2]*18+[1]:
        state=commit(state,fixture.try_step(state,ticks));states.append(clone(state))
    # Manual event boundary, independent of scheduler and try_event routine.
    state=clone(state);a=fixture.acceleration(state['y'],1).tolist()
    state.update(event_cursor=1,history_epoch=1,revision=20,output_cursor=20,phase='event-right',loads=BASE.copy(),
        algorithmic_acceleration=a,physical_acceleration=a)
    states.append(clone(state))
    for ticks in [2]*160:
        state=commit(state,fixture.try_step(state,ticks));states.append(clone(state))
    return states


def refused(call):
    try:call()
    except (ValueError,KeyError,TypeError):return True
    return False


def negative_checkpoints(fixture,left,right):
    records=[]
    for name,key,value in [('loads','loads',HOLD),('cursor','event_cursor',0),('phase','phase','event-left'),
        ('time','t',right['t']+fixture.unit),('step','n',right['n']+1),('revision','revision',999),
        ('output','output_cursor',999),('history_epoch','history_epoch',0),('array_length','v',[0.]),
        ('left_history_on_right','algorithmic_acceleration',left['algorithmic_acceleration'])]:
        state=clone(right);state[key]=value
        records.append(dict(name=name,detected=refused(lambda:validate(envelope(fixture,state),fixture))))
    e=envelope(fixture,right);e['checksum']='bad'
    records.append(dict(name='checksum',detected=refused(lambda:validate(e,fixture))))
    e=envelope(fixture,right);e['payload']['identity']='other';e['checksum']=digest(e['payload'])
    records.append(dict(name='model_or_policy_identity',detected=refused(lambda:validate(e,fixture))))
    state=clone(right);del state['algorithmic_acceleration']
    records.append(dict(name='missing_algorithmic_acceleration',detected=refused(lambda:validate(envelope(fixture,state),fixture))))
    return records


def main():
    fixture=Fixture();normal,_=run(fixture)
    reference=segmented_reference(fixture)
    hashes=[digest(s) for s in normal];keys=[s['output_cursor'] for s in normal]
    injected,guards=run(fixture,inject=True)
    restart_cases=[]
    for point in [(36,'holding'),(37,'event-left'),(37,'event-right'),(57,'released')]:
        resumed,metadata=run(fixture,inject=True,breakpoint=point)
        equal=bool([digest(s) for s in resumed]==hashes and [s['output_cursor'] for s in resumed]==keys)
        restart_cases.append(dict(tick=point[0],phase=point[1],all_states_and_output_equal=equal,**metadata))
        print(json.dumps(dict(stage='restart',**restart_cases[-1])),flush=True)
    left=next(s for s in normal if s['phase']=='event-left');right=next(s for s in normal if s['phase']=='event-right')
    s36=next(s for s in normal if s['tick']==36)
    guards['crossing_refused']='rejected' in fixture.try_step(s36,2)
    guards['pending_event_refused']='rejected' in fixture.try_step(left,2)
    guards['early_event_refused']='rejected' in fixture.try_event(fixture.initial())
    guards['repeat_event_refused']='rejected' in fixture.try_event(right)
    stale=fixture.try_step(fixture.initial(),2)
    guards['stale_parent_refused']=refused(lambda:commit(right,stale))
    guards['time_candidate_save_refused']=refused(lambda:save_load(fixture,stale['candidate'],'uncommitted-step'))
    event_candidate=fixture.try_event(left)['candidate']
    guards['event_candidate_save_refused']=refused(lambda:save_load(fixture,event_candidate,'uncommitted-event'))
    jump=np.array(right['physical_acceleration'])-np.array(left['physical_acceleration'])
    jump_error=float(np.linalg.norm(jump+fixture.hold)/np.linalg.norm(fixture.hold))
    bad=clone(right);bad['algorithmic_acceleration']=left['algorithmic_acceleration']
    wrong,_=run(fixture,start=bad)
    end=normal[-1];wrongend=wrong[-1]
    difference=diagnostic.errors(dict(y=np.array([wrongend['y']]),v=np.array([wrongend['v']])),
        dict(y=np.array([end['y']]),v=np.array([end['v']])),fixture.model)
    history_error=float(np.linalg.norm(np.array(wrong[1]['algorithmic_acceleration'])-np.array(normal[21]['algorithmic_acceleration']))/
        (fixture.model.w1**2*fixture.model.amplitude))
    negatives=negative_checkpoints(fixture,left,right)
    checks=dict(segmented_reference=[digest(s) for s in reference]==hashes,
        rollback_trajectory=[digest(s) for s in injected]==hashes,guards=all(v for k,v in guards.items() if k!='resumed'),
        restart_all=all(c['all_states_and_output_equal'] and c['resumed'] and c['rollback'] for c in restart_cases),
        event_once=end['event_cursor']==1 and end['n']==179 and end['revision']==180,
        output_once=keys==list(range(181)),
        event_position_velocity_continuous=left['y']==right['y'] and left['v']==right['v'],
        acceleration_jump=jump_error<=1e-12,invalid_checkpoints=all(n['detected'] for n in negatives),
        missing_history_reset_detected=difference['phase_space_max']>1e-8 or history_error>1e-8)
    archive=ROOT/'build/analysis/q9-event-transaction-states.npz'
    np.savez_compressed(archive,y=np.array([s['y'] for s in normal]),v=np.array([s['v'] for s in normal]),
        algorithmic_acceleration=np.array([s['algorithmic_acceleration'] for s in normal]),tick=np.array([s['tick'] for s in normal]))
    source=Path(__file__)
    result=dict(protocol='docs/verification/q9-event-transaction-experiment.md',source_sha256=release.sha(source),
        dependency_sha256={name:release.sha(source.with_name(name)) for name in ['check_q9_load_release.py','check_q9_transient_accuracy.py','check_q9_integrator_candidates.py']},
        tick_seconds=fixture.unit,release_seconds=EVENT*fixture.unit,end_seconds=END*fixture.unit,initialization=fixture.init,
        guards=guards,restarts=restart_cases,checks=checks,all_criteria_pass=all(checks.values()),
        event_left=left,event_right=right,acceleration_jump_relative=jump_error,
        invalid_checkpoints=negatives,wrong_history_final_difference=difference,wrong_history_next_algorithmic_scaled_difference=history_error,
        state_hashes=hashes,output_keys=keys,trajectory_file=str(archive.relative_to(ROOT)),trajectory_sha256=release.sha(archive),
        seconds=time.monotonic()-START,peak_process_working_bytes=plate.memory_bytes(),
        limitations='integer ticks, one known finite-force event, fixed supports, moving isotropic one-element 69-DOF trajectory; transition equivalence only, not independent time-accuracy or product checkpoint durability')
    (ROOT/'docs/verification/q9-event-transaction-result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],jump_error=jump_error,history_difference=difference)),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':
    try:main()
    except Exception as exc:
        (ROOT/'build/analysis/q9-event-transaction-stopped.json').write_text(json.dumps(dict(reason=repr(exc),source_sha256=release.sha(Path(__file__)),seconds=time.monotonic()-START),indent=2)+'\n',encoding='utf-8')
        raise
