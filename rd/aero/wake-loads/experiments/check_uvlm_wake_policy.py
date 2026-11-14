import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from dataclasses import replace
import json,hashlib,time,traceback,copy
import numpy as np
import uvlm_wake_assembly_v2 as aero
import uvlm_surface_bindings as binding
ROOT=rd_paths.ROOT
def main():
    start=time.monotonic();rows=[];failures=[]
    try:
        b=binding.q1('wing',1,1,[[0,-.5,0],[0,.5,0],[1,.5,0],[1,-.5,0]]);model=aero.Model([b],.025);flow=aero.Flow([1,0,.03]);s=aero.cold_state(model,b['q0']);store=aero.Store(model,s);ws=aero.Workspace()
        for i in range(2):
            t=store.begin(s.q,s.v,.1,flow,mode='prescribed');store.solve(t,flow,ws);store.accept(t)
        current=store.state;before=current.fingerprint();encoded=aero.checkpoint(model,current);restored=aero.restore(model,encoded,flow);resumed=aero.begin_trial(model,restored,restored.q,restored.v,.1,flow);assert resumed['mode']=='prescribed';rows.append(dict(case='prescribed_policy_and_both_circulations_restart',bitwise_state=restored.fingerprint()==before,policy_inherited=True))
        rejected=[]
        w=current.wakes[0];wrong_ids=list(w.vertex_ids);wrong_ids[1]=tuple((('other','row'),j) for j in range(len(wrong_ids[1])));badwake=aero.Wake(w.sid,w.vertices,w.mu,tuple(wrong_ids),w.row_ids,w.birth)
        for label,callback in [('window_mismatch',lambda:model.validate_state(replace(current,max_intervals=1))),('missing_accepted_history',lambda:model.validate_state(replace(current,gamma_previous=None,last_h=None))),('material_id_owner',lambda:model.validate_state(replace(current,wakes=(badwake,)))),('checkpoint_checksum',lambda:aero.restore(model,encoded.replace('wp1','broken'),flow))]:
            try:callback()
            except ValueError:rejected.append(label)
            else:raise AssertionError(label)
            assert store.state.fingerprint()==before
        assert len(rejected)==4;rows.append(dict(case='restart_semantics_and_atomic_rejection',negative_controls=rejected))
        # Qualification is attached to the actual load packet and physical
        # state, not just to the circulating coefficient array.
        good=store.begin(current.q,current.v,.1,flow);store.solve(good,flow,ws);bad=copy.deepcopy(good);bad['loads']['Q'][0,0]+=.1
        try:store.accept(bad)
        except ValueError:pass
        else:raise AssertionError('mutated load packet accepted')
        bad=copy.deepcopy(good);bad['state']=replace(bad['state'],v=bad['state'].v+.001)
        try:store.accept(bad)
        except ValueError:pass
        else:raise AssertionError('replaced state accepted')
        assert store.state.fingerprint()==before;rows.append(dict(case='qualification_covers_load_and_state_packet',negative_controls=2))
        n=ws.misses;g=model.geometry(current.q,current.v);ws.assemble(model,g,current.wakes);ws.assemble(model,g,current.wakes);hits=ws.hits;g2=model.geometry(current.q+np.array([0.,0.,.001]),current.v);ws.assemble(model,g2,current.wakes);assert ws.misses>=n+1 and hits>0
        try:model.eps=.1
        except AttributeError:pass
        else:raise AssertionError('model core changed without generation')
        rows.append(dict(case='geometry_cache_and_immutable_core',hits=hits,misses=ws.misses,core_change_rejected=True))
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    paths=[Path(__file__),Path(aero.__file__),ROOT/'rd/aero/wake-loads/experiments/uvlm_wake_assembly.py',Path(binding.__file__),ROOT/'rd/aero/wake-loads/records/uvlm-wake-assembly-result.json'];dest=ROOT/'rd/aero/wake-loads/records/uvlm-wake-policy-result.json';assert not dest.exists();dest.write_text(json.dumps(dict(passed=len(rows)==4 and not failures,records=rows,failures=failures,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8');print(json.dumps(rows,indent=2));print('WAKE_POLICY_RESULT',len(rows)==4 and not failures,flush=True)
if __name__=='__main__':main()
