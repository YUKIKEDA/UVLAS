import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,time,copy,numpy as np
from scipy.linalg import solve
from pg_uvlm_free_step import FreeProblem
import uvlm_wake_assembly_v2 as aero
from coupled_checkpoint_witness import packet,seal,restore,Store,physical_next
ROOT=rd_paths.ROOT
def main():
    start=time.monotonic();p=FreeProblem();s0=p.sstar.copy();s0[2:4]=[.002,-.0005];parent=p.decode(s0);g=p.model.geometry(parent.q,parent.v);ws=aero.Workspace();A,W=ws.assemble(p.model,g,parent.wakes);mu=np.concatenate([w.mu.ravel() for w in parent.wakes]);mu[:p.n]=0;rhs=np.sum(g['n']*(g['CV']-p.flow.value(g['P'],0.)),axis=1)-W@mu;s0[8:8+p.nb]=solve(A,rhs);physical0=p.decode(s0);s1=p.advance(s0)[0];physical1=physical_next(p,physical0,s1);encoded=seal(packet(p,s1,s0,physical1));r,phys,payload=restore(p,encoded);assert np.array_equal(r,s1);s2=p.advance(s1)[0];resumed=p.advance(r)[0];assert np.array_equal(s2,resumed);physical2=physical_next(p,physical1,s2);a=packet(p,s2,s1,physical2);b=packet(p,resumed,r,physical_next(p,phys,resumed));assert seal(a)==seal(b);legacy=phys.gamma_dot();difference=float(np.linalg.norm(legacy-np.array(payload['gamma_dot'])));assert difference>1e-8;rejected=[]
    for label,change in [('missing_history',lambda x:x['stencil']['gamma'].pop()),('BDF1_substitution',lambda x:x.update(gamma_dot=legacy.tolist())),('wrong_policy',lambda x:x.update(h=.5*p.h)),('wrong_model',lambda x:x.update(model_id='other')),('changed_load',lambda x:x['loads']['Q'][0].__setitem__(0,x['loads']['Q'][0][0]+.1))]:
        bad=copy.deepcopy(payload);change(bad)
        try:restore(p,seal(bad))
        except ValueError:rejected.append(label)
        else:raise AssertionError(label)
    store=Store(encoded);trial=store.qualify(p,a);bad=copy.deepcopy(trial);bad['encoded']=bad['encoded'].replace('physical-velocity','other-velocity');before=store.encoded
    try:store.commit(p,bad)
    except ValueError:rejected.append('changed_qualified_trial')
    else:raise AssertionError('changed_qualified_trial')
    assert store.encoded==before;store.commit(p,trial);after=store.encoded
    try:store.commit(p,trial)
    except ValueError:rejected.append('double_commit')
    else:raise AssertionError('double_commit')
    assert store.encoded==after;paths=[Path(__file__),ROOT/'rd/coupled/j1/experiments/coupled_checkpoint_witness.py',ROOT/'docs/specs/coupled/coupled-state-search-contract.md',ROOT/'rd/coupled/j1/records/pg-uvlm-free-step-result.json'];result=dict(passed=True,scope='actual free PG/UVLM physical IDs/time plus BDF2 stencil/load envelope and atomic commit',next_state_and_packet_bitwise=True,BDF1_would_change_current_derivative=difference,rejected=rejected,seconds=time.monotonic()-start,sources={q.relative_to(ROOT).as_posix():hashlib.sha256(q.read_bytes()).hexdigest() for q in paths});dest=ROOT/'rd/coupled/j1/records/coupled-checkpoint-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
