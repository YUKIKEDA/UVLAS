import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,copy,numpy as np
from scipy.linalg import subspace_angles
from flutter_search_witness import Search,canonical
ROOT=rd_paths.ROOT
def main():
    a=Search();a.run();records=[]
    for j,truth in enumerate([10.,14.]):
        b=a.state['brackets'][j];assert b['lo']<truth<b['hi'] and abs(b['estimate']-truth)<.02;records.append(dict(case='error_aware_critical_interval',mode=j,truth=truth,interval=[b['lo'],b['hi']],estimate=b['estimate'],frequency=b['frequency']))
    for stop in [3,8,12,a.state['actions']-1]:
        b=Search()
        for _ in range(stop):b.step()
        resumed=Search.restore(b.encode());resumed.run();assert resumed.encode()==a.encode()
    records.append(dict(case='search_restart_queue_modes_and_points_bitwise',restart_actions=[3,8,12,a.state['actions']-1],evaluations=a.state['evaluations'],unique_points=len(a.state['completed'])))
    p9=a.state['completed'][float(9).hex()];p12=a.state['completed'][float(12).hex()]
    def modes(p):
        x=np.array(p['modes']);return x[...,0]+1j*x[...,1]
    MAC=abs(modes(p9).conj().T@modes(p12))**2;assert MAC[0,0]>.95 and MAC[0,1]<.05 and p9['frequency'][0]<p9['frequency'][1] and p12['frequency'][0]>p12['frequency'][1];records.append(dict(case='frequency_order_crossing_keeps_physical_mode',MAC=MAC.tolist()))
    Q=np.eye(3)[:,:2];t=.6;internal=np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]]);angles=subspace_angles(Q,Q@internal);assert max(abs(angles))<1e-12;rank_rejected=np.linalg.matrix_rank(Q[:,:1])!=2;assert rank_rejected;records.append(dict(case='cluster_subspace_independent_of_internal_basis',max_angle=float(max(angles)),missing_cluster_dimension_rejected=True))
    p=a.state['completed'][float(10).hex()] if float(10).hex() in a.state['completed'] else a.evaluate(10.);assert a.sign(p,0)==0;records.append(dict(case='uncertain_growth_not_forced_to_stable_or_unstable'))
    rejected=[]
    for label,change in [('policy',lambda s:s.update(h=.03)),('missing_queue',lambda s:s.pop('queue')),('point_identity',lambda s:next(iter(s['completed'].values())).update(U=1.234)),('premature_complete',lambda s:s['brackets'][0].update(done=False))]:
        s=copy.deepcopy(a.state);change(s);encoded=canonical(dict(state=s,checksum=hashlib.sha256(canonical(s).encode()).hexdigest()))
        try:Search.restore(encoded)
        except ValueError:rejected.append(label)
        else:raise AssertionError(label)
    records.append(dict(case='invalid_search_state_rejections',rejected=rejected));paths=[Path(__file__),ROOT/'rd/coupled/j1/experiments/flutter_search_witness.py',ROOT/'docs/specs/coupled/coupled-state-search-contract.md'];result=dict(passed=True,scope='search control and mode identities on a fully known two-mode spectrum; not physical flutter validation',records=records,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest=ROOT/'rd/coupled/j1/records/flutter-search-control-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
