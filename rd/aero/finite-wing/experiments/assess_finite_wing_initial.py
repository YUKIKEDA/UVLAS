"""Preserve first finite campaign, including unmet reference convergence."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np,json,hashlib
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/finite-wing-2-3'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    records=[];artifacts={};failures=[]
    for d in [0,5]:
        data={}
        for side in ['candidate','reference']:
            for m,n in [(4,12),(4,24),(8,24),(8,48),(16,48),(16,96),(32,96),(16,192),(32,192)]:
                p=OUT/f'{side}-d{d}-m{m}-n{n}.json';r=json.loads(p.read_text());data[side,m,n]=r;artifacts[p.relative_to(ROOT).as_posix()]=sha(p);snap=p.with_suffix('.npz');assert sha(snap)==r['snapshot_sha256'];artifacts[snap.relative_to(ROOT).as_posix()]=sha(snap)
            a=data[side,32,96];b=data[side,32,192];lift=abs(a['F'][2]/b['F'][2]-1);moment=abs(a['M'][1]/b['M'][1]-1);bandmax=float(np.max(abs(np.array(a['bands'])/b['bands']-1)));bandl2=float(np.linalg.norm(np.array(a['bands'])-b['bands'])/np.linalg.norm(b['bands']))
            records.append(dict(case='finite_grid_uncertainty',side=side,dihedral=d,lift_adjacent_difference=lift,moment_adjacent_difference=moment,band_max_adjacent_difference=bandmax,band_l2_adjacent_difference=bandl2,space_budget=.003,passed=max(lift,moment)<=.003))
            if max(lift,moment)>.003:failures.append(f'{side} d{d}: spatial budget exceeded')
        a=data['candidate',32,192];b=data['reference',32,192];records.append(dict(case='independent_code_difference',dihedral=d,lift=abs(a['F'][2]/b['F'][2]-1),moment=abs(a['M'][1]/b['M'][1]-1),bands=float(np.max(abs(np.array(a['bands'])/b['bands']-1)))))
    sources={p.relative_to(ROOT).as_posix():sha(p) for p in [Path(__file__),ROOT/'rd/aero/finite-wing/experiments/run_finite_wing_campaign.py',ROOT/'rd/aero/finite-wing/experiments/uvlm_steady_reference.py',ROOT/'rd/aero/unsteady-reference/experiments/probe_aerosandbox_reference.py',ROOT/'rd/aero/unsteady-reference/records/uvlm-independent-reference-protocol.md']}
    result=dict(passed=not failures,failures=failures,records=records,sources=sources,artifacts=artifacts,reference_qualified=False)
    dest=ROOT/'rd/aero/finite-wing/records/finite-wing-initial-campaign-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(passed=result['passed'],records=records),indent=2))
if __name__=='__main__':main()
