"""Supplementary predeclared r=M*delta diagnosis; original failures retained."""
import hashlib
import json
import time
from pathlib import Path
import check_uvlm_harmonic_heave as trial

ROOT=Path(__file__).resolve().parents[2]


def main():
    rows=[];changes=[]
    for k in [.1,.5,1.]:
        cases={}
        for r in [.4,.2,.1]:
            for m in [16,32,64,128]:
                row,_,_=trial.harmonic(m,r/m,80.,k)
                row['cell_ratio']=r;rows.append(row);cases[(r,m)]=row
                print(json.dumps(row),flush=True)
        fine,_,_=trial.harmonic(128,.1/128,160.,k)
        fine['cell_ratio']=.1;rows.append(fine)
        changes.append(dict(k=k,joint_refinement=trial.relative(cases[(.1,64)],cases[(.1,128)]),
            time_refinement=trial.relative(cases[(.2,128)],cases[(.1,128)]),
            tail=trial.relative(cases[(.1,128)],fine)))
        trial.CACHE.clear()
    fine=[r for r in rows if r['length']==160.]
    checks=dict(complex_accuracy=all(r['complex_relative']<=.02 for r in fine),
        amplitude_accuracy=all(r['amplitude_relative']<=.02 for r in fine),
        phase_accuracy=all(abs(r['phase_degrees'])<=1 for r in fine),
        joint_change=all(c['joint_refinement']<=.005 for c in changes),
        time_change=all(c['time_refinement']<=.01 for c in changes),
        tail_change=all(c['tail']<=.002 for c in changes),
        nonpenetration=all(r['nonpenetration_relative']<=1e-11 for r in rows))
    source=Path(__file__)
    names=['check_uvlm_harmonic_heave.py',*json.loads((ROOT/'docs/verification/uvlm-harmonic-heave-result.json').read_text())['dependency_sha256']]
    result=dict(protocol='docs/verification/uvlm-harmonic-heave-experiment.md#7',
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={name:hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() for name in names},
        geometry_id=trial.geometry.GEOMETRY,harmonic_cases=rows,refinement_changes=changes,checks=checks,
        all_criteria_pass=bool(all(checks.values())),seconds=time.monotonic()-trial.START,
        peak_process_working_bytes=trial.memory.memory_bytes(),numpy=trial.np.__version__,scipy=trial.scipy.__version__,
        limits='supplementary joint refinement, not a product timestep criterion; original failed comparison retained')
    (ROOT/'docs/verification/uvlm-heave-resolution-result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,changes=changes,seconds=result['seconds'])),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':main()
