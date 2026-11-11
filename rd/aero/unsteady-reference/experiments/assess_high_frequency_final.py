import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,numpy as np
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/high-frequency-final-qualified'
def main():
    parent=ROOT/'rd/aero/unsteady-reference/records/ideal-pitch-frequency-result.json';old=json.loads(parent.read_text());meta={};values={};paths=[parent,Path(__file__),ROOT/'rd/aero/unsteady-reference/records/high-frequency-final-space-core-protocol.md',ROOT/'rd/aero/unsteady-reference/records/high-frequency-residual-evaluation-policy.md',ROOT/'rd/aero/unsteady-reference/records/high-frequency-final-symmetry-proof.json',ROOT/'rd/aero/unsteady-reference/records/high-frequency-symmetry-residual-diagnostic.json',ROOT/'rd/aero/unsteady-reference/experiments/run_high_frequency_final_qualified.py',ROOT/'rd/aero/unsteady-reference/experiments/run_high_frequency_final_bounded.py']
    for m in [128,256]:
        p=OUT/f'candidate-m{m}-n192-core1e-6.json';paths.extend([p,p.with_suffix('.npz')]);meta[m]=json.loads(p.read_text());values[m]=np.load(p.with_suffix('.npz'))['transfer'];assert meta[m]['seconds']<1200 and meta[m]['peak_bytes']<24*2**30 and meta[m]['residual']<1e-11 and max(meta[m]['uncondensed_sample_residuals'])<1e-11
    records=[]
    for r in old['records']:
        new=dict(r)
        if r['k']==.8:
            a=['heave','pitch'].index(r['motion']);b=['Fz','My'].index(r['observable']);limit=complex(*r['reference_limit']);value=values[256][a,b];difference=abs(value-limit)/abs(limit);combined=difference+r['reference_relative_uncertainty'];coarse_original=complex(*r['candidate']);new.update(candidate=[float(value.real),float(value.imag)],finite_candidate_to_reference_limit=float(difference),combined_relative_indicator=float(combined),passed=bool(combined<.01 and r['reference_relative_uncertainty']<.0025),candidate_chordwise_panels=256,candidate_spanwise_panels=192,candidate_core=1e-6,same_grid_core_change_at_m128=float(abs(values[128][a,b]-coarse_original)/abs(limit)),same_core_chord_change_128_to_256=float(abs(values[256][a,b]-values[128][a,b])/abs(limit)))
        else:new.update(candidate_chordwise_panels=128,candidate_spanwise_panels=192,candidate_core=1e-5)
        records.append(new)
    result=dict(passed=all(r['passed'] for r in records),scope='finite candidate values versus ideal-kernel independent extrapolated reference; three frequencies/two motions/two total loads; closed wake128m',records=records,initial_12_observation_result_remains_failed=True,rigorous_error_bound=False,resource_records={str(m):meta[m] for m in meta},sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest=ROOT/'rd/aero/unsteady-reference/records/high-frequency-final-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(passed=result['passed'],max_combined_indicator=max(r['combined_relative_indicator'] for r in records),high_frequency_records=[r for r in records if r['k']==.8]),indent=2))
if __name__=='__main__':main()
