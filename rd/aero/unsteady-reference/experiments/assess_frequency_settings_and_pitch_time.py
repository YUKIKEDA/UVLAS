import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,numpy as np
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/frequency-settings';OLD=ROOT/'build/analysis/pitch-frequency-2-3';USED=set()
def load(p):USED.add(p);return np.load(p)['transfer']
def main():
    records=[]
    for side in ['candidate','reference']:
        for k in [.05,.2,.8]:
            base=load(OUT/f'{side}-base-k{k}.npz');d={label:abs((load(OUT/f'{side}-{label}-k{k}.npz')-base)/base).tolist() for label in ['tail64','tail256','core1e-4','core1e-6','q24']};records.append(dict(side=side,k=k,relative_changes=d))
    def time(samples,cycles=3):
        p=OLD/f'pitch-time-s{samples}-c{cycles}.json';USED.add(p);USED.add(p.with_suffix('.npz'));a=np.array(json.loads(p.read_text())['transfer']);return a[:,0]+1j*a[:,1]
    ref=load(OLD/'reference-m2-n8-k0.2.npz')[1];ts=[time(n) for n in [256,512,1024]];period=abs((time(512,6)-ts[1])/ref);extra=[abs((2*ts[j+1]-ts[j]-ref)/ref) for j in range(2)];passed=max(extra[-1])<.001 and max(period)<.0001;paths=[Path(__file__),ROOT/'rd/aero/unsteady-reference/experiments/run_frequency_setting_controls.py',ROOT/'rd/aero/unsteady-reference/experiments/probe_ptera_pitch_time.py',ROOT/'rd/aero/unsteady-reference/records/frequency-setting-controls-protocol.md',ROOT/'rd/aero/unsteady-reference/records/uvlm-pitch-frequency-extension-protocol.md'];USED.update(paths);result=dict(passed=bool(passed),scope='setting sensitivities measured separately; native pitch input/axis/sign and continuous-time frequency provenance',settings_records=records,pitch_time_finite_errors=[abs((x-ref)/ref).tolist() for x in ts],pitch_time_extrapolated_errors=[x.tolist() for x in extra],period_3_to_6=period.tolist(),coarse_sensitivity_is_not_fine_grid_error_bound=True,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(USED)});dest=ROOT/'rd/aero/unsteady-reference/records/frequency-settings-and-pitch-time-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k!='sources'},indent=2))
if __name__=='__main__':main()
