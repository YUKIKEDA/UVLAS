import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,numpy as np
ROOT=rd_paths.ROOT
def main():
    records=[];paths=[]
    for d in [0,5]:
        for m,n in [(2,8),(4,16),(8,32)]:
            data={}
            for side in ['candidate','reference']:
                for q in [32,48]:
                    p=ROOT/f'build/analysis/multiple-surface-harmonic-ideal-q{q}/{side}-d{d}-m{m}-n{n}.npz';paths.append(p);data[side,q]=np.load(p)
            errors={key:float(np.linalg.norm(data['candidate',48][key]-data['reference',48][key])/np.linalg.norm(data['reference',48][key])) for key in ['A','gamma','dual']};changes={side:{key:float(np.linalg.norm(data[side,48][key]-data[side,32][key])/np.linalg.norm(data[side,48][key])) for key in ['A','gamma','dual']} for side in ['candidate','reference']};passed=max(errors.values())<1e-9 and max(v for r in changes.values() for v in r.values())<1e-9;records.append(dict(dihedral=d,m=m,n=n,errors=errors,quadrature_changes=changes,passed=bool(passed)))
    paths+=[Path(__file__),ROOT/'rd/aero/multi-surface/experiments/run_multiple_surface_harmonic_refined.py',ROOT/'rd/aero/multi-surface/records/multiple-surface-harmonic-quadrature-extension.md',ROOT/'rd/aero/multi-surface/records/multiple-surface-harmonic-ideal-result.json'];result=dict(passed=all(r['passed'] for r in records),scope='same ideal kernel and independently refined quadrature; harmonic four-input multiple lifting surfaces',records=records,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest=ROOT/'rd/aero/multi-surface/records/multiple-surface-harmonic-refined-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
