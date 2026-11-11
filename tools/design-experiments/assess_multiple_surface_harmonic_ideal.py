from pathlib import Path
import json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/multiple-surface-harmonic-ideal'
def main():
    paths=[];records=[]
    for d in [0,5]:
        for m,n in [(2,8),(4,16),(8,32)]:
            p=[OUT/f'{side}-d{d}-m{m}-n{n}.npz' for side in ['candidate','reference']];paths+=p;a,b=[np.load(q) for q in p];errors={key:float(np.linalg.norm(a[key]-b[key])/np.linalg.norm(b[key])) for key in ['A','gamma','dual']};records.append(dict(dihedral=d,m=m,n=n,errors=errors,passed=bool(max(errors.values())<1e-9)))
    paths+=[Path(__file__),ROOT/'tools/design-experiments/run_multiple_surface_harmonic_ideal.py',ROOT/'tools/design-experiments/run_multiple_surface_harmonic.py',ROOT/'docs/verification/ptera-core-model-diagnostic.md',ROOT/'docs/verification/multiple-surface-harmonic-result.json'];result=dict(passed=all(r['passed'] for r in records),scope='same ideal kernel; noncontact harmonic interaction; finite-epsilon cores remain different',records=records,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest=ROOT/'docs/verification/multiple-surface-harmonic-ideal-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
