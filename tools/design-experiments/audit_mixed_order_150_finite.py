from pathlib import Path
import hashlib,json,time
import numpy as np
import mixed_order_150_finite as f
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    path=ROOT/'docs/verification/mixed-order-150-finite-result.json';parent=json.loads(path.read_text());dest=ROOT/'docs/verification/mixed-order-150-finite-audit.json';assert not dest.exists()
    for p,h in parent['sources'].items():assert digest(ROOT/p)==h
    start=time.monotonic();records=[]
    for r in parent['records']:
        assert digest(ROOT/r['artifact'])==r['artifact_sha256'];s=np.load(ROOT/r['artifact'])
        mesh=dict(ref=s['ref'],cells=s['cells'],xy=[s[f'xy{i}'] for i in range(3)],free=[s[f'free{i}'] for i in range(3)])
        # Cardinal primary -> Vandermonde audit. Tensor stress norm and surface
        # metric are evaluated afresh, with no condensed matrix or saved QR.
        observed=f.observe(mesh,r['h'],s['u'],False)
        error_difference=abs(observed['error']-r['observed']['error']);h1_difference=abs(observed['h1']/r['observed']['h1']-1)
        records.append(dict(n=r['n'],h=r['h'],observed=observed,error_difference=error_difference,h1_difference=h1_difference,
            passed=error_difference<=1e-8 and h1_difference<=1e-8))
    sources=[Path(__file__),Path(f.__file__),path]
    dest.write_text(json.dumps(dict(records=records,all_reproduction_checks=bool(records) and all(r['passed'] for r in records),
        primary_complete=parent['all_numerical_checks'],independent_global_resolve=False,physical_adoption=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
    print(json.dumps(records,indent=2))
if __name__=='__main__':main()
