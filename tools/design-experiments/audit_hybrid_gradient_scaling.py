"""Verify the thickness scaling of the fixed, unchanged trace stabilization."""
from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    dest=ROOT/'docs/verification/hybrid-gradient-scaling-audit.json';assert not dest.exists()
    a=np.array([.2,-.3,.1]);b=np.array([-.1,.4,.2]);c=np.array([.6,.2,-.3]);p,w=np.polynomial.legendre.leggauss(9);records=[]
    for h in [.01,.001,.0001]:
        z=h*p/2;value=a+z[:,None]*b+.5*z[:,None]**2*c;observed=h/2*w@np.sum(value*value,axis=1)
        exact=h*np.linalg.norm(a+h*h*c/24)**2+h**3*np.linalg.norm(b)**2/12+h**5*np.linalg.norm(c)**2/720
        records.append(dict(h=h,integral=float(observed),closed_form=float(exact),relative_difference=float(abs(observed/exact-1))))
    result=json.loads((ROOT/'docs/verification/hybrid-gradient-refinement-result.json').read_text());r=result['records']
    physical=r[1]['observed']['elastic_error'][0]/r[0]['observed']['elastic_error'][0]
    sources=[Path(__file__),ROOT/'docs/verification/hybrid-gradient-refinement-result.json',ROOT/'docs/hybrid-gradient-ancf-formulation.md']
    dest.write_text(json.dumps(dict(records=records,all_checks_pass=all(r['relative_difference']<=1e-12 for r in records),
        total_error_amplification=result['thin_amplification'][0],physical_reconstructed_strain_error_amplification=physical,
        claim='For bounded c and scaled bending energy with fixed positive tau, r trace jump is O(h); no full ANCF impossibility or unique causal proof is claimed.',
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources}),indent=2),encoding='utf-8')
    print('scaling identities',records,'physical strain amplification',physical)
if __name__=='__main__':main()
