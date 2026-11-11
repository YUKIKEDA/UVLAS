"""Observation audit beyond the predeclared lift/pitch sensitivity criteria."""
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[2]
FOLDER=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    parent=FOLDER/'avl-source-measurement-corrected-result.json'
    result=json.loads(parent.read_text())
    cases=[dict(variant=c['variant'],case=c['case'],m=c['m'],n=c['n'],alpha=c['alpha'],
        side_coefficient=c['coefficients']['force'][1],yaw_coefficient=c['coefficients']['moment'][2]) for c in result['cases']]
    maxima={v:dict(side=max(abs(c['side_coefficient']) for c in cases if c['variant']==v),
        yaw=max(abs(c['yaw_coefficient']) for c in cases if c['variant']==v)) for v in ['baseline','measured','core0']}
    out=dict(kind='post_result_observation_not_revised_acceptance',source_sha256=digest(Path(__file__)),
        parent_result_sha256=digest(parent),cases=cases,maxima=maxima,
        initial_fourteen_checks_unchanged=result['checks'],
        core0_reference_adoption_eligible=False,
        unresolved='core-zero results break lateral reflection symmetry; near-collinear force-evaluation kernel is a hypothesis, not yet isolated')
    (FOLDER/'avl-core-symmetry-observation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(maxima))

if __name__=='__main__':main()
