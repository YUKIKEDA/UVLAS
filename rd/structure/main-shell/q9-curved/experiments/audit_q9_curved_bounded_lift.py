"""Check saved physical-coordinate H1 balls, residuals and KKT equations."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,time
import numpy as np
import check_q9_curved_bounded_lift as core
ROOT=core.ROOT;F=core.F
def main():
    target=F/'q9-curved-bounded-lift-audit.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-curved-bounded-lift-result.json';r=json.loads(parent.read_text());assert len(r['cases'])==22 and r['stop_reason'] is None
    assert r['source_sha256']==core.digest(core.__file__)
    for name,sha in r['dependencies'].items():assert core.digest(rd_paths.EXPERIMENTS/name)==sha
    start=time.monotonic();cases=[]
    for case in r['cases']:
        assert core.digest(ROOT/case['artifact'])==case['artifact_sha256'];a=np.load(ROOT/case['artifact']);A,F0,M=a['A'],a['F'],a['M'];maximum=dict(residual=0.,h1=0.,kkt=0.,constraint=0.)
        for index,p in enumerate(a['bounded']):
            k=index//2;fit=case['rows'][k]['bounded'][index%2];r0=A@p+F0[:,k];norm=float(np.linalg.norm(r0));P=p.reshape(len(M),6);MP=M@P;h1=float(np.sqrt(np.sum(P*MP)));grad=A.T@r0+fit['multiplier']*MP.reshape(-1);den=np.linalg.norm(A.T@F0[:,k])+fit['multiplier']*np.linalg.norm(MP)
            maximum['residual']=max(maximum['residual'],abs(norm-fit['residual'])/max(np.linalg.norm(F0[:,k]),1e-30));maximum['h1']=max(maximum['h1'],abs(h1-fit['h1'])/max(h1,fit['h1'],1e-30));maximum['kkt']=max(maximum['kkt'],float(np.linalg.norm(grad)/max(den,1e-30)));maximum['constraint']=max(maximum['constraint'],max(0.,h1/fit['bound']-1))
        checks=dict(residual=maximum['residual']<=1e-10,h1=maximum['h1']<=1e-9,kkt=maximum['kkt']<=1e-9,constraint=maximum['constraint']<=1e-10);cases.append(dict(geometry=case['geometry'],n=case['n'],h=case['h'],maximum=maximum,checks={k:bool(v) for k,v in checks.items()}))
    result=dict(source_sha256=core.digest(__file__),parent_sha256=core.digest(parent),protocol_sha256=core.digest(F/'q9-curved-bounded-lift-fine-experiment.md'),dependency_sha256=core.digest(core.__file__),cases=cases,seconds=time.monotonic()-start,all_criteria_pass=bool(len(cases)==22 and all(all(c['checks'].values()) for c in cases)),coarse_observation_accepted=False,product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(dict(all_criteria_pass=result['all_criteria_pass'],maximum={key:max(c['maximum'][key] for c in cases) for key in cases[0]['maximum']},seconds=result['seconds'])))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
