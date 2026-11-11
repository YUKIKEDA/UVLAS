"""Recalculate local gate from archived matrices using a separate solve path."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import eigvalsh
import check_q9_curved_load as original
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    start=time.monotonic();path=F/'q9-shear-projection-gate-result.json';data=json.loads(path.read_text());records=[]
    assert data['all_criteria_pass'] and data['stop_reason'] is None and len(data['cases'])==3
    assert data['source_sha256']==digest(ROOT/'tools/design-experiments/check_q9_shear_projection_gate.py')
    assert data['candidate_sha256']==digest(ROOT/'docs/q9-shear-projection-candidate.md')
    assert data['protocol_sha256']==digest(F/'q9-shear-projection-gate-experiment.md')
    for n,h in data['dependency_sha256'].items():assert digest(ROOT/'tools/design-experiments'/n)==h
    for c in data['cases']:
        assert all(c['checks'].values()) and digest(ROOT/c['artifact'])==c['artifact_sha256'];d=np.load(ROOT/c['artifact']);M,H,Kn,Kp,Kf=[d[n] for n in ['M','H','Kn','Kp','Kfull']]
        scale=1/np.sqrt(Kf.diagonal());transform=lambda k:scale[:,None]*k*scale[None,:];K=transform(Kp)
        ms=1/np.sqrt(M.diagonal());a=ms[:,None]*M*ms[None,:];ev=eigvalsh((a+a.T)/2);assert ev[0]>0 and ev[-1]/ev[0]<=1e8
        candidate=Kn+H.T@(ms[:,None]*np.linalg.solve(a,ms[:,None]*H));condensation=float(np.linalg.norm(transform(candidate-Kp))/np.linalg.norm(K));assert condensation<=1e-10
        e=eigvalsh((K+K.T)/2);singular=np.sort(abs(e))[::-1];ranks={str(t):int(np.count_nonzero(singular>t*singular[0])) for t in [1e-11,1e-10,1e-9]};assert ranks==c['rank'] and all(v==75 for v in ranks.values())
        delta=transform(Kf-Kp);de=eigvalsh((delta+delta.T)/2);assert e[0]>=-1e-10*max(1.,abs(e).max()) and de[0]>=-1e-10*max(1.,abs(de).max())
        assert np.linalg.norm(K-K.T)/np.linalg.norm(K)<=1e-12
        rigid=max(float(np.linalg.norm(K@(u/scale))/(np.linalg.norm(K,2)*np.linalg.norm(u/scale))) for u in d['rigid']);assert rigid<=1e-11
        patch=max(float(np.linalg.norm(scale*((Kp-Kf)@u))/np.linalg.norm(scale*(Kf@u))) for u in d['patch_fields']);assert patch<=1e-10
        reference=original.evaluate(d['q'],.1,9)['stiffness'];baseline=float(np.linalg.norm(transform(reference-Kf))/np.linalg.norm(transform(Kf)));assert baseline<=1e-10
        records.append(dict(case=c['case'],rank=ranks,condensation_difference=condensation,patch_maximum=patch,rigid=rigid,independent_full_matrix_difference=baseline))
        assert time.monotonic()-start<=120 and memory.memory_bytes()<=512*1024**2
    result=dict(source_sha256=digest(__file__),result_sha256=digest(path),audit_protocol_sha256=digest(F/'q9-shear-projection-audit-protocol.md'),
        dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [original,memory]},cases=records,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=True,product_shell_adopted=False)
    (F/'q9-shear-projection-gate-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_criteria_pass=True,cases=records)))
if __name__=='__main__':main()
