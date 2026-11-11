"""Independent rule/basis and numerical traction audit for local mixed states."""
from pathlib import Path
import hashlib,json,time
import numpy as np
import hybrid_gradient_ancf as f
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1e-30,np.linalg.norm(b)))
def main():
    parent_path=ROOT/'docs/verification/hybrid-gradient-local-result.json';parent=json.loads(parent_path.read_text());dest=ROOT/'docs/verification/hybrid-gradient-local-audit.json';assert not dest.exists()
    for p,h in parent['sources'].items():assert digest(ROOT/p)==h
    start=time.monotonic();records=[];stops=[]
    for r in parent['records']:
        try:
            assert digest(ROOT/r['artifact'])==r['artifact_sha256'];s=np.load(ROOT/r['artifact']);state=f.prepare(s['q'],r['h'],r['kind'],15,9,True)
            U,force,tangent,info=f.nonlinear(state,s['current'],s['v1'],s['v2'])
            A,_,_=f.generators(s['q'],r['h'],state['p'],r['kind'],True);primary_gradient=(A@s['coeff'])@s['D']
            metrics=dict(energy=rel(U-s['energy'],s['energy']),force=rel(force-s['force'],s['force']),tangent=rel(tangent-s['tangent'],s['tangent']),
                weak_gradient=rel(state['Gw']-primary_gradient,state['Gw']),constant_P_work=rel(np.einsum('p,pjn->nj',state['weight'],state['Gw'])-f.boundary_map(state,15,9,True),f.boundary_map(state,15,9,True)))
            fluxforce=np.zeros((135,3))
            for face,traction in zip(state['faces'],info['fluxes']):fluxforce[face['offset']:face['offset']+15]=np.einsum('p,pn,pa->na',face['measure'],face['hat'],traction)
            metrics['flux_dual']=rel(force.reshape(135,3)[75:]-fluxforce[75:],force.reshape(135,3)[75:])
            records.append(dict(kind=r['kind'],geometry=r['geometry'],h=r['h'],metrics=metrics,passed=all(v<=1e-8 for v in metrics.values())))
        except Exception as exc:stops.append(dict(kind=r['kind'],geometry=r['geometry'],h=r['h'],reason=repr(exc)))
    sources=[Path(__file__),Path(f.__file__),parent_path,ROOT/'docs/verification/hybrid-gradient-ancf-protocol.md']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_reproduction_checks=len(records)==len(parent['records']) and all(r['passed'] for r in records) and not stops,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,physical_adoption=False),indent=2),encoding='utf-8')
    print('audited',len(records),'failed',[(r['kind'],r['geometry'],r['h'],r['metrics']) for r in records if not r['passed']],'stops',stops)
if __name__=='__main__':main()
