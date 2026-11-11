"""Separate rule/basis audit including reconstruction-induced traction dual."""
from pathlib import Path
import hashlib,json,time
import numpy as np
import polynomial_trace_ancf as f
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1e-30,np.linalg.norm(b)))
def main():
    source=ROOT/'docs/verification/polynomial-trace-local-result.json';parent=json.loads(source.read_text());dest=ROOT/'docs/verification/polynomial-trace-local-audit.json';assert not dest.exists()
    start=time.monotonic();records=[];stops=[]
    for r in parent['records']:
        try:
            assert digest(ROOT/r['artifact'])==r['artifact_sha256'];s=np.load(ROOT/r['artifact']);state=f.prepare(s['q'],r['h'],r['kind'],15,9,True)
            primary=f.prepare(s['q'],r['h'],r['kind']);U,force,tangent,info=f.nonlinear(state,s['current'],s['v1'],s['v2'])
            metrics=dict(energy=rel(U-s['energy'],s['energy']),force=rel(force-s['force'],s['force']),tangent=rel(tangent-s['tangent'],s['tangent']),weak_gradient=rel(state['Gw']-f.evaluate_gradient(primary,state['p']),state['Gw']))
            flux=np.zeros((135,3))
            for face,t in zip(state['faces'],info['fluxes']):flux[face['offset']:face['offset']+15]=np.einsum('p,pn,pa->na',face['measure'],face['hat'],t)
            metrics['flux_dual']=rel(flux[75:]-force.reshape(135,3)[75:],force.reshape(135,3)[75:])
            jump_errors=[]
            for face,old,(axis,sign) in zip(state['faces'],primary['faces'],f.FACES):
                L=f.face_legendre(face['points'],axis);_,_,phi,_=f.old.face_data(s['q'],r['h'],axis,sign,15,9,True)
                J=np.zeros((len(phi),135));J[:,:75]=phi;J[:,face['offset']:face['offset']+15]=-phi[:,f.edge_nodes(axis,sign)]
                J+=L@old['lift_coeff']@primary['Hhi']@primary['D'];jump_errors.append(rel(J-face['jump'],face['jump']))
            metrics['stabilization_operator']=max(jump_errors)
            polynomial=f.polynomial_patch(state) if r['geometry']=='flat' else None
            records.append(dict(kind=r['kind'],geometry=r['geometry'],h=r['h'],metrics=metrics,polynomial=polynomial,passed=all(v<=1e-8 for v in metrics.values()) and (polynomial is None or polynomial['passed'])))
        except Exception as exc:stops.append(dict(kind=r['kind'],geometry=r['geometry'],h=r['h'],reason=repr(exc)))
    sources=[Path(__file__),Path(f.__file__),source,ROOT/'docs/verification/polynomial-trace-protocol.md']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_reproduction_checks=len(records)==len(parent['records']) and all(r['passed'] for r in records) and not stops,
        physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
    print('audit',len(records),'failed',[r for r in records if not r['passed']],'stops',stops)
if __name__=='__main__':main()
