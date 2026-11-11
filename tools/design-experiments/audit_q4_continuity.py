"""Re-evaluate stored negative controls, without certifying unsolved fits."""
from pathlib import Path
import hashlib,json,time
import numpy as np
import audit_tying_projection_v2 as independent
import diagnose_q4_continuity_v2 as primary
m=independent.m
ROOT=Path(__file__).resolve().parents[2]

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    parent_path=ROOT/'docs/verification/q4-continuity-v2-result.json';parent=json.loads(parent_path.read_text())
    dest=ROOT/'docs/verification/q4-continuity-v2-audit.json';assert not dest.exists()
    start=time.monotonic();records=[];stop=None
    try:
        for r in parent['records']:
            assert digest(ROOT/r['artifact'])==r['artifact_sha256'];saved=np.load(ROOT/r['artifact'])
            mesh=primary.relax(m.mesh(.4,2,.2),r['broken_fields']);p,w=m.quadrature(11,2);totals={kind:dict(error=0.,reference=0.,h1=0.) for kind in r['solutions']}
            area=0.;kernel2=matrix2=0.;gauge_error=0.
            for ids in mesh['cells']:
                q=mesh['ref'][ids].reshape(225);J,B0,B1=independent.fields(q,p,r['method']);weight=w*np.linalg.det(J)
                nodes=list(saved['free_rd']);lookup={int(node):i for i,node in enumerate(nodes)}
                local_kernel=np.zeros((225,saved['kernel'].shape[1]))
                for j,node in enumerate(ids):
                    if int(node) in lookup:local_kernel[3*j:3*j+3]=saved['kernel'][3*lookup[int(node)]:3*lookup[int(node)]+3]
                kernel2+=np.sum(weight[:,None,None]*(B0@local_kernel)**2);matrix2+=np.sum(weight[:,None,None]*B0**2)
                Nr,gr=m.shape(p,0,True);surface=np.einsum('na,pnb->pab',mesh['ref'][ids[:25]],gr);G=surface.transpose(0,2,1)@surface
                aw=w*np.sqrt(np.linalg.det(G));area+=aw.sum();x=(Nr@mesh['ref'][ids[:25]])[:,0]
                # The removed gauges must be pure element translations in both maps.
                if 0 in r['broken_fields']:
                    translations=np.zeros((225,3));translations[:75]=np.tile(np.eye(3),(25,1))
                    gauge_error=max(gauge_error,float(np.linalg.norm(B0@translations)),float(np.linalg.norm(B1@translations)))
                for kind,t in totals.items():
                    u=saved['u_'+kind][ids];observed=np.einsum('kij,pk->pij',m.MB,B1@u.reshape(225));target=np.einsum('kij,pk->pij',m.MB,m.target(kind,.4,x))
                    t['error']+=weight@independent.energy(observed-target);t['reference']+=weight@independent.energy(target)
                    for f in [1,2]:
                        N,g=m.shape(p,f,True);v=N@u[m.SLICES[f]];du=np.einsum('na,pnb->pab',u[m.SLICES[f]],g)
                        t['h1']+=aw@(np.sum(v*v,axis=1)+np.einsum('pai,pij,paj->p',du,np.linalg.inv(G),du))
            checks=dict(kernel=float(np.sqrt(kernel2/matrix2)/np.linalg.norm(saved['kernel']))<=1e-9,gauge=gauge_error<=1e-9);observations={}
            for kind,t in totals.items():
                e=float(np.sqrt(t['error']/t['reference']));h1=float(np.sqrt(t['h1']/area));de=abs(e-r['solutions'][kind]['error']);dh=abs(h1/r['solutions'][kind]['h1']-1)
                observations[kind]=dict(error=e,error_difference=de,h1_difference=dh);checks[kind]=de<=1e-8 and dh<=1e-8
            records.append(dict(broken_fields=r['broken_fields'],method=r['method'],solution_status='reproduced' if totals else 'not_solved_in_primary',
                                observations=observations,checks={k:bool(v) for k,v in checks.items()}))
            if not all(checks.values()):raise RuntimeError('independent_reproduction_failed')
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),parent_path,Path(independent.__file__),Path(primary.__file__)]
    dest.write_text(json.dumps(dict(records=records,all_available_fields_reproduced=len(records)==len(parent['records']) and stop is None,
        all_primary_solutions_available=all(bool(r['solutions']) for r in parent['records']),stop_reason=stop,
        independent_global_resolve=False,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
    print('audited',len(records),'stop',stop)

if __name__=='__main__':main()
