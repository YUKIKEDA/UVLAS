"""Independent central-mode gradients/Gram and boundary traction audit."""
from pathlib import Path
import hashlib,json,time
import numpy as np
import q9_curved_constrained_gradient as method
import audit_q9_curved_strain_screen as independent
import check_q9_curved_load as geometry
import check_q9_plate_convergence as memory
from q9_edge_moment_fields import quadrature
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def relative(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-30))
def main():
    target=F/'q9-curved-constrained-gradient-audit.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-curved-constrained-gradient-result.json';r=json.loads(parent.read_text());assert r['all_criteria_pass']
    assert r['source_sha256']==digest(ROOT/'tools/design-experiments/check_q9_curved_constrained_gradient.py')
    assert r['contract_sha256']==digest(ROOT/'docs/q9-curved-constrained-gradient-candidate.md')
    assert r['protocol_sha256']==digest(F/'q9-curved-constrained-gradient-experiment.md')
    for name,sha in r['dependencies'].items():assert digest(ROOT/'tools/design-experiments'/name)==sha
    for c in r['cases']:assert digest(ROOT/c['artifact'])==c['artifact_sha256']
    start=time.monotonic();records=[];stop=None
    try:
        for c in r['cases']:
            if c['rotated'] or c['chart']:continue
            name,h=c['geometry'],c['h'];q=geometry.geometry(name);state=method.prepare(q,h,0.,9);saved=np.load(ROOT/c['artifact']);points,w=quadrature(9);_,J,grad=independent.data(q,h,points);weight=w*np.linalg.det(J);V=weight.sum();psi=np.einsum('kaj,pj->pka',method.MB,grad[:,24,:]);mean=np.einsum('p,pki->ki',weight,psi)/V;b0=(1-points[:,0]**2)*(1-points[:,1]**2);bubble=b0/(weight@b0/V);phi=psi-bubble[:,None,None]*mean
            G=np.einsum('p,pki,pkj->ij',weight,phi,phi,optimize=True);_,E,R=method.base.evaluate(state,points);rhs=np.einsum('p,pki,pkn->in',weight,phi,E-R,optimize=True);alpha=np.linalg.solve(G,rhs)
            psi_error=relative(psi,E[:,:,72:75]);mean_error=float(np.linalg.norm(mean-saved['mode_mean'])/max(np.sqrt(np.einsum('p,pki,pki->',weight,psi,psi)/V),1e-30));gram=relative(G,saved['gram']);coefficient=relative(alpha,saved['mode_alpha'])
            ip,iw=quadrature(11);Ji,_,Bi=method.evaluate(state,ip);force=np.einsum('p,pkn->nk',iw*np.linalg.det(Ji),Bi[:,:,:81]);expected=independent.boundary(q,h)@np.diag([1,1,1,1/np.sqrt(2),1/np.sqrt(2),1/np.sqrt(2)]);boundary=relative(force,expected)
            checks=dict(psi=psi_error<=1e-10,mean=mean_error<=1e-10,gram=gram<=1e-10,coefficients=coefficient<=1e-10,boundary_work=boundary<=1e-10);records.append(dict(geometry=name,h=h,psi=psi_error,mean=mean_error,gram=gram,coefficients=coefficient,boundary_work=boundary,checks={k:bool(v) for k,v in checks.items()}));print(json.dumps(records[-1]),flush=True)
            if time.monotonic()-start>300 or memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
    except Exception as exc:stop=repr(exc)
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(parent),protocol_sha256=digest(F/'q9-curved-constrained-gradient-audit-protocol.md'),dependencies={Path(m.__file__).name:digest(m.__file__) for m in [method,method.base,independent,geometry,memory]},cases=records,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=bool(stop is None and len(records)==10 and all(all(c['checks'].values()) for c in records)),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['cases','dependencies']}))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
