"""Bounded diagnostic of the saved failed basis; never changes local gates."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.linalg import qr,solve_triangular,svdvals
import gsp_q3 as g

ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1.,np.linalg.norm(b)))

def main():
    dest=ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-basis-diagnostic.json'
    if dest.exists():raise FileExistsError(dest)
    primary=ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-local-result.json'
    parent=json.loads(primary.read_text());item=parent['states'][-1]
    assert not item['passed']
    artifact=ROOT/item['artifact'];assert digest(artifact)==item['artifact_sha256']
    data=np.load(artifact);state={k:data[k] for k in ['q','frame_scale','length','center','T','active','column_norms']}
    T=state['T'];selected=np.flatnonzero(np.linalg.norm(T,axis=1)>0)
    # Keep precisely the recorded independent subspace, no rank truncation.
    p,w=g.BASE.quadrature(9);weight,_=g.measure(state,0.,p,w)
    raw=g.stress_generator(state,p)
    def weighted(S,mu):return (np.einsum('ij,pjn->pin',g.CIH,S)*np.sqrt(mu)[:,None,None]).reshape(-1,T.shape[1])
    norms=state['column_norms'][selected]
    A=weighted(raw[:,:,selected]/norms,weight)
    Q,R=qr(A,mode='economic',check_finite=False)
    singular=svdvals(A,check_finite=False)
    records=[]
    for order in [9,11,13]:
        p,w=g.BASE.quadrature(order);mu,_=g.measure(state,0.,p,w)
        raw=g.stress_generator(state,p)
        inverse_field=raw@T
        stable_field=solve_triangular(R.T,(raw[:,:,selected]/norms).reshape(-1,len(selected)).T,
            lower=True,check_finite=False).T.reshape(len(p),6,-1)
        inverse_A=weighted(inverse_field,mu);stable_A=weighted(stable_field,mu)
        records.append(dict(order=order,
            recorded_inverse_gram=rel(inverse_A.T@inverse_A-np.eye(len(selected)),np.eye(len(selected))),
            triangular_gram=rel(stable_A.T@stable_A-np.eye(len(selected)),np.eye(len(selected))),
            triangular_q_difference=rel(stable_A-Q,Q) if order==9 else None,
            direct_q_gram=rel(Q.T@Q-np.eye(len(selected)),np.eye(len(selected))) if order==9 else None))
    result=dict(scope='failed_saved_state_only_no_local_rerun',rank=len(selected),
        selected_normalized_condition=float(singular[0]/singular[-1]),
        selected_sigma_min=float(singular[-1]),selected_sigma_max=float(singular[0]),
        qr_backward_residual=rel(A-Q@R,A),records=records,
        longdouble_mantissa=int(np.finfo(np.longdouble).nmant),
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in [Path(__file__),Path(g.__file__),primary,artifact]},
        changes_space=False,changes_threshold=False,local_gate_pass=False)
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
