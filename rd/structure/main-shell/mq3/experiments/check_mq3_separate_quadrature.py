"""Unchanged saved coefficients observed with separate quadrature nodes."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time,itertools
import numpy as np
from scipy.linalg import svdvals
import mq3_joint_moments as method
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def edge_residual(q,h,coeff):
    p,w=np.polynomial.legendre.leggauss(11);L=np.polynomial.legendre.legvander(p,2);errors=[];right=[]
    for tangent,sign in itertools.product([0,1],[-1.,1.]):
        points=np.zeros((11,3));points[:,tangent]=p;points[:,1-tangent]=sign;d=method.data(q,h,points)
        for i in range(3):
            rhs=np.einsum('p,pn->n',w*L[:,i],d['m'][:,tangent]);value=np.einsum('p,pn->n',w*L[:,i],d['M'][:,tangent])@coeff;errors.append(value-rhs);right.append(rhs)
        face=np.zeros((121,3));face[:,tangent]=np.repeat(p,11);face[:,1-tangent]=sign;face[:,2]=np.tile(p,11);d=method.data(q,h,face)
        for i,j in itertools.product(range(3),repeat=2):
            wf=np.outer(w*L[:,i],w*L[:,j]).reshape(-1);rhs=np.einsum('p,pn->n',wf,d['g'][:,tangent]);value=np.einsum('p,pn->n',wf,d['S'][:,tangent])@coeff;errors.append(value-rhs);right.append(rhs)
    return float(np.linalg.norm(errors)/max(1.,np.linalg.norm(right)))
def main():
    dest=ROOT/'rd/structure/main-shell/mq3/records/mq3-separate-quadrature-result.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();parent=ROOT/'rd/structure/main-shell/mq3/records/mq3-local-feasibility-result.json';data=json.loads(parent.read_text());assert data['all_checks_pass'];records=[];points,w=method.element.quadrature(11)
    for p,sha in data['sources'].items():assert digest(ROOT/p)==sha
    for record in data['records']:
        artifact=ROOT/record['artifact'];assert digest(artifact)==record['artifact_sha256'];saved=np.load(artifact);state={key:saved[key] for key in saved.files};q,h=state['q'],float(state['h']);d=method.data(q,h,points);R=method.evaluate(state,points);weight=w*np.linalg.det(d['J']);weight/=weight.sum();constant=float(np.max(abs(R[:,:,144:]-np.eye(6))));raw_work=np.einsum('p,pkn->kn',weight,d['E'][:,:,:144]);work=np.einsum('p,pkn->kn',weight,R[:,:,:144]);difference=float(np.linalg.norm(work-raw_work)/np.linalg.norm(raw_work));raw=(np.einsum('kl,pln->pkn',method.element.CH,d['E'][:,:,:144])*np.sqrt(weight)[:,None,None]).reshape(-1,144);mod=(np.einsum('kl,pln->pkn',method.element.CH,R[:,:,:144])*np.sqrt(weight)[:,None,None]).reshape(-1,144);s=svdvals(mod/np.linalg.norm(raw,axis=0),check_finite=False);ranks=method.spectrum_ranks(s);edge=edge_residual(q,h,state['coeff']);checks=dict(constant=constant<=1e-10,work=difference<=1e-10,moments=edge<=1e-10,rank=ranks==[138]*3);records.append(dict(kappa=record['kappa'],n=record['n'],cell=record['cell'],h=h,constant_error=constant,work_error=difference,moment_error=edge,ranks=ranks,checks={key:bool(v) for key,v in checks.items()}))
    sources=[Path(__file__),Path(method.__file__),parent,ROOT/'rd/structure/main-shell/mq3/records/mq3-separate-quadrature-protocol.md']
    result=dict(records=records,all_checks_pass=len(records)==20 and all(all(r['checks'].values()) for r in records),independent_implementation_audit=False,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start);dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print({key:max(r[key] for r in records) for key in ['constant_error','work_error','moment_error']});print('all_checks_pass',result['all_checks_pass'])
if __name__=='__main__':main()
