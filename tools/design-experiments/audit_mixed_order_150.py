"""Vandermonde / complex physical-z / tensor / fieldwise H1 route."""
from pathlib import Path
import hashlib,json,itertools,time
import numpy as np
from scipy.linalg import svd
import mixed_order_150 as m
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ranks(s):return [int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]]
def rel(a,b):return float(np.linalg.norm(a)/max(1.,np.linalg.norm(b)))
def energy(E):return .3/(1.3*.4)*np.trace(E,axis1=1,axis2=2)**2+2/2.6*np.einsum('pij,pij->p',E,E)

def main():
    path=ROOT/'docs/verification/mixed-order-150-result.json';parent=json.loads(path.read_text());dest=ROOT/'docs/verification/mixed-order-150-audit.json';assert not dest.exists()
    for f,h in parent['sources'].items():assert digest(ROOT/f)==h
    start=time.monotonic();locals=[];records=[];stop=None;p,w=m.quadrature(11,2);Nr,gr=m.shape(p,0,True)
    try:
        for item in parent['local_records']:
            assert digest(ROOT/item['artifact'])==item['artifact_sha256'];saved=np.load(ROOT/item['artifact']);q=saved['q'];finite=[]
            for index,old in enumerate(item['finite']):
                h=old['h'];pts,ws=m.quadrature(11);_,J,B=m.raw(q,pts,h*pts[:,2]/2,True)
                W=(np.einsum('ij,pjn->pin',m.CH,B)*np.sqrt(ws*np.linalg.det(J))[:,None,None]).reshape(-1,150)
                s=svd(W*saved[f'scale{index}'],compute_uv=False,lapack_driver='gesvd',check_finite=False)
                volume=h/2*np.einsum('p,pki->ik',ws*np.linalg.det(J),B);surface=m.boundary(q,h)
                metrics=dict(rank_reproduced=ranks(s)==old['ranks11'],surface_work=rel(volume-surface,surface))
                finite.append(dict(h=h,ranks=ranks(s),metrics=metrics,passed=bool(metrics['rank_reproduced'] and metrics['surface_work']<=1e-9)))
            locals.append(dict(kappa=item['kappa'],amplitude=item['amplitude'],cell=item['cell'],finite=finite,passed=all(f['passed'] for f in finite)))
            if not locals[-1]['passed']:raise RuntimeError('mixed_order_local_audit_failed')
        for item in parent['records']:
            assert digest(ROOT/item['artifact'])==item['artifact_sha256'];s=np.load(ROOT/item['artifact']);free=s['free_rd'];lookup={int(v):i for i,v in enumerate(free)}
            totals={kind:dict(error=0.,ref=0.,h1=0.,zero=0.) for kind in item['solutions']};area=kernel2=matrix2=0.
            for ids in s['cells']:
                q=s['ref'][ids].reshape(150);_,J,B0=m.raw(q,p,0.,True);_,_,Bc=m.raw(q,p,1e-25j,True);B1=Bc.imag/1e-25
                weight=w*np.linalg.det(J);K=np.zeros((150,s['kernel'].shape[1]))
                for local,node in enumerate(ids[:41]):
                    if int(node) in lookup:K[3*local:3*local+3]=s['kernel'][3*lookup[int(node)]:3*lookup[int(node)]+3]
                WB=np.einsum('ij,pjn->pin',m.CH,B0[:,:,:123]);WK=np.einsum('ij,pjn->pin',m.CH,B0@K)
                matrix2+=np.sum(weight[:,None,None]*WB**2);kernel2+=np.sum(weight[:,None,None]*WK**2)
                r=s['ref'][ids[:25]];Js=np.einsum('na,pnb->pab',r,gr);G=Js.transpose(0,2,1)@Js;aw=w*np.sqrt(np.linalg.det(G));area+=aw.sum();x=(Nr@r)[:,0]
                k=item['kappa'];a=np.column_stack((np.ones_like(x),np.zeros_like(x),k*x));D=1+k*k*x*x
                for kind in totals:
                    u=s['u_'+kind][ids];v={'quadratic':np.ones_like(x),'cubic':x,'quartic':x*x}[kind]
                    target=-v[:,None,None]*a[:,:,None]*a[:,None,:]/D[:,None,None]
                    E=np.einsum('kij,pk->pij',m.MB,B1@u.reshape(150))
                    totals[kind]['error']+=weight@energy(E-target);totals[kind]['ref']+=weight@energy(target)
                    totals[kind]['zero']+=np.sum(weight[:,None]*np.einsum('ij,pj->pi',m.CH,B0@u.reshape(150))**2)
                    for f in [1,2]:
                        N,g=m.shape(p,f,True);value=N@u[m.SLICES[f]];du=np.einsum('na,pnb->pab',u[m.SLICES[f]],g)
                        totals[kind]['h1']+=aw@(np.sum(value*value,axis=1)+np.einsum('pai,pij,paj->p',du,np.linalg.inv(G),du))
            kernel=float(np.sqrt(kernel2/matrix2)/np.linalg.norm(s['kernel']));checks=dict(kernel=kernel<=1e-9);observations={}
            for kind,t in totals.items():
                error=float(np.sqrt(t['error']/t['ref']));h1=float(np.sqrt(t['h1']/area));zero=float(np.sqrt(t['zero']/matrix2)/np.linalg.norm(s['u_'+kind][free]))
                difference=abs(error-item['solutions'][kind]['observed_error']);hdiff=abs(h1/item['solutions'][kind]['observed_h1']-1)
                checks[kind]=difference<=1e-8 and hdiff<=1e-8 and zero<=1e-9
                observations[kind]=dict(error=error,error_difference=difference,h1=h1,h1_difference=hdiff,zero=zero)
            records.append(dict(kappa=item['kappa'],amplitude=item['amplitude'],kernel_residual=kernel,observations=observations,checks={k:bool(v) for k,v in checks.items()}))
            print('audit',item['kappa'],item['amplitude'],checks,flush=True)
            if not all(checks.values()):raise RuntimeError('mixed_order_shared_audit_failed')
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(m.__file__),path,ROOT/'docs/verification/mixed-order-150-protocol.md']
    dest.write_text(json.dumps(dict(local_records=locals,records=records,
        all_reproduction_checks=bool(locals) and len(locals)==len(parent['local_records']) and len(records)==len(parent['records']) and stop is None,
        primary_complete=parent['all_numerical_checks'],independent_global_resolve=False,physical_adoption=False,stop_reason=stop,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
if __name__=='__main__':main()
