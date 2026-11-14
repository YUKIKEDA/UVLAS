"""QR-nullspace audit of M-Q3I1, with tensor patches and surface work."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from scipy.linalg import qr, solve, solve_triangular, svd, svdvals
import mq3i1_moments_v2 as candidate

ROOT=rd_paths.ROOT


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def relative(A,B):
    return float(np.linalg.norm(A)/max(1.,np.linalg.norm(B)))


def main():
    dest=ROOT/'rd/structure/main-shell/mq3/records/mq3i1-additional-local-audit.json'
    if dest.exists():
        raise FileExistsError(dest)
    start=time.monotonic()
    parent_path=ROOT/'rd/structure/main-shell/mq3/records/mq3i1-additional-local-result.json'
    parent=json.loads(parent_path.read_text())
    for path,sha in parent['sources'].items():
        assert digest(ROOT/path)==sha
    records,stop,preparation=[],None,None
    try:
        if parent['preparation_failure']:
            item=parent['preparation_failure']
            assert digest(ROOT/item['artifact'])==item['artifact_sha256']
            m=np.load(ROOT/item['artifact'])
            A,C=m['A'],m['C']
            _,s,V=svd(C,full_matrices=True,check_finite=False)
            rc=candidate.ranks(s)
            Z=V[rc[1]:].T
            _,sa,Va=svd(A@Z,full_matrices=True,check_finite=False)
            ra=candidate.ranks(sa)
            report=item['report']
            preparation=dict(C_ranks=rc,A_kernel_ranks=ra,kernel_dimension=Z.shape[1],
                agrees=rc==report['C_ranks'] and ra==report['A_kernel_ranks'])
            if ra[1]<Z.shape[1]:
                v=Z@Va[ra[1]]
                preparation.update(null_C=relative(C@v,C),null_A=relative(A@v,A),
                    target_coefficient_norm=float(np.linalg.norm(m['T']@v)))
                path=ROOT/item['artifact']
                witness=path.with_name('preparation-null-witness.npz')
                np.savez_compressed(witness,coefficient_direction=v,raw_target_direction=m['T']@v)
                preparation.update(witness=witness.relative_to(ROOT).as_posix(),witness_sha256=digest(witness))
            preparation['confirmed']=preparation['agrees'] and preparation.get('null_C',1.)<=1e-9 and preparation.get('null_A',1.)<=1e-9 and preparation.get('target_coefficient_norm',0.)>1e-12
        points,w=candidate.family.base.quadrature(11)
        N,_=candidate.family.base.geometry.shape(points)
        for record in parent['records']:
            assert digest(ROOT/record['state_artifact'])==record['state_sha256']
            assert digest(ROOT/record['artifact'])==record['artifact_sha256']
            m=np.load(ROOT/record['state_artifact'])
            state={key:m[key] for key in m.files}
            saved=np.load(ROOT/record['artifact'])
            h=record['h']
            C,d,_,_,_,_=candidate.gram(state,h)
            C=state['hard_U']@C;d=state['hard_U']@d
            Q,R=qr(C.T,mode='full',check_finite=False)
            n=len(C)
            particular=Q[:,:n]@solve_triangular(R[:n,:n].T,d,lower=True,check_finite=False)
            Z=Q[:,n:]
            AZ=state['A']@Z
            Qa,Ra=qr(AZ,mode='economic',check_finite=False)
            z=solve_triangular(Ra,Qa.T@(state['b']-state['A']@particular),check_finite=False)
            coefficients=particular+Z@z
            r=state['T'].shape[1]
            X=np.vstack((coefficients,np.zeros((len(C),150))))
            field=candidate.evaluate(state,h,X,points)
            primary=candidate.evaluate(state,h,saved['X'],points)
            raw=candidate.family.data(state['q'],h,points)
            weight=w*np.linalg.det(raw['J'])
            actual_weight=weight*h/2
            weight/=weight.sum()
            metrics=dict(field_difference=relative(field-primary,primary),constraints=relative(C@coefficients-d,d))
            zero=candidate.family.data(state['q'],0.,points)
            meanw=w*np.linalg.det(zero['J']);meanw/=meanw.sum()
            FF,gg=candidate.family.moments(state['q'],0.,zero,meanw)
            metrics['Qs_full_constraints']=relative(FF@state['T']@state['Qs'],FF@state['T'])
            metrics['Qs_membrane']=relative(zero['M']@state['T']@state['Qs'],zero['M']@state['T'])
            boundary=candidate.family.base.geometry.boundary(state['q'],h)
            work=np.einsum('p,pkn->nk',actual_weight,field[:,:,:144])
            metrics['surface_work']=relative(work-boundary,boundary)
            q=state['q'].reshape(16,3,3)
            zphys=h*points[:,2]/2
            xyz=N@q[:,0]+zphys[:,None]*(N@q[:,1])+zphys[:,None]**2/2*(N@q[:,2])
            E=np.zeros((len(points),3,3,4))
            for col,(component,variable) in enumerate([(0,0),(1,1),(0,1),(1,0)]):
                E[:,component,variable,col]+=xyz[:,variable]/2
                E[:,variable,component,col]+=xyz[:,variable]/2
            expected=np.einsum('kij,pijc->pkc',candidate.family.base.MB,E)
            observed=field[:,:,:144]@saved['patch_vectors']
            metrics['patch_error']=relative(observed-expected,expected)
            metrics['constant']=float(np.max(abs(field[:,:,144:]-np.eye(6))))
            eo=[];ee=[]
            for col in range(4):
                tensor=np.einsum('kij,pk->pij',candidate.family.base.MB,observed[:,:,col])
                target=E[:,:,:,col]
                value=lambda e:candidate.family.base.LAM*np.trace(e,axis1=1,axis2=2)**2+2*candidate.family.base.MU*np.einsum('pij,pij->p',e,e)
                eo.append(weight@value(tensor));ee.append(weight@value(target))
            metrics['patch_energy']=float(np.max(abs(np.array(eo)/ee-1)))
            S,b,_,_=candidate.system(state,h)
            direct=solve(S,b,assume_a='gen',check_finite=False)
            metrics['LU_vs_QR']=relative(direct[:r]-coefficients,coefficients)
            checks={key:bool(value<=1e-9) for key,value in metrics.items()}
            records.append(dict(kappa=record['kappa'],n=record['n'],cell=record['cell'],h=h,
                primary_passed=all(record['checks'].values()),metrics=metrics,checks=checks))
            print(record['kappa'],record['n'],record['cell'],h,checks,flush=True)
            if not all(checks.values()):
                raise RuntimeError('independent_local_audit_failed')
    except Exception as exc:
        stop=repr(exc)
    sources=[Path(__file__),Path(candidate.__file__),parent_path,
        ROOT/'rd/structure/main-shell/mq3/records/mq3i1-local-execution-protocol.md',ROOT/'rd/structure/main-shell/mq3/records/mq3i1-v2-pure-shear-protocol.md']
    result=dict(records=records,preparation_failure_audit=preparation,
        primary_complete=parent['all_checks_pass'],primary_conditions=len(parent['records']),
        all_reproduction_checks=(len(records)==len(parent['records']) and stop is None and
            (preparation is None or preparation['confirmed'])),
        physical_adoption=False,stop_reason=stop,independent_basis_implementation=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start)
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8')


if __name__=='__main__':
    main()
