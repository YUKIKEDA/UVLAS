"""QR stress fit, surface work and separate matrix checks of saved GSP states."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import qr,solve_triangular,svd
import gsp_q3_v2 as candidate
import q3_thin_expansion as geometry

ROOT=rd_paths.ROOT


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(A,B):return float(np.linalg.norm(A)/max(1.,np.linalg.norm(B)))


def main():
    dest=ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-local-v2-audit.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();path=ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-local-v2-result.json';parent=json.loads(path.read_text())
    for f,h in parent['sources'].items():assert digest(ROOT/f)==h
    records=[];edges=[];preparation=None;state_checks=[];stop=None
    try:
        if parent['preparation_failure']:
            failure=parent['preparation_failure'];assert digest(ROOT/failure['artifact'])==failure['artifact_sha256']
            witness=np.load(ROOT/failure['artifact']);s=svd(witness['normalized_generator'],compute_uv=False,lapack_driver='gesvd',check_finite=False)
            r=candidate.ranks(s);expected=failure['report']['stress_ranks_svd']
            preparation=dict(svd_driver='gesvd',ranks=r,primary_ranks=expected,reproduced=r==expected,
                unambiguous=r==[r[1]]*3)
        for item in parent['states']:
            assert digest(ROOT/item['artifact'])==item['artifact_sha256']
            saved=np.load(ROOT/item['artifact']);state={k:saved[k] for k in ['q','frame_scale','length','center','T','active','column_norms','coeff','denominator']}
            field_reports=[]
            for order in [9,11]:
                p,w=candidate.BASE.quadrature(order);weight,_=candidate.measure(state,0.,p,w)
                S=candidate.stress_generator(state,p)
                A=(np.einsum('ij,pjn->pin',candidate.CIH,S)*np.sqrt(weight)[:,None,None]).reshape(-1,263)
                active=state['active'];s=svd(A[:,active]/state['column_norms'][active],compute_uv=False,lapack_driver='gesvd',check_finite=False)
                field_reports.append(dict(order=order,ranks=candidate.ranks(s)))
            expected=item['basis_gates'];reproduced=all(v['ranks']==old['ranks'] for v,old in zip(field_reports,expected))
            state_checks.append(dict(kappa=item['kappa'],n=item['n'],cell=item['cell'],primary_passed=item['passed'],field_reports=field_reports,reproduced=reproduced))
        p9,w9=candidate.BASE.quadrature(9);p11,w11=candidate.BASE.quadrature(11)
        N11,_=candidate.BASE.geometry.shape(p11)
        for item in parent['records']:
            assert digest(ROOT/item['state_artifact'])==item['state_sha256']
            assert digest(ROOT/item['artifact'])==item['artifact_sha256']
            m=np.load(ROOT/item['state_artifact']);saved=np.load(ROOT/item['artifact'])
            state={k:m[k] for k in ['q','frame_scale','length','center','T','active','column_norms','coeff','denominator']};h=item['h']
            weight,volume=candidate.measure(state,h,p9,w9)
            S=candidate.stress_field(state,p9);_,E=candidate.linear_raw(state['q'],h,p9)
            A=(np.einsum('ij,pjn->pin',candidate.CIH,S)*np.sqrt(weight)[:,None,None]).reshape(-1,len(state['T'][0]))
            Y=(np.einsum('ij,jk,pkn->pin',candidate.CIH,candidate.C,E)*np.sqrt(weight)[:,None,None]).reshape(-1,150)
            Q,R=qr(A,mode='economic',check_finite=False);a=solve_triangular(R,Q.T@Y,check_finite=False)
            S11=candidate.stress_field(state,p11)
            field=np.einsum('ij,pjn->pin',candidate.CI,S11@a)
            primary=np.einsum('ij,pjn->pin',candidate.CI,S11@saved['a'])
            weight11,volume11=candidate.measure(state,h,p11,w11)
            _,raw11=candidate.linear_raw(state['q'],h,p11)
            metrics=dict(field_difference=rel(field-primary,primary),
                constant=float(np.max(abs(field[:,:,144:]-np.eye(6)))),
                work=rel(np.einsum('p,pkn->kn',weight11,field-raw11),np.einsum('p,pkn->kn',weight11,raw11)))
            boundary=candidate.BASE.geometry.boundary(state['q'],h)
            actual_volume=volume11*np.einsum('p,pkn->nk',weight11,field[:,:,:144])
            metrics['surface_work']=rel(actual_volume-boundary,boundary)
            q=state['q'].reshape(16,3,3);z=h*p11[:,2]/2
            X=N11@q[:,0]+z[:,None]*(N11@q[:,1])+z[:,None]**2/2*(N11@q[:,2])
            expected=np.zeros((len(p11),3,3,4))
            for col,(component,variable) in enumerate([(0,0),(1,1),(0,1),(1,0)]):
                expected[:,component,variable,col]+=X[:,variable]/2
                expected[:,variable,component,col]+=X[:,variable]/2
            observed=field[:,:,:144]@saved['patch_vectors']
            analytic=np.einsum('kij,pijc->pkc',candidate.MB,expected)
            metrics['patch']=rel(observed-analytic,analytic)
            energies=[];exact=[]
            for col in range(4):
                tensor=np.einsum('kij,pk->pij',candidate.MB,observed[:,:,col])
                energy=lambda e:candidate.BASE.LAM*np.trace(e,axis1=1,axis2=2)**2+2*candidate.BASE.MU*np.einsum('pij,pij->p',e,e)
                energies.append(weight11@energy(tensor));exact.append(weight11@energy(expected[:,:,:,col]))
            metrics['patch_energy']=float(np.max(abs(np.array(energies)/exact-1)))
            # Reconstruct local nonlinear force and both tangent contributions via QR coefficients.
            Ecur,Bcur,Fcur,grad,J=candidate.green(state,saved['nonlinear_current'],h,p9)
            Ycur=(np.einsum('ij,jk,pk->pi',candidate.CIH,candidate.C,Ecur)*np.sqrt(weight)[:,None]).reshape(-1)
            acur=solve_triangular(R,Q.T@Ycur,check_finite=False)
            sigma=S@acur
            v1,v2=saved['direction1'],saved['direction2']
            target_delta1=Bcur@v1;target_delta2=Bcur@v2
            fit1=Q.T@(np.einsum('ij,jk,pk->pi',candidate.CIH,candidate.C,target_delta1)*np.sqrt(weight)[:,None]).reshape(-1)
            fit2=Q.T@(np.einsum('ij,jk,pk->pi',candidate.CIH,candidate.C,target_delta2)*np.sqrt(weight)[:,None]).reshape(-1)
            material=fit1@fit2
            dF1=np.einsum('na,pnb->pab',v1.reshape(48,3),grad);dF2=np.einsum('na,pnb->pab',v2.reshape(48,3),grad)
            stress_tensor=np.einsum('kij,pk->pij',candidate.MB,sigma)
            geometric=np.einsum('p,pij,pai,paj->',weight,stress_tensor,dF1,dF2)
            metrics['material_tangent']=float(abs(material-saved['material_tangent'])/max(1.,abs(saved['material_tangent'])))
            metrics['geometric_tangent']=float(abs(geometric-saved['geometric_tangent'])/max(1.,abs(saved['geometric_tangent'])))
            checks={k:bool(v<=1e-9) for k,v in metrics.items()}
            records.append(dict(kappa=item['kappa'],n=item['n'],cell=item['cell'],h=h,primary_passed=all(item['checks'].values()),metrics=metrics,checks=checks))
            print(item['kappa'],item['n'],item['cell'],h,checks,flush=True)
            if not all(checks.values()):raise RuntimeError('QR_field_or_tensor_audit_failed')
        # Constant stress interface forces must cancel for matching C0 traces.
        if parent['all_checks_pass']:
            face_nodes=[(0,-1.,[0,4,8,12]),(0,1.,[3,7,11,15]),(1,-1.,[0,1,2,3]),(1,1.,[12,13,14,15])]
            for k in [0.,.4]:
                xy,cells,coords,_=geometry.element.mesh(k,2,.2)
                for h in [.01,.001,.0001]:
                    for i,(_,_,left) in enumerate(cells):
                        for j,(_,_,right) in enumerate(cells):
                            if j<=i or len(set(left)&set(right))!=4:continue
                            shared=set(left)&set(right)
                            lf=next((axis,sign) for axis,sign,nodes in face_nodes if set(left[nodes])==shared)
                            rf=next((axis,sign) for axis,sign,nodes in face_nodes if set(right[nodes])==shared)
                            a=candidate.BASE.geometry.boundary(coords[left].reshape(144),h,[lf])
                            b=candidate.BASE.geometry.boundary(coords[right].reshape(144),h,[rf])
                            global_force=np.zeros((len(xy),9,6))
                            global_force[left]+=a.reshape(16,9,6);global_force[right]+=b.reshape(16,9,6)
                            error=rel(global_force,np.concatenate((a,b)))
                            edges.append(dict(kappa=k,h=h,cells=[i,j],error=error,passed=error<=1e-9))
                            if not edges[-1]['passed']:raise RuntimeError('shared_constant_stress_face_work_failed')
    except Exception as exc:stop=repr(exc)
    reproduced=(len(records)==len(parent['records']) and all(all(r['checks'].values()) for r in records)
        and all(s['reproduced'] for s in state_checks) and (preparation is None or preparation['reproduced']) and stop is None)
    sources=[Path(__file__),Path(candidate.__file__),path,ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-local-v2-protocol.md']
    dest.write_text(json.dumps(dict(records=records,states=state_checks,shared_edge_checks=edges,preparation_failure_audit=preparation,
        all_reproduction_checks=reproduced,primary_complete=parent['all_checks_pass'],physical_adoption=False,
        stop_reason=stop,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')


if __name__=='__main__':main()
