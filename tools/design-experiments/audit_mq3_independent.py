"""Saved M-Q3 20-case audit with monomial/QR reconstruction and shared patches."""
from pathlib import Path
import hashlib,json,time,itertools
import numpy as np
from scipy.linalg import svdvals
import mq3_independent as audit
import mq3_joint_moments as primary
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rank(s):return [int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]]
def edge_check(state):
    p,w=np.polynomial.legendre.leggauss(11);F=[];b=[]
    for tangent,sign in itertools.product(range(2),[-1.,1.]):
        pts=np.zeros((11,3));pts[:,tangent]=p;pts[:,1-tangent]=sign;d=audit.data(state['q'],state['h'],pts)
        for i in range(3):F.append((w*p**i)@d['M'][:,tangent]);b.append((w*p**i)@d['m'][:,tangent])
        side=np.zeros((121,3));side[:,tangent]=np.repeat(p,11);side[:,1-tangent]=sign;side[:,2]=np.tile(p,11);d=audit.data(state['q'],state['h'],side)
        for i,j in itertools.product(range(3),range(3)):
            weights=np.outer(w*p**i,w*p**j).reshape(-1);F.append(weights@d['S'][:,tangent]);b.append(weights@d['g'][:,tangent])
    F=np.array(F);b=np.array(b);return float(np.linalg.norm(F@state['coeff']-b)/max(1.,np.linalg.norm(b)))
def covfield(state,points,u,modified=True):
    d=audit.data(state['q'],state['h'],points);E=audit.evaluate(state,points) if modified else d['E'];value=E[:,:,:144]@u.reshape(144);return audit.covariant(d['J'],value[:,:,None])[:,:,:,0]
def shared(states,k,h):
    side=7;coords=np.full((49,3,3),np.nan);cells={};volume=np.zeros((49,9,6));outside=np.zeros_like(volume);points,w=audit.quadrature(11);constant=0.
    for j,i in itertools.product(range(2),range(2)):
        cell=2*j+i;state=states[k,h,cell];ids=np.array([(3*j+b)*side+3*i+a for b,a in itertools.product(range(4),range(4))]);local=state['q'].reshape(16,3,3)
        for node,q in zip(ids,local):
            if np.all(np.isfinite(coords[node])):assert np.max(abs(coords[node]-q))<=1e-14
            coords[node]=q
        cells[i,j]=(ids,state);d=audit.data(state['q'],h,points);E=audit.evaluate(state,points);measure=w*np.linalg.det(d['J'])*h/2;volume[ids]+=np.einsum('p,pki->ik',measure,E[:,:,:144]).reshape(16,9,6);affine=np.einsum('kij,nj->nik',audit.MB,state['q'].reshape(48,3)).reshape(144,6);constant=max(constant,float(np.max(abs(E[:,:,:144]@affine-np.eye(6)))));faces=[(2,-1),(2,1)]
        for axis,index in [(0,i),(1,j)]:faces.append((axis,-1 if index==0 else 1))
        outside[ids]+=audit.geometry.boundary(state['q'],h,faces).reshape(16,9,6)
    u=np.random.default_rng(20261103).standard_normal((49,9));p,pw=np.polynomial.legendre.leggauss(11);worst_raw=worst_mod=0.;edge_count=0
    for axis in range(2):
        for j,i in itertools.product(range(2),range(2)):
            if (i if axis==0 else j)!=0:continue
            neighbor=(i+1,j) if axis==0 else (i,j+1);tangent=1-axis;pair=[cells[i,j],cells[neighbor]];moments={False:[],True:[]}
            for sign,(ids,state) in zip([1.,-1.],pair):
                middle=np.zeros((11,3));middle[:,axis]=sign;middle[:,tangent]=p;face=np.zeros((121,3));face[:,axis]=sign;face[:,tangent]=np.repeat(p,11);face[:,2]=np.tile(p,11)
                for modified in [False,True]:
                    cm=covfield(state,middle,u[ids],modified);cs=covfield(state,face,u[ids],modified);values=[float((pw*p**degree)@cm[:,tangent,tangent]) for degree in range(3)]
                    for a,b in itertools.product(range(3),range(3)):values.append(float(np.outer(pw*p**a,pw*p**b).reshape(-1)@cs[:,tangent,2]))
                    moments[modified].append(np.array(values))
            for modified in [False,True]:
                left,right=moments[modified];error=float(np.linalg.norm(left-right)/max(1.,np.linalg.norm(left),np.linalg.norm(right)))
                if modified:worst_mod=max(worst_mod,error)
                else:worst_raw=max(worst_raw,error)
            edge_count+=1
    work=float(np.linalg.norm(volume-outside)/np.linalg.norm(outside));checks=dict(raw_moments=worst_raw<=1e-9,modified_moments=worst_mod<=1e-9,global_work=work<=1e-9,constant_patch=constant<=1e-10)
    return dict(kappa=k,h=h,internal_edges=edge_count,raw_moment_difference=worst_raw,modified_moment_difference=worst_mod,global_work_difference=work,constant_patch_error=constant,checks=checks)
def main():
    target=ROOT/'docs/verification/mq3-independent-audit-result.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();parent=ROOT/'docs/verification/mq3-local-feasibility-result.json';prior=json.loads(parent.read_text());assert prior['all_checks_pass'];records=[];shared_records=[];states={};stop=None;out=ROOT/'build/analysis'/('mq3-independent-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();points,w=audit.quadrature(11)
    try:
        for path,sha in prior['sources'].items():assert digest(ROOT/path)==sha
        for record in prior['records']:
            artifact=ROOT/record['artifact'];assert digest(artifact)==record['artifact_sha256'];saved=np.load(artifact);old={key:saved[key] for key in saved.files};state=audit.prepare(old['q'],float(old['h']));d=audit.data(state['q'],state['h'],points);new=audit.evaluate(state,points);oldfield=primary.evaluate(old,points);weight=w*np.linalg.det(d['J'])*state['h']/2;W=np.sqrt(weight)[:,None,None];raw=np.einsum('ij,pjn->pin',audit.CH,d['E']);change=np.einsum('ij,pjn->pin',audit.CH,new-oldfield);difference=float(np.max(np.linalg.norm((change*W).reshape(-1,150),axis=0)/np.linalg.norm((raw*W).reshape(-1,150),axis=0)));constant=float(np.max(abs(new[:,:,144:]-np.eye(6))));work=np.einsum('p,pki->ik',weight,new[:,:,:144]);work_raw=np.einsum('p,pki->ik',weight,d['E'][:,:,:144]);workdiff=float(np.linalg.norm(work-work_raw)/np.linalg.norm(work_raw));boundary=audit.geometry.boundary(state['q'],state['h']);boundarydiff=float(np.linalg.norm(work-boundary)/np.linalg.norm(boundary));rawB=(raw[:,:,:144]*W).reshape(-1,144);mod=(np.einsum('ij,pjn->pin',audit.CH,new[:,:,:144])*W).reshape(-1,144);s=svdvals(mod/np.linalg.norm(rawB,axis=0),check_finite=False);branks=rank(s);edge=edge_check(state);checks=dict(basis_rank=state['basis_ranks']==record['space_ranks'],constraint_rank=state['constraint_ranks']==record['constraint_ranks'],constraint_residual=state['constraint_residual']<=1e-10,zero_rhs=state['zero_rhs']<=1e-10,field=difference<=1e-8,constant=constant<=1e-10,work=workdiff<=1e-10,boundary=boundarydiff<=1e-9,moments=edge<=1e-10,local_rank=branks==[138]*3,resources=time.monotonic()-start<=180 and memory.memory_bytes()<=2*1024**3);file=out/f'k{record["kappa"]:g}-n{record["n"]}-h{record["h"]:g}-c{record["cell"]}.npz';np.savez_compressed(file,**state);records.append(dict(kappa=record['kappa'],n=record['n'],h=record['h'],cell=record['cell'],basis_ranks=state['basis_ranks'],constraint_ranks=state['constraint_ranks'],B_ranks=branks,field_difference=difference,constant_error=constant,work_error=workdiff,boundary_error=boundarydiff,moment_error=edge,checks={key:bool(v) for key,v in checks.items()},artifact=file.relative_to(ROOT).as_posix(),artifact_sha256=digest(file)));print(record['kappa'],record['n'],record['h'],record['cell'],'difference',difference,'checks',checks,flush=True)
            if record['n']==2:states[record['kappa'],record['h'],record['cell']]=state
            if not all(checks.values()):stop='required_local_audit_gate_failed';break
        if stop is None:
            for k,h in itertools.product([0.,.4],[.01,.0001]):
                log=shared(states,k,h);shared_records.append(log);print('shared',log,flush=True)
                if not all(log['checks'].values()):stop='required_shared_gate_failed';break
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(audit.__file__),Path(audit.geometry.__file__),Path(primary.__file__),parent,ROOT/'docs/verification/mq3-independent-audit-protocol.md']
    target.write_text(json.dumps(dict(records=records,shared=shared_records,all_checks_pass=len(records)==20 and len(shared_records)==4 and stop is None,stop_reason=stop,independent_coefficient_construction=True,independent_mesh_generation=False,global_equilibrium_solved=False,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
