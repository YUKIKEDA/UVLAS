"""Leading shear constraints: local existence versus shared-node admissible space."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import solve_triangular,qr
import check_q9_distorted_plate as mesh_source
import check_q9_shear_projection_gate as projection
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification';THRESHOLDS=[1e-9,1e-10,1e-11]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=F/'q9-shear-gradient-compatibility-result.json'
    if output.exists():raise FileExistsError(output)
    out=ROOT/'build/analysis'/('q9-shear-gradient-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();start=time.monotonic();records=[];stop=None
    try:
        for amplitude in [0.,.2]:
            for n in [1,2,4,8]:
                vertices,nodes,cells=mesh_source.mesh(n,amplitude);side=2*n+1;wm={};dm={}
                for j in range(side):
                    for i in range(side):
                        node=j*side+i
                        if 0<i<side-1 and 0<j<side-1:wm[node]=len(wm)
                        if 0<j<side-1:dm[node,3]=len(dm)
                        if 0<i<side-1:dm[node,4]=len(dm)
                W=np.zeros((12*n*n,len(wm)));D=np.zeros((12*n*n,len(dm)));local=[];max_orthogonality=0.
                for index,(_,_,ids,_) in enumerate(cells):
                    data=projection.matrices(mesh_source.configuration(nodes[ids]),.1,9);p0=np.arange(0,36,3);other=np.setdiff1d(np.arange(36),p0)
                    M=data['M'][np.ix_(p0,p0)]/.1;H=data['H'][p0,:]/.1;L=np.linalg.cholesky(M);C=solve_triangular(L,H,lower=True)
                    max_orthogonality=max(max_orthogonality,float(np.linalg.norm(data['M'][np.ix_(p0,other)])/np.linalg.norm(data['M'])))
                    cw=C[:,2::9]/n;cd=C[:,np.array([[9*k+3,9*k+4] for k in range(9)]).reshape(-1)]
                    sc=np.linalg.svd(np.column_stack((cd,cw)),compute_uv=False);sd=np.linalg.svd(cd,compute_uv=False);common=max(sc[0],sd[0])
                    local.append({str(t):dict(direction_rank=int(np.count_nonzero(sd>t*common)),combined_rank=int(np.count_nonzero(sc>t*common))) for t in THRESHOLDS})
                    rows=slice(index*12,(index+1)*12)
                    for k,node in enumerate(ids):
                        node=int(node)
                        if node in wm:W[rows,wm[node]]+=cw[:,k]
                        for component in [3,4]:
                            if (node,component) in dm:D[rows,dm[node,component]]+=C[:,9*k+component]
                combined=np.column_stack((D,W));U,sd,_=np.linalg.svd(D,full_matrices=True);sc=np.linalg.svd(combined,compute_uv=False);common=max(sd[0],sc[0]);ranks={}
                for t in THRESHOLDS:
                    rd=int(np.count_nonzero(sd>t*common));rc=int(np.count_nonzero(sc>t*common));Z=U[:,rd:].T@W;sz=np.linalg.svd(Z,compute_uv=False);rz=int(np.count_nonzero(sz>t*common))
                    ranks[str(t)]=dict(direction_rank=rd,combined_rank=rc,defect_rank=rz,compatible_w_dimension=len(wm)-(rc-rd),nullspace_w_dimension=len(wm)-rz)
                def qr_rank(A):
                    _,r,_=qr(A,mode='economic',pivoting=True);return int(np.count_nonzero(abs(np.diag(r))>1e-10*common))
                qrd=qr_rank(D);qrc=qr_rank(combined);central=ranks['1e-10'];dims=[r['compatible_w_dimension'] for r in ranks.values()]
                checks=dict(thickness_orthogonal=max_orthogonality<=1e-12,local_range=all(all(r['direction_rank']==r['combined_rank']==12 for r in c.values()) for c in local),
                    rank_stable=len(set(dims))==1 and len({r['direction_rank'] for r in ranks.values()})==1 and len({r['combined_rank'] for r in ranks.values()})==1,
                    nullspace_agrees=all(r['compatible_w_dimension']==r['nullspace_w_dimension'] for r in ranks.values()),
                    qr_agrees=qrd==central['direction_rank'] and qrc==central['combined_rank'],rectangular_nontrivial=amplitude!=0 or central['compatible_w_dimension']>0)
                file=out/f'a{amplitude:g}-n{n}.npz';np.savez_compressed(file,D=D,W=W,nodes=nodes,vertices=vertices,w_nodes=np.array(list(wm)),d_node_components=np.array(list(dm)),direction_singular_values=sd,combined_singular_values=sc)
                records.append(dict(n=n,amplitude=amplitude,w_dofs=len(wm),direction_dofs=len(dm),constraint_rows=len(D),ranks=ranks,qr_ranks=dict(direction=qrd,combined=qrc),
                    local_rank_counts={str(t):dict(direction_min=min(c[str(t)]['direction_rank'] for c in local),combined_max=max(c[str(t)]['combined_rank'] for c in local)) for t in THRESHOLDS},
                    maximum_thickness_cross_term=max_orthogonality,checks=checks,artifact=file.relative_to(ROOT).as_posix(),artifact_sha256=digest(file)))
                print(json.dumps(dict(n=n,amplitude=amplitude,compatible_w_dimension=dims,local_range=checks['local_range'],checks=checks)),flush=True)
                if time.monotonic()-start>120 or memory.memory_bytes()>512*1024**2:raise RuntimeError('diagnostic_resource_limit')
        a=next(r for r in records if r['n']==1 and r['amplitude']==0);b=next(r for r in records if r['n']==1 and r['amplitude']==.2)
        A=np.load(ROOT/a['artifact']);B=np.load(ROOT/b['artifact']);same=bool(np.array_equal(A['nodes'],B['nodes']) and np.array_equal(A['D'],B['D']) and np.array_equal(A['W'],B['W']));assert same
        assert sum(p.stat().st_size for p in out.iterdir())<=16*1024**2
    except Exception as exc:stop=repr(exc);same=False
    result=dict(source_sha256=digest(__file__),protocol_sha256=digest(ROOT/'docs/q9-shear-gradient-compatibility-contract.md'),candidate_sha256=digest(ROOT/'docs/q9-shear-projection-candidate.md'),
        dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [mesh_source,projection,memory]},cases=records,n1_same_geometry_and_operator=same,
        seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),stop_reason=stop,all_criteria_pass=bool(stop is None and len(records)==8 and all(all(r['checks'].values()) for r in records)),product_shell_adopted=False)
    output.write_text(json.dumps(result,indent=2)+'\n')
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
