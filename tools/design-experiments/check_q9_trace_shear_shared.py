"""Shared traces and admissible w-space, before T-S36 global loads."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import qr
import q9_trace_shear as candidate
import check_q9_distorted_plate as mesh_source
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    target=F/'q9-trace-shear-shared-result.json'
    if target.exists():raise FileExistsError(target)
    gate=F/'q9-trace-shear-local-result.json';assert json.loads(gate.read_text())['all_criteria_pass'];out=ROOT/'build/analysis'/('q9-trace-shear-shared-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();start=time.monotonic();records=[];stop=None
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
                W=np.zeros((12*n*n,len(wm)));D=np.zeros((12*n*n,len(dm)));local_range=True;metric=0.;ties=[]
                for index,(_,_,ids,_) in enumerate(cells):
                    d=candidate.matrices(mesh_source.configuration(nodes[ids]),.1,9);C=d['C0'];ties.append(d['T0']);metric=max(metric,d['leading_metric_difference']);cw=C[:,2::9]/n;cd=C[:,np.array([[9*k+3,9*k+4] for k in range(9)]).reshape(-1)]
                    sc=np.linalg.svd(np.column_stack((cd,cw)),compute_uv=False);sd=np.linalg.svd(cd,compute_uv=False)
                    local_range=local_range and all(np.count_nonzero(sd>t*sc[0])==np.count_nonzero(sc>t*sc[0])==12 for t in [1e-9,1e-10,1e-11])
                    rows=slice(index*12,(index+1)*12)
                    for k,node in enumerate(ids):
                        node=int(node)
                        if node in wm:W[rows,wm[node]]+=cw[:,k]
                        for component in [3,4]:
                            if (node,component) in dm:D[rows,dm[node,component]]+=C[:,9*k+component]
                U,sd,_=np.linalg.svd(D,full_matrices=True);sc=np.linalg.svd(np.column_stack((D,W)),compute_uv=False);scale=max(sd[0],sc[0]);ranks={}
                for t in [1e-9,1e-10,1e-11]:
                    rd=int(np.count_nonzero(sd>t*scale));rc=int(np.count_nonzero(sc>t*scale));sz=np.linalg.svd(U[:,rd:].T@W,compute_uv=False);rz=int(np.count_nonzero(sz>t*scale))
                    ranks[str(t)]=dict(direction_rank=rd,combined_rank=rc,compatible_w_dimension=len(wm)-(rc-rd),nullspace_w_dimension=len(wm)-rz)
                qranks=[]
                for A in [D,np.column_stack((D,W))]:
                    _,R,_=qr(A,pivoting=True,mode='economic');qranks.append(int(np.count_nonzero(abs(np.diag(R))>1e-10*scale)))
                def embed(index,rows):
                    result=np.zeros((2,9*len(nodes)));result[:,cells[index][3]]=ties[index][rows];return result
                trace=0.;edges=0
                for j in range(n):
                    for i in range(n):
                        index=j*n+i
                        for other,a,b in ([(index+n,[2,5],[0,3])] if j+1<n else [])+([(index+1,[10,11],[6,7])] if i+1<n else []):
                            A=embed(index,a);B=embed(other,b);trace=max(trace,float(np.linalg.norm(A-B)/max(1.,np.linalg.norm(A),np.linalg.norm(B))));edges+=1
                q=mesh_source.configuration(nodes[cells[0][2]]);points=np.array([[-candidate.A,-1.,0.],[candidate.A,-1.,0.]])
                forward=candidate.sampled(q,.1,points,np.array([0,0]));reverse=candidate.sampled(q,.1,points[::-1],np.array([0,0]),[-1,-1]);wrong=candidate.sampled(q,.1,points[::-1],np.array([0,0]))
                orientation=float(np.linalg.norm(-reverse[::-1]-forward)/max(1.,np.linalg.norm(forward)));negative=float(np.linalg.norm(-wrong[::-1]-forward)/max(1.,np.linalg.norm(forward)))
                central=ranks['1e-10'];checks=dict(local_range=local_range,metric=metric<=1e-12,rank_stable=len({r['direction_rank'] for r in ranks.values()})==1 and len({r['combined_rank'] for r in ranks.values()})==1,
                    all_w_retained=all(r['compatible_w_dimension']==r['nullspace_w_dimension']==len(wm) for r in ranks.values()),qr=qranks==[central['direction_rank'],central['combined_rank']],trace=trace<=1e-12,orientation=orientation<=1e-12,negative=negative>1e-3)
                checks={k:bool(v) for k,v in checks.items()};file=out/f'a{amplitude:g}-n{n}.npz';np.savez_compressed(file,D=D,W=W,nodes=nodes,vertices=vertices,w_nodes=np.array(list(wm)),d_node_components=np.array(list(dm)),ties=np.array(ties))
                records.append(dict(n=n,amplitude=amplitude,w_dofs=len(wm),direction_dofs=len(dm),ranks=ranks,qr_ranks=qranks,internal_edges=edges,trace_difference=trace,orientation_difference=orientation,negative_difference=negative,
                    metric_difference=metric,checks=checks,artifact=file.relative_to(ROOT).as_posix(),artifact_sha256=digest(file)))
                print(json.dumps(dict(n=n,amplitude=amplitude,compatible_w_dimension=central['compatible_w_dimension'],checks=checks)),flush=True)
                if time.monotonic()-start>120 or memory.memory_bytes()>512*1024**2:raise RuntimeError('resource_limit')
        assert sum(p.stat().st_size for p in out.iterdir())<=16*1024**2
    except Exception as exc:stop=repr(exc)
    result=dict(source_sha256=digest(__file__),candidate_sha256=digest(ROOT/'docs/q9-trace-shear-candidate.md'),protocol_sha256=digest(F/'q9-trace-shear-gate-experiment.md'),local_gate_sha256=digest(gate),
        dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [candidate,mesh_source,memory]},cases=records,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=bool(stop is None and len(records)==8 and all(all(r['checks'].values()) for r in records)),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n')
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
