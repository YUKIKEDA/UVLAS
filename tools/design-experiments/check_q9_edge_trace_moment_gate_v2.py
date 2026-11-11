"""T-E24-M1 local stability before global-load adoption."""
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
import q9_edge_trace_moment as candidate
import q9_edge_trace24 as trace
import check_q9_edge_moment_shared_space_v2 as shared
import check_q9_curved_load as geometry
import check_q9_independent_load as material
from q9_edge_moment_fields import quadrature
from check_q9_trace_work_solvability import rank_metrics
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    target=F/'q9-edge-trace-moment-gate-v2-result.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-trace-work-solvability-result.json';assert json.loads(parent.read_text())['edge_work_all_w']
    start=time.monotonic();cases=[];spaces=[];stop=None;original=shared.method
    D=np.diag([1,1,1,np.sqrt(2),np.sqrt(2),np.sqrt(2)]);C=D@material.elastic@D
    out=ROOT/'build/analysis'/('q9-edge-trace-moment-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    try:
        for name,h in itertools.product(['RECT','SKEW','TRAP','SADDLE','NORMAL-C'],[.1,.01]):
            q=geometry.geometry(name);state=candidate.prepare(q,h,0.,5);points,w=quadrature(5);J,E,B=candidate.evaluate(state,points);measure=w*np.linalg.det(J);V=measure.sum();rawB=E[:,:,:81];red=B[:,:,:81]
            K=np.einsum('p,pik,ij,pjl->kl',measure,red,C,red,optimize=True);rawK=np.einsum('p,pik,ij,pjl->kl',measure,rawB,C,rawB,optimize=True);scale=1/np.sqrt(rawK.diagonal());ks=scale[:,None]*K*scale[None,:];s=np.linalg.svd(ks,compute_uv=False);ranks={str(t):int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]}
            value=float(np.max(np.linalg.norm(B[:,:,81:]-np.eye(6)[None],axis=1)));work=np.einsum('p,pkn->kn',measure,B-E)/V;work_error=float(np.linalg.norm(work)/np.linalg.norm(np.einsum('p,pkn->kn',measure,E)/V));sym=float(np.linalg.norm(ks-ks.T)/np.linalg.norm(ks));minimum=float(np.linalg.eigvalsh((ks+ks.T)/2)[0]);coords=q.reshape(9,3,3);modes=[]
            for axis in np.eye(3):
                u=np.zeros_like(coords);u[:,0,:]=axis;modes.extend([u.reshape(81),np.cross(np.broadcast_to(axis,coords.shape),coords).reshape(81)])
            rigid=max(float(np.linalg.norm(ks@(u/scale))/(s[0]*np.linalg.norm(u/scale))) for u in modes)
            checks=dict(value=value<=1e-10,work=work_error<=1e-10,rank=all(r==75 for r in ranks.values()),rigid=rigid<=1e-11,symmetry=sym<=1e-12,positive=minimum>=-1e-10*s[0])
            path=out/f'{name}-h{h:g}.npz';np.savez_compressed(path,K=K,rawK=rawK,singular_values=s)
            cases.append(dict(geometry=name,h=h,ranks=ranks,value=value,work=work_error,rigid=rigid,symmetry=sym,minimum_eigenvalue=minimum,checks={k:bool(v) for k,v in checks.items()},artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path)))
            print(json.dumps(dict(geometry=name,h=h,ranks=ranks,checks=cases[-1]['checks'])),flush=True)
        shared.method=candidate
        configs=[('TRAP',n,h) for h,n in itertools.product([.1,.01,.001],[2,4])]+[(name,n,.01) for name,n in itertools.product(['RECT','SKEW'],[2,4])]
        for name,n,h in configs:
            a=shared.assemble(name,n,h);r=rank_metrics(a['D'],a['W']);checks=all(v['compatible_w']==v['projected_w']==a['W'].shape[1] for v in r.values());spaces.append(dict(geometry=name,n=n,h=h,w_dofs=a['W'].shape[1],ranks=r,all_w_retained=bool(checks)))
            print(json.dumps(dict(geometry=name,n=n,h=h,all_w_retained=bool(checks))),flush=True)
            assert time.monotonic()-start<300 and shared.memory.memory_bytes()<1024**3
    except Exception as exc:stop=repr(exc)
    finally:shared.method=original
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(parent),protocol_sha256=digest(F/'q9-edge-trace-moment-gate-v2-experiment.md'),contract_sha256=digest(ROOT/'docs/q9-trace-work-solvability-contract.md'),dependencies={Path(m.__file__).name:digest(m.__file__) for m in [candidate,trace,shared,shared.raw,geometry,material]},local_cases=cases,shared_cases=spaces,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=shared.memory.memory_bytes(),all_criteria_pass=bool(stop is None and len(cases)==10 and len(spaces)==10 and all(all(c['checks'].values()) for c in cases) and all(c['all_w_retained'] for c in spaces)),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['local_cases','shared_cases','dependencies']}))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
