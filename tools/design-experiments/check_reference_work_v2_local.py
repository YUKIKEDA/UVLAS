"""Frozen RW-1 local/work necessary gates, including all six face contributions."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import svd
import reference_work_ancf_v2 as f
import check_polynomial_trace_local as prior
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1.,np.linalg.norm(b)))

def work(state):
    ids=[v for v in f.PID if sum(v)<=5]+[(5,5,k) for k in range(3)];phi,g=f.basis(state['p'],ids);u=f.reduce(state,lambda p:f.basis(p,ids)[0])
    truth=np.einsum('p,pjr,pjn->rn',state['mu'],state['S'],g,optimize=True);volume=state['V']@u[:75]-phi
    ev=np.einsum('p,pr,pn->rn',state['mu'],state['div'],volume,optimize=True);parts=[]
    for face in state['faces']:
        val,_=f.basis(face['points'],ids);err=face['hat']@u[face['offset']:face['offset']+face['n']]-val
        parts.append(np.einsum('p,pr,pn->rn',face['w']/8,face['sign']*face['S'][:,face['axis'],:],err,optimize=True))
    gw=state['D']@u;gp=state['Gwref']@u;moment=np.einsum('p,pjr,pjn->rn',state['mu'],state['S'],gp,optimize=True)
    return dict(fields=len(ids),volume_error_work=rel(ev,truth),vertical_error_work=max(rel(x,truth) for x in parts[2:]),top_bottom_error_work=max(rel(x,truth) for x in parts[:2]),
        weak_commutation=rel(gw-truth,truth),body_work=rel(moment-truth,truth),decomposition=rel(gw-truth+ev-sum(parts),truth))

def patch(state):
    ids=[v for v in f.PID if sum(v)<=5];_,g=f.basis(state['p'],ids);u=f.reduce(state,lambda p:f.basis(p,ids)[0]);true=np.einsum('pjn,pjk->pkn',g,state['invJ'])
    gradient=rel(state['Gw']@u-true,true);jump=max(rel(face['jump']@u,u) for face in state['faces'])
    return dict(fields=len(ids),gradient_error=gradient,stabilization_residual=jump,passed=bool(max(gradient,jump)<=1e-9))

def evaluate(q,name,h,out):
    s=f.prepare(q,h);ns=s['nscalar'];nc=3*ns;root=s['root'];sing=svd(root/np.linalg.norm(root,axis=0),compute_uv=False,check_finite=False)
    X0=s['X0'];ridx=np.array([0]+[face['offset'] for face in s['faces']]);rigid=[]
    for axis in np.eye(3):
        v=np.zeros_like(X0);v[ridx]=axis;rigid.extend([v.reshape(-1),np.cross(axis,X0).reshape(-1)])
    rigid=np.array(rigid).T;inside=np.einsum('p,pjn->nj',s['weight'],s['Gw']);outside=f.boundary_map(s)
    rng=np.random.default_rng(20261104);Fstar=np.array([[1.003,.002,0.],[0.,1.003,.001],[0.,0.,1.003]])
    current=X0@Fstar.T+h*1e-4*rng.normal(size=X0.shape);v1=rng.normal(size=nc);v2=rng.normal(size=nc)
    U,force,K,info=f.nonlinear(s,current,v1,v2);Uc,_,_,_=f.nonlinear(s,current.reshape(-1)+1e-25j*v1,v1,v2);_,fc,_,_=f.nonlinear(s,current.reshape(-1)+1e-25j*v2,v1,v2)
    axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis);A=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
    Q=np.eye(3)+np.sin(.71)*A+(1-np.cos(.71))*A@A;rot=current@Q.T;rot[ridx]+=[.2,-.1,.3];Ur,fr,_,_=f.nonlinear(s,rot,v1,v2)
    ff=np.zeros_like(current)
    for face,tr in zip(s['faces'],info['fluxes']):ff[face['offset']:face['offset']+face['n']]=np.einsum('p,pn,pa->na',face['measure'],face['hat'],tr)
    metrics=dict(reference_gradient=s['reference_gradient_error'],constant_P_work=rel(inside-outside,outside),rigid=rel(root@rigid,root)/np.linalg.norm(rigid),
        energy_force=float(abs(Uc.imag/1e-25-force@v1)/max(1e-30,abs(force@v1))),force_tangent=float(abs(v1@(fc.imag/1e-25)-K)/max(1e-30,abs(K))),
        rotation_energy=float(abs(Ur-U)/abs(U)),rotation_force=float(np.linalg.norm(fr-(force.reshape(ns,3)@Q.T).reshape(-1))/np.linalg.norm(force)),
        flux_dual=float(np.linalg.norm(force.reshape(ns,3)[75:]-ff[75:])/np.linalg.norm(force.reshape(ns,3)[75:])))
    wk=work(s);poly=patch(s) if name=='flat' else None
    checks=dict(rank=f.ranks(sing)==[nc-6]*3,identities=max(metrics.values())<=1e-9,work=max(v for k,v in wk.items() if k!='fields')<=1e-9,
        polynomial=poly is None or poly['passed'],potential_rank=s['potential_ranks']==[104]*3,gram=s['gram_error']<=1e-10,orientation=s['minimum_J']>0 and info['min_F']>0)
    path=out/f'{name}-{h:g}.npz';np.savez_compressed(path,q=q,current=current,v1=v1,v2=v2,energy=U,force=force,tangent=K,singular=sing,Gw=s['Gw'],D=s['D'],H=s['H'],jumps=np.vstack([x['jump'] for x in s['faces']]))
    return dict(geometry=name,h=h,nscalar=ns,ni=s['ni'],mixed_dimension=s['S'].shape[-1],mixed_completion_ranks=s['mixed_completion']['ranks'],face_dimensions=[x['n'] for x in s['faces']],
        face_completion_ranks=[x['completion']['ranks'] for x in s['faces']],ranks=f.ranks(sing),gram_error=s['gram_error'],metrics=metrics,work=wk,polynomial=poly,checks={k:bool(v) for k,v in checks.items()},artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path))

def main():
    dest=ROOT/'docs/verification/reference-work-v2-local-result.json';assert not dest.exists();start=time.monotonic();records=[];stops=[]
    out=ROOT/'build/analysis'/('reference-work-v2-local-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    geoms=[('flat',f.m.mesh(0.,1,0.)),('curved-distorted',f.m.mesh(.4,2,.2))]
    qs=[(name,grid['ref'][grid['cells'][0]].reshape(-1)) for name,grid in geoms]+[('quartic-reference',prior.generic())]
    for name,q in qs:
        for h in [.01,.001,.0001]:
            try:
                tick=time.monotonic();r=evaluate(q,name,h,out);r['checks']['resources']=time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3;records.append(r)
                print(name,h,r['nscalar'],r['ranks'],r['checks'],r['metrics'],r['work'],flush=True)
                (out/'checkpoint.json').write_text(json.dumps(dict(records=records,stops=stops),indent=2),encoding='utf-8')
            except Exception as exc:stops.append(dict(geometry=name,h=h,reason=repr(exc)));print('STOP',name,h,repr(exc),flush=True)
    sources=[Path(__file__),Path(f.__file__),ROOT/'docs/reference-work-ancf-formulation.md',ROOT/'docs/verification/reference-work-v2-protocol.md',ROOT/'docs/verification/reference-work-v2-local-preparation.json']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_local_checks=len(records)==9 and not stops and all(all(r['checks'].values()) for r in records),physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
