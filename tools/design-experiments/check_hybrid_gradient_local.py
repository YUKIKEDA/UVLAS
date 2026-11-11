"""Local necessary gates for the pre-registered mixed weak-gradient ANCF."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import svd
import hybrid_gradient_ancf as f
import check_q9_plate_convergence as memory
m=f.m;ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1e-30,np.linalg.norm(b)))

def generic():
    p,_=m.native_plane();x=(p[:,0]+1)/2;y=p[:,1]/2
    return np.vstack((np.column_stack((x,y,.2*x*x+.08*x**4*y**4)),np.column_stack((-.4*x-.32*x**3*y**4,-.32*x**4*y**3,np.ones(25))),np.zeros((25,3)))).reshape(-1)

def evaluate(q,name,h,kind,out):
    state=f.prepare(q,h,kind);root=state['root'];norm=np.linalg.norm(root,axis=0);singular=svd(root/norm,compute_uv=False,check_finite=False)
    X0=state['X0'];rigid=[];rindices=np.r_[np.arange(25),*[np.arange(75+15*i,80+15*i) for i in range(4)]]
    for axis in np.eye(3):
        translation=np.zeros((135,3));translation[rindices]=axis;rigid.extend([translation.reshape(-1),np.cross(axis,X0).reshape(-1)])
    rigid=np.array(rigid).T
    inside=np.einsum('p,pjn->nj',state['weight'],state['Gw']);outside=f.boundary_map(state)
    rng=np.random.default_rng(20261104);current=X0+.0001*rng.normal(size=X0.shape);v1=rng.normal(size=405);v2=rng.normal(size=405)
    U,force,tangent,info=f.nonlinear(state,current,v1,v2)
    Uc,_,_,_=f.nonlinear(state,current.reshape(-1)+1e-25j*v1,v1,v2);_,fc,_,_=f.nonlinear(state,current.reshape(-1)+1e-25j*v2,v1,v2)
    axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis);A=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
    Q=np.eye(3)+np.sin(.71)*A+(1-np.cos(.71))*A@A;rotated=current@Q.T;rotated[rindices]+=np.array([.2,-.1,.3])
    Ur,fr,_,_=f.nonlinear(state,rotated,v1,v2)
    fluxforce=np.zeros((135,3))
    for face,traction in zip(state['faces'],info['fluxes']):fluxforce[face['offset']:face['offset']+15]=np.einsum('p,pn,pa->na',face['measure'],face['hat'],traction)
    metrics=dict(reference_gradient=state['reference_gradient_error'],constant_P_work=rel(inside-outside,outside),rigid=rel(root@rigid,root)/np.linalg.norm(rigid),
        energy_force=float(abs(Uc.imag/1e-25-force@v1)/max(1e-30,abs(force@v1))),force_tangent=float(abs(v1@(fc.imag/1e-25)-tangent)/max(1e-30,abs(tangent))),
        rotation_energy=rel(Ur-U,U),rotation_force=rel(fr-(force.reshape(135,3)@Q.T).reshape(-1),force),
        flux_dual=rel(force.reshape(135,3)[75:]-fluxforce[75:],force.reshape(135,3)[75:]))
    checks=dict(basis=state['ranks']==[state['S'].shape[2]]*3 and max(state['gram_error'],state['span_error'])<=1e-10,
        rank=f.ranks(singular)==[399]*3,identities=all(v<=1e-9 for v in metrics.values()),orientation=state['minimum_J']>0 and info['min_F']>0)
    path=out/f'{kind}-{name}-{h:g}.npz';np.savez_compressed(path,q=q,current=current,v1=v1,v2=v2,energy=U,force=force,tangent=tangent,coeff=state['coeff'],D=state['D'],singular=singular)
    return dict(geometry=name,h=h,kind=kind,ranks=f.ranks(singular),basis_ranks=state['ranks'],gram=state['gram_error'],span=state['span_error'],metrics=metrics,
        min_J=state['minimum_J'],min_F=info['min_F'],checks={k:bool(v) for k,v in checks.items()},artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path))

def main():
    dest=ROOT/'docs/verification/hybrid-gradient-local-result.json';assert not dest.exists();start=time.monotonic();records=[];stops=[]
    out=ROOT/'build/analysis'/('hybrid-gradient-local-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    geometries=[('flat',m.mesh(0.,1,0.)['ref'][m.mesh(0.,1,0.)['cells'][0]].reshape(-1)),('curved-distorted',m.mesh(.4,2,.2)['ref'][m.mesh(.4,2,.2)['cells'][0]].reshape(-1)),('quartic-reference',generic())]
    for kind in ['S4','S3+']:
        try:
            for name,q in geometries:
                for h in [.01,.001,.0001]:
                    tick=time.monotonic();r=evaluate(q,name,h,kind,out);r['checks']['resources']=time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3
                    records.append(r);print(kind,name,h,r['ranks'],r['checks'],r['metrics'],flush=True)
                    if not r['checks']['resources']:raise RuntimeError('local_resource_limit')
        except Exception as exc:stops.append(dict(kind=kind,reason=repr(exc)));print('stop',kind,repr(exc),flush=True)
    sources=[Path(__file__),Path(f.__file__),Path(m.__file__),ROOT/'docs/hybrid-gradient-ancf-formulation.md',ROOT/'docs/verification/hybrid-gradient-ancf-protocol.md',ROOT/'docs/verification/hybrid-gradient-local-preparation.json']
    dest.write_text(json.dumps(dict(records=records,stops=stops,physical_adoption=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
