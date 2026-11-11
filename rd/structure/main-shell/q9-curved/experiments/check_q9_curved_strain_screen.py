"""Fixed G-COV/C-COV constant-field and dual-work screen; no load solve."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
from numpy.polynomial.legendre import leggauss
import q9_trace_shear as trace
import check_q9_curved_load as geometry
import check_q9_independent_load as common
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def tensor(v):
    e=np.zeros((len(v),3,3,v.shape[2]));e[:,0,0]=v[:,0];e[:,1,1]=v[:,1];e[:,2,2]=v[:,2]
    for row,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):e[:,i,j]=e[:,j,i]=v[:,row]/2
    return e
def voigt(e):return np.stack((e[:,0,0],e[:,1,1],e[:,2,2],2*e[:,1,2],2*e[:,0,2],2*e[:,0,1]),axis=1)
EB=tensor(np.eye(6)[None])[0].transpose(2,0,1)
SB=EB.copy();SB[3:]*=2
def tables(points,width,offset):
    scalar=[];gradient=[]
    for s,t,zeta in points:
        N,ds,dt=common.shape(s,t);z=offset+width*zeta/2
        scalar.append(np.column_stack((N,z*N,z*z*N/2)).reshape(27))
        gradient.append(np.column_stack((np.column_stack((ds,z*ds,z*z*ds/2)).reshape(27),np.column_stack((dt,z*dt,z*z*dt/2)).reshape(27),np.column_stack((np.zeros(9),N,z*N)).reshape(27)*width/2)))
    return np.array(scalar),np.array(gradient)
def kinematics(q,width,offset,points):
    N,g=tables(points,width,offset);J=np.einsum('ia,pib->pab',q.reshape(27,3),g);det=np.linalg.det(J)
    if np.min(det)<=0:raise ValueError('nonpositive reference Jacobian')
    grad=np.linalg.solve(J.transpose(0,2,1),g.transpose(0,2,1)).transpose(0,2,1);B=np.zeros((len(N),6,81))
    B[:,0,0::3]=grad[:,:,0];B[:,1,1::3]=grad[:,:,1];B[:,2,2::3]=grad[:,:,2];B[:,3,1::3]=grad[:,:,2];B[:,3,2::3]=grad[:,:,1];B[:,4,0::3]=grad[:,:,2];B[:,4,2::3]=grad[:,:,0];B[:,5,0::3]=grad[:,:,1];B[:,5,1::3]=grad[:,:,0]
    return J,B
def covariant(A,E):return np.einsum('pia,pijk,pjb->pabk',A,E,A,optimize=True)
def reduce(E,Es,A,As,S,components):
    cov=covariant(A,E);sample=covariant(As,Es);beta=np.array([sample[i,c,2] for i,c in enumerate(components)])
    replacement=np.einsum('pai,ik->pak',S,beta,optimize=True)
    cov[:,0,2]=cov[:,2,0]=replacement[:,0];cov[:,1,2]=cov[:,2,1]=replacement[:,1];inv=np.linalg.inv(A)
    return np.einsum('pai,pabk,pbj->pijk',inv,cov,inv,optimize=True)
def column_error(A,B):return float(np.max(np.linalg.norm(A-B,axis=0)/np.maximum(np.linalg.norm(B,axis=0),np.finfo(float).tiny)))
def field_error(O):
    delta=O-np.eye(6)[None];weights=np.array([1.,1.,1.,.5,.5,.5])
    return float(np.max(np.sqrt(np.sum(delta*delta*weights[None,:,None],axis=1)/weights[None,:])))
def rotation_maps(R):
    e=np.einsum('ia,kab,jb->kij',R,EB,R);s=np.einsum('ia,kab,jb->kij',R,SB,R)
    return voigt(e.transpose(1,2,0)[None])[0],np.array([[v[0,0],v[1,1],v[2,2],v[1,2],v[0,2],v[0,1]] for v in s]).T
def pull_force(f,K,R):
    value=np.einsum('ab,nbck->nack',K.T,f.reshape(9,3,3,6))
    return np.einsum('ab,ntbk->ntak',R.T,value).reshape(81,6)
def evaluate(q,width,offset,order,scheme):
    p,w=leggauss(order);points=np.array(list(itertools.product(p,p,p)));weights=np.array([a*b*c for a,b,c in itertools.product(w,w,w)]);ties,components=trace.sampling_points()
    J,B=kinematics(q,width,offset,points);Js,Bs=kinematics(q,width,offset,ties);Jc,_=kinematics(q,width,offset,[[0.,0.,0.]])
    A=J if scheme=='G-COV' else np.broadcast_to(Jc,J.shape);As=Js if scheme=='G-COV' else np.broadcast_to(Jc,Js.shape);S=trace.interpolation(points)
    O=voigt(reduce(np.broadcast_to(EB.transpose(1,2,0),(len(J),3,3,6)),np.broadcast_to(EB.transpose(1,2,0),(36,3,3,6)),A,As,S,components))
    Bb=voigt(reduce(tensor(B),tensor(Bs),A,As,S,components));measure=weights*np.linalg.det(J);raw=np.einsum('p,pik->ki',measure,B);work=np.einsum('p,pik->ki',measure,Bb)
    U=np.einsum('bij,nj->nib',EB,q.reshape(27,3)).reshape(81,6);affine=float(np.max(abs(np.einsum('pik,kj->pij',B,U)-np.eye(6)[None])))
    construction=float(np.max(abs(np.einsum('pik,kj->pij',Bb,U)-O)))
    extra=np.array([[-.71,.22,-.31],[.18,-.64,.73],[.62,.43,-.27]]);Je,_=kinematics(q,width,offset,extra);Ae=Je if scheme=='G-COV' else np.broadcast_to(Jc,Je.shape)
    extraO=voigt(reduce(np.broadcast_to(EB.transpose(1,2,0),(len(extra),3,3,6)),np.broadcast_to(EB.transpose(1,2,0),(36,3,3,6)),Ae,As,trace.interpolation(extra),components))
    return dict(O=O,work=work,raw_work=raw,J=J,value_error=max(field_error(O),field_error(extraO)),work_error=column_error(work,raw),affine_error=affine,construction_error=construction)
def main():
    target=F/'q9-curved-strain-screen-result.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();out=ROOT/'build/analysis'/('q9-curved-strain-screen-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();records=[];baselines={};stop=None
    try:
        for name in ['RECT','SKEW','TRAP','SADDLE','NORMAL-C']:
            for h in [.1,.01]:
                q0=geometry.geometry(name);J0,_=kinematics(q0,h,0.,np.array(list(itertools.product(leggauss(9)[0],repeat=3))))
                for rotated,reparameterized in itertools.product([False,True],repeat=2):
                    R=geometry.ROTATION if rotated else np.eye(3);q=geometry.transform(q0,R);K=np.eye(3);width=h;offset=0.
                    if reparameterized:
                        a=1.3;b=.1*h;K=np.array([[1.,b,b*b/2],[0.,a,a*b],[0.,0.,a*a]]);q=np.einsum('ab,nbc->nac',K,q.reshape(9,3,3)).reshape(81);width=h/a;offset=-b/a
                    Emap,Smap=rotation_maps(R);Einv,_=rotation_maps(R.T)
                    for scheme in ['G-COV','C-COV']:
                        values={n:evaluate(q,width,offset,n,scheme) for n in [7,9]};v=values[9];key=(name,h,scheme)
                        if not rotated and not reparameterized:baselines[key]=v
                        base=baselines[key];backO=np.einsum('ai,pij,jb->pab',Einv,v['O'],Emap,optimize=True);backforce=pull_force(v['work']@Smap,K,R);backraw=pull_force(v['raw_work']@Smap,K,R)
                        objective=float(np.max(abs(backO-base['O'])));workcov=column_error(backforce,base['work']);rawcov=column_error(backraw,base['raw_work']);geometry_error=float(np.max(abs(np.einsum('ab,pbc->pac',R.T,v['J'])-J0)))
                        quadrature=column_error(values[7]['work'],v['work']);checks=dict(geometry=geometry_error<=1e-12,affine=v['affine_error']<=1e-11,construction=v['construction_error']<=1e-11,objectivity=objective<=1e-10,work_covariance=workcov<=1e-10,raw_work_covariance=rawcov<=1e-10)
                        file=out/f'{name}-h{h:g}-R{int(rotated)}-P{int(reparameterized)}-{scheme}.npz';np.savez_compressed(file,q=q,width=width,offset=offset,O=v['O'],work=v['work'],raw_work=v['raw_work'])
                        records.append(dict(geometry=name,h=h,rotated=rotated,reparameterized=reparameterized,scheme=scheme,value_error=v['value_error'],work_error=v['work_error'],work_quadrature=quadrature,
                            value_gate=v['value_error']<=1e-10,work_gate=(v['work_error']<=1e-10 if quadrature<=1e-8 else None),integration_gate=quadrature<=1e-8,coordinate_defects=dict(geometry=geometry_error,value=objective,work=workcov,raw_work=rawcov),checks=checks,artifact=file.relative_to(ROOT).as_posix(),artifact_sha256=digest(file)))
                        if time.monotonic()-start>300 or memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
                print(json.dumps(dict(geometry=name,h=h,schemes=[dict(scheme=r['scheme'],value_error=r['value_error'],work_error=r['work_error'],work_quadrature=r['work_quadrature']) for r in records[-8:] if not r['rotated'] and not r['reparameterized']])),flush=True)
        if sum(p.stat().st_size for p in out.iterdir())>32*1024**2:raise RuntimeError('disk_limit')
    except Exception as exc:stop=repr(exc)
    result=dict(source_sha256=digest(__file__),protocol_sha256=digest(F/'q9-curved-strain-screen-experiment.md'),dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [trace,geometry,common,memory]},cases=records,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),diagnostic_checks_pass=bool(stop is None and len(records)==80 and all(all(r['checks'].values()) for r in records)),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(diagnostic_checks_pass=result['diagnostic_checks_pass'],cases=len(records),seconds=result['seconds'],stop_reason=stop)))
    if not result['diagnostic_checks_pass']:raise SystemExit(1)
if __name__=='__main__':main()
