"""G-EM1 simultaneous consistency/rank gates; no equilibrium solve."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
from numpy.polynomial.legendre import leggauss
import q9_edge_moment as method
import q9_trace_shear as flat
import check_q9_curved_strain_screen as raw
import check_q9_curved_strain_boundary as shared
import check_q9_curved_load as geometry
import check_q9_independent_load as material
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;CM=np.diag([1.,1.,1.,np.sqrt(2.),np.sqrt(2.),np.sqrt(2.)])@material.elastic@np.diag([1.,1.,1.,np.sqrt(2.),np.sqrt(2.),np.sqrt(2.)])
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def quadrature(order):
    p,w=leggauss(order);return np.array(list(itertools.product(p,p,p))),np.array([a*b*c for a,b,c in itertools.product(w,w,w)])
def fieldnorm(O):return float(np.max(np.linalg.norm(O-np.eye(6)[None],axis=1)))
def metrics(state,order,stiffness=False):
    points,w=quadrature(order);J,E,B=method.evaluate(state,points);measure=w*np.linalg.det(J);work=np.einsum('p,pik->ki',measure,B[:,:,:81]);original=np.einsum('p,pik->ki',measure,E[:,:,:81]);result=dict(J=J,O=B[:,:,81:],B=B[:,:,:81],raw=E[:,:,:81],work=work,raw_work=original,work_error=raw.column_error(work,original),value_error=fieldnorm(B[:,:,81:]))
    if stiffness:
        result['K']=np.einsum('p,pik,ij,pjl->kl',measure,B[:,:,:81],CM,B[:,:,:81],optimize=True);result['rawK']=np.einsum('p,pik,ij,pjl->kl',measure,E[:,:,:81],CM,E[:,:,:81],optimize=True)
    return result
def rotation(R):
    moved=np.einsum('ia,kab,jb->kij',R,method.MB,R);return np.einsum('aij,bij->ab',method.MB,moved)
def main():
    target=F/'q9-edge-moment-feasibility-result.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();out=ROOT/'build/analysis'/('q9-edge-moment-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();cases=[];pairs=[];baselines={};stop=None
    try:
        for name,h in itertools.product(['RECT','SKEW','TRAP','SADDLE','NORMAL-C'],[.1,.01]):
            q0=geometry.geometry(name)
            for rotated,reparameterized in itertools.product([False,True],repeat=2):
                R=geometry.ROTATION if rotated else np.eye(3);K3=np.eye(3);q=geometry.transform(q0,R);width=h;offset=0.
                if reparameterized:
                    a=1.3;b=.1*h;K3=np.array([[1.,b,b*b/2],[0.,a,a*b],[0.,0.,a*a]]);q=np.einsum('ab,nbc->nac',K3,q.reshape(9,3,3)).reshape(81);width=h/a;offset=-b/a
                states={o:method.prepare(q,width,offset,o) for o in [7,9]};d=metrics(states[9],9,True);check11=metrics(states[9],11);_,_,oldfield=method.evaluate(states[7],quadrature(9)[0]);field_difference=float(np.linalg.norm(oldfield[:,:,:81]-d['B'])/np.linalg.norm(d['raw']))
                scale=1/np.sqrt(d['rawK'].diagonal());ks=scale[:,None]*d['K']*scale[None,:];sv=np.linalg.svd(ks,compute_uv=False);ranks={str(t):int(np.count_nonzero(sv>t*sv[0])) for t in [1e-11,1e-10,1e-9]};symmetric=float(np.linalg.norm(ks-ks.T)/np.linalg.norm(ks));minimum=float(np.linalg.eigvalsh((ks+ks.T)/2)[0]);norm=float(np.linalg.norm(ks,2))
                coords=q.reshape(9,3,3);rigid=[]
                for axis in np.eye(3):
                    u=np.zeros_like(coords);u[:,0,:]=axis;rigid.append(u.reshape(81));rigid.append(np.cross(np.broadcast_to(axis,coords.shape),coords).reshape(81))
                rigid_error=max(float(np.linalg.norm(ks@(u/scale))/(norm*np.linalg.norm(u/scale))) for u in rigid)
                U=np.einsum('bij,nj->nib',method.MB,q.reshape(27,3)).reshape(81,6);affine=float(np.max(abs(np.einsum('pik,kj->pij',d['raw'],U)-np.eye(6)[None])))
                p=leggauss(9)[0];points=[[sign,t,z] if axis==0 else [t,sign,z] for axis,sign in itertools.product([0,1],[-1,1]) for t,z in itertools.product(p,[-.6,.1,.8])];points+= [[-.71,.22,-.31],[.18,-.64,.73],[.62,.43,-.27]]
                _,_,boundary=method.evaluate(states[9],points);value=max(d['value_error'],fieldnorm(boundary[:,:,81:]))
                Q=rotation(R);key=(name,h)
                if not rotated and not reparameterized:baselines[key]=d
                base=baselines[key];backO=np.einsum('ai,pij,jb->pab',Q.T,d['O'],Q,optimize=True);vdef=float(np.max(abs(backO-base['O'])));fdef=raw.column_error(raw.pull_force(d['work']@Q,K3,R),base['work']);T=np.kron(np.eye(9),np.kron(K3,R));s0=1/np.sqrt(base['rawK'].diagonal());kdef=float(np.linalg.norm(s0[:,None]*(T.T@d['K']@T-base['K'])*s0[None,:])/np.linalg.norm(s0[:,None]*base['K']*s0[None,:]));gdef=float(np.max(abs(np.einsum('ab,pbc->pac',R.T,d['J'])-base['J'])))
                flat_difference=None
                if name in ['RECT','SKEW','TRAP'] and not rotated and not reparameterized:
                    expected=flat.matrices(q,h,9)['Ktie'];flat_difference=float(np.linalg.norm(scale[:,None]*(d['K']-expected)*scale[None,:])/np.linalg.norm(scale[:,None]*expected*scale[None,:]))
                edgelog={str(key):dict(ranks=e['ranks'],condition=e['condition'],missing_mode_defect=e['missing_mode_defect']) for key,e in states[9]['edges'].items()}
                checks=dict(edge_rank=all(len(set(e['ranks']))==1 and e['condition']<=1e8 and e['missing_mode_defect']<=1e-10 for e in states[9]['edges'].values()),value=value<=1e-10,work=d['work_error']<=1e-10,independent_integration=check11['work_error']<=1e-10,
                    construction_convergence=field_difference<=1e-8,affine=affine<=1e-11,objective=vdef<=1e-10 and fdef<=1e-10 and kdef<=1e-10,geometry=gdef<=1e-12,
                    rank=all(r==75 for r in ranks.values()),symmetric=symmetric<=1e-12,positive=minimum>=-1e-10*max(1.,norm),rigid=rigid_error<=1e-11,flat_inheritance=flat_difference is None or flat_difference<=1e-10)
                checks={k:bool(v) for k,v in checks.items()};file=out/f'{name}-h{h:g}-R{int(rotated)}-P{int(reparameterized)}.npz';np.savez_compressed(file,q=q,width=width,offset=offset,K=d['K'],rawK=d['rawK'],O=d['O'],work=d['work'],raw_work=d['raw_work'],mean=states[9]['mean'],moment=states[9]['moment'],bubble_mean=states[9]['bubble_mean'],edge_coefficients=np.array([states[9]['edges'][key]['coefficients'] for key in itertools.product([0,1],[-1,1])]))
                cases.append(dict(geometry=name,h=h,rotated=rotated,reparameterized=reparameterized,value_error=value,work_error=d['work_error'],work_error_11=check11['work_error'],field_7_9_difference=field_difference,affine_error=affine,rank=ranks,singular_values=sv.tolist(),minimum_eigenvalue=minimum,rigid_error=rigid_error,symmetry=symmetric,flat_difference=flat_difference,coordinate_defects=dict(value=vdef,work=fdef,stiffness=kdef,geometry=gdef),edge_observability=edgelog,checks=checks,artifact=file.relative_to(ROOT).as_posix(),artifact_sha256=digest(file)))
                if time.monotonic()-start>300 or memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
            print(json.dumps(dict(geometry=name,h=h,checks=cases[-4]['checks'],value=cases[-4]['value_error'],work11=cases[-4]['work_error_11'],rank=cases[-4]['rank'])),flush=True)
        for name,h in itertools.product(['RECT','SKEW','TRAP','SADDLE','NORMAL-C'],[.1,.01]):
            coords,ids=shared.two_elements(name);left=np.array([[1,t,z] for t,z in itertools.product([-.8,-.1,.6],[-.6,0,.7])]);right=left.copy();right[:,0]=-1;operators=[];rawops=[];Js=[]
            for i,points in enumerate([left,right]):
                q=coords[ids[i]].reshape(81);state=method.prepare(q,h,0.,9);J,original,value=method.evaluate(state,points);S=np.einsum('pi,kij,pj->pk',J[:,:,1],method.MB,J[:,:,2]);local=np.einsum('pi,pij->pj',S,value[:,:,:81]);orig=np.einsum('pi,pij->pj',S,original[:,:,:81]);op=np.zeros((9,135));rop=op.copy();dofs=(9*ids[i][:,None]+np.arange(9)).reshape(-1);op[:,dofs]=local;rop[:,dofs]=orig;operators.append(op);rawops.append(rop);Js.append(J)
            diff=float(np.linalg.norm(operators[0]-operators[1])/max(np.linalg.norm(operators[0]),np.linalg.norm(operators[1])));base=float(np.linalg.norm(rawops[0]-rawops[1])/max(np.linalg.norm(rawops[0]),np.linalg.norm(rawops[1])));geo=float(np.max(abs(Js[0][:,:,[1,2]]-Js[1][:,:,[1,2]])));pairs.append(dict(geometry=name,h=h,trace_error=diff,raw_trace_error=base,geometry_error=geo,checks_pass=bool(diff<=1e-10 and base<=1e-10 and geo<=1e-12)))
            if time.monotonic()-start>300 or memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
        assert sum(p.stat().st_size for p in out.iterdir())<=32*1024**2
    except Exception as exc:stop=repr(exc)
    result=dict(source_sha256=digest(__file__),candidate_sha256=digest(ROOT/'rd/structure/main-shell/q9-edge-trace/q9-edge-moment-candidate.md'),protocol_sha256=digest(F/'q9-edge-moment-feasibility-experiment.md'),dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [method,flat,raw,shared,geometry,material,memory]},cases=cases,shared_pairs=pairs,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=bool(stop is None and len(cases)==40 and len(pairs)==10 and all(all(c['checks'].values()) for c in cases) and all(c['checks_pass'] for c in pairs)),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_criteria_pass=result['all_criteria_pass'],cases=len(cases),pairs=len(pairs),seconds=result['seconds'],stop_reason=stop)))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
