"""Predeclared curved constrained-gradient local/trace gates."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
import q9_curved_constrained_gradient as method
import q9_edge_trace_gradient as flat
import check_q9_curved_load as geometry
import check_q9_curved_strain_boundary as shared
import check_q9_curved_strain_screen as raw
import check_q9_independent_load as material
import check_q9_plate_convergence as memory
from q9_edge_moment_fields import quadrature
from check_q9_edge_moment_feasibility import rotation
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def relative(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-30))
def measure(state,order):
    points,w=quadrature(order);J,E,B=method.evaluate(state,points);weight=w*np.linalg.det(J);V=weight.sum();phi=method.modes(state,points,E);phi_norm=np.sqrt(np.einsum('p,pki,pki->i',weight,phi,phi));source_norm=np.sqrt(np.einsum('p,pkn,pkn->n',weight,E,E));projection=np.einsum('p,pki,pkn->in',weight,phi,B-E,optimize=True)
    return dict(points=points,J=J,E=E,B=B,weight=weight,phi=phi,mean=float(np.max(np.linalg.norm(np.einsum('p,pki->ki',weight,phi),axis=0)/(np.sqrt(V)*phi_norm))),restore=float(np.max(abs(projection)/(phi_norm[:,None]*source_norm[None,:]))),value=float(np.max(np.linalg.norm(B[:,:,81:]-np.eye(6)[None],axis=1))),work=relative(np.einsum('p,pkn->kn',weight,B),np.einsum('p,pkn->kn',weight,E)),raw_mode_mean=float(np.linalg.norm(state['mode_mean'])))
def main():
    target=F/'q9-curved-constrained-gradient-result.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-edge-trace-gradient-fine-audit.json';assert json.loads(parent.read_text())['all_criteria_pass']
    start=time.monotonic();out=ROOT/'build/analysis'/('q9-curved-constrained-gradient-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();cases=[];pairs=[];baselines={};stop=None
    D=np.diag([1,1,1,np.sqrt(2),np.sqrt(2),np.sqrt(2)]);C=D@material.elastic@D
    try:
        for name,h in itertools.product(['RECT','SKEW','TRAP','SADDLE','NORMAL-C'],[.1,.01]):
            for rotated,chart in itertools.product([False,True],repeat=2):
                q0=geometry.geometry(name);R=geometry.ROTATION if rotated else np.eye(3);q=geometry.transform(q0,R);K3=np.eye(3);width=h;offset=0.
                if chart:
                    a=1.3;b=.1*h;K3=np.array([[1,b,b*b/2],[0,a,a*b],[0,0,a*a]]);q=np.einsum('ab,nbc->nac',K3,q.reshape(9,3,3)).reshape(81);width=h/a;offset=-b/a
                state=method.prepare(q,width,offset,9);d=measure(state,9);ind=measure(state,11);B=d['B'][:,:,:81];E=d['E'][:,:,:81];weight=d['weight'];K=np.einsum('p,pik,ij,pjl->kl',weight,B,C,B,optimize=True);rawK=np.einsum('p,pik,ij,pjl->kl',weight,E,C,E,optimize=True);scale=1/np.sqrt(rawK.diagonal());ks=scale[:,None]*K*scale[None,:];s=np.linalg.svd(ks,compute_uv=False);ranks={str(t):int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]};sym=relative(ks,ks.T);minimum=float(np.linalg.eigvalsh((ks+ks.T)/2)[0]);coords=q.reshape(9,3,3);modes=[]
                for axis in np.eye(3):
                    u=np.zeros_like(coords);u[:,0,:]=axis;modes.extend([u.reshape(81),np.cross(np.broadcast_to(axis,coords.shape),coords).reshape(81)])
                rigid=max(float(np.linalg.norm(ks@(u/scale))/(s[0]*np.linalg.norm(u/scale))) for u in modes)
                ep=[];axes=[]
                for axis,sign in itertools.product(range(2),[-1,1]):
                    for t,z in itertools.product([-.8,-.1,.6],[-.7,0.,.4]):ep.append([t,sign,z] if axis==0 else [sign,t,z]);axes.append(axis)
                Je,Ee,_=method.evaluate(state,ep);phie=method.modes(state,ep,Ee);tan=Je[np.arange(len(ep)),:,axes];row=np.einsum('pi,kij,pj->pk',tan,method.MB,Je[:,:,2]);trace=float(np.linalg.norm(np.einsum('pk,pki->pi',row,phie))/max(np.linalg.norm(tan,axis=1).max()*np.linalg.norm(Je[:,:,2],axis=1).max()*np.linalg.norm(phie),1e-30))
                key=(name,h)
                if not rotated and not chart:baselines[key]=dict(K=K,rawK=rawK,B=B.copy(),O=d['B'][:,:,81:].copy())
                base=baselines[key];T=np.kron(np.eye(9),np.kron(K3,R));sb=1/np.sqrt(base['rawK'].diagonal());kdef=relative(sb[:,None]*(T.T@K@T)*sb[None,:],sb[:,None]*base['K']*sb[None,:]);Q=rotation(R);bdef=relative(np.einsum('pki,ij->pkj',B,T),np.einsum('ij,pjk->pik',Q,base['B']));odef=relative(np.einsum('ij,pjk,kl->pil',Q.T,d['B'][:,:,81:],Q),base['O'])
                flat_difference=None
                if name in ['RECT','SKEW','TRAP'] and not rotated and not chart:
                    expected=flat.evaluate(flat.prepare(q,width,offset,9),d['points'])[2];flat_difference=relative(d['B'],expected)
                checks=dict(gram=state['mode_gram_ranks']==[3,3,3] and state['mode_gram_condition']<=1e8,mean=ind['mean']<=1e-10,trace=trace<=1e-10,value=ind['value']<=1e-10,work11=ind['work']<=1e-10,restore11=ind['restore']<=1e-10,rank=all(v==75 for v in ranks.values()),rigid=rigid<=1e-11,symmetry=sym<=1e-12,positive=minimum>=-1e-10*s[0],coordinates=max(kdef,bdef,odef)<=1e-10,flat=flat_difference is None or flat_difference<=1e-10)
                path=out/f'{name}-h{h:g}-R{int(rotated)}-P{int(chart)}.npz';np.savez_compressed(path,K=K,rawK=rawK,gram=state['mode_gram'],mode_mean=state['mode_mean'],mode_alpha=state['mode_alpha'],singular_values=s)
                record=dict(geometry=name,h=h,rotated=rotated,chart=chart,gram_condition=state['mode_gram_condition'],mean11=ind['mean'],raw_mode_mean=ind['raw_mode_mean'],trace=trace,value11=ind['value'],work11=ind['work'],restore11=ind['restore'],ranks=ranks,rigid=rigid,symmetry=sym,minimum_eigenvalue=minimum,coordinate_defects=dict(stiffness=kdef,field=bdef,constant=odef),flat_difference=flat_difference,checks={k:bool(v) for k,v in checks.items()},artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path));cases.append(record)
                print(json.dumps(dict(geometry=name,h=h,rotated=rotated,chart=chart,checks=record['checks'],restore11=ind['restore'])),flush=True)
                if time.monotonic()-start>300 or memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
        for name,h in itertools.product(['RECT','SKEW','TRAP','SADDLE','NORMAL-C'],[.1,.01]):
            coords,ids=shared.two_elements(name);left=np.array([[1,t,z] for t,z in itertools.product([-.8,-.1,.6],[-.7,0.,.4])]);right=left.copy();right[:,0]=-1;operators=[]
            for i,points in enumerate([left,right]):
                q=coords[ids[i]].reshape(81);state=method.prepare(q,h,0.,9);J,_,B=method.evaluate(state,points);row=np.einsum('pi,kij,pj->pk',J[:,:,1],method.MB,J[:,:,2]);local=np.einsum('pk,pkn->pn',row,B[:,:,:81]);op=np.zeros((len(points),135));dofs=(9*ids[i][:,None]+np.arange(9)).reshape(-1);op[:,dofs]=local;operators.append(op)
            error=relative(operators[0],operators[1]);pairs.append(dict(geometry=name,h=h,trace=error,checks_pass=bool(error<=1e-10)))
        assert sum(p.stat().st_size for p in out.iterdir())<=32*1024**2
    except Exception as exc:stop=repr(exc)
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(parent),contract_sha256=digest(ROOT/'rd/structure/main-shell/q9-curved/q9-curved-constrained-gradient-candidate.md'),protocol_sha256=digest(F/'q9-curved-constrained-gradient-experiment.md'),dependencies={Path(m.__file__).name:digest(m.__file__) for m in [method,method.base,flat,shared,raw,geometry,material,memory]},cases=cases,shared_pairs=pairs,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=bool(stop is None and len(cases)==40 and len(pairs)==10 and all(all(c['checks'].values()) for c in cases) and all(c['checks_pass'] for c in pairs)),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['cases','shared_pairs','dependencies']}))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
