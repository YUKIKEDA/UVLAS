"""Physical outer-face loads, including independent top/bottom position traces."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import solve
import reference_work_fit as fit
import check_q9_plate_convergence as memory
f=fit.f;m=f.m;ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1e-30,np.linalg.norm(b)))
def setup(order=13,zorder=7,independent=False):
    original=m.mesh(.4,2,.2);states=[];registry={};trace_ref=[];fixed=[];connect=[];internals=[];X=[];offset=0
    for cell,ids in enumerate(original['cells']):
        s=f.prepare(original['ref'][ids].reshape(-1),.01,order=order,zorder=zorder,independent=independent);nodes,T=fit.register(s,fit.edge_keys(cell,2),registry,trace_ref,fixed)
        if np.linalg.norm(T-np.eye(len(T)))>1e-10:raise ValueError('equilibrium_case_face_orientation')
        connect.append(nodes);internals.append(np.arange(offset,offset+s['ni']));offset+=s['ni'];X.extend(s['X0'][:s['ni']]);stab=s['root'][s['elastic_rows']:];s['Kstab']=stab.T@stab;states.append(s)
    X.extend(trace_ref);maps=[np.r_[inside,offset+nodes] for inside,nodes in zip(internals,connect)]
    return np.array(X),np.array(maps),states,offset+np.array(fixed)
def local(s,current,tangent):
    nc=3*s['nscalar'];U,force,_,info=f.nonlinear(s,current,np.zeros(nc),np.zeros(nc));K=None
    if tangent:
        Gw=s['Gw'];F=np.eye(3)+np.einsum('na,pjn->paj',current-s['X0'],Gw);E=np.einsum('kij,pij->pk',m.MB,(F.transpose(0,2,1)@F-np.eye(3))/2);S=np.einsum('kij,pk->pij',m.MB,E@m.C.T)
        B=np.einsum('kij,pai,pjn->pkna',m.MB,F,Gw).reshape(-1,6,nc);material=(np.einsum('ij,pjn->pin',m.CH,B)*np.sqrt(s['weight'])[:,None,None]).reshape(-1,nc)
        geo=np.einsum('p,pin,pij,pjm->nm',s['weight'],Gw,S,Gw,optimize=True);K=material.T@material+np.kron(geo,np.eye(3))+s['Kstab']
    return U,force,K,info
def assemble(X,maps,states,tangent=True):
    force=np.zeros_like(X);K=np.zeros((X.size,X.size)) if tangent else None;U=0.;minimum=1e99;fluxes={};dual=0.
    for ids,s in zip(maps,states):
        u,fi,ki,info=local(s,X[ids],tangent);force[ids]+=fi.reshape(-1,3);U+=u;minimum=min(minimum,info['min_F'])
        if tangent:
            dofs=(3*ids[:,None]+np.arange(3)).reshape(-1);K[np.ix_(dofs,dofs)]+=ki
        for face,tr in zip(s['faces'],info['fluxes']):
            sl=slice(face['offset'],face['offset']+face['n']);nodes=ids[sl];fluxes.setdefault(int(nodes[0]),[]).append(tr)
            weak=np.einsum('p,pn,pa->na',face['measure'],face['hat'],tr);dual=max(dual,rel(weak-fi.reshape(-1,3)[sl],fi))
    return U,force,K,minimum,fluxes,dual
def loads(X0,maps,states,fixed,P):
    ext=np.zeros_like(X0);support=np.zeros_like(X0);counts={}
    for ids,s in zip(maps,states):
        for face in s['faces']:
            key=int(ids[face['offset']]);counts[key]=counts.get(key,0)+1
    for ids,s in zip(maps,states):
        for face in s['faces']:
            nodes=ids[face['offset']:face['offset']+face['n']]
            if counts[int(nodes[0])]!=1:continue
            force=np.einsum('p,pn,pa->na',face['measure'],face['hat'],face['normal']@P.T)
            if np.all(np.isin(nodes,fixed)):support[nodes]+=force
            else:ext[nodes]+=force
    return ext,support
def main():
    dest=ROOT/'docs/verification/reference-work-v3-equilibrium-result.json';assert not dest.exists();start=time.monotonic();history=[];stop=None;checks={};metrics={};artifact=None
    try:
        parent=json.loads((ROOT/'docs/verification/reference-work-v3-shared-result.json').read_text());assert parent['all_global_checks'] and parent['all_curved_improved']
        X0,maps,states,fixed=setup();free_nodes=np.setdiff1d(np.arange(len(X0)),fixed);free=(3*free_nodes[:,None]+np.arange(3)).reshape(-1)
        Fstar=np.eye(3)*1.003+np.array([[0,.0002,.0001],[.0002,0,.0001],[.0001,.0001,0]]);E=np.einsum('kij,ij->k',m.MB,(Fstar.T@Fstar-np.eye(3))/2);P=Fstar@np.einsum('kij,k->ij',m.MB,m.C@E)
        Xstar=X0@Fstar.T;ext,support=loads(X0,maps,states,fixed,P);X=X0.copy();X[fixed]=Xstar[fixed];norm=np.linalg.norm(ext[free_nodes])
        for iteration in range(16):
            U,force,K,minimum,fluxes,dual=assemble(X,maps,states);res=force-ext;rn=np.linalg.norm(res[free_nodes])/norm;history.append(dict(iteration=iteration,residual=float(rn),min_F=minimum))
            if rn<=1e-10:break
            if iteration==15:raise ValueError('newton_iteration_limit')
            T=K[np.ix_(free,free)];scale=1/np.sqrt(np.abs(np.diag(T)));dx=scale*solve(T*scale[:,None]*scale[None,:],-res.reshape(-1)[free]*scale,assume_a='sym',check_finite=False)
            accepted=False
            for backtrack in range(13):
                alpha=2.**(-backtrack);trial=X.copy();trial[free_nodes]+=(alpha*dx).reshape(-1,3);_,ft,_,det,_,_=assemble(trial,maps,states,False)
                if det>0 and np.linalg.norm((ft-ext)[free_nodes])/norm<rn:X=trial;history[-1]['alpha']=alpha;accepted=True;break
            if not accepted:raise ValueError('newton_line_search_failed')
        reaction=res.copy();reaction[free_nodes]=0;_,mi,si,fi=setup(15,9,True);ei,supi=loads(X0,mi,si,fi,P);_,force_i,_,_,_,_=assemble(X,mi,si,False)
        internal=[np.max(np.linalg.norm(v[0]+v[1],axis=1)) for v in fluxes.values() if len(v)==2]
        v=np.random.default_rng(204).normal(size=(2,3*states[0]['nscalar']));_,_,ke,_=local(states[0],X[maps[0]],True);_,_,bilinear,_=f.nonlinear(states[0],X[maps[0]],v[0],v[1])
        metrics=dict(displacement=rel(X-Xstar,Xstar-X0),equilibrium=rel(res[free_nodes],ext[free_nodes]),support_reaction=rel(reaction-supi,supi),independent_equilibrium=rel((force_i-ei)[free_nodes],ei[free_nodes]),
            internal_pointwise_flux=float(max(internal)/np.linalg.norm(P)),flux_dual=dual,tangent=float(abs(v[0]@ke@v[1]-bilinear)/max(1e-30,abs(bilinear))),min_F=minimum)
        checks=dict(newton=history[-1]['residual']<=1e-10,displacement=metrics['displacement']<=1e-7,identities=all(v<=1e-8 for k,v in metrics.items() if k not in ['displacement','min_F']),orientation=minimum>0,
            resources=time.monotonic()-start<=180 and memory.memory_bytes()<=2*1024**3)
        path=ROOT/'build/analysis'/('reference-work-v3-equilibrium-'+time.strftime('%Y%m%d-%H%M%S')+'.npz');np.savez_compressed(path,ref=X0,current=X,target=Xstar,maps=maps,force=force,external=ext,reaction=reaction,support=supi);artifact=dict(path=path.relative_to(ROOT).as_posix(),sha256=digest(path))
    except Exception as exc:stop=repr(exc)
    ps=[Path(__file__),Path(f.__file__),Path(fit.__file__),ROOT/'docs/verification/reference-work-v3-protocol.md',ROOT/'docs/verification/reference-work-v3-equilibrium-preparation.json']
    dest.write_text(json.dumps(dict(history=history,metrics=metrics,checks={k:bool(v) for k,v in checks.items()},all_checks_pass=bool(checks) and all(checks.values()) and stop is None,stop_reason=stop,artifact=artifact,manufactured_Ku_load=False,
        physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in ps},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8');print(json.dumps(dict(history=history,metrics=metrics,checks=checks,stop_reason=stop),indent=2))
if __name__=='__main__':main()
