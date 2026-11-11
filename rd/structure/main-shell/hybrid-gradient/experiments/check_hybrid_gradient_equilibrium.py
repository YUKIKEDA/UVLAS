"""Solve weak-gradient volume/trace equilibrium under physical affine tractions."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import solve
import hybrid_gradient_ancf as f
import check_q9_plate_convergence as memory
m=f.m;ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1e-30,np.linalg.norm(b)))

def setup(order=13,zorder=7,independent=False):
    grid=f.mesh(.4,2,.2);original=grid['original'];nc=len(original['cells']);nv=75*nc
    X0=np.vstack((np.vstack([original['ref'][ids] for ids in original['cells']]),grid['trace_ref']));maps=[];states=[]
    for c,ids in enumerate(original['cells']):
        indices=np.r_[np.arange(75*c,75*(c+1)),nv+grid['connect'][c]];maps.append(indices)
        state=f.prepare(original['ref'][ids].reshape(-1),.01,'S4',order,zorder,independent)
        stab=state['root'][state['elastic_rows']:];state['Kstab']=stab.T@stab;states.append(state)
    return grid,X0,np.array(maps),states,nv+grid['fixed']

def local(state,current,tangent):
    U,force,_,info=f.nonlinear(state,current,np.zeros(405),np.zeros(405));K=None
    if tangent:
        u=current-state['X0'];Gw=state['Gw'];F=np.eye(3)+np.einsum('na,pjn->paj',u,Gw)
        E=np.einsum('kij,pij->pk',m.MB,(F.transpose(0,2,1)@F-np.eye(3))/2);S=np.einsum('kij,pk->pij',m.MB,E@m.C.T)
        B=np.einsum('kij,pai,pjn->pkna',m.MB,F,Gw).reshape(-1,6,405)
        material=(np.einsum('ij,pjn->pin',m.CH,B)*np.sqrt(state['weight'])[:,None,None]).reshape(-1,405)
        geo=np.einsum('p,pin,pij,pjm->nm',state['weight'],Gw,S,Gw,optimize=True)
        K=material.T@material+np.kron(geo,np.eye(3))+state['Kstab']
    return U,force,K,info

def assemble(X,maps,states,tangent=True):
    force=np.zeros_like(X);K=np.zeros((X.size,X.size)) if tangent else None;U=0.;minimum=1e99;fluxes={};flux_residual=0.
    for ids,state in zip(maps,states):
        ui,fi,ki,info=local(state,X[ids],tangent);force[ids]+=fi.reshape(-1,3);U+=ui;minimum=min(minimum,info['min_F'])
        if tangent:
            dofs=(3*ids[:,None]+np.arange(3)).reshape(-1);K[np.ix_(dofs,dofs)]+=ki
        for face,t in zip(state['faces'],info['fluxes']):
            faceids=ids[face['offset']:face['offset']+15];key=int(faceids[0]);fluxes.setdefault(key,[]).append(t)
            weak=np.einsum('p,pn,pa->na',face['measure'],face['hat'],t)
            flux_residual=max(flux_residual,rel(weak-fi.reshape(135,3)[face['offset']:face['offset']+15],fi))
    return U,force,K,minimum,fluxes,flux_residual

def loads(X0,maps,states,fixed,P):
    external=np.zeros_like(X0);support=np.zeros_like(X0);counts={}
    for ids in maps:
        for face in range(4):key=int(ids[75+15*face]);counts[key]=counts.get(key,0)+1
    for ids,state in zip(maps,states):
        for face in state['faces']:
            nodes=ids[face['offset']:face['offset']+15]
            if counts[int(nodes[0])]!=1:continue
            force=np.einsum('p,pn,pa->na',face['measure'],face['hat'],face['normal']@P.T)
            if np.all(np.isin(nodes,fixed)):support[nodes]+=force
            else:external[nodes]+=force
        plane,w=m.quadrature(int(np.sqrt(len(state['p'])/len(np.unique(state['p'][:,2])))),2)
        for sign in [-1,1]:
            p=np.column_stack((plane,np.full(len(plane),sign)));phi,J,_=m.raw(state['q'],p,state['h']*p[:,2]/2,True)
            area=sign*np.cross(J[:,:,0],J[:,:,1]);external[ids[:75]]+=np.einsum('p,pn,pa->na',w,phi,area@P.T)
    return external,support

def main():
    dest=ROOT/'rd/structure/main-shell/hybrid-gradient/records/hybrid-gradient-equilibrium-result.json';assert not dest.exists();start=time.monotonic();history=[];stop=None;checks={};artifact=None;metrics={}
    try:
        grid,X0,maps,states,fixed=setup();free_nodes=np.setdiff1d(np.arange(len(X0)),fixed);free=(3*free_nodes[:,None]+np.arange(3)).reshape(-1)
        Fstar=np.eye(3)*1.003+np.array([[0,.0002,.0001],[.0002,0,.0001],[.0001,.0001,0]])
        E=np.einsum('kij,ij->k',m.MB,(Fstar.T@Fstar-np.eye(3))/2);S=np.einsum('kij,k->ij',m.MB,m.C@E);P=Fstar@S
        Xstar=X0@Fstar.T;external,support=loads(X0,maps,states,fixed,P);X=X0.copy();X[fixed]=Xstar[fixed];norm=np.linalg.norm(external[free_nodes])
        for iteration in range(16):
            U,force,K,minimum,fluxes,flux_dual=assemble(X,maps,states);res=force-external;rn=np.linalg.norm(res[free_nodes])/norm
            history.append(dict(iteration=iteration,residual=float(rn),min_F=minimum))
            if rn<=1e-10:break
            if iteration==15:raise ValueError('newton_iteration_limit')
            T=K[np.ix_(free,free)];scale=1/np.sqrt(np.abs(np.diag(T)));dx=scale*solve(T*scale[:,None]*scale[None,:],-res.reshape(-1)[free]*scale,assume_a='sym',check_finite=False)
            accepted=False
            for backtrack in range(13):
                alpha=2.**(-backtrack);trial=X.copy();trial[free_nodes]+=(alpha*dx).reshape(-1,3)
                _,ft,_,det,_,_=assemble(trial,maps,states,False)
                if det>0 and np.linalg.norm((ft-external)[free_nodes])/norm<rn:X=trial;history[-1]['alpha']=alpha;accepted=True;break
            if not accepted:raise ValueError('newton_line_search_failed')
        reaction=res.copy();reaction[free_nodes]=0
        _,_,maps_i,states_i,fixed_i=setup(15,9,True);ext_i,support_i=loads(X0,maps_i,states_i,fixed_i,P);_,force_i,_,_,_,_=assemble(X,maps_i,states_i,False)
        internal=[np.max(np.linalg.norm(v[0]+v[1],axis=1)) for v in fluxes.values() if len(v)==2]
        v=np.random.default_rng(204).normal(size=(2,405));_,_,ke,_=local(states[0],X[maps[0]],True);_,_,bilinear,_=f.nonlinear(states[0],X[maps[0]],v[0],v[1])
        metrics=dict(displacement=rel(X-Xstar,Xstar-X0),equilibrium=rel(res[free_nodes],external[free_nodes]),support_reaction=rel(reaction-support_i,support_i),
            independent_equilibrium=rel((force_i-ext_i)[free_nodes],ext_i[free_nodes]),internal_pointwise_flux=float(max(internal)/np.linalg.norm(P)),flux_dual=flux_dual,
            tangent=float(abs(v[0]@ke@v[1]-bilinear)/max(1e-30,abs(bilinear))),min_F=minimum)
        checks=dict(newton=history[-1]['residual']<=1e-10,displacement=metrics['displacement']<=1e-7,
            identities=all(v<=1e-8 for k,v in metrics.items() if k not in ['displacement','min_F']),orientation=minimum>0,
            resources=time.monotonic()-start<=180 and memory.memory_bytes()<=2*1024**3)
        path=ROOT/'build/analysis'/('hybrid-gradient-equilibrium-'+time.strftime('%Y%m%d-%H%M%S')+'.npz')
        np.savez_compressed(path,ref=X0,current=X,target=Xstar,maps=maps,force=force,external=external,reaction=reaction,support=support_i)
        artifact=dict(path=path.relative_to(ROOT).as_posix(),sha256=digest(path))
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(f.__file__),ROOT/'rd/structure/main-shell/hybrid-gradient/records/hybrid-gradient-equilibrium-protocol.md',ROOT/'rd/structure/main-shell/hybrid-gradient/records/hybrid-gradient-equilibrium-preparation.json',ROOT/'rd/structure/main-shell/hybrid-gradient/records/hybrid-gradient-local-result.json']
    dest.write_text(json.dumps(dict(history=history,metrics=metrics,checks={k:bool(v) for k,v in checks.items()},all_checks_pass=bool(checks) and all(checks.values()) and stop is None,
        stop_reason=stop,artifact=artifact,manufactured_Ku_load=False,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
    print(json.dumps(dict(history=history,metrics=metrics,checks=checks,stop_reason=stop,seconds=time.monotonic()-start),indent=2))
if __name__=='__main__':main()
