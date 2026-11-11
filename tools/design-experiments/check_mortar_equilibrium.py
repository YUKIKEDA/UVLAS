"""Nonlinear constrained equilibrium with physical outer tractions, no K*u load."""
from pathlib import Path
import hashlib,json,itertools,time
import numpy as np
from scipy.linalg import svd,solve
import mortar_q4 as f
import check_q9_plate_convergence as memory
m=f.m;ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1e-30,np.linalg.norm(b)))

def element(q0,q,h,tangent=True):
    plane,wp=m.quadrature(9,2);z,wz=np.polynomial.legendre.leggauss(5)
    p=np.array([[s,t,c] for (s,t),c in itertools.product(plane,z)]);w=np.array([a*b for a,b in itertools.product(wp,wz)])
    _,g,_=m.tables(p,h*p[:,2]/2);J=np.einsum('na,pnb->pab',q0.reshape(75,3),g);G=g@np.linalg.inv(J);mu=w*np.linalg.det(J)*h/2
    F=np.eye(3)+np.einsum('na,pnb->pab',(q-q0).reshape(75,3),G)
    E=np.einsum('kij,pij->pk',m.MB,(F.transpose(0,2,1)@F-np.eye(3))/2);stress=E@m.C.T
    B=np.einsum('kij,pai,pnj->pkna',m.MB,F,G).reshape(-1,6,225)
    force=np.einsum('p,pki,pk->i',mu,B,stress);U=np.einsum('p,pk,pk->',mu,E,stress)/2
    K=None
    if tangent:
        K=np.einsum('p,pki,kl,plj->ij',mu,B,m.C,B,optimize=True)
        tensor=np.einsum('kij,pk->pij',m.MB,stress);geo=np.einsum('p,pni,pij,pmj->nm',mu,G,tensor,G,optimize=True)
        K+=np.kron(geo,np.eye(3))
    return U,force,K,float(np.min(np.linalg.det(F)))

def assemble(mesh,q,h,tangent=True):
    U=0.;force=np.zeros(q.size);K=np.zeros((q.size,q.size)) if tangent else None;minimum=1e99
    for ids in mesh['cells']:
        indices=(3*ids[:,None]+np.arange(3)).reshape(-1);u,fe,ke,det=element(mesh['ref'][ids].reshape(-1),q[ids].reshape(-1),h,tangent)
        U+=u;force[indices]+=fe;minimum=min(minimum,det)
        if tangent:K[np.ix_(indices,indices)]+=ke
    return U,force,K,minimum

def loads(mesh,h,P,independent=False):
    external=np.zeros_like(mesh['ref']);support=np.zeros_like(external)
    for c,ids in enumerate(mesh['cells']):
        j,i=divmod(c,mesh['n']);faces=[(2,-1),(2,1)]
        if i==mesh['n']-1:faces.append((0,1))
        if j==0:faces.append((1,-1))
        if j==mesh['n']-1:faces.append((1,1))
        q=mesh['ref'][ids].reshape(-1);external[ids]+=f.boundary(q,h,P,faces,11 if independent else 9,independent)
        if i==0:support[ids]+=f.boundary(q,h,P,[(0,-1)],11 if independent else 9,independent)
    return external,support

def main():
    dest=ROOT/'docs/verification/mortar-q4-equilibrium-result.json';assert not dest.exists();start=time.monotonic();history=[];stop=None;checks={};artifact=None
    try:
        mesh=f.broken_mesh(.4,2,.2);h=.01;degree=1;G,meta=f.constraints(mesh,h,degree)
        norm=np.linalg.norm(G,axis=1);active=norm>1e-12*norm.max();A=G[active]/norm[active,None]
        U,s,V=svd(A,full_matrices=False,check_finite=False);r=f.ranks(s)[1]
        if f.ranks(s)!=[r]*3:raise ValueError('constraint_rank_uncertain')
        C=V[:r];free_nodes=np.concatenate(mesh['free']);fixed_nodes=np.setdiff1d(np.arange(len(mesh['ref'])),free_nodes)
        Cf=C[:,free_nodes];Z,zm=f.null_rows(Cf);Z3=np.kron(Z,np.eye(3));free=(3*free_nodes[:,None]+np.arange(3)).reshape(-1)
        if zm['rank']!=r:raise ValueError('supported_multiplier_rank_failed')
        Fstar=np.eye(3)*1.003+np.array([[0,.0002,.0001],[.0002,0,.0001],[.0001,.0001,0]])
        E=np.einsum('kij,ij->k',m.MB,(Fstar.T@Fstar-np.eye(3))/2);S=np.einsum('kij,k->ij',m.MB,m.C@E);P=Fstar@S
        qstar=mesh['ref']@Fstar.T;external,support=loads(mesh,h,P);external_i,support_i=loads(mesh,h,P,True)
        bc=qstar[fixed_nodes]-mesh['ref'][fixed_nodes];lift=np.linalg.lstsq(Cf,-C[:,fixed_nodes]@bc,rcond=1e-10)[0]
        q=mesh['ref'].copy();q[fixed_nodes]=qstar[fixed_nodes];q[free_nodes]+=lift
        fnorm=np.linalg.norm(external[free_nodes]);load=external.reshape(-1)
        for iteration in range(16):
            U,force,K,minimum=assemble(mesh,q,h);res=force-load;reduced=Z3.T@res[free];rn=np.linalg.norm(reduced)/fnorm
            history.append(dict(iteration=iteration,residual=float(rn),min_F=minimum))
            if rn<=1e-10:break
            if iteration==15:raise ValueError('newton_iterations_exceeded')
            Kr=Z3.T@K[np.ix_(free,free)]@Z3;scale=1/np.sqrt(np.abs(np.diag(Kr)))
            dz=scale*solve(Kr*scale[:,None]*scale[None,:],-reduced*scale,assume_a='sym',check_finite=False)
            step=Z3@dz;accepted=False
            for backtrack in range(13):
                alpha=2.**(-backtrack);trial=q.copy();trial[free_nodes]+=(alpha*step).reshape(-1,3)
                _,ft,_,det=assemble(mesh,trial,h,False);rt=np.linalg.norm(Z3.T@(ft-load)[free])/fnorm
                if det>0 and rt<rn:q=trial;accepted=True;history[-1]['alpha']=alpha;break
            if not accepted:raise ValueError('newton_line_search_failed')
        multiplier=np.linalg.solve(Cf@Cf.T,-Cf@res.reshape(-1,3)[free_nodes]);reaction=C.T@multiplier
        residual=force.reshape(-1,3)+reaction-external;support_computed=residual.copy();support_computed[free_nodes]=0
        force_i=np.zeros_like(mesh['ref'])
        for ids in mesh['cells']:
            _,fi,_=f.raw_nonlinear(mesh['ref'][ids].reshape(-1),q[ids].reshape(-1),h,np.zeros(225),np.zeros(225),11,7,True);force_i[ids]+=fi.reshape(75,3)
        Gi,mi=f.constraints(mesh,h,degree,11,7,True);lm=np.zeros((len(Gi),3))
        for a,entry in enumerate(mi):
            if entry['kind']=='piola':lm[a]=-h*P[:,entry['component']]
        expected_interface=Gi.T@lm
        v=np.random.default_rng(19).normal(size=(2,225));q0=mesh['ref'][mesh['cells'][0]].reshape(-1);qe=q[mesh['cells'][0]].reshape(-1)
        _,_,ke,_=element(q0,qe,h);_,_,analytic=f.raw_nonlinear(q0,qe,h,v[0],v[1]);tangent_error=abs(v[0]@ke@v[1]-analytic)/max(1e-30,abs(analytic))
        metrics=dict(displacement=rel(q-qstar,qstar-mesh['ref']),constraint=rel(Gi@q,mesh['ref']),free_equilibrium=rel(residual[free_nodes],external[free_nodes]),
            support_reaction=rel(support_computed-support_i,support_i),interface_reaction=rel(reaction-expected_interface,expected_interface),
            independent_equilibrium=rel((force_i+reaction-external_i)[free_nodes],external_i[free_nodes]),tangent=float(tangent_error),min_F=minimum)
        checks=dict(newton=history[-1]['residual']<=1e-10,displacement=metrics['displacement']<=1e-7,
            identities=all(v<=1e-8 for k,v in metrics.items() if k not in ['displacement','min_F']),orientation=minimum>0,
            resources=time.monotonic()-start<=180 and memory.memory_bytes()<=2*1024**3)
        path=ROOT/'build/analysis'/('mortar-equilibrium-'+time.strftime('%Y%m%d-%H%M%S')+'.npz')
        np.savez_compressed(path,ref=mesh['ref'],q=q,qstar=qstar,cells=mesh['cells'],G=G,C=C,multiplier=multiplier,reaction=reaction,external=external,support=support_i,force=force)
        artifact=dict(path=path.relative_to(ROOT).as_posix(),sha256=digest(path))
    except Exception as exc:stop=repr(exc);metrics={}
    sources=[Path(__file__),Path(f.__file__),ROOT/'docs/verification/mortar-q4-equilibrium-protocol.md',ROOT/'docs/verification/mortar-q4-equilibrium-preparation.json',ROOT/'docs/verification/mortar-q4-result.json']
    result=dict(history=history,metrics=metrics,checks={k:bool(v) for k,v in checks.items()},all_checks_pass=bool(checks) and all(checks.values()) and stop is None,
        stop_reason=stop,artifact=artifact,manufactured_Ku_load=False,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes())
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps({k:result[k] for k in ['history','metrics','checks','stop_reason','seconds']},indent=2))
if __name__=='__main__':main()
