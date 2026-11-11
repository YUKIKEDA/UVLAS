"""Same-space nonlinear, layer, inertia, signed spectrum and restart gates."""
from pathlib import Path
import copy,gc,hashlib,json,time,traceback
import numpy as np
from scipy.linalg import solve,eigh,svd,cholesky,solve_triangular,qr
import mass_paired_mechanics as m
from pg_ritz_state_reference import Reference,digest
ROOT=Path(__file__).resolve().parents[2]
def rel(a,b):return float(np.linalg.norm(a)/max(1e-30,np.linalg.norm(b)))

def main():
    dest=ROOT/'docs/verification/mass-paired-mechanics-result.json';assert not dest.exists()
    start=time.monotonic();records=[];stops=[];artifact=None
    try:
        model=m.shell.prepare(1,(.4,.2,.07,.1,-.06,.025),.05,True,'extended');q0=model['reference_coordinates']
        data=m.prepare(model,.05,order=45,zorder=5);Zs,ranks,trace=m.fixed_trace_kernel(model);Z=np.kron(Zs,np.eye(3))
        ua=q0@(m.FSTAR-np.eye(3)).T;external,fixed=m.external(model,.05,order=47,zorder=8)
        U,force,K,info=m.evaluate(data,ua,True);affine=rel(force-external-fixed,force)
        assert affine<=1e-8
        rng=np.random.default_rng(20261105);direction=rng.normal(size=ua.shape)
        direction*=.003/np.max(abs(np.einsum('na,pnb->pab',direction,data['G'])))
        u=ua+direction;U,f,K,state=m.evaluate(data,u,True);eps=1e-28
        Uc,fc,_,_=m.evaluate(data,u+1j*eps*direction)
        metrics=dict(energy_derivative=rel(np.array(Uc.imag/eps-np.sum(f*direction)),np.array(np.sum(f*direction))),
            force_derivative=rel(fc.imag/eps-(K@direction.ravel()).reshape(u.shape),fc.imag/eps),symmetric=rel(K-K.T,K))
        angle=.53;R=np.array([[np.cos(angle),0,np.sin(angle)],[0,1,0],[-np.sin(angle),0,np.cos(angle)]])
        Ur,fr,_,_=m.evaluate(data,(q0+u)@R.T-q0)
        metrics.update(objective_energy=abs(Ur-U)/abs(U),objective_force=rel(fr-f@R.T,fr))
        assert max(metrics.values())<=1e-8 and data['natural_metric_min']>0
        recovery=m.interface_recovery(model,.05,u)
        records.append(dict(case='nonlinear_layered_material',metrics=metrics,natural_metric_min=data['natural_metric_min'],recovery=recovery))
        print('MATERIAL',metrics,flush=True)
        split=m.prepare(model,.05,True,order=45,zorder=5);Us,fs,Ks,_=m.evaluate(split,u,True)
        metrics=dict(energy=abs(Us-U)/abs(U),force=rel(fs-f,f),tangent=rel(Ks-K,K),mass=rel(split['M']-data['M'],data['M']))
        assert max(metrics.values())<=1e-8
        records.append(dict(case='identical_layer_split',metrics=metrics));print('SPLIT',metrics,flush=True)
        del split,fs,Ks;gc.collect()
        x=rng.normal(size=Z.shape[1]);du=(Z@x).reshape(ua.shape)
        du*=.002/np.max(abs(np.einsum('na,pnb->pab',du,data['G'])));current=ua+du;history=[]
        for iteration in range(10):
            U,fi,Ki,st=m.evaluate(data,current,True);g=Z.T@(fi-external).ravel();Kr=Z.T@Ki@Z
            scale=1/np.sqrt(np.maximum(abs(np.diag(Kr)),1e-30));res=float(np.linalg.norm(scale*g)/np.linalg.norm(scale*(Z.T@external.ravel())))
            history.append(dict(iteration=iteration,residual=res,min_det=st['min_det']));print('NEWTON',history[-1],flush=True)
            if res<=1e-8:break
            dx=scale*solve(Kr*scale[:,None]*scale[None,:],-scale*g,assume_a='sym');dq=(Z@dx).reshape(current.shape)
            merit=float(U-np.sum(external*current));accepted=False
            for bt in range(14):
                trial=current+dq/2**bt;Ut,_,_,stt=m.evaluate(data,trial)
                if stt['min_det']>0 and Ut-np.sum(external*trial)<=merit+1e-14:current=trial;accepted=True;break
            assert accepted,'nonlinear_line_search'
        U,fi,K,st=m.evaluate(data,current,True)
        metrics=dict(free_residual=history[-1]['residual'],physical_displacement_error=float(np.max(abs(data['N']@(current-ua)))),affine_traction_balance=affine)
        assert max(metrics.values())<=1e-8
        records.append(dict(case='nonlinear_physical_traction_equilibrium',metrics=metrics,history=history,trace_ranks=ranks,trace_holdout=trace))
        fullM=np.kron(data['M'],np.eye(3));Mr=Z.T@fullM@Z;Kr=Z.T@K@Z
        values,vectors=eigh(Kr,Mr,check_finite=False)
        back=rel(Kr@vectors-(Mr@vectors)*values[None,:],Kr@vectors);orth=rel(vectors.T@Mr@vectors-np.eye(len(values)),np.eye(len(values)))
        assert values.min()>0 and max(back,orth)<=1e-8
        # Stream the unstressed strain operator; no giant assembled strain matrix.
        L=cholesky(fullM,lower=True);root=None;G=data['G']
        for begin in range(0,len(G),384):
            sl=slice(begin,begin+384);B=np.einsum('kij,pnj->pkni',m.MB,G[sl]).reshape(-1,6,3*model['ns'])
            CH=np.linalg.cholesky(data['C'][sl]).transpose(0,2,1)
            A=(CH@B*np.sqrt(data['weights'][sl])[:,None,None]).reshape(-1,3*model['ns'])
            root=qr(A if root is None else np.vstack((root,A)),mode='r',check_finite=False)[0][:3*model['ns']]
        sv=svd(solve_triangular(L,root.T,lower=True,check_finite=False).T,compute_uv=False)
        nullities=[int(sum(sv<=tol*sv[0])) for tol in [1e-9,1e-10,1e-11]];assert nullities==[6]*3
        trans=solve(data['M'],data['N'].T@(data['weights']*data['rho']),assume_a='pos')
        rigid=np.zeros((model['ns'],3,6));rigid[:,:,:3]=trans[:,None,None]*np.eye(3)[None]
        physical=np.zeros((len(data['X']),3,6));physical[:,:,:3]=np.eye(3)[None]
        for j in range(3):
            rigid[:,:,3+j]=np.cross(np.eye(3)[j],q0);physical[:,:,3+j]=np.cross(np.eye(3)[j],data['X'])
        Q=rigid.reshape(3*model['ns'],6);direct=np.einsum('p,pai,paj->ij',data['weights']*data['rho'],physical,physical)
        inertia=rel(Q.T@fullM@Q-direct,direct);assert inertia<=1e-9
        records.append(dict(case='mass_spectrum_inertia',eigen_backward=back,mass_orthogonality=orth,spatial_inertia=inertia,
            free_nullities=nullities,eigen_min=float(values.min()),eigen_max=float(values.max()),density_mass_bounds=np.linalg.eigvalsh(data['M'])[[0,-1]].tolist()))
        # Compression must survive the same material and support reduction.
        estar=np.einsum('kij,ij->k',m.MB,(m.FSTAR.T@m.FSTAR-np.eye(3))/2)
        negative=np.einsum('kij,ij->k',m.MB,np.diag([-.1,-.1,-.1]))
        neg=dict(data,e0=estar-np.linalg.solve(data['C'],np.broadcast_to(negative,(len(G),6))[:,:,None])[:,:,0])
        _,_,Kn,_=m.evaluate(neg,ua,True);negative_values=eigh(Z.T@Kn@Z,Mr,eigvals_only=True)
        assert negative_values.min()<0
        records.append(dict(case='negative_prestress_retained',negative_count=int(sum(negative_values<0)),minimum=float(negative_values.min()),filtered=False))
        out=ROOT/'build/analysis'/('mass-paired-mechanics-'+time.strftime('%Y%m%d-%H%M%S')+'.npz')
        np.savez_compressed(out,q0=q0,equilibrium=current,Z=Z,M=Mr,K=Kr,omega2=values,modes=vectors,fullM=fullM,
            Rcoeff=model['Rcoeff'],Dcoeff=model['Dcoeff'],Ccoeff=model['Ccoeff'],parameters=model['parameters'])
        artifact=dict(artifact=out.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(out.read_bytes()).hexdigest())
        ref=Reference(Mr,Kr,vectors,values,artifact['artifact_sha256']);k=np.arange(1,len(values)+1);amp=1e-5/k**2
        initial=ref.state(vectors@amp,vectors@(amp*np.sqrt(values)*.31));original=digest(initial)
        trial=ref.trial(initial,.0003);assert digest(initial)==original;accepted=ref.commit(initial,trial);rejected=0
        for state in [initial,accepted]:
            try:ref.commit(state,trial)
            except ValueError:rejected+=1
        stale=ref.trial(initial,.0003)
        try:ref.commit(accepted,stale)
        except ValueError:rejected+=1
        continuous=copy.deepcopy(initial);resumed=copy.deepcopy(initial)
        for _ in range(12):continuous=ref.commit(continuous,ref.trial(continuous,.0003))
        for _ in range(5):resumed=ref.commit(resumed,ref.trial(resumed,.0003))
        resumed=ref.load(ref.save(resumed))
        for _ in range(7):resumed=ref.commit(resumed,ref.trial(resumed,.0003))
        assert digest(continuous)==digest(resumed)
        def energy(s):
            x=np.array(s['x']);v=np.array(s['v']);return (x@Kr@x+v@Mr@v)/2
        error=abs(energy(continuous)-energy(initial))/energy(initial);assert error<1e-8
        for key,value in [('model','wrong'),('method','wrong'),('dtype','float32'),('dofs',0),('history_version',2),('history',{'a':[]}),('accepted',False),('step',-1),('cursor',99),('x',[]),('v',[float('inf')]),('time',float('nan'))]:
            bad=copy.deepcopy(initial);bad[key]=value
            try:ref.validate(bad)
            except ValueError:rejected+=1
            else:raise AssertionError('accepted corrupt '+key)
        wrapper=json.loads(ref.save(initial));wrapper['payload']['x'][0]+=.1
        try:ref.load(json.dumps(wrapper))
        except ValueError:rejected+=1
        else:raise AssertionError('checksum')
        assert digest(initial)==original and rejected==16
        records.append(dict(case='same_system_restart_reference',dofs=len(values),all_modes_excited=True,bitwise_replay=True,
            energy_error=error,negative_controls=rejected,product_integrator_adopted=False))
    except Exception as exc:stops.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',stops[-1],flush=True)
    paths=[Path(__file__),Path(m.__file__),Path(m.shell.__file__),Path(m.kernel.__file__),ROOT/'tools/design-experiments/pg_ritz_state_reference.py']
    result=dict(passed=len(records)==6 and not stops,records=records,stops=stops,state_artifact=artifact,seconds=time.monotonic()-start,R4_complete=False,
        sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('MECHANICS_RESULT',result['passed'],flush=True)

if __name__=='__main__':main()
