"""Independent uncondensed QR with a rebuilt global right-preconditioned operator."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import solve_triangular
import mixed_mean_stress_finite as fit
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2]

def solve_full(a):
    total=sum(c['ni'] for c in a['local'])
    lookup={int(v):i for i,v in enumerate(a['free'])}
    n=total+3*len(a['free']);maps=[];offset=0;counts=[]
    for index,cell in enumerate(a['local']):
        indices=list(range(offset,offset+cell['ni']));offset+=cell['ni']
        indices.extend(total+3*lookup[int(node)]+cart if int(node) in lookup else -1 for node in cell['trace_ids'] for cart in range(3))
        maps.append(np.array(indices));counts.append(len(fit.get_state(a,index)['root']))
    def build():
        A=np.zeros((sum(counts),n+fit.NT),order='F');offset=0
        for index,(cell,indices,count) in enumerate(zip(a['local'],maps,counts)):
            s=fit.get_state(a,index);W=s['root']/a['h']**1.5;ni=cell['ni']
            Wc=np.column_stack((W[:,:ni],W[:,ni:]@cell['trace_transform']))
            valid=indices>=0;rows=slice(offset,offset+count);offset+=count
            A[rows,indices[valid]]=Wc[:,valid]
            A[rows,n:]=fit.rhs(s,cell['q'],a['h'],a['grid']['kappa'])
        return A
    A=build();scale=1/np.linalg.norm(A[:,:n],axis=0)
    np.multiply(A[:,:n],scale[None,:],out=A[:,:n]);R0=fit.qrtools.rfactor(A);del A
    T0=R0[:n,:n].copy();del R0
    raw_rcond=fit.qrtools.rcond(T0)
    A=build()
    for start in range(0,len(A),256):
        block=A[start:start+256,:n]*scale[None,:]
        A[start:start+256,:n]=solve_triangular(T0.T,block.T,lower=True,check_finite=False).T
    R=fit.qrtools.rfactor(A);del A
    T=R[:n,:n];rhs=R[:n,n:n+fit.NT]
    map_relative=np.zeros(fit.NT);map_backward=np.zeros(fit.NT)
    normT0=np.linalg.norm(T0)
    def mapping(b,transpose=False):
        nonlocal map_relative,map_backward
        M=T0.T if transpose else T0
        x=solve_triangular(M,b,lower=transpose,check_finite=False)
        residual=np.linalg.norm(M@x-b,axis=0)
        map_relative=np.maximum(map_relative,residual/np.maximum(1e-30,np.linalg.norm(b,axis=0)))
        map_backward=np.maximum(map_backward,residual/np.maximum(1e-30,normT0*np.linalg.norm(x,axis=0)+np.linalg.norm(b,axis=0)))
        return x
    def unpack(x):
        local=[]
        for cell,indices in zip(a['local'],maps):
            u=np.zeros((len(indices),fit.NT));valid=indices>=0;u[valid]=x[indices[valid]]
            u[cell['ni']:]=cell['trace_transform']@u[cell['ni']:];local.append(u)
        trace=np.zeros((len(a['grid']['trace_ref']),3,fit.NT))
        trace[a['free']]=x[total:].reshape(-1,3,fit.NT)
        return dict(local=np.array(local),trace=trace)
    y=solve_triangular(T,rhs,check_finite=False);value=scale[:,None]*mapping(y)
    linear=np.linalg.norm(T@y-rhs,axis=0)/np.maximum(1e-30,np.linalg.norm(rhs,axis=0))
    sol=unpack(value);metrics,grad=fit.evaluate(a,sol);history=[dict(iteration=0,**metrics)];reason=None
    for iteration in range(1,5):
        if max(metrics['direct_kkt'])<=1e-8:break
        g=np.zeros((n,fit.NT))
        for cell,indices,gi in zip(a['local'],maps,grad):
            valid=indices>=0;canonical=np.vstack((gi[:cell['ni']],cell['trace_transform'].T@gi[cell['ni']:]))
            np.add.at(g,indices[valid],canonical[valid])
        b=mapping(scale[:,None]*g,True)
        z=solve_triangular(T.T,b,lower=True,check_finite=False);dy=solve_triangular(T,z,check_finite=False)
        linear=np.maximum(linear,np.linalg.norm(T.T@z-b,axis=0)/np.maximum(1e-30,np.linalg.norm(b,axis=0)))
        linear=np.maximum(linear,np.linalg.norm(T@dy-z,axis=0)/np.maximum(1e-30,np.linalg.norm(z,axis=0)))
        delta=scale[:,None]*mapping(dy);accepted=False
        for bt in range(4):
            alpha=2.**(-bt);trial_value=value-alpha*delta;trial=unpack(trial_value);tm,tg=fit.evaluate(a,trial)
            if max(tm['direct_kkt'])<max(metrics['direct_kkt']):
                value=trial_value;sol=trial;metrics,grad=tm,tg;history.append(dict(iteration=iteration,alpha=alpha,**metrics));accepted=True;break
        if not accepted:reason='original_stationarity_not_decreasing';break
    if max(metrics['direct_kkt'])>1e-8 and reason is None:reason='correction_iteration_limit'
    return sol,dict(**metrics,raw_full_rcond=raw_rcond,preconditioned_full_rcond=fit.qrtools.rcond(T),linear_residual=linear.tolist(),basis_map_rhs_relative_residual=map_relative.tolist(),basis_map_backward_error=map_backward.tolist(),history=history,stop_reason=reason,global_coordinates=n,local_condensation_factors_used=False,rebuilt_original_operator=True)

def main():
    source=ROOT/'docs/verification/mixed-mean-stress-finite-pair-result.json';parent=json.loads(source.read_text());assert parent['all_checks_pass']
    dest=ROOT/'docs/verification/mixed-mean-stress-finite-full-v2-result.json';assert not dest.exists()
    start=time.monotonic();records=[];stops=[];out=ROOT/'build/analysis'/('mixed-mean-stress-finite-full-v2-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    for r in parent['records']:
        try:
            tick=time.monotonic();a=fit.assemble(r['kappa'],r['n'],r['amplitude'],r['h']);sol,metrics=solve_full(a);obs=fit.observe(a,sol)
            old=np.load(ROOT/r['artifact']);field=np.linalg.norm(sol['local']-old['local'],axis=(0,1))/np.linalg.norm(old['local'],axis=(0,1))
            error=max(abs(x-y) for x,y in zip(obs['error'],r['observed']['error']))
            position=max(abs(x-y) for x,y in zip(obs['recovered_position_error'],r['observed']['recovered_position_error']))
            checks=dict(original_stationarity=max(metrics['direct_kkt'])<=1e-8,linear=max(metrics['linear_residual'])<=1e-8,preconditioned_condition=metrics['preconditioned_full_rcond']>=1e-10,basis_map_backward_error=max(metrics['basis_map_backward_error'])<=1e-12,fields=max(field)<=1e-8,error=error<=1e-8,position=position<=1e-8,resources=time.monotonic()-tick<=180 and memory.memory_bytes()<=2*1024**3)
            path=out/f"h{r['h']:g}.npz";np.savez_compressed(path,local=sol['local'],trace=sol['trace'])
            records.append(dict(h=r['h'],metrics=metrics,observed=obs,field_difference=field.tolist(),error_difference=error,position_difference=position,checks={k:bool(v) for k,v in checks.items()},seconds=time.monotonic()-tick,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            print('FULL_V2',r['h'],records[-1],flush=True)
            if not all(checks.values()):raise RuntimeError('MMS1_full_v2_reference_gate')
            if time.monotonic()-start>900:raise RuntimeError('MMS1_full_v2_reference_series_resources')
        except Exception as exc:stops.append(dict(h=r['h'],reason=repr(exc)));print('STOP',repr(exc),flush=True);break
    paths=[Path(__file__),Path(fit.__file__),source,ROOT/'docs/verification/mixed-mean-stress-finite-full-result.json',ROOT/'docs/verification/mixed-mean-stress-finite-full-v2-protocol.md']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_checks_pass=len(records)==2 and not stops,independent_global_resolve=True,main_shell_adopted=False,physical_PDE_acceptance=False,raw_basis_condition_failure_resolved=False,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
