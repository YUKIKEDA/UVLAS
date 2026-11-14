"""RW-3: QR-preconditioned defect correction on the original physical residual.

No FE space, stabilization or objective change. Triangular QR factors apply the
block inverse Hessian; no normal matrix is assembled. Corrections are accepted
only if the uncondensed scaled stationarity norm decreases.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import time
import numpy as np
from scipy.linalg import solve_triangular
import reference_work_fit as base
from reference_work_fit import f,m,KIND_TARGETS,target,assemble,observe,register,edge_keys

def initial(a):
    A=a['A'];b=a['b'];scale=1/np.linalg.norm(A,axis=0);n=A.shape[1];combined=np.empty((len(A),n+2),order='F');np.multiply(A,scale[None,:],out=combined[:,:n]);combined[:,n:]=b;a['A']=None;del A
    R=f.qrtools.rfactor(combined);del combined;T=R[:n,:n];rhs=R[:n,n:n+2];y=solve_triangular(T,rhs,check_finite=False);u=scale[:,None]*y
    trace=np.zeros((len(a['grid']['trace_ref']),3,2));trace[a['free']]=u.reshape(-1,3,2);local=[];error=np.zeros(2);gradient=np.zeros_like(trace)
    for cell in a['local']:
        ni,nlocal=cell['ni'],cell['n'];Rl=cell['R'];s=cell['scale'];canonical=trace[cell['trace_ids']].reshape(-1,2);ub=cell['trace_transform']@canonical;yb=ub/s[ni:,None]
        yi=solve_triangular(Rl[:ni,:ni],Rl[:ni,nlocal:nlocal+2]-Rl[:ni,ni:nlocal]@yb,check_finite=False);value=np.vstack((s[:ni,None]*yi,ub));local.append(value)
        res=Rl[:,:nlocal]@(value/s[:,None])-Rl[:,nlocal:nlocal+2];error+=np.sum(res*res,axis=0);gradient[cell['trace_ids']]+=(cell['tail'].T@(cell['tail']@canonical-cell['b'])).reshape(-1,3,2)
    kkt=np.linalg.norm(scale[:,None]*gradient[a['free']].reshape(-1,2),axis=0)/np.maximum(1e-30,np.sqrt(n)*np.linalg.norm(b,axis=0))
    linear=np.linalg.norm(T@y-rhs,axis=0)/np.maximum(1e-30,np.linalg.norm(rhs,axis=0))
    return dict(trace=trace,local=np.array(local),fit=dict(compressed_error=np.sqrt(error/a['reference']).tolist(),kkt=kkt.tolist(),linear_residual=linear.tolist(),local_rcond=a['min_local'],global_rcond=f.qrtools.rcond(T))),T,scale,rhs.copy()

def evaluate(a,sol):
    elastic=np.zeros(2);jump=np.zeros(2);ref=np.zeros(2);trace_g=np.zeros_like(sol['trace']);trace_n=np.zeros(sol['trace'].shape[:2]);inside=[];gradients=[];ncols=0
    for cell,values in zip(a['local'],sol['local']):
        s=f.prepare(cell['q'],a['h']);W=s['root']/a['h']**1.5;norm2=np.sum(W*W,axis=0);residuals=[]
        for j,kind in enumerate(KIND_TARGETS):
            star=target(cell['q'],s['p'],a['h'],a['grid']['kappa'],kind);E=s['B']@values[:,j]-star;ep=(E@m.CH.T*np.sqrt(s['weight'])[:,None]).reshape(-1)/a['h']**1.5;elastic[j]+=ep@ep;jr=[]
            for face in s['faces']:
                val=face['jump']@values[:,j].reshape(-1,3);jr.append((val*np.sqrt(s['tau']*face['measure'])[:,None]).reshape(-1)/a['h']**1.5)
            jp=np.concatenate(jr);jump[j]+=jp@jp;residuals.append(np.r_[ep,jp]);r=(star@m.CH.T*np.sqrt(s['weight'])[:,None]).reshape(-1)/a['h']**1.5;ref[j]+=r@r
        g=W.T@np.column_stack(residuals);gradients.append(g);ni=cell['ni'];inside.append(g[:ni]/np.sqrt(norm2[:ni,None]));ncols+=ni
        T=cell['trace_transform']
        if np.linalg.norm(T-np.eye(len(T)))>1e-10:raise ValueError('direct_stationarity_requires_canonical_benchmark_faces')
        trace_g[cell['trace_ids']]+=(T.T@g[ni:]).reshape(-1,3,2);trace_n[cell['trace_ids']]+=norm2[ni:].reshape(-1,3)
    grad=np.vstack(inside+[trace_g[a['free']].reshape(-1,2)/np.sqrt(trace_n[a['free']].reshape(-1,1))]);ncols+=3*len(a['free']);kkt=np.linalg.norm(grad,axis=0)/np.maximum(1e-30,np.sqrt(ncols)*np.sqrt(ref))
    return dict(error=np.sqrt((elastic+jump)/ref).tolist(),direct_kkt=kkt.tolist()),gradients

def correction(a,gradients,T,scale):
    condensed=np.zeros((len(a['grid']['trace_ref']),3,2));zs=[];solve_residual=np.zeros(2)
    def record(A,x,b):
        nonlocal solve_residual
        solve_residual=np.maximum(solve_residual,np.linalg.norm(A@x-b,axis=0)/np.maximum(1e-30,np.linalg.norm(b,axis=0)))
    for cell,g in zip(a['local'],gradients):
        ni=cell['ni'];n=cell['n'];R=cell['R'];s=cell['scale'];z=solve_triangular(R[:ni,:ni].T,s[:ni,None]*g[:ni],lower=True,check_finite=False);zs.append(z)
        record(np.triu(R[:ni,:ni]).T,z,s[:ni,None]*g[:ni])
        gb=g[ni:]-(R[:ni,ni:n].T@z)/s[ni:,None];condensed[cell['trace_ids']]+=(cell['trace_transform'].T@gb).reshape(-1,3,2)
    rhs=scale[:,None]*condensed[a['free']].reshape(-1,2);z=solve_triangular(T.T,rhs,lower=True,check_finite=False);dy=solve_triangular(T,z,check_finite=False)
    record(np.triu(T).T,z,rhs);record(np.triu(T),dy,z)
    dt=np.zeros_like(condensed);dt[a['free']]=(scale[:,None]*dy).reshape(-1,3,2);local=[]
    for cell,z in zip(a['local'],zs):
        ni,n=cell['ni'],cell['n'];R=cell['R'];s=cell['scale'];db=cell['trace_transform']@dt[cell['trace_ids']].reshape(-1,2)
        rhs=z-R[:ni,ni:n]@(db/s[ni:,None]);yi=solve_triangular(R[:ni,:ni],rhs,check_finite=False);record(np.triu(R[:ni,:ni]),yi,rhs)
        di=s[:ni,None]*yi;local.append(np.vstack((di,db)))
    return np.array(local),dt,solve_residual

def solve(a):
    start=time.monotonic();sol,T,scale,rhs=initial(a);metrics,grad=evaluate(a,sol);history=[dict(iteration=0,**metrics)];reason=None;solve_residual=np.array(sol['fit']['linear_residual'])
    for iteration in range(1,5):
        if max(metrics['direct_kkt'])<=1e-8:break
        if time.monotonic()-start>120:reason='correction_time_budget';break
        dl,dt,linear=correction(a,grad,T,scale);solve_residual=np.maximum(solve_residual,linear);accepted=False
        for backtrack in range(4):
            alpha=2.**(-backtrack);trial=dict(local=sol['local']-alpha*dl,trace=sol['trace']-alpha*dt);tm,tg=evaluate(a,trial)
            if max(tm['direct_kkt'])<max(metrics['direct_kkt']):
                sol['local']=trial['local'];sol['trace']=trial['trace'];metrics,grad=tm,tg;history.append(dict(iteration=iteration,alpha=alpha,correction_linear_residual=linear.tolist(),**metrics));accepted=True;break
        if not accepted:reason='original_stationarity_not_decreasing';break
    if max(metrics['direct_kkt'])>1e-8 and reason is None:reason='correction_iteration_limit'
    sol['fit']['initial_compressed_kkt']=sol['fit']['kkt'];sol['fit']['initial_linear_residual']=sol['fit']['linear_residual']
    gradient=np.zeros_like(sol['trace']);error=np.zeros(2)
    for cell,values in zip(a['local'],sol['local']):
        n=cell['n'];rr=cell['R'][:,:n]@(values/cell['scale'][:,None])-cell['R'][:,n:n+2];error+=np.sum(rr*rr,axis=0)
        u=sol['trace'][cell['trace_ids']].reshape(-1,2);gradient[cell['trace_ids']]+=(cell['tail'].T@(cell['tail']@u-cell['b'])).reshape(-1,3,2)
    final_y=sol['trace'][a['free']].reshape(-1,2)/scale[:,None]
    sol['fit']['compressed_error']=np.sqrt(error/a['reference']).tolist()
    sol['fit']['kkt']=(np.linalg.norm(scale[:,None]*gradient[a['free']].reshape(-1,2),axis=0)/np.maximum(1e-30,np.sqrt(len(scale))*np.linalg.norm(a['b'],axis=0))).tolist()
    sol['fit']['pre_correction_system_residual_at_final_solution']=(np.linalg.norm(T@final_y-rhs,axis=0)/np.maximum(1e-30,np.linalg.norm(rhs,axis=0))).tolist()
    sol['fit']['linear_residual']=solve_residual.tolist()
    sol['fit'].update(metrics);sol['fit'].update(correction_history=history,correction_stop_reason=reason,objective_evaluation='original_primary_physical_strain_and_stabilization',FE_operator_changed=False)
    return sol
