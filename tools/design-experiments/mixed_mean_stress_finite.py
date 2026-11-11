"""Full finite-thickness MMS-1 strain/stabilization projection diagnostic."""
import time
import numpy as np
from scipy.linalg import solve_triangular
import mixed_mean_stress_element_v2 as f
import reference_work_fit as topology
import moment_position_recovery_action as recovery
m=f.m;qrtools=f.base.qrtools;KINDS=['quartic','quintic','exponential'];NT=3
def reference(q,p,h,k):
    N,_=m.shape(p,0);x=(N@q.reshape(75,3)[:25])[:,0];z=h*p[:,2]/2;a=np.column_stack((np.ones_like(x),np.zeros_like(x),k*x));curvature=np.column_stack((x*x,20*x**3,np.exp(x)));den=1+k*k*x*x-k*z
    tensor=-z[:,None,None,None]*a[:,:,None,None]*a[:,None,:,None]*curvature[:,None,None,:]/den[:,None,None,None]
    E=np.einsum('kij,pijt->pkt',m.MB,tensor);b=np.column_stack((x**4/12,x**5,np.exp(x)-1-x));bp=np.column_stack((x**3/3,5*x**4,np.exp(x)-1));c=np.column_stack((-k*x**5/15,-5*k*x**6/6,-k*((x-1)*np.exp(x)+1-x*x/2)));u=np.zeros((len(p),3,NT));u[:,0]=c-z[:,None]*bp;u[:,2]=b-k*x[:,None]*z[:,None]*bp;return E,u
def strain_rows(s,E):
    _,high,low=f.project(s,E);return np.vstack(((np.einsum('ij,pjt->pit',m.CH,high)*np.sqrt(s['weight'])[:,None,None]).reshape(-1,E.shape[-1]),low))
def rhs(s,q,h,k):
    E,_=reference(q,s['p'],h,k);out=np.zeros((len(s['root']),NT));out[:s['elastic_rows']]=strain_rows(s,E);return out/h**1.5
def get_state(a,cell):
    if a['states'] is not None:return a['states'][cell]
    return f.prepare(a['local'][cell]['q'],a['h'],a['grid']['kappa'])
def assemble(k,n,amp,h):
    original=m.mesh(k,n,amp);registry={};trace_ref=[];fixed=[];local=[];norm=np.zeros(NT);minimum=1.;states=[] if n<=2 else None;light=[]
    for index,ids in enumerate(original['cells']):
        q=original['ref'][ids].reshape(-1);s=f.prepare(q,h,k);nodes,transform=topology.register(s,topology.edge_keys(index,n),registry,trace_ref,fixed);W=s['root']/h**1.5;bb=rhs(s,q,h,k);nc=W.shape[1];ni=3*s['ni'];scale=1/np.linalg.norm(W,axis=0);R=qrtools.rfactor(np.asfortranarray(np.column_stack((W*scale,bb))));T=np.kron(transform,np.eye(3));tail=(R[ni:,ni:nc]/scale[ni:][None,:])@T
        record=dict(q=q,R=R,scale=scale,n=nc,ni=ni,tail=tail,b=R[ni:,nc:nc+NT],trace_transform=T,trace_ids=nodes,reference=np.sum(bb*bb,axis=0),local_ranks=f.ranks(np.linalg.svd(R[:,:nc],compute_uv=False)),local_rcond=qrtools.rcond(R[:ni,:ni]));local.append(record);norm+=record['reference'];minimum=min(minimum,record['local_rcond']);light.append(recovery.light(s))
        if states is not None:states.append(s)
    free=np.setdiff1d(np.arange(len(trace_ref)),fixed);lookup={int(v):i for i,v in enumerate(free)};A=np.zeros((sum(len(c['tail']) for c in local),3*len(free)));b=np.zeros((len(A),NT));offset=0
    for cell in local:
        sl=slice(offset,offset+len(cell['tail']));offset=sl.stop;b[sl]=cell['b']
        for j,node in enumerate(cell['trace_ids']):
            if int(node) in lookup:A[sl,3*lookup[int(node)]:3*lookup[int(node)]+3]=cell['tail'][:,3*j:3*j+3]
    return dict(original=original,grid=dict(kappa=k,n=n,amplitude=amp,trace_ref=np.array(trace_ref),fixed=np.array(fixed)),free=free,A=A,b=b,local=local,reference=norm,h=h,min_local=minimum,states=states,recovery_states=light)
def initial(a):
    A=a['A'];b=a['b'];scale=1/np.linalg.norm(A,axis=0);n=A.shape[1];combined=np.empty((len(A),n+NT),order='F');np.multiply(A,scale[None,:],out=combined[:,:n]);combined[:,n:]=b;a['A']=None;del A
    R=qrtools.rfactor(combined);T=R[:n,:n];rhs0=R[:n,n:n+NT];y=solve_triangular(T,rhs0,check_finite=False);u=scale[:,None]*y;trace=np.zeros((len(a['grid']['trace_ref']),3,NT));trace[a['free']]=u.reshape(-1,3,NT);values=[]
    for cell in a['local']:
        ni,nc=cell['ni'],cell['n'];Rs=cell['R'];ss=cell['scale'];ub=cell['trace_transform']@trace[cell['trace_ids']].reshape(-1,NT);yi=solve_triangular(Rs[:ni,:ni],Rs[:ni,nc:nc+NT]-Rs[:ni,ni:nc]@(ub/ss[ni:,None]),check_finite=False);values.append(np.vstack((ss[:ni,None]*yi,ub)))
    return dict(local=np.array(values),trace=trace,fit=dict(local_rcond=a['min_local'],global_rcond=qrtools.rcond(T),initial_linear_residual=(np.linalg.norm(T@y-rhs0,axis=0)/np.maximum(1e-30,np.linalg.norm(rhs0,axis=0))).tolist())),T,scale,rhs0.copy()
def evaluate(a,sol):
    errors=np.zeros(NT);reference_norm=np.zeros(NT);trace_g=np.zeros_like(sol['trace']);trace_n=np.zeros(sol['trace'].shape[:2]);inside=[];gradients=[];ncols=0
    for index,(cell,values) in enumerate(zip(a['local'],sol['local'])):
        s=get_state(a,index);W=s['root']/a['h']**1.5;norm2=np.sum(W*W,axis=0);star,_=reference(cell['q'],s['p'],a['h'],a['grid']['kappa']);E=s['B']@values-star;S,high,low=f.project(s,E);r_el=strain_rows(s,E)/a['h']**1.5;errors+=np.sum(r_el*r_el,axis=0);b_ref=rhs(s,cell['q'],a['h'],a['grid']['kappa']);reference_norm+=np.sum(b_ref*b_ref,axis=0)
        g=np.einsum('p,pkn,pkt->nt',s['weight'],s['B'],S,optimize=True)
        u=values.reshape(s['nscalar'],3,NT)
        for face in s['faces']:
            jump=np.einsum('pn,nat->pat',face['jump'],u);ww=s['tau']*face['measure'];errors+=np.einsum('p,pat,pat->t',ww,jump,jump)/a['h']**3;g+=np.einsum('p,pn,pat->nat',ww,face['jump'],jump,optimize=True).reshape(cell['n'],NT)
        g/=a['h']**3;gradients.append(g);ni=cell['ni'];inside.append(g[:ni]/np.sqrt(norm2[:ni,None]));ncols+=ni;T=cell['trace_transform'];trace_g[cell['trace_ids']]+=(T.T@g[ni:]).reshape(-1,3,NT);trace_n[cell['trace_ids']]+=np.sum((W[:,ni:]@T)**2,axis=0).reshape(-1,3)
    grad=np.vstack(inside+[trace_g[a['free']].reshape(-1,NT)/np.sqrt(trace_n[a['free']].reshape(-1,1))]);ncols+=3*len(a['free']);kkt=np.linalg.norm(grad,axis=0)/np.maximum(1e-30,np.sqrt(ncols)*np.sqrt(reference_norm));return dict(error=np.sqrt(errors/reference_norm).tolist(),direct_kkt=kkt.tolist()),gradients
def correction(a,gradients,T,scale):
    condensed=np.zeros((len(a['grid']['trace_ref']),3,NT));zs=[];linear=np.zeros(NT)
    def record(A,x,b):
        nonlocal linear
        linear=np.maximum(linear,np.linalg.norm(A@x-b,axis=0)/np.maximum(1e-30,np.linalg.norm(b,axis=0)))
    for cell,g in zip(a['local'],gradients):
        ni,nc=cell['ni'],cell['n'];R=cell['R'];ss=cell['scale'];z=solve_triangular(R[:ni,:ni].T,ss[:ni,None]*g[:ni],lower=True,check_finite=False);zs.append(z);record(np.triu(R[:ni,:ni]).T,z,ss[:ni,None]*g[:ni]);gb=g[ni:]-(R[:ni,ni:nc].T@z)/ss[ni:,None];condensed[cell['trace_ids']]+=(cell['trace_transform'].T@gb).reshape(-1,3,NT)
    b=scale[:,None]*condensed[a['free']].reshape(-1,NT);z=solve_triangular(T.T,b,lower=True,check_finite=False);dy=solve_triangular(T,z,check_finite=False);record(np.triu(T).T,z,b);record(np.triu(T),dy,z);dt=np.zeros_like(condensed);dt[a['free']]=(scale[:,None]*dy).reshape(-1,3,NT);local=[]
    for cell,z in zip(a['local'],zs):
        ni,nc=cell['ni'],cell['n'];R=cell['R'];ss=cell['scale'];db=cell['trace_transform']@dt[cell['trace_ids']].reshape(-1,NT);b=z-R[:ni,ni:nc]@(db/ss[ni:,None]);yi=solve_triangular(R[:ni,:ni],b,check_finite=False);record(np.triu(R[:ni,:ni]),yi,b);local.append(np.vstack((ss[:ni,None]*yi,db)))
    return np.array(local),dt,linear
def solve(a):
    sol,T,scale,rhs0=initial(a);metrics,g=evaluate(a,sol);history=[dict(iteration=0,**metrics)];linear=np.array(sol['fit']['initial_linear_residual']);reason=None;start=time.monotonic()
    for iteration in range(1,5):
        if max(metrics['direct_kkt'])<=1e-8:break
        if time.monotonic()-start>120:reason='correction_time_budget';break
        dl,dt,lr=correction(a,g,T,scale);linear=np.maximum(linear,lr);accepted=False
        for bt in range(4):
            alpha=2.**(-bt);trial=dict(local=sol['local']-alpha*dl,trace=sol['trace']-alpha*dt);tm,tg=evaluate(a,trial)
            if max(tm['direct_kkt'])<max(metrics['direct_kkt']):sol['local']=trial['local'];sol['trace']=trial['trace'];metrics,g=tm,tg;history.append(dict(iteration=iteration,alpha=alpha,**metrics));accepted=True;break
        if not accepted:reason='original_stationarity_not_decreasing';break
    if max(metrics['direct_kkt'])>1e-8 and reason is None:reason='correction_iteration_limit'
    gradient=np.zeros_like(sol['trace']);error=np.zeros(NT)
    for cell,u in zip(a['local'],sol['local']):
        n=cell['n'];rr=cell['R'][:,:n]@(u/cell['scale'][:,None])-cell['R'][:,n:n+NT];error+=np.sum(rr*rr,axis=0);canonical=sol['trace'][cell['trace_ids']].reshape(-1,NT);gradient[cell['trace_ids']]+=(cell['tail'].T@(cell['tail']@canonical-cell['b'])).reshape(-1,3,NT)
    sol['fit'].update(metrics,compressed_error=np.sqrt(error/a['reference']).tolist(),compressed_kkt=(np.linalg.norm(scale[:,None]*gradient[a['free']].reshape(-1,NT),axis=0)/np.maximum(1e-30,np.sqrt(len(scale))*np.linalg.norm(a['b'],axis=0))).tolist(),linear_residual=linear.tolist(),pre_correction_system_residual_at_final_solution=(np.linalg.norm(T@(sol['trace'][a['free']].reshape(-1,NT)/scale[:,None])-rhs0,axis=0)/np.maximum(1e-30,np.linalg.norm(rhs0,axis=0))).tolist(),history=history,stop_reason=reason);return sol
def observe(a,sol,order=15,zorder=9):
    total=np.zeros(NT);refnorm=np.zeros(NT);states=[];values=[]
    for index,(cell,value) in enumerate(zip(a['local'],sol['local'])):
        s=f.prepare(cell['q'],a['h'],a['grid']['kappa'],order,zorder,True);s0=get_state(a,index);T=np.eye(s['nscalar'])
        for fm,fi in zip(s0['faces'],s['faces']):
            sl=slice(fi['offset'],fi['offset']+fi['n']);T[sl,sl]=fi['Mcoeff'].T@fm['Mcoeff']
        u=np.kron(T,np.eye(3))@value;star,_=reference(cell['q'],s['p'],a['h'],a['grid']['kappa']);res=strain_rows(s,s['B']@u-star);total+=np.sum(res*res,axis=0);rr=rhs(s,cell['q'],a['h'],a['grid']['kappa'])*a['h']**1.5;refnorm+=np.sum(rr*rr,axis=0)
        for face in s['faces']:
            jump=np.einsum('pn,nat->pat',face['jump'],u.reshape(s['nscalar'],3,NT));total+=np.einsum('p,pat,pat->t',s['tau']*face['measure'],jump,jump)
        states.append(recovery.light(s));values.append(u.reshape(s['nscalar'],3*NT))
    R=recovery.RecoveryAction(states,a['grid']['n'],values,order,zorder,clamp=True);position=np.zeros(NT);norm=np.zeros(NT)
    for index,cell in enumerate(a['local']):
        p,w=f.base.volume_points(order,zorder);_,_,det,_,_=f.base.geometry(cell['q'],a['h'],p,True);_,target=reference(cell['q'],p,a['h'],a['grid']['kappa']);v=R.evaluate(index,p).reshape(len(p),3,NT);position+=np.einsum('p,pat,pat->t',w*det,v-target,v-target);norm+=np.einsum('p,pat,pat->t',w*det,target,target)
    return dict(error=np.sqrt(total/refnorm).tolist(),recovered_position_error=np.sqrt(position/norm).tolist(),recovery_global_identity_allocated=hasattr(R,'eye'))
