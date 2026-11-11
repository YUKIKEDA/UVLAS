"""RW-3 square-root condensation of volume plus top/bottom face positions."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.linalg import solve_triangular
import reference_work_ancf_v3 as f
import polynomial_trace_fit as prior
m=f.m;KIND_TARGETS=prior.KIND_TARGETS;target=prior.target

def edge_keys(cell,n):
    j,i=divmod(cell,n);return [(0,i,j),(0,i+1,j),(1,j,i),(1,j+1,i)]

def register(state,keys,registry,trace_ref,fixed):
    indices=[];blocks=[]
    for face,key in zip(state['faces'][2:],keys):
        Mc=face['Mcoeff'];sl=slice(face['offset'],face['offset']+face['n'])
        if key not in registry:
            offset=len(trace_ref);registry[key]=(offset,Mc.copy());trace_ref.extend(state['X0'][sl])
            if key[0]==0 and key[1]==0:fixed.extend(range(offset,offset+face['n']))
        offset,canonical=registry[key]
        if canonical.shape!=Mc.shape:raise ValueError('shared_face_dimension')
        T=Mc.T@canonical
        if max(np.linalg.norm(T.T@T-np.eye(len(T))),np.linalg.norm(Mc@T-canonical))>1e-10:raise ValueError('shared_face_space')
        indices.extend(range(offset,offset+face['n']));blocks.append(T)
    nt=sum(len(x) for x in blocks);transform=np.zeros((nt,nt));offset=0
    for block in blocks:
        n=len(block);transform[offset:offset+n,offset:offset+n]=block;offset+=n
    return np.array(indices),transform

def make_local(q,h,k,state,transform):
    W=state['root']/h**1.5;n=W.shape[1];ni=3*state['ni'];targets=np.zeros((len(W),2));p=state['p'];weight=state['weight']
    for j,kind in enumerate(KIND_TARGETS):targets[:state['elastic_rows'],j]=(target(q,p,h,k,kind)@m.CH.T*np.sqrt(weight)[:,None]).reshape(-1)/h**1.5
    scale=1/np.linalg.norm(W,axis=0);R=f.qrtools.rfactor(np.column_stack((W*scale,targets)))
    T=np.kron(transform,np.eye(3));tail=(R[ni:,ni:n]/scale[ni:][None,:])@T
    return dict(q=q,R=R,scale=scale,n=n,ni=ni,tail=tail,b=R[ni:,n:n+2],trace_transform=T,local_ranks=f.ranks(np.linalg.svd(R[:,:n],compute_uv=False)),
        reference=np.sum(targets*targets,axis=0),local_rcond=f.qrtools.rcond(R[:ni,:ni]))

def assemble(k,n,amp,h,kind='RW3'):
    original=m.mesh(k,n,amp);registry={};trace_ref=[];fixed=[];local=[];norm=np.zeros(2);minimum=1.
    for cell,ids in enumerate(original['cells']):
        q=original['ref'][ids].reshape(-1);state=f.prepare(q,h,kind);nodes,transform=register(state,edge_keys(cell,n),registry,trace_ref,fixed)
        record=make_local(q,h,k,state,transform);record['trace_ids']=nodes;local.append(record);norm+=record['reference'];minimum=min(minimum,record['local_rcond']);del state
    free=np.setdiff1d(np.arange(len(trace_ref)),fixed);lookup={int(v):i for i,v in enumerate(free)}
    A=np.zeros((sum(len(c['tail']) for c in local),3*len(free)));b=np.zeros((len(A),2));offset=0
    for cell in local:
        sl=slice(offset,offset+len(cell['tail']));offset=sl.stop;b[sl]=cell['b']
        for j,node in enumerate(cell['trace_ids']):
            if int(node) in lookup:A[sl,3*lookup[int(node)]:3*lookup[int(node)]+3]=cell['tail'][:,3*j:3*j+3]
    return dict(original=original,grid=dict(kappa=k,n=n,amplitude=amp,trace_ref=np.array(trace_ref),fixed=np.array(fixed)),free=free,A=A,b=b,local=local,reference=norm,h=h,kind=kind,min_local=minimum)

def solve(a):
    A=a['A'];b=a['b'];scale=1/np.linalg.norm(A,axis=0);n=A.shape[1];combined=np.empty((len(A),n+2),order='F');np.multiply(A,scale[None,:],out=combined[:,:n]);combined[:,n:]=b;a['A']=None;del A
    R=f.qrtools.rfactor(combined);del combined;T=R[:n,:n];rhs=R[:n,n:n+2];y=solve_triangular(T,rhs,check_finite=False);u=scale[:,None]*y
    trace=np.zeros((len(a['grid']['trace_ref']),3,2));trace[a['free']]=u.reshape(-1,3,2);local=[];error=np.zeros(2);gradient=np.zeros_like(trace)
    for cell in a['local']:
        ni,nlocal=cell['ni'],cell['n'];Rl=cell['R'];s=cell['scale'];canonical=trace[cell['trace_ids']].reshape(-1,2);ub=cell['trace_transform']@canonical;yb=ub/s[ni:,None]
        yi=solve_triangular(Rl[:ni,:ni],Rl[:ni,nlocal:nlocal+2]-Rl[:ni,ni:nlocal]@yb,check_finite=False);value=np.vstack((s[:ni,None]*yi,ub))
        residual=Rl[:,:nlocal]@(value/s[:,None])-Rl[:,nlocal:nlocal+2];error+=np.sum(residual*residual,axis=0);local.append(value)
        condensed=cell['tail']@canonical-cell['b'];gradient[cell['trace_ids']]+=(cell['tail'].T@condensed).reshape(-1,3,2)
    kkt=np.linalg.norm(scale[:,None]*gradient[a['free']].reshape(-1,2),axis=0)/np.maximum(1e-30,np.sqrt(n)*np.linalg.norm(b,axis=0));linear=np.linalg.norm(T@y-rhs,axis=0)/np.maximum(1e-30,np.linalg.norm(rhs,axis=0))
    return dict(trace=trace,local=np.array(local),fit=dict(error=np.sqrt(error/a['reference']).tolist(),kkt=kkt.tolist(),linear_residual=linear.tolist(),local_rcond=a['min_local'],global_rcond=f.qrtools.rcond(T)))

def observe(a,solution,order=15,zorder=9):
    elastic=np.zeros(2);jump=np.zeros(2);ref=np.zeros(2);H=np.zeros(2);area=0.
    for cell,u in zip(a['local'],solution['local']):
        q=cell['q'];s=f.prepare(q,a['h'],order=order,zorder=zorder,independent=True)
        # Canonical coefficient spaces have the same orientation in this paired
        # evaluation; explicitly transform each facet for general geometries.
        sm=f.prepare(q,a['h']);transform=np.eye(s['nscalar'])
        for fm,fi in zip(sm['faces'],s['faces']):
            sl=slice(fi['offset'],fi['offset']+fi['n']);transform[sl,sl]=fi['Mcoeff'].T@fm['Mcoeff']
        ui=(np.kron(transform,np.eye(3))@u)
        for j,kind in enumerate(KIND_TARGETS):
            star=target(q,s['p'],a['h'],a['grid']['kappa'],kind);res=s['B']@ui[:,j]-star
            elastic[j]+=s['weight']@np.einsum('pk,kl,pl->p',res,m.C,res);ref[j]+=s['weight']@np.einsum('pk,kl,pl->p',star,m.C,star)
            for face in s['faces']:jump[j]+=s['tau']*face['measure']@np.sum((face['jump']@ui[:,j].reshape(-1,3))**2,axis=1)
            nodal=f.volume_nodal(s,ui[:225,j].reshape(75,3));hi,ar=prior.h1(q,nodal.reshape(-1),order);H[j]+=hi
            if j==0:area+=ar
    return dict(error=np.sqrt((elastic+jump)/ref).tolist(),elastic_error=np.sqrt(elastic/ref).tolist(),stabilization_error=np.sqrt(jump/ref).tolist(),h1=np.sqrt(H/area).tolist())
