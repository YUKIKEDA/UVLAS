"""Square-root condensation: volume ANCF variables to global position trace."""
import numpy as np
from scipy.linalg import solve_triangular
import hybrid_gradient_ancf as f
m=f.m;KIND_TARGETS=['quartic','exponential']
INS=np.array([25*a+5*j+i for a in range(3) for j in range(1,4) for i in range(1,4)])
BND=np.array([i for i in range(75) if i not in INS]);BASE_PERM=(3*np.r_[INS,BND][:,None]+np.arange(3)).reshape(-1)

def target(q,p,h,k,kind):
    N,_=m.shape(p,0);x=(N@q.reshape(75,3)[:25])[:,0];z=h*p[:,2]/2
    v=x*x if kind=='quartic' else np.exp(x);a=np.column_stack((np.ones_like(x),np.zeros_like(x),k*x));den=1+k*k*x*x-k*z
    return np.einsum('kij,pij->pk',m.MB,-z[:,None,None]*v[:,None,None]/den[:,None,None]*a[:,:,None]*a[:,None,:])

def h1(q,u,order=15):
    p,w=m.quadrature(order,2);N,g=m.shape(p,0,True);J=np.einsum('na,pnb->pab',q.reshape(75,3)[:25],g);G=J.transpose(0,2,1)@J
    weight=w*np.sqrt(np.linalg.det(G));total=0.
    for field in [1,2]:
        ns,ds=m.shape(p,field,True);values=ns@u.reshape(75,3)[m.SLICES[field]];du=np.einsum('na,pnb->pab',u.reshape(75,3)[m.SLICES[field]],ds)
        total+=weight@(np.sum(values*values,axis=1)+np.einsum('pai,pij,paj->p',du,np.linalg.inv(G),du))
    return float(total),float(weight.sum())

def make_local(q,h,k,kind,baseline=False):
    if baseline:
        p,w=f.volume_points(13,7);_,J,B=m.raw(q,p,h*p[:,2]/2);weight=w*np.linalg.det(J)*h/2
        W=(np.einsum('ij,pjn->pin',m.CH,B)*np.sqrt(weight)[:,None,None]).reshape(-1,225)
        perm=BASE_PERM;ni=81;W=W[:,perm]
    else:
        state=f.prepare(q,h,kind);p=state['p'];weight=state['weight'];W=state['root'];perm=np.arange(405);ni=225
    norm=h**1.5;W=W/norm;n=W.shape[1];targets=np.zeros((len(W),2))
    for j,ref in enumerate(KIND_TARGETS):
        E=target(q,p,h,k,ref);b=(E@m.CH.T*np.sqrt(weight)[:,None]).reshape(-1)/norm;targets[:len(b),j]=b
    scale=1/np.linalg.norm(W,axis=0);R=f.qrtools.rfactor(np.column_stack((W*scale,targets)))
    local_ranks=f.ranks(np.linalg.svd(R[:,:n],compute_uv=False))
    return dict(q=q,R=R,scale=scale,perm=perm,n=n,ni=ni,tail=R[ni:,ni:n]/scale[ni:][None,:],b=R[ni:,n:n+2],local_ranks=local_ranks,
        reference=np.sum(targets*targets,axis=0),local_rcond=f.qrtools.rcond(R[:ni,:ni]))

def assemble(k,n,amp,h,kind,baseline=False):
    grid=f.mesh(k,n,amp);original=grid['original'];local=[]
    if baseline:
        interior=set(int(node) for ids in original['cells'] for node in ids[INS]);free=np.array([node for node in np.concatenate(original['free']) if int(node) not in interior])
        trace_maps=[ids[BND] for ids in original['cells']]
    else:free=grid['free'];trace_maps=grid['connect']
    lookup={int(node):i for i,node in enumerate(free)};rows=(146 if baseline else 182)*n*n;A=np.zeros((rows,3*len(free)));b=np.zeros((rows,2));norm=np.zeros(2);minimum=1.;offset=0
    for c,ids in enumerate(original['cells']):
        q=original['ref'][ids].reshape(-1);cell=make_local(q,h,k,kind,baseline);count=len(cell['tail']);sl=slice(offset,offset+count);offset+=count
        b[sl]=cell['b'];nodes=trace_maps[c]
        for j,node in enumerate(nodes):
            if int(node) in lookup:A[sl,3*lookup[int(node)]:3*lookup[int(node)]+3]=cell['tail'][:,3*j:3*j+3]
        cell['trace_ids']=nodes;local.append(cell);norm+=cell['reference'];minimum=min(minimum,cell['local_rcond'])
    return dict(grid=grid,original=original,free=free,A=A,b=b,local=local,reference=norm,h=h,kind=kind,baseline=baseline,min_local=minimum)

def solve(a):
    A,b=a['A'],a['b'];scale=1/np.linalg.norm(A,axis=0);R=f.qrtools.rfactor(np.column_stack((A*scale,b)));n=A.shape[1]
    T=R[:n,:n];y=solve_triangular(T,R[:n,n:n+2],check_finite=False);u=scale[:,None]*y
    trace_nodes=len(a['original']['ref']) if a['baseline'] else len(a['grid']['trace_ref']);trace=np.zeros((trace_nodes,3,2));trace[a['free']]=u.reshape(-1,3,2)
    local_solutions=[];error=np.zeros(2)
    for cell in a['local']:
        ni,nlocal=cell['ni'],cell['n'];Rl=cell['R'];s=cell['scale'];ub=trace[cell['trace_ids']].reshape(-1,2);yb=ub/s[ni:,None]
        yi=solve_triangular(Rl[:ni,:ni],Rl[:ni,nlocal:nlocal+2]-Rl[:ni,ni:nlocal]@yb,check_finite=False)
        value=np.vstack((s[:ni,None]*yi,ub));value_original=np.zeros_like(value);value_original[cell['perm']]=value
        residual=Rl[:,:nlocal]@(value/s[:,None])-Rl[:,nlocal:nlocal+2];error+=np.sum(residual*residual,axis=0);local_solutions.append(value_original)
    residual=A@u-b;scaledA=A*scale
    kkt=np.linalg.norm(scaledA.T@residual,axis=0)/np.maximum(1e-30,np.linalg.norm(scaledA)*np.linalg.norm(b,axis=0))
    linear=np.linalg.norm(T@y-R[:n,n:n+2],axis=0)/np.maximum(1e-30,np.linalg.norm(R[:n,n:n+2],axis=0))
    return dict(trace=trace,local=np.array(local_solutions),fit=dict(error=np.sqrt(error/a['reference']).tolist(),kkt=kkt.tolist(),linear_residual=linear.tolist(),
        local_rcond=a['min_local'],global_rcond=f.qrtools.rcond(T)))

def observe(a,solution,order=15,zorder=9):
    error=np.zeros(2);elastic=np.zeros(2);jump=np.zeros(2);reference=np.zeros(2);H=np.zeros(2);area=0.
    for c,cell in enumerate(a['local']):
        q=cell['q'];h=a['h'];values=solution['local'][c]
        if a['baseline']:
            p,w=f.volume_points(order,zorder);_,J,B=m.raw(q,p,h*p[:,2]/2,True);weight=w*np.linalg.det(J)*h/2;state=None
        else:state=f.prepare(q,h,a['kind'],order,zorder,True);p=state['p'];B=state['B'];weight=state['weight']
        for j,kind in enumerate(KIND_TARGETS):
            star=target(q,p,h,a['grid']['kappa'],kind);res=B@values[:,j]-star
            energy=np.einsum('pk,kl,pl->p',res,m.C,res);ref=np.einsum('pk,kl,pl->p',star,m.C,star);elastic[j]+=weight@energy;reference[j]+=weight@ref
            if state is not None:
                u=values[:,j].reshape(135,3)
                for face in state['faces']:jump[j]+=state['tau']*face['measure']@np.sum((face['jump']@u)**2,axis=1)
            hi,ar=h1(q,values[:225,j],order);H[j]+=hi
            if j==0:area+=ar
    error=elastic+jump
    return dict(error=np.sqrt(error/reference).tolist(),elastic_error=np.sqrt(elastic/reference).tolist(),stabilization_error=np.sqrt(jump/reference).tolist(),h1=np.sqrt(H/area).tolist())

def rigid_graph(grid):
    nc=len(grid['original']['cells']);rows=[];seen={};supports=[]
    for c,ids in enumerate(grid['original']['cells']):
        X=f.reference_state(grid['original']['ref'][ids].reshape(-1));R=np.zeros((405,6));ridx=np.r_[np.arange(25),*[np.arange(75+15*i,80+15*i) for i in range(4)]]
        for a,axis in enumerate(np.eye(3)):
            v=np.zeros((135,3));v[ridx]=axis;R[:,a]=v.reshape(-1);R[:,a+3]=np.cross(axis,X).reshape(-1)
        for face in range(4):
            nodes=grid['connect'][c,15*face:15*(face+1)];key=int(nodes[0]);block=R[225+45*face:225+45*(face+1)]
            if key in seen:
                old,oldblock=seen[key];row=np.zeros((45,6*nc));row[:,6*c:6*c+6]=block;row[:,6*old:6*old+6]=-oldblock;rows.append(row)
            else:seen[key]=(c,block)
            if np.all(np.isin(nodes,grid['fixed'])):
                row=np.zeros((45,6*nc));row[:,6*c:6*c+6]=block;supports.append(row)
    G=np.vstack(rows);Gs=np.vstack(rows+supports);r=f.ranks(np.linalg.svd(G,compute_uv=False));rs=f.ranks(np.linalg.svd(Gs,compute_uv=False))
    return dict(free_ranks=r,supported_ranks=rs,expected_free=6*nc-6,expected_supported=6*nc,passed=r==[6*nc-6]*3 and rs==[6*nc]*3)
