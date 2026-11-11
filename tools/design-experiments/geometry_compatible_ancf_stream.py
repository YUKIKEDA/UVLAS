"""Same GC-ANCF space and quadrature; bounded-memory point and QR batches."""
import numpy as np
from scipy.linalg import svd,solve_triangular
import geometry_compatible_ancf_trial as previous
c=previous.c;old=previous.old;m=previous.m;f=previous.f
bending=previous.bending
BATCH=1024
def chunks(p,w):
    for offset in range(0,len(p),BATCH):yield p[offset:offset+BATCH],w[offset:offset+BATCH]
def shape_values(a,cell,q,p):
    fs=c.fields(a['n'],cell,q,p[:,:2],a['geo']);z=a['h']*p[:,2]/2
    return np.column_stack((fs[0][0],z[:,None]*fs[1][0],z[:,None]**2*fs[2][0]/2))
def initialize(k,n,amp,h,order=13,zorder=7,geo=None):
    geo=c.prepare_geometry(k,n,amp,order) if geo is None else geo;fs=c.fields(n,0,geo['refs'][0],np.zeros((1,2)),geo);nr,nd,nc=[v.shape[1] for v,gx,gy in fs];ns=nr+nd+nc
    a=dict(geo=geo,n=n,k=k,amp=amp,h=h,order=order,zorder=zorder,nr=nr,nd=nd,nc=nc,ns=ns);constraints=[];mass=np.zeros((ns,ns));reference=np.zeros(3)
    for cell,q in enumerate(geo['refs']):
        allp,allw=c.volume_rule(q,n,order,zorder)
        for p,w in chunks(allp,allw):
            N=shape_values(a,cell,q,p);_,_,det,_,_=f.base.geometry(q,h,p);E,_=old.reference(q,p,h,k)
            mass+=N.T@((w*det)[:,None]*N);reference+=np.einsum('p,pkt,kl,plt->t',w*det,E,m.C,E)/h**3
        if cell%n==0:
            edge=np.array([[-1.,t,z] for t in np.linspace(-1,1,33) for z in [-1.,0.,1.]])
            constraints.append(c.kinematics(n,cell,q,h,edge,geo)[0])
    C=np.vstack(constraints);_,s,V=svd(C,full_matrices=True,check_finite=False);ranks=[int(np.count_nonzero(s>tol)) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError('GC_boundary_rank:'+str(ranks))
    Zs=V[ranks[0]:].T;ms=1/np.sqrt(np.diag(mass));a.update(Zs=Zs,Z=np.kron(Zs,np.eye(3)),boundary_ranks=ranks,mass_scaled_min=float(np.linalg.eigvalsh((mass*ms)*ms[:,None])[0]),reference_norm=reference)
    return a
def anchor(a):
    R=None
    for q in a['geo']['refs']:
        allp,allw=c.volume_rule(q,a['n'],a['order'],a['zorder'])
        for p,w in chunks(allp,allw):
            _,_,det,_,_=f.base.geometry(q,a['h'],p);E,_=bending(q,p,a['h'],a['k'],a['n']);target,_=old.reference(q,p,a['h'],a['k']);weight=np.sqrt(w*det)[:,None,None]/a['h']**1.5
            A=(np.einsum('ij,pjn->pin',m.CH,E)*weight).reshape(-1,E.shape[-1]);b=(np.einsum('ij,pjt->pit',m.CH,target)*weight).reshape(-1,3);block=np.column_stack((A,b))
            R=old.qrtools.rfactor(np.asfortranarray(block if R is None else np.vstack((R,block))))
    n=R.shape[1]-3;scale=1/np.linalg.norm(R[:,:n],axis=0);R[:,:n]*=scale;R=old.qrtools.rfactor(np.asfortranarray(R));sv=svd(R[:n,:n],compute_uv=False);ranks=[int(np.count_nonzero(sv>tol*sv[0])) for tol in [1e-9,1e-10,1e-11]]
    a['anchor_ranks']=ranks;a['anchor_columns']=n
    if ranks!=[n]*3:raise ValueError('GC_anchor_rank:'+str(ranks))
    return scale[:,None]*solve_triangular(R[:n,:n],R[:n,n:],check_finite=False)
def residual_system(a,beta):
    rows=[];rhs=[]
    for cell,q in enumerate(a['geo']['refs']):
        allp,allw=c.volume_rule(q,a['n'],a['order'],a['zorder']);R=None
        for p,w in chunks(allp,allw):
            N,G,B,det,_=c.kinematics(a['n'],cell,q,a['h'],p,a['geo']);E,_=bending(q,p,a['h'],a['k'],a['n']);target,_=old.reference(q,p,a['h'],a['k']);target-=E@beta;weight=np.sqrt(w*det)[:,None,None]/a['h']**1.5
            A=(np.einsum('ij,pjn->pin',m.CH,B)*weight).reshape(-1,3*a['ns']);b=(np.einsum('ij,pjt->pit',m.CH,target)*weight).reshape(-1,3);block=np.column_stack((A,b))
            R=old.qrtools.rfactor(np.asfortranarray(block if R is None else np.vstack((R,block))))
        rows.append(R[:,:3*a['ns']]);rhs.append(R[:,3*a['ns']:])
    a['A']=np.vstack(rows);a['b']=np.vstack(rhs)
def solve(a,free_rank=False):
    beta=anchor(a);residual_system(a,beta);s=c.solve(a,free_rank)
    if s['passed_rank']:s['beta']=beta
    return s
def observe(a,sol,order,zorder):
    error=np.zeros(3);norm=np.zeros(3);pe=np.zeros(3);pn=np.zeros(3);gradient=np.zeros((3*a['ns'],3));coef=sol['coef'];beta=sol['beta']
    for cell,q in enumerate(a['geo']['refs']):
        allp,allw=c.volume_rule(q,a['n'],order,zorder)
        for p,w in chunks(allp,allw):
            N,G,B,det,_=c.kinematics(a['n'],cell,q,a['h'],p,a['geo']);E,U=bending(q,p,a['h'],a['k'],a['n']);target,pos=old.reference(q,p,a['h'],a['k']);e=B@coef+(E@beta-target);stress=np.einsum('ij,pjt->pit',m.C,e)
            error+=np.einsum('p,pkt,pkt->t',w*det,e,stress)/a['h']**3;norm+=np.einsum('p,pkt,kl,plt->t',w*det,target,m.C,target)/a['h']**3
            gradient+=np.einsum('p,pkn,pkt->nt',w*det,B,stress,optimize=True)/a['h']**3
            value=(N@coef.reshape(a['ns'],9)).reshape(len(p),3,3)+U@beta;pe+=np.einsum('p,pat,pat->t',w*det,value-pos,value-pos);pn+=np.einsum('p,pat,pat->t',w*det,pos,pos)
    scale=1/np.linalg.norm(a['A']@a['Z'],axis=0);g=scale[:,None]*(a['Z'].T@gradient);kkt=np.linalg.norm(g,axis=0)/np.maximum(1e-30,np.sqrt(len(scale))*np.sqrt(norm))
    return dict(error=np.sqrt(error/norm).tolist(),position=np.sqrt(pe/pn).tolist(),original_kkt=kkt.tolist())
