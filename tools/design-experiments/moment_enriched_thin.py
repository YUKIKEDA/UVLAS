"""All-DOF leading thin operator for the fixed RW-3 discretization."""
import numpy as np
from scipy.linalg import svd,qr,lstsq
import reference_work_fit as fit
import moment_enriched_ancf as f
c=f.c;m=f.m

def rank_info(s):return f.ranks(s)

def nullspace(A,normalize=False,independent=False):
    scale=1/np.maximum(np.linalg.norm(A,axis=0),1e-30) if normalize else np.ones(A.shape[1])
    B=A*scale
    # Compress tall constraints with an orthogonal row transform, retaining every column.
    if len(B)>B.shape[1]:B=f.qrtools.rfactor(np.asfortranarray(B))
    U,s,V=svd(B,full_matrices=True,lapack_driver='gesvd' if independent else 'gesdd',check_finite=False)
    ranks=rank_info(s)
    if len(set(ranks))!=1:raise ValueError('constraint_rank_ambiguous:'+str(ranks))
    Z=scale[:,None]*V[ranks[0]:].T
    if Z.shape[1]:Z=qr(Z,mode='economic',check_finite=False)[0]
    residual=float(np.linalg.norm(A@Z)/max(1.,np.linalg.norm(A)))
    return Z,dict(ranks=ranks,columns=A.shape[1],nullity=Z.shape[1],residual=residual,smallest_nonzero=float(s[ranks[0]-1]) if ranks[0] else None)

def operators(q,order=13,zorder=7,independent=False):
    s=f.prepare(q,.01,order=order,zorder=zorder,independent=independent);p=s['p'];N,g=m.shape(p,0,independent);coord=q.reshape(75,3)
    Jr=np.einsum('na,pnb->pab',coord[:25],g);director=N@coord[25:50];normal=np.cross(Jr[:,:,0],Jr[:,:,1]);d0=np.einsum('pa,pa->p',normal,director)/2
    if d0.min()<=0:raise ValueError('nonpositive_leading_volume')
    coeff,_=c.geometry_coefficients(np.asarray(q,dtype=np.float64).tobytes(),independent);L,_=c.scalar_basis(p,independent);C=np.einsum('pm,mabj->pabj',L,coeff.reshape(len(c.IDS),3,3,5),optimize=True)
    A=[]
    for power in range(3):
        G=np.einsum('pjn,pjk->pnk',s['Gwref'],C[:,:,:,power],optimize=True);A.append(np.einsum('kab,pnb->pkna',m.MB,G).reshape(len(p),6,3*s['nscalar']))
    return dict(state=s,p=p,w=s['w'],d0=d0,C=C,A=A)

def jump_root(op,which):
    s=op['state'];rows=[]
    for face in s['faces']:
        if which=='all' or (which=='tb')==(face['axis']==2):rows.append(np.kron(face['jump']*np.sqrt(face['w'])[:,None],np.eye(3)))
    return np.vstack(rows)

def constraints(op):
    s=op['state'];nc=3*s['nscalar'];w=np.sqrt(op['w'])[:,None,None];a0=(op['A'][0]*w).reshape(-1,nc);a1=(op['A'][1]*w).reshape(-1,nc)
    allj=jump_root(op,'all');tb=jump_root(op,'tb');C=np.zeros((2*len(a0)+len(allj)+len(tb),2*nc));offset=0
    C[:len(a0),:nc]=a0;offset+=len(a0);C[offset:offset+len(a0),:nc]=a1;C[offset:offset+len(a0),nc:]=a0;offset+=len(a0)
    C[offset:offset+len(allj),:nc]=allj;offset+=len(allj);C[offset:,nc:]=tb
    return C

def leading_rows(op,kappa):
    s=op['state'];nc=3*s['nscalar'];weighted=[]
    for A in op['A']:weighted.append((np.einsum('ij,pjn->pin',m.CH,A)*np.sqrt(op['w']/op['d0'])[:,None,None]).reshape(-1,nc))
    vrows=[]
    for face in s['faces'][2:]:
        Cf,_=c.geometry_coefficients(np.asarray(s['q'],dtype=np.float64).tobytes());Lf,_=c.scalar_basis(face['points']);C1=np.einsum('pm,mab->pab',Lf,Cf.reshape(len(c.IDS),3,3,5)[:,:,:,1],optimize=True)
        area=np.linalg.norm(C1[:,face['axis'],:],axis=1);vrows.append(np.kron(face['jump']*np.sqrt(s['tau']*face['w']*area)[:,None],np.eye(3)))
    vertical=np.vstack(vrows);nrows=len(weighted[0])+len(vertical);roots=[]
    for j in range(3):
        M=np.zeros((nrows,nc));M[:len(weighted[0])]=weighted[2-j]
        if j==1:M[len(weighted[0]):]=vertical
        roots.append(M)
    N,_=m.shape(op['p'],0);x=(N@s['q'].reshape(75,3)[:25])[:,0];a=np.column_stack((np.ones(len(x)),np.zeros(len(x)),kappa*x));targets=np.zeros((nrows,2))
    for j,kind in enumerate(fit.KIND_TARGETS):
        curvature=x*x if kind=='quartic' else np.exp(x);E=-op['p'][:,2,None,None]/2*curvature[:,None,None]*a[:,:,None]*a[:,None,:]/(1+kappa*kappa*x*x)[:,None,None]
        mandel=np.einsum('kij,pij->pk',m.MB,E);targets[:len(weighted[0]),j]=(mandel@m.CH.T*np.sqrt(op['w']*op['d0'])[:,None]).reshape(-1)
    return roots,targets

def local(q,kappa):
    op=operators(q);s=op['state'];nc=3*s['nscalar'];ni=3*s['ni']
    from moment_enriched_exact_isometry import local_basis
    Z,info=local_basis(op,kappa)
    roots,b=leading_rows(op,kappa);W01=np.column_stack(roots[:2])@Z;Pi=even_internal_basis(s);Wi=roots[2][:,:ni]@Pi;Wb=roots[2][:,ni:]
    U,sv,V=svd(Wi,full_matrices=False,check_finite=False);ranks=rank_info(sv)
    if len(set(ranks))!=1:raise ValueError('q2_internal_rank_ambiguous:'+str(ranks))
    Q=U[:,:ranks[0]];combined=np.column_stack((W01,Wb,b));combined-=Q@(Q.T@combined);R=f.qrtools.rfactor(np.asfortranarray(combined));del combined,U,Q
    return dict(q=q,Z=Z,info=info,nc=nc,ni=ni,R=R,ny=Z.shape[1],nb=nc-ni,faces=[dict(axis=x['axis'],sign=x['sign'],offset=x['offset'],n=x['n'],Mcoeff=x['Mcoeff']) for x in s['faces'][2:]],
        reference=np.sum(b*b,axis=0),q2_internal_ranks=ranks,nscalar=s['nscalar'])

def recover(cell,y,q2b,kappa):
    op=operators(cell['q']);roots,b=leading_rows(op,kappa);q01=cell['Z']@y;nc=cell['nc'];ni=cell['ni'];rhs=b-roots[0]@q01[:nc]-roots[1]@q01[nc:]-roots[2][:,ni:]@q2b
    Pi=even_internal_basis(op['state']);q2c,_,rank,sv=lstsq(roots[2][:,:ni]@Pi,rhs,cond=1e-10,lapack_driver='gelsd',check_finite=False);q2i=Pi@q2c
    if rank_info(sv)!=cell['q2_internal_ranks']:raise ValueError('q2_recovery_rank')
    return np.stack((q01[:nc],q01[nc:],np.vstack((q2i,q2b))))

def even_internal_basis(s):
    columns=[];ni=s['ni']
    for j,key in enumerate(f.VID):
        if key[2]%2==0:v=np.zeros(ni);v[j]=1;columns.append(v)
    bottom,top=s['faces'][:2]
    for j in range(25):
        v=np.zeros(ni);v[bottom['offset']+j]=1/np.sqrt(2);v[top['offset']+j]=1/np.sqrt(2);columns.append(v)
    return np.kron(np.column_stack(columns),np.eye(3))
