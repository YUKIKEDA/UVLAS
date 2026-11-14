"""Research ANCF EL2 with objective moment/constant-Piola mortar constraints."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import itertools
import numpy as np
from scipy.linalg import svd,qr
import qem_el2_q4 as m
MB,C,CH=m.MB,m.C,m.CH

def ranks(s):
    if not len(s) or not s[0]:return [0]*3
    return [int(np.count_nonzero(s>tol*s[0])) for tol in [1e-9,1e-10,1e-11]]

def broken_mesh(k,n,amp):
    original=m.mesh(k,n,amp);nc=n*n;fields=[];xy=[];ids=[];free=[];connect=[];offset=0
    for f in range(3):
        old=np.concatenate([cell[m.SLICES[f]] for cell in original['cells']])
        fields.append(original['ref'][old]);xy.append(original['xy'][f][old-original['field_ids'][f][0]])
        ids.append(np.arange(offset,offset+len(old)));free.append(ids[-1][abs(xy[-1][:,0])>1e-12])
        connect.append(np.arange(offset,offset+len(old)).reshape(nc,25));offset+=len(old)
    cells=np.array([np.concatenate([connect[f][c] for f in range(3)]) for c in range(nc)])
    faces=[]
    for j,i in itertools.product(range(n),repeat=2):
        c=j*n+i
        if i<n-1:faces.append((c,c+1,0))
        if j<n-1:faces.append((c,c+n,1))
    return dict(ref=np.concatenate(fields),xy=xy,field_ids=ids,free=free,cells=cells,faces=faces,n=n)

def constraints(mesh,h,degree,order=9,order_z=5,independent=False):
    t,wt=np.polynomial.legendre.leggauss(order);z,wz=np.polynomial.legendre.leggauss(order_z)
    total=len(mesh['ref']);rows=[];metadata=[]
    for face,(left,right,axis) in enumerate(mesh['faces']):
        pl=np.zeros((order,2));pr=pl.copy();pl[:,axis]=1;pr[:,axis]=-1;pl[:,1-axis]=pr[:,1-axis]=t
        for f in range(3):
            Nl,_=m.shape(pl,f,independent);Nr,_=m.shape(pr,f,independent)
            for power in range(degree+1):
                row=np.zeros(total);row[mesh['cells'][left][m.SLICES[f]]]=np.einsum('p,pi->i',wt*t**power,Nl)
                row[mesh['cells'][right][m.SLICES[f]]]-=np.einsum('p,pi->i',wt*t**power,Nr)
                rows.append(row);metadata.append(dict(face=face,kind='moment',field=f,power=power))
        points=np.array([[p[0],p[1],zz] for p,zz in itertools.product(pl,z)])
        points_r=points.copy();points_r[:,axis]=-1
        w=np.array([a*b for a,b in itertools.product(wt,wz)])
        q=mesh['ref'][mesh['cells'][left]].reshape(225)
        phi,J,_=m.raw(q,points,h*points[:,2]/2,independent)
        phir,Jr,_=m.raw(mesh['ref'][mesh['cells'][right]].reshape(225),points_r,h*points_r[:,2]/2,independent)
        # Physical cofactor with h/2 canceled analytically against /h.
        area=np.cross(J[:,:,(axis+1)%3],J[:,:,(axis+2)%3])/2
        for j in range(3):
            row=np.zeros(total);row[mesh['cells'][left]]=np.einsum('p,pn->n',w*area[:,j],phi)
            row[mesh['cells'][right]]-=np.einsum('p,pn->n',w*area[:,j],phir)
            rows.append(row);metadata.append(dict(face=face,kind='piola',component=j))
    return np.array(rows),metadata

def null_rows(G,independent=False):
    if not G.size:return np.eye(G.shape[1]),dict(rank=0,ranks=[0]*3,discarded=0,span=0.)
    norms=np.linalg.norm(G,axis=1);largest=norms.max()
    active=norms>1e-12*largest if largest else np.zeros(len(G),dtype=bool)
    if not np.any(active):return np.eye(G.shape[1]),dict(rank=0,ranks=[0]*3,discarded=len(G),span=0.)
    scaled=G[active]/norms[active,None]
    U,s,V=svd(scaled,full_matrices=True,lapack_driver='gesvd' if independent else 'gesdd',check_finite=False)
    rr=ranks(s)
    if rr!=[rr[1]]*3:raise ValueError('mortar_constraint_rank_uncertain:'+str(rr))
    r=rr[1];Z=V[r:].T
    span=float(np.linalg.norm(G@Z)/max(1e-30,np.linalg.norm(G)))
    if span>1e-10:raise ValueError('mortar_constraint_span_failed:'+str(span))
    return Z,dict(rank=r,ranks=rr,discarded=int(np.count_nonzero(~active)),span=span)

def rigid_scan(mesh,h,degree,independent=False):
    G,metadata=constraints(mesh,h,degree,11 if independent else 9,7 if independent else 5,independent)
    R=np.zeros((3*len(mesh['ref']),6*len(mesh['cells'])))
    for c,ids in enumerate(mesh['cells']):
        q=mesh['ref'][ids];center=q[:25].mean(axis=0)
        for a,axis in enumerate(np.eye(3)):
            u=np.zeros((75,3));u[:25]=axis;R[(ids[:,None]*3+np.arange(3)).reshape(-1),6*c+a]=u.reshape(-1)
            u=np.cross(axis,q);u[:25]-=np.cross(axis,center)
            R[(ids[:,None]*3+np.arange(3)).reshape(-1),6*c+3+a]=u.reshape(-1)
    _,meta=null_rows(G,independent)
    # Use all raw rows normalized; an independent SVD checks rigid connectivity.
    norms=np.linalg.norm(G,axis=1);active=norms>1e-12*norms.max();B=np.kron(G[active]/norms[active,None],np.eye(3))@R
    scales=1/np.linalg.norm(R,axis=0);B*=scales
    singular=svd(B,compute_uv=False,lapack_driver='gesvd' if independent else 'gesdd',check_finite=False);rr=ranks(singular)
    fixed=np.setdiff1d(np.arange(len(mesh['ref'])),np.concatenate(mesh['free']));fixed_dof=(3*fixed[:,None]+np.arange(3)).reshape(-1)
    supported=np.vstack((B,R[fixed_dof]*scales));rs=ranks(svd(supported,compute_uv=False,lapack_driver='gesvd' if independent else 'gesdd',check_finite=False))
    return dict(h=h,degree=degree,ranks=rr,supported_ranks=rs,expected_rank=R.shape[1]-6,expected_supported_rank=R.shape[1],
        spurious_free_modes=R.shape[1]-rr[1]-6,spurious_supported_modes=R.shape[1]-rs[1],constraint_metadata=meta,
        passed=rr==[R.shape[1]-6]*3 and rs==[R.shape[1]]*3)

def raw_nonlinear(q,current,h,v1,v2,order=9,order_z=5,independent=False):
    plane,wp=m.quadrature(order,2);z,wz=np.polynomial.legendre.leggauss(order_z)
    points=np.array([[s,t,c] for (s,t),c in itertools.product(plane,z)]);w=np.array([a*b for a,b in itertools.product(wp,wz)])
    _,g,_=m.tables(points,h*points[:,2]/2,independent);J=np.einsum('na,pnb->pab',q.reshape(75,3),g);G=g@np.linalg.inv(J)
    mu=w*np.linalg.det(J)*h/2;F=np.eye(3)+np.einsum('na,pnb->pab',(current-q).reshape(75,3),G)
    E=np.einsum('kij,pij->pk',MB,(F.transpose(0,2,1)@F-np.eye(3))/2)
    B=np.einsum('kij,pai,pnj->pkna',MB,F,G).reshape(-1,6,225);stress=E@C.T
    d1=np.einsum('na,pnb->pab',v1.reshape(75,3),G);d2=np.einsum('na,pnb->pab',v2.reshape(75,3),G)
    second=np.einsum('kij,pij->pk',MB,(d1.transpose(0,2,1)@d2+d2.transpose(0,2,1)@d1)/2)
    energy=np.einsum('p,pk,pk->',mu,E,stress)/2;force=np.einsum('p,pki,pk->i',mu,B,stress)
    tangent=np.einsum('p,pk,kl,pl->',mu,B@v1,C,B@v2)+np.einsum('p,pk,pk->',mu,stress,second)
    return energy,force,tangent

def boundary(q,h,P,faces,order=9,independent=False):
    plane,w=m.quadrature(order,2);force=np.zeros((75,3))
    for axis,sign in faces:
        points=np.zeros((len(plane),3));points[:,axis]=sign;points[:,[i for i in range(3) if i!=axis]]=plane
        phi,J,_=m.raw(q,points,h*points[:,2]/2,independent);J=J.copy();J[:,:,2]*=h/2
        area=sign*np.cross(J[:,:,(axis+1)%3],J[:,:,(axis+2)%3]);traction=area@P.T
        force+=np.einsum('p,pn,pa->na',w,phi,traction)
    return force

def patch(mesh,h,degree,F,independent=False):
    G,meta=constraints(mesh,h,degree,11 if independent else 9,7 if independent else 5,independent)
    E=np.einsum('kij,ij->k',MB,(F.T@F-np.eye(3))/2);S=np.einsum('kij,k->ij',MB,C@E);P=F@S
    internal=np.zeros_like(mesh['ref']);outer=np.zeros_like(internal)
    for c,ids in enumerate(mesh['cells']):
        q=mesh['ref'][ids].reshape(225);current=(mesh['ref'][ids]@F.T).reshape(225)
        _,f,_=raw_nonlinear(q,current,h,np.zeros(225),np.zeros(225),11 if independent else 9,7 if independent else 5,independent);internal[ids]+=f.reshape(75,3)
        j,i=divmod(c,mesh['n']);faces=[(2,-1),(2,1)]
        if i==0:faces.append((0,-1))
        if i==mesh['n']-1:faces.append((0,1))
        if j==0:faces.append((1,-1))
        if j==mesh['n']-1:faces.append((1,1))
        outer[ids]+=boundary(q,h,P,faces,11 if independent else 9,independent)
    multipliers=np.zeros((len(G),3))
    for a,entry in enumerate(meta):
        if entry['kind']=='piola':multipliers[a]=-h*P[:,entry['component']]
    reaction=G.T@multipliers
    denominator=max(1e-30,np.linalg.norm(outer));residual=float(np.linalg.norm(internal+reaction-outer)/denominator)
    gap=float(np.linalg.norm(G@(mesh['ref']@F.T))/max(1.,np.linalg.norm(mesh['ref'])))
    return dict(residual=residual,gap=gap,force_norm=float(np.linalg.norm(internal)),reaction_norm=float(np.linalg.norm(reaction)))

def shared(mesh,k,degree,order=9,independent=False):
    a=m.assemble(mesh,k,order,'raw');G,meta=constraints(mesh,0.,degree,order,7 if independent else 5,independent)
    Zrd,rdmeta=null_rows(G[:,a['free_rd']],independent);Zc,cmeta=null_rows(G[:,a['free_c']],independent)
    Zrd=np.kron(Zrd,np.eye(3));Zc=np.kron(Zc,np.eye(3))
    return a,Zrd,Zc,dict(rd=rdmeta,c=cmeta)
